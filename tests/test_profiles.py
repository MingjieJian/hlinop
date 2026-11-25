"""Line profiles from hlinop.f: STARK1, HPROFL, HLINOP and the cutoff weight HLINOP applies (HLCUT).

A line profile must integrate to one over frequency (per Hz).  HLINPROF: test_hlinprof.py.
STARK1 (default HLINOP_SET_STARK(1)) convolves the ion quasistatic profile with the electron
profile (impact core plus quasistatic electron wing); mode 0 is the earlier code, which added
the parts.
"""
import json
from pathlib import Path

import numpy as np
import pytest

import hltest as H
from cases import profile_cases

PROFILE_CASES = list(profile_cases())
DATA = Path(__file__).resolve().parent / "data"
# n_up >= 8: the Griem-type profile is meant for the lines HLINOP is used for (low lines such as
# H beta, with its central dip, and Lyman lines always come from the tables)
CDS = [e for e in json.loads((DATA / "stehle_cds_subset.json").read_text()) if e["nlo"] > 1 and e["nup"] >= 8]
VCS = json.loads((DATA / "lemke_vcs_brackett.json").read_text())

# Doppler-dominated cores: HPROFL takes the Doppler profile inside the core and ADDS the Stark and
# Lorentz profiles outside it, which adds area (to be fixed in HPROFL; STARK1 itself is normalised)
HPROFL_COMBINATION = {(2, 5, 5000.0, 1e12), (2, 5, 6300.0, 1e14), (2, 8, 5000.0, 1e12),
                      (2, 10, 5000.0, 1e12), (2, 14, 5000.0, 1e12), (3, 10, 5000.0, 1e12)}


def _ids(cases):
    return ["%d-%d_T%d_ne%.0e" % c[:4] for c in cases]


def _area(w, phi):
    with_tail = H.area_per_frequency(w, phi, tail=True)
    no_tail = H.area_per_frequency(w, phi, tail=False)
    # precondition on the integration grid, not on the code; pytest.fail is not an
    # AssertionError, so a bad grid shows as a real failure even under an xfail mark
    if abs(with_tail - no_tail) > 0.01 * with_tail:
        pytest.fail("integration grid too narrow: far wings %.3f of %.3f" % (with_tail - no_tail, with_tail))
    return with_tail


def _hwhm(x, y):
    i0 = int(np.argmax(y))
    half = y[i0] / 2
    r = i0 + int(np.argmax(y[i0:] < half))
    return float(np.interp(half, [y[r], y[r - 1]], [x[r], x[r - 1]]) - x[i0])


@pytest.mark.parametrize("case", PROFILE_CASES, ids=_ids(PROFILE_CASES))
def test_stark1_has_unit_area(case):
    nlo, nup, T, ne, nh, nhe = case
    w0 = H.line_centre(nlo, nup)
    w = H.detuning_grid(w0, rel_max=0.25)
    s = H.profile(nlo, nup, T, ne, nh, nhe, H.doppler(T), w, stark1=True)
    assert abs(_area(w, s) - 1.0) < 0.01


def _hprofl_params():
    out = []
    for c in PROFILE_CASES:
        marks = []
        if c[:4] in HPROFL_COMBINATION:
            marks = [pytest.mark.known_issue,
                     pytest.mark.xfail(reason="HPROFL uses the Doppler profile inside a Doppler-dominated core "
                                              "and adds the Stark and Lorentz profiles outside it; areas 1.02-1.39 "
                                              "here.  Fix: convolve in HPROFL.", strict=True, raises=AssertionError)]
        out.append(pytest.param(c, marks=marks, id=_ids([c])[0]))
    return out


@pytest.mark.parametrize("case", _hprofl_params())
def test_hprofl_has_unit_area(case):
    nlo, nup, T, ne, nh, nhe = case
    w0 = H.line_centre(nlo, nup)
    w = H.detuning_grid(w0, rel_max=0.25)
    p, _ = H.profile(nlo, nup, T, ne, nh, nhe, H.doppler(T), w)
    assert abs(_area(w, p) - 1.0) < 0.01


@pytest.mark.parametrize("entry", CDS, ids=["%d-%d_ne%.0e_T%g" % (e["nlo"], e["nup"], e["ne"], e["T"]) for e in CDS])
def test_stark1_matches_stehle_hutcheon(entry):
    """STARK1 against the Stark-only profiles of Stehle & Hutcheon (1999): HWHM and peak within
    10% (a benchmark over 171 Balmer and Paschen profiles gives medians of 2-5% and 95th
    percentiles of about 10%)."""
    nlo, nup, T, ne, f0 = entry["nlo"], entry["nup"], entry["T"], entry["ne"], entry["f0"]
    w0 = H.line_centre(nlo, nup)
    hw_ref = entry["hwhm_stark_alpha"] * f0
    x = np.linspace(0.0, 20.0 * hw_ref, 4001)
    s = H.profile(nlo, nup, T, ne, 0.0, 0.0, H.doppler(T), w0 + x, stark1=True) * H.C_ANGSTROM / (w0 + x) ** 2
    peak_ref = entry["sstark"][0] / f0
    assert abs(np.log(_hwhm(x, s) / hw_ref)) < np.log(1.10)
    assert abs(np.log(s[0] / peak_ref)) < np.log(1.10)


@pytest.mark.parametrize("ref", VCS, ids=["Br%d_T%g_ne%.0e" % (r["nup"], r["T"], r["ne"]) for r in VCS])
def test_stark1_matches_lemke_vcs_for_brackett(ref):
    """Brackett lines (no Stehle tables) against Lemke (1997) VCS Stark-only profiles: HWHM and
    peak within 15% (benchmark: medians 3%)."""
    nup, T, ne = ref["nup"], ref["T"], ref["ne"]
    w0 = H.line_centre(4, nup)
    x = np.linspace(0.0, 20.0 * ref["hwhm_A"], 4001)
    s = H.profile(4, nup, T, ne, 0.0, 0.0, H.doppler(T), w0 + x, stark1=True) * H.C_ANGSTROM / (w0 + x) ** 2
    assert abs(np.log(_hwhm(x, s) / ref["hwhm_A"])) < np.log(1.15)
    assert abs(np.log(s[0] / ref["peak_per_A"])) < np.log(1.15)


def test_lyman_alpha_keeps_the_earlier_stark_profile():
    """Lyman alpha uses the earlier construction in both modes (to be revisited)."""
    w0 = H.line_centre(1, 2)
    w = w0 + np.array([0.0, 0.1, 1.0, 10.0, -1.0, -10.0])
    a = H.profile(1, 2, 8000.0, 1e15, 1e17, 1e16, H.doppler(8000.0), w, stark1=True)
    b = H.profile(1, 2, 8000.0, 1e15, 1e17, 1e16, H.doppler(8000.0), w, stark1=True, stark=0)
    np.testing.assert_array_equal(a, b)


def test_hlinop_equals_hprofl_inside_the_cutoffs():
    """HLINOP is HPROFL with cutoffs at the series limit (blue) and a red cutoff."""
    nlo, nup, T, ne, nh, nhe = 2, 10, 6300.0, 1e14, 1e17, 1e16
    w0 = H.line_centre(nlo, nup)
    w = w0 + np.linspace(-20.0, 20.0, 81)
    p, q = H.profile(nlo, nup, T, ne, nh, nhe, H.doppler(T), w)
    np.testing.assert_allclose(q, p, rtol=1e-6)


@pytest.mark.parametrize("line", [(2, 10), (2, 3), (1, 3), (3, 4), (3, 12), (4, 5), (4, 11)],
                         ids=lambda l: "%d-%d" % l)
def test_hlinop_is_hprofl_times_hlcut(line):
    """HLCUT reproduces the cutoffs of HLINOP (blue cut at the series limit, red cut with the
    taper for alpha lines), so that HBOP can apply them to tabulated profiles."""
    nlo, nup = line
    T, ne, nh, nhe = 6300.0, 1e14, 1e17, 1e16
    w0 = H.line_centre(nlo, nup)
    limit = 1e8 / (H.IONH_CM - H.level_energy(nlo))
    w = np.concatenate([w0 + np.linspace(-30.0, 30.0, 61), limit + np.linspace(-20.0, 5.0, 26),
                        w0 * np.array([1.2, 1.5, 2.0])])
    p, q, c = H.profile(nlo, nup, T, ne, nh, nhe, H.doppler(T), w, with_cut=True)
    np.testing.assert_allclose(q, p * c, rtol=1e-6, atol=1e-35)
    assert np.any(c == 0.0)
