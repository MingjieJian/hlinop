"""HLINPROF (hlinprof.f) with the Stehle & Hutcheon (1999) table built from CDS (stehle_long.dat).

HLINPROF returns a profile per Angstrom.  Table: Lyman, Balmer and Paschen to n_up = 30,
T 2500 K - 1.259e6 K, log n_e 10 - 19.5 (fewer densities for higher n_up).
"""
import json
from pathlib import Path

import numpy as np
import pytest

import hltest as H
from cases import HLINPROF_CASES

CDS = json.loads((Path(__file__).resolve().parent / "data" / "stehle_cds_subset.json").read_text())
C_A = H.C_ANGSTROM


def _area(nlo, nup, T, ne, h1, he1, rel=0.1):
    """Area over +-rel of the line wavelength and the estimated far wings beyond (alpha^-2.5)."""
    w0 = H.line_centre(nlo, nup)
    w = H.detuning_grid(w0, rel_max=rel, npts=3000)
    p = H.hlinprof(nlo, nup, T, ne, h1, he1, H.doppler(T), w)
    assert np.all(np.isfinite(p)) and np.all(p >= 0)
    o = np.argsort(w)
    w, p = w[o], p[o]
    a = H.area_per_angstrom(w, p)
    tail = 2.0 / 3.0 * (p[0] * abs(w[0] - w0) + p[-1] * abs(w[-1] - w0))
    return a, tail


@pytest.mark.parametrize("entry", CDS, ids=["%d-%d_ne%.0e_T%g" % (e["nlo"], e["nup"], e["ne"], e["T"]) for e in CDS])
def test_hlinprof_reproduces_cds_at_table_nodes(entry):
    """At a tabulated temperature and density HLINPROF returns the CDS profile, S(alpha)/F0 with
    delta lambda = alpha * F0 (no self or helium broadening; points where the added radiative
    Lorentzian could matter at the 1e-4 level are skipped)."""
    nlo, nup, T, ne, f0 = entry["nlo"], entry["nup"], entry["T"], entry["ne"], entry["f0"]
    w0 = H.line_centre(nlo, nup)
    a = np.array(entry["alpha"])
    ref = np.array(entry["sdop"]) / f0
    dl = a * f0
    lor = 1e9 * w0 ** 2 / C_A / (2 * np.pi)                 # upper bound of the radiative half width (A)
    use = (dl <= 0.2) | (lor / (np.pi * np.maximum(dl, 1e-30) ** 2) < 1e-4 * ref)
    got = H.hlinprof(nlo, nup, T, ne, 0.0, 0.0, H.doppler(T), w0 + dl[use], wave0=w0)
    np.testing.assert_allclose(got, ref[use], rtol=2e-4)
    assert use.sum() >= 5


@pytest.mark.parametrize("case", HLINPROF_CASES, ids=["%d-%d_T%d_ne%.0e" % c[:4] for c in HLINPROF_CASES])
def test_hlinprof_has_unit_area(case):
    """Measured 0.955-1.004.  The published CDS profiles themselves integrate to 0.97-0.99 over
    their tabulated range (direct integration of the files), so HLINPROF follows the tables; the
    self-broadening Lorentzian added beyond 0.2 A adds up to a few per cent.  Tolerance 5%."""
    nlo, nup, T, ne, nh = case
    for h1, he1 in ((0.0, 0.0), (nh, nh / 10.0)):
        a, tail = _area(nlo, nup, T, ne, h1, he1)
        if tail > 0.01:
            pytest.fail("integration window too narrow: far wings %.3f" % tail)
        assert abs(a + tail - 1.0) < 0.05, (h1, a + tail)


@pytest.mark.parametrize("line,ne_last", [((2, 22), 1e13), ((2, 30), 10 ** 11.5), ((3, 14), 1e14)],
                         ids=["2-22", "2-30", "3-14"])
def test_density_above_the_table_uses_quasistatic_scaling(line, ne_last):
    """Above a line's last tabulated density the profile of that density is used in alpha units
    with F0 ~ n_e^(2/3): phi(dl; n_e) = phi(dl / s; n_last) / s with s = (n_e/n_last)^(2/3).
    Tolerance 1%: HLINPROF adds a radiative Lorentzian beyond 0.2 A that does not scale with
    density (0.1-0.5% of the profile at these detunings)."""
    nlo, nup = line
    T = 6300.0
    w0 = H.line_centre(nlo, nup)
    dl = np.array([0.0, 1.0, 3.0, 10.0, 30.0, -1.0, -10.0])
    base = H.hlinprof(nlo, nup, T, ne_last, 0.0, 0.0, H.doppler(T), w0 + dl, wave0=w0)
    for fac in (3.0, 30.0):
        s = fac ** (2.0 / 3.0)
        got = H.hlinprof(nlo, nup, T, ne_last * fac, 0.0, 0.0, H.doppler(T), w0 + dl * s, wave0=w0)
        np.testing.assert_allclose(got, base / s, rtol=1e-2)
    # no jump at the boundary
    just = H.hlinprof(nlo, nup, T, ne_last * 1.0001, 0.0, 0.0, H.doppler(T), w0 + dl, wave0=w0)
    np.testing.assert_allclose(just, base, rtol=1e-3)


def test_density_below_the_table_keeps_the_first_profile():
    """Below n_e = 1e10 the first tabulated profile is used in wavelength units."""
    nlo, nup, T = 2, 10, 5000.0
    w0 = H.line_centre(nlo, nup)
    w = w0 + np.array([0.0, 0.05, 0.2, 1.0, 5.0, -0.2])
    a = H.hlinprof(nlo, nup, T, 1e10, 0.0, 0.0, H.doppler(T), w, wave0=w0)
    b = H.hlinprof(nlo, nup, T, 1e8, 0.0, 0.0, H.doppler(T), w, wave0=w0)
    np.testing.assert_allclose(b, a, rtol=1e-6)


def test_missing_low_temperatures_at_high_density():
    """At high densities the tables omit the lowest temperatures; there the last density
    tabulated at both bracketing temperatures is scaled.  Result finite and normalised."""
    a, tail = _area(1, 3, 3000.0, 1e17, 0.0, 0.0)
    assert tail < 0.01 and abs(a + tail - 1.0) < 0.05


@pytest.mark.parametrize("line,T,ne", [((1, 2), 6000.0, 1e14), ((4, 11), 6000.0, 1e14), ((2, 31), 6000.0, 1e12),
                                       ((2, 10), 2000.0, 1e12)],
                         ids=["Lyman-alpha", "Brackett", "n_up>30", "T below table"])
def test_hlinop_is_used_where_there_is_no_table(line, T, ne):
    """Lyman alpha (by choice), lines without tables, and temperatures outside the table go to
    HLINOP: HLINPROF = HLINOP * c / lambda^2."""
    nlo, nup = line
    nh, nhe = 1e16, 1e15
    w0 = H.line_centre(nlo, nup)
    w = w0 + np.array([0.0, 0.5, 2.0, -2.0, 10.0])
    got = H.hlinprof(nlo, nup, T, ne, nh, nhe, H.doppler(T), w, wave0=w0)
    _, q = H.profile(nlo, nup, T, ne, nh, nhe, H.doppler(T), w, wave0=w0)
    np.testing.assert_allclose(got, q * C_A / w ** 2, rtol=1e-5)
