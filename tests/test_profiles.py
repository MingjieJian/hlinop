"""Line profiles: HPROFL/HLINOP (hlinop.f) and HLINPROF (hlinprof.f).

A line profile must integrate to one.  HPROFL integrates over frequency (per Hz),
HLINPROF over wavelength (per Angstrom).
"""
import numpy as np
import pytest

import hltest as H
from cases import HLINPROF_CASES, HLINPROF_VCS_CONDITIONS, profile_cases

PROFILE_CASES = list(profile_cases())


def _ids(cases):
    return ["%d-%d_T%d_ne%.0e" % c[:4] for c in cases]


def _hprofl_area(nlo, nup, T, ne, nh, nhe):
    w0 = H.line_centre(nlo, nup)
    w = H.detuning_grid(w0)
    p, _ = H.profile(nlo, nup, T, ne, nh, nhe, H.doppler(T), w)
    with_tail = H.area_per_frequency(w, p, tail=True)
    no_tail = H.area_per_frequency(w, p, tail=False)
    return with_tail, no_tail


@pytest.mark.known_issue
@pytest.mark.xfail(reason="STARK1 adds the ion quasistatic, electron quasistatic and electron impact parts "
                          "instead of convolving them; areas are 1.0-2.1 (1.74-2.14 wherever Stark broadening "
                          "dominates).  Fix: convolve the parts.", strict=True, raises=AssertionError)
@pytest.mark.parametrize("case", PROFILE_CASES, ids=_ids(PROFILE_CASES))
def test_hprofl_has_unit_area(case):
    with_tail, no_tail = _hprofl_area(*case)
    # precondition on the integration grid, not on the code: the estimated far wings beyond the
    # grid must be < 1% of the area.  pytest.fail is not an AssertionError, so a bad grid is
    # reported as a real failure even under the xfail mark.
    if abs(with_tail - no_tail) > 0.01 * with_tail:
        pytest.fail("integration grid too narrow: far wings %.3f of %.3f" % (with_tail - no_tail, with_tail))
    assert abs(with_tail - 1.0) < 0.01


@pytest.mark.parametrize("case", HLINPROF_CASES, ids=["2-%d_T%d_ne%.0e" % c[1:4] for c in HLINPROF_CASES])
def test_hlinprof_has_unit_area(case):
    """Stehle tables plus self, helium and radiative broadening.  Measured: 0.97-1.04.
    The self-broadening Lorentzian is added (not convolved) beyond 0.2 A from the centre, which
    adds up to 4% (H gamma, n_HI = 1e17); the Stehle-only area is 0.97-0.99 within +-5% of the
    line wavelength.  Tolerance 5%."""
    nlo, nup, T, ne, nh = case
    w0 = H.line_centre(nlo, nup)
    w = H.detuning_grid(w0, rel_max=0.05, npts=2500)
    for h1, he1 in ((0.0, 0.0), (nh, nh / 10.0)):
        a = H.area_per_angstrom(w, H.hlinprof(nlo, nup, T, ne, h1, he1, H.doppler(T), w))
        assert abs(a - 1.0) < 0.05, (h1, a)


@pytest.mark.parametrize("cond", HLINPROF_VCS_CONDITIONS,
                         ids=["T%d_ne%.0e" % c[:2] for c in HLINPROF_VCS_CONDITIONS])
def test_hlinprof_vcs_fallback_has_unit_area(cond):
    """Outside the Stehle grid of this repository (T 2500-10000 K, log n_e 10-16) HLINPROF uses
    the VCS tables (vcsbalmer.dat).  Measured: 0.96-1.02 for H alpha - H gamma.  Tolerance 5%."""
    T, ne, nh = cond
    for nup in (3, 4, 5):
        w0 = H.line_centre(2, nup)
        w = H.detuning_grid(w0, rel_max=0.05, npts=2500)
        for h1 in (0.0, nh):
            p = H.hlinprof(2, nup, T, ne, h1, h1 / 10.0, H.doppler(T), w)
            assert np.all(np.isfinite(p)) and np.all(p >= 0)
            a = H.area_per_angstrom(w, p)
            assert abs(a - 1.0) < 0.05, (nup, h1, a)


def test_hlinop_equals_hprofl_inside_the_cutoffs():
    """HLINOP is HPROFL with cutoffs at the series limit (blue) and a red cutoff."""
    nlo, nup, T, ne, nh, nhe = 2, 10, 6300.0, 1e14, 1e17, 1e16
    w0 = H.line_centre(nlo, nup)
    w = w0 + np.linspace(-20.0, 20.0, 81)
    p, q = H.profile(nlo, nup, T, ne, nh, nhe, H.doppler(T), w)
    np.testing.assert_allclose(q, p, rtol=1e-6)
