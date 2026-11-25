"""Line profiles from hlinop.f: HPROFL, HLINOP and the cutoff weight HLINOP applies (HLCUT).

A line profile must integrate to one over frequency (per Hz).  HLINPROF: test_hlinprof.py.
"""
import numpy as np
import pytest

import hltest as H
from cases import profile_cases

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
