"""Occupation probability WCALC (hbop.f).

The reference is an independent Python transcription of the published formulas:
- neutral perturbers: Hummer & Mihalas (1988) eq. 4.71, hard spheres with radii sqrt(<r^2>);
  for the H atom in level n this is r_n = (2.5 n^4 + 0.5 n^2)^(1/2) a0, the l = 0 value of
  <r^2> = n^2/2 [5 n^2 + 1 - 3 l(l+1)] a0^2; perturber radii 1.73 a0 (H, = sqrt(3) a0) and 1.02 a0 (He);
- charged perturbers: Hubeny, Hummer & Lanz (1994) Appendix A fit,
  w = f/(1+f), f = 0.1402 x beta^3 / (1 + 0.1285 x beta^(3/2)), x = (1+a)^3.15,
  a = 0.09 n_e^(1/6) T^(-1/2), beta = 8.3e14 n_e^(-2/3) K_n / n^4, K_n = 1 (n <= 3), 16/3 n/(n+1)^2 (n > 3).
"""
import numpy as np
import pytest

import hltest as H
from cases import WCALC_CASES

A0 = 5.29177e-9


def w_neutral(nh, nhe, n):
    r = np.sqrt(2.5 * n ** 4 + 0.5 * n ** 2) * A0
    return np.exp(-4.0 * np.pi / 3.0 * (nh * (r + 1.73 * A0) ** 3 + nhe * (r + 1.02 * A0) ** 3))


def w_charged(ne, n, T):
    if ne <= 10.0 or T <= 10.0:  # the code's guard
        return 1.0
    K = 1.0 if n <= 3 else 16.0 / 3.0 * n / (n + 1.0) ** 2
    a = 0.09 * ne ** (1.0 / 6.0) / np.sqrt(T)
    x = (1.0 + a) ** 3.15
    beta = 8.3e14 * ne ** (-2.0 / 3.0) * K / n ** 4
    f = 0.1402 * x * beta ** 3 / (1.0 + 0.1285 * x * beta ** 1.5)
    return f / (1.0 + f)


# The code rounds the exponents 1/6 and 2/3 to 0.16667 and 0.66667 and uses single-precision
# arguments; at n_e = 1e16 this gives up to ~3e-4 relative difference in w.
RTOL = 1e-3


@pytest.mark.legacy
def test_wcalc_matches_published_formulas_with_neutral_term():
    """Current default: w = w_charged * w_neutral.  A planned change makes the neutral term
    switchable (default off); this test then moves to the 'neutral on' setting."""
    got = H.wcalc(WCALC_CASES)
    ref = np.array([w_charged(ne, n, T) * w_neutral(nh, nhe, n) for nh, ne, nhe, n, T in WCALC_CASES])
    big = ref > 1e-30
    np.testing.assert_allclose(got[big], ref[big], rtol=RTOL)
    assert np.all(got[~big] < 1e-25)


def test_wcalc_charged_part_alone():
    """With no neutral perturbers WCALC is the HHL94 charged-particle factor."""
    cases = [(0.0, ne, 0.0, n, T) for ne in (1e11, 1e13, 1e15, 1e17)
             for n in (2, 3, 4, 6, 10, 20, 40) for T in (4000.0, 7000.0, 15000.0)]
    got = H.wcalc(cases)
    ref = np.array([w_charged(ne, n, T) for _, ne, _, n, T in cases])
    np.testing.assert_allclose(got, ref, rtol=RTOL)


def test_wcalc_limits_and_monotonicity():
    # no perturbers at all: w = 1
    np.testing.assert_allclose(H.wcalc([(0.0, 0.0, 0.0, n, 6000.0) for n in (1, 5, 30)]), 1.0)
    # 0 < w <= 1 and w falls with n at fixed conditions
    for nh, ne, nhe, T in [(1e16, 1e13, 1e15, 5500.0), (1e17, 1e16, 0.0, 12000.0)]:
        w = H.wcalc([(nh, ne, nhe, n, T) for n in range(1, 41)])
        assert np.all((w >= 0) & (w <= 1))
        assert np.all(np.diff(w) <= 0)
