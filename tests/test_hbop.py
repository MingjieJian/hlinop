"""HBOP (hbop.f): continuous (bound-free plus dissolved-level) and total opacity."""
import numpy as np
import pytest

import hltest as H
from cases import BALMER_LIMIT_WAVES, HBOP_CONDITIONS


@pytest.mark.parametrize("cond", HBOP_CONDITIONS, ids=[c[0] for c in HBOP_CONDITIONS])
def test_continuum_is_smooth_across_the_balmer_limit(cond):
    """The dissolved fraction D = 1 - w(n*)/w(n) goes to 1 at the series limit, so the
    pseudo-continuum joins the bound-free continuum without a jump: no step larger than 1%
    between 0.05-A neighbours from 3600 to 3720 A."""
    _, nh, nhe, ne, T = cond
    _, cont = H.hbop(nh, nhe, ne, T, H.doppler(T), [], BALMER_LIMIT_WAVES)
    assert np.all(cont > 0)
    step = np.abs(np.diff(np.log(cont)))
    assert step.max() < 0.01, BALMER_LIMIT_WAVES[np.argmax(step)]


@pytest.mark.parametrize("cond", HBOP_CONDITIONS, ids=[c[0] for c in HBOP_CONDITIONS])
def test_lines_add_to_the_continuum(cond):
    """TOTAL = CONTIN + line opacity, and the line opacity is positive in LTE (no lasing)."""
    _, nh, nhe, ne, T = cond
    lines = [(2, m, H.line_centre(2, m)) for m in range(3, 31)]
    waves = np.array([H.line_centre(2, m) for m in (4, 10, 20, 30)])
    tot, cont = H.hbop(nh, nhe, ne, T, H.doppler(T), lines, waves)
    _, cont0 = H.hbop(nh, nhe, ne, T, H.doppler(T), [], waves)
    np.testing.assert_allclose(cont, cont0, rtol=1e-6)
    assert np.all(tot > cont)
