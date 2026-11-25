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


def test_changing_the_excluded_volume_recomputes_the_populations():
    """HBOP caches level populations for unchanged (NH, NHE, NE, T).  Changing the excluded-volume
    setting between calls at the same conditions must give the same result as a fresh run."""
    _, nh, nhe, ne, T = HBOP_CONDITIONS[0]
    lines = [(2, m, H.line_centre(2, m)) for m in range(3, 31)]
    w = 3700.0
    seq = [1.0, 0.0, 1.0, 0.5, 0.0]
    tot, cont = H.hbop(nh, nhe, ne, T, H.doppler(T), lines, [w] * len(seq), scale=seq)
    for i, s in enumerate(seq):
        t1, c1 = H.hbop(nh, nhe, ne, T, H.doppler(T), lines, [w], scale=s)
        np.testing.assert_allclose([tot[i], cont[i]], [t1[0], c1[0]], rtol=1e-6)
    assert abs(tot[0] / tot[1] - 1.0) > 0.05  # the setting matters here (17% at 3700 A)


def _line_opacity(nlo, nup, waves, profile, cond=HBOP_CONDITIONS[0]):
    _, nh, nhe, ne, T = cond
    tot, cont = H.hbop(nh, nhe, ne, T, H.doppler(T), [(nlo, nup, H.line_centre(nlo, nup))], waves,
                       profile=profile)
    return tot - cont


@pytest.mark.parametrize("line", [(2, 4), (2, 10), (2, 22), (3, 14), (1, 3)], ids=lambda l: "%d-%d" % l)
def test_hbop_uses_the_tabulated_profiles(line):
    """By default the line opacity of a tabulated line follows HLINPROF (per A, converted to per
    Hz); HBOP_SET_LINE_PROFILE(0) follows HLINOP.  Their ratio equals the profile ratio."""
    nlo, nup = line
    _, nh, nhe, ne, T = HBOP_CONDITIONS[0]
    w0 = H.line_centre(nlo, nup)
    w = w0 + np.array([0.0, 0.5, 2.0, -2.0])
    l1 = _line_opacity(nlo, nup, w, None)
    l0 = _line_opacity(nlo, nup, w, 0)
    # HBOP passes the total H I density NH for the broadening by hydrogen (self broadening of
    # H alpha - H gamma in HLINPROF, resonance broadening in HLINOP)
    prof_tab = H.hlinprof(nlo, nup, T, ne, nh, nhe, H.doppler(T), w, wave0=w0) * w ** 2 / H.C_ANGSTROM
    _, prof_hlinop = H.profile(nlo, nup, T, ne, nh, nhe, H.doppler(T), w, wave0=w0)
    np.testing.assert_allclose(l1 / l0, prof_tab / prof_hlinop, rtol=1e-4)


def test_tabulated_lines_are_cut_at_the_series_limit():
    """Like HLINOP, the tabulated profiles in HBOP stop at the series limit: no line opacity of
    a Balmer line blueward of 3646 A (air)."""
    w = np.array([3600.0, 3620.0, 3640.0])
    np.testing.assert_array_equal(_line_opacity(2, 10, w, None), 0.0)
    assert np.all(_line_opacity(2, 10, np.array([3700.0, 3790.0]), None) > 0)


def test_dissolution_edge_is_the_series_limit():
    """HBNSTAR gives n* from the code's level energies, 1/n*^2 = (CONTH(n) - sigma_vac)/R_H, and 0
    blueward of the series limit; with HBOP's air-to-vacuum conversion the edge coincides with
    the blue cut that HLINOP/HLCUT apply to the lines (series limit in air)."""
    for nlo in (1, 2, 3, 4, 6, 10, 15):
        limit_vac = 1e8 / (H.IONH_CM - H.level_energy(nlo))
        waves = limit_vac * np.array([0.98, 0.99999, 1.00001, 1.01, 1.2, 1.5])
        out = H.nstar([(w, nlo) for w in waves])
        sig = 1e8 / waves
        q = H.IONH_CM - H.level_energy(nlo) - sig
        ref = np.where(q > 0, np.sqrt(H.RYDBERG_H / np.maximum(q, 1e-30)), 0.0)
        np.testing.assert_allclose(out[:, 3], ref, rtol=1e-6)
    # air: the edge is at VACAIR(limit) for the levels above 1 (as the line cut)
    for nlo in (2, 3, 4):
        limit_vac = 1e8 / (H.IONH_CM - H.level_energy(nlo))
        cut_air = H.nstar([(limit_vac, nlo)])[0, 2]
        out = H.nstar([(cut_air - 0.01, nlo), (cut_air + 0.01, nlo)])
        assert out[0, 4] == 0.0 and out[1, 4] > 100.0


def test_hbop_is_finite_in_the_infrared():
    """Far redward of an edge the extrapolated HBF overflows while the red-wing taper underflows;
    such levels are skipped (TOTAL and CONTIN were Inf at 6 micron and NaN at 8 micron)."""
    _, nh, nhe, ne, T = HBOP_CONDITIONS[0]
    lines = [(2, m, H.line_centre(2, m)) for m in range(3, 31)]
    tot, cont = H.hbop(nh, nhe, ne, T, H.doppler(T), lines, [40000.0, 60000.0, 80000.0, 120000.0])
    assert np.all(np.isfinite(tot)) and np.all(np.isfinite(cont)) and np.all(cont > 0)


def test_hbop_auto_reaches_20_micron():
    """With NBF = 15 HBOP_AUTO works up to the n = 15 bound-free edge (20.5 micron; it stopped
    at the n = 6 edge, 3.28 micron, before) and stops beyond it with a message."""
    total, cont, _ = H.hbop_auto(60, 1e16, 1e15, 1e13, 6000.0, H.doppler(6000.0), [40000.0, 150000.0, 200000.0])
    assert len(total) == 3 and np.all(cont > 0) and np.all(total >= cont)
    total, cont, msg = H.hbop_auto(60, 1e16, 1e15, 1e13, 6000.0, H.doppler(6000.0), [210000.0])
    assert len(total) == 0 and "20.5 micron" in msg
