"""Wavelength medium (HLINOP_SET_MEDIUM): air above 2000 A (default) or vacuum.  The medium sets
where series limits, red cuts and bound-free edges fall; profiles are relative to the line centre."""
import numpy as np

import hltest as H

LIMIT2_VAC = 1e8 / (H.IONH_CM - H.level_energy(2))       # Balmer series limit in vacuum
H10_VAC = H.line_centre(2, 10)
WAVES_VAC = np.array([3640.0, LIMIT2_VAC + 0.3, 3700.0, H10_VAC - 2.0, H10_VAC + 0.5, 4500.0, 8300.0,
                      16000.0])


def _calls(waves, centre):
    return [(w, 2, 10, centre) for w in waves]


def test_default_is_air():
    calls = _calls(H.vacair(WAVES_VAC), H.vacair(H10_VAC))
    np.testing.assert_array_equal(H.medium(None, calls), H.medium(1, calls))


def test_vacuum_medium_matches_air_at_the_same_physical_wavelengths():
    """Vacuum wavelengths with HLINOP_SET_MEDIUM(0) give what air wavelengths give in air: the same
    cuts, the same bound-free opacity (to the air/vacuum round trip) and the same profiles (to the
    frequency scale of air wavelengths, 3e-4)."""
    vac = H.medium(0, _calls(WAVES_VAC, H10_VAC))
    air = H.medium(1, _calls(H.vacair(WAVES_VAC), H.vacair(H10_VAC)))
    np.testing.assert_array_equal(vac[:, 0], air[:, 0])
    np.testing.assert_allclose(vac[:, 2], air[:, 2], rtol=1e-5)
    np.testing.assert_allclose(vac[:, 1], air[:, 1], rtol=1e-3, atol=1e-30)


def test_series_limit_in_vacuum():
    """In vacuum the Balmer lines are cut at the vacuum series limit (3647.05 A); in air at the air
    limit (3645.96 A), so a vacuum wavelength between the two is cut only in vacuum."""
    w = 0.5 * (LIMIT2_VAC + H.vacair(LIMIT2_VAC))
    c = _calls([w], H10_VAC)
    assert H.medium(0, c)[0, 0] == 0.0 and H.medium(1, c)[0, 0] == 1.0
