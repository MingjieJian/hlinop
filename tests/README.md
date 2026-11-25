# Tests

The tests run the Fortran routines through small driver programs (`drivers/`), without changes to the sources, and check the results from Python.

Run them with:

```
python -m pip install -r tests/requirements.txt
python -m pytest tests -rxX
```

Requirements:
- `gfortran` (another compiler can be set with `FC=...`);
- Python ≥ 3.9 with numpy and pytest.

The drivers are compiled into a temporary directory under `/tmp` for each run. `hlinprof.f` reads its table file names from `CHARACTER*80` parameters, so the path has to be short. Set `HLINOP_TEST_BUILD=/some/short/dir` to keep the build. A full run takes a few seconds.

GitHub Actions (`.github/workflows/tests.yml`) runs the same command on Linux for every pull request and every push to `develop`/`master`, and on macOS for pushes to `develop`/`master` only (macOS runners can queue for a long time).

## Three kinds of test

1. **Physics tests**: properties that correct code must have, whatever it returns today.
   - profiles have unit area;
   - `WCALC` agrees with an independent transcription of the published formulas (Hummer & Mihalas 1988 eq. 4.71; Hubeny, Hummer & Lanz 1994 App. A) for the default (no neutral term) and for excluded-volume scales 1, 0.5 and 0;
   - `HBOP` recomputes its cached populations when the excluded-volume setting changes;
   - `w` → 1 without perturbers, and `w` falls with n;
   - the continuum is smooth across the Balmer limit;
   - `STARK1` has unit area and matches the Stark-only profiles of Stehlé & Hutcheon (`data/stehle_cds_subset.json`) and of Lemke's VCS tables for Brackett (`data/lemke_vcs_brackett.json`, from `make_vcs_subset.py`); Lyman alpha keeps the earlier profile;
   - `HLINOP` equals `HPROFL` × `HLCUT` (the cutoff weight that HBOP applies to tabulated profiles);
   - `HLINPROF` returns the CDS profiles at table nodes (`data/stehle_cds_subset.json`, from `make_cds_subset.py`), has unit area within the tables' own normalisation, scales quasistatically above a line's last density, and returns HLINOP where there is no table;
   - `HBOP` uses the tabulated profiles by default and cuts them at the series limit.

   A physics test that the current code fails is marked `known_issue` and `xfail` with the reason. `xfail_strict` is on, so once a fix lands the test reports XPASS and fails. The `xfail` mark must then be removed in the same change. The suite therefore records which problems are known and when each one was fixed.

2. **Stable snapshots** (`test_regression.py`, kind `stable`): `HBF`, `HLINPROF` for tabulated lines, and `WCALC` by default and with the Hummer & Mihalas excluded volume. These are change detectors for parts that no planned change should alter. They do not certify the numbers.

3. **Legacy snapshots** (marker `legacy`): `HPROFL`/`HLINOP` Stark profiles and the `HBOP` opacities that use them: `hbop` (defaults), `hbop_hlinop_profiles` (`HBOP_SET_LINE_PROFILE(0)`) and `hbop_legacy_settings` (HLINOP profiles with `HLINOP_SET_CONVOLUTION(0)` and the Hummer & Mihalas excluded volume: the physics settings of the code before Nov 2026, but with the later bug fixes, so no longer equal to it bit for bit). These are parts that planned changes will alter on purpose. A change that alters them regenerates the snapshot in the same commit and says why:

   ```
   python tests/snapshots.py
   ```

   Do not regenerate snapshots to make an unexplained difference go away.

## Known issues

### Covered by `xfail` tests

None at present.

### Fixed (Nov 2026)

- **`HLINPROF` added the Lorentz broadening and ignored microturbulence.** The self-broadening Lorentzian was added beyond 0.2 Å (up to 4% extra area), and the Doppler width passed in was not used for tabulated lines. Both are now convolved with the tabulated profile (`HLPCONV`). The tests check that the area does not change with n_HI, that a Doppler width at or below the thermal one adds nothing, and that a large one sets the core width.
- **Cached convolutions.** They are now stored on 1200 instead of 240 logarithmic points (up to 2% interpolation error in steep shoulders), and the grid size no longer overflows for extremely narrow components (segmentation fault for Ly10 at 1e16 cm⁻³).
- **`NBF = 6`.** Bound–free absorption from n ≥ 7 was missing, although `HBF` has Gaunt factors up to n = 15 and the header said 15. It is now 15. `HBOP_AUTO` then works up to the n = 15 edge (20.5 µm); before, it stopped at 3.28 µm.
- **`HBOP` returned Inf/NaN beyond about 5 µm.** Far redward of an edge the extrapolated `HBF` overflowed while the taper had underflowed to 0. Such levels are now skipped.
- **The dissolution edge in `HBOP` was in the wrong place.** n\* was computed with `IONH` (the infinite-mass Rydberg) instead of the code's own level energies (R_H), which put the n = 2 dissolution edge 2 Å blueward of the series limit. The air wavelength was also used as if it were in vacuum. `HBNSTAR` now uses the code's level energies, and the continuum is computed in vacuum wavelengths (`WAVE` is in air above 2000 Å, as for the lines). The edge now coincides with the series limit at which `HLINOP` cuts the lines.
- **Self-broadening density.** `HBOP` passed the ground-state population `NP(1)` as the H I density for broadening by hydrogen; it now passes `NH`.
- **`HPROFL` added the broadening mechanisms outside the core.** Inside the core it used the dominant mechanism; outside it, it added all of them. Where Doppler broadening dominated, this added 2–39% to the area and left Stark broadening out of the core. The mechanisms are now convolved (`HPCONV`). Areas are 1 to 0.1%. The end-to-end test compares `HPROFL` with the Doppler-convolved Stehlé & Hutcheon profiles (within 10%).
- **`STARK1` added the ion quasistatic, electron quasistatic and electron impact parts instead of convolving them.** Areas were 1.74–2.14. The parts are now convolved (B: norm(ion) ⊛ norm(impact + electron wing)), with `HLINOP_SET_CONVOLUTION(0)` keeping the old code. The tests compare with Stehlé & Hutcheon (n_up ≥ 8; HWHM and peak within 10%) and with Lemke VCS for Brackett (within 15%).

### Found while writing the tests, not yet covered

3. **`STARK1` cache key.** The old `STARK1` caches `C1CON`/`C2CON`, which depend on `WAVEH`, keyed on (N, M) only, so a call with the same (N, M) and a different line centre reuses stale values. It is still so in `STARK1L`, kept unchanged. `STARK1PARTS`, used by the convolved profile, keys on (N, M, `WAVEH`).
4. **Undefined variable in `STARK1L`.** The old `STARK1`, now `STARK1L`, computes `T3NHE` from `HE1FRC`, which is not one of its arguments; the value is never used. It is removed in `STARK1PARTS`.
7. **`STEHLE` interpolation narrows Doppler-dominated cores and loses area.** STEHLE interpolates the tabulated profiles linearly in log S against log α, with nodes 0.2 dex apart. In Doppler-dominated cores (ln S ∝ −α²) the chord lies below the profile. The result:
   - half widths are 1–6% narrower than the widths given in the CDS headers (Hα, Hβ, H10, Pα at 5000 K, 1e10–1e14 cm⁻³);
   - integrated areas come out at 0.97–0.99.

   A monotone cubic (PCHIP) in the same log–log variables gives areas of 0.987–1.000 and widths of 0.98–0.99 of the header values. So the earlier conclusion that the published tables are not normalised was mostly an artefact of this interpolation. The same interpolation is in SMElib and in Turbospectrum's HLINPROF. Not yet changed. The `HLINPROF` area tolerance is 5%.
   - Lyman alpha is left to HLINOP on purpose (its far-wing satellites); to be revisited.
8. **`hlinprof.f` does not compile on its own.** It needs a `DATA.FILES` include that only `SYNTH/Makefile` generates. The tests generate it themselves.
