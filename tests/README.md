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
   - `STARK1` has unit area and matches the Stark-only profiles of Stehlé & Hutcheon (`data/stehle_cds_subset.json`) and of Lemke's VCS tables for Brackett (`data/lemke_vcs_brackett.json`, from `make_vcs_subset.py`); Lyman alpha keeps the earlier profile in HLINOP (HLINPROF has the tabulated one);
   - `HLINOP` equals `HPROFL` × `HLCUT` (the cutoff weight that HBOP applies to tabulated profiles);
   - `HLINPROF` returns the CDS profiles at table nodes (`data/stehle_cds_subset.json`, from `make_cds_subset.py`), has unit area within the tables' own normalisation, scales quasistatically above a line's last density, and returns HLINOP where there is no table;
   - `HBOP` uses the tabulated profiles by default and cuts them at the series limit;
   - the cached convolutions give the same numbers whatever the call order (grouped by line and depth, or wavelength by wavelength as a synthesis code calls them), and the cache index `HLCKEY` stores, finds and empties correctly, also on `HLINOP_CLEAR_CACHE` (`test_cache.py`);
   - interfaces for synthesis codes (`test_interface.py`): `HLINPROF_SET_FILES` uses the given table, `HLINPROF_SOURCE` reports the fallback to HLINOP, and `HBOP_BF` with `HBOP_LEVELS` equals the continuum of `HBOP` bit for bit.

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

- **The convolution caches thrashed in the call order of synthesis codes.** `STARK1C`, `HPCONV` and `HLPCONV` kept 128 states in turn. A code such as SME calls the profiles wavelength by wavelength, and at each wavelength for every depth and line, so with more than 128 (line, depth) states every call rebuilt the convolution (about 4.5 ms, against 0.1 µs for the earlier code). The caches now hold up to 6144 states each (hash index `HLCKEY` over densely stored entries, emptied when full), the cached profiles are stored as single-precision logarithms (changes below 2e-6; 9.6 kB per state, memory taken only for the states used; `HLINOP_CLEAR_CACHE` empties them, for a synthesis code to call at each new model), and builds are faster (real-input FFT with a twiddle table; `STEHLE` with precomputed logarithms and slopes, bisection, and both sides of the line in one call, bit for bit as before). For 72 depths × 28 Balmer lines called wavelength by wavelength: 0.7 ms per state to build (once), then 0.08 µs per call (0.05–0.1 µs with `HLINOP_SET_CONVOLUTION(0)`); before, 4.4 ms per call. The cached profiles are interpolated from the FFT grid in log space instead of linearly (closer to a 64 times finer reference by a factor of about 2 in the median); snapshots changed by up to 1.6e-4.
- **Lyman alpha was left to HLINOP.** The Griem-type Stark profile is 1.7-13 times too wide for the Lyman alpha core (no unshifted components) and the earlier HLINOP profile has half widths 0.1-0.9 of the Stehlé & Hutcheon ones and peaks up to 9 times too high at n_e = 1e15-1e18 cm⁻³. HLINPROF now takes the core from the table (convolved with the natural, resonance and He widths), multiplies the red side by the Boltzmann factor, and joins it linearly at -4000 cm⁻¹ to the Allard et al. H–H and H–H⁺ far wing (as in HLINOP) plus half the tabulated quasistatic wing for the electrons. Tests: blue side against the CDS profiles plus the natural Lorentzian (1.5%), half width (3%), continuity at -4000 cm⁻¹, far wing equal to HLINOP at 1300–3000 Å (1%).
- **`STEHLE` interpolation narrowed Doppler-dominated cores and lost area.** It interpolated linearly in log S against log α between nodes 0.2 dex apart. The chord lies below the profile in Doppler-dominated cores, so half widths were 1–6% narrower than the CDS header widths and areas came out at 0.96–0.99. This, not the published tables, was the cause of the earlier "tables not normalised" note. It now uses a monotone cubic (Fritsch–Carlson/PCHIP) in the same variables: areas are 0.98–1.00 and half widths are within 3% of the header widths. The tests also check monotonicity between nodes. `HLINOP_SET_CONVOLUTION(0)` keeps the linear form. The self-broadening grids (`HGRID`, `ISELF = 1`, not the default) still interpolate linearly.
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
8. **`hlinprof.f` does not compile on its own.** It needs a `DATA.FILES` include that only `SYNTH/Makefile` generates (the names can be overridden at run time with `HLINPROF_SET_FILES`, but the include is still needed to compile). The tests generate it themselves.
