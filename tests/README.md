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
   - `HLINOP` equals `HPROFL` × `HLCUT` (the cutoff weight that HBOP applies to tabulated profiles);
   - `HLINPROF` returns the CDS profiles at table nodes (`data/stehle_cds_subset.json`, from `make_cds_subset.py`), has unit area within the tables' own normalisation, scales quasistatically above a line's last density, and returns HLINOP where there is no table;
   - `HBOP` uses the tabulated profiles by default and cuts them at the series limit.

   A physics test that the current code fails is marked `known_issue` and `xfail` with the reason. `xfail_strict` is on, so once a fix lands the test reports XPASS and fails. The `xfail` mark must then be removed in the same change. The suite therefore records which problems are known and when each one was fixed.

2. **Stable snapshots** (`test_regression.py`, kind `stable`): `HBF`, `HLINPROF` for tabulated lines, and `WCALC` by default and with the Hummer & Mihalas excluded volume. These are change detectors for parts that no planned change should alter. They do not certify the numbers.

3. **Legacy snapshots** (marker `legacy`): `HPROFL`/`HLINOP` Stark profiles and the `HBOP` opacities that use them: `hbop` (defaults), `hbop_hlinop_profiles` (`HBOP_SET_LINE_PROFILE(0)`) and `hbop_original` (HLINOP profiles and the Hummer & Mihalas excluded volume, i.e. the code before Nov 2026; equal bit for bit to the snapshots taken before the changes). These are parts that planned changes will alter on purpose. A change that alters them regenerates the snapshot in the same commit and says why:

   ```
   python tests/snapshots.py
   ```

   Do not regenerate snapshots to make an unexplained difference go away.

## Known issues

### Covered by `xfail` tests

- **`STARK1` adds instead of convolving.** It adds the ion quasistatic, electron quasistatic and electron impact parts instead of convolving them. `HPROFL`/`HLINOP` areas are 1.74–2.14 wherever Stark broadening dominates (for example Balmer, Paschen and Brackett lines at T = 5000–10000 K, n_e = 1e12–1e16 cm⁻³). They are 1.02–1.13 only where Doppler broadening dominates the core.

### Found while writing the tests, not yet covered

1. **The dissolution edge in `HBOP` is in the wrong place.** The effective quantum number n\* in `HBOP` uses `IONH = 2.17991e-11` erg, the infinite-mass Rydberg (109 737 cm⁻¹). The level energies use R_H = 109 678.764 cm⁻¹. So the n = 2 dissolution edge sits at 3644.98 Å, while the series limit from the code's own levels is at 3647.02 Å (both vacuum). In addition, `HLINOP` cuts non-Lyman lines at the series limit in air, while the `HBOP` continuum is computed in vacuum.
2. **`NBF` disagrees with the header.** `NBF = 6`, so bound–free absorption from n ≥ 7 is omitted, while the header says the limit is 15.
3. **`STARK1` cache key.** `STARK1` caches `C1CON`/`C2CON`, which depend on `WAVEH`, keyed on (N, M) only. A call with the same (N, M) and a different line centre reuses stale values. SMElib added `WAVEH` to the key.
4. **Undefined variable in `STARK1`.** `STARK1` computes `T3NHE` from `HE1FRC`, which is not one of its arguments; the value is never used.
5. **Self-broadening density.** `HBOP` passes `NP(1)`, the ground-state population, as the H I density used for self-broadening, instead of the total `NH` (equal to < 0.1% in FGK photospheres).
6. **`HLINPROF` adds area.** The self-broadening Lorentzian is added (not convolved) beyond 0.2 Å from the centre. This adds up to 4% to the area (Hγ, n_HI = 1e17). The test tolerance is 5%.
7. **The published Stehle tables are not quite normalised.** Integrating the CDS profiles directly over their tabulated range gives 0.97–0.99 (Hβ, H10, Pa14, Lyβ at 1e14 cm⁻³), and HLINPROF reproduces this. The `HLINPROF` area tolerance is 5%.
   - Tabulated profiles contain thermal Doppler broadening only; the Doppler width passed in (with microturbulence) is used only for lines that go to HLINOP.
   - Lyman alpha is left to HLINOP on purpose (its far-wing satellites); to be revisited.
8. **`hlinprof.f` does not compile on its own.** It needs a `DATA.FILES` include that only `SYNTH/Makefile` generates. The tests generate it themselves.
