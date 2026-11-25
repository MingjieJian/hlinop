"""Shared test conditions.  Representative of FGK photospheres (T 5000-7000 K, n_e 1e12-1e15)
and of white-dwarf atmospheres (T ~1e4 K, n_e up to 1e16-1e17)."""
import numpy as np

# (NH, NE, NHE, n, T) for the occupation probability
WCALC_CASES = [(nh, ne, nhe, n, T)
               for nh in (1e14, 1e16, 1e17)
               for ne in (1e12, 1e14, 1e16)
               for nhe in (0.0, 1e16)
               for n in (1, 2, 3, 4, 5, 8, 12, 20, 30)
               for T in (5000.0, 10000.0)]

# lines for the profile tests: Balmer, Paschen and Brackett, low to high n
PROFILE_LINES = [(2, 5), (2, 8), (2, 10), (2, 14), (2, 22), (2, 30),
                 (3, 10), (3, 14), (3, 16), (3, 25),
                 (4, 11), (4, 15), (4, 20)]
# (T, n_e, n_HI, n_HeI); n_e 1e16 is used only for n_up <= 10, where +-15% of the line
# wavelength still contains the profile (far wings < 1%)
PROFILE_CONDITIONS = [(5000.0, 1e12, 1e16, 1e15), (6300.0, 1e14, 1e17, 1e16), (10000.0, 1e16, 1e16, 1e15)]


def profile_cases():
    for nlo, nup in PROFILE_LINES:
        for T, ne, nh, nhe in PROFILE_CONDITIONS:
            if ne > 1e15 and nup > 10:
                continue
            yield nlo, nup, T, ne, nh, nhe


# HLINPROF (Stehle tables in this repository cover H alpha - H gamma)
HLINPROF_CASES = [(2, nup, T, ne, nh) for nup in (3, 4, 5)
                  for T, ne, nh in ((5000.0, 1e12, 1e16), (6300.0, 1e14, 1e17), (10000.0, 1e15, 1e16))]

# HLINPROF outside its Stehle grid (T 2500-10000 K, log n_e 10-16): falls back to VCS
# (vcsbalmer.dat).  (T, n_e, n_HI)
HLINPROF_VCS_CONDITIONS = [(12000.0, 1e15, 1e16), (20000.0, 1e16, 1e15), (12000.0, 1e17, 1e16),
                           (5000.0, 1e17, 1e17), (5000.0, 1e9, 1e15)]

# HBOP: (label, NH, NHE, NE, T)
HBOP_CONDITIONS = [("FGK dwarf", 1e17, 1e16, 1e14, 6300.0),
                   ("FGK giant", 1e16, 1e15, 1e12, 4800.0),
                   ("DA white dwarf", 1e17, 0.0, 1e16, 12000.0)]
# wavelengths (A) around the Balmer limit
BALMER_LIMIT_WAVES = np.round(np.arange(3600.0, 3720.0001, 0.05), 3)
