"""Helpers for the HLINOP/HBOP test suite: build the Fortran drivers and call them.

The drivers in tests/drivers/ are small programs that read cases from stdin and write
results to stdout, so the Fortran sources are tested exactly as they are, without changes.
"""
import atexit
import os
import shutil
import subprocess
import tempfile
from functools import lru_cache
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DRIVERS = ROOT / "tests" / "drivers"
FC = os.environ.get("FC", "gfortran")
FFLAGS = ["-O2", "-std=legacy", "-w"]
C_ANGSTROM = 2.99792458e18  # speed of light in A/s
RYDBERG_H = 109677.576      # cm^-1, as used in hlinop.f/hbop.f for n >= 9
IONH_CM = 109678.764        # H ionisation energy in cm^-1, as in hlinop.f/hbop.f

# which Fortran sources each driver needs
PROGRAMS = {
    "drv_wcalc": ["hbop.f", "hlinprof.f", "hlinop.f"],
    "drv_profile": ["hlinop.f"],
    "drv_hbop": ["hbop.f", "hlinprof.f", "hlinop.f"],
    "drv_hbf": ["hbop.f", "hlinprof.f", "hlinop.f"],
    "drv_hlinprof": ["hlinprof.f", "hlinop.f"],
    "drv_nstar": ["hbop.f", "hlinprof.f", "hlinop.f"],
    "drv_hbop_auto": ["hbop_auto.f", "hbop.f", "hlinprof.f", "hlinop.f"],
    "drv_states": ["hlinprof.f", "hlinop.f"],
    "drv_hlckey": ["hlinop.f"],
    "drv_files": ["hlinprof.f", "hlinop.f"],
    "drv_hbopbf": ["hbop.f", "hlinprof.f", "hlinop.f"],
}
DATA = {"HFILE": "stehle_long.dat", "HSELFFILE": "bpo_self.grid.DEC", "HVCSFILE": "vcsbalmer.dat"}


@lru_cache(maxsize=None)
def build_dir():
    """Directory with the compiled drivers.  hlinprof.f reads its table file names from an
    INCLUDE file of fixed-form PARAMETER lines, so the data files are linked into a short path."""
    if os.environ.get("HLINOP_TEST_BUILD"):
        d = Path(os.environ["HLINOP_TEST_BUILD"])
    else:  # short path under /tmp, removed at exit
        d = Path(tempfile.mkdtemp(prefix="hlt", dir="/tmp"))
        atexit.register(shutil.rmtree, d, True)
    d.mkdir(parents=True, exist_ok=True)
    lines = []
    for key, name in DATA.items():
        dst = d / name
        if not dst.exists():
            os.symlink(ROOT / name, dst)
        if len(str(dst)) > 120:
            raise RuntimeError(f"data path too long for a fixed-form PARAMETER line: {dst}")
        lines += [f"      PARAMETER ({key}=", f"     ;'{dst}')"]
    (d / "DATA.FILES").write_text("\n".join(lines) + "\n")
    for prog, srcs in PROGRAMS.items():
        cmd = [FC, *FFLAGS, "-ffixed-line-length-132", f"-I{d}", str(DRIVERS / f"{prog}.f"),
               *[str(ROOT / s) for s in srcs], "-o", str(d / prog)]
        subprocess.run(cmd, check=True, cwd=d)
    return d


def run(prog, text):
    """Run a driver with the given stdin text; return the output as a 2-D float array."""
    out = subprocess.run([str(build_dir() / prog)], input=text, capture_output=True,
                         text=True, check=True).stdout
    rows = []
    for ln in out.splitlines():  # the codes print informational messages; keep numeric rows only
        try:
            rows.append([float(x) for x in ln.split()])
        except ValueError:
            continue
    rows = [r for r in rows if r]
    return np.array(rows, dtype=float)


def wcalc(cases, scale=None):
    """cases: iterable of (NH, NE, NHE, NS, T).  scale: excluded-volume scale passed to
    HBOP_SET_EXCLUDED_VOLUME before each case; None = never call it (the code's default)."""
    sc = -1.0 if scale is None else float(scale)
    txt = "".join("%.8e %.8e %.8e %.8e %.8e %.8e\n" % (*c, sc) for c in cases)
    return run("drv_wcalc", txt)[:, 0]


def hbf(cases):
    """cases: iterable of (wave_A, nlo)."""
    txt = "".join("%.10f %d\n" % (w, n) for w, n in cases)
    return run("drv_hbf", txt)[:, 2]


def level_energy(n):
    """Level energies (cm^-1) exactly as tabulated in hlinop.f/hbop.f."""
    tab = {1: 0.0, 2: 82259.105, 3: 97492.302, 4: 102823.893, 5: 105291.651,
           6: 106632.160, 7: 107440.444, 8: 107965.051}
    return tab[n] if n in tab else IONH_CM - RYDBERG_H / n ** 2


def line_centre(nlo, nup):
    """Vacuum line centre (A) from the code's level energies."""
    return 1e8 / (level_energy(nup) - level_energy(nlo))


def doppler(T, vturb_kms=1.5):
    """Reduced Doppler width (delta_lambda/lambda) for H at T with a microturbulence."""
    return np.sqrt(2 * 1.380649e-16 * T / 1.6735575e-24 + (vturb_kms * 1e5) ** 2) / 2.99792458e10


def detuning_grid(wave0, rel_max=0.15, dmin=1e-4, npts=3000):
    """Wavelengths symmetric about wave0: 0, then log-spaced |delta lambda| from dmin to rel_max*wave0."""
    d = np.geomspace(dmin, rel_max * wave0, npts)
    return np.concatenate([wave0 - d[::-1], [wave0], wave0 + d])


def profile(nlo, nup, T, ne, nh, nhe, dop, waves, wave0=None, with_cut=False, convolution=None, stark1=False):
    """HPROFL and HLINOP (per Hz) at the given wavelengths (and the HLCUT weight, or STARK1).
    convolution: mode for HLINOP_SET_CONVOLUTION (None = default)."""
    wave0 = line_centre(nlo, nup) if wave0 is None else wave0
    head = "%d %d %.10f %.6e %.6e %.6e %.6e %.6e %d\n" % (nlo, nup, wave0, T, ne, nh, nhe, dop,
                                                          -1 if convolution is None else convolution)
    out = run("drv_profile", head + "".join("%.10f\n" % w for w in waves))
    if stark1:
        return out[:, 4]
    if with_cut:
        return out[:, 1], out[:, 2], out[:, 3]
    return out[:, 1], out[:, 2]


def hlinprof(nlo, nup, T, ne, nh, nhe, dop, waves, wave0=None, convolution=None):
    """HLINPROF (per A) at the given wavelengths; convolution: HLINOP_SET_CONVOLUTION mode."""
    wave0 = line_centre(nlo, nup) if wave0 is None else wave0
    head = "%d %d %.10f %.6e %.6e %.6e %.6e %.6e %d\n" % (nlo, nup, wave0, T, ne, nh, nhe, dop,
                                                          -1 if convolution is None else convolution)
    return run("drv_hlinprof", head + "".join("%.10f\n" % w for w in waves))[:, 1]


def area_per_frequency(waves, phi_nu, tail=True):
    """Integral of a per-Hz profile over frequency (trapezoid), plus an estimate of the far wings
    beyond the grid assuming the quasistatic |delta nu|^-5/2 fall-off."""
    nu = C_ANGSTROM / waves
    o = np.argsort(nu)
    nu, p = nu[o], phi_nu[o]
    a = np.trapezoid(p, nu) if hasattr(np, "trapezoid") else np.trapz(p, nu)
    if tail:
        nu0 = nu[len(nu) // 2]
        for i in (0, -1):
            a += 2.0 / 3.0 * p[i] * abs(nu[i] - nu0)
    return float(a)


def area_per_angstrom(waves, phi_lam):
    o = np.argsort(waves)
    w, p = waves[o], phi_lam[o]
    return float(np.trapezoid(p, w) if hasattr(np, "trapezoid") else np.trapz(p, w))


def hbop(nh, nhe, ne, T, dop, lines, waves, scale=None, profile=None, convolution=None):
    """HBOP total and continuous opacity (cm^-1); lines: iterable of (nlo, nup, waveh).
    scale: excluded-volume scale (HBOP_SET_EXCLUDED_VOLUME); profile: line-profile mode
    (HBOP_SET_LINE_PROFILE).  Each is None (never set: the code's default), a number, or one
    number per wavelength (all calls run in one process, in order)."""
    lines = list(lines)
    waves = list(waves)

    def per_wave(x, default):
        if x is None or np.isscalar(x):
            return [default if x is None else x] * len(waves)
        return list(x)
    sc = per_wave(scale, -1.0)
    md = per_wave(profile, -1)
    sk = per_wave(convolution, -1)
    head = "%.6e %.6e %.6e %.6e %.6e %d\n" % (nh, nhe, ne, T, dop, len(lines))
    head += "".join("%d %d %.10f\n" % l for l in lines)
    out = run("drv_hbop", head + "".join("%.10f %.8e %d %d\n" % (w, float(s), int(m), int(k))
                                         for w, s, m, k in zip(waves, sc, md, sk)))
    return out[:, 1], out[:, 2]


def nstar(cases):
    """cases: iterable of (wave, nlo); returns array with columns wave, AIRVAC, VACAIR,
    HBNSTAR(wave), HBNSTAR(AIRVAC(wave))."""
    return run("drv_nstar", "".join("%.10f %d\n" % c for c in cases))


def hbop_auto(nmax, nh, nhe, ne, T, dop, waves):
    """Run HBOP_AUTO; returns (total, contin, stderr)."""
    head = "%d %.6e %.6e %.6e %.6e %.6e\n" % (nmax, nh, nhe, ne, T, dop)
    r = subprocess.run([str(build_dir() / "drv_hbop_auto")], input=head + "".join("%.10f\n" % w for w in waves),
                       capture_output=True, text=True)
    rows = []
    for ln in r.stdout.splitlines():
        try:
            rows.append([float(x) for x in ln.split()])
        except ValueError:
            continue
    a = np.array([x for x in rows if x]).reshape(-1, 3)
    return a[:, 1], a[:, 2], r.stderr + r.stdout


def states(calls):
    """calls: iterable of (nlo, nup, wave0, T, ne, nh, nhe, dop, wave), evaluated in this order in
    one process.  Returns HLINPROF (per A) and HLINOP (per Hz)."""
    txt = "".join("%d %d %.10f %.6e %.6e %.6e %.6e %.6e %.10f\n" % c for c in calls)
    out = run("drv_states", txt)
    return out[:, 0], out[:, 1]


def hlckey(keys):
    """The cache index HLCKEY (16 hash slots, 12 entries): (entry, hit, number used) for each key."""
    return run("drv_hlckey", "".join("%.17g\n" % k for k in keys)).astype(int)


def files(stehle_path, calls, messages=None, raw=False):
    """HLINPROF with the Stehle table set by HLINPROF_SET_FILES (None: the DATA.FILES default) and
    the message mode set by HLINPROF_SET_MESSAGES (None: not called).  calls: (nlo, nup, wave0, T,
    ne, nh, nhe, dop, wave).  Returns PROF (per A) and HLINPROF_SOURCE, or the raw output."""
    head = ("-" if stehle_path is None else str(stehle_path)) + "\n"
    head += "%d\n" % (-1 if messages is None else messages)
    txt = head + "".join("%d %d %.10f %.6e %.6e %.6e %.6e %.6e %.10f\n" % c for c in calls)
    if raw:
        return subprocess.run([str(build_dir() / "drv_files")], input=txt, capture_output=True, text=True)
    out = run("drv_files", txt)
    return out[:, 0], out[:, 1].astype(int)


def hbop_bf(nh, nhe, ne, T, waves):
    """HBOP_LEVELS, HBOP_BF and HBOP_BF_DISSOLVED against HBOP's continuum: (sum NLTE / NH, W(1),
    W(15)) and an array of (CONTIN from HBOP, HBOP_BF, HBOP_BF_DISSOLVED to n = 15, to n = 4) per
    wavelength."""
    txt = "%.8e %.8e %.8e %.8e\n" % (nh, nhe, ne, T) + "".join("%.10f\n" % w for w in waves)
    out = run("drv_hbopbf", txt)
    return out[0][:3], out[1:, 1:]

