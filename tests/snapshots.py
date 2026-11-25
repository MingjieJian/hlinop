"""Regression snapshots: what the code returns now for a fixed set of inputs.

A snapshot is a change detector, not a statement that the numbers are right.  Two kinds:
- stable: parts that no planned change should alter (HBF, HLINPROF for tabulated lines,
  WCALC by default and with the Hummer & Mihalas excluded volume);
- legacy: parts that planned changes will alter on purpose (HPROFL/HLINOP Stark profiles and
  the HBOP opacities that use them: by default for the lines without tables, with
  HBOP_SET_LINE_PROFILE(0) for all lines).  A change that alters them must
  regenerate the snapshot (python tests/snapshots.py) in the same commit and say why.
"""
import json
from pathlib import Path

import numpy as np

import hltest as H
from cases import HBOP_CONDITIONS, WCALC_CASES

FILE = Path(__file__).resolve().parent / "data" / "snapshots.json"

HBF_INPUTS = [(float(w), n) for n in range(1, 7)
              for w in np.geomspace(300.0, 1.2 * 911.27 * n * n, 12).round(3)]
PROFILE_SETS = [(2, 10, 6300.0, 1e14, 1e17, 1e16), (3, 14, 5000.0, 1e13, 1e16, 1e15),
                (4, 15, 4800.0, 1e12, 1e16, 1e15)]
PROFILE_OFFSETS = [0.0, 0.05, 0.3, 1.0, 3.0, 10.0, 30.0, -0.3, -3.0, -30.0]  # A from the centre
HLINPROF_SETS = [(2, 3, 5800.0, 1e13, 1e17, 1e16), (2, 4, 6300.0, 1e14, 1e17, 1e16),
                 (2, 5, 10000.0, 1e15, 1e16, 1e15), (1, 3, 7000.0, 1e14, 1e16, 1e15),
                 (2, 22, 6300.0, 1e14, 1e17, 1e16), (3, 14, 5000.0, 1e13, 1e16, 1e15)]
HLINPROF_OFFSETS = [0.0, 0.1, 0.5, 2.0, 8.0, 30.0, -0.5, -8.0]
HBOP_WAVES = [3640.0, 3655.0, 3670.0, 3700.0, 3750.0, 3800.0, 3889.0, 3970.0, 4102.0, 4300.0,
              4861.0, 6563.0, 8000.0, 8600.0, 16000.0]


def compute():
    snap = {}
    snap["hbf"] = {"kind": "stable", "inputs": HBF_INPUTS, "values": H.hbf(HBF_INPUTS).tolist()}
    snap["wcalc"] = {"kind": "stable", "inputs": WCALC_CASES, "values": H.wcalc(WCALC_CASES).tolist()}
    snap["wcalc_hm88"] = {"kind": "stable", "excluded_volume_scale": 1.0, "inputs": WCALC_CASES,
                          "values": H.wcalc(WCALC_CASES, scale=1.0).tolist()}
    prof = {}
    for nlo, nup, T, ne, nh, nhe in PROFILE_SETS:
        w0 = H.line_centre(nlo, nup)
        p, q = H.profile(nlo, nup, T, ne, nh, nhe, H.doppler(T), w0 + np.array(PROFILE_OFFSETS))
        prof["%d-%d_T%g_ne%g" % (nlo, nup, T, ne)] = {"hprofl": p.tolist(), "hlinop": q.tolist()}
    snap["profile"] = {"kind": "legacy", "offsets_A": PROFILE_OFFSETS, "values": prof}
    hlp = {}
    for nlo, nup, T, ne, nh, nhe in HLINPROF_SETS:
        w0 = H.line_centre(nlo, nup)
        hlp["%d-%d_T%g_ne%g" % (nlo, nup, T, ne)] = H.hlinprof(
            nlo, nup, T, ne, nh, nhe, H.doppler(T), w0 + np.array(HLINPROF_OFFSETS)).tolist()
    snap["hlinprof"] = {"kind": "stable", "offsets_A": HLINPROF_OFFSETS, "values": hlp}
    lines = [(2, m, H.line_centre(2, m)) for m in range(3, 31)] + \
            [(3, m, H.line_centre(3, m)) for m in range(4, 31)] + \
            [(4, m, H.line_centre(4, m)) for m in range(5, 31)]
    # hbop: defaults (tabulated profiles, no excluded volume); hbop_hlinop_profiles: HLINOP for
    # all lines; hbop_original: HLINOP for all lines and the Hummer & Mihalas excluded volume,
    # i.e. the behaviour of the code before Nov 2026
    for key, scale, mode in (("hbop", None, None), ("hbop_hlinop_profiles", None, 0),
                             ("hbop_original", 1.0, 0)):
        hb = {}
        for label, nh, nhe, ne, T in HBOP_CONDITIONS:
            tot, cont = H.hbop(nh, nhe, ne, T, H.doppler(T), lines, HBOP_WAVES, scale=scale, profile=mode)
            hb[label] = {"total": tot.tolist(), "contin": cont.tolist()}
        snap[key] = {"kind": "legacy", "excluded_volume_scale": scale, "line_profile_mode": mode,
                     "waves_A": HBOP_WAVES,
                     "lines": "Balmer, Paschen, Brackett to n_up = 30, vacuum centres", "values": hb}
    return snap


if __name__ == "__main__":
    FILE.parent.mkdir(exist_ok=True)
    FILE.write_text(json.dumps(compute(), indent=1) + "\n")
    print("wrote", FILE)
