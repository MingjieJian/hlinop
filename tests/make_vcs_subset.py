"""Write tests/data/lemke_vcs_brackett.json: reference HWHM and peak of Stark-only Brackett profiles
from the VCS tables of Lemke (1997, A&AS 122, 285; CDS J/A+AS/122/285, file 'bracknd', without
Doppler broadening), used to check STARK1 for Brackett lines (no Stehle tables exist for them).

    python tests/make_vcs_subset.py PATH_TO_bracknd

Br11 and Br12 are not used: their VCS widths are about 4 times narrower than the neighbouring
lines and Griem theory (a defect of the table).
"""
import json
import sys
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent / "data" / "lemke_vcs_brackett.json"
CASES = [(nup, T, ne) for nup in (13, 16, 20) for T in (5000.0, 10000.0) for ne in (1e12, 1e13, 1e14)]


def read(fn):
    tok = open(fn).read().split()
    p = 0
    nline = int(tok[p])
    p += 1
    hdr = []
    for _ in range(nline):
        v = tok[p:p + 11]
        p += 11
        hdr.append(dict(nl=int(v[0]), nu=int(v[1]), la0=float(v[2]), lne0=float(v[3]), lt0=float(v[4]),
                        dla=float(v[5]), dne=float(v[6]), dt=float(v[7]), mp=int(v[8]), mne=int(v[9]),
                        mt=int(v[10])))
    out = {}
    for h in hdr:
        p += 2
        n = (h["mp"] + 1) * h["mt"] * h["mne"]
        a = np.array(tok[p:p + n], float).reshape(h["mne"], h["mt"], h["mp"] + 1)
        p += n
        h["logS"] = a[:, :, 1:]
        h["logalpha"] = h["la0"] + h["dla"] * np.arange(h["mp"])
        h["logne"] = h["lne0"] + h["dne"] * np.arange(h["mne"])
        h["logT"] = h["lt0"] + h["dt"] * np.arange(h["mt"])
        out[(h["nl"], h["nu"])] = h
    return out


def profile(h, T, ne, dl):
    """Profile per A (bilinear in log ne, log T on log S; alpha^-2.5 beyond the table)."""
    def idx(g, x):
        i = int(np.clip(np.searchsorted(g, x) - 1, 0, len(g) - 2))
        return i, (x - g[i]) / (g[i + 1] - g[i])
    j, u = idx(h["logne"], np.log10(ne))
    i, t = idx(h["logT"], np.log10(T))
    L = h["logS"]
    ls = (1 - u) * ((1 - t) * L[j, i] + t * L[j, i + 1]) + u * ((1 - t) * L[j + 1, i] + t * L[j + 1, i + 1])
    F0 = 1.25e-9 * ne ** (2 / 3)
    x = np.log10(np.maximum(np.abs(dl) / F0, 1e-30))
    y = np.interp(x, h["logalpha"], ls)
    hi = x > h["logalpha"][-1]
    y[hi] = ls[-1] - 2.5 * (x[hi] - h["logalpha"][-1])
    return 10 ** y / F0


def main(fn):
    tab = read(fn)
    out = []
    for nup, T, ne in CASES:
        x = np.linspace(0, 2000, 2000001)
        p = profile(tab[(4, nup)], T, ne, x)
        half = p[0] / 2
        k = int(np.argmax(p < half))
        hw = float(np.interp(half, [p[k], p[k - 1]], [x[k], x[k - 1]]))
        out.append(dict(nlo=4, nup=nup, T=T, ne=ne, hwhm_A=hw, peak_per_A=float(p[0])))
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print("wrote", OUT, len(out), "cases")


if __name__ == "__main__":
    main(sys.argv[1])
