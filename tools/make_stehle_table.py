#!/usr/bin/env python3
"""Build the binary Stehle & Hutcheon (1999) table read by HLINPROF (STEHLE) from the CDS catalogue.

Source: CDS J/A+AS/140/93 (VizieR VI/98A), "Extensive tabulations of Stark broadened hydrogen line
profiles", Stehle C., Hutcheon R., A&AS 140, 93 (1999).  Lyman (n_up 2-30), Balmer (3-30) and
Paschen (4-30) lines, 10 temperatures (2500 K - 1.259e6 K), up to 20 electron densities
(1e10 - 3.162e19 cm^-3; fewer for higher n_up because of the Inglis-Teller limit).  The Paschen
files are the corrected ones of 23 Nov 2007 (directory pa/, not pa.ori/).

Usage:
    python tools/make_stehle_table.py CDS_DIR OUTPUT
where CDS_DIR contains ly/, ba/ and pa/ with the profilN.dat and nraie.dat files of each line, e.g.
downloaded from https://cdsarc.cds.unistra.fr/ftp/J/A+AS/140/93/ (or a VizieR mirror).

Output: Fortran sequential unformatted file, little-endian, 4-byte record markers, format version 2:
  1  NWL, NTEMP, NNE, NLINE, IVERS(=2)                      int32
  2  (N1(L), M1(L), L=1,NLINE)                              int32
  3  (T(I), I=1,NTEMP)                         K            float32
  4  (XNE(J), J=1,NNE)                         log10 cm^-3  float32
  5  ((F0(J,L), J=1,NNE), L=1,NLINE)  normal Holtsmark field = A per alpha unit   float32
  6  ((WING(J,L), J=1,NNE), L=1,NLINE) asymptotic iStark * alpha^2.5              float32
  7  (((WIDTH(I,J,L), I=1,NTEMP), J=1,NNE), L=1,NLINE) HWHM incl. Doppler (alpha) float32
  8  (((NALPHA(I,J,L), I=1,NTEMP), J=1,NNE), L=1,NLINE)                            int32
  then for each line L: (((ALPHA(K,I,J), K=1,NWL), I=1,NTEMP), J=1,NNE)             float32
                        (((SPROF(K,I,J), K=1,NWL), I=1,NTEMP), J=1,NNE)             float32
SPROF is the Doppler-convolved profile (column iDoppler, units 1/alpha), as in the earlier
tables of HLINOP and SME.  NALPHA(I,J,L) = 0 marks a profile that is not tabulated: densities
above the line's Inglis-Teller limit, and at high densities the lowest temperatures (the tables
require r0/Debye < 1).  Missing profiles and points beyond NALPHA are zero-filled.  The format
version (IVERS = 2) differs from the earlier 4-integer header so that old readers stop with an
error instead of misreading the file.
"""
import re
import struct
import sys
from pathlib import Path

import numpy as np

SERIES = [("ly", 1, range(2, 31)), ("ba", 2, range(3, 31)), ("pa", 3, range(4, 31))]
NUM = re.compile(r"[-+]?\d+\.\d+E[-+]\d+")
TRIMMED = []   # (n_lo, n_up, density index, temperature index) of profiles cut at a padding row


def parse_profile_file(path):
    """Return dict with ne, f0, wing, T (list), width_dop (list), and per-T arrays (alpha, sdop)."""
    lines = path.read_text().splitlines()
    hdr = {}
    hdr["ne"] = float(NUM.findall(lines[1])[0])
    hdr["f0"] = float(NUM.findall(lines[2])[0])
    hdr["wing"] = float(NUM.findall(lines[3])[0])
    temps, widths, profs = [], [], []
    i = 4
    while i < len(lines):
        if lines[i].startswith("T(K)"):
            ts = [float(x) for x in NUM.findall(lines[i])]
            wd = [float(x) for x in NUM.findall(lines[i + 3])]   # "width incl.Dop."
            assert lines[i + 3].startswith("width incl"), (path, lines[i + 3])
            nt = len(ts)
            cols = [([], []) for _ in range(nt)]
            j = i + 5
            while j < len(lines) and lines[j].strip() and not lines[j].startswith("T(K)"):
                row = lines[j]
                alpha = float(row[:11])
                for k in range(nt):
                    field = row[11 + 25 * k: 11 + 25 * (k + 1)]
                    v = NUM.findall(field)
                    if v:
                        cols[k][0].append(alpha)
                        cols[k][1].append(float(v[0]))   # iDoppler (first number in the field)
                j += 1
            temps += ts
            widths += wd
            profs += [(np.array(a), np.array(s)) for a, s in cols]
            i = j
        else:
            i += 1
    return hdr, temps, widths, profs


def build(cds):
    lines = []
    for d, nlo, nups in SERIES:
        for nup in nups:
            tdir = cds / d / ("%s%02d" % (d, nup))
            nraie = int((tdir / "nraie.dat").read_text().split()[0])
            dens = [parse_profile_file(tdir / ("profil%d.dat" % k)) for k in range(1, nraie + 1)]
            lines.append((nlo, nup, dens))
    temps = max((dd[0][1] for _, _, dd in lines), key=len)   # full temperature grid
    nt = len(temps)
    nne = max(len(dd) for _, _, dd in lines)
    # common log10 density grid from the line with the most densities
    ref = max(lines, key=lambda x: len(x[2]))[2]
    xne = np.log10([h["ne"] for h, *_ in ref])
    nwl = max(len(p[0]) for _, _, dd in lines for _, _, _, pp in dd for p in pp)
    nl = len(lines)
    f0 = np.zeros((nl, nne), np.float32)
    wing = np.zeros((nl, nne), np.float32)
    width = np.zeros((nl, nne, nt), np.float32)
    nalpha = np.zeros((nl, nne, nt), np.int32)
    alpha = np.zeros((nl, nne, nt, nwl), np.float32)
    sprof = np.zeros((nl, nne, nt, nwl), np.float32)
    for L, (nlo, nup, dens) in enumerate(lines):
        for J, (h, ts, wd, pp) in enumerate(dens):
            assert abs(np.log10(h["ne"]) - xne[J]) < 1e-3, (nlo, nup, J, h["ne"])
            f0[L, J], wing[L, J] = h["f0"], h["wing"]
            for t, w, (a, s) in zip(ts, wd, pp):      # high densities omit the lowest temperatures
                I = int(np.argmin(np.abs(np.log10(temps) - np.log10(t))))
                assert abs(temps[I] / t - 1) < 1e-3, (nlo, nup, J, t)
                width[L, J, I] = w
                # Some CDS profiles end with a padding row printed with alpha = 0 but a non-zero
                # intensity; cut each profile at the first non-increasing alpha (and drop zero
                # intensities, which the log interpolation cannot use).
                bad = np.where(np.diff(a) <= 0)[0]
                if len(bad):
                    a, s = a[:bad[0] + 1], s[:bad[0] + 1]
                    TRIMMED.append((nlo, nup, J, I))
                keep = s > 0
                a, s = a[keep], s[keep]
                nalpha[L, J, I] = len(a)
                alpha[L, J, I, :len(a)] = a
                sprof[L, J, I, :len(a)] = s
    ids = np.array([(nlo, nup) for nlo, nup, _ in lines], np.int32)
    return dict(nwl=nwl, nt=nt, nne=nne, nl=nl, ids=ids, temps=np.array(temps, np.float32),
                xne=xne.astype(np.float32), f0=f0, wing=wing, width=width, nalpha=nalpha,
                alpha=alpha, sprof=sprof)


def write(tab, out):
    def rec(fh, arr):
        b = np.ascontiguousarray(arr).tobytes()
        fh.write(struct.pack("<i", len(b)) + b + struct.pack("<i", len(b)))
    with open(out, "wb") as fh:
        rec(fh, np.array([tab["nwl"], tab["nt"], tab["nne"], tab["nl"], 2], "<i4"))
        rec(fh, tab["ids"].astype("<i4"))                       # (N1,M1) pairs, L order
        rec(fh, tab["temps"].astype("<f4"))
        rec(fh, tab["xne"].astype("<f4"))
        rec(fh, tab["f0"].astype("<f4"))                        # C order (L,J) = Fortran (J,L)
        rec(fh, tab["wing"].astype("<f4"))
        rec(fh, tab["width"].astype("<f4"))                     # (L,J,I) = Fortran (I,J,L)
        rec(fh, tab["nalpha"].astype("<i4"))
        for L in range(tab["nl"]):
            rec(fh, tab["alpha"][L].astype("<f4"))              # (J,I,K) = Fortran (K,I,J)
            rec(fh, tab["sprof"][L].astype("<f4"))


if __name__ == "__main__":
    cds, out = Path(sys.argv[1]), Path(sys.argv[2])
    t = build(cds)
    write(t, out)
    print("profiles cut at a padding row:", len(TRIMMED))
    print("wrote %s: %d lines, NWL %d, NTEMP %d, NNE %d, log ne %.1f-%.1f, T %g-%g K"
          % (out, t["nl"], t["nwl"], t["nt"], t["nne"], t["xne"][0], t["xne"][-1], t["temps"][0], t["temps"][-1]))
