"""Write tests/data/stehle_cds_subset.json: a few profiles read directly from the CDS text files
of Stehle & Hutcheon (1999), J/A+AS/140/93, used to check HLINPROF at table nodes.

    python tests/make_cds_subset.py CDS_DIR

Each entry: line (n_lo, n_up), electron density, temperature, F0 (A per alpha unit), the
Doppler-convolved profile (alpha, iDoppler in 1/alpha) and the Stark-only profile (sstark, iStark
in 1/alpha) as printed in the CDS file, and the Stark-only HWHM (alpha units, file header).
"""
import re
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import make_stehle_table as M  # noqa: E402

# (series dir, n_lo, n_up, density file index, temperature) - table nodes
NODES = [("ly", 1, 2, 9, 10000.0), ("ly", 1, 2, 13, 10000.0), ("ly", 1, 3, 9, 10000.0), ("ly", 1, 10, 5, 5000.0),
         ("ba", 2, 4, 9, 5000.0), ("ba", 2, 10, 7, 10000.0), ("ba", 2, 22, 5, 5000.0),
         ("ba", 2, 30, 3, 10000.0), ("ba", 2, 5, 15, 19950.0),
         ("pa", 3, 10, 9, 5000.0), ("pa", 3, 14, 7, 10000.0), ("pa", 3, 25, 4, 5000.0)]
OUT = Path(__file__).resolve().parent / "data" / "stehle_cds_subset.json"


NUM = re.compile(r"[-+]?\d+\.\d+E[-+]\d+")


def header_width_dop(path, T):
    """Header HWHM including Doppler broadening ('width incl.Dop.', alpha units) at temperature T."""
    lines = path.read_text().splitlines()
    for i, ln in enumerate(lines):
        if ln.startswith("T(K)"):
            ts = [float(x) for x in NUM.findall(ln)]
            ws = [float(x) for x in NUM.findall(lines[i + 3])]
            for t, w in zip(ts, ws):
                if abs(t / T - 1) < 1e-3:
                    return w
    raise ValueError(path)


def stark_only(path, T):
    """Stark-only column and header HWHM ('width Stark') for temperature T."""
    lines = path.read_text().splitlines()
    i = 4
    while i < len(lines):
        if lines[i].startswith("T(K)"):
            ts = [float(x) for x in NUM.findall(lines[i])]
            ws = [float(x) for x in NUM.findall(lines[i + 2])]
            j, rows = i + 5, []
            while j < len(lines) and lines[j].strip() and not lines[j].startswith("T(K)"):
                rows.append(lines[j]); j += 1
            for c, t in enumerate(ts):
                if abs(t / T - 1) < 1e-3:
                    a, st = [], []
                    for r in rows:
                        v = NUM.findall(r[11 + 25 * c: 11 + 25 * (c + 1)])
                        if len(v) == 2:
                            al = float(r[:11])
                            if a and al <= a[-1]:
                                break
                            a.append(al); st.append(float(v[1]))
                    return ws[c], a, st
            i = j
        else:
            i += 1
    raise ValueError(path)


def main(cds):
    out = []
    for d, nlo, nup, k, T in NODES:
        h, ts, wd, pp = M.parse_profile_file(cds / d / ("%s%02d" % (d, nup)) / ("profil%d.dat" % k))
        i = min(range(len(ts)), key=lambda j: abs(ts[j] - T))
        assert abs(ts[i] / T - 1) < 1e-3
        a, s = pp[i]
        keep = [0]
        for j in range(1, len(a)):           # stop at a padding row (non-increasing alpha)
            if a[j] <= a[j - 1]:
                break
            keep.append(j)
        hw, ast, sst = stark_only(cds / d / ("%s%02d" % (d, nup)) / ("profil%d.dat" % k), T)
        out.append(dict(nlo=nlo, nup=nup, ne=h["ne"], T=ts[i], f0=h["f0"], wave0=None,
                        alpha=[float(a[j]) for j in keep], sdop=[float(s[j]) for j in keep],
                        hwhm_stark_alpha=hw, alpha_stark=ast, sstark=sst,
                        hwhm_dop_alpha=header_width_dop(cds / d / ("%s%02d" % (d, nup)) / ("profil%d.dat" % k), T)))
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print("wrote", OUT, len(out), "profiles")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
