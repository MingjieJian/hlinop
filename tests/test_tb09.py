"""Stehle -> TB09 blending of the tabulated Balmer profiles (HLINOP_SET_TB09; hlinprof.f STEHL2/HTB09P).

HTB09P reads the TB09 table as distributed with SYNSPEC (data/tremblay.dat, Lemke format).  The table is not in
this repository, so most tests use a synthetic table in the same format with a known profile shape (two blocks,
headers over several lines, a zero-padded label as in tremblay.dat).  They check: default off (STEHL2 equals
the Stehle evaluator STEHLS bit for bit), the weight w = 3t^2 - 2t^3 across the window, pure TB09 above it
(against an independent evaluation of the synthetic profile), the mirror point, lines outside H10-H22, cache
clearing on a setting change, the argument checks, and that a file in another format stops with a message.
With HLINOP_TB09_FILE set to SYNSPEC's tremblay.dat, pure TB09 is also checked against an independent reading
of that file.
"""
from pathlib import Path
import os

import numpy as np
import pytest

import hltest as H

LA0, DLA, MP = -3.0, 0.1, 55          # log10 alpha grid
LNE0, DNE, MNE = 10.0, 0.5, 17        # log10 n_e grid
LT0, DLT, MT = np.log10(2500.0), np.log10(2.0), 7


def shape(alpha, nup):
    """Synthetic S(alpha): Gaussian core plus an alpha^-2.5 wing (not normalised)."""
    a0 = 1.0 + 0.05 * (nup - 10)
    return np.exp(-(alpha / a0) ** 2) + 0.05 * (1.0 + alpha) ** -2.5


def write_lemke(path, blocks):
    """Write a TB09 table in the Lemke/SYNSPEC format.  blocks: list of lists of (nl, nu, header, rows),
    header = (la0, lne0, lt0, dla, dne, dlt, mp, mne, mt), rows = [(flag, logS[mp])] * mne*mt (n_e outer)."""
    out = [f"{len(blocks):6d}"]
    for block in blocks:
        out.append(f"{len(block):12d}")
        for nl, nu, (la0, lne0, lt0, dla, dne, dlt, mp, mne, mt), _ in block:
            out += [f"{nl:12d}{nu:12d}{la0:24.15f}{lne0:24.15f}", f"{lt0:24.15f}{dla:24.15f}{dne:24.15f}",
                    f"{dlt:24.15f}{mp:16d}{mne:12d}{mt:12d}"]
        for nl, nu, (*_, mp, mne, mt), rows in block:
            out.append(f"nl={nl:02d} nu={nu:2d}")
            for flag, vals in rows:
                v = [f"{float(flag):8.4f}"] + [f"{x:8.4f}" for x in vals]
                out += ["".join(v[i:i + 9]) for i in range(0, len(v), 9)]
    Path(path).write_text("\n".join(out) + "\n")


@pytest.fixture(scope="module")
def table(tmp_path_factory):
    path = tmp_path_factory.mktemp("tb09") / "tb09_synthetic.dat"
    la = LA0 + DLA * np.arange(MP)
    head = (LA0, LNE0, LT0, DLA, DNE, DLT, MP, MNE, MT)
    lyman = [(1, 3, head, [(0, np.full(MP, -9.0))] * (MNE * MT))]
    balmer = [(2, nup, head, [(0, np.round(np.log10(shape(10 ** la, nup)), 4))] * (MNE * MT)) for nup in range(10, 23)]
    write_lemke(path, [lyman, balmer])
    return path


def tb09_expected(nup, ne, dlam):
    """Independent value of the synthetic profile per A: S(alpha)/F0 / (2 int_0^inf S dalpha), with the
    table's interpolation (linear in log alpha, alpha^-2.5 beyond the last point)."""
    la = LA0 + DLA * np.arange(MP); ls = np.round(np.log10(shape(10 ** la, nup)), 4)
    x = np.linspace(LA0, la[-1], 200001); s = 10 ** np.interp(x, la, ls)
    area = 10 ** LA0 * 10 ** ls[0] + np.trapezoid(s * 10 ** x * np.log(10), x) + 10 ** ls[-1] * 10 ** la[-1] / 1.5
    f0 = 1.25e-9 * ne ** (2 / 3); xa = np.log10(max(abs(dlam) / f0, 10 ** LA0))
    val = 10 ** np.interp(xa, la, ls) if xa <= la[-1] else 10 ** (ls[-1] - 2.5 * (xa - la[-1]))
    return val / f0 / (2 * area)


def call(nup, T, lne, dlam, iboth=0):
    return ("stehl2", nup, H.line_centre(2, nup), T, 10.0 ** lne, dlam, iboth)


def test_default_off_equals_stehle():
    r = H.tb09([call(n, 7000.0, lne, d) for n in (12, 18, 22) for lne in (12.0, 14.0) for d in (0.0, 0.7, 5.0)])
    np.testing.assert_array_equal(r[:, 0], r[:, 2])


@pytest.mark.parametrize("nup", [12, 18, 22])
def test_weight_across_the_window(table, nup):
    """Below the window: Stehle exactly; above: TB09; inside: (1-w) Stehle + w TB09, w = 3t^2-2t^3."""
    T, d = 7000.0, 0.6
    lim = H.tb09([call(nup, T, 11.0, d)])[0, 4]                   # last tabulated Stehle log n_e
    lnes = [lim - 2.6, lim - 2.0, lim - 1.9, lim - 1.75, lim - 1.6, lim - 1.5, lim - 1.0]
    blend = H.tb09([("set", 1, 2.0, 1.5, table)] + [call(nup, T, x, d) for x in lnes])
    pure = H.tb09([("set", 1, 20.0, 19.0, table)] + [call(nup, T, x, d) for x in lnes])      # w = 1
    for x, b, p in zip(lnes, blend, pure):
        t = min(1.0, max(0.0, (x - (lim - 2.0)) / 0.5)); w = t * t * (3 - 2 * t)
        np.testing.assert_allclose(b[0], (1 - w) * b[2] + w * p[0], rtol=1e-5)   # REAL*4 arithmetic
        if w == 0:
            assert b[0] == b[2]


def test_pure_tb09_matches_the_table(table):
    nup, T = 16, 7000.0
    cases = [(lne, d) for lne in (12.0, 13.5, 15.0) for d in (0.0, 0.3, 2.0, 30.0)]
    r = H.tb09([("set", 1, 20.0, 19.0, table)] + [call(nup, T, lne, d) for lne, d in cases])
    exp = [tb09_expected(nup, 10 ** lne, d) for lne, d in cases]
    np.testing.assert_allclose(r[:, 0], exp, rtol=2e-3)


def test_mirror_point(table):
    """IBOTH = 1: the mirror value is blended with the same TB09 value (TB09 depends on |dlambda| only)."""
    nup, T = 18, 7000.0
    lim = H.tb09([call(nup, T, 11.0, 0.5)])[0, 4]
    r = H.tb09([("set", 1, 2.0, 1.5, table), call(nup, T, lim - 1.75, 3.0, 1), call(nup, T, lim - 1.75, -3.0, 0)])
    np.testing.assert_allclose(r[0, 1], r[1, 0], rtol=1e-6)


def test_lines_outside_h10_h22_unchanged(table):
    cs = [call(n, 7000.0, lne, 0.5) for n in (8, 9, 23, 26) for lne in (12.0, 14.0)]
    off = H.tb09(cs)
    on = H.tb09([("set", 1, 20.0, 19.0, table)] + cs)
    np.testing.assert_array_equal(on[:, 0], off[:, 0])


def test_setting_change_clears_the_profile_cache(table):
    """HLINPROF caches convolved profiles; a later HLINOP_SET_TB09 must not return the cached Stehle one."""
    nup = 18; w0 = H.line_centre(2, nup); c = ("hlinprof", nup, w0, 7000.0, 1e14, H.doppler(7000.0), 0.4)
    r = H.tb09([c, ("set", 1, 20.0, 19.0, table), c])
    fresh = H.tb09([("set", 1, 20.0, 19.0, table), c])
    assert r[0, 0] != r[1, 0]
    assert r[1, 0] == fresh[0, 0]


def test_argument_checks(tmp_path):
    bad = H.tb09([("set", 1, 1.5, 2.0, tmp_path / "x"), call(18, 7000.0, 13.0, 0.5)], raw=True)
    assert "need DLO > DHI" in bad.stdout + bad.stderr
    missing = H.tb09([("set", 1, 2.0, 1.5, tmp_path / "missing.dat"), call(18, 7000.0, 13.0, 0.5)], raw=True)
    assert "TB09 table not found" in missing.stdout + missing.stderr


def test_other_format_stops(tmp_path):
    """A file that is not in the SYNSPEC/Lemke format (here the earlier compact copy) stops with a message."""
    bad = tmp_path / "compact.dat"
    bad.write_text("1\n2 18 3 2 2 -3 0.1 10 0.5 3.4 0.3\n" + "0 -1 -2 -3\n" * 4)
    r = H.tb09([("set", 1, 20.0, 19.0, bad), call(18, 7000.0, 13.0, 0.5)], raw=True)
    assert "cannot read the TB09 table" in r.stdout + r.stderr


def read_lemke(path):
    """Independent reader of the SYNSPEC table: {(nl, nu): (header, logS[mne, mt, mp])}."""
    tok = Path(path).read_text().split("\n")
    i, out = 1, {}
    def nums(i, n):
        v = []
        while len(v) < n:
            v += [float(x.replace("D", "E")) for x in tok[i].split()]; i += 1
        return v, i
    nblk = int(tok[0])
    for _ in range(nblk):
        while not tok[i].strip(): i += 1
        nb = int(tok[i]); i += 1; heads = []
        for _ in range(nb):
            v, i = nums(i, 11); heads.append(v)
        for v in heads:
            nl, nu, la0, lne0, lt0, dla, dne, dlt, mp, mne, mt = v; mp, mne, mt = int(mp), int(mne), int(mt)
            while not tok[i].strip(): i += 1
            i += 1
            rows = []
            for _ in range(mne * mt):
                r, i = nums(i, mp + 1); rows.append(r[1:])
            out[(int(nl), int(nu))] = ((la0, lne0, lt0, dla, dne, dlt), np.array(rows).reshape(mne, mt, mp))
    return out


@pytest.mark.skipif(not os.environ.get("HLINOP_TB09_FILE"), reason="HLINOP_TB09_FILE (SYNSPEC tremblay.dat) not set")
def test_real_table_pure_tb09():
    """Pure TB09 from SYNSPEC's tremblay.dat at table nodes (T, n_e) against an independent evaluation."""
    path = os.environ["HLINOP_TB09_FILE"]
    tab = read_lemke(path)
    assert len(tab) == 69 and (2, 22) in tab and (4, 14) in tab
    cases = [(nup, it, ine, d) for nup in (10, 14, 18, 22) for it, ine in ((2, 8), (4, 10)) for d in (0.0, 0.4, 3.0)]
    cmds = [("set", 1, 20.0, 19.0, path)]
    exp = []
    for nup, it, ine, d in cases:
        (la0, lne0, lt0, dla, dne, dlt), ls3 = tab[(2, nup)]
        T, ne = 10 ** (lt0 + dlt * it), 10 ** (lne0 + dne * ine)
        cmds.append(call(nup, T, np.log10(ne), d))
        ls = ls3[ine, it]; la = la0 + dla * np.arange(len(ls))
        x = np.linspace(la0, la[-1], 200001); s = 10 ** np.interp(x, la, ls)
        area = 10 ** la0 * 10 ** ls[0] + np.trapezoid(s * 10 ** x * np.log(10), x) + 10 ** ls[-1] * 10 ** la[-1] / 1.5
        f0 = 1.25e-9 * ne ** (2 / 3); xa = np.log10(max(abs(d) / f0, 10 ** la0))
        val = 10 ** np.interp(xa, la, ls) if xa <= la[-1] else 10 ** (ls[-1] - 2.5 * (xa - la[-1]))
        exp.append(val / f0 / (2 * area))
    r = H.tb09(cmds)
    np.testing.assert_allclose(r[:, 0], exp, rtol=2e-3)
