"""Interfaces for synthesis codes: table files at run time (HLINPROF_SET_FILES), which profile
HLINPROF returned (HLINPROF_SOURCE), informational messages (HLINPROF_SET_MESSAGES), and HBOP's
populations and bound-free opacity on their own (HBOP_LEVELS, HBOP_BF, HBOP_BF_DISSOLVED)."""
import os
import subprocess

import numpy as np
import pytest

import hltest as H
from cases import HBOP_CONDITIONS

LINES = [(2, 4), (2, 10), (3, 8), (1, 2), (4, 10), (2, 40)]


def _calls(T=6000.0, ne=1e14):
    out = []
    for nlo, nup in LINES:
        w0 = H.line_centre(nlo, nup)
        for dl in (0.0, 0.3, -2.0):
            out.append((nlo, nup, w0, T, ne, 1e16, 1e15, float(H.doppler(T)), w0 + dl))
    return out


def test_set_files_uses_the_given_table(tmp_path):
    """The same table under another name gives the same profiles; a missing file stops the
    program (so the name set at run time is the one used)."""
    other = tmp_path / "another_name.dat"
    os.symlink(H.ROOT / "stehle_long.dat", other)
    p0, _ = H.files(None, _calls())
    p1, _ = H.files(other, _calls())
    np.testing.assert_array_equal(p0, p1)
    r = H.files(tmp_path / "missing.dat", _calls()[:1], raw=True)
    assert "did not find" in r.stdout + r.stderr


def test_source_reports_the_fallback_to_hlinop():
    """1 for tabulated lines (Lyman alpha included), 0 where HLINOP is used: no table for the line
    (Brackett, n_up > 30) or a temperature below the table."""
    _, src = H.files(None, _calls())
    expect = {(2, 4): 1, (2, 10): 1, (3, 8): 1, (1, 2): 1, (4, 10): 0, (2, 40): 0}
    assert list(src) == [expect[l] for l in LINES for _ in range(3)]
    _, cold = H.files(None, _calls(T=1500.0))
    assert np.all(cold == 0)


def test_messages_can_be_switched_off():
    """By default HLINPROF prints which tables it uses; HLINPROF_SET_MESSAGES(0) silences that
    (results unchanged)."""
    calls = _calls()
    on = H.files(None, calls, raw=True)
    off = H.files(None, calls, messages=0, raw=True)
    assert "Will use Stehle" in on.stdout and "Using Hlinop" in on.stdout
    assert "Will use" not in off.stdout and "Using Hlinop" not in off.stdout
    num = lambda r: [ln for ln in r.stdout.splitlines() if ln.strip()[:1].isdigit() or ln.strip()[:1] == "-"]
    assert num(on) == num(off)


@pytest.mark.parametrize("cond", HBOP_CONDITIONS, ids=[c[0] for c in HBOP_CONDITIONS])
def test_hbop_bf_equals_the_continuum_of_hbop(cond):
    """HBOP_BF with the populations of HBOP_LEVELS is HBOP's continuous opacity, bit for bit; the
    LTE populations add up to NH."""
    _, nh, nhe, ne, T = cond
    waves = [900.0, 3640.0, 3660.0, 3700.0, 4861.0, 8200.0, 8300.0, 14600.0, 16000.0, 40000.0]
    (s, w1, w15), out = H.hbop_bf(nh, nhe, ne, T, waves)
    np.testing.assert_array_equal(out[:, 0], out[:, 1])
    assert abs(s - 1.0) < 1e-5
    assert 0.0 < w15 < w1 <= 1.0


@pytest.mark.parametrize("cond", HBOP_CONDITIONS, ids=[c[0] for c in HBOP_CONDITIONS])
def test_dissolved_part_of_the_bound_free_opacity(cond):
    """HBOP_BF_DISSOLVED counts only levels below their edge: nothing blueward of the Lyman edge,
    the n = 2 pseudo-continuum redward of the Balmer limit (part of HBOP_BF there), and at 1.6
    micron (levels n >= 5 still above their edges) the same with n <= 4 as with n <= 15."""
    _, nh, nhe, ne, T = cond
    _, out = H.hbop_bf(nh, nhe, ne, T, [900.0, 3700.0, 16000.0])
    cont, bf, d15, d4 = out.T
    assert d15[0] == 0.0
    assert 0.0 < d15[1] < bf[1]
    assert d15[2] == d4[2]


@pytest.mark.parametrize("src", ["hlinop.f", "hlinprof.f", "hbop.f"])
def test_compiles_without_legacy_options(src, tmp_path):
    """Codes such as SME compile these files without -std=legacy (newer gfortran then rejects,
    e.g., a scalar passed for an array argument)."""
    r = subprocess.run([H.FC, "-c", "-O2", "-I%s" % H.build_dir(), str(H.ROOT / src),
                        "-o", str(tmp_path / "x.o")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-2000:]
