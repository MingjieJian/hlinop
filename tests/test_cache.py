"""The convolution caches of STARK1C, HPCONV and HLPCONV (hlinop.f, hlinprof.f).

A synthesis code calls the profiles wavelength by wavelength, and at each wavelength for every
depth and line, so consecutive calls belong to different states (line, T, n_e, ...).  Each state is
built once and kept (hash index HLCKEY, up to 6144 states stored densely); the result must not
depend on the order of the calls.
"""

import numpy as np

import hltest as H

LINES = [(2, 4), (2, 8), (2, 17), (3, 12), (2, 35), (4, 14)]   # tabulated and HLINOP lines
DEPTHS = [(4000.0 + 400.0 * i, 10.0 ** (10.0 + 0.4 * i)) for i in range(15)]
OFFSETS = np.concatenate([-np.geomspace(1e-3, 30.0, 12), [0.0], np.geomspace(1e-3, 30.0, 12)])


def _call(line, depth, dl):
    nlo, nup = line
    T, ne = depth
    w0 = H.line_centre(nlo, nup)
    return (nlo, nup, w0, T, ne, 1e16, 1e15, float(H.doppler(T)), w0 + dl)


def test_results_do_not_depend_on_the_call_order():
    """Grouped by state (each state built once and then used) and in the order of a synthesis
    code (wavelength, then depth, then line: the state changes at every call): identical."""
    grouped = [(l, d, x) for l in LINES for d in DEPTHS for x in OFFSETS]
    synthesis = [(l, d, x) for x in OFFSETS for d in DEPTHS for l in LINES]
    p1, q1 = H.states([_call(*c) for c in grouped])
    p2, q2 = H.states([_call(*c) for c in synthesis])
    i2 = {c: i for i, c in enumerate(synthesis)}
    order = [i2[c] for c in grouped]
    np.testing.assert_array_equal(p1, p2[order])
    np.testing.assert_array_equal(q1, q2[order])


def test_cache_index():
    """HLCKEY (16 hash slots, 12 entries): new keys miss and are stored in the next entry (dense
    storage), repeated keys hit their entry; when all entries are used the next new key empties
    the cache."""
    keys = [1000.0 + 37.5 * i for i in range(12)]
    first = H.hlckey(keys)
    assert np.all(first[:, 1] == 0) and list(first[:, 2]) == list(range(1, 13))
    assert list(first[:, 0]) == list(range(1, 13))
    both = H.hlckey(keys + keys)
    assert np.all(both[12:, 1] == 1)
    np.testing.assert_array_equal(both[12:, 0], first[:, 0])
    over = H.hlckey(keys + [5000.0] + keys[:1])
    assert over[12, 1] == 0 and over[12, 2] == 1          # emptied, then the new key stored
    assert over[13, 1] == 0 and over[13, 2] == 2          # the earlier key is gone


def test_clear_cache():
    """HLINOP_CLEAR_CACHE empties the cache at its next use (entries restart at 1)."""
    keys = [1000.0 + 37.5 * i for i in range(5)]
    out = H.hlckey(keys + [-1.0] + keys[:2])
    assert list(out[:5, 0]) == [1, 2, 3, 4, 5]
    assert list(out[5:, 1]) == [0, 0] and list(out[5:, 2]) == [1, 2]


def test_clearing_the_cache_does_not_change_results():
    """Profiles computed after HLINOP_CLEAR_CACHE (rebuilt) equal those from the cache."""
    calls = [_call(l, d, x) for x in OFFSETS[::4] for d in DEPTHS[::5] for l in LINES]
    clear = (0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    p, q = H.states(calls + [clear] + calls)
    n = len(calls)
    np.testing.assert_array_equal(p[:n], p[n:])
    np.testing.assert_array_equal(q[:n], q[n:])

