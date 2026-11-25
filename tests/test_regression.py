"""Compare the code with the stored snapshots (see snapshots.py for what they mean)."""
import json

import numpy as np
import pytest

import snapshots as S

# Outputs are single precision; allow for compiler and libm differences between platforms.
RTOL = 2e-5
ATOL = 1e-37

STORED = json.loads(S.FILE.read_text())
NOW = S.compute()


def _flat(x):
    if isinstance(x, dict):
        return [v for k in sorted(x) for v in _flat(x[k])]
    return list(np.ravel(x))


def _check(name):
    a, b = np.array(_flat(NOW[name]["values"])), np.array(_flat(STORED[name]["values"]))
    np.testing.assert_allclose(a, b, rtol=RTOL, atol=ATOL)


@pytest.mark.parametrize("name", ["hbf", "hlinprof", "wcalc", "wcalc_hm88"])
def test_stable_snapshot(name):
    _check(name)


@pytest.mark.legacy
@pytest.mark.parametrize("name", ["profile", "hbop", "hbop_hlinop_profiles", "hbop_original"])
def test_legacy_snapshot(name):
    _check(name)
