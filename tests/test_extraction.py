"""Regression test pinning the synthetic example shipped with the original
release scripts (both the Python script and the Mathematica notebook produce
these values)."""

import numpy as np
import pytest

from rectcond import extract_tensor

# Example inputs from "extract conductivity tensor - python.py"
D1, D2 = 2.7, 1.2
R1 = 0.1636153846153846
R2 = 0.01676923076923077
R3 = -0.055923076923076916
R5 = 0.04815384615384615


@pytest.fixture(scope="module")
def result():
    return extract_tensor(D1, D2, R1, R2, R3, R5)


def test_matches_original_script_outputs(result):
    # Values printed by the original script (and the Mathematica notebook)
    assert result.r == pytest.approx(0.9124725968112841, rel=1e-8)
    assert result.sigma_hall == pytest.approx(0.9995256774753994, rel=1e-8)
    assert result.geometric_mean == pytest.approx(3.1613351020566602, rel=1e-8)
    assert result.z5 == pytest.approx(0.7095695295796434, rel=1e-8)
    assert result.a == pytest.approx(0.5000678518624317, rel=1e-6)
    assert result.alpha == pytest.approx(3.141368029833166, rel=1e-6)
    assert result.anisotropy_ratio == pytest.approx(0.6323714416917322, rel=1e-6)
    assert result.sigma_minus == pytest.approx(1.9991380361582496, rel=1e-6)
    assert result.sigma_plus == pytest.approx(4.999174367519501, rel=1e-6)


def test_recovers_round_number_tensor(result):
    # The shipped example is a synthetic test whose recovered tensor is
    # consistent with sigma_+ = 5, sigma_- = 2, sigma_H = 1 and principal
    # axes aligned with the sample edges (alpha = 0 mod pi).
    assert result.sigma_plus == pytest.approx(5.0, rel=1e-3)
    assert result.sigma_minus == pytest.approx(2.0, rel=1e-3)
    assert result.sigma_hall == pytest.approx(1.0, rel=1e-3)
    assert result.alpha % np.pi == pytest.approx(0.0, abs=1e-3) or \
        result.alpha % np.pi == pytest.approx(np.pi, abs=1e-3)


def test_internal_consistency(result):
    assert result.sigma_plus >= result.sigma_minus
    assert result.geometric_mean == pytest.approx(
        np.sqrt(result.sigma_plus * result.sigma_minus), rel=1e-9)
    assert 0.0 < result.r < 1.0
    assert 0.0 < result.a < 1.0
    assert 0.0 < result.z5 < 1.0
    assert 0.0 <= result.alpha < np.pi
