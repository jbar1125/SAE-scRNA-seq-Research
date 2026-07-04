"""Unit tests for the participation-ratio baseline math.
Run: python tests/test_baselines_logic.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import baselines_nmf_pca as b


def test_participation_ratio_known_cases():
    assert abs(b.participation_ratio(np.ones(50)) - 50) < 1e-6      # isotropic -> D
    assert abs(b.participation_ratio([9.0, 0, 0, 0]) - 1) < 1e-6    # rank-1 -> 1
    assert abs(b.participation_ratio([2, 2, 2, 2, 0, 0]) - 4) < 1e-6
    assert np.isnan(b.participation_ratio([0, 0, 0]))
    print("PASS test_participation_ratio_known_cases")


def test_pr_of_representation_rank():
    rng = np.random.default_rng(0)
    # a rank-3 representation embedded in 20 dims -> PR close to 3
    Z = rng.standard_normal((500, 3)) @ rng.standard_normal((3, 20))
    pr = b.pr_of_representation(Z)
    assert 2.0 < pr < 4.0, pr
    print("PASS test_pr_of_representation_rank")


if __name__ == "__main__":
    test_participation_ratio_known_cases()
    test_pr_of_representation_rank()
    print("\nALL BASELINE LOGIC TESTS PASSED")
