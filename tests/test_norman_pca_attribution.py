#!/usr/bin/env python3
"""Projection can rise when all selected marker genes fall; strata must be retained."""
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from norman_pca_attribution import contributions, weighted_score_difference


def main():
    # A control-correlated marker decreases, but a negative-weight other gene
    # decreases further. This establishes a possible algebraic explanation only.
    x = np.array([[2., 4.], [2., 4.], [1., 1.], [1., 1.],
                  [12., 14.], [12., 14.], [11., 11.], [11., 11.]])
    w = np.array([1., -1.])
    strata = np.array(['a']*4 + ['b']*4)
    pert = np.array([2, 3, 6, 7]);ctrl = np.array([0, 1, 4, 5])
    c = contributions(x, w, pert, ctrl, strata)
    assert np.allclose(c, [-1., 3.])
    assert np.isclose(c.sum(), weighted_score_difference(x @ w, pert, ctrl, strata))
    # Cell-count imbalance can create a pooled shift from batch offsets alone.
    y = np.array([[0., 0.]]*8 + [[10., 0.]]*8)
    unbalanced_strata = np.array(['a']*8 + ['b']*8)
    unbalanced_pert = np.array([0, 1, 8, 9, 10, 11]);controls = np.array([2, 3, 4, 5, 12, 13])
    assert not np.isclose((y @ w)[unbalanced_pert].mean(), (y @ w)[controls].mean())
    assert np.allclose(contributions(y, w, unbalanced_pert, controls, unbalanced_strata), 0)
    print('ALL PCA-ATTRIBUTION TESTS PASSED (signed mean identity and stratum imbalance)')


if __name__ == '__main__':
    main()
