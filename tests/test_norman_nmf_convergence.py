#!/usr/bin/env python3
"""Actual NMF warm optimization plus complete seven-arm correction behavior."""
from pathlib import Path
import sys
import warnings
import numpy as np
from sklearn.decomposition import NMF

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from norman_nmf_convergence import warm_restart, joint_followup, NEW_ARM
from norman_shared_inference import ARMS


def main():
    rng = np.random.default_rng(7)
    x = (rng.uniform(.1, 1., (30, 2)) @ rng.uniform(.1, 1., (2, 8))).astype(np.float32)
    first = NMF(n_components=2, init='random', solver='mu', max_iter=2, random_state=0)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        w = first.fit_transform(x)
    old_w, old_h = w.copy(), first.components_.copy()
    settings = {'n_components': 2, 'init': 'custom', 'solver': 'mu', 'beta_loss': 'frobenius',
                'tol': .0001, 'max_iter': 500, 'random_state': 0}
    scores, h, info = warm_restart(x, w, first.components_, settings)
    assert np.array_equal(w, old_w) and np.array_equal(first.components_, old_h)
    assert np.linalg.norm(x - scores @ h) < np.linalg.norm(x - old_w @ old_h)
    assert info['iterations'] > 0 and np.isfinite(scores).all() and (scores >= 0).all()
    def fixture(p):
        return {'programs': {str(i): {'low_power': i != 0, 'per_perturbation': [
            {'pert': 'A', 'mw_p': p if i == 0 else 1., 'auc': .6, 'status': 'tested' if i == 0 else 'unmatched_program'}]}
            for i in range(9)}, 'identity_heterogeneity': {'programs': {
                str(i): {'empirical_p': p if i == 0 else 1., 'status': 'tested' if i == 0 else 'no_marker_match'} for i in range(9)}}}
    original = {arm: fixture(1.) for arm in ARMS}
    records, a, i = joint_followup(original, fixture(.0008))
    assert (a, i) == (63, 63)
    row = records[NEW_ARM]['programs']['0']['per_perturbation'][0]
    assert not row['seven_arm_modulates'] and abs(row['seven_arm_joint_q'] - .0504) < 1e-12
    assert 'seven_arm_joint_q' not in original[ARMS[0]]['programs']['0']['per_perturbation'][0]
    del original[ARMS[-1]]
    try:
        joint_followup(original, fixture(.0008))
        raise AssertionError('selected-arm correction accepted')
    except ValueError:
        pass
    print('ALL NMF-CONVERGENCE TESTS PASSED (actual warm fit, initial artifacts, seven-arm family)')


if __name__ == '__main__':
    main()
