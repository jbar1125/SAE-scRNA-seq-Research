#!/usr/bin/env python3
"""Guard against selected-arm inference and incompatible recovered factor scores."""
import copy
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from norman_shared_inference import ARMS, correct_complete_family, require_reproduced


def fixture():
    return {arm: {'programs': {str(i): {'per_perturbation': [
        {'pert': 'A', 'mw_p': .001 if arm == ARMS[0] and i == 0 else 1.,
         'mw_q': .009, 'auc': .6, 'status': 'tested' if i == 0 else 'unmatched_program'}]}
        for i in range(9)}, 'identity_heterogeneity': {'programs': {
        str(i): {'empirical_p': .001 if arm == ARMS[0] and i == 0 else 1.,
                 'status': 'tested' if i == 0 else 'low_power_marker_panel'} for i in range(9)}}}
        for arm in ARMS}


def main():
    x = fixture()
    result, associations, identities = correct_complete_family(x)
    assert (associations, identities) == (54, 54)
    first = result[ARMS[0]]['programs']['0']['per_perturbation'][0]
    assert abs(first['joint_bh_q'] - .054) < 1e-12 and not first['joint_modulates']
    assert not result[ARMS[0]]['identity_heterogeneity']['programs']['0']['joint_heterogeneity']
    assert 'joint_bh_q' not in x[ARMS[0]]['programs']['0']['per_perturbation'][0]
    high_resolution = copy.deepcopy(x)
    high_resolution[ARMS[0]]['identity_heterogeneity']['programs']['0']['empirical_p'] = .0001
    corrected, _, _ = correct_complete_family(high_resolution)
    assert corrected[ARMS[0]]['identity_heterogeneity']['programs']['0']['joint_heterogeneity']
    broken = copy.deepcopy(x)
    del broken[ARMS[-1]]
    try:
        correct_complete_family(broken)
        raise AssertionError('selected arm family accepted')
    except ValueError:
        pass
    broken = copy.deepcopy(x)
    broken[ARMS[-1]]['programs']['0']['per_perturbation'] = []
    try:
        correct_complete_family(broken)
        raise AssertionError('missing hypotheses accepted')
    except ValueError:
        pass
    old = {'control_split': {'ids': [1, 2]}, 'programs': {'P': {
        'feature': 0, 'component_index': 0, 'orientation': 'positive', 'low_power': False,
        'top_associated_genes': ['G'], 'per_perturbation': [
            {'pert': 'A', 'status': 'tested', 'auc': .6, 'mw_p': .01}]}}}
    require_reproduced(copy.deepcopy(old), old)
    incompatible = copy.deepcopy(old)
    incompatible['programs']['P']['per_perturbation'][0]['auc'] = .601
    try:
        require_reproduced(incompatible, old)
        raise AssertionError('incompatible score recovery accepted')
    except ValueError:
        pass
    print('ALL SHARED-INFERENCE TESTS PASSED (complete families, resolution, score recovery)')


if __name__ == '__main__':
    main()
