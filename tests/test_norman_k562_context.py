#!/usr/bin/env python3
"""Preserve all models, failed panels and perturbations in the expanded family."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from norman_k562_context import ALL_ARMS, correct_expanded


def main():
    records = {arm: {'programs': {str(i): {'per_perturbation': [
        {'pert': 'A', 'mw_p': .001 if arm == ALL_ARMS[0] and i == 0 else 1., 'auc': .6,
         'status': 'tested' if i == 0 else 'unmatched_program'}]} for i in range(10)},
        'identity': {str(i): {'empirical_p': .0001 if arm == ALL_ARMS[0] and i == 0 else 1.} for i in range(10)}}
        for arm in ALL_ARMS}
    out, a, i = correct_expanded(records)
    assert (a, i) == (80, 80)
    r = out[ALL_ARMS[0]]['programs']['0']['per_perturbation'][0]
    assert abs(r['context_bh_q'] - .08) < 1e-12 and not r['context_called']
    assert out[ALL_ARMS[0]]['identity']['0']['context_holm_p'] == .008
    assert 'context_bh_q' not in records[ALL_ARMS[0]]['programs']['0']['per_perturbation'][0]
    del records['direct_markers']
    try:
        correct_expanded(records)
        raise AssertionError('direct-score arm was silently dropped')
    except ValueError:
        pass
    print('ALL K562-CONTEXT TESTS PASSED (complete arms, unscored panels, expanded correction)')


if __name__ == '__main__':
    main()
