#!/usr/bin/env python3
"""Target RNA induction and alternative names must not erase absent/negative cases."""
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from norman_target_rna import analyze_targets


def main():
    labels = np.array(['control']*40 + ['A']*40 + ['B']*40)
    blocks = np.array(['a', 'b']*60)
    x = np.zeros((120, 1))
    x[40:80] = 1
    result = analyze_targets(x, ['A'], labels, blocks, {'B': 'A'})
    assert len(result) == 3
    by = {(r['pert'], r['comparison']): r for r in result}
    assert by['A', 'direct_name']['target_rna_increase']
    assert by['B', 'direct_name']['mw_p'] == 1
    assert by['B', 'direct_name']['status'] == 'transcript_symbol_absent'
    assert not by['B', 'author_alternative']['target_rna_increase']
    assert all('joint_bh_q' in r for r in result)
    print('ALL TARGET-RNA TESTS PASSED (induction, missing transcript, negative alternative retained)')


if __name__ == '__main__':
    main()
