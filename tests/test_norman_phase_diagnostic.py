"""Retain unmatched PCA hypotheses and both phase directions when joining families."""
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from norman_phase_diagnostic import correct_joint


def main():
    pca={'programs':{f'panel{i}':{'per_perturbation':[{'mw_p':1.}]} for i in range(9)}}
    rows=[{'primary':{'mw_p':.001,'status':'tested','auc':.7},
           'target_omission':{'mw_p':1.,'status':'tested','auc':.5}},
          {'primary':{'mw_p':1.,'status':'tested','auc':.3},
           'target_omission':{'mw_p':1.,'status':'tested','auc':.5}}]
    assert correct_joint(rows,pca)==11
    assert np.isclose(rows[0]['primary']['joint_holm_p'],.011)
    assert rows[0]['primary']['flagged'] and not rows[0]['target_omission']['flagged']
    assert not rows[1]['primary']['flagged']
    assert all(p['per_perturbation'][0]['joint_holm_primary']==1 for p in pca['programs'].values())
    print('ALL PHASE-DIAGNOSTIC FAMILY TESTS PASSED')


if __name__=='__main__':main()
