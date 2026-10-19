"""Reference bias and depth imbalance must not masquerade as robust associations."""
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from marker_program_baseline import marker_program_baseline
from marker_baseline_sensitivity import analyze


def main():
    n=400
    labels=np.array(['control']*n+['T']*160)
    guides=np.array(['C'+str(i//100) for i in range(n)]+['Tguide']*160)
    groups=np.array(['g0']*len(labels))
    ids=np.array(['cell'+str(i) for i in range(len(labels))])
    genes=np.array(['a','b','c','d'])
    panels={'P':['a','b','c','d']}
    # Perturbed cells match three control references; only the fourth reference differs.
    x=np.repeat(np.where(guides=='C0',1.,3.)[:,None],4,axis=1)
    p=marker_program_baseline(x,genes,labels,panels,strata=groups,cell_ids=ids,n_identity_shuffle=0)
    assert p['programs']['P']['modulators']==['T']
    s=analyze(x,genes,labels,ids,groups,guides,np.arange(len(labels))+100,p)
    assert s['control_reference_family_size']==4
    assert s['summary']['P']['all_control_references']==[]
    assert s['summary']['P']['control_reference_target_omission_robust']['C0']==['T']
    assert s['summary']['P']['control_reference_target_omission_robust']['C1']==[]
    # Score depends entirely on library depth; perturbation identity adds nothing within bins.
    umi=np.r_[np.tile([100.,1000.],n//2),[1000.]*160]
    x=np.repeat(np.where(umi==1000.,3.,1.)[:,None],4,axis=1)
    p=marker_program_baseline(x,genes,labels,panels,strata=groups,cell_ids=ids,n_identity_shuffle=0)
    assert p['programs']['P']['modulators']==['T']
    s=analyze(x,genes,labels,ids,groups,guides,umi,p)
    assert s['summary']['P']['depth_matched_target_omission_robust']==[]
    assert s['comparisons']['depth_matched'][0]['auc']==.5
    # A real within-bin shift survives both checks and the original heldout split is retained.
    x[labels=='T']+=3
    p=marker_program_baseline(x,genes,labels,panels,strata=groups,cell_ids=ids,n_identity_shuffle=0)
    s=analyze(x,genes,labels,ids,groups,guides,umi,p)
    assert s['summary']['P']['primary_and_depth_and_every_reference']==['T']
    assert len(s['per_gemgroup_auc_descriptive'])==1
    print('ALL MARKER-BASELINE SENSITIVITY TESTS PASSED')


if __name__=='__main__': main()
