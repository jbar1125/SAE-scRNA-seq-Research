"""Interaction scale, technical composition and missing arms must be distinguished."""
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from norman_dual_marker import interaction,score_panel


def main():
    residual=np.tile(np.linspace(-1,1,200),2)
    labels=np.repeat(['control','A','B','A+B'],400)
    strata=np.tile(np.repeat(['s0','s1'],200),4);control=np.flatnonzero(labels=='control')
    additive=np.concatenate([residual+offset for offset in [0,1,1,2]])
    null=interaction(additive,labels,strata,'A+B',control)
    assert null['status']=='tested' and abs(null['contrast'])<1e-12 and null['p']>.99
    shifted=additive.copy();shifted[labels=='A+B']+=1
    positive=interaction(shifted,labels,strata,'A+B',control)
    assert np.isclose(positive['contrast'],1) and positive['p']<1e-10 and positive['interval_95'][0]>.2
    # Unequal batch composition creates a pooled apparent interaction, absent within strata.
    labs=[];blocks=[];scores=[]
    for name,counts in [('control',(100,300)),('A',(300,100)),('B',(300,100)),('A+B',(100,300))]:
        for j,n in enumerate(counts):
            labs.extend([name]*n);blocks.extend([str(j)]*n);scores.extend(np.linspace(-1,1,n)+4*j)
    labs=np.array(labs);blocks=np.array(blocks);scores=np.array(scores);controls=np.flatnonzero(labs=='control')
    pooled=interaction(scores,labs,np.repeat('all',len(labs)),'A+B',controls)
    conditioned=interaction(scores,labs,blocks,'A+B',controls)
    assert abs(pooled['contrast'])>1 and abs(conditioned['contrast'])<1e-12
    missing=interaction(scores,labs,blocks,'A+missing',controls)
    assert missing['status']=='insufficient_comparison_cells' and missing['p']==1
    assert score_panel(np.ones((100,2)),[0,1],{},np.arange(50)) is None
    print('ALL DUAL-MARKER INTERACTION TESTS PASSED')


if __name__=='__main__':main()
