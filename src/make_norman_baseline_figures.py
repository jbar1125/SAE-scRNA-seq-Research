#!/usr/bin/env python3
"""Render the archived real Norman marker baseline without retraining or selection."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['svg.hashsalt']='norman-marker-baseline-v1'
matplotlib.rcParams['svg.fonttype']='none'
import matplotlib.pyplot as plt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,default=Path('data_norman_continuation/marker_baseline/gemgroup_seed0.json'))
    p.add_argument('--out-dir',type=Path,default=Path('docs/figures'))
    a=p.parse_args()
    r=json.loads(a.input.read_text())
    names=list(r['programs'])
    fig,(left,right)=plt.subplots(1,2,figsize=(13,5.8),gridspec_kw={'width_ratios':[1.4,1]})
    y=list(range(len(names)))
    n=[len(r['programs'][name]['modulators']) for name in names]
    robust=[len(r['programs'][name]['target_omission_robust_modulators']) for name in names]
    left.barh([v-.18 for v in y],n,height=.33,color='#A9BCCF',label='Primary score')
    left.barh([v+.18 for v in y],robust,height=.33,color='#1A85FF',label='Also passes target omission')
    left.set_yticks(y,[name+f" ({len(r['programs'][name]['eligible_markers'])} eligible markers)" for name in names])
    left.invert_yaxis()
    left.set_xlabel('Perturbations meeting q < .05 and AUC ≥ .55 (of 105)')
    left.set_title('A. External marker scores; held-out controls',loc='left',fontsize=11)
    left.legend(fontsize=8,loc='lower right')
    for i,(x0,x1) in enumerate(zip(n,robust)):
        left.text(max(x0,x1)+.5,i,f'{x0} / {x1}',va='center',fontsize=8)
    left.set_xlim(0,max(n)+11)
    subsets=['Ery_TF','Ery_Heme','Ery_Membrane']
    rows=[next(v for v in r['programs'][name]['per_perturbation'] if v['pert']=='KLF1') for name in subsets]
    right.bar([i-.18 for i in range(3)],[v['auc'] for v in rows],width=.33,color='#A9BCCF',label='Primary score')
    right.bar([i+.18 for i in range(3)],[v['target_omission']['auc'] for v in rows],width=.33,color='#E66100',label='Target omitted')
    right.axhline(.5,color='#555',linewidth=.8,linestyle=':')
    right.axhline(.55,color='#555',linewidth=.8,linestyle='--',label='Declared effect floor')
    right.set_xticks(range(3),['TF panel','Heme panel','Membrane panel'])
    right.set_ylabel('Gemgroup-stratified AUC: KLF1 vs held-out controls')
    right.set_ylim(.35,.71)
    right.set_title('B. KLF1: direct target vs downstream scores',loc='left',fontsize=11)
    right.legend(fontsize=8,loc='upper right')
    for ax in (left,right):
        ax.spines[['top','right']].set_visible(False)
    fig.suptitle('Norman CRISPRa: marker-score associations are not SAE-specific evidence',x=.05,ha='left',fontsize=14,y=.98)
    fig.text(.05,.885,'111,445 cells; 105 single perturbations; 5,926 discovery / 5,929 inference controls; eight technical gemgroups.\nBH across all 945 program–perturbation comparisons, separately for primary and omission sensitivity.',fontsize=9,color='#444')
    fig.text(.05,.065,'Panel B is a post hoc diagnostic of an existing candidate, not a separately selected hypothesis test.\nThree low-coverage panels are unscored. Technical partitions and control splitting do not establish biological replication.\nAssociations can reflect generic severity, cell state or compositional effects; no mechanism, novelty or SAE advantage is established.',fontsize=8.5,color='#444',linespacing=1.5)
    fig.subplots_adjust(left=.22,right=.98,top=.77,bottom=.23,wspace=.5)
    a.out_dir.mkdir(parents=True,exist_ok=True)
    for extension in ('png','svg'):
        path=a.out_dir/f'fig_norman_marker_baseline.{extension}'
        fig.savefig(path,dpi=160,bbox_inches='tight',metadata={'Date':None} if extension=='svg' else {})
        if extension=='svg':
            path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
    plt.close(fig)
    print('Rendered archived baseline counts and KLF1 diagnostic')


if __name__=='__main__':main()
