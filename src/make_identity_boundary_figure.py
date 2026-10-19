#!/usr/bin/env python3
"""Render the archived fixed boundary audit; do not rerun or select simulations."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def main():
    summary=json.loads((ROOT/'data_program_modulation_audit/boundary_audit/summary.json').read_text())
    order=[('identical_generic','Identical generic response\n(valid cell null)'),
           ('graded_generic_severity','Unequal generic severity\n(rank heterogeneity is real)'),
           ('confounded_batch_pooled','Identity confounded with batch\n(pooled cell shuffling)'),
           ('confounded_batch_conditioned','Same confounded data\n(condition on batch)'),
           ('clustered_null_cell_shuffle','Clustered population null\n(invalid cell shuffling)'),
           ('clustered_null_unit_means','Same clustered data\n(independent cluster means)')]
    with plt.rc_context({'font.size':10,'svg.hashsalt':'identity-boundaries-v1'}):
        fig,ax=plt.subplots(figsize=(10,5.8))
        for i,(name,label) in enumerate(order):
            r=summary['analyses'][name];x=100*r['any_rejection_rate'];lo,hi=np.array(r['exact_95_percent_interval'])*100
            if name=='confounded_batch_conditioned':
                ax.text(2,i,'No exchangeable identities; no test',va='center',fontsize=9,color='#555555')
            else:
                color='#176b64' if name in ['identical_generic','clustered_null_unit_means'] else '#a24438'
                ax.errorbar(x,i,xerr=[[x-lo],[hi-x]],fmt='o',color=color,capsize=3,markersize=6)
            ax.text(124,i,f"{r['any_rejection_count']}/100",va='center',ha='right',fontsize=9)
        ax.axvline(5,color='#777777',linestyle=':',linewidth=1,zorder=0)
        ax.set_yticks(range(len(order)),[label for _,label in order]);ax.invert_yaxis()
        ax.set_xlim(-3,128);ax.set_xticks([0,20,40,60,80,100])
        ax.set_xlabel('Simulations with any Holm rejection (%)')
        ax.set_title('Identity heterogeneity can reflect severity, batch or cell dependence',pad=16,fontsize=12)
        ax.spines[['top','right','left']].set_visible(False);ax.tick_params(axis='y',length=0,pad=12)
        ax.grid(axis='x',alpha=.15)
        fig.text(.02,.035,'100 seeds per analysis; 9 stipulated matches; 999 permutations when identifiable; exact 95% intervals.\nDotted 5% line is a nominal reference for valid nulls only. Paired analyses are not independent trials.',fontsize=9)
        fig.subplots_adjust(left=.36,right=.97,bottom=.19,top=.87)
        out=ROOT/'docs/figures/fig_identity_boundaries'
        fig.savefig(out.with_suffix('.png'),dpi=180)
        fig.savefig(out.with_suffix('.svg'),metadata={'Date':None})
        svg=out.with_suffix('.svg');svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
        plt.close(fig)


if __name__=='__main__':main()
