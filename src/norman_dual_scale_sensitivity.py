#!/usr/bin/env python3
"""Declared log/linear score comparison; transformation and standardization both vary."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import subprocess

import numpy as np
from statsmodels.stats.multitest import multipletests

from marker_program_baseline import normalize_counts_to_markers,split_controls,_file_sha256
from norman_dual_marker import analyze

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--adata',required=True,type=Path)
    parser.add_argument('--primary',required=True,type=Path)
    parser.add_argument('--out',required=True,type=Path)
    args=parser.parse_args()
    if args.out.exists():parser.error('choose a fresh output path')
    plan=ROOT/'config/norman_dual_scale_sensitivity_plan.json'
    base=ROOT/'data_norman_continuation/marker_baseline/gemgroup_seed0.json'
    sources=[Path(__file__).resolve(),plan,base,ROOT/'src/norman_dual_marker.py',ROOT/'src/marker_program_baseline.py',
             ROOT/'config/human_markers.json',args.primary]
    before={str(p):_file_sha256(p) for p in sources}
    primary=json.loads(args.primary.read_text());baseline=json.loads(base.read_text())
    import anndata as ad
    a=ad.read_h5ad(args.adata);labels=a.obs['perturbation'].astype(str).to_numpy();ids=a.obs_names.astype(str).to_numpy()
    if a.obs['gemgroup'].isna().any():raise ValueError('missing gemgroup')
    strata=a.obs['gemgroup'].astype(str).to_numpy()
    panels=json.loads((ROOT/'config/human_markers.json').read_text())['marker_sets']
    log_expression,genes,_=normalize_counts_to_markers(a.X,a.var_names.to_numpy(),panels);del a
    # Reproduce every original primary/omission estimate before checking another scale.
    reproduced=analyze(log_expression,genes,labels,strata,ids,baseline)
    for name,r in primary['programs'].items():
        if reproduced['programs'][name]!=r:raise ValueError('primary dual output failed exact reproduction: '+name)
    linear=np.expm1(log_expression);discovery,_=split_controls(labels,strata=strata,seed=0)
    for j,g in enumerate(genes):
        r=baseline['discovery_gene_parameters'][g];sd=float(linear[discovery,j].std(ddof=0))
        r.update(mean=float(linear[discovery,j].mean()),population_sd=sd,scale=max(sd,.1),
                 n_nonzero=int(np.count_nonzero(linear[discovery,j])))
    linear_result=analyze(linear,genes,labels,strata,ids,baseline)
    linear_result['baseline_gene_parameters_reproduced']='linear discovery-control parameters, with original marker eligibility retained'
    keys=[(name,i) for name in primary['programs'] for i in range(len(primary['pairs']))]
    adjusted={}
    for which in ['primary','target_omission']:
        values=[result['programs'][name]['comparisons'][i][which]['p'] for result in [primary,linear_result] for name,i in keys]
        q=multipletests(values,method='holm')[1];offset=len(keys)
        for j,(name,i) in enumerate(keys):
            row=adjusted.setdefault((name,i),{})
            row[which]={'log_joint_holm_p':float(q[j]),'linear_joint_holm_p':float(q[offset+j])}
    comparisons=[];counts={name:0 for name in primary['programs']}
    for name,i in keys:
        log_row=primary['programs'][name]['comparisons'][i];linear_row=linear_result['programs'][name]['comparisons'][i]
        assert log_row['pair']==linear_row['pair']
        values=[log_row[w] for w in ['primary','target_omission']]+[linear_row[w] for w in ['primary','target_omission']]
        ps=[adjusted[(name,i)][w][scale+'_joint_holm_p'] for scale in ['log','linear'] for w in ['primary','target_omission']]
        good=all(r['status']=='tested' and abs(r['contrast'])>=.2 for r in values) and all(p<.05 for p in ps)
        if good:good=len({int(np.sign(r['contrast'])) for r in values})==1
        counts[name]+=bool(good)
        comparisons.append({'program':name,'pair':log_row['pair'],'joint':adjusted[(name,i)],
            'log_primary_contrast':log_row['primary']['contrast'],'linear_primary_contrast':linear_row['primary']['contrast'],
            'two_scale_same_sign_consistent':bool(good)})
    if before!={str(p):_file_sha256(p) for p in sources}:raise RuntimeError('source changed during execution')
    output={'analysis':'exploratory_dual_marker_log_linear_sensitivity','joint_family_size':len(keys)*2,
        'primary_reproduced_exactly':True,'linear_result':linear_result,'joint_comparisons':comparisons,
        'two_scale_consistent_counts':counts,'scope':json.loads(plan.read_text())['reporting'],
        'provenance':{'completed_utc':datetime.now(timezone.utc).isoformat(),'source_sha256':before,
            'input_path':str(args.adata),'input_sha256':_file_sha256(args.adata),
            'git_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
            'git_dirty':bool(subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()),
            'versions':primary['provenance']['versions']}}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with args.out.open('x') as stream:json.dump(output,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps({'out':str(args.out),'joint_family_size':len(keys)*2,'two_scale_consistent_counts':counts}),flush=True)


if __name__=='__main__':main()
