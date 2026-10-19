#!/usr/bin/env python3
"""External phase-score directions compared with archived control-annotated PCA effects."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import subprocess

import numpy as np
from statsmodels.stats.multitest import multipletests

from marker_program_baseline import normalize_counts_to_markers,marker_program_baseline,_compare_up,split_controls,_file_sha256

ROOT=Path(__file__).resolve().parents[1]


def correct_joint(phase_rows,pca):
    """Keep every PCA panel and both phase directions in each correction family."""
    pca_rows=[r for p in pca['programs'].values() for r in p['per_perturbation']]
    for which in ['primary','target_omission']:
        values=[r['mw_p'] for r in pca_rows]+[r[which]['mw_p'] for r in phase_rows]
        q=multipletests(values,method='holm')[1]
        for row,value in zip(pca_rows,q[:len(pca_rows)]):
            row['joint_holm_'+which]=float(value)
        for row,value in zip(phase_rows,q[len(pca_rows):]):
            r=row[which];r['joint_holm_p']=float(value)
            r['flagged']=bool(r['status']=='tested' and value<.05 and r['auc']>=.55)
    return len(values)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--adata',required=True,type=Path)
    parser.add_argument('--out',required=True,type=Path)
    args=parser.parse_args()
    if args.out.exists():parser.error('choose a fresh output path')
    plan=ROOT/'config/norman_phase_diagnostic_plan.json';gene_file=ROOT/'config/cell_cycle/regev_lab_cell_cycle_genes.txt'
    pca_path=ROOT/'data_norman_continuation/direct_pca/seed0.json'
    old_path=ROOT/'data_norman_continuation/marker_baseline/gemgroup_seed0.json'
    sources=[Path(__file__).resolve(),plan,gene_file,pca_path,old_path,ROOT/'config/cell_cycle/source.json',ROOT/'src/marker_program_baseline.py']
    before={str(p.relative_to(ROOT)):_file_sha256(p) for p in sources}
    genes_requested=gene_file.read_text().splitlines()
    if len(genes_requested)!=97 or len(set(genes_requested))!=97:raise ValueError('declared97-gene list required')
    panels={'S':genes_requested[:43],'G2M':genes_requested[43:]}
    import anndata as ad
    a=ad.read_h5ad(args.adata)
    if a.obs['gemgroup'].isna().any():raise ValueError('missing gemgroup')
    labels=a.obs['perturbation'].astype(str).to_numpy();ids=a.obs_names.astype(str).to_numpy()
    strata=a.obs['gemgroup'].astype(str).to_numpy()
    x,genes,normalization=normalize_counts_to_markers(a.X,a.var_names.to_numpy(),panels);del a
    phase=marker_program_baseline(x,genes,labels,panels,strata=strata,cell_ids=ids,seed=0,n_identity_shuffle=0)
    original=json.loads(old_path.read_text());pca=json.loads(pca_path.read_text())
    if phase['control_split']!=original['control_split']:raise ValueError('original control split differs')
    discovery,inference=split_controls(labels,strata=strata,seed=0)
    gidx={g:i for i,g in enumerate(genes)};pars=phase['discovery_gene_parameters']
    standardized=(x-np.array([pars[g]['mean'] for g in genes]))/np.array([pars[g]['scale'] for g in genes])
    rows=[]
    for name,panel in phase['programs'].items():
        indices=[gidx[g] for g in panel['eligible_markers']]
        for record in panel['per_perturbation']:
            target=record['pert'];pert_rows=np.flatnonzero(labels==target)
            omitted=[i for i in indices if genes[i]!=target]
            for direction,sign in [('up',1),('down',-1)]:
                def compare(columns):
                    if len(columns)<3:return {'status':'insufficient_marker_coverage','auc':.5,'mw_p':1.}
                    return _compare_up(sign*standardized[:,columns].mean(axis=1),pert_rows,inference,strata,20)
                rows.append({'phase':name,'pert':target,'direction':direction,'primary':compare(indices),
                             'target_omission':compare(omitted),'target_removed':len(omitted)!=len(indices)})
    family=correct_joint(rows,pca)
    if family!=1365 or len(rows)!=420:raise ValueError('full planned comparison family differs')
    cycling={r['pert']:r for r in pca['programs']['Cycling']['per_perturbation']}
    marker={r['pert']:r for r in original['programs']['Cycling']['per_perturbation']}
    output_rows=[]
    for target in sorted(cycling):
        p=cycling[target];r={'pert':target,'pca_auc':p['auc'],'pca_original_up_flag':p['modulates'],
            'pca_joint_holm_p':p['joint_holm_primary'],'pca_joint_up_flag':bool(p['joint_holm_primary']<.05 and p['auc']>=.55),
            'broad_marker_auc':marker[target]['auc'],'phases':{}}
        for name in panels:
            chosen=[v for v in rows if v['pert']==target and v['phase']==name]
            r['phases'][name]={v['direction']:v for v in chosen}
        output_rows.append(r)
    counts={}
    for name in panels:
        counts[name]={}
        for direction in ['up','down']:
            chosen=[r for r in rows if r['phase']==name and r['direction']==direction]
            counts[name][direction]={'primary':sum(r['primary']['flagged'] for r in chosen),
                'omission_consistent':sum(r['primary']['flagged'] and r['target_omission']['flagged'] for r in chosen)}
    if before!={str(p.relative_to(ROOT)):_file_sha256(p) for p in sources}:raise RuntimeError('source changed')
    output={'analysis':'exploratory_cell_cycle_phase_score_diagnostic','family_size':family,'phase_direction_hypotheses':420,
        'control_split_reproduced':True,'phase_fit_record':phase,
        'fit_record_note':'The fitting helper retains its auxiliary upward BH210 fields; all diagnostic flags use the joint Holm1365 families below, not those fields.',
        'directional_results':rows,'comparisons':output_rows,'counts':counts,
        'all_pca_programs_with_joint_correction':pca['programs'],'normalization':normalization,
        'scope':json.loads(plan.read_text())['limits'],
        'provenance':{'completed_utc':datetime.now(timezone.utc).isoformat(),'source_sha256':before,
            'input_path':str(args.adata),'input_sha256':_file_sha256(args.adata),
            'git_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
            'git_dirty':bool(subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()),
            'versions':original['provenance']['versions']}}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with args.out.open('x') as stream:json.dump(output,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps({'out':str(args.out),'family_size':family,'counts':counts}),flush=True)


if __name__=='__main__':main()
