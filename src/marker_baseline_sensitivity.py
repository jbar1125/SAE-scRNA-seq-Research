#!/usr/bin/env python3
"""Audit marker-score associations against every control reference and depth strata.

This is the explicit exploratory plan in config/marker_baseline_sensitivity_plan.json.
The primary fit/split is preserved. No favorable reference or technical group is chosen.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
from statsmodels.stats.multitest import multipletests

from marker_program_baseline import (_compare_up, _file_sha256,
    normalize_counts_to_markers, marker_program_baseline)


def analyze(expression, genes, labels, cell_ids, gemgroups, guide_ids, umi_counts, primary):
    x = np.asarray(expression, dtype=float)
    genes, labels, ids, groups, guides = [np.asarray(v).astype(str)
        for v in (genes, labels, cell_ids, gemgroups, guide_ids)]
    umi = np.asarray(umi_counts, dtype=float)
    if x.shape != (len(labels), len(genes)) or any(len(v) != len(labels) for v in (ids,groups,guides,umi)):
        raise ValueError('unaligned sensitivity inputs')
    if len(set(ids)) != len(ids) or not np.isfinite(umi).all() or (umi <= 0).any():
        raise ValueError('cell IDs must be unique; UMI counts must be finite and positive')
    by_id = {v:i for i,v in enumerate(ids)}
    discovery = np.array([by_id[v] for v in primary['control_split']['discovery_cell_ids']])
    inference = np.array([by_id[v] for v in primary['control_split']['inference_cell_ids']])
    if np.intersect1d(discovery,inference).size or np.any(labels[np.r_[discovery,inference]]!='control'):
        raise ValueError('invalid primary control split')
    pars = primary['discovery_gene_parameters']
    z = (x - np.array([pars[g]['mean'] for g in genes])) / np.array([pars[g]['scale'] for g in genes])
    gidx = {g:i for i,g in enumerate(genes)}
    depth = np.empty(len(labels),dtype=object)
    boundaries = {}
    for block in np.unique(groups):
        d = discovery[groups[discovery]==block]
        if len(d)<4:
            raise ValueError('insufficient discovery controls for depth quantiles')
        cuts = np.unique(np.quantile(np.log1p(umi[d]),[.25,.5,.75]))
        ix = np.flatnonzero(groups==block)
        bins = np.searchsorted(cuts,np.log1p(umi[ix]),side='right')
        depth[ix] = [block+'|depth'+str(b) for b in bins]
        boundaries[block] = cuts.tolist()
    references = sorted(set(guides[inference]))
    rows_by_condition = {'depth_matched':[]}
    rows_by_condition.update({'reference:'+r:[] for r in references})
    descriptive = []
    for name, program in primary['programs'].items():
        indices = [gidx[g] for g in program['eligible_markers']]
        score = z[:,indices].mean(1) if not program['low_power'] else np.zeros(len(labels))
        for pert in primary['perturbations']:
            pert_rows = np.flatnonzero(labels==pert)
            omitted = [i for i in indices if genes[i]!=pert]
            omission_score = z[:,omitted].mean(1) if len(omitted)>=3 else np.zeros(len(labels))
            for condition in rows_by_condition:
                controls, strata = inference, groups
                if condition=='depth_matched': strata=depth.astype(str)
                else: controls=inference[guides[inference]==condition.split(':',1)[1]]
                p = (_compare_up(score,pert_rows,controls,strata,20) if not program['low_power'] else
                     dict(auc=.5,mw_p=1.,n_cells=0,n_controls=0,status='low_power_marker_panel'))
                o = (_compare_up(omission_score,pert_rows,controls,strata,20) if len(omitted)>=3 else
                     dict(auc=.5,mw_p=1.,n_cells=0,n_controls=0,status='low_power_after_target_omission'))
                rows_by_condition[condition].append(dict(program=name,pert=pert,**p,target_omission=o))
            for block in np.unique(groups):
                result = (_compare_up(score,pert_rows[groups[pert_rows]==block],
                      inference[groups[inference]==block],None,20) if not program['low_power'] else
                      dict(auc=.5,mw_p=1.,n_cells=0,n_controls=0,status='low_power_marker_panel'))
                descriptive.append(dict(program=name,pert=pert,gemgroup=block,
                    **{k:result[k] for k in ('auc','n_cells','n_controls','status')}))
    def adjust(rows):
        p = multipletests([r['mw_p'] for r in rows],method='fdr_bh')[1]
        o = multipletests([r['target_omission']['mw_p'] for r in rows],method='fdr_bh')[1]
        for row,pq,oq in zip(rows,p,o):
            row.update(mw_q=float(pq),modulates=bool(pq<.05 and row['auc']>=.55 and row['status']=='tested'))
            sensitivity = row['target_omission']
            sensitivity.update(mw_q=float(oq),modulates=bool(oq<.05 and sensitivity['auc']>=.55 and sensitivity['status']=='tested'))
            row['survives_target_omission']=bool(row['modulates'] and sensitivity['modulates'])
    adjust(rows_by_condition['depth_matched'])
    adjust([r for name,rows in rows_by_condition.items() if name.startswith('reference:') for r in rows])
    def candidates(rows,name):
        return {r['pert'] for r in rows if r['program']==name and r['survives_target_omission']}
    summary={}
    for name,program in primary['programs'].items():
        by_ref={r:sorted(candidates(rows_by_condition['reference:'+r],name)) for r in references}
        stable=set.intersection(*(set(v) for v in by_ref.values())) if by_ref else set()
        original=set(program['target_omission_robust_modulators'])
        matched=candidates(rows_by_condition['depth_matched'],name)
        summary[name]={
            'primary_target_omission_robust':sorted(original),
            'depth_matched_target_omission_robust':sorted(matched),
            'control_reference_target_omission_robust':by_ref,
            'all_control_references':sorted(stable),
            'primary_and_depth_and_every_reference':sorted(original & matched & stable)}
    return {'analysis':'exploratory_marker_baseline_sensitivity','summary':summary,
        'control_references':references,'control_reference_family_size':len(references)*primary['family_size'],
        'depth_family_size':primary['family_size'],'discovery_log1p_umi_quartiles':boundaries,
        'comparisons':rows_by_condition,'per_gemgroup_auc_descriptive':descriptive,
        'scope':'Conditional cell-level sensitivity; depth may be a mediator; guide references and gemgroups are not biological replicates'}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--adata',type=Path,required=True)
    p.add_argument('--primary',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--plan',type=Path,default=Path('config/marker_baseline_sensitivity_plan.json'))
    args=p.parse_args()
    if args.out.exists(): p.error('choose a fresh output path')
    source_paths={'source':Path(__file__),'baseline_source':Path(__file__).with_name('marker_program_baseline.py'),'plan':args.plan}
    source_hashes={name:_file_sha256(path) for name,path in source_paths.items()}
    primary=json.loads(args.primary.read_text())
    if _file_sha256(args.adata)!=primary['provenance']['input']['sha256']:
        raise ValueError('input does not match primary baseline')
    import anndata as ad
    a=ad.read_h5ad(args.adata)
    panels={name:rec['markers_requested'] for name,rec in primary['programs'].items()}
    x,genes,_=normalize_counts_to_markers(a.X,a.var_names.to_numpy(),panels,duplicate_genes='sum')
    labels=a.obs['perturbation'].astype(str).to_numpy()
    groups=a.obs['gemgroup'].astype(str).to_numpy()
    ids=a.obs_names.astype(str).to_numpy()
    guides=a.obs['guide_identity'].astype(str).to_numpy()
    umi=a.obs['UMI_count'].to_numpy()
    del a
    check=marker_program_baseline(x,genes,labels,panels,strata=groups,cell_ids=ids,n_identity_shuffle=0)
    for name,record in primary['programs'].items():
        if record!=check['programs'][name]:
            raise ValueError('primary marker-score computation did not reproduce exactly')
    result=analyze(x,genes,labels,ids,groups,guides,umi,primary)
    for name,path in source_paths.items():
        if _file_sha256(path)!=source_hashes[name]: raise RuntimeError('source changed during run')
    result['provenance']={'primary_sha256':_file_sha256(args.primary),
        'input_sha256':primary['provenance']['input']['sha256'],'source_sha256':source_hashes,
        'completed_utc':datetime.now(timezone.utc).isoformat(),'primary_reproduced_exactly':True}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with args.out.open('x') as f: json.dump(result,f,indent=2,allow_nan=False)
    print(json.dumps({'summary':result['summary'],'out':str(args.out)}),flush=True)


if __name__=='__main__': main()
