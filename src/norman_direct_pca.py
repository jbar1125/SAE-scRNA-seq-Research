#!/usr/bin/env python3
"""Exploratory ordinary-PCA baseline with discovery-only annotation and heldout controls.

Both signs of every component are searched and included in the annotation
correction. No SAE is trained. The fixed plan states comparison limitations.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess

import numpy as np
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests

from causal_grounding import gene_association
from marker_program_baseline import _compare_up, split_controls
from program_modulation import select_program_feature
from program_modulation_specificity import program_modulation_specificity


def analyze(scores, expression, genes, labels, panels, strata, *, seed=0,
            top_k=50, n_shuffle=999, min_active=50):
    x=np.asarray(expression);z=np.asarray(scores);labels=np.asarray(labels).astype(str)
    genes=list(map(str,genes));strata=np.asarray(strata).astype(str)
    if (z.ndim!=2 or x.ndim!=2 or len(x)!=len(z) or len(labels)!=len(z)
            or len(strata)!=len(z) or x.shape[1]!=len(genes) or len(set(genes))!=len(genes)
            or not np.isfinite(z).all() or not np.isfinite(x).all()):
        raise ValueError('finite aligned matrices and unique genes are required')
    discovery,inference=split_controls(labels,strata=strata,seed=seed)
    if min(len(discovery),len(inference))<20:raise ValueError('insufficient split controls')
    # No orientation is chosen from perturbation outcomes or inference controls.
    signed=np.column_stack([z,-z])
    mask=np.zeros(len(z),dtype=bool);mask[discovery]=True
    associations=gene_association(signed,x,mask)
    eligible=np.flatnonzero((signed[discovery]>0).sum(axis=0)>=min_active)
    gidx={g:i for i,g in enumerate(genes)}
    pidx={p:np.flatnonzero(labels==p) for p in np.unique(labels)
          if p!='control' and '+' not in p and np.sum(labels==p)>=20}
    out={'analysis':'direct_signed_pca_heldout_control_baseline',
         'n_cells':len(z),'n_genes':len(genes),'n_components':z.shape[1],
         'n_signed_orientations':signed.shape[1],'n_eligible_orientations':len(eligible),
         'n_perturbations':len(pidx),'family_size':len(pidx)*len(panels),
         'control_split':{'n_discovery':len(discovery),'n_inference':len(inference),
                          'discovery_rows':discovery.tolist(),'inference_rows':inference.tolist()},
         'multiplicity':'BH over all named panels x eligible single perturbations, no-match p=1',
         'programs':{}}
    flat=[]
    for name,markers in panels.items():
        indices=sorted({gidx[g] for g in markers if g in gidx})
        f,raw,adjusted=select_program_feature(associations,indices,len(genes),min(top_k,len(genes)),eligible)
        matched=len(indices)>=3 and f>=0 and adjusted<=.05
        rec={'feature':int(f) if matched else None,'enrich_p':raw,'enrich_p_bonf':adjusted,
             'low_power':not matched,'markers_present':[genes[i] for i in indices],
             'reason':None if matched else ('insufficient_marker_coverage' if len(indices)<3 else 'no_significant_program_match'),
             'component_index':int(f%z.shape[1]) if matched else None,
             'orientation':('positive' if f<z.shape[1] else 'negative') if matched else None,
             'top_associated_genes':([genes[i] for i in np.argsort(associations[f])[::-1][:12]] if matched else []),
             'per_perturbation':[]}
        for pert,rows in pidx.items():
            r=(_compare_up(signed[:,f],rows,inference,strata,20) if matched else
               {'auc':.5,'mw_p':1.,'n_cells':0,'n_controls':0,'status':'unmatched_program'})
            r['pert']=pert;rec['per_perturbation'].append(r);flat.append(r)
        out['programs'][name]=rec
    q=multipletests([r['mw_p'] for r in flat],method='fdr_bh')[1] if flat else []
    for r,qq in zip(flat,q):r.update(mw_q=float(qq),modulates=bool(qq<.05 and r['auc']>=.55 and r['status']=='tested'))
    for rec in out['programs'].values():
        rec['modulators']=[r['pert'] for r in rec['per_perturbation'] if r['modulates']]
    keep=labels!='control';keep[inference]=True
    out['identity_heterogeneity']=program_modulation_specificity(signed[keep],labels[keep],out['programs'],
        strata=strata[keep],n_shuffle=n_shuffle,seed=seed)
    return out


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--adata',required=True,type=Path)
    parser.add_argument('--out',required=True,type=Path)
    parser.add_argument('--model-out',required=True,type=Path)
    args=parser.parse_args()
    if args.out.exists() or args.model_out.exists():parser.error('choose fresh result/model paths')
    root=Path(__file__).resolve().parents[1]
    plan=root/'config/norman_direct_pca_plan.json'
    declared=json.loads(plan.read_text());settings=declared['fit']
    sources=[Path(__file__).resolve(),plan,root/'config/human_markers.json']+[root/'src'/n for n in
        ['causal_pipeline.py','causal_grounding.py','marker_program_baseline.py','program_modulation.py','program_modulation_specificity.py']]
    before={str(p.relative_to(root)):sha(p) for p in sources}
    from causal_pipeline import load_perturbseq
    import anndata as ad
    panels=json.loads((root/'config/human_markers.json').read_text())['marker_sets']
    x,genes,labels,_=load_perturbseq(args.adata,'perturbation','control',n_hvg=4000,
                                   keep_genes={g for panel in panels.values() for g in panel})
    backed=ad.read_h5ad(args.adata,backed='r')
    try:
        if backed.obs['gemgroup'].isna().any():raise ValueError('missing gemgroup')
        strata=backed.obs['gemgroup'].astype(str).to_numpy()
        ids=backed.obs_names.astype(str).to_numpy()
    finally:backed.file.close()
    print(json.dumps({'stage':'loaded','shape':list(x.shape),'n_components':settings['n_components']}),flush=True)
    model=PCA(n_components=settings['n_components'],svd_solver=settings['svd_solver'],
              iterated_power=settings['iterated_power'],random_state=settings['random_state'],whiten=False)
    scores=model.fit_transform(x)
    args.model_out.parent.mkdir(parents=True,exist_ok=True)
    with args.model_out.open('xb') as stream:
        np.savez_compressed(stream,components=model.components_,mean=model.mean_,
                            explained_variance_ratio=model.explained_variance_ratio_)
    print(json.dumps({'stage':'pca_fit','explained_variance_sum':float(model.explained_variance_ratio_.sum())}),flush=True)
    result=analyze(scores,x,genes,labels,panels,strata)
    result['fit']={'settings':settings,'explained_variance_ratio':model.explained_variance_ratio_.tolist(),
                   'explained_variance_sum':float(model.explained_variance_ratio_.sum())}
    for which in ['discovery','inference']:
        result['control_split'][which+'_cell_ids']=ids[result['control_split'].pop(which+'_rows')].tolist()
    after={str(p.relative_to(root)):sha(p) for p in sources}
    if before!=after:raise RuntimeError('analysis source changed during execution')
    result['provenance']={'completed_utc':datetime.now(timezone.utc).isoformat(),
        'input_path':str(args.adata),'input_sha256':sha(args.adata),'model_path':str(args.model_out),
        'model_sha256':sha(args.model_out),'source_sha256':before,
        'git_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
        'git_dirty':bool(subprocess.check_output(['git','status','--porcelain'],cwd=root,text=True).strip()),
        'versions':{n:importlib.metadata.version(n) for n in ['numpy','scipy','scikit-learn','scanpy','anndata','statsmodels']}}
    result['scope']=declared['comparison_limits']
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with args.out.open('x') as stream:json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps({'stage':'complete','out':str(args.out),'modulator_counts':{n:len(r['modulators']) for n,r in result['programs'].items()}}),flush=True)


if __name__=='__main__':main()
