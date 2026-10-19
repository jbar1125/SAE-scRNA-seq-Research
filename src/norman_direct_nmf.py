#!/usr/bin/env python3
"""Fixed direct NMF continuation; retain warnings and the actual joint-fit factors."""
import argparse
from datetime import datetime,timezone
import importlib.metadata
import json
from pathlib import Path
import subprocess
import time
import warnings

import numpy as np
from sklearn.decomposition import NMF
from sklearn.exceptions import ConvergenceWarning

from causal_pipeline import load_perturbseq
from norman_direct_pca import analyze,sha

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--adata',required=True,type=Path)
    parser.add_argument('--out',required=True,type=Path)
    parser.add_argument('--model-out',required=True,type=Path)
    parser.add_argument('--scores-out',required=True,type=Path)
    args=parser.parse_args()
    if any(p.exists() for p in [args.out,args.model_out,args.scores_out]):parser.error('choose fresh output paths')
    plan=ROOT/'config/norman_direct_nmf_plan.json';declared=json.loads(plan.read_text());fit=declared['fit']
    sources=[Path(__file__).resolve(),plan,ROOT/'config/human_markers.json']+[ROOT/'src'/n for n in
        ['norman_direct_pca.py','causal_pipeline.py','causal_grounding.py','marker_program_baseline.py','program_modulation.py','program_modulation_specificity.py']]
    before={str(p.relative_to(ROOT)):sha(p) for p in sources}
    panels=json.loads((ROOT/'config/human_markers.json').read_text())['marker_sets']
    x,genes,labels,_=load_perturbseq(args.adata,'perturbation','control',n_hvg=4000,
        keep_genes={g for panel in panels.values() for g in panel})
    if x.shape!=(111445,4081) or np.any(x<0):raise ValueError('declared nonnegative full-size input required')
    import anndata as ad
    backed=ad.read_h5ad(args.adata,backed='r')
    try:
        if backed.obs['gemgroup'].isna().any():raise ValueError('missing gemgroup')
        strata=backed.obs['gemgroup'].astype(str).to_numpy();ids=backed.obs_names.astype(str).to_numpy()
    finally:backed.file.close()
    settings={k:fit[k] for k in ['n_components','init','solver','beta_loss','tol','max_iter','random_state']}
    print(json.dumps({'stage':'loaded','shape':list(x.shape),'settings':settings}),flush=True)
    model=NMF(**settings);started=time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        scores=model.fit_transform(x)
    fit_seconds=time.perf_counter()-started
    notices=[{'category':w.category.__name__,'message':str(w.message)} for w in caught]
    converged=not any(issubclass(w.category,ConvergenceWarning) for w in caught)
    print(json.dumps({'stage':'fit_complete','iterations':int(model.n_iter_),'converged':converged,
        'elapsed_seconds':fit_seconds,'warnings':notices}),flush=True)
    for p in [args.model_out,args.scores_out]:p.parent.mkdir(parents=True,exist_ok=True)
    with args.model_out.open('xb') as stream:np.savez_compressed(stream,components=model.components_)
    with args.scores_out.open('xb') as stream:np.save(stream,scores,allow_pickle=False)
    result=analyze(scores,x,genes,labels,panels,strata)
    result['analysis']='direct_nmf_heldout_control_baseline'
    result['fit']={'settings':settings,'iterations':int(model.n_iter_),'convergence_warning_absent':converged,
        'warnings':notices,'reconstruction_frobenius_error':float(model.reconstruction_err_),
        'relative_reconstruction_error':float(model.reconstruction_err_/np.linalg.norm(x)),
        'fit_elapsed_seconds':fit_seconds,
        'orientation_note':'Nonnegative NMF scores; negative orientations have no positive discovery controls and are ineligible.'}
    for which in ['discovery','inference']:
        result['control_split'][which+'_cell_ids']=ids[result['control_split'].pop(which+'_rows')].tolist()
    if before!={str(p.relative_to(ROOT)):sha(p) for p in sources}:raise RuntimeError('source changed during execution')
    result['provenance']={'completed_utc':datetime.now(timezone.utc).isoformat(),'input_path':str(args.adata),
        'input_sha256':sha(args.adata),'model_path':str(args.model_out),'model_sha256':sha(args.model_out),
        'actual_fitted_scores_path':str(args.scores_out),'scores_sha256':sha(args.scores_out),'source_sha256':before,
        'git_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'git_dirty':bool(subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()),
        'versions':{n:importlib.metadata.version(n) for n in ['numpy','scipy','scikit-learn','scanpy','anndata','statsmodels']}}
    result['scope']=declared['limits']
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with args.out.open('x') as stream:json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps({'stage':'complete','out':str(args.out),'counts':{n:len(r['modulators']) if not r['low_power'] else None for n,r in result['programs'].items()}}),flush=True)


if __name__=='__main__':main()
