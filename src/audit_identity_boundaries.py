#!/usr/bin/env python3
"""Execute fixed synthetic severity, confounding and clustered-cell boundary worlds."""
import argparse
from datetime import datetime,timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy.stats import binomtest

import program_modulation_specificity as identity

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def worlds(seed,p):
    rng=np.random.default_rng(np.random.SeedSequence([seed,p['generator_salt'],1]))
    ng=p['n_perturbations'];nc=p['cells_per_perturbation'];nctrl=p['n_controls'];nf=p['n_programs']
    labels=np.array(['control']*nctrl+[f'P{g}' for g in range(ng) for _ in range(nc)])
    baseline=rng.normal(size=(len(labels),nf))
    generic=baseline.copy();generic[nctrl:]+=3
    severity=np.repeat(np.linspace(0,4,ng),nc)
    loadings=np.linspace(.5,1.5,nf)
    graded=baseline.copy();graded[nctrl:]+=severity[:,None]*loadings[None,:]
    # Rank-one group offsets encode one generic factor, not nine selective pathways.
    effects=np.linspace(0,4,ng)[:,None]*loadings[None,:]
    assert np.linalg.matrix_rank(effects)==1
    confounded=baseline.copy();confounded[nctrl:]+=severity[:,None]
    blocks=np.array([f'B{i%ng}' for i in range(nctrl)]+[f'B{g}' for g in range(ng) for _ in range(nc)])
    crng=np.random.default_rng(np.random.SeedSequence([seed,p['generator_salt'],2]))
    offsets=crng.normal(size=(ng,4,nf))
    cells=offsets[:,:,None,:]+crng.normal(0,.15,(ng,4,10,nf))
    controls=crng.normal(size=(nctrl,nf))
    clustered=np.vstack([controls,cells.reshape(ng*nc,nf)])
    unit_scores=np.vstack([controls,cells.mean(axis=2).reshape(ng*4,nf)])
    unit_labels=np.array(['control']*nctrl+[f'P{g}' for g in range(ng) for _ in range(4)])
    return {
        'identical_generic':(generic,labels,None,20),
        'graded_generic_severity':(graded,labels,None,20),
        'confounded_batch_pooled':(confounded,labels,None,20),
        'confounded_batch_conditioned':(confounded,labels,blocks,20),
        'clustered_null_cell_shuffle':(clustered,labels,None,20),
        'clustered_null_unit_means':(unit_scores,unit_labels,None,2)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out-dir',required=True,type=Path)
    args=parser.parse_args()
    if args.out_dir.exists():parser.error('refusing to overwrite prior experiment')
    plan=ROOT/'config/identity_boundary_audit_plan.json'
    declared=json.loads(plan.read_text());p=declared['fixed']
    sources=[Path(__file__).resolve(),Path(identity.__file__).resolve(),plan]
    before={str(s.relative_to(ROOT)):sha(s) for s in sources}
    args.out_dir.mkdir(parents=True,exist_ok=False)
    provenance={'started_utc':datetime.now(timezone.utc).isoformat(),'plan':declared,'source_sha256':before,
        'git_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'git_dirty':bool(subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()),
        'versions':{n:importlib.metadata.version(n) for n in ['numpy','scipy','statsmodels']}}
    (args.out_dir/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    matches={f'program{j}':{'feature':j,'enrich_p_bonf':1e-6} for j in range(p['n_programs'])}
    counts={};started=time.perf_counter()
    with (args.out_dir/'results.jsonl').open('x',buffering=1) as stream:
        for seed in range(p['seed_start'],p['seed_start']+p['simulation_seeds']):
            for world,(scores,labels,blocks,minimum) in worlds(seed,p).items():
                result=identity.program_modulation_specificity(scores,labels,matches,strata=blocks,
                    min_cells=minimum,n_shuffle=p['permutations'],seed=seed,alpha=p['alpha'])
                n=sum(r['conditional_heterogeneity'] for r in result['programs'].values())
                row={'seed':seed,'world':world,'any_rejection':bool(n),'n_rejections':n,
                     'n_permutations':result['n_permutations'],'n_movable_strata':result['n_movable_strata'],
                     'n_exchangeable_rows':result['n_exchangeable_cells'],'minimum_rows_per_identity':minimum,
                     'programs':result['programs']}
                if world=='confounded_batch_conditioned':
                    assert result['n_movable_strata']==0 and result['n_permutations']==0 and n==0
                    assert all(r['status']=='identity_confounded_with_strata' and r['empirical_p']==1 for r in result['programs'].values())
                stream.write(json.dumps(row,allow_nan=False)+'\n')
                tally=counts.setdefault(world,{'any_rejection_count':0,'all_nine_rejection_count':0,'total_rejected_hypotheses':0})
                tally['any_rejection_count']+=bool(n);tally['all_nine_rejection_count']+=n==9;tally['total_rejected_hypotheses']+=n
            if (seed+1)%20==0:print(json.dumps({'completed_seeds':seed+1,'elapsed_seconds':time.perf_counter()-started}),flush=True)
    if before!={str(s.relative_to(ROOT)):sha(s) for s in sources}:raise RuntimeError('source changed; partial evidence retained')
    for tally in counts.values():
        interval=binomtest(tally['any_rejection_count'],p['simulation_seeds']).proportion_ci(.95,method='exact')
        tally.update(n_simulations=p['simulation_seeds'],any_rejection_rate=tally['any_rejection_count']/p['simulation_seeds'],
                     exact_95_percent_interval=[float(interval.low),float(interval.high)])
    summary={'completed_utc':datetime.now(timezone.utc).isoformat(),'scope':declared['scope'],
             'results_sha256':sha(args.out_dir/'results.jsonl'),'source_sha256':before,
             'analyses':counts,'elapsed_seconds':time.perf_counter()-started,
             'interpretation':'Graded generic severity is an alternative to the identity null, not a false rejection under it. Confounded conditioning is uninformative, not evidence of absence. Cluster-cell permutations violate independence; unit means are a separate simulated-unit analysis. These paired worlds do not estimate independent biological replication.'}
    (args.out_dir/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary),flush=True)


if __name__=='__main__':main()
