#!/usr/bin/env python3
"""Fixed-score, scale-dependent dual-minus-single marker contrasts; exploratory only."""
import argparse
from datetime import datetime,timezone
import importlib.metadata
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.stats import t
from statsmodels.stats.multitest import multipletests

from marker_program_baseline import normalize_counts_to_markers,split_controls,_file_sha256

ROOT=Path(__file__).resolve().parents[1]


def interaction(score,labels,strata,pair,controls,minimum=20):
    """Weighted within-stratum four-independent-group Welch contrast."""
    score=np.asarray(score,float);labels=np.asarray(labels).astype(str);strata=np.asarray(strata).astype(str)
    if score.ndim!=1 or len(score)!=len(labels) or len(strata)!=len(labels) or not np.isfinite(score).all():
        raise ValueError('finite aligned score, label and stratum arrays required')
    genes=pair.split('+')
    if len(genes)!=2 or len(set(genes))!=2:raise ValueError('two distinct targets required')
    groups={'dual':np.flatnonzero(labels==pair),'A':np.flatnonzero(labels==genes[0]),
            'B':np.flatnonzero(labels==genes[1]),'control':np.asarray(controls)}
    out={'status':'insufficient_comparison_cells','contrast':None,'standard_error':None,'df':None,
         'interval_95':None,'p':1.,'group_cells':{g:len(v) for g,v in groups.items()},'strata':[]}
    if min(map(len,groups.values()))<minimum:return out
    for block in np.unique(strata[groups['dual']]):
        rows={g:idx[strata[idx]==block] for g,idx in groups.items()}
        if min(map(len,rows.values()))<2:continue
        means={g:float(score[idx].mean()) for g,idx in rows.items()}
        contributions={g:float(score[idx].var(ddof=1)/len(idx)) for g,idx in rows.items()}
        out['strata'].append({'stratum':str(block),'n':{g:len(v) for g,v in rows.items()},
            'contrast':means['dual']-means['A']-means['B']+means['control'],
            'variance_contributions':contributions})
    used=sum(r['n']['dual'] for r in out['strata']);coverage=used/len(groups['dual'])
    out['dual_cell_coverage']=coverage
    if coverage<.8:out['status']='insufficient_stratum_overlap';return out
    delta=variance=denominator=0.
    for r in out['strata']:
        w=r['n']['dual']/used;r['weight']=w;delta+=w*r['contrast']
        for g,v in r['variance_contributions'].items():
            term=w*w*v;variance+=term;denominator+=term*term/(r['n'][g]-1)
    if variance<=0 or denominator<=0:out['status']='zero_estimated_variance';return out
    se=float(np.sqrt(variance));df=float(variance*variance/denominator)
    probability=float(2*t.sf(abs(delta/se),df));radius=float(t.ppf(.975,df)*se)
    out.update(status='tested',contrast=float(delta),standard_error=se,df=df,p=probability,
               interval_95=[float(delta-radius),float(delta+radius)])
    return out


def score_panel(expression,indices,fit,discovery):
    if len(indices)<3:return None
    standardized=(expression[:,indices]-fit['means'][indices])/fit['scales'][indices]
    score=standardized.mean(axis=1);sd=float(score[discovery].std(ddof=0))
    return score/sd if sd>1e-8 else None


def analyze(expression,genes,labels,strata,ids,baseline):
    genes=list(map(str,genes));gidx={g:i for i,g in enumerate(genes)}
    discovery,inference=split_controls(labels,strata=strata,seed=0)
    for which,rows in [('discovery',discovery),('inference',inference)]:
        if ids[rows].tolist()!=baseline['control_split'][which+'_cell_ids']:raise ValueError('control split differs')
    parameters=baseline['discovery_gene_parameters']
    for i,g in enumerate(genes):
        r=parameters[g]
        if (not np.isclose(expression[discovery,i].mean(),r['mean'],rtol=0,atol=1e-12)
            or not np.isclose(expression[discovery,i].std(ddof=0),r['population_sd'],rtol=0,atol=1e-12)
            or np.count_nonzero(expression[discovery,i])!=r['n_nonzero']):raise ValueError('baseline fit failed reproduction: '+g)
    fit={'means':np.array([parameters[g]['mean'] for g in genes]),
         'scales':np.array([parameters[g]['scale'] for g in genes])}
    pairs=sorted(p for p in set(labels) if '+' in p)
    out={'analysis':'conditional_cell_level_dual_marker_additivity_contrasts','pairs':pairs,
         'family_size':len(pairs)*len(baseline['programs']),'n_dual_cells':int(sum('+' in p for p in labels)),
         'score_units':'discovery-control SD of the fixed marker average, separately for target omission',
         'multiplicity':'Holm over all panels x observed dual pairs, separately for primary and target omission',
         'baseline_gene_parameters_reproduced':True,'control_split_reproduced':True,'programs':{}}
    flat=[]
    for name,panel in baseline['programs'].items():
        indices=[gidx[g] for g in panel['eligible_markers']]
        primary=score_panel(expression,indices,fit,discovery)
        rec={'eligible_markers':panel['eligible_markers'],'low_power':primary is None,'comparisons':[]}
        for pair in pairs:
            omitted=[i for i in indices if genes[i] not in pair.split('+')]
            sensitivity=score_panel(expression,omitted,fit,discovery) if len(omitted)!=len(indices) else primary
            def compare(values):
                return interaction(values,labels,strata,pair,inference) if values is not None else {
                    'status':'insufficient_marker_coverage','contrast':None,'p':1.}
            row={'pair':pair,'primary':compare(primary),'target_omission':compare(sensitivity),
                 'omitted_targets':[g for g in pair.split('+') if g in panel['eligible_markers']],
                 'markers_remaining':len(omitted)}
            rec['comparisons'].append(row);flat.append(row)
        out['programs'][name]=rec
    for which in ['primary','target_omission']:
        q=multipletests([r[which]['p'] for r in flat],method='holm')[1]
        for row,adjusted in zip(flat,q):
            r=row[which];r['holm_p']=float(adjusted)
            r['flagged']=bool(r['status']=='tested' and adjusted<.05 and abs(r['contrast'])>=.2)
    for row in flat:
        a,b=row['primary'],row['target_omission']
        row['same_sign_omission_robust']=bool(a['flagged'] and b['flagged'] and np.sign(a['contrast'])==np.sign(b['contrast']))
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--adata',required=True,type=Path)
    parser.add_argument('--out',required=True,type=Path)
    args=parser.parse_args()
    if args.out.exists():parser.error('choose a fresh output path')
    base=ROOT/'data_norman_continuation/marker_baseline/gemgroup_seed0.json'
    plan=ROOT/'config/norman_dual_marker_plan.json'
    sources=[Path(__file__).resolve(),plan,base,ROOT/'config/human_markers.json',
             ROOT/'src/marker_program_baseline.py']
    before={str(p.relative_to(ROOT)):_file_sha256(p) for p in sources}
    import anndata as ad
    a=ad.read_h5ad(args.adata)
    labels=a.obs['perturbation'].astype(str).to_numpy();ids=a.obs_names.astype(str).to_numpy()
    if a.obs['gemgroup'].isna().any():raise ValueError('missing gemgroup')
    strata=a.obs['gemgroup'].astype(str).to_numpy()
    panels=json.loads((ROOT/'config/human_markers.json').read_text())['marker_sets']
    expression,genes,normalization=normalize_counts_to_markers(a.X,a.var_names.to_numpy(),panels)
    del a
    result=analyze(expression,genes,labels,strata,ids,json.loads(base.read_text()))
    if before!={str(p.relative_to(ROOT)):_file_sha256(p) for p in sources}:raise RuntimeError('source changed')
    result['normalization']=normalization
    result['scope']=json.loads(plan.read_text())['reporting']
    result['provenance']={'completed_utc':datetime.now(timezone.utc).isoformat(),'source_sha256':before,
        'input_path':str(args.adata),'input_sha256':_file_sha256(args.adata),
        'git_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'git_dirty':bool(subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()),
        'versions':{n:importlib.metadata.version(n) for n in ['numpy','scipy','anndata','scanpy','statsmodels']}}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with args.out.open('x') as stream:json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps({'out':str(args.out),'family_size':result['family_size'],'omission_robust':{
        n:sum(r['same_sign_omission_robust'] for r in p['comparisons']) for n,p in result['programs'].items()}}),flush=True)


if __name__=='__main__':main()
