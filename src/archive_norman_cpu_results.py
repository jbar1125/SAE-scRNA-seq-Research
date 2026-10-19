#!/usr/bin/env python3
"""Validate/archive completed full-size CPU seeds and publish scoped evidence commits.

The user authorized GitHub updates. Only generated evidence paths are committed;
unrelated staged changes are preserved. No result is invented for a pending seed.
"""
import argparse
from datetime import datetime,timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time

import numpy as np
from statsmodels.stats.multitest import multipletests

ROOT=Path(__file__).resolve().parents[1]
RUN=Path('/workspace/research-runs/norman-sae-primary')
ARCHIVE=ROOT/'data_norman_continuation/sae_primary'
INPUT_SHA='a09744916427b08d75dade996d059f717e4812463292a1e5a4e5450589cf2c6a'
PANELS=list(json.loads((ROOT/'config/human_markers.json').read_text())['marker_sets'])


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def validate_record(r,seed):
    s=r['sae'];p=r['provenance'];a=p['arguments']
    if not (s['epochs']==s['epochs_completed']==200 and s['latent_dim']==2048
            and s['k']==32 and s['batch']==1024 and s['trained_on']=='all'
            and s['n_train_cells']==111445 and s['d']==4081):
        raise ValueError('not a completed declared full-size seed')
    if a['seed']!=seed or p['input_sha256']!=INPUT_SHA:
        raise ValueError('seed or input identity mismatch')
    if set(r['programs'])!=set(PANELS) or set(r['identity_heterogeneity']['programs'])!=set(PANELS):
        raise ValueError('the complete nine-panel family is required')
    identity=r['identity_heterogeneity']
    if (identity['n_planned_programs']!=9 or identity['n_permutations_requested']!=999
            or not identity['stratified'] or a['strata_col']!='gemgroup'):
        raise ValueError('identity-inference plan mismatch')
    for rec in identity['programs'].values():
        if not .001<=rec['empirical_p']<=1:raise ValueError('invalid finite-permutation p-value')
        if not rec['matched'] and rec['empirical_p']!=1:
            raise ValueError('unmatched hypotheses must remain p=1')
    if any(rec['matched'] for rec in identity['programs'].values()) and identity['n_permutations']!=999:
        raise ValueError('the requested permutations did not execute')


def verify_large_artifacts(r,seed):
    import torch
    model=RUN/f'seed{seed}.pt';acts=RUN/f'seed{seed}_activations.npy'
    if sha(model)!=r['provenance']['checkpoint_sha256'] or sha(acts)!=r['provenance']['activations_sha256']:
        raise ValueError('model or activation hash mismatch')
    c=torch.load(model,map_location='cpu',weights_only=True)
    if (c['completed_epochs']!=200 or c['input']['shape']!=[111445,4081]
            or c['training_params']['latent_dim']!=2048 or c['training_params']['seed']!=seed):
        raise ValueError('checkpoint is incomplete or has the wrong input')
    a=np.load(acts,mmap_mode='r',allow_pickle=False)
    if a.shape!=(111445,2048) or a.dtype!=np.float32:
        raise ValueError('activation artifact shape or dtype mismatch')


def copy_exclusive(source,destination):
    """A killed copy must not leave a partial result at its final archive path."""
    temporary=None
    try:
        with tempfile.NamedTemporaryFile(dir=destination.parent,delete=False) as dst:
            temporary=Path(dst.name)
            with source.open('rb') as src:shutil.copyfileobj(src,dst)
            dst.flush();os.fsync(dst.fileno())
        os.link(temporary,destination)  # Atomic creation; never replace an existing file.
    finally:
        if temporary is not None:temporary.unlink(missing_ok=True)


def joint_summary(records):
    if len(records)!=3:raise ValueError('all three complete seeds are required')
    by_seed={r['provenance']['arguments']['seed']:r for r in records}
    if set(by_seed)!={0,1,2}:raise ValueError('need distinct seeds 0,1,2')
    for seed,r in by_seed.items():validate_record(r,seed)
    fingerprints=[r['provenance']['source_sha256'] for r in records]
    if any(v!=fingerprints[0] for v in fingerprints):
        raise ValueError('source versions differ across seeds; review before combining')
    pairs=[(seed,name) for seed in range(3) for name in PANELS]
    raw=[by_seed[seed]['identity_heterogeneity']['programs'][name]['empirical_p'] for seed,name in pairs]
    corrected=multipletests(raw,method='holm')[1]
    adjusted={(seed,name):float(p) for (seed,name),p in zip(pairs,corrected)}
    programs={}
    for name in PANELS:
        seen={}
        for seed in range(3):
            for gene in by_seed[seed]['programs'][name].get('modulators',[]):seen[gene]=seen.get(gene,0)+1
        programs[name]={
            'matched_seeds':sum(not by_seed[s]['programs'][name].get('low_power',True) for s in range(3)),
            'joint_holm_p_by_seed':{str(s):adjusted[(s,name)] for s in range(3)},
            'joint_identity_rejections':sum(adjusted[(s,name)]<=.05 for s in range(3)),
            'descriptive_modulator_counts_by_seed':{str(s):by_seed[s]['programs'][name].get('n_modulators') for s in range(3)},
            'descriptive_candidates_in_at_least_two_seeds':sorted(g for g,n in seen.items() if n>=2)}
    return {'complete_seeds':[0,1,2],'family_size':27,'multiplicity':'Holm over all 27 program-by-seed identity tests',
        'scope':'Conditional cell-level identity heterogeneity. Seed agreement is computational robustness, not biological replication. Legacy per-program modulator lists remain exploratory: this correction does not repair their control-selection or per-gene inference.',
        'programs':programs}


def publish(paths):
    relative=[str(p.relative_to(ROOT)) for p in paths]
    branch=subprocess.check_output(['git','branch','--show-current'],cwd=ROOT,text=True).strip()
    if branch!='work':raise RuntimeError('checkout branch changed; archive retained, automatic publication stopped')
    subprocess.run(['git','add','--',*relative],cwd=ROOT,check=True)
    changes=subprocess.run(['git','diff','--cached','--quiet','--',*relative],cwd=ROOT)
    if changes.returncode==0:return
    if changes.returncode!=1:raise RuntimeError('cannot inspect staged evidence')
    subprocess.run(['git','commit','--only','-m','Archive validated completed Norman CPU evidence','--',*relative],cwd=ROOT,check=True)
    env=os.environ.copy();env['GIT_TERMINAL_PROMPT']='0'
    result=subprocess.run(['git','push','origin','HEAD:refs/heads/claude/code-handoff-doc-2v3jt3'],cwd=ROOT,env=env)
    if result.returncode:print('Git push failed; local validated evidence retained for later publication',flush=True)


def collect(publish_git=False):
    ARCHIVE.mkdir(parents=True,exist_ok=True)
    new=[];records=[];owned=[]
    for seed in range(3):
        source=RUN/f'seed{seed}.json'
        if not source.exists() or time.time()-source.stat().st_mtime<3:continue
        r=json.loads(source.read_text());validate_record(r,seed)
        destination=ARCHIVE/source.name
        if destination.exists():
            if sha(destination)!=sha(source):raise ValueError('archived result differs; refusing to overwrite')
        else:
            verify_large_artifacts(r,seed)
            copy_exclusive(source,destination);new.append(destination)
        log=RUN/f'seed{seed}.log';destlog=ARCHIVE/log.name
        if destlog.exists():
            if sha(destlog)!=sha(log):raise ValueError('archived log differs; refusing to overwrite')
        else:copy_exclusive(log,destlog);new.append(destlog)
        owned.extend([destination,destlog])
        records.append(r)
    if len(records)==3:
        summary=joint_summary(records)
        target=ARCHIVE/'joint_summary.json'
        payload=json.dumps(summary,indent=2)+'\n'
        if target.exists() and target.read_text()!=payload:
            raise ValueError('existing joint summary differs; preserving it for review')
        if not target.exists():target.write_text(payload);new.append(target)
        owned.append(target)
        report=ROOT/'docs/NORMAN_SAE_RESULTS.md'
        marker='<!-- Generated by archive_norman_cpu_results.py; edit the source template, not this generated report. -->'
        lines=[marker,'# Completed Norman CPU SAE runs','',summary['scope'],'',
            'All three declared seeds completed 200 epochs; input/model/activation hashes were checked.',
            'These results require scientific interpretation and comparison with the marker baseline before manuscript conclusions.','',
            '| Program | Matched seeds | Identity rejections after joint Holm | Exploratory candidate counts (seeds 0/1/2) |',
            '|---|---:|---:|---|']
        for name,p in summary['programs'].items():
            counts='/'.join(str(p['descriptive_modulator_counts_by_seed'][str(s)]) for s in range(3))
            lines.append(f"| {name} | {p['matched_seeds']}/3 | {p['joint_identity_rejections']}/3 | {counts} |")
        lines.extend(['','Full records: [joint summary](../data_norman_continuation/sae_primary/joint_summary.json).',
            'No new mechanism, pathway specificity, independent replication, novelty or SAE advantage is established by these tables alone.',''])
        content='\n'.join(lines)
        if report.exists() and not report.read_text().startswith(marker):
            raise ValueError('existing report is not owned by this generator; preserving it')
        if not report.exists() or report.read_text()!=content:report.write_text(content);new.append(report)
        owned.append(report)
    manifest=ARCHIVE/'manifest.json'
    if records and (new or not manifest.exists()):
        if manifest.exists() and json.loads(manifest.read_text()).get('generator')!='archive_norman_cpu_results.py':
            raise ValueError('existing manifest is not owned by this generator; preserving it')
        m={'generator':'archive_norman_cpu_results.py','updated_utc':datetime.now(timezone.utc).isoformat(),'complete_seeds':sorted(r['provenance']['arguments']['seed'] for r in records),
            'artifacts':{p.name:sha(p) for p in ARCHIVE.iterdir() if p.is_file() and p.name!='manifest.json'}}
        manifest.write_text(json.dumps(m,indent=2)+'\n');new.append(manifest)
        print(json.dumps({'new_archive_paths':[str(p.relative_to(ROOT)) for p in new],'complete_seeds':m['complete_seeds']}),flush=True)
    if records and publish_git:
        # Also recover files copied before a prior process died before committing.
        publish([*owned,manifest])
    return len(records)==3


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--watch',action='store_true')
    parser.add_argument('--publish',action='store_true',help='publish only validated generated evidence; user-authorized branch')
    args=parser.parse_args()
    while True:
        complete=collect(args.publish)
        if complete or not args.watch:break
        time.sleep(30)
