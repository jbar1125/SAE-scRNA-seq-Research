#!/usr/bin/env python3
"""Complete direct and author-alternative target-transcript induction diagnostic."""
import argparse
import csv
from datetime import datetime, timezone
import importlib.metadata
import json
from pathlib import Path
import subprocess

import numpy as np
from statsmodels.stats.multitest import multipletests

from marker_program_baseline import normalize_counts_to_markers, split_controls, _compare_up
from norman_direct_pca import sha

ROOT = Path(__file__).resolve().parents[1]


def analyze_targets(x, genes, labels, strata, alternatives):
    _, inference = split_controls(labels, strata=strata, seed=0)
    targets = sorted(p for p in np.unique(labels) if p != 'control' and '+' not in p and np.sum(labels == p) >= 20)
    gidx = {str(g): i for i, g in enumerate(genes)}
    requests = [(p, p, 'direct_name') for p in targets]
    if not set(alternatives) <= set(targets):
        raise ValueError('alternative labels absent from intervention library')
    requests += [(p, gene, 'author_alternative') for p, gene in alternatives.items()]
    records = []
    for target, gene, kind in requests:
        rows = np.flatnonzero(labels == target)
        if gene in gidx:
            r = _compare_up(x[:, gidx[gene]], rows, inference, strata, 20)
        else:
            r = {'status': 'transcript_symbol_absent', 'auc': .5, 'mw_p': 1., 'n_cells': len(rows),
                 'n_controls': len(inference)}
        records.append({'pert': target, 'transcript': gene, 'comparison': kind, **r})
    q = multipletests([r['mw_p'] for r in records], method='fdr_bh')[1]
    for row, value in zip(records, q):
        row['joint_bh_q'] = float(value)
        row['target_rna_increase'] = bool(row['status'] == 'tested' and value < .05 and row['auc'] >= .55)
    return records


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--adata', type=Path, required=True)
    p.add_argument('--out-dir', type=Path, required=True)
    args = p.parse_args()
    if args.out_dir.exists():
        p.error('choose a fresh output directory')
    plan_path = ROOT / 'config/norman_target_rna_plan.json'
    plan = json.loads(plan_path.read_text())
    primary_path = ROOT / 'data_norman_continuation/marker_baseline/gemgroup_seed0.json'
    sources = [Path(__file__).resolve(), plan_path, primary_path,
               ROOT / 'src/marker_program_baseline.py', ROOT / 'src/norman_direct_pca.py']
    before = {str(a.relative_to(ROOT)): sha(a) for a in sources}
    if sha(args.adata) != plan['input_sha256']:
        raise ValueError('source SHA differs')
    import anndata as ad
    a = ad.read_h5ad(args.adata)
    if a.obs['gemgroup'].isna().any():
        raise ValueError('missing gemgroup')
    labels = a.obs['perturbation'].astype(str).to_numpy()
    ids = a.obs_names.astype(str).to_numpy()
    strata = a.obs['gemgroup'].astype(str).to_numpy()
    targets = sorted(p for p in np.unique(labels) if p != 'control' and '+' not in p)
    requested = set(targets) | set(plan['alternative_transcripts'].values())
    x, genes, normalization = normalize_counts_to_markers(a.X, a.var_names.to_numpy(), {'targets': sorted(requested)})
    del a
    rows = analyze_targets(x, genes, labels, strata, plan['alternative_transcripts'])
    if len(rows) != 109 or len({(r['pert'], r['comparison']) for r in rows}) != 109:
        raise ValueError('complete planned family differs')
    discovery, inference = split_controls(labels, strata=strata, seed=0)
    split = {'n_discovery': len(discovery), 'n_inference': len(inference),
             'discovery_cell_ids': ids[discovery].tolist(), 'inference_cell_ids': ids[inference].tolist()}
    old_split = json.loads(primary_path.read_text())['control_split']
    if split != old_split:
        raise ValueError('original control split differs')
    summary = {kind: {'n_planned': sum(r['comparison'] == kind for r in rows),
                     'n_present': sum(r['comparison'] == kind and r['status'] == 'tested' for r in rows),
                     'n_target_rna_increases': sum(r['comparison'] == kind and r['target_rna_increase'] for r in rows)}
               for kind in ['direct_name', 'author_alternative']}
    if before != {str(a.relative_to(ROOT)): sha(a) for a in sources}:
        raise RuntimeError('source changed during execution')
    result = {'analysis': 'exploratory_target_rna_diagnostic', 'family_size': len(rows),
        'summary': summary, 'comparisons': rows, 'normalization': normalization, 'control_split': split,
        'scope': plan['limits'], 'provenance': {'completed_utc': datetime.now(timezone.utc).isoformat(),
            'input_sha256': plan['input_sha256'], 'source_sha256': before,
            'git_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
            'versions': {n: importlib.metadata.version(n) for n in ['numpy', 'scipy', 'anndata', 'statsmodels']}}}
    args.out_dir.mkdir(parents=True)
    with (args.out_dir / 'results.json').open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    with (args.out_dir / 'comparisons.csv').open('x', newline='') as stream:
        keys = ['pert', 'transcript', 'comparison', 'status', 'auc', 'mw_p', 'joint_bh_q', 'target_rna_increase',
                'n_cells', 'n_controls']
        writer = csv.DictWriter(stream, fieldnames=keys, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
