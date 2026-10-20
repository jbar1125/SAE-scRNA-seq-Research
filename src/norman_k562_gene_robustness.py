#!/usr/bin/env python3
"""Individual RNA and single-marker omissions for the fixed historical program."""
import argparse
import copy
import csv
from datetime import datetime, timezone
import gzip
import importlib.metadata
import json
from pathlib import Path
import subprocess

import numpy as np
from statsmodels.stats.multitest import multipletests

from marker_program_baseline import normalize_counts_to_markers, marker_program_baseline, _compare_up, split_controls
from norman_direct_pca import sha

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--adata', type=Path, required=True)
    p.add_argument('--out-dir', type=Path, required=True)
    args = p.parse_args()
    if args.out_dir.exists():
        p.error('choose a fresh output directory')
    plan_path = ROOT / 'config/norman_k562_gene_robustness_plan.json'
    parent_path = ROOT / 'data_norman_continuation/k562_context/results.json.gz'
    sources = [Path(__file__).resolve(), plan_path, parent_path, ROOT / 'src/marker_program_baseline.py',
               ROOT / 'src/norman_direct_pca.py']
    before = {str(v.relative_to(ROOT)): sha(v) for v in sources}
    with gzip.open(parent_path, 'rt') as f:
        parent = json.load(f)
    if sha(args.adata) != parent['provenance']['input_sha256']:
        raise ValueError('exact original count source required')
    primary = parent['direct_fit']
    name = parent['context_program']
    markers = primary['programs'][name]['eligible_markers']
    if len(markers) != 8 or markers != parent['context_markers']:
        raise ValueError('all eight fixed markers must retain original eligibility/order')
    import anndata as ad
    a = ad.read_h5ad(args.adata)
    panels = {n: r['markers_requested'] for n, r in primary['programs'].items()}
    labels = a.obs['perturbation'].astype(str).to_numpy()
    strata = a.obs['gemgroup'].astype(str).to_numpy()
    ids = a.obs_names.astype(str).to_numpy()
    x, genes, normalization = normalize_counts_to_markers(a.X, a.var_names.to_numpy(), panels)
    del a
    check = marker_program_baseline(x, genes, labels, panels, strata=strata, cell_ids=ids, n_identity_shuffle=0)
    if check['programs'] != primary['programs'] or check['control_split'] != primary['control_split']:
        raise ValueError('all original primary direct records did not reproduce exactly')
    _, inference = split_controls(labels, strata=strata, seed=0)
    gidx = {g: i for i, g in enumerate(genes)}
    pars = primary['discovery_gene_parameters']
    columns = [gidx[g] for g in markers]
    z = (x[:, columns] - np.array([pars[g]['mean'] for g in markers])) / np.array([pars[g]['scale'] for g in markers])
    new = []
    for i, gene in enumerate(markers):
        scores = {'individual_transcript': x[:, columns[i]],
                  'leave_one_out': z[:, [j for j in range(8) if j != i]].mean(axis=1)}
        for kind, score in scores.items():
            for target in primary['perturbations']:
                values = _compare_up(score, np.flatnonzero(labels == target), inference, strata, 20)
                new.append({'kind': kind, 'gene': gene, 'pert': target, **values})
    base = copy.deepcopy(parent['representations'])
    old_rows = [row for v in base.values() for panel in v['programs'].values() for row in panel['per_perturbation']]
    if len(old_rows) != 8400 or len(new) != 1680:
        raise ValueError('planned original and new families differ')
    all_rows = old_rows + new
    q = multipletests([row['mw_p'] for row in all_rows], method='fdr_bh')[1]
    for row, value in zip(all_rows, q):
        row['gene_robustness_bh_q'] = float(value)
        row['gene_robustness_called'] = bool(value < .05 and row['auc'] >= .55 and row['status'] == 'tested')
    summary = []
    for target in primary['perturbations']:
        gene_rows = [r for r in new if r['pert'] == target and r['kind'] == 'individual_transcript']
        omission = [r for r in new if r['pert'] == target and r['kind'] == 'leave_one_out']
        if len(gene_rows) != 8 or len(omission) != 8:
            raise ValueError('missing gene/omission diagnostics')
        primary_row = next(r for r in base['direct_markers']['programs'][name]['per_perturbation'] if r['pert'] == target)
        summary.append({'pert': target, 'updated_primary_call': primary_row['gene_robustness_called'],
            'n_individual_transcript_increases': sum(r['gene_robustness_called'] for r in gene_rows),
            'n_passing_omissions': sum(r['gene_robustness_called'] for r in omission),
            'all_omissions_pass': all(r['gene_robustness_called'] for r in omission),
            'minimum_individual_gene_auc': min(r['auc'] for r in gene_rows),
            'minimum_omission_auc': min(r['auc'] for r in omission),
            'individual_auc': {r['gene']: r['auc'] for r in gene_rows},
            'omission_auc': {r['gene']: r['auc'] for r in omission}})
    if before != {str(v.relative_to(ROOT)): sha(v) for v in sources}:
        raise RuntimeError('source changed during execution')
    output = {'analysis': 'exploratory_fixed_signature_gene_and_omission_diagnostic',
        'family_size': len(all_rows), 'new_comparisons': new, 'original_models_with_updated_correction': base,
        'summary': summary, 'normalization': normalization,
        'scope': json.loads(plan_path.read_text())['limits'],
        'provenance': {'completed_utc': datetime.now(timezone.utc).isoformat(), 'source_sha256': before,
            'input_sha256': parent['provenance']['input_sha256'], 'primary_reproduced_exactly': True,
            'git_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
            'versions': {n: importlib.metadata.version(n) for n in ['numpy', 'scipy', 'anndata', 'statsmodels']}}}
    args.out_dir.mkdir(parents=True)
    with (args.out_dir / 'results.json').open('x') as f:
        json.dump(output, f, indent=2, allow_nan=False)
        f.write('\n')
    with (args.out_dir / 'comparisons.csv').open('x', newline='') as f:
        fields = ['source', 'name', 'pert', 'status', 'auc', 'mw_p', 'gene_robustness_bh_q', 'gene_robustness_called']
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for arm, v in base.items():
            for program, panel in v['programs'].items():
                for row in panel['per_perturbation']:
                    writer.writerow({'source': arm, 'name': program, **{key: row[key] for key in fields[2:]}})
        for row in new:
            writer.writerow({'source': row['kind'], 'name': row['gene'], **{key: row[key] for key in fields[2:]}})
    print(json.dumps({'family_size': len(all_rows), 'updated_primary_calls': sum(r['updated_primary_call'] for r in summary),
        'primary_and_all_omissions': sum(r['updated_primary_call'] and r['all_omissions_pass'] for r in summary),
        'preidentified_cases': [r for r in summary if r['pert'] in ['CBL', 'UBASH3B', 'PTPN12', 'KLF1']]}), flush=True)


if __name__ == '__main__':
    main()
