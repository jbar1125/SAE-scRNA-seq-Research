#!/usr/bin/env python3
"""Fixed alternative library denominator; do not infer absolute RNA abundance."""
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

from marker_program_baseline import normalize_counts_to_markers, split_controls, _compare_up
from norman_direct_pca import sha

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--adata', type=Path, required=True)
    p.add_argument('--out-dir', type=Path, required=True)
    args = p.parse_args()
    if args.out_dir.exists():
        p.error('choose a fresh output directory')
    plan_path = ROOT / 'config/norman_k562_composition_plan.json'
    context_path = ROOT / 'data_norman_continuation/k562_context/results.json.gz'
    gene_path = ROOT / 'data_norman_continuation/k562_gene_robustness/results.json.gz'
    sources = [Path(__file__).resolve(), plan_path, context_path, gene_path,
               ROOT / 'src/marker_program_baseline.py', ROOT / 'src/norman_direct_pca.py']
    before = {str(v.relative_to(ROOT)): sha(v) for v in sources}
    with gzip.open(context_path, 'rt') as f:
        context = json.load(f)
    with gzip.open(gene_path, 'rt') as f:
        genes_result = json.load(f)
    if sha(args.adata) != context['provenance']['input_sha256']:
        raise ValueError('original unique raw-count source required')
    primary = context['direct_fit']
    name, markers = context['context_program'], context['context_markers']
    import anndata as ad
    a = ad.read_h5ad(args.adata)
    labels = a.obs['perturbation'].astype(str).to_numpy()
    strata = a.obs['gemgroup'].astype(str).to_numpy()
    ids = a.obs_names.astype(str).to_numpy()
    all_genes = list(map(str, a.var_names))
    counts = np.asarray(a.X[:, [all_genes.index(g) for g in markers]].toarray(), dtype=np.float64)
    totals = np.asarray(a.X.sum(axis=1, dtype=np.float64)).ravel()
    outside = totals - counts.sum(axis=1)
    if np.any(outside <= 0) or np.any(totals <= 0) or not np.isfinite(outside).all():
        raise ValueError('outside-signature normalization undefined; no cells will be silently dropped')
    panels = {n: r['markers_requested'] for n, r in primary['programs'].items()}
    old_x, old_genes, _ = normalize_counts_to_markers(a.X, a.var_names.to_numpy(), panels)
    del a
    _, inference = split_controls(labels, strata=strata, seed=0)
    if ids[inference].tolist() != primary['control_split']['inference_cell_ids']:
        raise ValueError('original inference controls differ')
    pars = primary['discovery_gene_parameters']
    old_z = (old_x - np.array([pars[g]['mean'] for g in old_genes])) / np.array([pars[g]['scale'] for g in old_genes])
    indices = [list(old_genes).index(g) for g in markers]
    old_program = old_z[:, indices].mean(axis=1)
    new_x = np.log1p(counts / (outside / 1e4)[:, None])
    new_z = np.array((new_x - np.array([pars[g]['mean'] for g in markers])) /
                     np.array([pars[g]['scale'] for g in markers]), order='F')
    new_program = new_z.mean(axis=1)
    old_gene_rows = {(r['gene'], r['pert']): r for r in genes_result['new_comparisons'] if r['kind'] == 'individual_transcript'}
    old_program_rows = {r['pert']: r for r in primary['programs'][name]['per_perturbation']}
    new, recovered = [], 0
    for target in primary['perturbations']:
        rows = np.flatnonzero(labels == target)
        for i, gene in enumerate(markers + [name]):
            is_program = i == 8
            original = _compare_up(old_program if is_program else old_x[:, indices[i]], rows, inference, strata, 20)
            expected = old_program_rows[target] if is_program else old_gene_rows[gene, target]
            if (original['status'] != expected['status'] or abs(original['auc'] - expected['auc']) > 1e-7
                    or abs(original['mw_p'] - expected['mw_p']) > 1e-10):
                raise ValueError('original normalization comparison failed declared recovery bounds')
            recovered += 1
            values = _compare_up(new_program if is_program else new_x[:, i], rows, inference, strata, 20)
            new.append({'origin': 'outside_signature_program' if is_program else 'outside_signature_transcript',
                        'label': gene, 'pert': target, **values})
    old = []
    for arm, v in genes_result['original_models_with_updated_correction'].items():
        for program, panel in v['programs'].items():
            old.extend({'origin': arm, 'label': program, **copy.deepcopy(row)} for row in panel['per_perturbation'])
    old.extend({'origin': row['kind'], 'label': row['gene'], **copy.deepcopy(row)} for row in genes_result['new_comparisons'])
    all_rows = old + new
    if (len(old), len(new), recovered, len(all_rows)) != (10080, 945, 945, 11025):
        raise ValueError('complete declared correction families differ')
    q = multipletests([r['mw_p'] for r in all_rows], method='fdr_bh')[1]
    for row, value in zip(all_rows, q):
        row['normalization_bh_q'] = float(value)
        row['normalization_called'] = bool(value < .05 and row['auc'] >= .55 and row['status'] == 'tested')
    fraction = counts.sum(axis=1) / totals
    summaries = []
    for target in primary['perturbations']:
        old_genes_rows = [r for r in old if r['origin'] == 'individual_transcript' and r['pert'] == target]
        new_genes_rows = [r for r in new if r['origin'] == 'outside_signature_transcript' and r['pert'] == target]
        new_p = next(r for r in new if r['origin'] == 'outside_signature_program' and r['pert'] == target)
        old_p = next(r for r in old if r['origin'] == 'direct_markers' and r['label'] == name and r['pert'] == target)
        group = np.flatnonzero(labels == target)
        summaries.append({'pert': target, 'original_program_auc': old_p['auc'], 'outside_program_auc': new_p['auc'],
            'original_program_called': old_p['normalization_called'], 'outside_program_called': new_p['normalization_called'],
            'original_individual_calls': sum(r['normalization_called'] for r in old_genes_rows),
            'outside_individual_calls': sum(r['normalization_called'] for r in new_genes_rows),
            'original_gene_auc': {r['label']: r['auc'] for r in old_genes_rows},
            'outside_gene_auc': {r['label']: r['auc'] for r in new_genes_rows},
            'signature_fraction_quartiles': np.quantile(fraction[group], [.25, .5, .75]).tolist()})
    if before != {str(v.relative_to(ROOT)): sha(v) for v in sources}:
        raise RuntimeError('source changed during execution')
    output = {'analysis': 'exploratory_fixed_signature_denominator_sensitivity', 'family_size': len(all_rows),
        'comparisons': all_rows, 'summary': summaries,
        'inference_control_signature_fraction_quartiles': np.quantile(fraction[inference], [.25, .5, .75]).tolist(),
        'scope': json.loads(plan_path.read_text())['limits'], 'provenance': {
            'completed_utc': datetime.now(timezone.utc).isoformat(), 'source_sha256': before,
            'input_sha256': context['provenance']['input_sha256'], 'original_comparisons_recovered': recovered,
            'frozen_discovery_parameters_used': True,
            'git_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
            'versions': {n: importlib.metadata.version(n) for n in ['numpy', 'scipy', 'scanpy', 'anndata', 'statsmodels']}}}
    args.out_dir.mkdir(parents=True)
    with (args.out_dir / 'results.json').open('x') as f:
        json.dump(output, f, indent=2, allow_nan=False)
        f.write('\n')
    with (args.out_dir / 'comparisons.csv').open('x', newline='') as f:
        keys = ['origin', 'label', 'pert', 'status', 'auc', 'mw_p', 'normalization_bh_q', 'normalization_called']
        writer = csv.DictWriter(f, fieldnames=keys, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(all_rows)
    print(json.dumps({'family_size': len(all_rows), 'original_comparisons_recovered': recovered,
        'inference_control_fraction_quartiles': output['inference_control_signature_fraction_quartiles'],
        'preidentified_cases': [r for r in summaries if r['pert'] in ['CBL', 'UBASH3B', 'PTPN12', 'KLF1']]}), flush=True)


if __name__ == '__main__':
    main()
