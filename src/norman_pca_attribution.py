#!/usr/bin/env python3
"""Attribute projected PCA score means, with explicit randomized-SVD fidelity checks."""
import argparse
import csv
from datetime import datetime, timezone
import gzip
import io
import importlib.metadata
import json
from pathlib import Path
import subprocess

import numpy as np

from marker_program_baseline import _compare_up, split_controls
from norman_direct_pca import sha

ROOT = Path(__file__).resolve().parents[1]


def contributions(expression, weights, pert_rows, control_rows, strata):
    """Exact weighted mean attribution; this is not a rank/AUC decomposition."""
    total_pairs = 0
    effect = np.zeros(len(weights), dtype=np.float64)
    for block in np.unique(strata[pert_rows]):
        left = pert_rows[strata[pert_rows] == block]
        right = control_rows[strata[control_rows] == block]
        if len(left) < 2 or len(right) < 2:
            continue
        pairs = len(left) * len(right)
        difference = (expression[left].mean(axis=0, dtype=np.float64)
                      - expression[right].mean(axis=0, dtype=np.float64))
        effect += pairs * weights * difference
        total_pairs += pairs
    if not total_pairs:
        raise ValueError('no usable within-stratum comparisons')
    return effect / total_pairs


def weighted_score_difference(scores, pert_rows, control_rows, strata):
    total = 0.
    pairs_total = 0
    for block in np.unique(strata[pert_rows]):
        left = pert_rows[strata[pert_rows] == block]
        right = control_rows[strata[control_rows] == block]
        if len(left) < 2 or len(right) < 2:
            continue
        pairs = len(left) * len(right)
        total += pairs * (scores[left].mean() - scores[right].mean())
        pairs_total += pairs
    return float(total / pairs_total)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--adata', type=Path, required=True)
    p.add_argument('--common-run', type=Path, required=True)
    p.add_argument('--out-dir', type=Path, required=True)
    args = p.parse_args()
    if args.out_dir.exists():
        p.error('choose a fresh output directory')
    common_path = args.common_run / 'results.json'
    common = json.loads(common_path.read_text())
    if not common['provenance']['model_artifacts']['direct_pca_seed0']['archived_association_inference_reproduced']:
        raise ValueError('original PCA inference was not reproduced')
    plan_path = ROOT / 'config/norman_pca_attribution_plan.json'
    plan = json.loads(plan_path.read_text())
    original_path = ROOT / 'data_norman_continuation/direct_pca/seed0.json'
    phase_path = ROOT / 'data_norman_continuation/phase_diagnostic/comparisons.csv'
    gene_path = ROOT / 'config/cell_cycle/regev_lab_cell_cycle_genes.txt'
    sources = [Path(__file__).resolve(), plan_path, original_path, phase_path, gene_path,
               ROOT / 'config/human_markers.json'] + [ROOT / 'src' / n for n in
               ['causal_pipeline.py', 'causal_grounding.py', 'marker_program_baseline.py', 'norman_direct_pca.py']]
    before = {str(v.relative_to(ROOT)): sha(v) for v in sources}
    if sha(phase_path) != common['provenance']['source_sha256'][str(phase_path.relative_to(ROOT))]:
        raise ValueError('phase archive changed')
    if sha(args.adata) != common['provenance']['input_sha256']:
        raise ValueError('input SHA differs')
    old = json.loads(original_path.read_text())
    record = old['programs']['Cycling']
    if (record['component_index'], record['orientation']) != (3, 'positive'):
        raise ValueError('original component selection differs')
    model_path = Path(old['provenance']['model_path'])
    if sha(model_path) != old['provenance']['model_sha256']:
        raise ValueError('original model SHA differs')
    original_model = np.load(model_path, allow_pickle=False)
    weights = original_model['components'][3].astype(np.float64)
    mean = original_model['mean'].astype(np.float64)
    recovered = common['provenance']['model_artifacts']['direct_pca_seed0']
    score_path = Path(recovered['reproduced_score_path'])
    if sha(score_path) != recovered['reproduced_score_sha256']:
        raise ValueError('recovered score SHA differs')
    original_scores = np.load(score_path, mmap_mode='r', allow_pickle=False)[:, 3].astype(np.float64)
    from causal_pipeline import load_perturbseq, _topk_input_metadata
    import anndata as ad
    panels = json.loads((ROOT / 'config/human_markers.json').read_text())['marker_sets']
    x, genes, labels, _ = load_perturbseq(args.adata, 'perturbation', 'control', n_hvg=4000,
        keep_genes={g for panel in panels.values() for g in panel})
    if _topk_input_metadata(x) != common['provenance']['preprocessed_input']:
        raise ValueError('preprocessed expression bytes differ')
    a = ad.read_h5ad(args.adata, backed='r')
    try:
        strata = a.obs['gemgroup'].astype(str).to_numpy()
        ids = a.obs_names.astype(str).to_numpy()
    finally:
        a.file.close()
    _, control_rows = split_controls(labels, strata=strata, seed=0)
    if ids[control_rows].tolist() != old['control_split']['inference_cell_ids']:
        raise ValueError('inference controls differ')
    projected = np.empty(len(x), dtype=np.float64)
    offset = mean @ weights
    for start in range(0, len(x), 2048):
        projected[start:start + 2048] = x[start:start + 2048].astype(np.float64) @ weights - offset
    normalized_rmse = float(np.sqrt(np.mean((projected - original_scores)**2)) / original_scores.std())
    with phase_path.open() as stream:
        phase = {r['pert']: r for r in csv.DictReader(stream)}
    marker_names = set(gene_path.read_text().splitlines())
    partition = [('phase' if gene in marker_names else 'other') +
                 ('_positive' if w > 0 else '_negative' if w < 0 else '_zero') for gene, w in zip(genes, weights)]
    rows, all_contributions = [], []
    corrected = {r['pert']: r for r in common['representations']['direct_pca_seed0']['programs']['Cycling']['per_perturbation']}
    original = {r['pert']: r for r in record['per_perturbation']}
    for target in sorted(original):
        indices = np.flatnonzero(labels == target)
        c = contributions(x, weights, indices, control_rows, strata)
        estimated = weighted_score_difference(projected, indices, control_rows, strata)
        if not np.isclose(c.sum(), estimated, atol=1e-8, rtol=1e-8):
            raise ValueError('gene contributions fail to recover projected score mean')
        auc = _compare_up(projected, indices, control_rows, strata, 20)['auc']
        totals = {name: float(c[np.array(partition) == name].sum())
                  for name in [group + '_' + sign for group in ['phase', 'other']
                               for sign in ['positive', 'negative', 'zero']]}
        row = {'pert': target, 'original_auc': original[target]['auc'], 'projected_auc': auc,
               'original_up_call': original[target]['modulates'],
               'common_joint_up_call': corrected[target]['joint_modulates'],
               'S_auc': float(phase[target]['S_auc']), 'G2M_auc': float(phase[target]['G2M_auc']),
               'projected_mean_difference': estimated,
               'original_mean_difference': weighted_score_difference(original_scores, indices, control_rows, strata),
               'contribution_sum_error': float(c.sum() - estimated),
               'negative_weight_contribution': float(c[weights < 0].sum()), 'partitions': totals}
        rows.append(row)
        all_contributions.append((target, c))
    maximum_auc_difference = max(abs(r['projected_auc'] - r['original_auc']) for r in rows)
    passed = normalized_rmse <= .01 and maximum_auc_difference <= .005
    discordant = [r for r in rows if r['original_auc'] >= .55 and r['S_auc'] <= .45 and r['G2M_auc'] <= .45]
    summary = {'n_perturbations': len(rows), 'n_input_genes': len(genes),
        'requested_phase_genes': len(marker_names), 'phase_genes_in_pca_input': sum(g in marker_names for g in genes),
        'normalized_score_rmse': normalized_rmse, 'maximum_absolute_auc_difference': maximum_auc_difference,
        'original_score_fidelity_gate_passed': passed,
        'interpretation_allowed': 'Projected mean attribution can describe the original score approximately within the declared fidelity bounds.' if passed else
                                  'Different projection: do not explain the original score using these contributions.',
        'direction_discordant_cases': len(discordant),
        'discordant_positive_projected_mean': sum(r['projected_mean_difference'] > 0 for r in discordant),
        'discordant_positive_mean_with_positive_negative_weight_sum': sum(r['projected_mean_difference'] > 0 and r['negative_weight_contribution'] > 0 for r in discordant),
        'all_105_positive_projected_mean': sum(r['projected_mean_difference'] > 0 for r in rows),
        'all_105_positive_negative_weight_sum': sum(r['negative_weight_contribution'] > 0 for r in rows),
        'maximum_absolute_contribution_sum_error': max(abs(r['contribution_sum_error']) for r in rows)}
    args.out_dir.mkdir(parents=True)
    with (args.out_dir / 'gene_contributions.csv.gz').open('xb') as raw_stream:
        with gzip.GzipFile(fileobj=raw_stream, mode='wb', mtime=0) as compressed:
            with io.TextIOWrapper(compressed, newline='') as text:
                writer = csv.writer(text)
                writer.writerow(['pert', 'gene', 'component_weight', 'partition', 'weighted_mean_contribution'])
                for target, c in all_contributions:
                    writer.writerows(zip([target]*len(genes), genes, weights, partition, c))
    if before != {str(v.relative_to(ROOT)): sha(v) for v in sources}:
        raise RuntimeError('source changed during execution')
    result = {'analysis': 'exploratory_pca_mean_attribution', 'summary': summary, 'comparisons': rows,
        'scope': plan['limits'], 'provenance': {'completed_utc': datetime.now(timezone.utc).isoformat(),
            'source_sha256': before, 'input_sha256': common['provenance']['input_sha256'],
            'common_result_sha256': sha(common_path), 'original_model_sha256': sha(model_path),
            'recovered_scores_sha256': sha(score_path),
            'gene_contributions_sha256': sha(args.out_dir / 'gene_contributions.csv.gz'),
            'versions': {name: importlib.metadata.version(name) for name in
                         ['numpy', 'scipy', 'scanpy', 'anndata', 'scikit-learn']},
            'git_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()}}
    with (args.out_dir / 'results.json').open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    with (args.out_dir / 'comparisons.csv').open('x', newline='') as stream:
        columns = [key for key in rows[0] if key != 'partitions']
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows({key: row[key] for key in columns} for row in rows)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
