#!/usr/bin/env python3
"""Harmonize inference for six completed Norman representations without SAE retraining."""
from __future__ import annotations

import argparse
import copy
import csv
from datetime import datetime, timezone
import gc
import importlib.metadata
import json
from pathlib import Path
import subprocess
import time

import numpy as np
from statsmodels.stats.multitest import multipletests

from norman_direct_pca import analyze, sha

ROOT = Path(__file__).resolve().parents[1]
ARMS = ['all_cell_sae_seed0', 'all_cell_sae_seed1', 'all_cell_sae_seed2',
        'control_only_sae_seed0', 'direct_pca_seed0', 'capped_direct_nmf_seed0']


def correct_complete_family(records):
    """Correct every arm and panel, keeping unscored hypotheses and original results."""
    if set(records) != set(ARMS):
        raise ValueError('all six declared representations are required')
    result = copy.deepcopy(records)
    panels = list(result[ARMS[0]]['programs'])
    expected = None
    tests, identity = [], []
    for arm in ARMS:
        if set(result[arm]['programs']) != set(panels):
            raise ValueError('panel family differs across representations')
        if set(result[arm]['identity_heterogeneity']['programs']) != set(panels):
            raise ValueError('identity family is incomplete')
        for name in panels:
            rows = result[arm]['programs'][name]['per_perturbation']
            names = [r['pert'] for r in rows]
            if len(set(names)) != len(names):
                raise ValueError('duplicate perturbation hypothesis')
            if expected is None:
                expected = set(names)
            if set(names) != expected:
                raise ValueError('perturbation family differs across arms or panels')
            tests.extend(rows)
            identity.append(result[arm]['identity_heterogeneity']['programs'][name])
    q = multipletests([r['mw_p'] for r in tests], method='fdr_bh')[1]
    for row, value in zip(tests, q):
        row['joint_bh_q'] = float(value)
        row['joint_modulates'] = bool(value < .05 and row['auc'] >= .55 and row['status'] == 'tested')
    adjusted = multipletests([r['empirical_p'] for r in identity], method='holm')[1]
    for row, value in zip(identity, adjusted):
        row['joint_holm_p'] = float(value)
        row['joint_heterogeneity'] = bool(value < .05 and row['status'] == 'tested')
    return result, len(tests), len(identity)


def require_reproduced(reproduced, archived):
    """Abort if regenerated PCA/NMF scores no longer reproduce archived inference."""
    if reproduced['control_split'] != archived['control_split']:
        raise ValueError('control split was not reproduced')
    for name, old in archived['programs'].items():
        new = reproduced['programs'][name]
        for key in ['feature', 'component_index', 'orientation', 'low_power', 'top_associated_genes']:
            if new[key] != old[key]:
                raise ValueError(f'{name}: {key} was not reproduced')
        pairs = {r['pert']: r for r in new['per_perturbation']}
        if len(pairs) != len(old['per_perturbation']):
            raise ValueError('comparison family differs')
        for row in old['per_perturbation']:
            r = pairs[row['pert']]
            if (r['status'] != row['status'] or abs(r['auc'] - row['auc']) > 1e-7
                    or abs(r['mw_p'] - row['mw_p']) > 1e-10):
                raise ValueError(f'{name}/{row["pert"]}: score inference was not reproduced')


def summarize(records, phase):
    summary = {'association_family_size': 5670, 'identity_family_size': 54,
               'programs': {}, 'computational_candidates_in_at_least_two_all_cell_seeds': {},
               'cycling_direction_comparisons': {}}
    panels = list(records[ARMS[0]]['programs'])
    for name in panels:
        summary['programs'][name] = {}
        for arm in ARMS:
            r = records[arm]['programs'][name]
            i = records[arm]['identity_heterogeneity']['programs'][name]
            candidates = [p['pert'] for p in r['per_perturbation'] if p['joint_modulates']]
            summary['programs'][name][arm] = {'matched': not r['low_power'],
                'feature': r['feature'], 'match_p_bonf': r['enrich_p_bonf'],
                'n_associations': len(candidates) if not r['low_power'] else None,
                'candidates': candidates, 'identity_raw_p': i['empirical_p'],
                'identity_joint_holm_p': i['joint_holm_p']}
        counts = {}
        for arm in ARMS[:3]:
            for candidate in summary['programs'][name][arm]['candidates']:
                counts[candidate] = counts.get(candidate, 0) + 1
        summary['computational_candidates_in_at_least_two_all_cell_seeds'][name] = sorted(
            c for c, n in counts.items() if n >= 2)
    for arm in ARMS:
        cycling = records[arm]['programs']['Cycling']
        if cycling['low_power']:
            summary['cycling_direction_comparisons'][arm] = {'matched': False}
            continue
        descriptive_up = [r for r in cycling['per_perturbation'] if r['auc'] >= .55]
        called = [r for r in cycling['per_perturbation'] if r['joint_modulates']]
        def count(rows):
            return {'n': len(rows), 'S_auc_le_045': sum(float(phase[r['pert']]['S_auc']) <= .45 for r in rows),
                    'G2M_auc_le_045': sum(float(phase[r['pert']]['G2M_auc']) <= .45 for r in rows),
                    'both_auc_le_045': sum(float(phase[r['pert']]['S_auc']) <= .45 and
                                          float(phase[r['pert']]['G2M_auc']) <= .45 for r in rows)}
        summary['cycling_direction_comparisons'][arm] = {
            'matched': True, 'auc_floor_only': count(descriptive_up), 'joint_association_calls': count(called)}
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--adata', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    args = parser.parse_args()
    if args.out_dir.exists():
        parser.error('choose a fresh run directory')
    args.out_dir.mkdir(parents=True)
    plan_path = ROOT / 'config/norman_shared_inference_plan.json'
    plan = json.loads(plan_path.read_text())
    sources = [Path(__file__).resolve(), plan_path, ROOT / 'config/human_markers.json'] + [ROOT / 'src' / n for n in
        ['norman_direct_pca.py', 'causal_pipeline.py', 'causal_grounding.py', 'marker_program_baseline.py',
         'program_modulation.py', 'program_modulation_specificity.py']]
    inputs = [ROOT / 'data_norman_continuation' / arm / 'seed0.json' for arm in ['direct_pca', 'direct_nmf']]
    phase_path = ROOT / 'data_norman_continuation/phase_diagnostic/comparisons.csv'
    sources += inputs + [phase_path]
    before = {str(p.relative_to(ROOT)): sha(p) for p in sources}
    if sha(args.adata) != plan['input_sha256']:
        raise ValueError('input SHA mismatch')
    from causal_pipeline import load_perturbseq, _topk_input_metadata
    import anndata as ad
    panels = json.loads((ROOT / 'config/human_markers.json').read_text())['marker_sets']
    x, genes, labels, _ = load_perturbseq(args.adata, 'perturbation', 'control', n_hvg=4000,
        keep_genes={g for panel in panels.values() for g in panel})
    if x.shape != (111445, 4081):
        raise ValueError('input expression shape differs')
    metadata = _topk_input_metadata(x)
    import torch
    checkpoint = torch.load('/workspace/research-runs/norman-sae-primary/seed0.pt',
                            map_location='cpu', weights_only=True)
    if metadata != checkpoint['input']:
        raise ValueError('preprocessed input bytes differ from the actual primary fit')
    del checkpoint
    a = ad.read_h5ad(args.adata, backed='r')
    try:
        if a.obs['gemgroup'].isna().any():
            raise ValueError('missing gemgroup')
        strata = a.obs['gemgroup'].astype(str).to_numpy()
        ids = a.obs_names.astype(str).to_numpy()
    finally:
        a.file.close()
    baseline = json.loads(inputs[0].read_text())
    records, artifacts = {}, {}
    print(json.dumps({'stage': 'loaded', 'shape': list(x.shape)}), flush=True)
    for index, arm in enumerate(ARMS):
        start = time.perf_counter()
        if index < 4:
            directory = ('norman-sae-primary' if index < 3 else 'norman-sae-control-sensitivity')
            seed = index if index < 3 else 0
            run = Path('/workspace/research-runs') / directory
            old_path = run / f'seed{seed}.json'
            old = json.loads(old_path.read_text())
            score_path = run / f'seed{seed}_activations.npy'
            if sha(score_path) != old['provenance']['activations_sha256']:
                raise ValueError('activation SHA mismatch')
            scores = np.load(score_path, mmap_mode='r', allow_pickle=False)
            artifacts[arm] = {'activation_path': str(score_path), 'activation_sha256': sha(score_path),
                              'original_result_sha256': sha(old_path), 'fit': old['sae']}
        elif index == 4:
            from sklearn.decomposition import PCA
            old = baseline
            original_model = Path(old['provenance']['model_path'])
            if sha(original_model) != old['provenance']['model_sha256']:
                raise ValueError('original PCA model SHA mismatch')
            settings = old['fit']['settings']
            model = PCA(n_components=settings['n_components'], svd_solver=settings['svd_solver'],
                        iterated_power=settings['iterated_power'], random_state=settings['random_state'], whiten=False)
            scores = model.fit_transform(x)
            score_path = args.out_dir / 'reproduced_pca_fit_transform.npy'
            with score_path.open('xb') as stream:
                np.save(stream, scores, allow_pickle=False)
            artifacts[arm] = {'reproduced_score_path': str(score_path), 'reproduced_score_sha256': sha(score_path),
                              'original_model_sha256': sha(original_model), 'fit': old['fit']}
        else:
            old = json.loads(inputs[1].read_text())
            score_path = Path(old['provenance']['actual_fitted_scores_path'])
            if sha(score_path) != old['provenance']['scores_sha256']:
                raise ValueError('actual NMF factor SHA mismatch')
            scores = np.load(score_path, mmap_mode='r', allow_pickle=False)
            artifacts[arm] = {'activation_path': str(score_path), 'activation_sha256': sha(score_path), 'fit': old['fit']}
        r = analyze(scores, x, genes, labels, panels, strata, seed=0, n_shuffle=9999)
        r['analysis'] = 'common_conditional_downstream_inference'
        for which in ['discovery', 'inference']:
            r['control_split'][which + '_cell_ids'] = ids[r['control_split'].pop(which + '_rows')].tolist()
        if r['control_split'] != baseline['control_split']:
            raise ValueError('the declared control split differs')
        if index >= 4:
            require_reproduced(r, old)
            artifacts[arm]['archived_association_inference_reproduced'] = True
        records[arm] = r
        with (args.out_dir / f'{arm}.json').open('x') as stream:
            json.dump(r, stream, indent=2, allow_nan=False)
            stream.write('\n')
        print(json.dumps({'stage': 'arm_complete', 'arm': arm,
                          'elapsed_seconds': time.perf_counter() - start}), flush=True)
        del scores
        gc.collect()
    combined, association_n, identity_n = correct_complete_family(records)
    if (association_n, identity_n) != (5670, 54):
        raise ValueError('complete joint families differ')
    with phase_path.open() as stream:
        phase = {r['pert']: r for r in csv.DictReader(stream)}
    summary = summarize(combined, phase)
    if before != {str(p.relative_to(ROOT)): sha(p) for p in sources}:
        raise RuntimeError('source or input archive changed during execution')
    result = {'analysis': 'exploratory_six_representation_common_inference',
        'association_family_size': association_n, 'identity_family_size': identity_n,
        'representations': combined, 'summary': summary, 'scope': plan['limits'],
        'provenance': {'completed_utc': datetime.now(timezone.utc).isoformat(),
            'input_sha256': plan['input_sha256'], 'preprocessed_input': metadata,
            'source_sha256': before, 'model_artifacts': artifacts,
            'git_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
            'versions': {n: importlib.metadata.version(n) for n in
                         ['numpy', 'scipy', 'scikit-learn', 'torch', 'scanpy', 'anndata', 'statsmodels']}}}
    with (args.out_dir / 'results.json').open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    with (args.out_dir / 'summary.json').open('x') as stream:
        json.dump(summary, stream, indent=2, allow_nan=False)
        stream.write('\n')
    columns = ['arm', 'program', 'pert', 'status', 'auc', 'mw_p', 'mw_q', 'joint_bh_q', 'joint_modulates',
               'S_auc', 'G2M_auc']
    with (args.out_dir / 'comparisons.csv').open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for arm in ARMS:
            for name, panel in combined[arm]['programs'].items():
                for row in panel['per_perturbation']:
                    writer.writerow({'arm': arm, 'program': name, **{k: row[k] for k in columns[2:9]},
                                     'S_auc': phase[row['pert']]['S_auc'], 'G2M_auc': phase[row['pert']]['G2M_auc']})
    print(json.dumps({'stage': 'complete', 'association_family': association_n, 'identity_family': identity_n}), flush=True)


if __name__ == '__main__':
    main()
