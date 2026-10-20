#!/usr/bin/env python3
"""Bounded warm-restart sensitivity, preserving the original capped NMF evidence."""
import argparse
import copy
import csv
from datetime import datetime, timezone
import gzip
import importlib.metadata
import json
from pathlib import Path
import subprocess
import time
import warnings

import numpy as np
from sklearn.decomposition import NMF
from sklearn.exceptions import ConvergenceWarning
from statsmodels.stats.multitest import multipletests

from causal_pipeline import load_perturbseq, _topk_input_metadata
from norman_direct_pca import analyze, sha
from norman_shared_inference import ARMS

ROOT = Path(__file__).resolve().parents[1]
NEW_ARM = 'warm_restarted_nmf_seed0'


def warm_restart(x, w, h, settings):
    if (x.ndim != 2 or w.shape != (len(x), settings['n_components']) or
            h.shape != (settings['n_components'], x.shape[1])):
        raise ValueError('actual factor and input dimensions must align')
    if any(not np.isfinite(a).all() or np.any(a < 0) for a in [x, w, h]):
        raise ValueError('finite nonnegative input and factors required')
    model = NMF(**settings, verbose=1)
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        scores = model.fit_transform(x, W=w.copy(), H=h.copy())
    info = {'settings': settings, 'iterations': int(model.n_iter_),
        'fit_elapsed_seconds': time.perf_counter() - started,
        'convergence_warning_absent': not any(issubclass(a.category, ConvergenceWarning) for a in caught),
        'warnings': [{'category': a.category.__name__, 'message': str(a.message)} for a in caught],
        'reconstruction_frobenius_error': float(model.reconstruction_err_),
        'stopping_reference': 'New custom-initialization error; not bit-identical uninterrupted continuation.'}
    return scores, model.components_, info


def joint_followup(original, new):
    if set(original) != set(ARMS) or NEW_ARM in original:
        raise ValueError('the complete original six arms are required')
    records = copy.deepcopy({**original, NEW_ARM: new})
    names = list(records[ARMS[0]]['programs'])
    all_rows, all_identity = [], []
    expected = None
    for arm, r in records.items():
        if set(r['programs']) != set(names) or set(r['identity_heterogeneity']['programs']) != set(names):
            raise ValueError('all planned panels are required')
        for name in names:
            rows = r['programs'][name]['per_perturbation']
            perts = [a['pert'] for a in rows]
            if len(set(perts)) != len(perts):
                raise ValueError('duplicate comparison')
            if expected is None:
                expected = set(perts)
            if set(perts) != expected:
                raise ValueError('complete perturbation family required')
            all_rows.extend(rows)
            all_identity.append(r['identity_heterogeneity']['programs'][name])
    q = multipletests([a['mw_p'] for a in all_rows], method='fdr_bh')[1]
    for a, value in zip(all_rows, q):
        a['seven_arm_joint_q'] = float(value)
        a['seven_arm_modulates'] = bool(value < .05 and a['auc'] >= .55 and a['status'] == 'tested')
    q = multipletests([a['empirical_p'] for a in all_identity], method='holm')[1]
    for a, value in zip(all_identity, q):
        a['seven_arm_holm_p'] = float(value)
    return records, len(all_rows), len(all_identity)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--adata', type=Path, required=True)
    p.add_argument('--out-dir', type=Path, required=True)
    args = p.parse_args()
    if args.out_dir.exists():
        p.error('choose a fresh output directory')
    args.out_dir.mkdir(parents=True)
    plan_path = ROOT / 'config/norman_nmf_convergence_plan.json'
    plan = json.loads(plan_path.read_text())
    old_path = ROOT / 'data_norman_continuation/direct_nmf/seed0.json'
    common_path = ROOT / 'data_norman_continuation/shared_inference/results.json.gz'
    sources = [Path(__file__).resolve(), plan_path, old_path, common_path, ROOT / 'config/human_markers.json'] + [ROOT / 'src' / n for n in
        ['norman_direct_pca.py', 'norman_shared_inference.py', 'causal_pipeline.py', 'causal_grounding.py',
         'marker_program_baseline.py', 'program_modulation.py', 'program_modulation_specificity.py']]
    before = {str(a.relative_to(ROOT)): sha(a) for a in sources}
    old = json.loads(old_path.read_text())
    with gzip.open(common_path, 'rt') as stream:
        common = json.load(stream)
    model_path = Path(old['provenance']['model_path'])
    score_path = Path(old['provenance']['actual_fitted_scores_path'])
    if sha(model_path) != plan['initialization']['components_sha256'] or sha(score_path) != plan['initialization']['actual_scores_sha256']:
        raise ValueError('original actual NMF artifact SHA differs')
    if sha(args.adata) != plan['input_sha256']:
        raise ValueError('source input SHA differs')
    panels = json.loads((ROOT / 'config/human_markers.json').read_text())['marker_sets']
    x, genes, labels, _ = load_perturbseq(args.adata, 'perturbation', 'control', n_hvg=4000,
        keep_genes={g for panel in panels.values() for g in panel})
    if _topk_input_metadata(x) != common['provenance']['preprocessed_input']:
        raise ValueError('preprocessed input bytes differ')
    w = np.load(score_path, mmap_mode='r', allow_pickle=False)
    h = np.load(model_path, allow_pickle=False)['components']
    original_w, original_h = sha(score_path), sha(model_path)
    print(json.dumps({'stage': 'loaded', 'shape': list(x.shape), 'initial_error': old['fit']['reconstruction_frobenius_error'],
                      'settings': plan['fit'], 'threads': 4}), flush=True)
    scores, components, info = warm_restart(x, w, h, plan['fit'])
    info['original_capped_error'] = old['fit']['reconstruction_frobenius_error']
    info['relative_reconstruction_error'] = float(info['reconstruction_frobenius_error'] / np.linalg.norm(x))
    if info['reconstruction_frobenius_error'] > info['original_capped_error'] * (1 + 1e-5):
        raise ValueError('warm restart increased reconstruction error')
    if (sha(score_path), sha(model_path)) != (original_w, original_h):
        raise ValueError('original artifacts changed')
    print(json.dumps({'stage': 'fit_complete', 'fit': info}), flush=True)
    new_model = args.out_dir / 'model.npz'
    new_scores = args.out_dir / 'fitted_scores.npy'
    with new_model.open('xb') as stream:
        np.savez_compressed(stream, components=components)
    with new_scores.open('xb') as stream:
        np.save(stream, scores, allow_pickle=False)
    import anndata as ad
    a = ad.read_h5ad(args.adata, backed='r')
    try:
        strata = a.obs['gemgroup'].astype(str).to_numpy()
        ids = a.obs_names.astype(str).to_numpy()
    finally:
        a.file.close()
    result = analyze(scores, x, genes, labels, panels, strata, seed=0, n_shuffle=9999)
    for which in ['discovery', 'inference']:
        result['control_split'][which + '_cell_ids'] = ids[result['control_split'].pop(which + '_rows')].tolist()
    if result['control_split'] != old['control_split']:
        raise ValueError('original control split differs')
    records, n_tests, n_identity = joint_followup(common['representations'], result)
    if (n_tests, n_identity) != (6615, 63):
        raise ValueError('declared joint families differ')
    if before != {str(a.relative_to(ROOT)): sha(a) for a in sources}:
        raise RuntimeError('scientific source changed during fit')
    summary = {name: {arm: {'matched': not r['programs'][name]['low_power'],
        'n_associations': sum(v['seven_arm_modulates'] for v in r['programs'][name]['per_perturbation'])
                          if not r['programs'][name]['low_power'] else None,
        'identity_holm_p': r['identity_heterogeneity']['programs'][name]['seven_arm_holm_p']}
        for arm, r in records.items()} for name in panels}
    output = {'analysis': 'exploratory_warm_restart_nmf_convergence', 'fit': info,
        'association_family_size': n_tests, 'identity_family_size': n_identity,
        'representations': records, 'summary': summary, 'scope': plan['interpretation'],
        'provenance': {'completed_utc': datetime.now(timezone.utc).isoformat(), 'input_sha256': plan['input_sha256'],
            'source_sha256': before, 'initial_model_sha256': original_h, 'initial_scores_sha256': original_w,
            'new_model_path': str(new_model), 'new_model_sha256': sha(new_model),
            'new_scores_path': str(new_scores), 'new_scores_sha256': sha(new_scores),
            'git_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
            'versions': {n: importlib.metadata.version(n) for n in ['numpy', 'scipy', 'scikit-learn', 'scanpy', 'anndata', 'statsmodels']}}}
    with (args.out_dir / 'results.json').open('x') as stream:
        json.dump(output, stream, indent=2, allow_nan=False)
        stream.write('\n')
    with (args.out_dir / 'comparisons.csv').open('x', newline='') as stream:
        fields = ['arm', 'program', 'pert', 'status', 'auc', 'mw_p', 'seven_arm_joint_q', 'seven_arm_modulates']
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for arm, r in records.items():
            for name, panel in r['programs'].items():
                for row in panel['per_perturbation']:
                    writer.writerow({'arm': arm, 'program': name, **{k: row[k] for k in fields[2:]}})
    print(json.dumps({'stage': 'complete', 'fit': info, 'new_arm_counts': {n: r[NEW_ARM] for n, r in summary.items()}}), flush=True)


if __name__ == '__main__':
    main()
