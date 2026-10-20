#!/usr/bin/env python3
"""Test the exact historical K562 program across all completed reference models."""
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

from causal_pipeline import load_perturbseq, _topk_input_metadata
from marker_program_baseline import normalize_counts_to_markers, marker_program_baseline
from norman_direct_pca import analyze, sha
from norman_nmf_convergence import ARMS, NEW_ARM

ROOT = Path(__file__).resolve().parents[1]
ALL_ARMS = ARMS + [NEW_ARM, 'direct_markers']


def correct_expanded(records):
    if set(records) != set(ALL_ARMS):
        raise ValueError('all seven learned models and direct markers are required')
    out = copy.deepcopy(records)
    panels = list(out[ALL_ARMS[0]]['programs'])
    expected = None
    rows, identity = [], []
    for arm in ALL_ARMS:
        r = out[arm]
        if set(r['programs']) != set(panels) or set(r['identity']) != set(panels):
            raise ValueError('complete expanded panel families are required')
        for name in panels:
            these = r['programs'][name]['per_perturbation']
            names = [v['pert'] for v in these]
            if len(names) != len(set(names)):
                raise ValueError('duplicate perturbation hypothesis')
            if expected is None:
                expected = set(names)
            if set(names) != expected:
                raise ValueError('all eligible perturbations must remain')
            rows.extend(these)
            identity.append(r['identity'][name])
    q = multipletests([v['mw_p'] for v in rows], method='fdr_bh')[1]
    for r, value in zip(rows, q):
        r['context_bh_q'] = float(value)
        r['context_called'] = bool(value < .05 and r['auc'] >= .55 and r['status'] == 'tested')
    q = multipletests([v['empirical_p'] for v in identity], method='holm')[1]
    for r, value in zip(identity, q):
        r['context_holm_p'] = float(value)
    return out, len(rows), len(identity)


def require_direct_recovery(new, old):
    if new['control_split'] != old['control_split']:
        raise ValueError('direct-marker split differs')
    for gene, parameters in old['discovery_gene_parameters'].items():
        actual = new['discovery_gene_parameters'][gene]
        for key, expected in parameters.items():
            if isinstance(expected, (float, int)) and not isinstance(expected, bool):
                if not np.isclose(actual[key], expected, rtol=1e-10, atol=1e-10):
                    raise ValueError('discovery parameter was not reproduced')
            elif actual[key] != expected:
                raise ValueError('discovery eligibility was not reproduced')
    for name, panel in old['programs'].items():
        current = new['programs'][name]
        if current['eligible_markers'] != panel['eligible_markers'] or current['low_power'] != panel['low_power']:
            raise ValueError('direct-marker coverage differs')
        by = {r['pert']: r for r in current['per_perturbation']}
        for r in panel['per_perturbation']:
            value = by[r['pert']]
            if (value['status'] != r['status'] or abs(value['auc'] - r['auc']) > 1e-7
                    or abs(value['mw_p'] - r['mw_p']) > 1e-10):
                raise ValueError('original direct-marker inference not reproduced')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--adata', type=Path, required=True)
    p.add_argument('--out-dir', type=Path, required=True)
    args = p.parse_args()
    if args.out_dir.exists():
        p.error('choose a fresh output directory')
    plan_path = ROOT / 'config/norman_k562_context_plan.json'
    plan = json.loads(plan_path.read_text())
    follow_path = ROOT / 'data_norman_continuation/nmf_convergence/results.json.gz'
    validated = follow_path.with_name('validation.json')
    if not validated.exists() or json.loads(validated.read_text())['status'] != 'PASS':
        p.error('wait for the complete, validated bounded NMF follow-up')
    common_path = ROOT / 'data_norman_continuation/shared_inference/results.json.gz'
    old_marker_path = ROOT / 'data_norman_continuation/marker_baseline/gemgroup_seed0.json'
    sources = [Path(__file__).resolve(), plan_path, follow_path, validated, common_path,
               old_marker_path, ROOT / 'config/human_markers.json'] + [ROOT / 'src' / n for n in
               ['norman_direct_pca.py', 'norman_nmf_convergence.py', 'norman_shared_inference.py',
                'causal_pipeline.py', 'causal_grounding.py', 'marker_program_baseline.py',
                'program_modulation.py', 'program_modulation_specificity.py']]
    before = {str(a.relative_to(ROOT)): sha(a) for a in sources}
    if sha(args.adata) != plan['input_sha256']:
        raise ValueError('source input SHA differs')
    with gzip.open(follow_path, 'rt') as f:
        follow = json.load(f)
    with gzip.open(common_path, 'rt') as f:
        common = json.load(f)
    old = json.loads(old_marker_path.read_text())
    original_panels = json.loads((ROOT / 'config/human_markers.json').read_text())['marker_sets']
    name, markers = plan['program_name'], plan['genes']
    panels = {**original_panels, name: markers}
    # Do not silently change the fixed model input by force-keeping the new markers.
    x, genes, labels, _ = load_perturbseq(args.adata, 'perturbation', 'control', n_hvg=4000,
        keep_genes={g for panel in original_panels.values() for g in panel})
    if _topk_input_metadata(x) != common['provenance']['preprocessed_input']:
        raise ValueError('actual fitted expression bytes differ')
    import anndata as ad
    a = ad.read_h5ad(args.adata)
    strata = a.obs['gemgroup'].astype(str).to_numpy()
    ids = a.obs_names.astype(str).to_numpy()
    marker_x, marker_genes, normalization = normalize_counts_to_markers(a.X, a.var_names.to_numpy(), panels)
    del a
    direct = marker_program_baseline(marker_x, marker_genes, labels, panels,
        strata=strata, cell_ids=ids, seed=0, n_identity_shuffle=9999)
    require_direct_recovery(direct, old)
    print(json.dumps({'stage': 'direct_complete', 'original_945_comparisons_reproduced': True}), flush=True)
    records, artifacts = {}, {}
    for arm in ARMS + [NEW_ARM]:
        if arm == NEW_ARM:
            path = Path(follow['provenance']['new_scores_path'])
            expected = follow['provenance']['new_scores_sha256']
        else:
            info = common['provenance']['model_artifacts'][arm]
            key = 'reproduced_score' if arm == 'direct_pca_seed0' else 'activation'
            path = Path(info[key + '_path'])
            expected = info[key + '_sha256']
        if sha(path) != expected:
            raise ValueError('actual score SHA differs')
        z = np.load(path, mmap_mode='r', allow_pickle=False)
        context = analyze(z, x, genes, labels, {name: markers}, strata, seed=0, n_shuffle=9999)
        for which in ['discovery', 'inference']:
            context['control_split'][which + '_cell_ids'] = ids[context['control_split'].pop(which + '_rows')].tolist()
        if context['control_split'] != old['control_split']:
            raise ValueError('context annotation/test split differs')
        base = follow['representations'][arm]
        records[arm] = {'programs': {**base['programs'], **context['programs']},
                        'identity': {**base['identity_heterogeneity']['programs'],
                                     **context['identity_heterogeneity']['programs']},
                        'context_annotation': {k: v for k, v in context.items() if k not in ['programs', 'identity_heterogeneity', 'control_split']}}
        artifacts[arm] = {'path': str(path), 'sha256': expected}
        print(json.dumps({'stage': 'learned_context_complete', 'arm': arm}), flush=True)
        del z
    records['direct_markers'] = {'programs': direct['programs'], 'identity': direct['identity_heterogeneity']['programs']}
    corrected, a_n, i_n = correct_expanded(records)
    if (a_n, i_n) != (8400, 80):
        raise ValueError('complete expanded correction families differ')
    summary = {program: {arm: {'matched': not v['programs'][program]['low_power'],
        'n_associations': sum(r['context_called'] for r in v['programs'][program]['per_perturbation'])
                          if not v['programs'][program]['low_power'] else None,
        'candidates': [r['pert'] for r in v['programs'][program]['per_perturbation'] if r['context_called']],
        'identity_holm_p': v['identity'][program]['context_holm_p']}
        for arm, v in corrected.items()} for program in panels}
    if before != {str(a.relative_to(ROOT)): sha(a) for a in sources}:
        raise RuntimeError('scientific source changed during context execution')
    output = {'analysis': 'exploratory_historical_k562_context_extension', 'association_family_size': a_n,
        'identity_family_size': i_n, 'representations': corrected, 'summary': summary,
        'context_program': name, 'context_markers': markers,
        'present_in_fixed_model_input': [g for g in markers if g in genes],
        'absent_from_fixed_model_input': [g for g in markers if g not in genes],
        'direct_fit': direct, 'normalization': normalization, 'original_direct_marker_recovery': True,
        'scope': plan['limits'], 'provenance': {'completed_utc': datetime.now(timezone.utc).isoformat(),
            'source_sha256': before, 'input_sha256': plan['input_sha256'], 'actual_score_files': artifacts,
            'git_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
            'versions': {n: importlib.metadata.version(n) for n in ['numpy', 'scipy', 'scanpy', 'anndata', 'statsmodels']}}}
    args.out_dir.mkdir(parents=True)
    with (args.out_dir / 'results.json').open('x') as f:
        json.dump(output, f, indent=2, allow_nan=False)
        f.write('\n')
    with (args.out_dir / 'comparisons.csv').open('x', newline='') as f:
        keys = ['arm', 'program', 'pert', 'status', 'auc', 'mw_p', 'context_bh_q', 'context_called']
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for arm, v in corrected.items():
            for program, panel in v['programs'].items():
                for row in panel['per_perturbation']:
                    writer.writerow({'arm': arm, 'program': program, **{k: row[k] for k in keys[2:]}})
    print(json.dumps({'stage': 'complete', 'context_program_summary': summary[name]}), flush=True)


if __name__ == '__main__':
    main()
