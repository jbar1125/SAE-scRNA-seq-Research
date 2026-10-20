"""Read and validate the completed declared Norman model artifacts; archive control seed.

Run from the repository root with its src directory importable. Large files stay
outside Git. Existing evidence is compared byte for byte, never silently replaced.
"""
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sys

import numpy as np
import torch
from statsmodels.stats.multitest import multipletests

ROOT = Path.cwd()
sys.path.insert(0, str(ROOT / 'src'))
from archive_norman_cpu_results import sha, validate_record, verify_large_artifacts, joint_summary


def copy_identical(source, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        assert sha(source) == sha(destination), str(destination)
    else:
        with destination.open('xb') as stream:
            with source.open('rb') as incoming:
                shutil.copyfileobj(incoming, stream)


def main():
    primary = []
    checked = []
    expected_input = 'a09744916427b08d75dade996d059f717e4812463292a1e5a4e5450589cf2c6a'
    assert sha('/workspace/datasets/sae-scrna/norman_pert_unique.h5ad') == expected_input
    for seed in range(3):
        path = ROOT / f'data_norman_continuation/sae_primary/seed{seed}.json'
        r = json.loads(path.read_text())
        validate_record(r, seed)
        verify_large_artifacts(r, seed)
        assert sha(path) == sha(f'/workspace/research-runs/norman-sae-primary/seed{seed}.json')
        for name, fingerprint in r['provenance']['source_sha256'].items():
            assert sha(ROOT / 'src' / name) == fingerprint, name
        primary.append(r)
        checked.append({'arm': 'all', 'seed': seed, 'epochs': 200,
                        'training_cells': 111445, 'activated_features': r['sae']['features_ever_active'],
                        'checkpoint_sha256': r['provenance']['checkpoint_sha256'],
                        'activations_sha256': r['provenance']['activations_sha256']})
    assert joint_summary(primary) == json.loads((ROOT / 'data_norman_continuation/sae_primary/joint_summary.json').read_text())
    run = Path('/workspace/research-runs/norman-sae-control-sensitivity')
    control = json.loads((run / 'seed0.json').read_text())
    s, p = control['sae'], control['provenance']
    assert s['epochs'] == s['epochs_completed'] == 200
    assert (s['d'], s['latent_dim'], s['k'], s['batch']) == (4081, 2048, 32, 1024)
    assert (s['trained_on'], s['n_train_cells']) == ('control', 11855)
    assert p['arguments']['seed'] == 0 and p['input_sha256'] == expected_input
    assert p['source_sha256'] == primary[0]['provenance']['source_sha256']
    for kind, filename in [('checkpoint', 'seed0.pt'), ('activations', 'seed0_activations.npy')]:
        assert sha(run / filename) == p[kind + '_sha256']
    checkpoint = torch.load(run / 'seed0.pt', map_location='cpu', weights_only=True)
    assert checkpoint['completed_epochs'] == checkpoint['target_epochs'] == 200
    assert checkpoint['input']['shape'] == [11855, 4081]
    assert (checkpoint['training_params']['seed'], checkpoint['training_params']['latent_dim']) == (0, 2048)
    acts = np.load(run / 'seed0_activations.npy', mmap_mode='r', allow_pickle=False)
    assert acts.shape == (111445, 2048) and acts.dtype == np.float32
    active = np.zeros(acts.shape[1], dtype=bool)
    for start in range(0, len(acts), 4096):
        block = acts[start:start + 4096]
        assert np.isfinite(block).all() and (block >= 0).all()
        assert np.all((block > 0).sum(axis=1) <= 32)
        active |= (block > 0).any(axis=0)
    assert int(active.sum()) == s['features_ever_active']
    identity = control['identity_heterogeneity']
    assert identity['n_planned_programs'] == 9 and identity['n_permutations'] == 999
    assert identity['stratified'] and p['arguments']['strata_col'] == 'gemgroup'
    panels = list(primary[0]['programs'])
    assert set(panels) == set(control['programs']) == set(identity['programs'])
    for r in identity['programs'].values():
        assert .001 <= r['empirical_p'] <= 1
        if not r['matched']:
            assert r['empirical_p'] == 1
    checked.append({'arm': 'control', 'seed': 0, 'epochs': 200,
                    'training_cells': 11855, 'activated_features': int(active.sum()),
                    'checkpoint_sha256': p['checkpoint_sha256'],
                    'activations_sha256': p['activations_sha256']})
    archive = ROOT / 'data_norman_continuation/sae_control_sensitivity'
    for suffix in ['json', 'log']:
        copy_identical(run / f'seed0.{suffix}', archive / f'seed0.{suffix}')
    summary = {'scope': 'One exploratory control-only seed. No biological replication or three-seed sensitivity.',
               'programs': {n: {'matched': not r['low_power'], 'candidate_count': r.get('n_modulators'),
                                'candidates': r.get('modulators', []),
                                'feature_background_percentile': r.get('feature_background_percentile'),
                                'identity_holm_9': identity['programs'][n]['holm_p']}
                            for n, r in control['programs'].items()}}
    summary_path = archive / 'summary.json'
    payload = json.dumps(summary, indent=2, allow_nan=False) + '\n'
    if summary_path.exists():
        assert summary_path.read_text() == payload
    else:
        summary_path.write_text(payload)
    # Explicitly exploratory combination of the complete four-fit, nine-panel family.
    all_records = primary + [control]
    raw = [r['identity_heterogeneity']['programs'][n]['empirical_p'] for r in all_records for n in panels]
    joint = multipletests(raw, method='holm')[1]
    validation = {'completed_utc': datetime.now(timezone.utc).isoformat(), 'status': 'PASS',
                  'input_sha256': expected_input, 'models': checked,
                  'primary_joint_27_reproduced': True,
                  'cross_arm_identity_family': {'size': 36, 'method': 'Holm',
                    'scope': 'Exploratory four-fit combination; conditional cell-level heterogeneity only.',
                    'records': [{'arm': 'all' if i < 3 else 'control', 'seed': i if i < 3 else 0,
                                 'program': n, 'holm_p': float(joint[9*i+j])}
                                for i in range(4) for j, n in enumerate(panels)]},
                  'validator_sha256': sha(Path(__file__))}
    output = ROOT / 'data_norman_continuation/provenance/completed_models_validation.json'
    with output.open('x') as stream:
        json.dump(validation, stream, indent=2, allow_nan=False)
        stream.write('\n')
    copy_identical(Path(__file__), output.with_name('validate_completed_models.py'))
    print(json.dumps({'status': 'PASS', 'models': checked, 'control_summary': summary}, indent=2))


if __name__ == '__main__':
    main()
