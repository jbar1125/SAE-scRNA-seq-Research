"""Adapt a verified direct fit to the unchanged sensitivity CLI, preserving provenance."""
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path('/workspace/SAE-scRNA-seq-Research')
RUN = Path('/workspace/research-runs/norman-k562-context-sensitivity-v1')
PRIMARY = ROOT / 'data_norman_continuation/k562_context/results.json.gz'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    if RUN.exists():
        raise FileExistsError(RUN)
    assert json.loads(PRIMARY.with_name('validation.json').read_text())['status'] == 'PASS'
    with gzip.open(PRIMARY, 'rt') as stream:
        parent = json.load(stream)
    primary = parent['direct_fit']
    primary['provenance'] = {'input': {'path': '/workspace/datasets/sae-scrna/norman_pert_unique.h5ad',
        'sha256': parent['provenance']['input_sha256']},
        'extraction': {'parent_archive_sha256': sha(PRIMARY), 'parent_field': 'direct_fit',
                       'adapter_sha256': sha(Path(__file__)),
                       'note': 'All fit fields copied losslessly. This new provenance describes extraction, not a new primary execution.'}}
    RUN.mkdir()
    adapted = RUN / 'primary.json'
    adapted.write_text(json.dumps(primary, indent=2) + '\n')
    with (RUN / 'core_run.log').open('x') as log:
        subprocess.run([sys.executable, str(ROOT / 'src/marker_baseline_sensitivity.py'),
            '--adata', '/workspace/datasets/sae-scrna/norman_pert_unique.h5ad',
            '--primary', str(adapted), '--out', str(RUN / 'results.json'),
            '--plan', str(ROOT / 'config/norman_k562_context_sensitivity_plan.json')],
            cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
    r = json.loads((RUN / 'results.json').read_text())
    assert (r['control_reference_family_size'], r['depth_family_size']) == (4200, 1050)
    assert len(r['control_references']) == 4 and r['provenance']['primary_reproduced_exactly']
    summary = {}
    for name, original in parent['representations']['direct_markers']['programs'].items():
        main_calls = {v['pert'] for v in original['per_perturbation'] if v['context_called']}
        passed = set(r['summary'][name]['primary_and_depth_and_every_reference'])
        intersection = sorted(main_calls & passed)
        summary[name] = {'eight_arm_primary_calls': len(main_calls),
            'primary_and_omission_and_depth_and_every_reference': len(intersection), 'candidates': intersection,
            'sensitivity_details': r['summary'][name]}
    name = parent['context_program']
    cases = {}
    for target in ['CBL', 'UBASH3B', 'PTPN12', 'KLF1']:
        rows = [row for condition, values in r['comparisons'].items()
                if condition.startswith('reference:') for row in values if row['program'] == name and row['pert'] == target]
        depth = next(row for row in r['comparisons']['depth_matched'] if row['program'] == name and row['pert'] == target)
        cases[target] = {'survives_full_intersection': target in summary[name]['candidates'],
                         'minimum_reference_auc': min(row['auc'] for row in rows), 'depth_matched_auc': depth['auc']}
    payload = {'summary': summary, 'preidentified_context_cases': cases,
        'context_program': name, 'scope': 'Post-primary descriptive robustness; same cells and scores. Depth can overadjust a mediator. References/technical partitions are not biological replication.',
        'provenance': {'created_utc': datetime.now(timezone.utc).isoformat(), 'parent_context_sha256': sha(PRIMARY),
                       'adapted_fit_sha256': sha(adapted), 'sensitivity_sha256': sha(RUN / 'results.json'),
                       'adapter_sha256': sha(Path(__file__))}}
    (RUN / 'summary.json').write_text(json.dumps(payload, indent=2) + '\n')
    with (RUN / 'comparisons.csv').open('x', newline='') as stream:
        keys = ['condition', 'program', 'pert', 'status', 'auc', 'mw_p', 'mw_q', 'modulates', 'survives_target_omission']
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        for condition, values in r['comparisons'].items():
            for row in values:
                writer.writerow({'condition': condition, **{key: row[key] for key in keys[1:]}})
    print(json.dumps({'program_counts': {n: [v['eight_arm_primary_calls'], v['primary_and_omission_and_depth_and_every_reference']] for n, v in summary.items()},
                      'context_cases': cases}), flush=True)


if __name__ == '__main__':
    main()
