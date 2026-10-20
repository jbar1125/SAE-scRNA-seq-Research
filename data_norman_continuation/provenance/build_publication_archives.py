"""Create checked review archives without changing source data or model files."""
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile

ROOT = Path('/workspace/SAE-scRNA-seq-Research')
OUT = Path('/workspace/publication_review_2026-10-20')
OUT.mkdir(exist_ok=True)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def archive(name, files):
    target = OUT / (name + '.tar.gz')
    assert not target.exists(), target
    expected = {str(p.relative_to('/workspace')): {'sha256': digest(p), 'bytes': p.stat().st_size}
                for p in files}
    for p in files:
        assert p.is_file() and not p.is_symlink()
    with target.open('xb') as raw:
        with gzip.GzipFile(fileobj=raw, mode='wb', compresslevel=1, mtime=0) as zipped:
            with tarfile.open(fileobj=zipped, mode='w|') as tar:
                for p in files:
                    tar.add(p, arcname=str(p.relative_to('/workspace')), recursive=False)
    checked = {}
    with tarfile.open(target, 'r|gz') as tar:
        for member in tar:
            assert member.isfile() and member.name in expected
            h = hashlib.sha256()
            stream = tar.extractfile(member)
            for block in iter(lambda: stream.read(8*1024*1024), b''):
                h.update(block)
            checked[member.name] = {'sha256': h.hexdigest(), 'bytes': member.size}
    assert checked == expected
    record = {'filename': target.name, 'bytes': target.stat().st_size,
              'sha256': digest(target), 'members': expected, 'roundtrip_validation': 'PASS'}
    print(json.dumps({'archive': target.name, 'bytes': record['bytes'], 'status': 'PASS'}), flush=True)
    return record


def main():
    records = []
    for arm, seeds in [('norman-sae-primary', range(3)), ('norman-sae-control-sensitivity', [0])]:
        for seed in seeds:
            run = Path('/workspace/research-runs') / arm
            files = [run / f'seed{seed}.{suffix}' for suffix in ['pt', 'json', 'log']]
            files.append(run / f'seed{seed}_activations.npy')
            records.append(archive(arm + f'-seed{seed}', files))
    files = [Path('/workspace/research-runs/norman-direct-pca-v1/seed0_model.npz'),
             Path('/workspace/research-runs/norman-direct-nmf-v2/seed0_model.npz'),
             Path('/workspace/research-runs/norman-direct-nmf-v2/seed0_fitted_scores.npy'),
             Path('/workspace/research-runs/norman-shared-inference-v1/reproduced_pca_fit_transform.npy')]
    records.append(archive('norman-direct-factor-models', files))
    source = Path('/workspace/datasets/sae-scrna/norman_pert_unique.h5ad')
    assert digest(source) == 'a09744916427b08d75dade996d059f717e4812463292a1e5a4e5450589cf2c6a'
    records.append(archive('norman-prepared-unique-counts', [source]))
    manifest = {'created_utc': datetime.now(timezone.utc).isoformat(),
                'git_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                'purpose': 'Review archives of actual large artifacts. Original files unchanged. These are not biological replications or a new experiment.',
                'archives': records,
                'restore': 'Inspect members first. Relative paths start research-runs/ or datasets/. Restore only into a fresh workspace; do not overwrite existing data. Recorded absolute source paths are cloud-instance provenance.',
                'publication_status': 'Local verified archives; public release or external deposition not implied.'}
    (OUT / 'large_artifact_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print('COMPLETE', flush=True)


if __name__ == '__main__':
    main()
