# Reproducing the fresh Norman continuation

The cloud executions and all compact evidence are archived. This is a CLI research
workflow. No web server or GPU is required. A clean-machine rebuild of every
historical experiment has not been demonstrated; original causal JSONs and some
old structural weights remain missing. New Norman artifacts do not recreate them.

## Inputs and exact archived models

Use the Python environment and versions in [REPRODUCIBILITY.md](REPRODUCIBILITY.md).
All source-data normalization starts from raw counts. The actual public raw source
URL and observed fingerprint are in
[`norman_download_manifest.json`](../data_norman_continuation/provenance/norman_download_manifest.json).
The supported preparation is `src/prep_norman.py --in <raw> --out <new-prepared>`
followed by `src/prep_norman_unique.py --in <prepared> --out <new-unique>`.
The latter sums duplicate counts before normalization, preserves cell totals,
and retains source paths in file metadata. Do not overwrite original inputs.

The exact unique-count file used here has SHA-256
`a09744916427b08d75dade996d059f717e4812463292a1e5a4e5450589cf2c6a`.
Prepared HDF5 bytes include provenance metadata, so a different source path or
serialization version can change the file hash without changing count values.
Do not alter a hash merely to bypass a guard. Compare logical matrices, ordered
genes, cell IDs, labels and strata, and explicitly record any revised provenance.

Large verified review archives contain the actual four SAE checkpoints/activation
matrices, original PCA/NMF factors, recovered original PCA scores and the exact
prepared unique-count input. They currently reside under
`/workspace/publication_review_2026-10-20/`. GitHub release-asset upload is blocked
at the network proxy; these files have not been publicly deposited. Inspect the
archive manifest and member paths, and restore only into a fresh workspace.
The GitHub draft release is a review placeholder, not proof of uploaded assets.

## Fresh output paths

From the existing checkout and activated environment, choose one new run directory:

```bash
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export XDG_CACHE_HOME=/workspace/.cache
export MPLCONFIGDIR=/workspace/.cache/matplotlib
export NUMBA_CACHE_DIR=/workspace/.cache/numba
NORMAN_REPRO_RUN=$(mktemp -d /workspace/norman-reproduction.XXXXXX)
```

The following CLI paths were executed on the real dataset; the fresh directory
names below are examples, not a claim that this entire block ran as one command.
The data path refers to this cloud instance's exact cached source.

```bash
python src/marker_program_baseline.py \
  --adata /workspace/datasets/sae-scrna/norman_pert_unique.h5ad \
  --out "$NORMAN_REPRO_RUN/marker.json"
python src/marker_baseline_sensitivity.py \
  --adata /workspace/datasets/sae-scrna/norman_pert_unique.h5ad \
  --primary "$NORMAN_REPRO_RUN/marker.json" \
  --out "$NORMAN_REPRO_RUN/sensitivity.json"
python src/norman_direct_pca.py \
  --adata /workspace/datasets/sae-scrna/norman_pert_unique.h5ad \
  --out "$NORMAN_REPRO_RUN/pca.json" \
  --model-out "$NORMAN_REPRO_RUN/pca-model.npz"
python src/norman_target_rna.py \
  --adata /workspace/datasets/sae-scrna/norman_pert_unique.h5ad \
  --out-dir "$NORMAN_REPRO_RUN/target-rna"
```

Marker analysis on the unique-count input uses the same count-conserving policy;
the original marker run instead requested duplicate summation on the preceding
prepared input. Their input/policy metadata differ. Compare numeric estimates,
control IDs and fitted parameters rather than expecting identical provenance JSON.

The specialized common/attribution runners intentionally refer to the declared
archived models and comparison records, including canonical `/workspace` paths.
They do not automatically discover or substitute newly fitted models. With those
actual reference artifacts present:

```bash
python src/norman_shared_inference.py \
  --adata /workspace/datasets/sae-scrna/norman_pert_unique.h5ad \
  --out-dir "$NORMAN_REPRO_RUN/common"
python src/norman_pca_attribution.py \
  --adata /workspace/datasets/sae-scrna/norman_pert_unique.h5ad \
  --common-run "$NORMAN_REPRO_RUN/common" \
  --out-dir "$NORMAN_REPRO_RUN/attribution"
```

The common runner verifies the exact preprocessed input, actual activation hashes,
original PCA/NMF inference and the control split. These are conditional fixed-model
reanalyses. They are not new biological replicates or independent training holdouts.
For new SAE fits, use [the exact run configurations](NORMAN_CPU_RUNS.md), separate
paths and actual source/model/activation hashes; do not rerun finished jobs as a
routine readiness check. The warm-restarted NMF is a separate bounded numerical
sensitivity with its original result preserved.

Inspect all primary and low-power records and compare complete correction
families. JSON timestamps and paths need not match. CPU/GPU or library changes
can alter weights and scores; a common seed does not ensure byte equality.
The [verification table](CONTINUATION_VERIFICATION.md) distinguishes actual checks
from unavailable independent biological validation and publication requirements.
