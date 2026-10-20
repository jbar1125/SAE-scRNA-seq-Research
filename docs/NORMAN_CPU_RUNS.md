# Running and resuming Norman on CPU

**2026-10-20 status:** all three primary seeds and the separate control-only seed
completed 200 epochs. The sequence and archive workers exited after completion.
The six-model common-inference and signed-weight diagnostics also completed;
their compact full evidence is archived. Original PCA scores were successfully
recovered, with the complete old inference reproduced before their use.
Read [the completed-fit interpretation](NORMAN_COMPLETED_MODELS.md); do not rerun
these commands merely because no trainer is alive. Commands below preserve the
executed configurations and support recovery if artifacts are independently lost.

Use the existing isolated checkout. Do not create a Git worktree unless requested.
The user selected this CPU environment. Cached raw/prepared inputs are under
`/workspace/datasets/sae-scrna/`; models and activations are under
`/workspace/research-runs/norman-sae-primary/`. Large files remain outside Git.
Keep the archived input manifests and preserve original data files.

```bash
cd /workspace/SAE-scRNA-seq-Research
source /workspace/.venvs/sae-scrna/bin/activate
export XDG_CACHE_HOME=/workspace/.cache
export MPLCONFIGDIR=/workspace/.cache/matplotlib
export NUMBA_CACHE_DIR=/workspace/.cache/numba
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
```

Before launching, inspect both training processes and any
`src/continue_norman_cpu.sh` sequence worker to avoid duplicate training. A worker
was launched behind the existing seed-0 process and will validate/continue seeds
1 and 2 automatically. If the environment is restored with no live worker or
trainer, `bash src/continue_norman_cpu.sh` resumes the incomplete sequence. Do not
start a second worker while the first is alive. Its full parameters and resume
selection were dry-run checked; the underlying checkpointed CLI was executed on
synthetic data before the full launch. The primary configuration is recorded in
`config/norman_cpu_execution_plan.json`. Seed 0 was launched with:

```bash
python -u src/program_modulation.py \
  --adata /workspace/datasets/sae-scrna/norman_pert_unique.h5ad \
  --pert-col perturbation --control-value control --train-on all \
  --latent 2048 --k 32 --n-hvg 4000 --epochs 200 --batch 1024 --seed 0 \
  --n-shuffle 0 --n-identity-shuffle 999 --strata-col gemgroup \
  --checkpoint-every 5 \
  --checkpoint /workspace/research-runs/norman-sae-primary/seed0.pt \
  --activations-out /workspace/research-runs/norman-sae-primary/seed0_activations.npy \
  --out /workspace/research-runs/norman-sae-primary/seed0.json
```

If that process exits before completion, run the same command with `--resume`.
`--epochs 200` means total epochs, not 200 additional epochs. The checkpoint
rejects different data bytes, shapes, dtypes, training parameters or Torch version.
The resumed command must use a fresh JSON/activation path if a partial prior
analysis already produced one. Never overwrite archived analysis results.
Seeds 1 and 2 use their corresponding seed number and distinct artifact paths.

The command was exercised with a synthetic AnnData file, real two-epoch training,
saved activations and checkpoint provenance before the full-data launch. The
checkpoint suite also verifies exact five-epoch weights and RNG equality after
two epochs plus a three-epoch resume. Full 200-epoch real-data completion is now
verified from actual checkpoints, complete outputs and model/activation hashes.

The old full-label shuffle is disabled because it does not calibrate specificity.
The separate identity test uses 999 within-gemgroup permutations and Holm over all
nine planned panels. Claims across seeds need the complete 27-hypothesis family,
and agreement across seeds is not biological replication. Match failures and
negative outcomes must remain in the report.

A separate `src/archive_norman_cpu_results.py --watch --publish` worker completed.
It validates actual completion, configuration, the input fingerprint, the entire
nine-panel family, 999-permutation execution, and large model/activation hashes
before copying compact outputs. After all three seeds complete, it applies Holm
to all 27 identity tests and generates `docs/NORMAN_SAE_RESULTS.md`. Its scoped Git
commits preserve unrelated staged work and push the user-authorized branch.
Inspect `archive-worker.log` and the existing worker before starting another copy.
If absent, start it from the activated environment with the same command. A copied
result interrupted before publication is recovered on the next pass; differing
existing outputs are preserved for review. Large weights remain outside Git.
This automation does not replace scientific interpretation or biological validation.

The completed marker-baseline command and its exploratory sensitivity command are:

```bash
python src/marker_program_baseline.py \
  --adata /workspace/datasets/sae-scrna/norman_pert.h5ad \
  --out /tmp/norman-marker-new-result.json \
  --strata-col gemgroup --duplicate-genes sum --n-identity-shuffle 999
python src/marker_baseline_sensitivity.py \
  --adata /workspace/datasets/sae-scrna/norman_pert.h5ad \
  --primary /tmp/norman-marker-new-result.json \
  --out /tmp/norman-marker-new-sensitivity.json
```

Equivalent commands completed on the full dataset. Choose fresh output paths.
Their plans, raw outputs, source fingerprints and comparison table are archived
in `data_norman_continuation/`; results are summarized in
[NORMAN_CONTINUATION_RESULTS.md](NORMAN_CONTINUATION_RESULTS.md).

## Separate control-only sensitivity

The [concrete sensitivity plan](../config/norman_control_sae_sensitivity_plan.json)
retains the architecture and 200-epoch target but trains seed 0 on 11,855 controls,
then encodes all 111,445 cells. One CPU thread is used, and it is queued behind
the ordinary NMF process. This does not replace the three all-cell primary seeds
or complete a three-seed control sensitivity. Full-data preprocessing still
defines the gene panel; it is not an independent preprocessing holdout.

Artifacts are under `/workspace/research-runs/norman-sae-control-sensitivity/`.
Check its waiting worker/trainer before starting another copy. Its command uses
the primary CLI flags with these explicit differences:

```bash
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python -u src/program_modulation.py \
  --adata /workspace/datasets/sae-scrna/norman_pert_unique.h5ad \
  --pert-col perturbation --control-value control --train-on control \
  --latent 2048 --k 32 --n-hvg 4000 --epochs 200 --batch 1024 --seed 0 \
  --n-shuffle 0 --n-identity-shuffle 999 --strata-col gemgroup \
  --checkpoint-every 5 \
  --checkpoint /workspace/research-runs/norman-sae-control-sensitivity/seed0.pt \
  --activations-out /workspace/research-runs/norman-sae-control-sensitivity/seed0_activations.npy \
  --out /workspace/research-runs/norman-sae-control-sensitivity/seed0.json
```

Use the same one-thread settings and `--resume` for a partial checkpoint, after
checking existing output artifacts. The original synthetic control-only test
found 0/8 planted regulators; that is a limitation to report alongside real-data
results, not a reason to suppress this arm or call it a validated replacement.
