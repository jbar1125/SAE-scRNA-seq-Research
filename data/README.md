# Data inputs (not in git)

**CURRENT STATUS, 2026-10-18.** The files below describe the legacy Phase 1 input
bundle, not the executed Gate-0 datasets. The public source datasets and current
artifact registry are in `docs/PROJECT_MASTER_SUMMARY.md` sections 12 and 14.

| Executed Gate-0 arm | Directory | Matrix | Matrix MD5 |
|---|---|---|---|
| Mouse, Paul15 | `data_g0/` | 2,730 x 2,012 | `184ea5ddfca23db2bbb343b311688c45` |
| Human, Setty 2019 CD34+ marrow | `data_g0_human/` | 4,142 x 2,032 | `627961897d52a50c7a81f5bbf7e04097` |

Both arms use 512-latent, 10-seed Gate-0 SAEs. Small reports and result artifacts are
committed; expression matrices and checkpoints are ignored. The human run is not the
legacy 25,000-cell CELLxGENE sample. See the preprocessing reports in each directory
and `docs/REPRODUCIBILITY.md` for the supported Gate-0 workflow.

The causal experiments use separate Replogle/Norman inputs and all-cell TopK SAEs.
Norman preparation requires raw counts (`norman_2019_raw()`, not the processed
`norman_2019()` download). Historical causal raw JSONs from rented compute are not
archived in this checkout; reported summaries cannot substitute for that archive.

## Legacy Phase 1 bundle

The legacy V0b pipeline reads these from `--data-dir` (default `./data`). They are not
committed (too large; private). Place them here locally, or point `--data-dir` at
the Colab Drive folder that holds them.

## Required files

| File | Shape / form | Size | Checksum |
|------|--------------|------|----------|
| `expression_matrix.npy` | 2730 x 2000 float32, scaled | ~20.8 MB | MD5 `60183a17983c8b977d036e0f3a58da61` (verified at load; OPTIONAL — only used for MD5 check + activation dynamics) |
| `gene_names.csv` | 2000 mouse HVGs (MGI symbols), first column | ~13 KB | — |
| `cell_metadata_palantir.csv` | V0a output, `prob_<term>_W<id>` branch-prob columns | ~485 KB | — |
| `sae_seed0.pt` … `sae_seed4.pt` | 5 mouse SAE checkpoints (encoder/decoder Linear 2000<->128) | ~2.06 MB each | — |

The MD5 above applies only to the legacy matrix and is asserted by the legacy
loader. Do not substitute it for either Gate-0 hash. If you have only the small files (gene names, palantir CSV,
checkpoints), the core module-assignment + asymmetry decision still runs; the MD5
check and `submodule_dynamics.csv` are skipped with a warning.

## Where these live in Drive (as of this audit)

- `expression_matrix.npy`, `gene_names.csv`: folder id `1svskB4IltG8BEbnwZhPyLF2AEqX-fQhC`
- `cell_metadata_palantir.csv`, `gene_names.csv`, `sae_seed0-4.pt`: folder id `1uzrYWEw9O5q9yd9JJ68nf_m9CPC6ZNxV` (a general My-Drive dump, not a clean working dir)
- mouse checkpoints also in `SAE_stability` (folder id `103HfmZxJP78M4cfZ4elXdm1_bZKHnwRr`)
- human checkpoints `sae_human_seed0-2.pt` in `SAE_human` (folder id `1Dw-McC38RoYmSkZF3fAaZMIs0nJbN5DX`)
- At the Phase 1 audit, `v0b_outputs/` (folder id `1qCFiUlfiuBi4ByQ7x4BCV9NsDBqXZZ8F`) contained only `cell_metadata_palantir.csv`, with no frozen `module_assignments_v0b.csv`. Gate 0 was subsequently completed and frozen in the separate repository directories above; this old Drive observation does not describe its status.

## Notebook trap (important)

Three V0b notebooks exist in Drive (`v0b_module_definitions.ipynb`, `_1`, `_2`).
The canonically-named, most-recently-modified one is **V1** (see PROJECT_AUDIT.md
for the five-point evidence). Do not run the notebooks. Use
`src/v0b_module_definitions.py`, the legacy canonical V2. The current Gate-0 decision
uses `src/v0b_v3_1_decision.py`; do not mix the generations of inputs or metrics.
