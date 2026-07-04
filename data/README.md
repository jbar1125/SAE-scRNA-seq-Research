# Data inputs (not in git)

The V0b pipeline reads these from `--data-dir` (default `./data`). They are not
committed (too large; private). Place them here locally, or point `--data-dir` at
the Colab Drive folder that holds them.

## Required files

| File | Shape / form | Size | Checksum |
|------|--------------|------|----------|
| `expression_matrix.npy` | 2730 x 2000 float32, scaled | ~20.8 MB | MD5 `60183a17983c8b977d036e0f3a58da61` (verified at load; OPTIONAL — only used for MD5 check + activation dynamics) |
| `gene_names.csv` | 2000 mouse HVGs (MGI symbols), first column | ~13 KB | — |
| `cell_metadata_palantir.csv` | V0a output, `prob_<term>_W<id>` branch-prob columns | ~485 KB | — |
| `sae_seed0.pt` … `sae_seed4.pt` | 5 mouse SAE checkpoints (encoder/decoder Linear 2000<->128) | ~2.06 MB each | — |

The MD5 of `expression_matrix.npy` is a project non-negotiable and is asserted at
load time. If you have only the small files (gene names, palantir CSV,
checkpoints), the core module-assignment + asymmetry decision still runs; the MD5
check and `submodule_dynamics.csv` are skipped with a warning.

## Where these live in Drive (as of this audit)

- `expression_matrix.npy`, `gene_names.csv`: folder id `1svskB4IltG8BEbnwZhPyLF2AEqX-fQhC`
- `cell_metadata_palantir.csv`, `gene_names.csv`, `sae_seed0-4.pt`: folder id `1uzrYWEw9O5q9yd9JJ68nf_m9CPC6ZNxV` (a general My-Drive dump, not a clean working dir)
- mouse checkpoints also in `SAE_stability` (folder id `103HfmZxJP78M4cfZ4elXdm1_bZKHnwRr`)
- human checkpoints `sae_human_seed0-2.pt` in `SAE_human` (folder id `1Dw-McC38RoYmSkZF3fAaZMIs0nJbN5DX`)
- `v0b_outputs/` (folder id `1qCFiUlfiuBi4ByQ7x4BCV9NsDBqXZZ8F`) currently contains ONLY a copy of `cell_metadata_palantir.csv` — there is NO frozen `module_assignments_v0b.csv`, confirming V0b has never been completed.

## Notebook trap (important)

Three V0b notebooks exist in Drive (`v0b_module_definitions.ipynb`, `_1`, `_2`).
The canonically-named, most-recently-modified one is **V1** (see PROJECT_AUDIT.md
for the five-point evidence). Do not run the notebooks. Use
`v0b_module_definitions.py` in the repo root, which is the canonical V2.
