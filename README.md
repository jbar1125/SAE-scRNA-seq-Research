# SAE / scRNA-seq Hematopoiesis Research

A fully computational research project (Jacob Barzideh) testing whether sparse
autoencoders (SAEs) can discover **real** gene programs in blood-cell development,
and building a rigorous framework to tell a genuine biological program from a
decomposition artifact.

**One line:** the original "asymmetric modularity" headline is not supported; what
replaces it is a pre-registered, cross-species, cross-method demonstration that
gene-program modularity is method-dependent, motivating a causal test as the real
arbiter.

Start here → **[`LAB_NOTEBOOK.md`](./LAB_NOTEBOOK.md)** (dated day-to-day log) and
**[`docs/COMPONENT0_RESULTS.md`](./docs/COMPONENT0_RESULTS.md)** (the executed
results + frozen SHAs). For agents: **[`CLAUDE.md`](./CLAUDE.md)** is the operating
manual, read it first.

---

## Repository map

```
README.md            <- you are here (project map)
CLAUDE.md            <- operating manual for Claude Code (rules, status)
LAB_NOTEBOOK.md      <- dated research log (thoughts + what was done each session)
requirements.txt     <- pinned environment

src/                 <- all executable code
  v0b_module_definitions.py   canonical v2 module logic (markers, coverage, loaders,
                              V0B_MARKERS_JSON override for the human arm)
  v0b_v3_loading.py           continuous decoder-loading strength (no gate)
  v0b_v3_1_decision.py        THE pre-registered metric (control-referenced,
                              abundance-matched null, effect-size floor, 60% rule)
  train_sae.py                overcomplete SAE trainer (L0 QC band reporting)
  preprocess_paul15.py        mouse marker-aware preprocessing (+ h5 fallback loader)
  preprocess_human_marrow.py  human CELLxGENE preprocessing (ortholog markers)
  baselines_nmf_pca.py        PCA/NMF participation ratio + same v3.1 modularity
  marker_sensitivity.py       leave-one-marker-out robustness

tests/               <- synthetic logic tests (no data/torch needed); add src/ to path
config/              <- frozen specs + marker panels
  preregistration_spec.json / .sha256   frozen metric spec (SHA fc342829)
  human_markers.json                    1:1 human orthologs of the mouse submodules
docs/                <- narrative + reference documents (see index below)

data/                <- manifest for the large Drive inputs (gitignored)
data_g0/             <- EXECUTED mouse run (matrices/checkpoints gitignored;
                        small result artifacts committed as evidence)
data_g0_human/       <- EXECUTED human run (same convention)
```

## Documentation index (`docs/`)

| Doc | What it is |
|-----|-----------|
| `PROJECT_HANDOFF.md` | Phase-1-closeout snapshot (historical record). |
| `PROJECT_AUDIT.md` | Authoritative corrections + rigor log; sections J/K are current findings. |
| `STRATEGY_AND_POSITIONING.md` | Field situating, novelty audit vs Kendiukhov et al., the reframe. |
| `PREREGISTRATION.md` | The frozen analysis spec (what was locked before running). |
| `COMPONENT0_STATUS.md` | Gate-0 state + the Colab/GPU reproduction runner. |
| `COMPONENT0_RESULTS.md` | The EXECUTED Gate-0 numbers (mouse + human) + frozen SHAs. |
| `COMPONENT0_HARDENING.md` | Robustness battery (R1-R4 + TRRUST + GRN + ICA) stress-testing the frozen verdict. |
| `ABSTRACT.md` | One-page competition/mentor summary (draft), traced to committed numbers. |
| `phase2_research_plan_v6.md` | Reframed plan; the causal head-to-head is the spine. |
| `VERSION_HISTORY.md` | Metric/pipeline version lineage (v2 -> v3 -> v3.1). |
| `COMPUTE.md` | What runs on CPU vs GPU and why (incl. why the AMD 5700 XT isn't practical). |
| `COMPONENT2_PLAN.md` | Executable design for the causal CRISPRi head-to-head vs the 6.2% null (GPU-blocked; specified, not run). |

## How to run (from repo root)

```bash
pip install -r requirements.txt

# logic tests (no data needed)
python3 tests/test_v0b_logic.py
python3 tests/test_v0b_v3_1_logic.py
python3 tests/test_baselines_logic.py

# mouse pipeline (needs Paul15 -> data_g0/expression_matrix.npy)
python3 src/preprocess_paul15.py --n-hvg 2000 --out-dir ./data_g0
python3 src/train_sae.py --matrix ./data_g0/expression_matrix.npy --latent-dim 512 --l1 0.8 --seeds 10 --epochs 300 --out-dir ./data_g0
python3 src/v0b_v3_1_decision.py --data-dir ./data_g0 --output-dir ./data_g0/v0b_outputs
python3 src/baselines_nmf_pca.py --data-dir ./data_g0 --rank 512 --seeds 10 --out ./data_g0/baselines
python3 src/marker_sensitivity.py --data-dir ./data_g0 --out ./data_g0/sensitivity

# human replication arm (needs Setty marrow -> data_g0_human/)
python3 src/preprocess_human_marrow.py --out-dir ./data_g0_human
V0B_MARKERS_JSON=config/human_markers.json python3 src/train_sae.py --matrix ./data_g0_human/expression_matrix.npy --latent-dim 512 --l1 0.8 --seeds 10 --epochs 300 --out-dir ./data_g0_human
V0B_MARKERS_JSON=config/human_markers.json python3 src/v0b_v3_1_decision.py --data-dir ./data_g0_human --output-dir ./data_g0_human/v0b_outputs
```

Runs are CPU-friendly (~50-75 s/seed; a full 10-seed run is ~8-12 min).

## Current status

**Gate 0 (Component 0): COMPLETE and frozen** (mouse + human). Original claim NOT
SUPPORTED in both species; the modularity verdict is method-dependent (PCA/NMF/SAE
disagree); marker leave-one-out is robust (0/34 mouse, 0/36 human); nulls calibrated.
Frozen bundles: mouse `172861...`, human `cc8159c8...`.

**Next: Component 2** — Replogle CRISPRi causal head-to-head vs Kendiukhov's 6.2%
embedding-space null (the centerpiece). Optional cheap hardening first: SCENIC as a
5th baseline, TRRUST/ChIP-Atlas cross-reference.

## Data provenance

- Mouse: Paul15 myeloid atlas, obtained via git-LFS mirror, hash-verified.
- Human: Setty 2019 CD34+ bone-marrow (distributed with Palantir), sha256
  `be4f8603...`, hash-verified.
- Large binaries (matrices, checkpoints, raw h5ad) are gitignored and regenerable;
  the committed small JSON/CSV artifacts are the reproducible evidence.
