# SAE / scRNA-seq Hematopoiesis Research

A fully computational research project (Jacob Barzideh) testing whether sparse
autoencoders (SAEs) can discover **real** gene programs in blood-cell development,
and building a rigorous framework to tell a genuine biological program from a
decomposition artifact.

**One line:** two acts. (1) A rigor foundation — a pre-registered, cross-species,
cross-method demonstration that SAE gene-program modularity is method-dependent, so no
single-method structural claim is trustworthy. (2) The centerpiece it motivates — a
causal head-to-head asking whether SAEs trained on gene EXPRESSION recover the
regulatory logic that single-cell foundation models discard (only ~6% causally
grounded), tested on Replogle CRISPRi. The metric + pipeline for (2) are built and
CPU-verified; the GPU run is the open step. See `docs/ELEVATION_PLAN.md`.

Start here → **[`LAB_NOTEBOOK.md`](./LAB_NOTEBOOK.md)** (dated day-to-day log) and
**[`docs/COMPONENT0_RESULTS.md`](./docs/COMPONENT0_RESULTS.md)** (the executed
results + frozen SHAs). For agents: **[`CLAUDE.md`](./CLAUDE.md)** is the operating
manual, read it first.

---

## Status at a glance

| Item | State | Where |
|------|-------|-------|
| Gate 0 executed + frozen (mouse + human) | ✅ done, in `main` | `docs/COMPONENT0_RESULTS.md` |
| Robustness battery R1–R4 (both species) | ✅ done | `docs/COMPONENT0_HARDENING.md` |
| 5 decomposition families × 2 species | ✅ done | `docs/COMPONENT0_HARDENING.md` |
| TRRUST regulon cross-reference | ✅ done | `data_g0*/trrust/` |
| Result figures (deterministic) | ✅ done | `docs/figures/`, `src/make_figures.py` |
| Competition abstract + visual brief | ✅ done | `docs/ABSTRACT.md`, `docs/web/result_brief.html` |
| Reproducibility verification | ✅ passes | `docs/REPRODUCIBILITY.md` |
| **Causal head-to-head — metric + pipeline** | ✅ **built, CPU-verified** | `src/causal_grounding.py`, `src/causal_pipeline.py`, `tests/test_causal_grounding.py` |
| Causal metric pre-registered (frozen spec) | ✅ v6 (effect-size-controlled; M1) | `config/causal_grounding_spec.json` (SHA `81b44e04`) |
| Causal head-to-head — Arm B (expression) run | 🔴 **honest negative** (grounds housekeeping/effect-size, NOT lineage regulators; no master TF robust) | `docs/COMPONENT2_RESULTS.md` |
| Causal head-to-head — Arm A (embedding) + specificity fixes | 🟡 next | `docs/UPGRADE_BACKLOG.md` (S.1-S.3, Tier 1.1) |
| Causally-supervised SAE (a method idea) | ❌ prototyped, did NOT validate — shelved | `docs/ELEVATION_PLAN.md` §4 |
| Second human dataset / full SCENIC / mouse TRRUST | ⛔ blocked (no reachable data/DB here) | logged in `LAB_NOTEBOOK.md` |

**Act 1 finding (done):** the original "asymmetric modularity" claim is NOT SUPPORTED
in either species; the modularity verdict is method-dependent (survives in 1 of 10
method × species cells); every stress test confirms it. Structural, not yet causal.

**Act 2 (the winnable result, built and ready):** does an expression-space SAE beat the
~6-10% causal-grounding ceiling of foundation-model SAEs on Replogle CRISPRi? Metric and
pipeline are built and verified on synthetic (recovers planted regulators, FDR-calibrated,
defeats the triviality trap); the GPU run is the open step.

Follow-up hardening lives on PR #3 (draft); Gate 0 is already merged to `main`.

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
| `DIFFERENTIATION.md` | Verified competitor landscape + differentiation table + rehearsed novelty answer (all competitors are embedding-space; this is expression-space). |
| `PREREGISTRATION.md` | The frozen analysis spec (what was locked before running). |
| `COMPONENT0_STATUS.md` | Gate-0 state + the Colab/GPU reproduction runner. |
| `COMPONENT0_RESULTS.md` | The EXECUTED Gate-0 numbers (mouse + human) + frozen SHAs. |
| `COMPONENT2_RESULTS.md` | Causal Arm-B (expression) preliminary run: seed table, sensitivity sweep, the essential-gene confound (NOT frozen). |
| `COMPONENT2_NEXT_DIRECTIONS.md` | Results-driven execution map: the v4 re-run, the head-to-head grid redesign, and the exhaustive prioritized direction tree. |
| `COMPONENT0_HARDENING.md` | Robustness battery (R1-R4 + TRRUST + GRN + ICA) stress-testing the frozen verdict. |
| `ABSTRACT.md` | One-page competition/mentor summary (draft), traced to committed numbers. |
| `REPRODUCIBILITY.md` | Scripted end-to-end verification + expected hashes/outputs. |
| `ELEVATION_PLAN.md` | **The pivot**: how to turn the null into a positive, causal, benchmark-beating result. |
| `UPGRADE_BACKLOG.md` | Exhaustive, prioritized backlog of upgrades + redirections (tiers 0-6 + big reframes). |
| `IMPACT_STRATEGY.md` | How to make the project impressive/influential AFTER the causal negative: the auditor reframe, the effect-size-residual metric, the benchmark, the 3 project shapes. |
| `RUNPOD_EXECUTION.md` | Turnkey GPU steps for the expression-vs-embedding causal head-to-head. |
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
