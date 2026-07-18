# COMPONENT 2 — CAUSAL GROUNDING RESULTS (Arm B, expression-space)

> **STATUS: PRELIMINARY. NOT FROZEN.** Arm B (expression-space) only — the
> embedding-space Arm A that makes this a head-to-head is NOT yet run. Single cell
> line (K562). The grounded set is essential-gene-confounded (see below). Do NOT
> report these as a final finding or SHA-freeze them. Recorded 2026-07-18 as an
> honest interim checkpoint.

## Run configuration

- Data: Replogle 2022 K562 essential Perturb-seq (`pertpy.data.replogle_2022_k562_essential`), 310,385 cells.
- Genes: top-2000 HVG UNION perturbed-TF genes = 2,136.
- Scored perturbations: 162 TFs (`config/tf_lists/hs_hgnc_tfs.txt`, >= min_cells).
- SAE: TopK, trained on control (non-targeting) cells only, GPU (cuda).
- Metric: cross-fit Mann-Whitney grounding, spec SHA `b766d21c` (`config/causal_grounding_spec.json`).

## Seed stability (latent 2048, k 32)

| seed | grounded | rate |
|-----:|:--------:|-----:|
| 0 | 22/162 | 0.136 |
| 1 | 31/162 | 0.191 |
| 2 | 42/162 | 0.259 |
| 3 | 31/162 | 0.191 |
| 4 | 32/162 | 0.198 |
| **mean** | | **0.195 (sd 0.039, range 0.136-0.259)** |

All five seeds sit far above the label-shuffle null (~0.00). The finding "expression-
space SAE features show real perturbation-specific grounding" is SEED-ROBUST. The exact
rate is noisy (~20% coefficient of variation); the effect is not.

## Seed-stable grounded set (grounded in >= 3/5 seeds): 26

CNOT3, CSNK2B, CXXC1, DMAP1, E4F1, **GATA1**, GTF2A2, HINFP, HSF1, ILF2, MAX, NCBP2,
NELFB, POLD2, PSMD12, RFC2, RFC3, SFPQ, SMUG1, SRP9, TBP, TFDP1, THOC2, TIMELESS, ZMAT2,
ZNF335.

Honest categorization:
- **~24/26 are essential / general-transcription machinery / non-TF:** TBP, GTF2A2
  (basal apparatus), POLD2/RFC2/RFC3 (replication), SRP9 (signal recognition particle),
  PSMD12 (proteasome), THOC2/NCBP2/SFPQ (RNA processing/export), CSNK2B (kinase), SMUG1
  (base-excision repair), CNOT3, DMAP1, TIMELESS, ZMAT2, NELFB, ILF2, HINFP.
- **Notable real regulators that DO ground robustly:** **GATA1** — the master erythroid
  TF — grounds in >=3/5 seeds in an erythroleukemia line. This is a positive control
  PASSING: the metric catches a genuine lineage regulator. Plus MAX (MYC network), HSF1
  (heat-shock), E4F1, CXXC1, ZNF335 as defensible sequence-specific TFs.
- **Interpretation:** regulatory signal is PRESENT but BURIED under essential-gene
  transcriptional collapse. The three specificity fixes (S.1-S.3 in `UPGRADE_BACKLOG.md`)
  are required before any regulatory-specificity claim.

## Hyperparameter sensitivity (seed 0)

| config | grounded | rate |
|:------|:--------:|-----:|
| latent 512, k 32 | 64/162 | 0.395 |
| latent 1024, k 32 | 46/162 | 0.284 |
| latent 2048, k 32 | 22/162 | 0.136 |
| latent 2048, k 16 | 31/162 | 0.191 |
| latent 2048, k 64 | 58/162 | 0.358 |

- **Dictionary size (at k=32) is inversely related to the rate:** 512 -> 0.395, 1024 ->
  0.284, 2048 -> 0.136. Smaller/broader dictionaries ground MORE, consistent with the
  confound: broad features capture global-state axes that essential knockdowns suppress.
  A higher rate is therefore NOT "better" — it can mean less specific features.
- **k (at latent 2048) trends up with sparsity budget** (k64 = 0.358 highest), noisy at
  single seed.
- **The rate spans 0.136-0.395 from hyperparameters ALONE** (~3x). The absolute grounding
  rate is not a standalone quotable quantity.

## Honest conclusions

1. Grounding is real and seed-robust (all seeds >> null), but the ABSOLUTE rate is
   hyperparameter-dominated. The expression-vs-embedding head-to-head MUST fix (latent, k)
   identical across arms and report the rate across a grid, not a single number. Comparing
   the raw rate to Kendiukhov's 6.2% is invalid until Arm A is run at matched settings.
2. The grounded set is essential-gene-dominated, so the current rate largely reflects
   essential-knockdown transcriptional collapse, not regulatory specificity — BUT GATA1
   (and a few other real TFs) ground robustly, so the metric does detect genuine
   regulators; they are swamped, not absent.
3. Smaller dictionaries grounding more reinforces conclusion 2: the rate is inflated by
   broad, global-state features.

## Not done (the real next steps)

- **Arm A** (scGPT / Geneformer embedding-space SAE) at MATCHED hyperparameters — the
  actual head-to-head. This is the deliverable; the number above is not.
- The three specificity fixes S.1-S.3 (Lambert 2018 sequence-specific TF list; a
  program-coherence specificity filter that actually bites; an essential-gene control
  split via DepMap) — see `UPGRADE_BACKLOG.md`.
- Second cell line (RPE1), effect-size CIs, cross-fit K-repeat stability.

## Reproduce

Regenerable from the pinned code + seeds (raw JSONs were on an ephemeral GPU box and are
not committed; they regenerate exactly):

```bash
for s in 0 1 2 3 4; do
  python3 src/causal_pipeline.py --adata replogle.h5ad --pert-col gene \
    --control-value non-targeting --rep expression --latent 2048 --k 32 --seed $s \
    --n-hvg 2000 --tf-list config/tf_lists/hs_hgnc_tfs.txt --n-shuffle 50 \
    --out causal_out/expr_grounding_seed$s.json
done
```
