# COMPONENT 2 PLAN — the causal head-to-head (GPU-blocked; specified here)

Component 2 is the centerpiece and the project's genuine novelty: test whether
**expression-space** SAE programs are causally grounded, and compare the grounding
rate head-to-head against Kendiukhov 2026's **embedding-space** SAE null of **6.2%**.
Gate 0 established that structural modularity is method-dependent; Component 2 asks
the question structure cannot answer: **are these programs causally real?**

This unit is NOT runnable in the current CPU container (needs GPU + the large
perturbation atlas). This document is the frozen, executable design so it can be run
unchanged on Colab/GPU. Nothing here is executed yet; every number is a target, not
a result. Do not cite anything in this file as a finding.

## Objective

For each SAE feature (candidate gene program), decide "causally grounded: yes/no" by
a pre-specified rule using CRISPRi perturbation data, then report the fraction of
features grounded. Compare that fraction to 6.2% (Kendiukhov, embedding-space) as the
null to beat. The scientific claim under test: expression-space SAE programs are
causally grounded at a rate meaningfully different from (hypothesis: higher than) the
embedding-space rate.

## Data

- **Replogle et al. 2022 genome-wide CRISPRi** (K562 and/or RPE1): each perturbation
  = knock-down of one gene, with the resulting single-cell / pseudobulk expression
  profile. This provides, per knocked-down gene g, a differential-expression vector
  delta_g over genes.
- The SAE decoder from Gate 0 (frozen), retrained on a matrix whose gene space
  overlaps the Replogle gene space (a Replogle-space SAE; the Paul15 SAE is not on
  the same genes, so Component 2 trains its own SAE on the perturbation-atlas control
  cells, reusing the exact frozen architecture and QC).

## Grounding criterion (to be pre-registered before running)

For SAE feature f with decoder column w_f (gene loadings):
1. Define the feature's **putative drivers** = top-T genes by |w_f| that are
   themselves perturbed in Replogle (T fixed in advance, e.g. 5).
2. Feature f is **causally grounded** iff knocking down at least one putative driver
   produces a differential-expression vector delta_g that aligns with the feature's
   loading w_f beyond a null: cosine(delta_g, w_f) exceeds the 95th percentile of an
   abundance/degree-matched null distribution of cosines (same matching philosophy as
   the Gate-0 abundance-matched null), with BH-FDR across features.
3. Grounding rate = fraction of features (with >=1 perturbed putative driver) that are
   grounded. Report per seed; require the standard >=60% seed majority for a feature
   to count as grounded.

Rationale: a real program should respond coherently when its own driver is silenced.
This mirrors the Gate-0 rigor (matched null, FDR, seed majority) so the causal test
inherits the same guardrails.

## The head-to-head

- **Expression-space arm (ours):** the grounding rate above, on the Replogle-space
  SAE.
- **Embedding-space arm (baseline):** reproduce Kendiukhov's protocol on a
  foundation-model (scGPT-scale) embedding SAE and recompute the grounding rate under
  the SAME criterion, so the 6.2% is measured, not just cited. If reproduction is out
  of scope, use 6.2% as the reported literature null and state that it was not
  re-measured under our exact criterion (an honesty caveat, not a silent assumption).

## Metric + decision

- Primary: expression-space grounding rate vs the embedding-space rate (6.2%),
  with a bootstrap CI over features and a two-proportion test.
- The claim is supported only if the expression-space rate's CI excludes 6.2% in the
  hypothesized direction; otherwise report the honest null/opposite.

## Compute

- GPU required: the perturbation atlas is large, and the embedding arm needs a
  single-cell foundation model (scGPT-scale). The SAE training itself is still small,
  but the data handling and the embedding arm are not CPU-feasible.
- Run on Colab (T4+) or NVIDIA hardware via Remote Control. See `docs/COMPUTE.md`.

## Runner outline (to fill in when on GPU)

```
# 1. acquire Replogle CRISPRi atlas (control cells + per-perturbation profiles)
# 2. preprocess control cells -> Replogle-space expression matrix (reuse
#    src/preprocess_* conventions: normalize+log1p+marker-aware HVG, MD5 record)
# 3. train the Gate-0 SAE on control cells (src/train_sae.py, l1 retuned to L0 band)
# 4. for each feature: putative drivers -> cosine(delta_g, w_f) vs matched null
#    (new src/component2_causal.py, mirroring the v3.1 null machinery)
# 5. grounding rate + head-to-head vs 6.2%; freeze the decision bundle
```

## Honest limitations (state up front)

- Cross-cell-type transfer: Replogle is K562/RPE1, not primary hematopoiesis. A
  program grounded there is grounded in that system; generalization to marrow is a
  separate claim.
- The grounding criterion (cosine + matched null) is one reasonable operationalization
  of "causal"; alternatives exist and should be sensitivity-tested (as Gate 0 did for
  its thresholds).
- If the embedding arm is not re-measured, the 6.2% comparison is literature-anchored,
  not apples-to-apples under our exact criterion.
