# Abstract / one-page results summary

Draft for competition (ISEF / Regeneron STS) and mentor review. Written from the
executed, version-controlled results only; every number traces to a committed
artifact. Honest by design: it leads with the negative and labels every limit.

## Title (working)

Gene-program modularity in hematopoiesis is decomposition-dependent: a
pre-registered, cross-species stress test of sparse-autoencoder programs.

## Background

Sparse autoencoders (SAEs) are increasingly used to extract interpretable "gene
programs" from single-cell RNA-seq by disentangling superimposed signals into sparse,
ideally monosemantic features. A prevailing narrative in hematopoiesis is "asymmetric
modularity": that granulocyte commitment is a unified program while erythroid
commitment is distributed across many. The deeper and rarely-tested question is
whether such program structure is a property of the biology or an artifact of the
chosen decomposition.

## Question

(1) Does the asymmetric-modularity claim survive a rigorous, pre-registered test on
sparsity-corrected, multi-seed SAEs? (2) Is any modularity verdict invariant to the
decomposition method, or is it method-dependent?

## Methods

- **Data.** Mouse Paul15 myeloid atlas (2730 x 2012) and human Setty 2019 CD34+ marrow
  (4142 x 2032), both hash-verified; marker-aware preprocessing with corrected,
  lineage-complete panels (mouse; human orthologs).
- **Model.** Overcomplete SAEs (512 latents), 10 seeds/species, sparsity tuned so mean
  L0 lands in a pre-set 20-50 quality band (mouse 35.8, human 32.5).
- **Metric (pre-registered, SHA-frozen before running).** A control-referenced
  modularity test: a submodule counts as a real program only if its decoder-mass
  concentration exceeds a Progenitor-control baseline by an effect-size floor AND beats
  an abundance-matched permutation null (BH-FDR). Decision rule: >= 60% of seeds.
- **Rigor.** PCA/NMF/GRN/ICA baselines under the identical metric; leave-one-marker-out
  sensitivity; decision-threshold, L0-band, seed-count (to 20), and HVG-count sweeps;
  mouse->human replication; TRRUST regulatory-database cross-reference.

## Results

- **The original asymmetric-modularity claim is NOT SUPPORTED** in mouse (0/10 seeds)
  or human (0/10). Where a direction exists, the data show the REVERSE of the headline:
  granulocyte carries more real above-control programs than erythroid.
- **The verdict is method-dependent.** Across five decomposition families and two
  species, the original claim survives in exactly **1 of 10** method x species
  combinations (mouse NMF, the most parts-based/artifact-prone). PCA resolves little;
  ICA (at its convergent rank) and the SAE agree on the granulocyte direction; NMF
  disagrees. No unsupervised decomposition yields a method-invariant program structure.
- **The negative is robust.** It does not flip under any decision-threshold choice
  (0/12 configs), across the L0 band, at 20 seeds, or across HVG counts (1000-3000), in
  either species. No single marker drives it (0/34 mouse, 0/36 human leave-one-out).
- **The SAE is not hallucinating gene sets.** A TRRUST cross-reference shows the SAE
  feature that captures a lineage's transcription factors also captures their curated
  targets (human erythroid: 100% of seeds, p = 6.4e-5).
- **Null calibration holds:** the Progenitor control is non-significant in both species.

## Interpretation and contribution

The contribution is not "SAEs on single-cell data" (now a published wave, e.g.
Kendiukhov 2026, CytoSAE). It is a pre-registered **artifact-vs-signal framework** that
demonstrates, across five decompositions and two species, that gene-program modularity
is co-determined by the chosen method. This reframes a class of single-method "gene
program" claims as method-contingent, and motivates a causal test as the true arbiter.

## Limitations

Structural (decoder-mass concentration), not causal. n=10 seeds, one dataset per
species; the human erythroid globin program is under-represented by CD34+ selection.
ICA is reported at rank 50 (its convergent regime), not the 512 of the other families.
Mouse TRRUST and full motif-pruned SCENIC were not reachable in the run environment.

## Next (the centerpiece)

A causal head-to-head: test whether expression-space SAE programs are validated by
Replogle CRISPRi perturbation at a higher rate than the ~6.2% reported for
embedding-space SAE features. This is specified in `COMPONENT2_PLAN.md` and requires a
GPU environment.

---

Provenance: numbers from `COMPONENT0_RESULTS.md` (frozen bundles mouse `172861...`,
human `cc8159c8...`) and `COMPONENT0_HARDENING.md`; figures in `docs/figures/`.
