# Abstract draft

Manuscript draft based on committed Gate-0 evidence. Later causal work is described
as exploratory context because its original result JSONs were not preserved.
This is not a claim of a new biological discovery.

## Working title

Measured hematopoietic gene-program modularity depends on decomposition:
a cross-species audit of sparse-autoencoder interpretation

## Abstract

Sparse autoencoders (SAEs) produce interpretable features from single-cell RNA
sequencing, but a plausible gene set does not establish that the inferred program
structure is intrinsic to the cells. This study tests an initial hypothesis that
erythroid commitment is more distributed than granulocyte commitment and examines
whether its verdict persists across decomposition methods.

A version-locked final analysis was applied to mouse Paul15 myeloid data
(2,730 cells, 2,012 genes) and human Setty 2019 CD34+ marrow data
(4,142 cells, 2,032 genes), after exploratory metric and marker development.
Ten SAE seeds per dataset were evaluated with a control-referenced decoder-loading
metric, an abundance-matched permutation null, a fixed effect floor, and
false-discovery-rate control. Structural comparisons included PCA, NMF,
co-expression regulons, and ICA; ICA used a smaller convergent rank. Robustness
analyses varied thresholds, sparsity, seed count, and gene-panel size.

The original hypothesis was supported in none of the ten SAE seeds in either
dataset. Median erythroid/granulocyte above-control program counts were 1.5/3
in mouse and 1/2 in human. Across five decomposition families and two datasets,
only mouse NMF supported the original decision rule. No single-marker omission
changed the SAE decision among 34 mouse and 36 human omissions. The original
hypothesis remained unsupported in the reported threshold, sparsity, seed-count,
and gene-panel sweeps, although some conditions produced ties rather than the
opposite direction.

These results demonstrate method dependence of this operational program-modularity
measure in two datasets. They do not establish a method-independent biological
partition or a causal lineage mechanism. Exploratory Perturb-seq analyses motivated
further tests of program-level modulation, but do not yet establish an SAE advantage
over foundation-model representations or a novel lineage-regulatory discovery.

## Evidence and scope

- Structural numbers trace to `data_g0/` and `data_g0_human/`: the frozen
  `v0b_outputs/` bundles, `baselines/`, `grn/`, `ica/`, `sensitivity/`, and
  `robustness/` summaries.
- One dataset per species limits generalization. Human CD34+ selection and the
  smaller ICA rank limit direct comparability. Marker panels and the chosen
  metric define the operational quantity measured.
- SHA locking records an analysis version; exploratory development preceded the
  final Gate-0 runs. Later causal amendments followed inspected outcomes.
- Original Component 2 JSONs were lost with ephemeral compute instances.
  [The causal results record](COMPONENT2_RESULTS.md) distinguishes historical
  reports from preserved machine-readable evidence. A publication extending this
  abstract to causal results requires archived reruns and appropriate uncertainty.
