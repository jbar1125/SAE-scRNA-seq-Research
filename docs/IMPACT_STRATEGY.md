# IMPACT STRATEGY: evidence and publication priorities

**CURRENT STATUS, 2026-10-18.** This replaces the 2026-07-18 strategy whose central
effect-size explanation was contradicted by the next real-data experiment. Historical
proposals and their outcomes are recorded below; they are not findings. Read
[PROJECT_MASTER_SUMMARY.md](PROJECT_MASTER_SUMMARY.md), especially sections 8-10 and 18,
for the experiment history and [COMPONENT2_RESULTS.md](COMPONENT2_RESULTS.md) for the
reported Replogle results.

## What the evidence supports

- Gate 0: the original asymmetric-modularity claim fails in both species under the
  frozen SAE analysis. Across five decomposition families and two species, only mouse
  NMF supports it. The defensible structural result is **method-dependence**, not an
  invariant reversal of erythroid and granulocyte organization.
- Replogle CRISPRi: the reported grounded sets are mostly housekeeping/chromatin
  programs. The effect-size control retained all 22/162 seed-0 hits; the measured
  effect-size/suppression correlation was rho = 0.28. This does not support the former
  claim that grounding simply tracks perturbation effect size.
- The reported representation ladder compares expression -> SAE, PCA-128 -> SAE and
  NMF-128 -> SAE. Its similar rates do **not** establish whether an SAE adds value over
  direct PCA/NMF factors. Foundation-model embeddings were never run.
- Norman CRISPRa: the aggregate rate did not exceed its perturbation-label null
  (reported p = 0.902). A seed-stable erythroid program is an observation, not yet a
  validated program-specific causal result. GATA1 is absent from that library.
- Regulon recovery: 8/42 TF-by-seed cases exceed the random-set chance mean and 0/42
  exceed its 95th percentile. The existing runs do not establish regulon recovery.

**Evidence limit:** the causal-run raw JSON archive from rented compute is missing;
the figures and summaries preserve reported results, not an independently auditable
run archive. Recover or regenerate that archive before presenting exact causal
numbers as publication-ready evidence. The Gate-0 small result artifacts are committed.

## The thesis that survives

> An audit of sparse-autoencoder interpretability in single-cell expression data:
> lineage modularity depends on the decomposition, and biological annotation alone
> does not establish perturbation-specific regulatory meaning. Explicit nulls and
> provenance checks expose several ways a promising metric can mislead.

This is a bounded methods claim about the tested datasets and implementations.
It does not establish that all SAEs, all expression representations, or all
foundation-model features fail to encode regulation.

## Which July proposals survived testing

| Original direction | Observed outcome or current status |
|---|---|
| M1: diagnose effect-size confounding and add a residual gate | The implemented gate is an AND-filter: it can remove hits, never rescue them. It removed 0/22 Replogle seed-0 hits. The proposed effect-size explanation was not supported. |
| M2: representation ladder and embedding head-to-head | PCA/NMF input transformations were run before SAE training. The foundation-model arm remains unrun; the literature's 6.2% uses a different metric and panel. |
| M3: artifact taxonomy | Supported as a documented audit: method-dependence, wrong enrichment background, failed column specificity, dominant-axis/gate saturation, hub-target selection, TopK competition and training-mode mismatch. Separate demonstrated diagnostics from proposed mechanisms. |
| M4: GATA1 boundary case | GATA1 was a Replogle v3 hit but lost seed-majority support under v4. It is not a validated universal positive control and was not perturbed in Norman. |
| M5: standalone benchmark | A possible future deliverable. Packaging does not repair missing run provenance or validate the metric. |
| M6: preprint | Appropriate once claims, archived results, methods and citations are auditable. Publication or competition success cannot be promised. |
| M7: communication | Use the supported thesis above; retire the effect-size scatter as a centerpiece and the unexecuted expression-vs-embedding claim. |

The earlier "~10x essential enrichment" calculation used all 1,839 TFs rather than
the 162 tested. In the reported v4 split, core-essential-tagged TFs ground at 6/51
(11.8%) and the other tested TFs at 14/111 (12.6%). The latter are **not** a validated
nonessential reference set. S.2 column specificity retained housekeeping hits and
reduced GATA1 support; it was rolled back. Neither is a solved confound contribution.

## Publication priorities

1. **Make the existing record auditable.** Recover causal JSONs, or rerun the exact
   reported configurations with dataset checksums, code/spec hashes, environment,
   seeds, training-cell membership and complete outputs. Distinguish new runs from
   reconstructed historical reports.
2. **Test the erythroid program against its declared nulls.** Follow the program-
   modulation spec and report every seed, the feature-background comparison and
   specificity controls. Archive negative or mixed outcomes too. The all-cell SAE is
   the prior executed mode; control-only training failed the synthetic self-test and
   is a sensitivity analysis, not a retrospectively implemented safeguard.
3. **Write the audit alongside the experiment.** Gate 0 and the metric failure modes
   remain reportable whichever way program modulation goes. Keep causal attribution
   proportional to the experiment and distinguish feature selection from testing.
4. **Extend only after the current claims are secure.** A matched embedding arm,
   direct PCA/NMF causal baselines, held-out training and independent replication
   answer different questions; none has been established by the existing ladder.

## A usable one-sentence pitch

"I audited whether sparse-autoencoder gene programs in single-cell data survive
structural and perturbation-based tests, found that modularity depends on the method,
and documented the controls needed before calling an annotated feature regulatory."
