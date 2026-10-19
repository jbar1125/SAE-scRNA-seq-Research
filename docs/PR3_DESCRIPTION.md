# Suggested PR title

Audit SAE program interpretation and preserve reproducible publication evidence

# Description

The historical program-modulation quantitative gates accept a shared response to
every perturbation as though it established perturbation-identity specificity.
An absent marker match also returned feature -1 and silently scored the last
activation column. This branch now reproduces those failures, leaves unmatched
programs unscored, and separately tests conditional identity-associated rank
heterogeneity with controls fixed, optional exchangeability strata, and Holm
correction across every planned program.

The branch includes the earlier structural hardening, causal pipeline and dated
errata, plus the continuation audit. Documentation now distinguishes preserved
Gate-0 evidence from historical causal outputs that were not archived. It corrects
all-cell training, the wrong-background essentiality figure, the PCA/NMF-to-SAE
ladder interpretation, and unrun embedding comparisons. Frozen specifications
and Gate-0 decision bundles remain unchanged.

## New evidence and checks

- Four synthetic worlds, 20 seeds each, 99 permutations: historical quantitative
  gates pass 20/20 generic one-axis worlds and 20/20 absent-marker worlds. The
  corrected identity test rejects 0/20 generic worlds and detects 20/20 distinct-
  responder worlds. Manual biological coherence was not scored.
- Separate deployed-default null calibration: nine stipulated matched features,
  999 permutations, 200 independent simulations; 6/200 familywise Holm rejections,
  3.0% (exact 95% interval 1.11–6.42%). Finite simulated calibration is not proof
  for arbitrary experimental designs.
- Raw simulation outputs, executed generator, source hashes, versions, parameters
  and summaries are retained under `data_program_modulation_audit/`.
- All original standalone suites and new counterexample/specificity suites,
  syntax checks, frozen hashes, and synthetic actual-training CLI checks have
  been exercised. See the continuation audit and reproducibility guide.
- Corrected historical figures include descriptive uncertainty and are
  deterministic across consecutive generations.

## Fresh real-data continuation

The public Norman raw file is now downloaded, fingerprinted and prepared. A
full-data external-marker baseline uses 5,926 discovery and 5,929 heldout controls,
eight gemgroup strata, all nine panels and all 105 sufficiently sampled single
perturbations (BH over 945 comparisons). Heme associations number 37; 29 survive
target omission, discovery-defined RNA-depth matching and all four negative-control
references. The reference sensitivity corrects all 3,780 comparisons jointly.

KLF1's TF-panel AUC falls from 0.615 to 0.425 after its own transcript is omitted,
while its heme and membrane associations persist. CBL and PTPN12 survive all
sensitivities; UBASH3B falls below the fixed effect floor against one reference.
These are conditional marker-score associations, not biological novelty or an
SAE advantage. Complete JSONs, source/data hashes, a 945-row table and reproducible
figures are archived in `data_norman_continuation/`.

The declared ordinary PCA comparison also completed: 64 components, both signed
orientations, discovery-only annotation and the same inference controls, with
joint BH over 945 tests. Ery_TF and Cycling matched (28 and 52 associations);
the other seven panels were unmatched and remain unscored. Different model
capacities and score definitions prevent a superiority claim from these counts.

All 15 suites passed together with no skips; the completed-result archiver and
direct-PCA suites passed separately (17 distinct suites exercised). Full-size CPU seed 0 is
training toward 200 epochs with resumable model/optimizer/RNG checkpoints; seeds
1 and 2 are queued sequentially. A separate copied epoch-20 checkpoint completed
full-data inference; its complete nine-panel outputs are archived as an interim
readiness check, not a completed primary seed. A tested watcher will validate and
publish compact outputs after actual full-seed completion.

Crossref verified nine bibliographic references. Targeted primary-source review
documents known CBL/erythroid and UBASH3B/erythroleukemia context and distinguishes
model-internal circuits from cellular causality. Novelty remains unestablished;
the original Norman full text and supplements still require review.

## Scope and outstanding work

No new biological discovery is established. GitHub API access recovered and this
PR is updated through REST; no GPU is attached and the user chose CPU execution.
The manuscript now includes preserved structural evidence, the simulation audit
and the fresh marker/PCA baselines. Full SAE completion, remaining primary-text review, confirmed
authorship/disclosures, independent validation and submission remain outstanding.

Conditional cell-level identity heterogeneity does not establish pathway
specificity, biological replication, novelty, or SAE superiority. Full Norman
runs, direct NMF and appropriate guide/batch validation remain separate experiments.
