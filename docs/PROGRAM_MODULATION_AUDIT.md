# Program modulation: continuation audit, 2026-10-18

The historical quantitative program-modulation gates do not establish
perturbation-identity specificity. A shared response to every perturbation can pass
both gates when it is concentrated in one feature. An independent implementation
audit also found that zero marker overlap returned feature `-1`, after which the
caller scored the last activation column as the named program.

This is a reproducible methods finding. It is not a biological discovery, an
estimate of performance on Norman data, or a demonstration that SAEs cannot learn
useful programs. The historical specification and all four frozen specifications
and both Gate-0 bundles remain unchanged.

## Executed counterexamples

Historical source was loaded from full Git commit
`a08db6d5ef4440adb9f1e1613e96856ab6515533`. Each simulated world has 500 controls,
20 perturbations with 60 cells each, 40 directly generated activation features,
205 genes and five external program markers. These are downstream-metric tests;
no SAE was trained for this benchmark. World construction was fixed before the
20-seed run, after inspection of an initial counterexample. Worlds share baseline
draws within a seed; comparisons are paired.

The historical gates are full-label shuffle p < 0.05 and feature-background
percentile at least 90. The specification's manual biological-coherence criterion
was not operationalized or scored.

| World | Historical quantitative gates pass | Current conditional identity test rejects |
|---|---:|---:|
| Four distinct responders among 20 perturbations | 20/20 | 20/20 |
| Every perturbation identically activates one feature | 20/20 | 0/20 |
| No perturbation effect | 0/20 | 0/20 |
| No expressed marker match; unrelated last feature responds | 20/20 | 0/20; all unmatched |

This benchmark uses one planned program and 99 permutations. It does not evaluate
the deployed nine-program multiplicity family. Exact 95% binomial intervals for
the simulation frequencies are 83.16–100% for 20/20 and 0–16.84% for 0/20. These
wide intervals must not be presented as proof of perfect sensitivity or zero
false positives. All four planted responders were recovered in every specific
world; two seeds also produced one additional exploratory modulator call.

Full outputs, the executed generator, source hashes, summaries and reanalysis
are retained in [the evidence archive](../data_program_modulation_audit/README.md).

## Corrected scope and implementation

`program_modulation.py` now leaves absent or insufficiently enriched marker
matches unscored. Every planned program remains in the output, including
low-power programs. Duplicate gene symbols are rejected before training rather
than silently assigning the same symbol to its last occurrence. A future data
aggregation policy must be explicit; the original large Norman file has not
been reprocessed in this environment.

The full-label shuffle is exploratory because it holds fixed a feature selected
using the original controls rather than repeating the selection for permuted
controls. A feature-background percentile is descriptive: different features
have different activation prevalence and power, and a selected named feature is
not exchangeable with random features. Neither readout warrants a specificity
claim.

The separate `program_modulation_specificity.py` tests whether eligible
single-perturbation identities have different rank locations on the selected
feature. It holds controls, feature selection, eligible groups and the fitted
label-blind representation fixed. The statistic is the weighted between-group
sum of squares of pooled midranks. Identities are permuted among eligible
perturbed cells, within supplied exchangeability strata. Each p-value includes
the observed assignment via `(1 + exceedances)/(1 + permutations)`.

All nine planned programs, including unmatched programs at p=1, enter Holm
correction. The default 999 permutations allow sufficiently small p-values for
this family; 50 permutations have a minimum p of 1/51, too coarse for one Holm
rejection among nine hypotheses. Duplicate selected features reuse the same
randomization evidence and remain separate planned hypotheses.

A rejection establishes conditional cell-level identity-associated rank
heterogeneity. Unequal generic response severity can also reject. It does not
establish a pathway-specific mechanism, biological replication, SAE superiority,
or novelty. Batch/guide structure must justify the chosen exchangeability unit.
Completely confounded strata are reported as uninformative.

## Deployed-default null calibration

A separate fixed simulation ran 200 independent null worlds, each with 100
controls, ten perturbations with 40 cells each, nine distinct stipulated matched
features and 999 permutations. The familywise outcome is any Holm rejection
among nine programs. It produced **6/200 rejections, 3.0%**, with exact 95%
binomial interval **1.11–6.42%**. This interval includes the nominal 5% rate.
There was no detected excess rejection under this simulated null; the experiment
does not prove validity for arbitrary experimental designs. Five paired uniform
generic-offset checks preserved rank-test results exactly; they are not five
additional independent null trials.

The complete 200 per-seed outcomes, parameter choices, source hashes, package
versions and summary are under
[`data_program_modulation_audit/calibration/`](../data_program_modulation_audit/calibration/).
These stipulated-feature simulations assess the conditional rank test, not the
whole biological feature-discovery process.

## Workflow and remaining experiment

The original nine suites passed during cloud setup. New tests cover the
counterexamples, marker-match failures, full planned-family retention, identical
features, batch confounding, permutation reproducibility, input validation, and
protection against overwriting evidence. The updated CLI was executed on
synthetic raw-count AnnData with actual CPU SAE training, optional strata, and
hashed provenance. Its unmatched programs were correctly left unscored. That
smoke test provides no biological result.

The [new validation plan](../config/program_modulation_validation_plan.json)
preserves v1.1 as historical evidence and states the scope for the next real-data
run. Seed agreement is computational robustness; it does not create biological
replication. If claiming findings across three seeds, correct over all 27 planned
program-by-seed hypotheses rather than selecting promising seeds afterward.

Required next steps are the full Norman run, transparent direct marker-score
baselines, target-transcript omission sensitivity, design-supported guide/batch
checks, and independent validation of any proposed mechanism. Direct PCA/NMF
feature baselines remain unrun: the earlier causal ladder instead trained an SAE
on each reduced representation.

The Norman source and GitHub API initially returned network-policy 403 responses.
Access later recovered: the raw file was downloaded and fingerprinted, preparation
completed, and PR #3 was updated through the REST API. A fresh full-data marker
baseline and explicit sensitivities are now archived; see
[NORMAN_CONTINUATION_RESULTS.md](NORMAN_CONTINUATION_RESULTS.md). The user selected
CPU execution and full-size SAE training is running with checkpoints. No completed
SAE result, biological novelty review or publication submission is claimed. See
the [manuscript draft](MANUSCRIPT_DRAFT.md) for the current bounded paper scope.
