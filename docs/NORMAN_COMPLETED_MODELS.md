# Interpretation of the completed Norman SAE fits

Updated 2026-10-20. All three declared all-cell seeds and the separate control-only
seed finished 200 epochs. The full-size experiment is no longer pending. The
archived results do not establish a new regulator or an SAE advantage.

## Original analysis, complete three-seed family

Every all-cell fit used 111,445 cells, 4,081 genes, 2,048 latent features, k=32 and
batch size 1,024. The original all-control annotation matched Ery_TF,
Ery_Effector and Cycling in every seed. Heme, membrane and the other four panels
failed the annotation gate in every seed; they remain unscored rather than negative
biological results. Each matched panel had raw identity p=.001 and joint Holm
p=.027 over all 27 panel-by-seed hypotheses. This tests conditional cell-level
rank heterogeneity, which can reflect unequal generic perturbation severity.

| Program | Legacy candidate counts, seeds 0/1/2 | Candidates in at least two seeds | Feature-background percentiles, seeds 0/1/2 |
|---|---|---|---|
| Ery_TF | 0/0/11 | None | 0/0/98 |
| Ery_Effector | 1/1/1 | KLF1 | 63/73/68 |
| Cycling | 0/0/0 | None | 0/0/0 |

The only repeated legacy candidate, KLF1, does not meet the old 90th-percentile
background gate in any seed. The TF candidate list appears in seed 2 alone.
Consequently there is no candidate that combines the old background criterion
with at-least-two-seed agreement. The background gate is itself descriptive and
was shown insufficient to establish specificity. Passing it would not establish
a discovery. Legacy per-panel candidate inference also reused annotation controls
and differs from the newer stratified, heldout-control baseline.

The epoch-20 heme annotation does not survive to the fixed 200-epoch endpoints.
This is a descriptive trajectory observation; the interim checkpoint was a
readiness diagnostic, and no epoch was selected to rescue the claim. Feature
indices are specific to their models and cannot identify a shared feature across
seeds. Counts of features ever active were 474, 461 and 466 of 2,048; these are
encoding diagnostics, not a calibrated model-quality or biological criterion.

See the [generated all-seed table](NORMAN_SAE_RESULTS.md),
[complete joint summary](../data_norman_continuation/sae_primary/joint_summary.json)
and individual full JSON/log files in the same archive.

## Separate control-only sensitivity

One seed trained on 11,855 controls and encoded all 111,445 cells. It matched heme,
effector and Cycling; each had raw identity p=.001 and within-arm Holm p=.009.
Heme had KLF1 as its only legacy candidate and background percentile 82; the other
two panels had no candidates. The remaining six panels were unmatched. The
complete original four-fit family gives Holm p=.036 for each matched identity
test across all 36 hypotheses, including failures. This retrospective combination
is exploratory; the original declared primary 27-test results are preserved.

This real control-only model does recover named features despite the earlier
synthetic control-only failure. It does not establish general success or failure
of control-only training: only one real seed was run. Its 2,014 ever-active
features and different annotations show sensitivity to the training population,
without identifying the best population. Its gene selection used all cells, and
annotation used every control, including controls used in model training.

Complete [control-only records](../data_norman_continuation/sae_control_sensitivity/summary.json)
and a [four-model validation record](../data_norman_continuation/provenance/completed_models_validation.json)
are archived. The validator checks exact model/activation hashes, checkpoint
epochs and training dimensions, original input SHA, unchanged scientific sources
and all planned hypotheses. Large checkpoints and activations remain outside Git.

## Common-inference follow-up

A [separate exploratory plan](../config/norman_shared_inference_plan.json) was
committed after these original outcomes and before new downstream outcomes. It
applies discovery-only annotation, the same inference controls and stratified
tests to every saved SAE model, alongside direct PCA and capped NMF. All six arms
enter joint BH over 5,670 perturbation associations and Holm over 54 identity
tests. A pre-execution amendment increases identity permutations to 9,999 because
999 permutations cannot resolve a first Holm rejection in a 54-test family.
Original results are preserved. The follow-up has now completed with full
score/data checks: see [common inference](NORMAN_SHARED_INFERENCE.md) and the
[signed-weight diagnostic](NORMAN_PCA_ATTRIBUTION.md).

The common protocol holds out annotation/reference controls conditional on fitted
models; it does not hold cells out of all-cell representation training. Unequal
rank, marker eligibility, shared cells, nonconverged NMF and the one-seed control
arm remain limitations. Neither the original nor the follow-up establishes
independent biological replication, pathway specificity, novelty or superiority.
