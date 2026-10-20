# Exploratory dual-perturbation marker contrasts

The 131 observed dual labels contain 41,759 cells, with at least 54 cells per pair
and sampled single arms for every pair. This screen was declared before examining
dual expression outcomes, after the single-perturbation baselines and the abstract
of the original Norman paper had been read. It is an exploratory continuation,
not an independent prospective or biological replication experiment.

## Fixed score and contrast

All original discovery-control marker parameters and the exact control split
reproduced before analysis. The nine original panels retain their eligibility
rules. Scores average fixed standardized log-normalized markers and are expressed
in discovery-control score SD units. Within each usable gemgroup:

`interaction = mean(AB) - mean(A) - mean(B) + mean(inference controls)`

Stratum weights depend on dual-cell counts; each group needs at least two cells
per used stratum, at least 20 pooled cells, and at least 80% dual-cell coverage.
The independent-cell Welch–Satterthwaite approximation supplies two-sided p-values
and pointwise 95% intervals. Holm correction covers all 9 × 131=1,179 comparisons,
including unscorable p=1 hypotheses. A flag also requires absolute contrast ≥ 0.2 SD.
Omitting both target transcripts is a separate same-sized Holm family, requiring
at least three remaining markers. Same-sign intersection is the reported
target-omission consistency rule. See [the fixed plan](../config/norman_dual_marker_plan.json).

## All-panel outcomes

| Panel | Primary flags | Same-sign target-omission flags | Above additive | Below additive |
|---|---:|---:|---:|---:|
| Ery_TF | 1 | 0 | 0 | 0 |
| Ery_Heme | 25 | 25 | 17 | 8 |
| Ery_Membrane | 8 | 8 | 6 | 2 |
| Ery_Effector | Unscored | Unscored | Unscored | Unscored |
| Gran_TF | 14 | 6 | 1 | 5 |
| Gran_Primary | Unscored | Unscored | Unscored | Unscored |
| Gran_Secondary | Unscored | Unscored | Unscored | Unscored |
| Progenitor | 11 | 11 | 8 | 3 |
| Cycling | 10 | 10 | 9 | 1 |

The direction columns concern the same-sign target-omission set. For heme, the
targets are not members of the marker panel, so omission does not change the
score. An unchanged result is not extra independent evidence.

CBL+CNN1, reported as an erythroid synergy in the original paper's abstract, has
heme contrast 0.877 control-score SD (pointwise interval 0.633–1.120; Holm
p=4.61×10^-9). Its membrane contrast 0.191 does not meet the effect floor and has
Holm p=1. Its TF-panel contrast is also unflagged. The known heme association is
contextual recovery using the same source dataset, not a novel mechanism or
independent validation of our analysis.

Complete estimates, stratum counts, primary/omission results and failures are
retained as [compressed exact JSON](../data_norman_continuation/dual_marker/results.json.gz),
with a [portable 1,179-row table](../data_norman_continuation/dual_marker/comparisons.csv)
and [summary](../data_norman_continuation/dual_marker/summary.json). Decompress with
Python's `gzip` or `gzip -dc`; the archive manifest records the exact uncompressed
source hash as well as the compressed artifact hash.

## Log/linear sensitivity

The declared post-primary sensitivity completed after reproducing every original
primary and omission estimate exactly. Original marker eligibility and the split
remained fixed; linear normalized expression received its own discovery-control
standardization. Thus transformation and standardization both change. Joint Holm
covered all 2,358 log/linear hypotheses, separately for primary and target omission.
Consistency requires every estimate to pass that correction and the 0.2 SD floor
with the same sign.

| Panel | Consistent on both scales | Above additive | Below additive |
|---|---:|---:|---:|
| Ery_TF | 0 | 0 | 0 |
| Ery_Heme | 21 | 14 | 7 |
| Ery_Membrane | 6 | 4 | 2 |
| Ery_Effector | Unscored | Unscored | Unscored |
| Gran_TF | 6 | 1 | 5 |
| Gran_Primary | Unscored | Unscored | Unscored |
| Gran_Secondary | Unscored | Unscored | Unscored |
| Progenitor | 4 | 4 | 0 |
| Cycling | 6 | 5 | 1 |

Known CBL+CNN1 retains its heme contrast: 0.877 SD on the log score and 1.995 SD on
the linear score, with joint corrected p<1e-8 on each. Its membrane comparison
remains unflagged. Larger numerical SD effects on one scale do not mean stronger
biological synergy; each scale defines a different score and control variance.

Complete linear outputs and all joint p-values are in the
[compressed sensitivity JSON](../data_norman_continuation/dual_marker/scale_sensitivity.json.gz)
and [portable comparison table](../data_norman_continuation/dual_marker/scale_comparisons.csv).
These intersections are post-primary computational sensitivity summaries, not
independent confirmation. The three low-coverage panels remain unscored.

## Interpretation

These are departures from additivity on this particular marker-score scale.
Log transformation, library normalization, saturation, cell-state mixtures and
construct effects can produce deviations even without a molecular interaction.
Different marker panels also measure different outcomes. Independent-cell standard
errors can be too small when cells share unmodeled experimental units; gemgroups
are technical strata, not independent biological samples. Pointwise intervals
are not simultaneous familywise intervals. No dual-specific control-reference
or depth sensitivity has yet been completed.

The [log/linear sensitivity](../config/norman_dual_scale_sensitivity_plan.json)
does not turn these score contrasts into validated biological genetic interactions.
The original full text/supplements and independent experimental validation remain
necessary before a new biological claim.
