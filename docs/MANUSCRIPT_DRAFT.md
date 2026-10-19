# Measured hematopoietic gene-program modularity depends on decomposition: a computational audit of sparse-autoencoder interpretation

**Internal manuscript draft, 2026-10-18. Not ready for submission.** Primary numerical
results below come from preserved Gate-0 artifacts and a fresh, archived Norman
external-marker analysis. Historical causal analyses lack their original result
JSONs and are excluded from the primary results. The program-modulation design audit uses synthetic data and tests statistical
interpretation, not biological discovery. No novel biological mechanism is claimed.

**Authors and affiliations:** TODO: author-confirmed author list, order, affiliations,
and corresponding-author contact. **Mentorship:** TODO: confirm scientific supervision
and each mentor's role. No mentor approval is implied by this draft.

## Abstract

Sparse autoencoders (SAEs) can describe single-cell expression as overlapping gene
programs, but biologically recognizable features need not imply a method-independent
program structure. We evaluated an initial hypothesis that erythroid commitment is
more distributed than granulocyte commitment. Following exploratory metric and marker
development, a version-locked analysis assessed mouse Paul15 data (2,730 cells,
2,012 genes) and human Setty CD34+ marrow data (4,142 cells, 2,032 genes), with ten SAE
seeds per dataset. A control-referenced decoder-loading measure combined a fixed
effect floor with an abundance-matched permutation null. Comparisons included PCA,
NMF, co-expression regulons, and ICA. No SAE seed supported the original hypothesis
in either dataset. Median erythroid/granulocyte above-control program counts were
1.5/3 in mouse and 1/2 in human. Only mouse NMF supported the decision rule among
ten method-by-dataset combinations. Marker omissions changed none of 34 mouse or
36 human decisions; additional sensitivity analyses also failed to restore support.
These findings establish dependence of an operational modularity measure on the
analysis representation and preprocessing choices. They do not identify an invariant
biological partition or a causal lineage mechanism. In a separate synthetic audit, legacy quantitative modulation gates passed in
20/20 simulations of an identical generic response across all perturbations. A
conditional identity test rejected in 0/20 generic-response simulations and 20/20
specific-response simulations. This counterexample shows why program responses
and feature-background rankings alone do not establish identity-specific effects.
A fresh heldout-control marker baseline found 37 heme-panel associations among
105 Norman perturbations; 29 survived target omission, RNA-depth matching and
every negative-control reference. KLF1's TF-panel AUC fell from 0.615 to 0.425
after its own transcript was removed, while heme and membrane associations
persisted. These are conditional marker-score associations, not an SAE advantage
or a new mechanism. Precise nulls and comparison baselines are necessary before
descriptive features become causal biological claims.

## Introduction

Erythroid and granulocyte differentiation involve recognizable transcription factors,
biosynthetic machinery, and terminal effector genes, providing a test of whether
features represent reproducible organization rather than appealing annotations. Single-cell studies have
characterized heterogeneous myeloid progenitors and human differentiation
trajectories [1,2]. A decomposition of those data, however, measures expression
through its own constraints: variance maximization, nonnegative factors, independent
components, or sparse reconstruction.

An SAE expresses observations through a dictionary of learned features with sparse
activity. Matching a feature to familiar genes demonstrates a descriptive
relationship. It does not determine how many biological programs exist, whether
another method would recover the same partition, or which perturbations regulate
the feature. These questions require different controls. The present project began with a directional claim: granulocyte
commitment was treated as unified while erythroid commitment was considered
distributed. Exploratory metric and marker corrections changed its interpretation,
motivating an explicit final decision rule and comparisons across methods.

Related foundation-model work distinguishes model-internal SAE ablations from
cellular causal structure and quantifies annotation-selection and sign-bias
effects in validation [9]. Our expression-space analysis asks different questions;
it does not claim to introduce SAE interpretation or biological causal discovery.

We test that direction across two datasets and five decomposition families, then
audit a perturbation-validation design. The contribution is an inspectable account
of what these measurements support, without claiming that SAE use in single-cell
biology is original or that all biological modules must be representation-invariant.

## Methods

### Evidence set and analysis chronology

The primary evidence set comprises committed preprocessing and training reports,
decision JSONs, per-seed CSVs, baseline summaries, and robustness summaries under
[`data_g0/`](../data_g0/) and [`data_g0_human/`](../data_g0_human/). Analyses were
developed on these biological systems before the final Gate-0 specification was
version-locked. The lock preceded the final runs; it did not precede all exploration.
SHA-256 checks establish preservation of the specified files, not independent
prospective registration or protection against every researcher choice.

Later causal metric amendments followed inspected outcomes. Historical Replogle
Perturb-seq [5] and
CRISPRa records are useful context, but original machine-readable outputs were lost
with ephemeral compute instances. Their numerical claims are consequently excluded
from the primary evidence presented here. This draft reports preserved-artifact analysis and a design audit, not repeated
full-scale biological training.

### Data and preprocessing

Mouse data derive from the Paul et al. myeloid progenitor dataset [1]. The supported
preprocessing applies a source-specific log policy to the depth-adjusted expression
matrix, selects highly variable genes, retains available predefined lineage markers,
and scales genes. The resulting matrix contains 2,730 cells and 2,012 genes. The
human arm uses Setty et al. CD34+ marrow data distributed with Palantir [2]. Counts
are normalized to a library-size target of 10,000, log-transformed, subjected to
marker-aware gene selection, and scaled, producing 4,142 cells and 2,032 genes.
Scanpy tooling [4] supports mouse data access. Assay and cell-selection differences require source-specific preprocessing.

Retained genes and marker coverage are recorded in each
`preprocess_report.json` and `gene_names.csv`. Marker panels contain four erythroid submodules and three
granulocyte submodules, with Progenitor and Cycling panels as comparators. Human
symbols are specified in [`config/human_markers.json`](../config/human_markers.json).
These panels define the operational units tested; they are not a complete census
of all lineage biology.

### Sparse autoencoders and structural comparisons

Gate-0 SAEs use 512 nonnegative ReLU features, a linear decoder, reconstruction
loss with an L1 activity penalty, and ten random seeds per dataset. The recorded
fits used 300 epochs. Mean active-feature counts were approximately 35.8 in mouse
and 32.5 in human, within the target range of 20–50. Feature indices were seed-specific.

PCA, NMF, an absolute-Pearson co-expression regulon representation, and FastICA
provided structural comparisons. PCA and NMF used rank 512; NMF required the
nonnegative part of the scaled matrix. FastICA was reported at its convergent
rank of 50 with three seeds. The co-expression representation used correlations
with panel transcription factors and was deterministic; it was not motif-pruned
SCENIC. PCA copies of one fitted basis are not independent random-seed replicates.
Accordingly, these comparisons interrogate a practical analysis pipeline, including
representation-specific transformations, rather than isolating a single algorithm
under identical rank, input values, and effective capacity.

### Operational program measure and decision

For a marker panel of size K among G retained genes, each feature's enrichment is
the fraction of its absolute decoder loading assigned to those markers, multiplied
by G/K. Panel strength is the mean of the three strongest feature enrichments.
For each seed, nonnegative excess strength is measured relative to the Progenitor
panel. A panel contributes an above-control program count only when excess
exceeds 0.10 and its abundance-matched permutation q-value is below 0.05.

The null draws marker-size-matched genes from 20 bins defined by total absolute
decoder mass, using 1,000 draws and Benjamini–Hochberg adjustment across tested
panels [3]. This matches a loading-mass proxy, not measured transcript abundance
directly. The implementation samples within bins with replacement. Permutation
significance uses strength averaged across seeds; the resulting panel significance
is then shared across seed-specific excess decisions. Seeds therefore describe
algorithmic variability, not independent biological replication.

The original directional claim is supported for a seed only if erythroid has at
least two above-control panels, at least as many as granulocyte, and greater
normalized entropy of excess strength. Dataset-level support requires a 60% seed
majority. The rule is more specific than a comparison of median counts alone.
Implementation and frozen constants are in
[`src/v0b_v3_1_decision.py`](../src/v0b_v3_1_decision.py) and
[`config/preregistration_spec.json`](../config/preregistration_spec.json).

### Robustness and perturbation-validation design

Structural sensitivity analyses varied effect thresholds, significance cutoffs,
sparsity, seed count, and gene-panel size. Marker leave-one-out analyses examined
dependence on an individual marker. These evaluate stability within the tested
design and do not replace additional biological datasets.

The subsequent program-modulation design selects a feature using predefined marker
enrichment among control cells, tests its activation under perturbations, and
compares modulator counts with full-label shuffles and other eligible features.
These diagnostics answer different questions. Full-label shuffling removes the
overall perturbation-versus-control response, while a feature-background ranking
compares responsiveness with other learned features. Neither by itself tests
whether perturbation identities have distinct effects. Selection on controls also
uses control membership; describing it as entirely label-free is inaccurate.

A separate conditional identity test therefore retains controls and permutes
perturbation identities among eligible perturbed cells, optionally within nuisance
strata. A rank-location heterogeneity statistic assesses whether feature activity
differs among identities, with Holm correction across all planned programs.
Feature selection must remain invariant under these conditional permutations.
This test addresses cell-level identity heterogeneity, not pathway specificity,
mechanistic mediation, or replication across biological samples. Unequal generic
response severity can also produce heterogeneity.

The synthetic audit compared four worlds over seeds 0–19, with 500 controls,
20 perturbation identities, 60 cells per identity, and 40 planted features. Five
markers identified the target feature. Four identities shifted that feature in
the specific-response world; all twenty shifted it identically in the generic-response
world. A null world had no shifts. In an absent-marker world, constant markers
matched no feature while the last feature had a generic response. Activations were
supplied directly; no SAE was trained. Each analysis used 99 permutations. The
legacy quantitative gates were shuffle p < 0.05 and feature-background percentile
at least 90; the specification's manual biological-coherence criterion was not
scored. The conditional reanalysis included a marker-match gate and one planned
program, so it does not benchmark multiplicity across many programs. Exact binomial
intervals summarize frequencies over simulation seeds only. Corresponding worlds
share baseline draws within seed; each scenario has twenty independent seeds, not
eighty independent replicates across the four scenarios.

A separate fixed calibration used 200 independent null simulations, nine distinct
planned features with stipulated marker matches, and 999 permutations per simulation.
The generator supplied independent normal features. The outcome was any Holm-adjusted
rejection at familywise alpha 0.05. Five paired uniform-offset checks assessed rank
invariance; they were not additional independent simulations.

### Fresh external-marker baseline and sensitivities

A separate real-data analysis used the raw Norman K562 CRISPRa counts [6]. Its plan
was declared before examining these new baseline outcomes, informed by prior
project history; this remains an exploratory continuation rather than independent
prospective confirmation. The prepared dataset contains 111,445 cells, including
11,855 controls, 57,831 single-perturbation cells and 41,759 dual-perturbation cells.
All 105 adequately sampled single-perturbation labels were analyzed, without
requiring that their own target transcript be observed.

Duplicate symbols were summed in raw counts, reducing 33,694 columns to 33,660
unique symbols while preserving exact per-cell library totals. All-gene
normalization to 10,000 and log1p preceded marker subsetting. Within each of eight
technical gemgroups, an outcome-blind split assigned half the controls to score
fitting, yielding 5,926 discovery and 5,929 inference controls. Eligible markers
had nonzero expression in at least ten discovery controls and positive variance.
Scores averaged discovery-control-standardized marker expression with a 0.1
standard-deviation floor, requiring at least three eligible markers.

One-sided, tie-corrected Mann–Whitney U statistics and their null variances were
summed within gemgroups; AUC was weighted by within-group comparison pairs. All
nine programs crossed with 105 perturbations entered a single BH family of 945,
including unscorable p=1 records. Activation required q < 0.05 and AUC ≥ 0.55.
The target transcript was omitted where present, with at least three remaining
markers and a separate same-sized correction family. A controls-fixed identity
omnibus used 999 within-gemgroup permutations and Holm over all nine panels.

Additional sensitivities were declared after primary outcomes and before their
own outcomes. They retained the primary scores and control split. Comparisons
used each of four negative-control construct identities, with joint BH across
all 3,780 program–perturbation–reference hypotheses. Depth matching used
within-gemgroup log-UMI quartiles learned from discovery controls, with BH over
945 comparisons. Original and target-omitted scores were corrected separately.
Depth may be a perturbation mediator, so this sensitivity can overadjust genuine
biology. Per-gemgroup effects are descriptive, not biological replication.

### Direct PCA features

After observing the marker baseline and an explicitly interim epoch-20 SAE
inference check, we fixed a separate exploratory PCA analysis before its outcomes.
The identical normalized 111,445 × 4,081 matrix supplied 64 ordinary principal
components (randomized SVD, three power iterations, seed 0). No SAE was trained.
Both signs of each component entered the annotation search. Discovery controls
alone supplied gene associations and marker matching, with Bonferroni over all
eligible orientations; associations used the same inference controls and the
stratified test, with BH over all 945 hypotheses. Identity inference used 999
within-gemgroup permutations and Holm across all nine panels. All-cell fitting
still saw inference-control expression, so this split is not an independent
representation-training holdout. Neither capacity nor all analysis choices are
matched to the 2,048-feature primary SAE; no superiority test is asserted.

## Results

### The original direction was unsupported in both SAE analyses

The frozen SAE decision supported the original hypothesis in 0/10 mouse seeds
and 0/10 human seeds. Median above-control panel counts were 1.5 erythroid versus
3 granulocyte in mouse, and 1 versus 2 in human. The mouse decision JSON's prose
rounds the erythroid median to 2; the numeric field and per-seed data show 1.5,
which is the value reported here. These results reject support under the defined
rule; they do not prove that erythroid differentiation lacks distributed programs.

### Structural interpretation changed across decomposition pipelines

| Representation | Mouse erythroid/granulocyte | Human erythroid/granulocyte |
|---|---:|---:|
| PCA | 1/1 | 1/3 |
| ICA, rank 50 | 1/3 | 2/3 |
| NMF | 2.5/2 | 0/0 |
| SAE | 1.5/3 | 1/2 |
| Co-expression regulons | 2/2 | 1/2 |

Only mouse NMF satisfied the original dataset-level decision, in eight of ten
seeds. Thus one of ten method-by-dataset combinations supported the claim.
The SAE-associated opposite direction was echoed by ICA but was not universal:
mouse PCA and co-expression counts tied, while human NMF detected neither lineage
above the chosen baseline. The ten comparisons are descriptive, not independent population samples. Their disagreement
demonstrates that this operational answer depends on the representation and its
associated preprocessing and capacity choices.

### Tested sensitivity analyses did not restore the original claim

None of twelve threshold combinations per dataset supported the original direction.
Marker omission changed no verdict among 34 mouse and 36 human omissions. Extending
SAE fits to twenty seeds retained the unsupported decision in both datasets; mouse
granulocyte count remained three in all twenty seeds, while human granulocyte count
remained two. Sparsity sweeps also preserved the decision. Gene-panel sweeps covered
approximately 1,000–3,000 selected genes and did not restore support. At the smallest
human panel, counts tied at one rather than reversing the unsupported hypothesis.

These sweeps have qualifications. Fixed regularization pushed active-feature counts
outside the target quality band at gene-panel extremes. Repeated seeds assess a
fixed sample of cells. Stability of a decision can coexist with an incomplete panel,
a biased metric, or weak biological identifiability.

A secondary TRRUSTv2 [7] comparison in the human arm found curated erythroid target
enrichment in all ten SAE seeds, with median unadjusted p = 6.44 × 10^-5. This is
a structural association based on a curated database and feature selection rule,
not a perturbation experiment or proof that the identified feature implements a
particular transcription-factor mechanism.

### A generic response passed both legacy quantitative gates

| Synthetic world | Legacy quantitative gates passed | Current conditional identity rejection |
|---|---:|---:|
| Four specific responders | 20/20 | 20/20 |
| Identical response in all twenty identities | 20/20 | 0/20 |
| No response | 0/20 | 0/20 |
| No marker-feature match; last feature responds | 20/20 | 0/20, all unmatched |

The generic-response world had no identity-specific effects by construction, yet
both legacy quantitative gates passed in every seed. Consequently those two gates
are insufficient evidence for identity specificity; this does not evaluate the
unscored manual-coherence gate. The specific-response world detected all four true
responders in every seed, with one additional false-positive modulator in each of
two seeds. In the absent-marker world, historical selection returned index −1 and
silently analyzed the last feature. The current marker gate instead marked all
twenty runs unmatched, preventing that invalid interpretation.

Exact 95% binomial intervals are 0.832–1.000 for 20/20 frequencies and 0–0.168 for
0/20 frequencies. Thus zero observed generic-response rejections does not establish
universal calibration or a precise low error rate. These are simulation-seed
intervals, not biological uncertainty. The retained
[legacy results](../data_program_modulation_audit/legacy_counterexamples_20seeds.json),
[current conditional reanalysis](../data_program_modulation_audit/specificity_reanalysis_20seeds.json),
and [intervals](../data_program_modulation_audit/rates_and_intervals.json)
contain the numerical evidence. The inference failure concerns the interpretation
of a diagnostic, not the validity of an experimentally observed biological program.

The separate nine-program null calibration produced 6/200 familywise rejections
(3.0%; exact 95% interval 1.11–6.42%), without evidence of excess false positives
under that generator. All five paired offset-invariance checks passed. The
[calibration summary](../data_program_modulation_audit/calibration/summary.json)
and its provenance record identify the exact executed sources. This finite
calibration does not validate arbitrary real-data exchangeability, marker matching,
or trained-SAE inference.

### Raw marker responses distinguish self-target signal from broader association

The fresh baseline found seven TF-panel, 37 heme-panel, 25 membrane-panel, six
granulocyte-TF and three progenitor associations among 105 single perturbations.
Cycling had no positive association under the fixed AUC criterion. Erythroid
effector, granulocyte-primary and granulocyte-secondary panels were unscorable
because only two, one and two markers, respectively, met discovery eligibility.
They are not negative biological findings.

Target omission retained six TF-panel, 37 heme, 25 membrane, four granulocyte-TF
and three progenitor associations. Requiring the original association also to
survive depth matching and every control reference retained three, 29, 20, four
and two, respectively. These intersections are descriptive robustness, not
independent validation. All comparisons and references are archived in
[`data_norman_continuation/`](../data_norman_continuation/).

KLF1 illustrates the direct-target trap. Its TF-panel AUC was 0.615, but fell to
0.425 after excluding KLF1 from its own score and no longer met the activation
criterion. Its heme and membrane AUCs were 0.626 and 0.580; those signatures do
not contain KLF1, and both survived the sensitivities. Removing CEBPA or CEBPE
left only two granulocyte-TF markers, so those omission comparisons were
unscorable rather than clean negatives.

CBL and PTPN12 retained heme and membrane associations across the sensitivities.
UBASH3B's heme AUC was 0.571 against pooled inference controls, but its lowest
control-reference AUC was 0.548, below the 0.55 floor. The earlier candidate trio
therefore cannot be endorsed without qualification. All six adequately covered
panels, including Cycling, showed conditional identity heterogeneity at the
permutation floor (raw p=0.001, Holm p=0.009). This broad pattern does not establish
erythroid pathway specificity, and no SAE advantage follows from a marker baseline.

### Ordinary PCA recovered different annotated responses

The 64 components explained 35.77% of expression variance. Ery_TF matched the
positive orientation of component 12 (zero-based), with 28 positive associations;
Cycling matched component 3, with 52. Both rejected the identity null (Holm
p=0.009). The other seven panels were unmatched and unscored, not biological
negatives. Ery_Effector's adjusted annotation p=0.0561 failed the fixed threshold.

The TF-matched feature's leading associated genes include globins, mitochondrial
transcripts and LMO2. Its annotation does not make it equivalent to the direct
TF score. Similarly, 52 associations for the Cycling PCA feature coexist with
zero for the direct Cycling marker score because these are different measurements.
These results show that ordinary PCA can recover annotated identity-associated
features; they do not rank its biological fidelity or demonstrate an SAE advantage.
Complete comparisons and failures are archived alongside the marker baseline.

The full-size 2,048-feature, 200-epoch SAE run is ongoing on CPU; no completed
new SAE result is included. Its model, optimizer, means, input fingerprints and
RNG states are checkpointed. Small tests verified exact uninterrupted/resumed
training equivalence. See [the continuation record](NORMAN_CONTINUATION_RESULTS.md)
for complete method details and actual run status.

## Discussion

The preserved evidence supports a narrower conclusion than either the original
asymmetry story or its reversal: program counts produced by this metric depend on
the chosen decomposition pipeline. Robustness within an SAE does not establish
invariance across models. Conversely, disagreement does not identify which method
is closest to a true biological organization. Methods can emphasize complementary
levels of a hierarchy, and rank or preprocessing differences can change resolution.
These results support operational dependence, not universal artifact claims or
absence of biological modules.

A diagnostic can measure a real response while answering the wrong question.
Concentration of a generic response in an erythroid-like feature would not establish
a distinct set of erythroid regulators. Conditional
identity heterogeneity is a more targeted question, but still admits differences
in stress severity, guide quality, or cell-state composition. Establishing pathway
specificity requires matched generic-response controls, appropriate gene-expression
baselines, and independent evidence.

Historical causal work in this project cannot settle those questions. All-cell
SAE training included perturbed cells even when earlier text specified control-only
training. Matching/test cross-fitting did not create an independent model-training
holdout. The causal representation ladder also trained an SAE on PCA/NMF inputs;
it did not directly compare SAE features with plain PCA or NMF features. Together
with missing original outputs, these issues prevent the historical experiments
from supporting a publication claim of causal superiority. A matched embedding arm,
held-out validation, and durable evidence capture remain required experiments.

## Limitations and next experiments

There is one biological dataset per species, and human CD34+ selection restricts
the states observed. Predetermined, unequal lineage panel collections constrain the
counting measure. Matching on decoder mass is not a universal biological null, and
shared significance across seeds limits independence. NMF rectification, unequal
ICA rank, and deterministic comparator fits complicate attribution to algorithms
alone. Biological effect sizes and uncertainty across independent samples are not
available from these summaries.

Priority work is to broaden the synthetic audit to graded generic severity,
multiple programs and experimental dependence; finish the full-size SAE runs,
compare them with the now-archived external-marker baseline and direct PCA/NMF
features, and evaluate effects with defensible guide or sample structure. Norman
construct strings often reflect vector orientation rather than independent
sgRNA reagents; they must not be treated as biological replicates.
Independent dataset replication and primary-literature review are needed before
any potentially interesting modulator can be characterized as new biology. None
of these future outcomes is assumed in this draft.

Known biology also limits novelty: the original Norman abstract explicitly reports
CBL/CNN1 synergy driving erythroid differentiation [6]. Prior work links FLI1 and
UBASH3B in erythroleukemia, using interventions and endpoints different from our
CRISPRa marker-score analysis [8]. Neither source independently validates our
scores. The [targeted literature review](LITERATURE_REVIEW.md) documents what was
read and what remains unchecked, including the original Norman full text and
supplements. We do not characterize PTPN12 or UBASH3B as newly discovered regulators.

## Data and code availability

Code and small structural evidence artifacts are available in
[the project repository](https://github.com/jbar1125/SAE-scRNA-seq-Research).
Primary-source data are described in [1,2,6]. Large matrices and checkpoints are
excluded from Git; their preprocessing reports record provenance. The frozen mouse
decision bundle has SHA-256 `172861417a923a83705812a96e5ad1e555055e968323289879447e53beb84396`;
the human bundle has `cc8159c882cf326bca5bbbc7f74d110d4a81338ec0ba07e208eb935bc8f729d3`.
Each bundle concatenates its decision JSON and per-seed CSV. Verification is
documented in [REPRODUCIBILITY.md](REPRODUCIBILITY.md). The synthetic audit retains
its exact executed generator source and source hashes in
[`data_program_modulation_audit/`](../data_program_modulation_audit/). The historical
benchmark implementation is pinned to commit `a08db6d`; the subsequent conditional
reanalysis is identified by its source hash. The initial generator did not record
a working-tree dirty flag, a provenance limitation disclosed in its manifest.

Relevant artifacts are each dataset's `v0b_outputs/v0b_v3_1_decision.json`,
`v0b_outputs/v0b_v3_1_per_seed.csv`, `baselines/baselines_results.json`,
`grn/grn_baseline.json`, `ica/ica_baseline.json`, and the `robustness/` and
`sensitivity/` summaries. The TRRUST result is
[`data_g0_human/trrust/trrust_crossref.json`](../data_g0_human/trrust/trrust_crossref.json).
Historical causal evidence limitations are detailed in
[COMPONENT2_RESULTS.md](COMPONENT2_RESULTS.md).

Fresh Norman baseline and sensitivity JSONs, a 945-row table, data fingerprints,
source hashes and technical-overlap metadata are retained under
[`data_norman_continuation/`](../data_norman_continuation/). Large source files and
resumable SAE weights remain outside Git. Their source fingerprints and preparation
policies are recorded; no unavailable full-SAE result is asserted.

**TODO:** deposit durable run archives with machine-readable outputs, exact code
revisions, environment manifests, and source-data access instructions; add archive
identifiers and a release tag. No archive identifier is invented here.

## Contributions, funding, and disclosures

**TODO, author confirmation required:** record contributions for conception,
software, analysis, validation, writing, and supervision; disclose computational
and writing assistance, including AI assistance where required. Confirm funding,
compute support, conflicts of interest, and acknowledgments. Confirm applicable
public-data ethics and competition or journal disclosure requirements; no ethics
determination, institutional approval, or exemption is asserted in this draft.

## References

Bibliographic fields were checked against Crossref, with retrieval provenance in
[`data_literature_review/`](../data_literature_review/). **TODO before submission:**
complete claim verification against primary texts and supplements. Metadata checks
do not imply that every cited paper was read in full.

1. Paul F, et al. Transcriptional heterogeneity and lineage commitment in myeloid
   progenitors. *Cell*. 2015;163:1663–1677.
   https://doi.org/10.1016/j.cell.2015.11.013
2. Setty M, et al. Characterization of cell fate probabilities in single-cell data
   with Palantir. *Nature Biotechnology*. 2019;37:451–460.
   https://doi.org/10.1038/s41587-019-0068-4
3. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and
   powerful approach to multiple testing. *Journal of the Royal Statistical
   Society: Series B*. 1995;57:289–300.
   https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
4. Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression
   data analysis. *Genome Biology*. 2018;19:15.
   https://doi.org/10.1186/s13059-017-1382-0
5. Replogle JM, et al. Mapping information-rich genotype-phenotype landscapes with
   genome-scale Perturb-seq. *Cell*. 2022;185:2559–2575.e28.
   https://doi.org/10.1016/j.cell.2022.05.013
6. Norman TM, et al. Exploring genetic interaction manifolds constructed from
   rich single-cell phenotypes. *Science*. 2019;365(6455):786–793.
   https://doi.org/10.1126/science.aax4438
7. Han H, et al. TRRUST v2: an expanded reference database of human and mouse
   transcriptional regulatory interactions. *Nucleic Acids Research*.
   2018;46(D1):D380–D386. https://doi.org/10.1093/nar/gkx1013
8. Wang J, et al. FLI1 induces erythroleukemia through opposing effects on
   UBASH3A and UBASH3B expression. *BMC Cancer*. 2024;24:326.
   https://doi.org/10.1186/s12885-024-12075-2
9. Kendiukhov I. Causal circuit tracing reveals distinct computational
   architectures in single-cell foundation models: inhibitory dominance,
   biological coherence, and cross-model convergence. *Bioinformatics*.
   2026;42(7):btag379. https://doi.org/10.1093/bioinformatics/btag379
