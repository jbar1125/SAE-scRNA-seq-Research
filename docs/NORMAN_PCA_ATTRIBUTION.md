# Why a control-annotated Cycling PCA score can increase as phase scores decrease

The [declared diagnostic](../config/norman_pca_attribution_plan.json) completed on
all 105 single perturbations. It supplies an algebraic explanation of the observed
score means. It does not identify a molecular mechanism or establish PCA failure.

Randomized-SVD fit_transform scores were first recovered with the original
settings and their complete 945 comparison estimates verified by the
[common-inference run](NORMAN_SHARED_INFERENCE.md). The saved component-3 linear
projection was then compared against those recovered scores rather than presumed
identical. Its RMSE was 0.000154 original-score standard deviations and its maximum
absolute AUC difference across all 105 perturbations was 0.0000378. Both pass the
declared .01/.005 fidelity bounds. The projected means therefore approximate this
original component closely enough for the declared descriptive interpretation.

For every perturbation and input gene, the contribution is its saved component
weight times the perturbation-minus-inference-control mean expression difference,
computed within gemgroup. Fixed perturbed-cell × control-cell pair counts weight
the partitions, as in the stratified AUC. All 105 × 4,081 = 428,505 contributions
were retained. They sum to the projected mean difference, with maximum absolute
error 2.22e-16. This is a mean decomposition; it does not decompose ranks or AUCs.

Of 97 requested external phase genes, 32 enter the PCA input. Their positive,
negative and zero-weight contributions are separated from all other genes.
Twenty perturbations have original PCA AUC≥.55 while both external phase AUCs
are ≤.45. All 20 have positive projected mean changes and a net positive
contribution from genes with negative component weights. A negative weight can
raise the score when its gene expression decreases.

A post hoc descriptive breakdown of the retained partition sums sharpens that
explanation: the net contribution of the 32 represented phase genes is negative
in every one of these 20 cases. Other genes with positive weights contribute
positively in all 20, as do other genes with negative weights. Median contributions
are approximately −0.261 for the net phase class, +0.455 for other positive-weight
genes and +0.585 for other negative-weight genes. Medians need not sum to the
median total (+0.775). The increase is therefore supported by other weighted
expression changes despite a falling phase-gene contribution.

Negative-weight increases are not unique to this selected set: among all 105
perturbations, 94 have positive sums from negative-weight genes and 85 have positive
projected mean changes. These are descriptive counts in the observed data, not
independent successes or a new significance family. The 20-case subset is defined
by fixed AUC thresholds; the separate directional Holm analysis has 19 doubly
decreased cases among its 50 retained PCA calls. Do not interchange those counts.

The result shows why a control-associated label does not fix a projection's
intervention response. It does not prove phase arrest, true proliferation,
regulatory causality, biological novelty or superiority of another decomposition.
Gene coverage, covariance and normalization differ between this projection and
the direct signature scores. Independent biological evidence remains absent.

The [full result](../data_norman_continuation/pca_attribution/results.json),
[105 comparison rows](../data_norman_continuation/pca_attribution/comparisons.csv),
complete compressed gene table, log and
[validation record](../data_norman_continuation/pca_attribution/validation.json)
are archived. The input, scientific sources, original model, recovered scores and
common-analysis output are fingerprinted. Every gene and partition sum was checked.
