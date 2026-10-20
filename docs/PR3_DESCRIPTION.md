# Suggested PR title

Audit gene-program interpretation with completed CPU fits and archived evidence

# Description

This branch tests what biologically annotated features support in structural and
perturbation analyses. It preserves the failed original asymmetry hypothesis,
repairs an unmatched-feature index bug, demonstrates limits of legacy modulation
gates and conditional identity inference, and now archives the completed Norman
CPU experiments alongside direct baselines and publication materials.

Frozen specifications and both Gate-0 bundles remain unchanged. Historical causal
JSONs were lost; their old numbers are excluded from the primary manuscript.
Dated corrections cover all-cell training, the wrong-background essentiality
claim, the PCA/NMF-to-SAE ladder and unrun embedding comparisons.

## Completed evidence

- Frozen structural analyses:0/10 SAE seeds support the original hypothesis in
  each dataset; only mouse NMF supports it across five methods and two datasets.
- Synthetic audits: legacy gates accept identical generic responses 20/20;
  conditional identity inference still rejects all nine panels in 100/100
  variable-strength rank-one worlds. Invalid clustered-cell shuffling rejects
 100/100 versus 8/100 for a separate cluster-mean analysis. Complete outputs,
  source hashes, exact intervals and untestable cases are retained.
- Fresh heldout-control marker baseline: all 105 singles, all nine panels and
 945 tests. Heme 37 primary associations,29 surviving target omission, depth
  matching and all four control references. KLF 1 TF AUC.615 becomes.425 after
  own-transcript omission; downstream heme/membrane responses persist.
- All three 2048-feature primary SAE seeds finished 200 epochs. TF, effector and
  Cycling match in every seed (joint Holm 27 p=.027). KLF 1 is the sole repeated
  legacy candidate but background percentiles 63/73/68 fail the old gate. No
  candidate combines that gate and two-seed agreement. One control-only seed
  also finished and is separately archived; it is not biological replication.
- Direct rank 64 PCA and NMF completed. NMF reached its fixed 200-iteration cap
  without convergence; this limitation is preserved in every comparison.
- Common downstream inference across all six representations completed with
 5670 association hypotheses and 54 identity hypotheses,9999 permutations and
  identical discovery/inference controls. Every original PCA/NMF comparison
  reproduced before interpreting the combined family. All six match Cycling,
  but only PCA has upward associations(52). No superiority follows from counts.
- External phase diagnostic: no upward S/G 2 M calls;52 S and 37 G 2 M decreases.
  Among 50 jointly retained PCA calls,19 decrease in both phases. In the common
  family,20/52 meet both phase AUC<=.45 thresholds, a descriptive effect rule.
- Signed PCA attribution retains all 428505 gene contributions. Its projection
  closely matches recovered original scores(RMSE/SD=.000154,max AUC difference
  .0000378). In all 20 direction-discordant cases, represented phase genes lower
  the mean while other weighted genes raise it. This explains score means, not
  ranks, cell-cycle progression or molecular causality.
- All 131 double constructs and nine panels are retained in 1179-test marker
  contrasts, plus joint 2358-test log/linear sensitivity. Twenty-one heme contrasts
  persist across scales. Known CBL/CNN 1 biology is contextual recovery, not novelty.

## Publication and validation

The manuscript, four main figures, supplementary figure, claim-to-evidence ledger,
full result records and targeted literature review are versioned. Ten citations
have verified metadata. Two primary full texts and original Norman analysis code
were reviewed selectively; the original Norman paper/supplements remain unchecked.

Twenty-one distinct standalone suites were exercised(fifteen together and six
separately), including complete-family, score-recovery and mean-attribution tests.
Syntax, actual end-to-end research executions, exact completed-model artifacts,
frozen hashes and deterministic figure checks pass. Full command/status scopes
are in docs/CONTINUATION_VERIFICATION.md.

No new biological mechanism, independent replication or SAE advantage is
established. The paper reports representation dependence and score-response
discordance under the tested conditions. Authorship/affiliations, disclosures,
primary-text claim review and durable large-artifact deposition remain required
before submission. This remains a draft PR.
