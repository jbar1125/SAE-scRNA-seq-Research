# Fresh Norman continuation evidence

This archive contains a completed full-data external-marker baseline and explicit
exploratory sensitivities, three completed full-size SAE seeds and a separate
one-seed control-only sensitivity.
See [the result record](../docs/NORMAN_CONTINUATION_RESULTS.md) and
[run/resume instructions](../docs/NORMAN_CPU_RUNS.md).

- `marker_baseline/gemgroup_seed0.json`: all nine panels, all 105 single
  perturbations, heldout-control split, fitted parameters, all 945 comparisons,
  separate target-omission family, identity omnibus and hashed provenance.
- `marker_baseline/sensitivity_seed0.json`: original computations reproduced
  exactly; all four negative-control references (3,780-comparison correction),
  discovery-defined depth matching (945 comparisons), target omission, and
  descriptive technical-partition AUCs. These are post-primary exploratory checks.
- `marker_baseline/comparisons.csv`: a 945-row portable comparison table.
- `marker_baseline/summary.json`: all-panel counts and candidates, including
  low-power and negative outcomes. No selected favorable reference is omitted.
- `direct_pca/`: a completed ordinary 64-component signed-PCA analysis, all 945
  comparisons, annotation failures and full source/model provenance. No SAE is
  trained in this baseline; unequal model capacities prevent a superiority claim.
- `direct_nmf/`: a completed 200-iteration capped fit with an explicit convergence
  warning, 945 comparisons and fitted-score/model hashes. Its nonconvergence and
  unequal capacity preclude a definitive method ranking.
- `sae_primary/`: all three completed 200-epoch primary JSON/log pairs and joint
  Holm over all 27 panel-by-seed identity tests. Its own manifest was created by
  the validated archive worker. No candidate combines the legacy background gate
  with at-least-two-seed agreement.
- `sae_control_sensitivity/`: one completed control-only 200-epoch seed, all nine
  outcomes and a compact summary. It does not complete a three-seed sensitivity.
- `phase_diagnostic/`: all 105 external S/G2M score comparisons with the Cycling
  PCA feature, 420 directional hypotheses and all 945 original PCA hypotheses in
  joint correction families. This post-primary diagnostic measures score directions,
  not cell-cycle progression, arrest or independent biological replication.
- `shared_inference/`: six representations using common discovery/inference
  controls, all 5,670 association hypotheses and 54 identity hypotheses, 9,999
  permutations, complete outputs and source/actual-score validation. Original
  PCA/NMF inference was reproduced before interpreting the combined family.
- `pca_attribution/`: all 105 comparison estimates and 428,505 gene contributions
  to the saved component's projected mean. Explicit original-score fidelity
  gates pass. This explains a score mean, not ranks or cellular causality.
- `nmf_convergence/`: actual bounded warm restart met tolerance after 380 further
  iterations; original capped model preserved. All 6,615 comparisons and 63
  identity hypotheses, model/score hashes and warnings are retained.
- `target_rna/`: all 109 direct/alternative transcript comparisons and a separate
  105-row post hoc pairing with the phase diagnostic. No old label is excluded.
- `k562_context/`: the exact historical eight-gene program across seven learned
  models and direct scores, all 8,400 association and 80 identity hypotheses;
  original direct scores reproduced first. The gene list was selected historically
  from these same data, so nominal values do not calibrate its full selection.
- `k562_context_sensitivity/`: exact primary reproduced, every control reference
  (4,200 comparisons), depth strata (1,050) and technical-partition estimates.
- `k562_gene_robustness/`: all 10,080 old/new hypotheses, each individual transcript
  and single-gene omission; all 105 labels retained. Thirty-eight primary calls
  also survive references, depth and all omissions, as a descriptive intersection.
- `k562_composition/`: all 11,025 hypotheses for the fixed-weight alternative
  denominator, after recovering 945 original gene/program comparisons. Neither
  denominator establishes absolute RNA abundance.
- `dual_marker/`: all 1,179 planned dual-marker contrasts and target-omission
  estimates, with exact compressed JSON, stratum counts, complete CSV and source
  fingerprints. These are scale-dependent additive departures, not validated
  molecular genetic interactions. See [the result record](../docs/NORMAN_DUAL_MARKER_RESULTS.md).
- `provenance/`: observed source fingerprints, exact count-conserving symbol
  aggregation, metadata-only technical overlap and timing preflights. The initial
  metadata heuristic omitted `gemgroup`; its explicit cross-tab verifies eight
  technical partitions with controls and all interventions. They are not
  biological replicates. Timing estimates are variable preflight measurements.
- `validation.json`: 15 standalone suites passed together and nine additional
  suites separately, for 24 distinct suites exercised. Actual CLI/training,
  frozen-integrity and review-rendering checks are distinguished from biological
  claims. The completed-model validator and its complete four-model record are
  in `provenance/`.
- `sae_interim_epoch20/`: a separate copied partial model's full-data inference
  diagnostic, with all nine outcomes. Twenty of the planned 200 epochs had
  completed; this is not a completed primary seed or an independent replicate.
- `sae_progress.json`: timestamped observation of a partial, mutable local
  checkpoint. It is not a completed result or a permanent weight archive.

Raw/prepared counts, model checkpoints and activations live outside Git under
`/workspace/datasets/sae-scrna/` and `/workspace/research-runs/` in this environment.
The raw file's measured checksum is an observed source fingerprint, not an
independently authenticated expected checksum. No model seed or technical
partition constitutes independent biological replication. Robust intersections
are descriptive sensitivity summaries. These associations do not establish a
new mechanism, novelty, pathway specificity or an SAE advantage.

`provenance/large_artifact_manifest.json` lists seven actual large review archives
and exact roundtrip-checked member hashes. They are local and outside Git.
GitHub's empty draft release has no model/data assets because the upload host is
blocked at the network proxy. Prepared archives are not public deposition.
