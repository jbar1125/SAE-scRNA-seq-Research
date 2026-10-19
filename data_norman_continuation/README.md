# Fresh Norman continuation evidence

This archive contains a completed full-data external-marker baseline and explicit
exploratory sensitivities. It contains no completed full-size SAE result yet.
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
- `dual_marker/`: all 1,179 planned dual-marker contrasts and target-omission
  estimates, with exact compressed JSON, stratum counts, complete CSV and source
  fingerprints. These are scale-dependent additive departures, not validated
  molecular genetic interactions. See [the result record](../docs/NORMAN_DUAL_MARKER_RESULTS.md).
- `provenance/`: observed source fingerprints, exact count-conserving symbol
  aggregation, metadata-only technical overlap and timing preflights. The initial
  metadata heuristic omitted `gemgroup`; its explicit cross-tab verifies eight
  technical partitions with controls and all interventions. They are not
  biological replicates. Timing estimates are variable preflight measurements.
- `validation.json`: 15 standalone suites passed together, plus the completed-
  result archiver, direct-PCA and dual-marker suites passed separately; frozen integrity,
  actual CLI/training checks and rendering checks.
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
