# Program-modulation design audit

These are simulations of the downstream statistic, not biological discoveries or
trained-SAE performance estimates. Activations are generated directly. The four
worlds share baseline draws within each simulation seed; each world's 20 seeds
provide its own independent simulation-frequency estimate.

| Simulated world | Historical quantitative gates pass | Current conditional identity test rejects |
|---|---:|---:|
| Four specific responders among 20 perturbations | 20/20 | 20/20 |
| All 20 perturbations identically activate one axis | 20/20 | 0/20 |
| No perturbation effect | 0/20 | 0/20 |
| No expressed marker match; unrelated last feature responds | 20/20 | 0/20, all unmatched |

The historical gates are full-label shuffle p < 0.05 and feature-background
percentile >= 90. The specification's manual biological-coherence criterion is
not scored. Passing these two gates does not establish perturbation-identity
specificity. In the absent-marker world, the historical selector returns -1 and
silently scores the last feature; the current implementation leaves the program
untested. The specific world recovers all four planted responders in every seed;
two seeds also contain one extra modulator call.

This benchmark uses **one planned program and 99 permutations**, not the deployed
nine-program, 999-permutation analysis. Minimum permutation p is 0.01. Exact 95%
Clopper-Pearson intervals are [0.8316, 1] for 20/20 and [0, 0.1684] for 0/20.
They quantify simulation frequencies, not biological replication or the precise
false-positive rate of the deployed procedure. Its separate null calibration is
under `calibration/`.

- `legacy_run/metadata.json`, `results.jsonl`, and `run.log` retain the original
  run byte for byte. The JSONL contains all 80 complete historical outputs.
- `legacy_counterexamples_20seeds.json` is the original summary, also unaltered.
  Its auxiliary identity calculation is an ungated rank diagnostic. In particular,
  it does not protect the invalid -1 feature. Use the current reanalysis below
  for marker-gated inference.
- `executed_generator_source.txt` is the exact generator source copied during
  the original run before subsequent output-safety and provenance enhancements.
  Historical `program_modulation.py` and `causal_grounding.py` are loaded in full
  from Git commit `a08db6d5ef4440adb9f1e1613e96856ab6515533`; their SHA-256 values
  appear in the original metadata.
- `specificity_reanalysis_20seeds.json` reconstructs those recorded worlds and
  exercises the actual current `program_modulation_specificity` implementation.
  It records the SHA-256 of its input summary, the current generator, and the
  diagnostic source. The summary hash identifies the unaltered legacy summary
  above; the generator hash differs from the original executed source because
  the reanalysis entry point and integrity safeguards were added subsequently.
- `rates_and_intervals.json` combines the historical quantitative gate counts
  and current diagnostic counts with exact intervals. `source_manifest.json`
  supplies hashes for the retained artifacts. The original generator did not
  record its working-tree dirty flag; this is explicitly unknown. The current
  CLI records that flag and its own source fingerprint.

Reproduce with a fresh output directory. Equivalent commands were executed:

```bash
python3 src/benchmark_program_modulation.py \
  --out-dir /tmp/program-modulation-reproduction --seeds 20 --n-shuffle 99
python3 src/benchmark_program_modulation.py \
  --reanalyse-summary data_program_modulation_audit/legacy_counterexamples_20seeds.json \
  --out-dir /tmp/program-modulation-specificity-reproduction --n-shuffle 99
python3 tests/test_program_modulation_counterexamples.py
```

The current generator preserves the original world construction and historical
metric calculations. Runtime, output-safety metadata, and provenance fields can
differ across reruns; the original artifacts remain preserved here.

The later [`boundary_audit/`](boundary_audit/) retains 100 seeds for six analyses
of generic severity, complete batch confounding and correlated cells. Graded
severity yielded all-nine identity rejections in 100/100 worlds despite rank-one
group effects. Cell shuffling under a clustered population null rejected in
100/100 worlds, while a separate independent-cluster-mean analysis rejected in
8/100 (exact 95% interval 3.52–15.16%). These are scope/dependence examples, not
proof of pathway specificity, universal calibration or biological discovery.
See the [full audit](../docs/PROGRAM_MODULATION_AUDIT.md) for all six outcomes.
