# CPU continuation verification record

Snapshot: 2026-10-18 ET. PASS means the stated scope was executed; it does not
establish biological novelty or validate a different experimental design.
Run commands from the checkout with `/workspace/.venvs/sae-scrna/bin/python`
and the cache/thread settings in [NORMAN_CPU_RUNS.md](NORMAN_CPU_RUNS.md).

| Command or check | Status | Evidence and scope |
|---|---|---|
| `python -m py_compile` on changed analysis/archiver/test files | PASS | Current sources compile. Syntax alone is not an end-to-end execution. |
| `python tests/test_v0b_logic.py` | PASS | Rerun after the continuation edits; original synthetic logic checks pass. Low-coverage fixture warnings are expected and retained in the combined log. |
| Every original suite plus checkpoint, unique-count, marker, sensitivity and identity tests | PASS | Fifteen standalone scripts executed together with zero skipped suites; [complete log](../data_norman_continuation/validation_logs/fifteen_suites.log). |
| `python tests/test_norman_result_archive.py` | PASS | [Separate log](../data_norman_continuation/validation_logs/result_archive.log): partial/incompatible runs rejected, complete family retained, interrupted publication recovered, existing evidence preserved. |
| `python tests/test_norman_direct_pca.py` | PASS | [Separate log](../data_norman_continuation/validation_logs/direct_pca.log): both signs, planted responses, no-match retention and discovery-only annotation conditional on fixed scores. |
| Actual checkpointed small SAE/AnnData CLI | PASS | Real small training, saved factors, model/optimizer/RNG continuation and output provenance; see [validation record](../data_norman_continuation/validation.json). Not a completed full-size SAE. |
| Full-data marker baseline and specified sensitivities | PASS | All 111,445 cells, 945 primary comparisons, 3,780 reference comparisons; original baseline reproduced before sensitivity. Complete outputs and hashes archived. |
| Full-data ordinary PCA command | PASS | 64 components, all 945 comparisons and annotation failures; source/model fingerprints and exact existing control split verified. No SAE trained. |
| `python src/audit_identity_boundaries.py --out-dir <fresh-directory>` | PASS | 100 seeds × six analyses, 600 unique complete records, all nine hypotheses; summary counts recomputed and source/output hashes verified. No biological model fitted. |
| `python src/make_identity_boundary_figure.py` twice | PASS | PNG and SVG byte-identical after layout correction; rendered figure inspected. |
| Four frozen spec hashes and both Gate-0 bundle hashes | PASS | Rechecked after new experiments; all match, with original evidence preserved. |
| Full-data direct NMF command | UNVERIFIED | Running under the fixed plan. Actual iterations, convergence warnings and outputs must be checked before a result is reported. |
| Three full-size 200-epoch SAE seeds | UNVERIFIED | Seed 0 is partial/running; seeds 1 and 2 are queued. Copied epoch 20 inference is an interim readiness check. The archive watcher validates actual completion rather than inferring it from a process. |
| Independent biological replication, novel mechanism or SAE superiority | UNVERIFIED | Computational seeds, gemgroups, annotations and cell-level p-values do not establish these claims. |

The two separately executed suites make 17 distinct standalone suites exercised;
this is not a claim that a single 17-suite invocation was run. The 100-seed boundary
experiment is a research execution, not an additional unit-test suite. Source
fingerprints in run provenance identify the executed code revision; later
documentation commits do not change the underlying scientific outputs.
