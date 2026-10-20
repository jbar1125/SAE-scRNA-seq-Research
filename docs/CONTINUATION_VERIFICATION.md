# CPU continuation verification record

Snapshot: 2026-10-20 UTC. PASS means the stated scope was executed; it does not
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
| `python tests/test_norman_dual_marker.py` | PASS | [Separate log](../data_norman_continuation/validation_logs/dual_marker.log): additive null, planted interaction, technical-composition artifact and unavailable-arm handling. |
| `python tests/test_norman_phase_diagnostic.py` | PASS | [Separate log](../data_norman_continuation/validation_logs/phase_diagnostic.log): original PCA plus both directions of phase hypotheses retained, separate omission family. |
| `python tests/test_norman_shared_inference.py` | PASS | Complete six-arm families and no-match hypotheses retained; selected-arm and incompatible score-recovery inputs rejected. Actual six-arm execution is recorded separately below. |
| `python tests/test_norman_pca_attribution.py` | PASS | Signed mean decomposition and a technical-composition counterexample exercised. |
| `python tests/test_norman_nmf_convergence.py` | PASS | Actual warm fit reduces error without changing initial artifacts; complete seven-arm correction retained. |
| `python tests/test_norman_target_rna.py` | PASS | Induced transcript detected; absent direct target and noninduced alternative retained. |
| `python tests/test_norman_k562_context.py` | PASS | All eight arms and unscored panels remain in expanded families; selected-arm analysis rejected. |
| Actual checkpointed small SAE/AnnData CLI | PASS | Real small training, saved factors, model/optimizer/RNG continuation and output provenance; see [validation record](../data_norman_continuation/validation.json). Not a completed full-size SAE. |
| Full-data marker baseline and specified sensitivities | PASS | All 111,445 cells, 945 primary comparisons, 3,780 reference comparisons; original baseline reproduced before sensitivity. Complete outputs and hashes archived. |
| Full-data ordinary PCA command | PASS | 64 components, all 945 comparisons and annotation failures; source/model fingerprints and exact existing control split verified. No SAE trained. |
| Full-data dual-marker screen and log/linear sensitivity | PASS | 1,179 pair-panel contrasts; original fits/split and all original estimates reproduced. Both scales enter joint 2,358-test Holm families. Complete outputs and failures are archived. No validated molecular mechanism follows. |
| `python src/audit_identity_boundaries.py --out-dir <fresh-directory>` | PASS | 100 seeds × six analyses, 600 unique complete records, all nine hypotheses; summary counts recomputed and source/output hashes verified. No biological model fitted. |
| `python src/make_identity_boundary_figure.py` twice | PASS | PNG and SVG byte-identical after layout correction; rendered figure inspected. |
| Four frozen spec hashes and both Gate-0 bundle hashes | PASS | Rechecked after new experiments; all match, with original evidence preserved. |
| Full-data direct NMF command | PASS, capped | Actual 200-iteration execution completed; convergence tolerance was not met. The warning, complete 945 comparisons and actual fitted-score/model/source hashes are retained. No definitive method ranking follows. |
| Full-data external S/G2M diagnostic | PASS | All 420 phase/direction tests plus all 945 original PCA tests retained in joint Holm families. Scores measure signature direction, not cell-cycle progression. |
| `python data_norman_continuation/provenance/validate_completed_models.py` | PASS | [Validation record](../data_norman_continuation/provenance/completed_models_validation.json): source/input/checkpoint/activation hashes, exact completed training dimensions, original 27-test summary and all four complete nine-panel families verified. |
| Control-only full-size SAE sensitivity | PASS, one seed | Actual 200 epochs on 11,855 controls; all 111,445 cells encoded. Complete output/log archived. This is not a three-seed sensitivity or an independent holdout. |
| Three full-size 200-epoch SAE seeds | PASS | All three completed; actual artifacts checked twice, compact JSON/logs and joint Holm 27 archived and pushed. No candidate combines the old background gate with two-seed agreement. |
| Six-representation common inference | PASS | All 5,670 comparisons and 54 identity hypotheses, 9,999 permutations per arm; original PCA/NMF inference and all splits reproduced. All source/score hashes verified; complete outputs archived. |
| Signed PCA attribution | PASS | All 428,505 unique gene contributions and 105 partition sums recover score means; fixed fidelity gate passes (RMSE/SD=.000154, max AUC difference=.0000378). This is descriptive algebra, not mechanism or rank attribution. |
| Bounded full-size NMF warm restart | PASS | Numerical tolerance met after 380 further iterations; actual factors and all 6,615/63 hypotheses verified, archived and pushed. Stopping reference resets; no global-optimum claim. |
| Full target-RNA diagnostic | PASS | All 109 direct/alternative comparisons retained; three missing symbols left unscored. All 105 cross-diagnostic rows retained separately. |
| Historical K562 program | PASS | All 8,400 associations and 80 identity hypotheses retained. Original 945 direct-marker comparisons and parameters reproduced; actual scores/input/source hashes verified. |
| Context reference/depth sensitivity | PASS | Exact primary reproduced; all 4,200 reference, 1,050 depth and 8,400 technical-stratum records retained. Descriptive robust intersection, not replication. |
| Signature transcript/omission diagnostic | PASS | All 10,080 hypotheses and 1,680 new comparisons; original primary reproduced, joint BH independently recomputed. |
| Signature denominator sensitivity | PASS | All 945 original comparisons reproduced, all 11,025 old/new hypotheses retained with fixed weights; joint BH independently recomputed. No absolute-RNA claim. |
| Large-archive build and streamed roundtrip checks | PASS | Seven local archives, 1,696,336,750 bytes, exact member checksums verified. Upload remains blocked; this is not public deposition. |
| GitHub release-asset upload | FAIL, network | CONNECT to uploads.github.com returns HTTP403 before authentication. Empty release draft created; required host saved in configuration draft, which does not apply runtime policy. |
| `python src/make_norman_comparison_figure.py` twice | PASS | Both PNG/SVG byte-identical; rendered figure inspected. Counts and paired directions use complete archived outputs. |
| Robustness figure title and rendering | PASS | Title now describes displayed sensitivity settings; count calculations unchanged. Actual PNG/SVG rendering reproduced byte-identically twice. |
| `python src/build_manuscript_review.py --out-dir <fresh-directory>` | PASS | Actual PDF/Word build, source/generator/figure/output hashes, PDF content and Word CRC/five embedded figures checked. Final title, DOI links and figure pages inspected; review copies retain author-confirmation placeholders. [Render manifest](review/render_manifest.json). |
| Independent biological replication, novel mechanism or SAE superiority | UNVERIFIED | Computational seeds, gemgroups, annotations and cell-level p-values do not establish these claims. |

Nine separately executed additional suites make 24 distinct standalone suites exercised;
this is not a claim that a single 24-suite invocation was run. The 100-seed boundary
experiment is a research execution, not an additional unit-test suite. Source
fingerprints in run provenance identify the executed code revision; later
documentation commits do not change the underlying scientific outputs.
