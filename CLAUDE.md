# CLAUDE.md

Operating manual for this repository. Read this file before making changes.

## Start here

Read [README.md](README.md), [NORMAN_CONTINUATION_RESULTS.md](docs/NORMAN_CONTINUATION_RESULTS.md),
and [NORMAN_CPU_RUNS.md](docs/NORMAN_CPU_RUNS.md). Check existing CPU trainers and
`continue_norman_cpu.sh` and `archive_norman_cpu_results.py` workers before
starting another long job. All three 200-epoch primary SAE fits and one control-only
fit completed; read [NORMAN_COMPLETED_MODELS.md](docs/NORMAN_COMPLETED_MODELS.md).
Ordinary PCA and the capped NMF completed; a separate bounded NMF warm restart
met tolerance after 380 further iterations. The historical K562 program,
reference/depth, single-gene omission and composition checks also completed.
Read [NORMAN_K562_CONTEXT.md](docs/NORMAN_K562_CONTEXT.md) for positive conditional
associations and selection limits. Do not duplicate completed fits. Read
[LITERATURE_REVIEW.md](docs/LITERATURE_REVIEW.md) before proposing novelty.
Then read the
current-state and limitations sections of
[PROJECT_MASTER_SUMMARY.md](docs/PROJECT_MASTER_SUMMARY.md), and the latest
[LAB_NOTEBOOK.md](LAB_NOTEBOOK.md) entries. Use
[COMPONENT0_RESULTS.md](docs/COMPONENT0_RESULTS.md) and committed JSON/CSV artifacts
for Gate-0 numbers. Use [COMPONENT2_RESULTS.md](docs/COMPONENT2_RESULTS.md) for the
current causal interpretation and evidence gaps. Earlier strategy documents are
historical plans, not authority for current findings.

Cloud tasks already have isolated checkouts. Use the existing checkout; do not
create a Git worktree unless the user explicitly requests one. Preserve user files
and pre-existing changes. Keep new runs outside frozen result directories.

## Code delivery

Run code before presenting it as working. Before handing over a Python change:

1. Run `python3 -m py_compile` on changed Python files.
2. Run `python3 tests/test_v0b_logic.py`.
3. Run relevant tests and an actual synthetic or real-data execution of the changed
   path. A parser check or import is not an end-to-end run.
4. Report a verification table with commands and PASS, FAIL, or UNVERIFIED.

If full training or private data is unavailable, state exactly what was verified.
Never call an unexecuted command or a stand-in activation matrix a verified
biological analysis. Do not promise identical regenerated GPU outputs merely
because the random seed is fixed; record package versions, device, and determinism
settings and compare the regenerated results.

For code downloaded into another environment, pin the commit SHA and include a
version guard. Avoid commands that leave a stale copy under the expected filename.
Keep credentials out of code, outputs, saved instructions, and chat.

## Scientific reporting

- Lead with failed hypotheses and important limitations. Do not invent a discovery
  to satisfy a desired outcome.
- Distinguish no detected effect from an underpowered test or an analysis whose
  eligibility gates never fired.
- State uncertainty, seed counts, effect sizes, multiple-testing families, and
  what the permutation null actually randomizes.
- A frozen hash establishes file integrity. It is not independent registration,
  replication, or a reason to call exploratory analyses confirmatory.
- Gate 0 followed exploratory metric/marker development and then a version-locked
  final analysis. Causal amendments were informed by inspected results. Preserve
  that chronology in publications.
- Never silently revise frozen specifications or decision bundles. New analyses
  need dated amendments and separate outputs. Preserve historical notebook entries
  and append explicit corrections.
- Original Component 2 result JSONs disappeared with ephemeral GPU instances.
  Historical tables and notebook reports remain, but exact causal regeneration
  has not been demonstrated. Archive reruns before treating these as submission
  evidence.
- Separate association, perturbation responsiveness, direct regulatory mechanism,
  and biological novelty. One does not establish the next.

## Data and implementation facts

- The project is computational; there is no wet-bench validation.
- Paul15's depth-adjusted source matrix has raw-scale values. The supported
  `preprocess_paul15.py` policy applies `log1p` when appropriate; do not copy the
  superseded blanket instruction that Paul15 is already log transformed. Human
  Setty marrow counts are library-size normalized and then logged by their own
  preprocessing script. Follow each source-specific pipeline.
- Gate-0 matrices verify against `preprocess_report.json`: mouse MD5
  `184ea5ddfca23db2bbb343b311688c45`, human
  `627961897d52a50c7a81f5bbf7e04097`. MD5
  `60183a17983c8b977d036e0f3a58da61` belongs to superseded Phase 1 data.
- The canonical V0b implementation is `src/v0b_module_definitions.py`. Do not run
  the historical Drive notebooks; their names do not reliably identify versions.
- Core V0b logic deliberately works without Torch. Do not introduce top-level
  Torch calls into that module's pure-logic paths.
- SAE feature indices have no shared identity across seeds. Match or annotate
  features separately within each seed; never compare feature index alone.
- Use `V0B_MARKERS_JSON=config/human_markers.json` for human structural analyses.
- All historical causal SAEs trained on **all cells**, without perturbation labels.
  `--train-on all` preserves that mode. Control-only training is a sensitivity
  analysis whose current synthetic training test recovers 0/8 planted regulators,
  compared with 8/8 in the all-cell mode. This is a known limitation, not evidence
  that every possible control-only method fails.
- In the causal representation ladder, PCA/NMF transform the **input to an SAE**.
  Every arm still trains an SAE. Do not describe that ladder as direct comparison
  against plain PCA/NMF features.
- Causal cross-fitting separates matching and testing cells; the SAE itself still
  saw all cells. Any statement about generalization must disclose this scope.
- For Norman data use raw counts. Pertpy's `norman_2019()` is processed;
  `norman_2019_raw()` is the intended source. `prep_norman.py` guards against
  non-count input. GATA1 is absent from the Norman perturbation library; it was
  not tested and failed there.
- Apply the analysis-specific seed rule. Frozen Gate 0 uses a 60% decision
  majority over 10 seeds; program modulation specifies seeds 0,1,2 and stable
  modulators in at least 2/3. Do not substitute a generic seed convention.

## Current evidence

Gate 0 is executed and frozen in two datasets: mouse Paul15, 2,730 cells by 2,012
genes; human Setty CD34+ marrow, 4,142 by 2,032. Each has 10 SAE seeds. The original
asymmetric-modularity hypothesis is supported in 0/10 SAE seeds in each dataset.
Across five decomposition families and two datasets, it is supported only for
mouse NMF. Marker leave-one-out changes 0/34 mouse and 0/36 human verdicts.
The reversal seen in SAE counts is method-contingent, not established biological
ground truth. Do not alter either frozen bundle to fix prose rounding.

Historical Replogle K562/RPE1 results identify primarily housekeeping/chromatin
perturbation responses. The recorded effect-size control removed no seed-0 calls;
the former 10-fold essential-enrichment claim used the wrong background. The
recorded v4 comparison is 6/51 Hart-core-essential-tagged TFs versus 14/111 other
tested TFs; the latter group is not a validated nonessential reference set.
Column specificity failed to establish the proposed fix. The expression-versus-
foundation-model embedding comparison was never run, and a published 6.2% rate is
not a matched null for this pipeline.

Norman grounding and single-TF regulon recovery did not support their intended
claims. Fresh full-size training and historical-program follow-ups are complete.
The historically selected eight-gene K562 program has 63/61/61 all-cell SAE calls,
57 direct-score calls and 38 surviving references, depth and every gene omission;
PCA also recovers it. Historical selection used the same cells, so these are
conditional exploratory associations, not independent confirmation, a new
regulator mechanism or an SAE advantage. Read the complete archived outputs and
NORMAN_K562_CONTEXT.md before proposing additional work. Model/data archives are
verified locally; GitHub asset upload remains blocked by the network proxy.

## Working style

Be concise and direct. Use no em-dashes. Research prose should be formal and
evidence-limited. Explain corrections plainly. Execute the work, then distinguish
what passed, what failed, and what remains unverified.
