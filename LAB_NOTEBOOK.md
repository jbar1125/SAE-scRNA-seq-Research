# LAB NOTEBOOK

Day-to-day research log for the SAE / scRNA-seq hematopoiesis project (Jacob
Barzideh). Newest entries at the BOTTOM within each date; dates newest-first would
break the narrative, so this reads top-to-bottom chronologically. Each entry: what
I set out to do, what I actually did, the decisions and reasoning, bugs, results,
and the commit SHAs that captured it. Semi-brief on purpose. The git history is the
byte-level record; this is the human-readable why.

Conventions: "SUPPORTED / NOT SUPPORTED" always refers to the ORIGINAL asymmetric-
modularity hypothesis (erythroid distributed, granulocyte unified). Commit SHAs are
short. All code lives in `src/`, docs in `docs/`, specs/markers in `config/`.

---

## 2026-06-29 — Project seeded from the Phase-1 handoff

**Goal.** Get the project into version control and make it legible to Claude Code.

- Added `PROJECT_HANDOFF.md` (551-line Phase-1-closeout snapshot: thesis, the
  asymmetric-modularity finding + the 6.2% causal caveat, full V0-V15 plan, ISEF
  paperwork, references).
- Decision: keep the handoff intact as a historical record; corrections go in a
  separate audit rather than editing the handoff in place.

Commits: `4ac09e9`.

---

## 2026-06-30 — Canonical V0b, the notebook trap, first honest run

**Goal.** Replace the ambiguous Drive notebooks with one trustworthy implementation
and start testing the headline claim.

- Wrote `v0b_module_definitions.py` (canonical v2) to SUPERSEDE the three Drive
  notebooks. Audit finding: the canonically-named, most-recently-modified
  `v0b_module_definitions.ipynb` is actually **v1** (uses `CROSS_SEED_THRESHOLD=3`
  instead of per-seed winner-take-all, pads markers with genes absent from the
  2000-HVG set, omits `prob_19Lymph`, scores with the biased `H_ery/H_gran` ratio),
  and no `module_assignments_v0b.csv` exists anywhere — so V0b had never actually
  been completed. This is the "notebook trap."
- Added `PROJECT_AUDIT.md` (authoritative correction layer) and synthetic logic
  tests (`tests/test_v0b_logic.py`).
- **Bug shipped and fixed same day:** a lazy-import refactor removed the top-level
  `torch` import but `main()` still called `torch.manual_seed` → NameError. It went
  out without being run. Fix `ca1a746`. This is the incident that motivated the
  hard "run before you send" rule.
- Rewrote `CLAUDE.md` with non-negotiable code-verification + honesty rules and a
  smart summary.
- **First real V0b (v2) run: asymmetric modularity NOT SUPPORTED.** But this is an
  ARTIFACT, not a clean negative: the hypergeometric + winner-take-all gate is
  structurally blind to distributed signal (diluted markers never clear q<0.05
  across ~896 mostly-null tests). Recorded as such — do not report "erythroid is
  unified" from a blind test.
- Wrote `v0b_v3_loading.py`: continuous decoder-loading test, no significance gate,
  so distributed signal can register. First v3 run came back "SUPPORTED" — but the
  diagnostic showed that verdict is a count/abundance artifact, not a finding.

Commits: `b320bb9`, `ca1a746`, `ed79c6a`, `b1aaa5e`, `70315d4`, `874f244`.

**Thought.** Two "results" in one day (v2 NOT SUPPORTED, v3 SUPPORTED) that are both
untrustworthy for opposite reasons. The lesson crystallizing: the metric itself is
the experiment. Need a pre-registered, control-referenced, null-calibrated metric
before any verdict means anything.

---

## 2026-07-01 — The marker bug, the reversal, the reframe, Component 0 built

**Goal.** Fix the data/compute confounds, build the real metric, and situate the
project honestly in the field.

Data/compute upgrade:
- Overcomplete SAEs (512 latents) + marker-aware preprocessing (`train_sae.py`,
  `preprocess_paul15.py`): force-include the full marker union so no lineage is
  analyzed with its effector program removed.
- Log-transform inconsistency found: the handoff calls Paul15 "pre-log-transformed",
  but the source loads RAW (max ~168). Replaced the hard guard with a `--log1p auto`
  policy; recorded as audit inconsistency A11. (This contradicts a CLAUDE.md
  "non-negotiable"; the audit wins on facts.)

**The marker-curation bug (my error, owned).** Earlier runs used wrong globin
symbols (`Hba-a1`/`Hbb-bs`, absent from Paul15) and omitted `Alas2`. So the
erythroid effector program was invisible for a bookkeeping reason and erythroid
looked falsely weak/distributed. Fixed to Paul15's real symbols (`Hba-a2`,
`Hbb-b1`) + `Alas2` + `Ermap`. After the fix the hemoglobin program is the
STRONGEST feature in the data.

**The reversal.** With corrected markers, control-referenced counts give erythroid
1 dominant program (unified), granulocyte 3 (distributed) — the REVERSE of the v1
headline. Recorded in audit section J. Explicitly NOT frozen yet (still L0 too high,
n=5, mouse-only, metric not pre-registered).

The metric, finalized:
- `v0b_v3_1_decision.py`: Progenitor-referenced baseline, abundance-matched
  permutation null + BH-FDR, effect-size floor, 60%-majority decision rule.
- `baselines_nmf_pca.py` (PCA/NMF participation ratio + same v3.1 modularity),
  `marker_sensitivity.py` (leave-one-marker-out).
- Pre-registration freeze: `config/preregistration_spec.json` + `.sha256`
  (`fc342829`) + `PREREGISTRATION.md`. Froze the metric/markers, NOT any verdict.
- Loader hardening: seed auto-detection (any number of `sae_seed*.pt`), latent-dim
  inferred from the checkpoint, dual-source MD5 (locked original OR the marker-aware
  matrix's own recorded hash).

Positioning:
- `STRATEGY_AND_POSITIONING.md`: after verifying the Kendiukhov 2026 abstract, RETRACTED
  the overstated "SAEs on single-cell is novel" claim (it is now published:
  Kendiukhov 2026 arXiv 2603.02952, CytoSAE, the 2025-26 wave). Reframed the novelty
  to the artifact-vs-signal framework + the causal CRISPRi test.
- `phase2_research_plan_v6.md`: replan around the corrected direction and the causal
  head-to-head vs the 6.2% null.
- Added the h5 fallback loader so real Paul15 can be obtained where `sc.datasets`
  is blocked.

Commits: `1f39d29`, `84f2b22`, `f5fd742`, `f4ca4c1`, `a43dc32`, `bae8c7a`,
`edd9502`, `312ea5b`, `91a023c`, `f3bc7c8`, `ea259ff`, `35173e4`, `15d00e4`,
`1a0866a`, `41e4b81`.

**Thought.** The project's real contribution is now clear: not "we ran SAEs on
cells" but "here is how to tell whether an SAE program is real." Every prior number
is retired; the only ones that count come from the pre-registered metric on a
properly-regularized run.

---

## 2026-07-02 — Gate 0 EXECUTED end to end (mouse + human)

**Goal.** Actually run Component 0 on real data, in-container, and freeze.

Environment reality: no GPU in the container. SAEs are ~2M params on 2-4k cells, so
CPU trains a 10-seed run in ~8-12 min. Data obtained via git-LFS mirrors (Paul15
from a scanpy mirror; human Setty CD34+ marrow from the Palantir repo), both
hash-verified at load.

Sparsity retune (the QC lesson): L0 depends on training LENGTH, not just L1. At 150
epochs l1=3.0 gave L0 34.6; at the real 300 epochs the same L1 collapsed to L0 ~17
(too sparse). Retuned at 300 epochs → **l1=0.8 → L0 37 (band center)** for mouse.

**Mouse primary result (frozen).** 10-seed 512-latent SAE, L0 mean 35.75 (all in
20-50 band). Original claim NOT SUPPORTED (0/10 seeds). Reverse is seed-robust:
granulocyte 3 real above-control programs in 10/10 seeds; erythroid 1 dominant
hemoglobin program (Ery_Effector strength 5.07). Progenitor control perm_q 0.085
(non-significant → null calibrated). Leave-one-marker-out: 0/34 flips.

**Baselines — the real finding: method-dependence.** Participation ratio SAE 196 vs
PCA 14 vs NMF 5. Running the SAME v3.1 metric on each decomposition gives DIFFERENT
answers: PCA resolves ~nothing, NMF says erythroid-distributed ("SUPPORTED", the
original direction, but parts-based artifact), SAE says granulocyte-distributed. No
unsupervised decomposition recovers a method-invariant program structure.

**Human replication (frozen).** Setty 2019 CD34+ marrow, 4142x2032, CELLxGENE
pipeline (normalize_total + log1p) with 1:1 human ortholog markers
(`config/human_markers.json`) injected via a new `V0B_MARKERS_JSON` override in
`v0b_module_definitions.py` (defaults to mouse; byte-for-byte identical when unset).
Retuned to l1=0.8 → L0 32.5. Claim NOT SUPPORTED (0/10). Granulocyte-more-distributed
direction REPLICATES. Honest species difference: the human globin program is present
but NOT significant vs the abundance-matched null (Setty is CD34+ progenitor-selected,
so terminal erythroblasts are underrepresented). Human baselines: PCA agrees with SAE
(granulocyte), NMF detects nothing. Leave-one-marker-out: 0/36 flips. Control perm_q
0.748.

**Cross-species x cross-method synthesis.** The original claim survives in exactly 1
of 6 method x species combinations (mouse NMF, the most artifact-prone). The SAE
granulocyte-distributed direction is the only claim holding in both species, but NMF
disagrees in both — so even the reversal is method-contingent. This is the
artifact-vs-signal thesis, demonstrated across 2 species and 3 decompositions.

Frozen decision bundles (properly-regularized 10-seed runs, negative-plus-method-
dependence, NOT a positive claim from a broken run):
- mouse `172861417a923a83705812a96e5ad1e555055e968323289879447e53beb84396`
- human `cc8159c882cf326bca5bbbc7f74d110d4a81338ec0ba07e208eb935bc8f729d3`

Recorded in `docs/COMPONENT0_RESULTS.md` (every number) and audit section K.
**Gate 0 is complete.** PR #1 updated to reflect the executed scope.

Commits: `b36fc0b`, `2181170`, `f562c66`, `1b8495a`.

**Thought.** The headline the project started with is dead, and that is fine — what
replaces it (rigorous method-dependence, cross-species) is stronger and more honest
for a competition. Everything structural is now settled; whether the programs are
CAUSALLY real is Component 2's job, not Gate 0's.

---

## 2026-07-04 — Repository reorganization + this notebook

**Goal.** Make the project a clean, competition- and mentor-legible repo: everything
accessible, prior versions preserved, nothing lost, and a real dated lab notebook.

- Reorganized the flat root into: `src/` (8 code files), `tests/`, `docs/` (7
  narrative/reference docs), `config/` (frozen specs + `human_markers.json`).
  Executed run dirs `data_g0/` (mouse) and `data_g0_human/` (human) kept as-is so
  the results record stays path-accurate. All moves via `git mv` (history preserved).
- Updated every execution-affecting reference: tests now add `src/` to the path
  (plus repo-root for the one cross-test import), `preprocess_human_marrow.py`
  default markers → `config/human_markers.json`, CLAUDE.md links → `docs/`, the
  Colab runner download URLs → `src/`.
- Added `README.md` (entry point / project map), this `LAB_NOTEBOOK.md`, and
  `docs/VERSION_HISTORY.md` (the metric-lineage record: v2 gated → v3 continuous →
  v3.1 pre-registered).
- Verified nothing broke: `py_compile src/*.py` clean; all 4 logic tests PASS from
  the new layout; a live `src/v0b_v3_1_decision.py` run on `data_g0/` reproduces the
  NOT-SUPPORTED verdict.

**Open threads / next.**
- Optional cheap hardening: SCENIC/pySCENIC as a 5th baseline; TRRUST/ChIP-Atlas
  cross-reference of the top SAE programs.
- Consider a second human dataset (non-CD34-selected) to test whether the human
  globin program's non-significance is purely a composition effect.
- Component 2 (centerpiece): Replogle CRISPRi causal head-to-head vs the 6.2%
  embedding-space null. Needs a GPU environment (Remote Control on own hardware, or
  Colab) for the larger models/atlases.

### Gate-0 hardening: robustness battery (human-gated autonomous pipeline)

Operating mode agreed: I self-drive the mechanical work and stop-and-ask only at
real forks; phone push on finish/decision. Queue (all CPU-tractable; Component 2
waits for GPU): R1 decision-parameter sensitivity, R2 L0-band sensitivity, R3 more
seeds, R4 HVG-count sensitivity. Then gated: SCENIC 5th baseline, TRRUST/ChIP-Atlas
cross-ref, second human dataset.

**R1 decision-parameter sensitivity — DONE (`src/robustness_sweeps.py`).** Sweeps the
two "arbitrary" knobs of the frozen metric (EXCESS_FLOOR in {0.05,0.10,0.15,0.20} x
PERM_Q in {0.01,0.05,0.10}, 12 configs) on the existing 10-seed checkpoints; the
expensive permutation null is computed once and the grid evaluated on top. Result:
in BOTH species, **0/12 configs support the original claim and granulocyte >=
erythroid real-program count at every point** -> the NOT-SUPPORTED verdict is not a
knob-tuning artifact. Artifacts: `data_g0{,_human}/robustness/`.

**Autonomy update.** Per `/goal`, mode switched to fully autonomous: no halting to
ask; forks are resolved with a conservative default and logged here as "Decision
(auto)". Only irreversible/outward-facing actions or hard capability limits block.

**Decision (auto) — compute documentation.** Added `docs/COMPUTE.md` recording the
CPU-vs-GPU reasoning (Gate 0 is CPU-fine; GPU starts at Component 2; the AMD RX 5700
XT is not practical for the PyTorch/ROCm stack). Reason: the user asked how GPU work
is handled; documenting it prevents repeated confusion for mentors.

**Decision (auto) — external DB reachability (for later units).** Probed: TRRUST
direct host grnpedia.org is proxy-BLOCKED (000); pySCENIC installs from PyPI but its
cisTarget motif-ranking databases (aertslab) are large and expected to be blocked.
Plan when those units come up: try a GitHub mirror for TRRUST; if blocked, defer the
DB cross-reference and note it as an external-resource limit rather than fake it. For
SCENIC: if cisTarget DBs are unreachable, run a GRNBoost2-only co-expression baseline
labeled honestly as such (NOT full SCENIC), or defer. No fabricated regulons.

**R2 L0-band sensitivity — DONE.** Retrained mouse (5 seeds/point) at L1 spanning the
QC band: l1=1.0 -> L0 27.7 (gran3/ery2), l1=0.7 -> L0 42.2 (gran3/ery2), l1=0.6 ->
L0 52.0 (just outside band, gran3/ery1). **NOT SUPPORTED at every point and
granulocyte >= erythroid throughout** -> the verdict is not specific to L0=37, it
holds across the whole band and at the boundary. Artifact:
`data_g0/robustness/r2_l0_band_summary.json` (retrain dirs gitignored).

**Event: PR #1 MERGED** (main = commit 9a05f53). Per repo protocol, restarted the
designated branch from merged main (clean; e30ab93 is an ancestor of main, no force
needed) and continue robustness work as a NEW PR. All Gate-0 + reorg + R1 work is now
in main.

Next: R3 (20 seeds, stability of n=10) then R4 (HVG-count sensitivity, uses the raw
Paul15 h5 on disk).

**Decision (auto) — TRRUST ground-truth cross-reference (gated unit, resolved).**
Direct host grnpedia.org is blocked, but the canonical TRRUST v2 human file is
mirrored across many GitHub repos (same sha a701172d); fetched
`config/trrust/trrust_rawdata.human.tsv` (795 TFs; all lineage TFs GATA1/KLF1/TAL1/
CEBPA/CEBPE/RUNX1 present). Mouse TRRUST mirrors 404'd; decision: use human TRRUST
and uppercase mouse symbols (edges ~conserved), and report the mouse cross-ref only
as an underpowered cross-species check, not a finding.

**TRRUST cross-ref — DONE (`src/trrust_crossref.py`).** Pre-specified test: for each
lineage, pick the single SAE feature most enriched for that lineage's TF submodule,
then hypergeometric-test whether the lineage TFs' TRRUST targets are over-represented
in that feature's top-30 genes. **Human: erythroid TF targets enriched in the
erythroid-TF feature in 100% of seeds (median p=6.4e-5); granulocyte in 50% (median
p=0.13).** So SAE features recover real regulatory structure, and, consistently with
the main finding, the UNIFIED erythroid program captures its regulon tightly while
the DISTRIBUTED granulocyte program does so less from any single feature. Mouse
cross-ref (uppercased vs human TRRUST) is null (0-10%), reported as an underpowered
cross-species artifact (small HVG panel x human DB), not evidence. Artifacts:
`data_g0{,_human}/trrust/trrust_crossref.json`.

**Decision (auto) — SCENIC 5th baseline: DEFERRED (external-resource limit).** Probed:
pySCENIC/arboreto/ctxcore install from PyPI and the aertslab TF list is on GitHub,
but the cisTarget motif-ranking databases (resources.aertslab.org, ~GB .feather
files) are proxy-BLOCKED (000). Full SCENIC (motif pruning + AUCell) is therefore
infeasible in-container. Rejected a GRNBoost2-only surrogate because it yields a
TF x target adjacency, not a cell x component decomposition, so it does not slot into
the PCA/NMF/SAE modularity comparison and would over-claim the "SCENIC" label. The
TRRUST cross-reference already provides the curated-regulatory-database validation
SCENIC would contribute, so deferring loses little. To run when unblocked (Colab):
`pip install pyscenic; pyscenic grn`, then `pyscenic ctx` with the hg38 mc9nr feather
DBs, then `pyscenic aucell`, then run v3.1 modularity on the AUCell cell x regulon
matrix as a 5th decomposition.

**Decision (auto) — second human dataset: DEFERRED (documented next step).** Purpose
would be to test whether the human globin program's non-significance is a CD34-
selection artifact, by adding a non-CD34-selected marrow (total BMMC incl. mature
erythroblasts). Deferred because it needs both a reachable non-CD34 h5ad (uncertain
via the proxy; GEO/figshare are blocked) and a full CPU train cycle (fragile given
the container's ~10-15 min restart cadence). The mechanism is already stated in
COMPONENT0_RESULTS and the cross-species DIRECTION already replicated on the one human
set, so this is a hardening add, not a gap in the core result. Acquisition approach
when possible: a GitHub-LFS-mirrored BMMC h5ad (like the Setty acquisition) or scanpy
built-ins on a networked box, then reuse `preprocess_human_marrow.py`.

**R3 (efficient) — running.** After finding the 20-seed retrain re-trained the 10
frozen seeds needlessly and was too long for the restart cadence, added `--seed-start`
to train_sae and launched only the NEW seeds 10-19 (combined with the frozen 0-9 for
a 20-seed v3.1).

**R3 seed stability — DONE.** 20 seeds (10 frozen + 10 new; all L0 in-band 35-37.5).
20-seed v3.1 reproduces the 10-seed verdict EXACTLY: **gran_n_real = 3 in 20/20
seeds; ery_n_real = 1 in ten seeds and 2 in ten (median 1.5); granulocyte > erythroid
in 20/20; 0/20 support the original pattern; NOT SUPPORTED.** Doubling the seed count
changes nothing -> n=10 was sufficient and the granulocyte-more-distributed direction
is highly stable. Artifact: `data_g0/robustness/r3_20seed_summary.json`.

**R4 HVG-count sensitivity — DONE.** Re-preprocessed mouse at n-hvg 1000 (1023 genes)
and 3000 (3006 genes), 5 seeds each, l1=0.8. **NOT SUPPORTED at 1000/2000/3000 HVGs;
granulocyte 3 / erythroid 1-2 throughout.** Honest caveat: at fixed l1=0.8, L0 leaves
the band at 1000 genes (~90, too dense) and 3000 genes (~18, too sparse) because L0
depends on gene count; a per-count l1 retune would recenter L0, but the verdict is
stable without it. Artifact: `data_g0/robustness/r4_hvg_summary.json`.

**GRN baseline — DONE (`src/grn_baseline.py`), replacing the deferred SCENIC.** Full
SCENIC stays infeasible (cisTarget DBs blocked), but rather than only deferring, built
an executable co-expression GRN as a 4th decomposition family: for each aertslab TF in
the panel, regulon-loading = |Pearson corr| with each gene -> a genes x TF matrix run
through the IDENTICAL v3.1 modularity pipeline (labeled honestly as a correlation
regulon, NOT motif-pruned SCENIC). Result: **mouse ery=2/gran=2 (a TIE) - a fourth,
distinct answer** vs PCA (1/1), NMF (2.5/2 "supported"), SAE (1.5/3); human ery=1/gran=2
(agrees with SAE). So FOUR decomposition families give FOUR different mouse modularity
answers - the method-dependence thesis is now demonstrated across 4 method families,
not 3. Artifacts: `data_g0{,_human}/grn/grn_baseline.json`; TF lists in
`config/tf_lists/`.

**Autonomy continuation (later same day).** Pipeline pushed further per `/goal`.

**Decision (auto) — genuine resource blocks, re-probed and documented (not skipped).**
- Second human dataset: api.github.com is gated to the one in-scope repo ("access not
  enabled"); the Palantir repo exposes only `marrow_sample_scseq_counts.h5ad` (all other
  candidate filenames 404). No reachable non-CD34 human marrow set. BLOCKED until the
  user enables another data repo or link-shares a file. Not faked.
- Mouse TRRUST: re-probed 3 GitHub mirrors for `trrust_rawdata.mouse.tsv`; all 404
  (consistent with the earlier probe). Mouse curated-regulon cross-ref stays deferred.
- Full SCENIC (motif-pruned): cisTarget DBs remain blocked; the co-expression GRN
  (`src/grn_baseline.py`, committed a04d1ec) is the honest tractable substitute and is
  the 4th decomposition family in the method-dependence table.

**PR #3 check-in (cron 2fb86799).** State open/draft, mergeable_state clean, no CI
configured, no review comments. No action; cron stays armed until merge/close.

**Human-arm robustness (R2h/R3h/R4h) — DONE.** Extended L0-band, 20-seed stability, and
HVG-count robustness to the human Setty arm. All confirm the human verdict: NOT
SUPPORTED and granulocyte >= erythroid across the L0 band (21.4/32.5/41.5), across 20
seeds, and across HVG counts (1000/2000/3000). Human R4 note: at 1000 genes both
lineages collapse to 1 real program (L0 ~98, far out of band) - a tie, still not the
claim. Summaries in `data_g0_human/robustness/`; hardening doc R2-R4 rows now report
both species. The robustness battery is now symmetric across mouse and human.

**ICA 5th decomposition family — DONE, with an honesty save.** Added `src/ica_baseline.py`
(FastICA components -> same v3.1 modularity). First run at rank 512 (SAE-matched) did
NOT converge (mouse 0/3, human 1/3 seeds) and the non-converged human run FALSELY read
"SUPPORTED" (ery=3/gran=2). Caught it: non-converged ICA components are unreliable, so I
would not present them as a result. Diagnosed 512-component ICA as ill-posed; ICA
converges 3/3 at rank 50. Re-ran at rank 50 (converged): **mouse ery=1/gran=3, human
ery=2/gran=3, both NOT SUPPORTED**, agreeing with the SAE granulocyte direction. The
rank differs from the other families (512), stated openly. Method-dependence now spans
5 families x 2 species; the original claim survives in exactly 1 of 10 combos (mouse
NMF). This is the honesty rule working as intended: a convergence artifact nearly
entered the record and was rejected.

**Result figures — DONE (`src/make_figures.py`).** Two publication-quality figures
generated from the committed JSON artifacts (regenerable, no retraining): Fig 1 the
5-family x 2-species method-dependence (diverging bars; the lone warm bar = mouse NMF);
Fig 2 the mouse robustness verdict stable across every stress test. CVD-safe diverging
palette, sign encoded by position + color, direct value labels. PNG+SVG under
`docs/figures/`, referenced from `COMPONENT0_HARDENING.md`.

**Reproducibility pass + a real self-correction — DONE.** Ran an end-to-end verification
(compile, all 4 logic tests, frozen-bundle SHAs, figure determinism) and wrote
`docs/REPRODUCIBILITY.md` with the verification table.
- **Bug caught: pre-registration hash mismatch.** `config/preregistration_spec.sha256`
  and every doc record `fc342829...`, but `sha256sum config/preregistration_spec.json`
  returned `6398c5ed...`. Diagnosed: `fc342829` is the hash of the canonical
  `json.dumps(spec, indent=2)`; the committed file had a trailing newline, so its
  raw-byte hash differed and the documented verification FAILED. This has been true
  since the freeze commit `35173e4`. Fixed by normalizing the file to the exact
  canonical bytes; spec CONTENT is byte-identical to the freeze (verified `json.load`
  equality), so nothing frozen changed, only the trailing byte. `sha256sum` now prints
  `fc342829...`. PREREGISTRATION.md verification block corrected (was a placeholder).
- **Figures made deterministic** (`svg.hashsalt` + drop embedded Date): two runs now
  produce byte-identical PNG/SVG.
- Frozen decision bundles re-verified: mouse `172861...`, human `cc8159c8...` still match.

## 2026-07-05 — THE PIVOT: from a rigorous null to a positive causal head-to-head

User's honest call: "It is not a strong result as of now and will not win anything...
make this a million times better... with a stronger result." Correct. Ran an
autonomous literature + reasoning pass and rebuilt the project's centerpiece.

**Diagnosis (docs/ELEVATION_PLAN.md).** The current result is a negative + "methods
disagree" — hygiene, not a discovery; no causal/biological payload; defensive framing;
crowded substrate. It is a foundation, not a result.

**The reframe, grounded in the literature.** The key fact: SAEs on single-cell
FOUNDATION MODELS (Geneformer/scGPT) encode organized knowledge but "minimal
regulatory logic" - only 3/48 TFs (6.2%, ->10.4% multi-tissue) are causally grounded
vs Replogle CRISPRi (arXiv 2603.02952). Nobody has asked: is that a limit of SAEs or of
EMBEDDINGS? New central question + falsifiable hypothesis: SAEs trained directly on
EXPRESSION recover causal regulatory logic that foundation models discard, beating the
6-10% ceiling. Novelty confirmed by search (only tangential prior work). K562 (Replogle)
keeps the myeloid framing and gives real causal ground truth.

**Triviality trap identified + defeated (ELEVATION_PLAN 3b).** "Expression sees the
genes, so of course it wins" - defeated by: train on control only; match by program
with the KD gene EXCLUDED (program-level, not target self-drop); require perturbation-
specific suppression AND a confident single-target-program match. The residual
expression advantage then IS the scientific claim (causal logic is legible in
expression, compressed away in embeddings).

**Built + validated in-container (CPU, no GPU):**
- `src/causal_grounding.py` - the metric. Unit-tested (`tests/test_causal_grounding.py`):
  6/6 planted regulators grounded, null + global(non-specific) + background rejected,
  0 false positives under no causal structure. Two real design bugs caught and fixed by
  testing (response-based matching was circular; column-specificity alone couldn't
  reject a global stressor -> added program-match confidence).
- `src/causal_pipeline.py` - the RunPod head-to-head (TopK SAE on control cells ->
  encode -> grounding; --rep expression vs --rep embedding). `--synthetic` self-test
  trains a real SAE at realistic scale and recovers 6-7/8 regulators with ~1/250 false
  positives -> full pipeline (learning + grounding) validated end to end.
- Pre-registered the causal metric: `config/causal_grounding_spec.json` (+ .sha256
  `b7d28fae`), frozen BEFORE any real run.
- `docs/RUNPOD_EXECUTION.md` - turnkey steps for the GPU run (data via pertpy, both
  arms, compare, budget). Honest about what is CPU-verified vs what runs on RunPod.

**Negative result, reported plainly.** Prototyped a causally-supervised SAE (add a
decoder-side reconstruction of held-out perturbation DE signatures to the objective).
It did NOT improve held-out causal grounding on synthetic - it hurt (plain 4/4 vs
causal 0/4 on separable programs; ~tie-to-worse on entangled). Removed the unvalidated
code (`src/causal_sae.py`) rather than ship a broken "method"; shelved pending a better
formulation. The head-to-head does not depend on it. This is the honesty rule working:
I tried to make it "more," it did not validate, so I say so.

**Next (needs RunPod GPU):** run Arm B (expression) + Arm A (scGPT embedding) on
Replogle; the head-to-head number; if positive, build the causally-grounded ATLAS (a
low-risk formatting of the grounded programs + their perturbation certificates) and
reformat the writeup around the causal result.

**Hardening battery COMPLETE.** R1-R4 + TRRUST + GRN + ICA all confirm the frozen Gate-0
verdict (no method resurrects the original claim) and extend the method-dependence
result.
Summary written to `docs/COMPONENT0_HARDENING.md`. Terminal state for the CPU
container: SCENIC + second human dataset deferred (documented), Component 2 specified
(`docs/COMPONENT2_PLAN.md`, GPU-blocked). PR #3 carries all of it.

---

## 2026-07-18 — First real Replogle run returned 0/1789; deep-dive found TWO real metric bugs

**Goal.** Run Arm B (expression-space causal grounding) on the real Replogle K562
Perturb-seq on a rented GPU (vast.ai, RTX A5000). It came back **0/1789 grounded =
0.000, BELOW the shuffled null (0.009)**. Per the honesty rule the first sentence is:
this was NOT a clean negative, it was a broken metric. I did the requested deep dive
and found two independent bugs, both now fixed and validated on synthetic before any
re-run on real data. Nothing was frozen or reported from the broken run.

**Setup problems fixed first (not the science).**
- Repo was private; clone on the rented box failed auth. A GitHub PAT was pasted in
  plaintext -> told the user to REVOKE it immediately; made the repo public so no token
  ever touches a rented machine. (Security: never put a token on rented compute.)
- OOM thrash: dense 310385 x 8563 = ~10 GB blew up the 24 GB box (memory red, GPU + CPU
  idle). Added `--n-hvg` to `load_perturbseq`: keep the top-N HVGs UNION the tested
  perturbation genes, subset BEFORE densify. Verified on a tiny AnnData.

**Bug 1 — SD-standardization is invalid for sparse TopK activations.** The v1 metric
required a standardized activation drop >= 0.25 SD. A TopK feature fires in only ~k/L of
cells, so its activation SD is large vs its mean; even a FULL shutoff is only ~0.1 SD --
unreachable by any floor -> forced 0 grounded, mechanically. Diagnosed by printing the
per-feature drop distribution. Fix: the suppression test is now a one-sided Mann-Whitney
U (KD cells vs control cells) reported as AUC = U/(n_kd*n_ctrl); the floor is on AUC
(<= 0.45), a rank statistic with no scale assumption -> robust to sparse/zero-inflated
activations.

**Bug 2 — double-dipping.** The same cells both selected the matched feature and tested
its suppression, so a diffuse/background perturbation could be "grounded" by whichever
feature happened to dip in exactly those cells. Fix: **cross-fitting** -- split each
perturbation's cells A|B, match on A, test suppression on held-out B. This is the fix
that made the oracle test reject background/global perturbations cleanly.

**Wrong gene set (contributing).** The first run tested ALL ~1789 perturbations,
including housekeeping/essential genes whose knockdown is non-specific. Added
`--tf-list config/tf_lists/hs_hgnc_tfs.txt` (1,839 human TFs) to restrict scoring to
transcription factors, where a program-level causal signature is meaningful.

**Validation (CPU, synthetic, before touching real data again).**
- Oracle unit test `tests/test_causal_grounding.py`: PASS -- 6/6 planted regulators
  grounded, null + global + background all rejected, null-world calibration 1/13.
- Pipeline self-test `--synthetic`: 8/8 true regulators grounded, 0 false positives,
  shuffled null mean 0.000 (p=0.048). PIPELINE SELF-TEST PASSED.
- Real-vs-shuffle separation holds across background sizes (n_bg=40: 5/8 true, 0 false;
  n_bg=250: 8/8 true, 0 false) vs shuffle 0.000.

**Pre-registration amended, not silently overwritten.** The metric changed materially
(SD-standardized drop -> cross-fit Mann-Whitney AUC), so `config/causal_grounding_spec.json`
now carries a dated `amendments` block recording exactly what changed and why, and the
frozen SHA moved `b7d28fae -> b766d21c`. Legitimate because NO valid real result existed
under v1 -- the amendment is pre-registered before the first valid real run, per the
no-frozen-artifact-from-a-broken-run rule.

**Next.** User re-runs Arm B on vast.ai off the fixed metric with
`--tf-list config/tf_lists/hs_hgnc_tfs.txt --n-hvg 2000` after `git pull`. Then the real
expression-space grounding number, then Arm A (scGPT embeddings) for the head-to-head.

### Later 2026-07-18 — first real Component 2 number (single seed) + an essential-gene confound

**Result (PRELIMINARY, single seed, NOT frozen).** Arm B ran on real Replogle K562 (GPU,
cuda): loaded 310385 cells x 2136 genes (top-2000 HVG UNION perturbed TF genes), 162
testable TFs. **grounded 22/162 = 0.136**, shuffled null mean 0.001, p=0.048 (p is at the
20-shuffle resolution floor 1/21; 0/20 shuffles reached the observed value -> clean
separation, not weak signal). SAE latent 2048, k 32. diag: pass floor(auc<=0.45) 35,
mw_q<0.05 24, pass match 162/162.

**Honest read (result vs artifact).** The 13.6% is REAL (survives label shuffle) but the
grounded SET is the wrong flavor. The 22 are dominated by essential / general-transcription
machinery, replication/repair, RNA-processing, and even non-TFs: TBP, TAF7, GTF2A2 (basal
apparatus), TFDP1, NELFB, HINFP, CXXC1, E4F1, POLD2, RFC2, SMUG1, DMAP1, TIMELESS, SFPQ,
THOC2, NCBP2, ILF2, SRP9 (SRP, not a TF), CSNK2B (kinase), PTPMT1 (phosphatase); only ~2
(HSF1, MAX) are clean sequence-specific TFs with a defined program. NO lineage regulators
(GATA1, SPI1, TAL1, RUNX1, KLF1, CEBPA) grounded. Interpretation: knocking down an
essential gene collapses transcription broadly, a "global cell health" feature drops
specifically in those cells -> causal in a trivial sense, NOT regulatory-specific. The
triviality trap in a new coat (KD gene's own mRNA is excluded, but "the KD makes the cell
sick" is not). The match gate passing 162/162 confirms it is not enforcing program
specificity. So the honest headline is "13.6% perturbation-specific grounding, but
essential-gene-dominated; regulatory specificity NOT yet shown," NOT "beat the 6.2%
benchmark."

**Three fixes queued (implement next session; do NOT ship a specificity claim without
them):**
1. **Stricter TF list.** `hs_hgnc_tfs.txt` includes basal machinery (TBP, TAFs, GTFs). Swap
   to the Lambert et al. 2018 "The Human Transcription Factors" sequence-specific-DBD list
   (~1600 regulators) and re-measure. Needs the curated Lambert list committed to
   `config/tf_lists/`.
2. **A specificity filter that bites.** Require the matched feature's program to be a
   coherent, GO-enriched module (not a diffuse global-lowness feature), and tighten the
   match gate so it actually filters (currently 162/162 pass -> non-binding). Add a
   program-coherence score to `causal_grounding.py`.
3. **Essential-gene control.** Score grounding separately for essential vs non-essential
   (DepMap common-essential list). If the rate is high only for essentials, THAT is the
   finding (and a different one). Needs the DepMap essential list committed.

**Running overnight (2026-07-18):** multi-seed (seeds 0-4, n-shuffle 50) for seed-stability
+ the >=3/5 rule; then a latent x k sensitivity sweep (512/1024/2048 x 16/32/64, n-shuffle 0)
to show 13.6% is not a knob-tuning artifact. Arm A (scGPT embeddings) is the real
head-to-head and is a separate, hands-on setup (not an unattended job).

### Overnight results in (2026-07-18) -> `docs/COMPONENT2_RESULTS.md`

- **Seeds (2048/32):** mean 0.195, sd 0.039, range 0.136-0.259 across seeds 0-4. All >>
  the ~0 shuffle null -> grounding is SEED-ROBUST; the exact rate is noisy (~20% CV).
- **Seed-stable set (>=3/5): 26 perturbations.** Still ~24/26 essential/machinery/non-TF,
  BUT **GATA1 grounds robustly** (master erythroid TF, in an erythroleukemia line) -> the
  metric DOES catch a real lineage regulator; regulatory signal is buried, not absent.
  Positive-control (Tier 0.4) partially satisfied.
- **Sensitivity sweep (seed 0):** rate spans **0.136-0.395 from hyperparameters alone**.
  Dictionary size is INVERSELY related (512->0.395, 1024->0.284, 2048->0.136): smaller/
  broader dictionaries ground MORE, consistent with the essential-gene/global-state
  confound. Higher rate != better. -> the absolute rate is not quotable standalone; only
  the matched-hyperparameter head-to-head is.
- Logged PRELIMINARY, NOT frozen (Arm A missing, essential-gene-confounded, one cell line).
  GPU box then destroyed; raw JSONs regenerable from pinned seeds+code.
