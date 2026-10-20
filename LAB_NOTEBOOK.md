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


### Even later 2026-07-18 — the causal fix FAILED on real data (honest negative)

Ran v4 (column specificity) on real Replogle, 5 seeds. It did NOT work and actively hurt:
seed-stable set (22) is entirely housekeeping/machinery, GATA1 fell from >=3/5 (v3) to 2/5
(v4 removed the one lineage regulator), all other master TFs (SPI1/TAL1/RUNX1/KLF1/CEBPA/
GFI1/GATA2/MYB) 0/5. Root cause: distinct essential KDs hit distinct stress features (gate
passes each), and GATA1's K562 erythroid program has multiple regulators (gate cuts GATA1
as "non-specific"). Rolled back to default OFF (spec v5, SHA 98104b2f).

Also corrected my own stat: the "~10x essential enrichment (p~4e-7)" used the wrong
background. The tag-breakdown shows essential 11.8% vs non-essential 12.6% grounding -- no
preference. The pool is ~31% essential because Replogle-essential IS an essential-gene
library. Grounded set is housekeeping-heavy because the POOL is.

Honest bottom line: on this dataset the metric tracks effect-DETECTABILITY, not regulatory
specificity; no master hematopoietic regulator robustly grounds. Real negative for the
"expression recovers regulatory logic" hypothesis here. Not ruled out: the head-to-head
(does embedding do worse?) and a dataset that actually perturbs lineage TFs. Third
self-correction of the day; recorded, not buried.

### 2026-07-18 (M1 result) — effect-size is NOT the confound; grounding is real but housekeeping

Ran raw vs effect-size-controlled (M1) on real Replogle. The control changed NOTHING: raw
22/162 == controlled 22/162, 0 dropped, 0 rescued. (Own error: framed M1 as a "rescue" -- but
an AND-filter can only remove, never add.) The dropped-zero IS the finding: all 22 hits pass
the effect-size control (they suppress more than effect-matched peers), so they are NOT
effect-size artifacts. Diagnostic weak: rho(effect,suppression)=0.28, mean effect grounded
3.81 vs 3.07. This refutes my earlier "grounding just tracks effect-detectability" claim --
overstated. GATA1 effsize_p=0.435 (not even effect-relative-specific) -> the lineage signal
isn't in this dataset, not hidden by a confound. Honest reframe: expression-space SAE
genuinely grounds ~22 causally-real, effect-controlled programs, but they're proliferation/
housekeeping (the biology present in K562-essential), not lineage. Method works; dataset lacks
the lineage biology. Next: (1) head-to-head vs embeddings, (2) a lineage-perturbation dataset.

### 2026-07-18 (representation ladder) — SAE not special; grounding rate flat across expression decompositions

Ran the metric with the SAE trained on expression vs PCA-128 vs NMF-128 of the same Replogle
data. Rate is FLAT: expression 0.121, PCA 0.130, NMF 0.111 -- statistically the same (PCA even
marginally higher). So the SAE recovers no more causal grounding than plain matrix
factorization; it is not special for this metric. The grounded SETS are representation-
dependent though: seed0 expression∩PCA overlap = 6/22 (~2x chance), expression∩NMF = 2 -- a
weak reproducible core + a large method-dependent tail. Convergent with Gate 0's structural
method-dependence: the causal-grounding SET is representation-dependent while the RATE is a
~13% near-constant. Kills "expression SAEs are special"; refines the claim to "expression
SPACE grounds ~13% robust to decomposition -- does that beat the embedding ~6%?" -> the
embedding head-to-head is now the pivotal remaining experiment.

### 2026-07-18 (RPE1 replication) — ladder-flat replicates; rate cell-line-dependent; still housekeeping

Ran RPE1 (second cell line) expression seeds 0-2 + pca/nmf seed 0. (1) Ladder FLAT again:
expression 0.076 ~ pca 0.070 ~ nmf 0.076 -> "SAE not special vs matrix factorization"
replicates in a second cell line (robust). (2) Rate is LOWER than K562 (7.6% vs 13%) ->
the grounding rate is cell-line-dependent, not a universal constant. (3) Grounded set again
housekeeping/chromatin/ribosomal (RPL35, HMGB1, HMGA1, SMARCC1, GABPB1, THAP11, ZNF720),
with ZERO overlap with K562's housekeeping set and NO lineage TFs (all 0/3). Consistent
meta-finding: the metric grounds whichever essential/proliferation programs dominate a given
dataset -- real but dataset-specific, not-special, never hematopoietic. Expression side now
thoroughly characterized (2 cell lines x 3 representations x effect-controlled). Diminishing
returns on more expression runs; the pivotal unrun experiment is the embedding arm.

### 2026-08-10 (Norman CRISPRa) — first lineage-relevant grounding; PROVISIONAL, null still missing

Pivoted off the Replogle CRISPRi screens (essential-gene libraries that do not perturb lineage
TFs -- a wrong-substrate negative, see 07-18 entries) to Norman 2019 CRISPRa (K562, ~111k
cells, dual-guide overexpression via pertpy). This is the RIGHT substrate: it actually
overexpresses hematopoietic master TFs (KLF1, CEBPA, SPI1, ETS2, ...). Ran the metric with
`--direction up` (activation: match +DE program, MW "greater", floor auc>=0.55). Data prep:
parse `guide_identity` -> perturbation label (drop NegCtrl tokens; single/dual guides), map
Ensembl var_names to gene symbols so perturbations are testable (the earlier ZeroDivisionError
was 0 testable perts from the Ensembl/symbol namespace mismatch; fixed).

FIRST genuinely encouraging signal, and it is biologically coherent -- but READ THE CAVEATS:
- UP-direction grounding rate mean ~0.121 across seeds. KLF1 (erythroid master) grounded 2/3
  seeds; ETS2 grounded 3/3. Seed-stable grounded set = KLF1, SLC4A1, ETS2, RUNX1T1, CBFA2T3 --
  erythroid / hematopoietic, NOT the housekeeping/ribosomal set that dominated every Replogle
  run. That qualitative shift (lineage TFs instead of machinery) is the actual news.
- DOWN-direction on the same data gave ~0.363 -- this is the documented TopK feature-competition
  artifact (overexpression pushes competitor features out of top-k so they look "suppressed");
  `down` is the WRONG sign for OE and must not be reported as the result.
- Myeloid TFs (SPI1/CEBPA) did NOT ground consistently -- consistent with K562's erythroid/
  megakaryocytic bias (they are not strongly myeloid).

WHAT IS NOT YET ESTABLISHED (do not overclaim, do not freeze):
- NO shuffle null yet (this run used --n-shuffle 0). The rate ~0.121 is the same MAGNITUDE as
  Replogle's ~13%; without a label-shuffle null on this dataset we cannot yet say the up-rate is
  above chance. The claim rests on the SET being lineage-relevant, not on the rate.
- Have NOT yet confirmed that KLF1's grounded feature IS the erythroid/hemoglobin program. Added
  a `top_program_genes` readout (commits 4b335e6, 566c334) to surface the matched feature's top
  gene weights inline; the re-run must show KLF1 -> HBB/HBA/GYPA/ALAS2/SLC4A1 for this to be a
  clean, hand-verifiable positive rather than a coincidental match.
- Single dataset, single cell line. Not frozen. Provisional pending: (1) --n-shuffle 50 null,
  (2) top_program_genes confirmation, (3) replication.

Next run: Norman `--direction up --n-shuffle 50` reusing the cached h5ad, reading out
top_program_genes for KLF1/ETS2. If the null collapses below 0.121 AND KLF1's feature is
hemoglobin, this is the project's first real positive discovery.

### 2026-08-10 (Norman null IS IN) — CORRECTION: up-rate does NOT beat shuffle (p=0.90). Negative.

The shuffle null from the prior entry is now run (--n-shuffle 50, direction up, latent2048/k32,
n_hvg 4000, 102 testable perturbations, on GPU). It came back NEGATIVE and I am correcting the
provisional-positive framing of the previous entry:
- seed 0: grounded 15/102 = 0.147 | shuffle null mean 0.175 | empirical p = 0.902. The real
  up-rate is BELOW the label-shuffle null. Seeds 1/2 (no null): 0.127, 0.196 -- all inside the
  null band. On the pre-registered "rate above chance" test, Norman CRISPRa FAILS, same verdict
  as Replogle CRISPRi (different mechanism).

Distinguish result from artifact (both are present, per CLAUDE.md):
- REAL (null-independent): the SAE recovers a pristine, SEED-STABLE erythroid/hemoglobin feature
  -- top genes HBZ, ALAS2, HBG1/HBG2, HBA1, SLC25A37 (mitoferrin), GYPB, BLVRB, FTH1 -- identical
  program across all 3 seeds (feat 537/630/1982). Genuine, hand-verifiable biology. The SAE
  decomposition is finding true structure.
- ARTIFACT (why the null is not beaten): grounding is PROMISCUOUS. The hemoglobin feature is
  grounded by CBL/PTPN12/UBASH3B/ZBTB25 (ubiquitin ligases / phosphatases -- signaling brakes),
  NOT by GATA1 (doesn't ground) or KLF1 (grounds a separate, self-dominated KLF1+HBZ feature). A
  second cluster (MAP2K3/MAP2K6/CSRNP1) grounds one stress/proteostasis feature
  (PRDX1/LDHB/PSMA4/UBB/HSP90AA1) -- the same housekeeping grounding as Replogle.

Mechanism of the null (itself the finding): K562 has ~2 dominant transcriptional axes (erythroid,
stress); the SAE isolates them cleanly; the matching step then grounds almost ANY grouping --
real or shuffled -- to one of those dominant features. So the metric's per-dataset false-positive
rate is ~17%, swamping the ~15% real rate. The binary "grounded rate vs shuffle" is the wrong
instrument once a couple of axes dominate the variance.

Honest bottom line: TWO nulls now -- CRISPRi (wrong substrate: essential screens, no lineage TFs)
and CRISPRa (right substrate, but dominant-axis promiscuity). The broad claim "expression-space
SAE grounding beats chance" is NOT supported in either. NOT a positive discovery. Do not freeze.

Interesting honest lead (NOT yet tested, must pre-register before running): the brake genes that
DO ground the erythroid feature (CBL/UBASH3B/PTPN12) are negative regulators of KIT/EPOR signaling
whose overexpression is known to push K562 toward erythroid differentiation -- so their activation
of the erythroid feature may be REAL directed biology, not noise. A DIRECTED, pre-registered test
-- "do erythroid-promoting perturbations activate the erythroid SAE feature above a matched null,
while non-erythroid ones do not" -- is the one remaining honest shot at a positive, on the ONE
feature we can name. It may also come back null; register it first, report it as exploratory.

### 2026-08-10 (regulon recovery v2) — STOP. Real negative: SAE features do not match curated regulons

v2 (spec sha c653f6d1) fixed both defects the v1 diagnostic found. Both fixes worked
MECHANICALLY: 133 hub targets dropped, panel 4047 -> 5791 genes, tested pool 10 -> 14 TFs.
The result got WORSE: 0/14 recovered, shuffle null 0.00, p=1.000, all 3 seeds, and every
enrichQ = 1.0 (v1 at least reached 0.086). Reason is visible in the numbers: force-keeping
targets made regulons BIGGER (JUN 57->102, SPI1 32->54, CEBPA 22->35) while hits stayed 2-4,
and a bigger regulon needs proportionally more hits to clear the hypergeometric bar.

Per the v2 spec's own pre-commitment ("if v2 again returns an uninformative or negative
headline, report it as such and stop"), this line of inquiry is CLOSED. No v3.

But v2 is NOT merely uninformative -- it upgrades to a REAL, QUANTIFIED negative. Simulating
the actual selection procedure (best of 2048 random top-50 gene sets, 5791-gene panel,
src/regulon_chance_baseline.py) gives the chance expectation, and the observed overlaps sit AT
OR BELOW it for every regulon size:

  regulon targets:   13    18    32    35    40    54   102
  observed overlap:   2     2     2     3     3     4     4
  chance best-of-2048: 2.28  2.67  3.20  3.29  3.58  4.12  5.34

FINDING: SAE features do not concentrate curated single-TF regulons any better than the best of
2048 random gene sets of the same size. The causal question was unanswerable with this ground
truth because the feature/regulon alignment does not exist at single-TF resolution.

Scope it correctly (do NOT overstate):
- This is about SINGLE-TF CURATED-REGULON alignment, not about whether SAE features are
  biologically meaningful. The same SAEs recover a textbook-clean erythroid/hemoglobin program
  (HBZ, ALAS2, HBG1/2, HBA1, SLC25A37, GYPB, BLVRB), seed-stable across all runs.
- It does NOT contradict Gate 0's TRRUST positive (human erythroid targets enriched, 100% of
  seeds, p=6.4e-5). That test was at LINEAGE granularity (union of a lineage's TF targets,
  512 features, marrow data). Precise synthesis: regulon alignment is detectable at coarse
  lineage level and ABSENT at single-TF resolution with per-feature multiple-testing control.
- Likely cause (unverified, stated as hypothesis): TRRUST edges are curated across many cell
  types/conditions, while SAE features are co-expression modules in one cell line. Co-expression
  modules need not equal curated direct-target sets. This is a known gap in the field, not a bug.

TALLY of the Component-2 causal program, all executed, all honest:
  1. Replogle CRISPRi (K562-essential, RPE1): grounds ~13%/7.6%, but housekeeping only, no
     lineage TFs. WRONG SUBSTRATE (essential-gene libraries do not perturb lineage TFs).
  2. Norman CRISPRa rate test: up-rate 0.147 vs shuffle null 0.175, p=0.902. NOT above chance;
     dominant-axis promiscuity (K562 has ~2 axes; ~any labeling grounds).
  3. Norman regulon recovery v1: uninformative (test never fired; panel truncation + hub hijack).
  4. Norman regulon recovery v2: REAL NEGATIVE, quantified above.
The broad hypothesis "expression-space SAE features are causally grounded, single-perturbation
resolution" is NOT SUPPORTED across 2 screens, 2 cell lines, 2 metric families. That is the
honest state of the project. Do not freeze any result table from these runs.

WHAT SURVIVES (the asset worth consolidating): a seed-stable, hand-verifiable erythroid program
recovered unsupervised, reproducible across seeds AND representations, which a specific
non-obvious set of perturbations (CBL, UBASH3B, PTPN12 -- all negative regulators of KIT/EPOR
signaling) causally activates. Plus a rigorous negative-result methodology with a diagnostic
that distinguishes "no signal" from "test never fired" -- which is exactly what caught v1.

### 2026-08-10 (self-correction) — the shuffle null in causal_pipeline is CONSERVATIVE, not neutral

Building the program-modulation oracle surfaced a null-design issue that also applies to
`shuffle_control` in src/causal_pipeline.py, so it is recorded here rather than buried.

Both permute labels among PERTURBED cells only, controls fixed. That destroys perturbation
IDENTITY but PRESERVES the perturbed-vs-control contrast. In the program-modulation oracle this
inflated the null badly: null mean 14.4 modulators vs 5 real, p=1.000, on a synthetic world where
the ground truth is 4 real modulators. Mechanism: spreading strongly-activating cells across every
shuffled group lifts every group above the fixed control.

Consequence for the reported Norman rate result (up-rate 0.147 vs null 0.175, p=0.902): that null
retains generic perturbation-vs-control signal, so it is a STRINGENT null, not a neutral one. What
that p actually establishes: the real rate does not exceed a null that itself contains generic
perturbation response. It does NOT establish that grounding is at or below true chance. The
earlier entry's headline ("NOT above chance") should be read with that caveat. I am NOT retracting
the negative -- the regulon-recovery chance-level finding is independent of any shuffle and stands
on its own -- but the p=0.902 is weaker evidence than it was presented as, and that is my error.

Fix going forward: program_modulation defaults to shuffle_mode="full" (permute ALL labels incl.
control, so group membership is independent of biology), with mode="perturbed" retained to ask the
different, stricter question. Oracle-verified: full-mode gives p=0.032 on the specific world where
perturbed-mode gave p=1.000.

---

## 2026-10-18: full-project audit for the master summary; four latent defects found

**Goal.** Compile `docs/PROJECT_MASTER_SUMMARY.md` (everything, start to finish). Doing it from the
repo rather than from memory surfaced four real problems. Leading with them.

**1. The SAE was never trained on control cells only.** `config/causal_grounding_spec.json` says
"trained on control (non-targeting) cells only; SAE never sees a perturbation" and ELEVATION_PLAN 3b
lists it as triviality-trap defense #1. Every caller (`causal_pipeline.run`, `regulon_recovery`,
`program_modulation`) passed the full matrix to `train_topk_sae`, so every real run fit on ALL
cells (labels unseen). First in-container run of the pipeline self-test (torch from PyPI): all-cell
training grounds 8/8 planted regulators with 0 false positives; control-only grounds 0/8 (signal
present, matched AUC ~0.40 vs 0.31, best mw_q ~0.078, just over the FDR cut). So the pre-registered
design was never validated, and in this synthetic world it is underpowered. Handling, without
changing any past number: `train_rows()` + `--train-on {all,control}` (default `all`, results
byte-identical to before); result JSONs now record `sae.trained_on`; dated erratum in the causal
spec (`ee5dacb7` -> `aae1a9d9`); program-modulation spec v1.1 (`b396834d` -> `ddb29bcd`) makes the
all-cell primary explicit and pre-declares control-only as a sensitivity run; erratum notes in
`causal_grounding.py`, `causal_pipeline.py`, COMPONENT2_RESULTS, ELEVATION_PLAN,
METHODOLOGY_ADDENDUM_OE; `tests/test_pipeline_training.py`. Interpretation of the Replogle
positive narrows: features learned from perturbed + control cells respond specifically; NOT
"programs learned from natural variation behave causally".

**2. My 2026-08-20 re-download instruction was wrong and unrun.** I gave `pt.dt.norman_2019()`
without running it. Reading both pertpy files' schemas over HTTP: `norman_2019` is PROCESSED
(111,255 cells, log-normalized X, raw counts only in `layers['counts']`, symbol var_names,
`perturbation_name`); `norman_2019_raw` is raw (111,445 cells x 33,694, integer counts, Ensembl
var_names + `gene_symbols`, `guide_identity`, no `perturbation_name`). The August schema matches
the raw file, so every reported Norman result was normalized exactly once. Following my August
instruction would have double-normalized the program-modulation run. Fix: `ensure_raw_counts` in
`prep_norman.py` (passes raw, repairs from `layers['counts']`, refuses otherwise) +
`tests/test_prep_norman.py`. Verified on the real raw file: 11,855 controls, 57,831 single-gene
cells, 41,759 dual-gene cells, 102 testable (identical to August).

**3. GATA1 is not in the Norman library.** 0 cells (RUNX1 too; RUNX1T1 is present). The
2026-08-10 entry's "NOT [grounded] by GATA1 (doesn't ground)" implied GATA1 was tested and failed.
It was never perturbed.

**4. The regulon chance comparison was worded too strongly.** Per TF and seed from the pasted v2
output, against the simulated best-of-2048 distribution: 8/42 cells above its mean, 0/42 above its
95th percentile (SPI1 seed 0 exactly at it). "At or below chance for every regulon size" was too
strong; "within the chance distribution everywhere" is right. Conclusion unchanged. Also: the v2
amendment expected the force-keep to make KLF1/LYL1/CBFA2T3/MEIS1 testable; it did not, because
their TRRUST regulons are tiny (KLF1: 2 targets).

**Verified today.** All 9 test suites pass in a fresh container; all 4 spec hashes and both Gate-0
bundles verify; the Norman data path, loader and all three Perturb-seq CLIs run end to end on real
data (subsample, tiny SAE; plumbing only, numbers not interpreted). Timing on this 4-core CPU: an
all-cell 2048-latent fit is 2.0 h/seed, control-only 0.22 h/seed, so the decisive run needs a GPU.

**Also established.** Regeneron STS 2027 closes 2026-11-05, 8:00 PM ET (18 days). The
program-modulation test remains unrun; the verified recipe is in the master summary section 18.1.

---

## 2026-10-18 — continuation: environment, publication claims, and modulation design audit

The requested continuation began from full local history (93 commits at `a08db6d`),
not an additional shallow-history recovery. The CPU environment was rebuilt in an
isolated virtualenv outside the checkout. All nine original standalone suites,
all four specification hashes, and both frozen Gate-0 bundle hashes passed. The
installation was rerun successfully and reusable install/start instructions were
saved in the cloud configuration draft. No original frozen artifact was changed.

**Additional claim correction.** `build_representation` supplies expression/PCA/NMF
inputs, but every causal arm then trains a TopK SAE. The historical ladder cannot
support “SAE no better than plain PCA/NMF.” It reports similar observed rates for
SAEs trained on different inputs. Original causal JSONs from destroyed instances
are not present; notebook tables cannot establish exact regeneration. Current
README, operating instructions, abstract, strategies and historical-plan notices
now disclose those limits. The historical 10-fold essential-enrichment figure was
replaced with the recorded v4 6/51 core-essential-tagged vs 14/111 other-tested-TF
comparison, with descriptive Wilson intervals and an unarchived-source warning.
The 111 comparator is not the Hart NEGv1 validated nonessential set.

**New methods finding, not biological discovery.** In 20 independent simulated
seeds, all 20 perturbations identically activating one axis passed the two old
quantitative program-modulation gates in 20/20 cases. A distinct four-responder
world also passed in 20/20, and a no-effect world passed in 0/20. The benchmark
used one planned program and 99 permutations; it did not score the manual
biological-coherence criterion. Direct synthetic activations were used, not
trained SAEs. An absent-marker world exposed a separate bug: every feature has
marker-enrichment p=1, the selector returns -1, and the caller scores the last
column; all 20 seeds passed the quantitative gates despite no marker match.
The current code explicitly leaves absent/weak matches unscored and retains all
planned hypotheses. Exact 95% simulation-frequency intervals are 83.16–100% for
20/20 and 0–16.84% for 0/20. Full raw outputs and source fingerprints are archived
under `data_program_modulation_audit/`, including the original executed generator.

**Separate inference.** A controls-fixed, optionally stratified permutation test
of between-identity rank variation now applies Holm correction over all planned
programs. It rejected all 20 distinct-responder worlds and none of the 20 generic
or 20 no-effect worlds; absent matches receive p=1. A deployed-default calibration
with nine stipulated matched features, 999 permutations and 200 independent null
worlds produced 6/200 familywise rejections (3.0%; exact 95% interval 1.11–6.42%).
This is finite synthetic calibration, not proof for arbitrary cell/guide designs.
A rejection means conditional observed-cell identity heterogeneity; variable
stress severity can also produce it. It does not establish pathway specificity,
biological replication, novelty or SAE superiority. The old spec remains intact;
`config/program_modulation_validation_plan.json` records the next-run scope.

**Executed software path.** The revised modulation CLI ran on synthetic raw-count
AnnData with actual CPU SAE fitting, strata and hashed provenance. It found no
significant named-program matches and correctly left them unscored. Invalid
zero-epoch runs and existing output paths are rejected before data loading. This
is an implementation smoke test only, not the planned full Norman analysis.

**External prerequisites remain.** The Norman public-data hostname and GitHub API
returned proxy-policy 403 responses. Both additions were saved in the environment
configuration draft; draft saving does not apply runtime access. Existing native
Git authentication works and correction commit `c7fe97e` was pushed to the same
branch. No NVIDIA GPU is attached. No full real-data modulation result or new
biological discovery is claimed. See `docs/PROGRAM_MODULATION_AUDIT.md` and
`docs/MANUSCRIPT_DRAFT.md`; the manuscript has explicit citation, evidence,
independent-validation, author/disclosure and submission tasks remaining.

## 2026-10-18 — CPU continuation: archived real Norman baseline

After the usage interruption, the user explicitly chose this CPU environment and
asked to continue. The public Norman file had been downloaded and prepared.
Checkpoint/marker-baseline work left in the workspace was inspected and tested,
then committed as `d1ab30c`; checkpoint-enabled CLI and explicit unique-count
preparation were committed as `ebdc0b0`. The raw source observed SHA is
`32cc993abad41336eb0f3635f97e802e13fe5721dd9277c1c55e7c13853782e0`.
Original files were preserved. All 34 duplicated symbol columns were summed before
normalization with exact per-cell count conservation; none overlaps the 59 markers
or 105 single-target symbols. The unique input SHA is
`a09744916427b08d75dade996d059f717e4812463292a1e5a4e5450589cf2c6a`.

The full-data marker baseline completed: 111,445 cells, 5,926 discovery and 5,929
inference controls split within eight gemgroups, all nine panels, all 105 single
perturbations, BH over 945 comparisons. Target omission left six TF-panel, 37
heme, 25 membrane, four granulocyte-TF and three progenitor associations. Three
panels were unscorable (<3 eligible markers); Cycling had no positive association.
All six adequately covered panels nonetheless had conditional identity
heterogeneity (raw permutation floor .001, Holm .009), so this does not isolate
an erythroid-specific mechanism.

A subsequent, explicitly exploratory plan retained the original scores/split and
examined discovery-defined RNA-depth bins plus every negative-control reference.
Joint reference correction covered 3,780 comparisons. Twenty-nine heme, 20
membrane, three TF-panel, four granulocyte-TF and two progenitor associations
survived the original and all specified sensitivities. Reference intersections
are descriptive robustness, not independent replication; depth adjustment can
remove perturbation-mediated effects. Original primary computations reproduced
exactly before sensitivity testing.

KLF1's TF-panel AUC .615 becomes .425 after its own transcript is excluded, while
heme .626 and membrane .580 associations remain. CEBPA/CEBPE omission tests fall
below minimum marker coverage and are unscorable, not clean negatives. CBL and
PTPN12 survive the heme/membrane sensitivities. UBASH3B's heme AUC .571 falls to
.548 against one negative-control reference, below the fixed .55 floor, so the
former three-candidate story is not confirmed wholesale. These are marker-score
associations, not novel mechanisms or an SAE advantage.

Full JSONs, source/data fingerprints, exact plans, a 945-row comparison table and
figures are archived under `data_norman_continuation/`. All 15 standalone suites
passed, with zero skips; all four original spec hashes and both frozen Gate-0
bundles still match. Figure outputs are byte-identical across consecutive runs
and were visually checked.

The 2048-feature, k32, 200-epoch all-cell SAE seed 0 is running on four CPU cores,
with complete epoch-boundary model/optimizer/RNG/input-hash checkpoints. Seeds 1
and 2 are queued sequentially by `src/continue_norman_cpu.sh`; no completed new
SAE result is claimed. Exact uninterrupted/resumed small-model equality and a
real synthetic checkpoint/activation CLI were tested. An explicit interim plan
allows a separate copied checkpoint at or after epoch 20 solely to validate
full-data inference, without selecting a favorable training epoch.

GitHub REST access and dataset access work. The PR title/body were updated and
verified after the old gh-pr GraphQL route hit a deprecated Projects field.
Crossref/EuropePMC requests still returned403; their domains were added to the
configuration draft for citation review. Updated startup instructions preserve
CPU resume behavior and distinguish completed marker results from unfinished
SAE training. Primary-text citation review, independent biological validation,
confirmed disclosures and publication submission remain outstanding.

## 2026-10-18 — Interim full-shape inference and durable CPU continuation

The predeclared readiness check copied the epoch-20 seed-0 checkpoint from the
200-epoch parent without stopping or modifying the primary trainer. Full-data
inference reconstructed the copied model without optimizer steps. Four panels
matched (Ery_TF, Ery_Heme, Ery_Effector, Cycling), each with per-seed identity Holm
p=.009; the other five retained p=1. The heme feature had seven exploratory
modulators and background percentile 86. No matched feature reached percentile 90.
All outcomes and snapshot hashes are archived under
`data_norman_continuation/sae_interim_epoch 20/`; none is a completed primary result.
This observation does not change the 200-epoch target or select a favorable epoch.

HBB was detected in only three discovery controls, failing the direct baseline's
minimum-ten rule. SAE matching checks presence rather than that detection rule.
Thus matching Ery_Effector while its direct score is unscorable cannot establish
an SAE advantage. This eligibility mismatch is explicitly recorded.

A separate completed-result watcher validates configuration, all nine hypotheses,
permutation execution and model/activation hashes before archiving compact seed
outputs and publishing scoped evidence commits. All three completed seeds are
required for the joint 27-test Holm table. Tests reject partial or incompatible
runs and exercise recovery after interrupted archive publication, without claiming
the still-future 200-epoch real outputs have completed. This sixteenth suite passed
separately after the full fifteen-suite run. The sequential training worker and
primary fit remain active; no duplicate trainer was launched.

Citation APIs recovered. Crossref metadata and the EuropePMC Norman abstract were
retrieved. The original Norman abstract already reports CBL/CNN1 synergy driving
erythroid differentiation, so CBL-associated erythroid recovery cannot alone be
presented as novel biology. Legitimate free full-text hosts were added to the
network draft; access and a fuller claim-by-claim literature check remain pending.

## 2026-10-18 — Ordinary PCA comparison and primary-source citation checks

The ordinary PCA plan and sign/heldout-selection tests were committed before
viewing full-data results (`5fbf6e8`). The fit uses 64 components, randomized SVD,
three power iterations, seed 0, all cells and the identical 4,081-gene SAE input.
Both signed orientations enter the annotation search and its correction; no SAE
is trained. Annotation uses the existing discovery-control split, and associations
use heldout controls, gemgroup strata and joint BH over 945 tests. This split is for annotation
and association: all-cell representation fitting still sees every expression row.
The small test's inference-control independence holds conditional on fixed scores,
not as an independent representation-training holdout. Unequal model capacities
and differing procedures preclude a superiority claim from the response counts.
The real-data command was launched with one CPU thread alongside the main SAE.

Bibliographic metadata for nine references was verified with Crossref. The Norman
abstract was read through EuropePMC; its free HTML remained inaccessible and its
non-OA REST full-text route was unavailable. Two other open-access XMLs were
retrieved: Wang et al. (2024) on FLI1/UBASH3A/B and Kendiukhov (2026) on SAE model-internal
circuits. Targeted methods/results/limitations sections support the statements in
`docs/LITERATURE_REVIEW.md`; the experiments were not reproduced. Seven first-page
query logs are archived as bibliographic metadata, without treating hit counts or
missing hits as novelty evidence. The manuscript now includes verified Norman and
TRRUST references and explicitly distinguishes prior biology from this analysis.

The declared direct PCA analysis completed. All 128 signed orientations were
eligible, explaining35.77% of expression variance across 64 components. Ery_TF
matched component 12 with 28 positive associations and Cycling matched component 3
with 52; both had identity Holm p=.009. Seven panels were unmatched and stayed
unscored, including Ery_Effector at adjusted annotation p=0.0561. The TF-matched
feature's leading genes include globins and mitochondrial transcripts as well
as LMO2, illustrating why the panel name does not make different scores equivalent.
Cycling's 52 PCA associations versus zero for its direct marker score reflect
different measurements, not a method ranking. All 945 comparisons and full outputs
are archived; the exact existing split and fitted-model/source hashes verified.
The direct-PCA suite passed separately, bringing distinct exercised suites to 17.

## 2026-10-18 — Fixed identity-test severity and dependence audit

Plan and generator were committed before outcomes (`7a3c220`). One hundred seeds
crossed with six analyses retain600 unique complete records, nine hypotheses each.
With 999 permutations, identical generic responses gave 4/100 familywise rejections
(exact 95% interval 1.10–9.93%). A single shared factor with unequal generic severity
gave all-nine rejections in 100/100: correct identity rank heterogeneity does not
establish separate pathway regulation. Pooled shuffling under complete identity/
batch confounding rejected in 100/100; within-batch conditioning correctly found
no exchangeable identities and ran no permutations in every case.

A separate population-null world had four independent clusters per perturbation
and ten correlated cells each. Invalid cell shuffling rejected in 100/100 worlds;
analysis of 40 independent cluster means (minimum-two-unit eligibility) rejected
in 8/100, exact 95% interval 3.52–15.16%. This wide interval includes 5%; it does not
prove calibration or establish Norman's experimental units. Paired worlds share
some draws, and stipulated matches do not calibrate feature selection. The test
source, running SAE settings and thresholds remain unchanged.

All raw outcomes, exact plans, sources, versions, logs and summaries are archived
under `data_program_modulation_audit/boundary_audit/`. Counts were recomputed from
all 600 records and artifact/source hashes verified. The standalone PNG/SVG figure
renders archived data and is byte-identical across two generations; its initially
overlapping count labels were corrected before committing. The manuscript now
reports these limitations quantitatively rather than assuming them away.

Later nominally successful HTML requests still did not retrieve the Norman text:
PMC served a browser challenge, EuropePMC an application shell, and its explicit
no-JavaScript route failed. Those are not full-text review. The abstract and the
two successfully retrieved open-access comparison papers remain the verified
review sources; no challenge or access restriction was bypassed.

## 2026-10-18 — Direct NMF launch and verification index

The direct NMF settings were committed before inspecting outcomes: rank 64,
NNDSVDa initialization, multiplicative updates, Frobenius loss, tolerance 1e-4,
maximum 200 iterations and seed 0. It uses the identical nonnegative matrix,
existing discovery/inference split and unchanged factor-analysis routine. Negative
orientations have no positive discovery controls and are ineligible for this
nonnegative score representation. Warnings and actual iterations will be retained;
capacity/search differences preclude a superiority claim.

The initial launch was stopped before its loaded-data log or any fit outcomes,
to correct a timing field that would otherwise include downstream analysis in
fit duration. The log/path was preserved. The corrected run uses a fresh v2
directory, `/workspace/research-runs/norman-direct-nmf-v2/`, and one CPU thread.
The primary SAE was not interrupted. Actual fitted cell factors and gene components
will be retained outside Git with hashes, rather than replacing joint-fit factors
with a potentially different later NMF transform.

The original V0b logic suite and four spec/two bundle hashes passed again after
the continuation changes. The 15-suite raw log and separate archiver/PCA logs are
now retained in Git. `docs/CONTINUATION_VERIFICATION.md` provides commands and scoped
PASS/UNVERIFIED statuses, explicitly distinguishing the two still-running fits
from completed evidence and from biological claims.

## 2026-10-18 — Dual-marker additivity screen and scale follow-up

Metadata showed 131 dual labels (41,759 cells, minimum 54 per pair), with both
single arms adequately sampled for every pair. The plan and contrast tests were
committed before dual expression outcomes (`93d7f2c`). The original marker-fit
parameters and exact control split reproduced; all 9 × 131 comparisons remained
in the Holm family. The four-group, within-gemgroup Welch contrast uses fixed
dual-count weights and discovery-control score SD units. Tests recover a planted
interaction, avoid a batch-composition artifact and retain unavailable arms at p=1.

Primary/target-omission same-sign flags were 0 TF, 25 heme (17 above/8 below), 8 membrane
(6 above/2 below), 6 granulocyte TF (1 above/5 below), 11 progenitor (8 above/3 below), 10 Cycling
(9 above/1 below); three low-coverage panels remained unscored. Heme target omission
does not change its score because the perturbation targets are not panel members.
CBL+CNN1 has heme contrast .877 control SD, pointwise interval .633–1.120 and Holm
p=4.61e-9; this recovers prior published context in the same data, not novel biology.
Its membrane contrast .191 is unflagged. Full JSONs are archived with lossless
compression, exact raw/compressed hashes, all 1,179 portable rows and summaries.

The next sensitivity was fixed after those primary results and before linear
outcomes (`7a808ed`): compare log and linear normalized expression using the
same markers/split and respective discovery-control standardization; joint Holm
covers all 2,358 scale-specific tests. This changes both transformation and fitted
standardization, so it cannot isolate a pure log-link effect. All pairs remain,
including the known-context pair regardless of its outcome. The sensitivity is
running; neither stable results nor independent experimental units are assumed.

The log/linear sensitivity then completed, reproducing every original panel/pair
primary and omission estimate exactly before the new scale was analyzed. Under
joint 2,358-test correction and same-sign/effect-floor rules on both scales and
their target omissions, 21 heme contrasts remained (14 above/7 below), 6 membrane
(4 above/2 below), 6 granulocyte TF (1 above/5 below), 4 progenitor (all above) and 6 Cycling
(5 above/1 below); TF had none and three panels remained unscored. The known
CBL+CNN1 heme contrast was .877 SD on the log score and 1.995 SD on the linear score;
both joint p-values are <1e-8. This is scale robustness in the same data, not
independent confirmation or evidence of a new mechanism. All 1,179 joint records,
linear outputs and exact raw/compressed hashes are archived. The dual suite passed
separately, bringing distinct exercised standalone suites to 18.

## 2026-10-18 — Capped NMF output and phase-diagnostic launch

The NMF fit completed its fixed 200 iterations in 1,799.97 seconds of concurrent
one-thread wall time, retaining the explicit convergence warning. Relative
uncentered Frobenius error is .33394, not a PCA explained-variance ratio. It matches
TF, heme, effector and Cycling panels with6,21, 54,0 upward associations; other panels
are unscored. All four matched identity tests have Holm p=.009. This unconverged
cap cannot support a definitive method ranking, and effector marker-presence vs
direct-score detection eligibility remains a mismatch. Actual cell-factor/model
hashes, dimensions, nonnegativity and source fingerprints verified; full outputs,
all 945 rows and warnings are archived. The queued control-only SAE then started
without interrupting the all-cell primary. Both full-size fits remain unfinished.

An explicitly exploratory archived-output comparison found 23 of 52 PCA Cycling-up
calls with direct Cycling-marker AUC<=.45; 48 of 52 had marker AUC<.5. This need not
mean PCA failure: broad scores may combine different cell-cycle programs. A phase
diagnostic was fixed before its outcomes (`241b95c`), using 97 external S/G2M genes,
the existing split and all 105 singles. Both directions plus all 945 original PCA
hypotheses enter joint Holm 1,365 families, including separate target omission.
The new family-retention test passed; full phase scores are running.

The maintained Scanpy tutorial moved, and its current notebook points to a
checksum-pinned Dropbox ZIP. That host was added to the network draft but returned
403. A Broad Institute workshop copy is pinned to its commit and exact Git blob;
ETH/NBIS indexed copies have the same blob. All 97 symbols are distinct and the
current tutorial defines first 43 as S, remaining 54 as G2M. The canonical ZIP itself
was not retrieved or independently verified. Source provenance records that limit.
Crossref verifies the cited Tirosh paper's 2016 print year despite the tutorial's
2015 wording; this is bibliographic verification, not full primary-text review.

The phase diagnostic completed. Forty-two of 43 S markers and all 54 G2M markers
were eligible (only MLF1IP was absent; no alias mapping). Across all 105 singles,
joint correction gave zero upward calls in either phase, 52 S decreases and 37
G2M decreases, all retained after target omission. Of 50 PCA Cycling-up calls
retained under joint correction, 28 also had S decreases, 21 G2M decreases, and 19
both. The original 52-case set gives 29, 22 and 20 respectively. The finer-phase
upregulation explanation is unsupported under these score/threshold choices;
neither actual phase arrest nor PCA failure is established. All directional
records, 105 paired comparisons and the complete 1,365-test families are archived.

## 2026-10-20: full-size SAE completion and common-inference continuation

All three primary seeds finished their declared 200 epochs during the pause. The
archive watcher validated and pushed complete JSON/log pairs and joint Holm27.
Today a second artifact check confirmed actual source/input/model/activation hashes,
checkpoint epochs and dimensions. TF, effector and Cycling matched in every seed
(raw identity p=.001, joint Holm p=.027); heme and all other panels remained unscored.
The sole repeated legacy candidate is effector-associated KLF1, with background
percentiles63/73/68. The TF list had11 candidates in seed2 only (percentile98).
No candidate combines the old background gate and at-least-two-seed agreement.
These gates would not by themselves establish specificity if met. The epoch20
heme match does not persist to the declared endpoints; no favorable epoch is selected.

The separate control-only seed also finished200 epochs on11855 controls, encoding
all111445 cells. Heme, effector and Cycling matched; KLF1 was the only heme candidate
(background percentile82), with no effector/Cycling candidates. All four fits share
scientific source hashes. Their complete36-hypothesis exploratory identity family
gives Holm p=.036 for matched tests. The real control-only fit recovers named
features despite the earlier synthetic0/8 result; one seed does not establish
general success or failure. Full compact outputs/log and validator are archived.

A separate common-inference plan was committed before new annotation outcomes
(30af531), after the original model/baseline/phase results. A pre-execution precision
correction (5d19cc9) fixes9999 identity permutations for all6 arms:999 would have
minimum p=.001 and could never reject first in a54-test Holm family. Saved SAE/NMF
scores will enter discovery-only annotation and identical inference controls;
PCA fit_transform must be exactly recovered and all945 old comparison estimates
reproduced before use. Joint BH covers5670 associations and Holm54 identity tests.
All panels, failures and105 singles are retained; no new SAE optimization occurs.
The new family/recovery suite and required V0b checks passed before execution.
Twenty distinct standalone suites have now been exercised,15 together and5 separately.
The full common-inference run is in progress; unit tests are not full-data evidence.

### 2026-10-20: completed common inference and signed-weight diagnostic

All six arms completed common inference. All 5,670 association hypotheses and 54
identity hypotheses were retained; 9,999 permutations executed per arm. Exact
original PCA/NMF annotation, AUC and p-value estimates and every control split
reproduced before combined interpretation. All matched identity tests have raw
p=.0001 and joint Holm p=.0054. All six models match Cycling, but only PCA has
upward associations (52); the other five have none. All-cell SAE TF counts are
0/0/5, effector 1/1/1; KLF1 remains the sole repeated candidate. Control-only heme
no longer matches (adjusted p=.0990). Complete outputs and all failures archived.

Twenty of the 52 PCA-up associations have both external phase AUCs<=.45; this is
a descriptive effect rule, distinct from19 jointly significant doubly decreased
cases among50 PCA calls in the phase diagnostic. A separately declared signed-
weight follow-up (a0402ff) passes its fixed score-fidelity gate: RMSE/SD=.000154
and maximum AUC difference=.0000378. All428505 contributions across105 singles
and4081 genes recover projected means (max absolute error2.22e-16). Only32 of97
external phase genes enter the PCA input. In all20 direction-discordant cases,
the represented phase genes lower the mean while other positive-weight and
negative-weight classes raise it. This is deterministic mean algebra, not rank
attribution or biological mechanism. Among all105,94 have positive negative-
weight sums, so the behavior is not unique to the selected cases. Complete
contributions, partition estimates, source/model/score hashes and checks retained.

The four-panel manuscript now includes the completed fits and these diagnostics,
plus a claim-to-evidence ledger. The new comparison figure is visually inspected
and deterministic as PNG/SVG. The attribution test passes separately, for21
distinct suites (15 together and6 separately); syntax, required V0b and all four
frozen specs/two bundles pass again. No new regulatory mechanism or SAE advantage
is established, and the manuscript is not submitted.

Original authors' Git blobs were retrieved and verified. Their guide-performance
code distinguishes the single Perturb-seq guide from two-guide averages in the
separate fitness screen, and explicitly studies cassette orientation. Their
examples include PTPN12/SNAI1 and CBL/CNN1 CD235a/CD71 cDNA plotting. This is code
context, not full original-article/supplement review or independent validation.
Ten citation records now include Tirosh2016; no original-paper gene-list audit
is claimed.

## 2026-10-20: completed numerical and historical-program follow-ups

The bounded NMF warm restart met its declared tolerance after380 further
iterations,2339.6seconds on four CPU threads. Frobenius error3077.269531 became
3050.078613 (0.88% lower); the stopping reference resets, so this is not an
uninterrupted580-iteration trajectory or a global-optimum claim. Joint seven-fit
BH6615/Holm63 retains the original capped fit. New TF/heme/membrane/effector/
Cycling counts are10/35/27/52/0. A membrane annotation appears after a small
error change. Actual factors/source hashes, complete raw comparisons and warnings
were checked and pushed automatically (ea054d5).

The exact historical K562 list was fixed before new context outcomes (54b1633):
HBZ,ALAS2,HBG1,HBG2,HBA1,SLC25A37,GYPB,BLVRB. All eight are in the unchanged
4081-gene input and eligible for direct scores. Original945 direct comparisons
and parameters reproduced. Complete eight-arm BH8400/Holm80 includes all ten
panels, all105 singles, seven learned fits and direct markers. Historical-program
counts are63/61/61 all-cell SAE,8 control-only SAE,65 PCA,54/52 NMF and57 direct.
Sixty candidates recur in all three primary seeds,62 in at least two. This is
computational agreement on the same cells. CBL,UBASH3B,PTPN12 recur, but PCA/NMF
also recover the program. The list changes context and granularity together, and
its historical selection on the same data is not independent confirmation.
Nominal p/q calculations do not replay that selection process.

The unchanged reference/depth helper reproduced the full primary exactly. All
4200 reference comparisons,1050 depth comparisons and8400 technical-partition
records are retained. Thirty-nine of57 direct calls survive all references/depth;
CBL,UBASH3B,PTPN12,KLF1 remain. Minimum reference AUCs are.679,.616,.683,.814 and
depth AUCs.750,.701,.753,.862. Technical references are not biological replicates.
The extraction adapter records that its added provenance describes lossless
extraction of the actual direct_fit, not a new primary execution.

A separately fixed gene diagnostic (1b4437b) retains10080 hypotheses, including
840 individual transcripts and840 single-marker omissions. Forty-six primary
calls survive every omission;38 also survive references/depth. All three historic
leads meet the criterion for all eight individual genes and every omission.
KLF1 has seven individual calls because GYPB AUC=.519, despite all omissions
passing. No absence or mechanism is inferred. This prompted a fixed denominator
sensitivity (ec6fa81), retaining all11025 old/new hypotheses after recovering945
original estimates. Excluding the eight signature genes from the library sum,
while fixing score weights, retains57 aggregate calls. KLF1 GYPB becomes.521
and remains below the.55 floor. This adjustment does not explain the threshold
result; neither denominator estimates absolute RNA abundance. No genes or
thresholds were revised after outcomes.

The target-RNA screen retains109 comparisons:90/105 direct labels pass,three
symbols are absent andtwo of four author alternatives pass. BAK1/BCL2L11 do show
RNA induction, so an unexplained old-code exclusion cannot mean absent RNA
induction. A descriptive pairing retains all105 labels:19/20 PCA-up/both-phase-
down cases have target-RNA increases; C3orf72 is the exception. No cases are
dropped and no new significance test is claimed.

Seven actual model/count archives (1696336750 compressed bytes) passed complete
streamed checksum roundtrips. Original files remain unchanged. A GitHub draft
release was created, but uploads.github.com is denied by CONNECT HTTP403 before
authentication. It has no assets. The host was saved in the environment draft;
saving does not apply runtime policy. Authorship/affiliation/contact/funding facts
were requested asynchronously and remain unconfirmed. No paper was submitted.
Twenty-four distinct standalone suites have passed (15 together,9 separately).


### 2026-10-20: publication review copies and final consolidation

A 25-page PDF and editable Word review draft were rendered with system
Pandoc/XeLaTeX. All five figures are embedded; source/figure/output hashes, PDF
content checks and Word ZIP CRC are recorded in docs/review/render_manifest.json.
The abstract starts on its own page, DOI links wrap without clipping, and the
supplementary figure/caption stay together. The legacy robustness image's title
was narrowed to displayed sensitivity settings; its count calculations were
unchanged, and two PNG/SVG renderings were byte-identical. Required V0b logic
was rerun successfully after the presentation edits. Author facts remain pending.

The PR description now covers completed primary/control fits, the positive
historical-program result and its selection limits, numerical NMF sensitivity,
signed PCA attribution, target RNA and all 24 exercised suites. The cloud startup
draft records all fits as complete and preserves the tested installation and
network settings. Saving the draft does not apply or publish it. No submission,
merge, new mechanism, independent replication or public archive deposition is
claimed. Large-asset upload remains denied at the network CONNECT stage.


### 2026-10-20: final candidate-membership and literature access checks

An exact post hoc set comparison confirmed that all 57 direct-score candidates
are identical under the original context family and both count denominators in
the expanded family. Candidate lists, CSV source hashes and comparison scope are
archived separately; no new hypothesis test or threshold was introduced.
All 150 local Markdown links in the consolidated review/evidence commit resolve.

Additional original-Norman full-text routes failed. A supplementary API returned
HTTP200 carrying an explicit not-open-access XML error, not supplements. Crossref
and EuropePMC link the published article to public preprint 10.1101/601096, but
its PDF was denied by proxy CONNECT before retrieval. Metadata/abstract review
is distinct from full-text review; the latter remains incomplete. No access
restriction was bypassed and no known regulator was relabeled as novel.
