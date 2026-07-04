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
knob-tuning artifact. Artifacts: `data_g0{,_human}/robustness/`. Next: R2 (retrain
at L0 spanning the 20-50 band) and R3 (20 seeds).
