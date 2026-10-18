# PROJECT MASTER SUMMARY

> **Continuation audit, 2026-10-18:** the record below predates additional same-day
> checks. Read [PROGRAM_MODULATION_AUDIT.md](PROGRAM_MODULATION_AUDIT.md) and
> [MANUSCRIPT_DRAFT.md](MANUSCRIPT_DRAFT.md) for the continuation. The old two
> quantitative modulation gates pass a generic-response counterexample in 20/20
> simulations. A no-marker match silently selected feature -1; this is now guarded.
> A separate controls-fixed identity-heterogeneity test is implemented and calibrated
> on synthetic nulls. No new real-data biological discovery has been established.
> The causal PCA/NMF ladder trained an SAE on each reduced input; it did not compare
> against direct PCA/NMF features. Historical original causal JSONs remain missing.
> Current entrypoints and strategy docs have been corrected, and the false 10-fold
> enrichment panel has been replaced. Preserve the frozen specs and bundles.

The complete record of the SAE / scRNA-seq hematopoiesis project, start to finish, as of
**2026-10-18**. Compiled by reading every tracked file (the 21 reference docs in `docs/`,
`README.md`, `CLAUDE.md`, `LAB_NOTEBOOK.md`, 22 scripts, 9 test suites, the frozen specs and result
bundles), all 90 commits, and the surviving session transcripts, and by re-running the full
verification battery in a fresh container. This document indexes and
narrates the authoritative docs; it does not replace them. Where another doc is stale or wrong,
section 17 says so explicitly.

How to read it: section 0 is the state in one screen. Sections 1-9 are the narrative in order.
Sections 10-16 are registries (scorecard, every self-correction, frozen hashes, every file, every
dataset, compute, competitors). Section 17 lists stale or wrong statements elsewhere. Section 18 is
future directions, ordered around the STS deadline. Section 19 is the resume checklist. Section 20
is a glossary.

---

## 0. The state of the project in one screen

**Deadline first.** Regeneron STS 2027 applications close **November 5, 2026, 8:00 PM ET**
(support-request cutoff November 4), per Society for Science. Today is 2026-10-18: **18 days**.
Scholars (top 300) are announced January 7, 2027; finalists (top 40) January 21; finals March
11-17. ISEF 2027 goes through an affiliated regional fair first. Section 18 is ordered around this.

**What the project is.** A fully computational high-school project (Jacob Barzideh) that trains
sparse autoencoders (SAEs) on single-cell RNA-seq of blood development and asks whether the "gene
programs" they find are real, and whether they are causal.

**What it set out to show, and what happened.**
- The original headline, "asymmetric modularity" (granulocyte commitment is one unified program,
  erythroid commitment is distributed): **dead**. Not supported in mouse (0/10 seeds) or human
  (0/10). The SAE points the opposite way, and even that reversal is method-dependent: across five
  decomposition families and two species the original claim survives in 1 of 10 cells.
- That rigor work (Gate 0) is complete, pre-registered, frozen, and stress-tested. It is the floor.
- The causal program (Component 2) asked whether SAE programs respond causally to CRISPR
  perturbations. Four executed tests, **no positive result at single-perturbation resolution**:
  1. Replogle CRISPRi (K562, RPE1): grounding is real and beats a strict null (K562 22/162 =
     13.6%), but only for housekeeping programs; no lineage TFs; SAE no better than PCA/NMF.
  2. Norman CRISPRa rate test: 0.147 vs shuffle null 0.175 (p = 0.902). Which TF was overexpressed
     adds nothing beyond a generic perturbation response.
  3. Regulon recovery v1: uninformative (the test never fired).
  4. Regulon recovery v2: a real negative. SAE features match curated TF regulons no better than
     the best of 2048 random gene sets.
- **One positive asset survives.** A seed-stable, hand-checkable erythroid/hemoglobin program in
  K562 (HBZ, ALAS2, HBG1/2, HBA1, SLC25A37, GYPB, BLVRB), activated under CRISPRa by CBL, UBASH3B
  and PTPN12, negative regulators of KIT/EPOR signaling. Not yet tested with a dedicated,
  non-circular null.

**The experiment that decides whether the project has a positive finding** is built,
oracle-tested, pre-registered (spec v1.1 `ddb29bcd`), smoke-tested on the real data today, and
**not yet run**: `src/program_modulation.py`. Minutes per seed on a rented GPU; 2 hours per seed on
this CPU container. The verified recipe is in section 18.1.

**Found today while compiling this summary** (all fixed or recorded; section 11, items 36-45):
- Every Component 2 run fit the SAE on ALL cells, although the frozen spec says control cells
  only. Only the all-cell mode passes the pipeline's own synthetic self-test (8/8 planted
  regulators); control-only gets 0/8. The spec now carries a dated erratum; the code exposes
  `--train-on {all,control}` with the default unchanged, so every past result is reproducible.
- My 2026-08-20 re-download instruction (`pt.dt.norman_2019()`) pointed at pertpy's
  log-normalized file. The correct call is `pt.dt.norman_2019_raw()`. The August results used the
  raw file (verified against both files' schemas). `prep_norman.py` now refuses or repairs
  non-raw input.
- GATA1 is **not in the Norman library** (0 cells). An August note implied it was tested and
  failed to ground; it was never perturbed.
- The regulon-recovery chance comparison was worded too strongly ("at or below chance for every
  regulon size"). Exact: 8 of 42 TF x seed cells sit above the chance mean, 0 above its 95th
  percentile. The conclusion is unchanged.

**Repo state.** Branch `claude/code-handoff-doc-2v3jt3`; draft PR #3 open (mergeable). Gate 0 was
merged to `main` via PR #1 on 2026-07-04. All 9 test suites pass in a fresh container and every
frozen hash verifies (2026-10-18).

---

## 1. Identity, people, constraints

| | |
|---|---|
| Researcher | Jacob Barzideh (GitHub `jbar1125`), high school |
| Program | School science-research program. `PROJECT_HANDOFF.md` names it "Jericho Science Research"; a later session summary wrote "Rocho", almost certainly a slip |
| Mentor | Mrs. Han (faculty mentor); remote affiliation with Cold Spring Harbor Laboratory (mentored, not employed) |
| Competitions | Regeneron STS 2027 (deadline 2026-11-05) and ISEF 2027 |
| Project type | Fully computational; public, de-identified data; no wet bench (mentor confirmed the student cannot do independent bench work) |
| Related work | Co-author on a team ICML 2026 mechanistic-interpretability workshop paper ("Semantic Drift Recovery (SDR) Pipeline"; wrote part of the intro and results; decision was due around 2026-06-12; outcome not recorded in this repo). Splunk Agentic Ops hackathon with Aaron Yang |
| Lab notebook | JSR-compliant bound carbon copy, plus the electronic `LAB_NOTEBOOK.md` |
| Repository | `github.com/jbar1125/SAE-scRNA-seq-Research`, **public** since 2026-07-18 (made public so no access token ever touches rented compute) |
| Outside the repo | ISEF paperwork drafts (`isef_paperwork_synopsis.docx`, `isef_all_forms.docx`, `isef_addendums.docx`) and Phase 2 plan v1-v5 (`.docx`), all in Google Drive. A weekly summer reflection form for the research program was filled in during an earlier session and then shortened on request; its text was never committed and is not recoverable from the repo or the surviving transcript |

**Working style** (CLAUDE.md section 6, handoff section 15): extremely concise; no em-dashes;
per-bullet formatting; code ready to run; review of literature in formal third person;
methodology in imperative future tense; next steps phrased "I must ..."; brutal, constructive
honesty ("BE AS HONEST AS POSSIBLE, CONSTRUCTIVE CRITICISM TO ALL HEAVENS").

**Non-negotiable rules** (CLAUDE.md sections 1-3):
- Run code before presenting it: `py_compile`, `tests/test_v0b_logic.py`, a real or synthetic run
  of the changed path, then a verification table with every row PASS. Never hand over a command
  that has not been run; if it cannot be run, say so and label it UNVERIFIED. (Origin: the
  2026-06-30 `torch.manual_seed` NameError, shipped without running.)
- Pastes into rented or Colab machines: pin the commit SHA, `rm -f` then `wget -qO`, and a grep
  guard that prints OK or STALE before running.
- Honesty: never sugarcoat; lead with the bad news; distinguish a result from an artifact (a null
  from a blind or underpowered method is not a clean negative); label provisional numbers every
  time; never SHA-freeze output from a broken or inconclusive run; self-correct out loud.
- Science conventions: per-seed winner-take-all, never cross-seed feature-index matching (SAEs do
  not preserve feature identity across seeds); BH-FDR per test family; Cohen's d with p-values; a
  finding needs >= 3/5 mouse seeds or >= 2/3 human seeds (a 60% majority at 10 seeds).
- Git: develop only on `claude/code-handoff-doc-2v3jt3`; open PRs as drafts; never put
  credentials on rented machines.

---

## 2. The arc in six acts

| Act | Dates | Question | Outcome |
|---|---|---|---|
| 0. Phase 1 (pre-repo, Colab) | spring 2026 to 06-29 | Do SAEs on scRNA-seq show "asymmetric modularity"? | Claim formed on buggy pipelines; never validated |
| 1. Gate 0 | 06-29 to 07-02 | Does the claim survive a pre-registered, control-referenced test? | No (0/10 mouse, 0/10 human); modularity is method-dependent; frozen |
| 2. Hardening | 07-04 to 07-05 | Is that negative robust? | Yes, across thresholds, sparsity, seeds, gene counts, 5 method families, TRRUST |
| 3. Causal pivot + Replogle | 07-05 to 07-20 | Are SAE programs causally grounded on CRISPRi? | Real but housekeeping-only; SAE no better than PCA/NMF; wrong substrate |
| 4. Overexpression + Norman | 07-21 to 08-20 | Do overexpressed lineage TFs ground their programs? | Rate at null; regulon alignment at chance; one erythroid program survives; decisive test built |
| 5. Audit + consolidation | 10-18 | What is actually true, and what is next? | This document; four latent defects found, fixed or recorded |

---

## 3. Master timeline (dated, with commits)

**Before the repository (spring 2026, Google Colab, free T4).** Phase 1 notebooks
(`SAE_rebuilt.ipynb` and variants, in Drive): mouse SAEs (5 seeds, 128 latents), human SAEs (3
seeds, 25k CELLxGENE cells), NMF baseline, V0a pseudotime fix, ISEF paperwork, Phase 2 plan v1-v5,
the ICML workshop paper. 2026-05-29: the canonically named Drive notebook
`v0b_module_definitions.ipynb` was last modified; it is actually version 1 (the "notebook trap").

| Date | Commits | What happened |
|---|---|---|
| 06-29 | `4ac09e9` | Phase 1 handoff committed (551-line snapshot) |
| 06-30 | `4fdb988` `b320bb9` `ca1a746` `ed79c6a` `b1aaa5e` `70315d4` `874f244` | Initial commit; canonical V0b v2 script, audit, logic tests; the shipped-untested `torch.manual_seed` NameError fixed; CLAUDE.md rules written; first V0b v2 run (NOT SUPPORTED, but a blind test); V0b v3 continuous loading; v3's SUPPORTED shown to be an artifact. PR #1 opened |
| 07-01 | `1f39d29` `84f2b22` `f5fd742` `f4ca4c1` `a43dc32` `bae8c7a` `edd9502` `312ea5b` `91a023c` `f3bc7c8` `ea259ff` `35173e4` `15d00e4` `1a0866a` `41e4b81` | Overcomplete SAEs + marker-aware preprocessing; Paul15 log-transform fix; v3.1 metric; marker symbols corrected; the reversal; strategy doc; novelty claim retracted; plan v6; pre-registration frozen (`fc342829`); PCA/NMF baselines; leave-one-marker-out; seed auto-detection; h5 fallback loader |
| 07-02 | `b36fc0b` `2181170` `f562c66` `1b8495a` | Gate 0 executed and frozen, mouse then human. Gate 0 complete |
| 07-04 | `e3440d7` `f3d21eb` `c949ed2` `e30ab93` `9a05f53` `5288cfb` `e4b2af6` `485792a` `31587de` `92fdc3c` `5dc751c` `739dc4f` `068f1f9` `4a8da2a` `a04d1ec` | Repo reorganized; LAB_NOTEBOOK, README, VERSION_HISTORY; R1; COMPUTE.md; **PR #1 merged to main**; PR #2 (main sync); R2; TRRUST; SCENIC deferral; COMPONENT2_PLAN; `--seed-start`; R3; R4; GRN family. PR #3 opened |
| 07-05 | `21b0e49` `4978a56` `7c0e5d1` `ff4a101` `6a0d328` `5f5f1be` `f323bca` `1a48eeb` `a1483f8` `a73c0a5` `396ee6f` `7bfd18e` `de18ac2` `cf5ee09` `1cedbb0` `bfb2fd9` `762bb2b` | Human robustness; ICA (non-convergence artifact caught); figures; abstract; reproducibility pass (hash bug fixed); visual brief; ELEVATION_PLAN pivot; causal metric; pipeline; RunPod guide + causal spec v1 (`b7d28fae`); causally-supervised SAE failed and was shelved; label-shuffle null |
| 07-17 | `c00e54d` | `--n-hvg` to bound memory |
| 07-18 | `d91f5c2` `254a9f7` `5f6643c` `3eef866` `d580218` `7fb9765` `6474e7c` `757d6fe` `b70b248` `de632c1` `c701f7c` `ed1b7ab` `0c3e65d` `5bd9a0a` | GPU day on vast.ai. Broken first run (0/1789) and two metric bugs fixed (spec `b766d21c`); backlog; first real number (13.6%); 5 seeds + sweep; essential-gene tooling; column specificity (v4 `029ab066`); next-directions map; falsifiable predictions; DIFFERENTIATION; PCA/NMF rungs; diagnostic figure; v4 failed and rolled back (v5 `98104b2f`); IMPACT_STRATEGY; effect-size control M1 (v6 `81b44e04`) |
| 07-19 | `5d2ff84` | M1 result: effect size is not the confound |
| 07-20 | `f0a28da` `617def2` | Representation ladder flat; RPE1 replication |
| 07-21 | `cf07251` `627b50e` | Overexpression direction (v7 `ee5dacb7`); h5ad schema inspector |
| late 07 to early 08 | (no commits) | Joung TF Atlas judged too heavy (16.7 GB raw TAR); switch to Norman 2019 CRISPRa via pertpy; three data-prep failures on vast.ai |
| 08-04 | `a356749` | Guard against 0 testable perturbations (ZeroDivisionError) |
| 08-10 | `4b335e6` `566c334` `9ab480c` `39c5533` `08b5e50` `2f5578a` `743cccc` `8bc8734` `853758b` | Top-program-genes readout; provisional positive logged; `prep_norman.py`; null run negative, positive retracted; regulon recovery v1 (`bbfac449`); regulon diagnostic; v2 (`c653f6d1`); v2 verdict, line closed |
| 08-20 | `4a3728e` | Program-modulation test (`b396834d`) + shuffle-null self-correction. Last exchange before today: "where do I run this" |
| 10-18 | (this session) | This audit: train-cell erratum, raw-data guard, GATA1 correction, chance-wording correction, master summary |

90 commits are reachable from the branch head; 58 of them are on PR #3 and not yet in `main`.

---

## 4. Act 0: Phase 1 (before the repository)

### 4.1 Data and model
- **Mouse:** Paul et al. 2015 myeloid progenitors via `scanpy.datasets.paul15()`: 2,730 cells x
  2,000 HVGs, 19 clusters (1Ery-6Ery, 7MEP, 8Mk, 9GMP, 10GMP, 11DC, 12Baso, 13Baso, 14Mo, 15Mo,
  16Neu, 17Neu, 18Eos, 19Lymph). Matrix MD5 `60183a17983c8b977d036e0f3a58da61` (superseded).
- **Human:** a 25,000-cell stratified subsample of 759,000 healthy bone-marrow cells from the
  CELLxGENE Census (the Census release date was never pinned).
- **SAE:** `Linear(2000->128)`, ReLU, `Linear(128->2000)`; AdamW lr 5e-4; L1 0.1; 300 epochs;
  unit-norm decoder columns; dead-latent resampling every 500 steps. Mouse 5 seeds, human 3. Mouse
  reconstruction MSE 0.6157 (sd 0.0055); L0 about 45-48 active features per cell.
- **NMF baseline:** sklearn, NNDSVD init, rank 128, 5 seeds.

### 4.2 V0a: the first real bug (pseudotime root)
The original diffusion pseudotime used `iroot=0`, cell W31105, a committed 7MEP: every pseudotime
value was rooted inside the erythroid lineage. Replaced with Palantir rooted at the 10GMP centroid
(W39126, index 2690) with 7 manual terminal states. Spearman vs canonical cluster order 0.733;
1Ery self-probability 0.90. Output `cell_metadata_palantir.csv`. (The frozen Gate-0 metric does not
use pseudotime; it is decoder-loading based.)

### 4.3 The Phase 1 claims (all provisional; none survived)
- Asymmetric modularity: granulocyte unified module (mean z -5.53 mouse, -10.68 human);
  erythroid distributed with sub-additive double knockouts; "replicated" in 3/3 human seeds.
- Commitment crossover at pseudotime ~0.37 (derived quantity); 49 "bimodal" features (sigma=3
  peak-finding heuristic, not a statistic).
- Participation ratio SAE 24.2 vs PCA 43.4 or 42.3 (documents disagreed); NMF "1.8x SAE" vs
  "34.94" (34.94/24.2 = 1.44x; internally inconsistent).
- NMF recovers the sub-additive erythroid hierarchy in 1/5 seeds vs SAE 4/5; Feature 73 (Id1)
  r = -0.828 with pseudotime; unsupervised recovery of GATA1, TAL1, FECH, Aqp1, Erfe, GFI1B, Id1.
- Feature-to-TF annotation table: 640 entries, SHA
  `bd5310926e2b0d735856cb0f4317ad8f401917267fdae71dcd254f6b4bd06373`, with a known bug (GATA1/TAL1
  absorbed into GFI1B by marker overlap).

### 4.4 The V0-V15 validation plan, and what became of each item

| Item | Plan | Fate |
|---|---|---|
| V0a | Palantir root fix | Done (above) |
| V0b | Gene-content module definitions | Superseded by the repo scripts and Gate 0 |
| V0c | Hartigan dip test for bimodality | Never built; the 49-feature count is retired |
| V0d, V14, V15 | Non-circular in-silico perturbation; OOD threshold bug (self-distance = 0); decoder/classifier scale bug | Never fixed. In-silico perturbation was abandoned for real CRISPR data (Act 3) |
| V1-V3 | Data, preprocessing, training reproducibility | Superseded by MD5-verified Gate-0 pipelines |
| V4 | Marker recovery across top-K | Superseded (v3.1 uses the full decoder, no top-K) |
| V5 | Cross-seed feature stability | Not done (per-seed winner-take-all instead) |
| V7, V8 | Root and pseudotime-method sensitivity | Partial; not continued |
| V9 | Logistic trajectory fits | Never done; the 0.37 crossover is retired |
| V10 | Modularity robustness | Done inside Gate 0 |
| V11, V12 | Participation ratio; NMF modularity | Done in Gate 0 at 512 latents |
| V13 | Cross-species concordance | Done in Gate 0 on Setty 2019 (not the CELLxGENE subsample) |

### 4.5 Phase 2 plan (v5 `.docx`, then v6 in the repo)
Five parts, 12 weeks, 4 decision gates: Part 1 van Galen 2019 AML (GSE116256); Part 2 Replogle
CRISPRi vs the 6.2% null (24 canonical TFs; Gate 2 strong pass >= 8/24); Part 3 Hegde 2025 tumor
marrow (GSE270148); Part 4 scGPT/Geneformer embedding SAEs; Part 5 Park 2020 thymus (optional).
Only Part 2's question was executed, in redesigned form (Act 3). Parts 1, 3, 4 and 5 were never run.

### 4.6 ISEF/STS paperwork (from the handoff; status of signatures unknown here)
Needed: Forms 1, 1A, 1B, 3 (risk assessment, computer use), Official Abstract; Form 2 defensively
if the mentor signs. Not needed: 1C, 4, 5A, 5B, 6A, 6B, 7 (7 pending fair confirmation). Rules
Wizard: "None of the Above". The AI Usage Policy addendum applies (prompt logbook + citation). Human
participants, vertebrate animals, hazardous agents: N/A (public de-identified data; the mouse data
were generated under Weizmann and Mount Sinai IACUC approvals).

---

## 5. Act 1: Gate 0 (Component 0), executed and frozen

### 5.1 The metric lineage (why there are four versions)

| Version | File | Method | Why it failed |
|---|---|---|---|
| v1 | Drive notebook (never in repo) | cross-seed index matching (`CROSS_SEED_THRESHOLD=3`); `H_ery/H_gran` entropy ratio; markers absent from the HVG set; `prob_19Lymph` omitted | The version that "supported" the claim. Cross-seed matching is invalid for SAEs; the ratio is biased by 3 vs 2 submodules; ~273 high-Lymph cells mislabeled (uncommitted 973 vs ~700) |
| v2 | `src/v0b_module_definitions.py` | hypergeometric enrichment of top-30 decoder genes, BH-FDR, per-seed winner-take-all (>= 2 marker overlap), normalized entropy + largest-submodule fraction | First run 2026-06-30: NOT SUPPORTED, but structurally blind. Diluted markers never clear q < 0.05 across ~896 mostly-null tests per seed; only Gran_Primary (29 significant) and Ery_Heme (2) passed. Small sets have an FDR floor (Gran_TF best-case p ~9e-4). Assignments SHA `f18b088a95e8e4f3248561d72a711f2bd6d1302127afd60b753a4e9d211d96b5`, NOT frozen |
| v3 | `src/v0b_v3_loading.py` | continuous decoder-mass enrichment, no gate, gene-label permutation null | SUPPORTED 5/5, but an artifact: the largest-fraction floor is 1/n (about 0.17 of the 0.244 gap is bin count); the count-fair entropy gap is 0.048; the null was mis-calibrated (Cycling control p = 0.001); Ery_TF and Ery_Membrane sat below the Progenitor control |
| **v3.1** | `src/v0b_v3_1_decision.py` | enrichment = (sum of abs decoder weight on markers / total) x (N_genes / K); strength = mean of the top 3 features; Progenitor-only baseline (Cycling secondary); abundance-matched permutation null (20 bins, 1,000 permutations, BH-FDR); REAL iff excess > 0.10 AND q < 0.05; SUPPORTED iff in >= 60% of seeds ery_n_real >= 2, ery_n_real >= gran_n_real, and ery excess-entropy > gran | THE metric. Pre-registered |

### 5.2 The marker-curation bug and the reversal (my error, owned in the record)
Earlier runs used globin symbols absent from Paul15 (`Hba-a1`, `Hbb-bs`) and omitted `Alas2` and
`Ermap`, so the erythroid effector program was invisible for a bookkeeping reason. Corrected to
`Hba-a2`, `Hbb-b1`, `Alas2`, `Ermap`: hemoglobin became the strongest feature in the data
(strength 5.07) and the direction flipped (erythroid 1 dominant program, granulocyte 3). Also
found: "Paul15 is pre-log-transformed" is false for `sc.datasets.paul15()` (raw maximum 168);
`preprocess_paul15.py --log1p auto` handles it (audit A11). Unfrozen intermediate runs: the
marker-aware 512-latent run at L1 0.1 had L0 ~213 (not sparse); the corrected-marker run at L1 0.5
had L0 ~63, just above the band.

### 5.3 The pre-registration (frozen 2026-07-01, commit `35173e4`)
`config/preregistration_spec.json`, SHA `fc342829191da5a0bb5244975c3f61c8934943cb512818446e5b5f0984f6011a`:
the v3.1 constants (EXCESS_FLOOR 0.10, PERM_Q 0.05, N_PERM 1000, N_ABUND_BINS 20, FEAT_TOP 3,
SEED 42), the decision rule, nine marker sets, cross-seed rules, MD5 at load, and the Phase 2 causal
null (6.2%, one-sided binomial, Gate 2 thresholds). Deliberately not frozen: the L1 strength and
dictionary size (tuned to the L0 20-50 quality band) and the seed count.

Mouse marker sets (frozen): Ery_TF (8) Gata1 Klf1 Tal1 Lmo2 Zfpm1 Stat5a Bcl11a Myb; Ery_Heme (7)
Alas2 Fech Hmbs Ppox Cpox Urod Alad; Ery_Membrane (7) Gypc Ank1 Rhag Aqp1 Epor Tspo Ermap;
Ery_Effector (2) Hba-a2 Hbb-b1; Gran_TF (3) Cebpa Cebpe Runx1; Gran_Primary (4) Mpo Elane Prtn3
Ctsg; Gran_Secondary (8) Ltf Lcn2 Camp Ngp S100a8 S100a9 Mmp8 Mmp9; Progenitor (4, the control) Kit
Cd34 Mllt3 Eif4ebp1; Cycling (16, secondary control). Human (`config/human_markers.json`) uses 1:1
orthologs: Ery_Effector is HBB HBA1 HBA2 (3 genes); Gran_Secondary has 7 (no NGP ortholog).

### 5.4 The executed runs (2026-07-02, in-container, CPU)

| | Mouse (primary) | Human (replication) |
|---|---|---|
| Data | Paul15 via a git-LFS mirror (direct download proxy-blocked), sha256 `bd75eb90...` | Setty 2019 CD34+ marrow, `marrow_sample_scseq_counts.h5ad` from the Palantir repo, sha256 `be4f8603...` |
| Matrix | 2,730 x 2,012; log1p + marker-aware HVG; z-scaled; MD5 `184ea5dd...` | 4,142 x 2,032; normalize_total(1e4) + log1p + marker-aware HVG; MD5 `62796189...` |
| L1 sweep at 300 epochs | 0.5 -> L0 67.4; 0.75 -> 38.9; **0.8 -> 37.0**; 1.0 -> 28.0; 1.5 -> 19.4; 2.0 -> 17.4 | 0.4 -> 99.5; 0.6 -> 54.4; **0.8 -> 32.5**; 1.0 -> 21.4 |
| Final SAE | 512 latents, L1 0.8, 10 seeds, L0 35.75 (34.6-37.0), MSE 0.4642 | 512 latents, L1 0.8, 10 seeds, L0 32.5, MSE ~0.666 |
| Verdict | NOT SUPPORTED, 0/10 seeds (6 needed) | NOT SUPPORTED, 0/10 |
| Counts | gran_n_real = 3 in 10/10 seeds; ery_n_real 1 or 2 (median 1.5) | gran median 2, ery median 1 |
| Control | Progenitor perm q 0.085 (correctly null) | Progenitor perm q 0.748 |

Sparsity lesson: L0 depends on training length, not only L1 (L1 3.0 gave L0 34.6 at 150 epochs but
about 17 at 300).

Mouse submodule strength (1.0 = null) and abundance-matched q: Ery_Effector 5.07 (0.0018, the one
dominant erythroid program); Gran_Primary 4.30 (0.0018); Gran_Secondary 3.36 (0.0018); Gran_TF 2.93
(0.0154); Ery_Heme 2.40 (0.0018, borderline vs the floor); Progenitor 2.33 (0.085, control);
Cycling 2.01 (0.0018, a real program as expected); Ery_Membrane 1.99 and Ery_TF 1.90 (below the
control).

Human: Gran_Primary 5.44 (0.0045); Gran_Secondary 4.02 (0.043); Ery_Effector 3.75 (0.216, present
but not significant: CD34+ selection under-samples hemoglobin-high cells); Cycling 3.67; Ery_TF
2.74 (0.018); Gran_TF 2.74 (0.298); Ery_Membrane 2.51; Progenitor 2.22 (0.748); Ery_Heme 2.21.

### 5.5 The real finding: method-dependence
Participation ratio: mouse SAE 195.9, PCA 14.0, NMF 5.2; human SAE 278.6, PCA 270.7, NMF 5.9.
Same metric, same matrix, different decompositions (median ery/gran real programs):

| Method | Mouse | Human | Supports the original claim? |
|---|:-:|:-:|---|
| PCA | 1 / 1 | 1 / 3 | no |
| ICA (rank 50) | 1 / 3 | 2 / 3 | no |
| NMF | 2.5 / 2 | 0 / 0 | mouse **yes** (the lone outlier); human detects nothing |
| SAE | 1.5 / 3 | 1 / 2 | no |
| GRN (co-expression) | 2 / 2 | 1 / 2 | no |

The original claim survives in **1 of 10** method x species cells (mouse NMF, the most
parts-based, artifact-prone method). The SAE's granulocyte-more-distributed direction is echoed by
ICA in both species and by human PCA/GRN, but NMF disagrees in both species, so even the reversal
is method-contingent. No unsupervised decomposition yields a method-invariant program structure.

### 5.6 Marker robustness
Leave-one-marker-out: **0/34** mouse and **0/36** human drops flip the verdict. In mouse, dropping
`Alas2` pulls Ery_Heme below the floor, and dropping either globin retires the 2-gene Ery_Effector,
but erythroid never reaches the claim.

### 5.7 Frozen decision bundles
Mouse `172861417a923a83705812a96e5ad1e555055e968323289879447e53beb84396`; human
`cc8159c882cf326bca5bbbc7f74d110d4a81338ec0ba07e208eb935bc8f729d3` (decision JSON `28301fec...`,
per-seed CSV `5fdc2e8c...`). Bundle = `cat decision.json per_seed.csv | sha256sum`. Both verified
2026-10-18.

### 5.8 Limitations (as recorded)
Structural (decoder-mass concentration), not causal; one dataset per species; Ery_Effector is two
genes in mouse; the human globin program is not significant because of CD34+ selection; ICA is at
rank 50, not 512; the 20-50 L0 band and the 0.10 floor are project conventions (frozen before the
runs). Minor text inconsistency: the mouse decision JSON's `plain_finding` string says "erythroid
(2)" while its own `median_ery_n_real` is 1.5 (rounding in the text only).

---

## 6. Act 2: the hardening battery (2026-07-04 to 07-05)

Every unit ran on CPU in-container; artifacts are the JSON summaries under `data_g0*/robustness/`,
`trrust/`, `grn/`, `ica/`. Bottom line: no stress test resurrects the original claim.

| Unit | What varies | Mouse | Human |
|---|---|---|---|
| R1 decision knobs | excess floor {0.05, 0.10, 0.15, 0.20} x q {0.01, 0.05, 0.10} = 12 configs, 500 permutations | 0/12 support; gran >= ery everywhere | 0/12; same |
| R2 sparsity band | L1 -> mean L0 | L0 27.7 (gran 3 / ery 2), 42.2 (3/2), 52.0 (3/1, just out of band) | L0 21.4, 32.5, 41.5 (gran 2 / ery 1 throughout) |
| R3 seed count | 10 -> 20 seeds (`--seed-start`) | gran 3 in 20/20; ery 1 in 10, 2 in 10 | gran 2 in 20/20; ery 1 in 19, 0 in 1 |
| R4 gene count | 1000 / 2000 / 3000 HVGs | 1023 / 2012 / 3006 genes; L0 90.1 / 35.75 / 18.1; gran 3 throughout | 1036 / 2032 / 3026; L0 97.6 / 32.5 / 12.9; at 1036 genes both lineages 1 (a tie) |
| TRRUST cross-ref | curated TF -> target edges | null (0-10% of seeds; human DB on uppercased mouse symbols, underpowered, not evidence) | erythroid TF targets (31 in panel) overlap the erythroid-TF feature's top-30 by 4-7 genes; significant in 10/10 seeds, median p 6.4e-5; granulocyte 50%, median p 0.13 |
| GRN (4th family) | abs-Pearson regulons, NOT motif-pruned SCENIC | 176 TFs; ery 2 / gran 2 (a tie) | 206 TFs; ery 1 / gran 2 |
| ICA (5th family) | FastICA | rank 512 does not converge (0/3); rank 50 converges 3/3: ery 1 / gran 3 | rank 512 converged 1/3 and falsely read "supported" (ery 3 / gran 2); rank 50: ery 2 / gran 3 |

Deferred or blocked, documented rather than faked: full SCENIC (cisTarget motif databases
proxy-blocked); a second, non-CD34 human marrow dataset (no reachable source); mouse TRRUST (all
mirrors 404). Also produced: deterministic figures (`fig_method_dependence`, `fig_robustness`;
`svg.hashsalt` plus no embedded date, so reruns are byte-identical), `docs/ABSTRACT.md`,
`docs/web/result_brief.html` ("Gene-program modularity in blood-cell development is an artifact of
the method, not a property of the data"), and `docs/REPRODUCIBILITY.md`. The reproducibility pass
caught that the pre-registration file had a trailing newline (raw-byte hash `6398c5ed...` instead
of the recorded `fc342829...`); the bytes were normalized and the content verified identical.

---

## 7. Act 3a: the causal pivot (design, 2026-07-05)

**Trigger.** The user's verdict on Gate 0: "It is not a strong result as of now and will not win
anything ... make this a million times better." Correct: a negative plus "methods disagree" is
hygiene, not a discovery (`docs/ELEVATION_PLAN.md` section 1).

**The question.** Kendiukhov 2026 (arXiv 2603.02952) trained SAEs on Geneformer and scGPT
residual streams and found only 3 of 48 TFs (6.2%; 10.4% with a multi-tissue control) show
target-specific responses to Replogle CRISPRi. Is that ceiling a limit of SAEs or of
foundation-model embeddings? Hypothesis: SAEs trained directly on expression ground more.

**The triviality trap and its defenses** (ELEVATION_PLAN 3b). "Expression sees the genes, so of
course it wins." Planned defenses: (1) train the SAE on control cells only; (2) match by program
with the perturbed gene excluded from both the program and the DE signature; (3) demand
specificity. **Erratum (2026-10-18): defense 1 was never implemented** (section 9.10). Defenses 2
and 3 and cross-fitting were.

**The metric as finally built** (`src/causal_grounding.py`, spec v7 + erratum):
1. The SAE's dictionary: `train_topk_sae` in `src/causal_pipeline.py`: TopK (`k` largest
   activations kept, ReLU), `Linear(d->L)` encoder, bias-free `Linear(L->d)` decoder with
   unit-norm atoms, AdamW lr 4e-4, 200 epochs, batch 1024, mean-centered input. Fit on all cells
   in every real run (labels unseen).
2. Program of feature f = Pearson correlation of f's activation with each gene over **control**
   cells. Matchable only if f fires in >= 50 control cells.
3. Cross-fitting: each perturbation's cells are split A|B (seeded). Match on A, test on held-out B.
4. Match: argmax over features of the dot product of the normalized program with the normalized
   -DE (knockdown) or +DE (overexpression) signature, perturbed gene excluded. Match confidence:
   the best alignment must be a robust-z (median/MAD) upper-tail outlier among all features,
   p < 0.10.
5. Test: one-sided Mann-Whitney of the matched feature's activation in split-B cells vs all
   control cells; AUC = U / (n_pert x n_ctrl). Floor AUC <= 0.45 (down) or >= 0.55 (up). BH-FDR
   q < 0.05 over perturbations.
6. Column specificity (optional; default OFF since v5): the perturbation must move the feature more
   than other perturbations do.
7. Effect-size control M1 (default ON since v6, alpha 0.10): the matched-feature movement must be
   an upper-tail outlier among the k nearest perturbations in total transcriptional effect size;
   a no-op below 30 perturbations.
8. Pipeline-level label-shuffle null (`shuffle_control`): permute labels among perturbed cells
   with controls fixed, n relabelings, empirical p = (1 + #null >= real) / (1 + n).

`load_perturbseq` applies normalize_total(1e4) + log1p, keeps the top-N HVGs UNION every tested
perturbation's gene (plus any `keep_genes`), and densifies after subsetting. `--rep` selects the SAE
input (`expression`, `embedding` .npy, `pca`, `nmf`), so the representation ladder reuses the same
metric. Validation: an oracle unit test (planted regulators grounded; null, global and background
perturbations rejected) and a pipeline self-test with a real SAE (8/8 planted regulators, 0 false
positives, null 0.000; all-cell mode).

**A method idea that failed, honestly shelved.** A causally-supervised SAE (decoder-side
reconstruction of held-out perturbation signatures added to the loss) did not improve held-out
grounding on synthetic data and sometimes hurt it (plain 4/4 vs supervised 0/4 on separable
programs). `src/causal_sae.py` was deleted rather than shipped.

---

## 8. Act 3b: Replogle CRISPRi (2026-07-17 to 07-20)

### 8.1 Operations
Rented vast.ai GPU (RTX A5000, 24 GB box). The private repo would not clone; a GitHub token was
pasted in plaintext, so the user was told to revoke it and the repo was made public (no token
ever goes on rented compute). Out-of-memory thrash: the dense 310,385 x 8,563 matrix (~10 GB)
blew up the box; fixed with `--n-hvg`.

### 8.2 The broken first run and two real bugs
The first run grounded **0/1789**, below the shuffle null (0.009): a broken metric, not a negative.
(1) An SD-standardized drop is unreachable for sparse TopK activations (a feature fires in ~k/L of
cells, so even a full shutoff is ~0.1 SD); replaced by Mann-Whitney AUC. (2) Double dipping (the
same cells selected and validated the feature); fixed by cross-fitting. Scoring was also restricted
to TFs (`config/tf_lists/hs_hgnc_tfs.txt`, 1,839 pySCENIC TFs).

### 8.3 The causal-grounding spec chain

| Version | SHA | Change | Fate |
|---|---|---|---|
| v1 | `b7d28fae` | SD-standardized drop, no cross-fitting | broken (0/1789); nothing reported |
| v3 | `b766d21c` | cross-fit Mann-Whitney AUC (the internal name skips v2) | first valid real runs |
| v4 | `029ab066` | + column specificity | failed on real data (cut GATA1) |
| v5 | `98104b2f` | column specificity default OFF | rollback |
| v6 | `81b44e04` | + effect-size control (M1) | a no-op on real data |
| v7 | `ee5dacb7` | + direction (down/up) | used for all Norman runs |
| v7 + erratum | `aae1a9d9` | 2026-10-18: `trained_on` corrected to what ran (all cells) | current |

### 8.4 The first real numbers
Replogle 2022 K562 essential library (`pertpy.data.replogle_2022_k562_essential`): 310,385 cells
x 2,136 genes (top-2000 HVG UNION perturbed TFs), 162 testable TFs, latent 2048, k 32.
- Seed 0: **22/162 = 0.136**, shuffle null mean 0.001, p = 0.048 (the 20-shuffle floor). Gates:
  match 162/162, floor 35, mw_q < 0.05 24.
- Seeds 0-4: 22 / 31 / 42 / 31 / 32 -> mean **0.195**, sd 0.039. Every seed far above the null.
- Seed-stable set (>= 3/5 seeds), 26 genes: CNOT3, CSNK2B, CXXC1, DMAP1, E4F1, **GATA1**, GTF2A2,
  HINFP, HSF1, ILF2, MAX, NCBP2, NELFB, POLD2, PSMD12, RFC2, RFC3, SFPQ, SMUG1, SRP9, TBP, TFDP1,
  THOC2, TIMELESS, ZMAT2, ZNF335. About 20 are housekeeping machinery (basal transcription,
  replication, RNA processing, proteasome); 6 are sequence-specific TFs; GATA1 is the only
  lineage TF.
- Hyperparameter sweep (seed 0): latent 512 -> 0.395, 1024 -> 0.284, 2048 -> 0.136 (k 32); k 16
  -> 0.191, k 64 -> 0.358 (latent 2048). A 3x range from knobs alone: the absolute rate is not a
  quotable number; only a matched-hyperparameter comparison is.

### 8.5 The essential-gene confound, claimed and then corrected
Claimed: the seed-stable set is ~10x enriched for Hart CEGv2 core-essential genes (8/26 vs a 3.0%
base rate, hypergeometric p ~4e-7) and depleted of NEGv1 nonessentials (0/26). Wrong background:
the base rate used all 1,839 TFs instead of the 162 tested. Against the right background (the v4
tag split), essential 6/51 = 11.8% vs non-essential 14/111 = 12.6%: **no essential preference**.
The Replogle "essential" library is essential-gene-dominated by construction. Two falsifiable
predictions (P1, sweep flattening; P2, higher seed stability) were recorded before v4 but never
cleanly tested, because v4 failed.

### 8.6 Column specificity (v4) failed
v4 grounded 20-34/162 (mean ~0.157), the same housekeeping set. GATA1 fell from >= 3/5 seeds to
2/5; SPI1, TAL1, RUNX1, KLF1, CEBPA, GFI1, GATA2, MYB were 0/5. Cause: each essential knockdown
suppresses its own distinct stress feature (so the gate passes it), while GATA1's erythroid program
in K562 has several regulators (so GATA1 looks non-specific and is cut). Rolled back.

### 8.7 The effect-size control (M1) changed nothing, which was the finding
Raw 22/162 = controlled 22/162 (0 dropped; the gate is an AND-filter, so "rescue" was never
possible, a framing error of mine). Spearman rho(effect size, suppression) = 0.28; mean effect
grounded 3.81 vs 3.07. GATA1 effect-size p = 0.435; MYC 0.018 but fails the floor. **Effect size is
not the confound.** This refutes the pitch written into `docs/IMPACT_STRATEGY.md` the same day.

### 8.8 The representation ladder: the SAE is not special
SAE on expression 0.121 (seeds 0.136 / 0.093 / 0.136) vs on PCA-128 0.130 (0.130 / 0.111 / 0.148)
vs on NMF-128 0.111 (seed 0): flat. The grounded sets differ: expression and PCA share 6 (chance
~2.85), expression and NMF 2. A weak shared core plus a method-dependent tail, mirroring Gate 0.

### 8.9 RPE1, a second cell line
Expression 0.076, PCA 0.070, NMF 0.076: flat again, and lower than K562. Grounded set again
housekeeping and chromatin (RPL35, HMGB1, HMGA1, SMARCC1, GABPB1, THAP11, ZNF720), zero overlap with
K562's set, no lineage TFs (0/3 seeds).

### 8.10 Synthesis of the CRISPRi arm
Expression-space SAE grounding is real (above a stringent null), effect-size controlled, flat
across decompositions, cell-line dependent, and always lands on whichever housekeeping programs
dominate the dataset. It never recovered lineage regulation. Read with the 2026-10-18 erratum: in
every one of these runs the SAE was fit on perturbed as well as control cells, so "programs learned
from natural variation behave causally" was never what was tested; "features learned from a mix of
perturbed and control cells respond to specific perturbations" was. The embedding head-to-head
(Arm A), the centerpiece of the original plan, was never run.

---

## 9. Act 4: overexpression and Norman CRISPRa (2026-07-21 to 08-20)

### 9.1 Why the pivot
The Replogle libraries knock down essential genes and do not perturb lineage TFs: a
wrong-substrate null, not proof the method cannot find regulatory logic. The intended substrate was
Joung et al. 2023, the TF Atlas (GSE216481): ~3,550 TF ORFs overexpressed in human stem cells,
~254,519 cells, with independent proof that TF-specific signal exists (Jain and Sharma 2026, arXiv
2604.02511: 59 of 61 testable TFs recovered, vs 27 by naive one-vs-rest). Its raw GEO deposit is a
16.7 GB TAR, too heavy for the setup, so the project switched to **Norman et al. 2019** (K562
CRISPRa, dual-guide, via pertpy), which overexpresses KLF1, CEBPA, SPI1, ETS2 and others.

### 9.2 Directional metric and the TopK competition artifact
`direction="up"` matches the +DE program, uses Mann-Whitney "greater", floor AUC >= 0.55, and the
upper tail for column specificity; everything else is unchanged. Synthetic: on planted activations
`up` grounds 6/6 and wrong-sign `down` 0/6. Subtlety: in a TopK SAE an overexpressed program wins a
top-k slot and pushes competitor features out, so they look "suppressed"; the wrong-sign `down`
then grounds competitors. On Norman, `down` gave ~0.363 vs `up` ~0.121; `down` is the artifact.
Documented in `docs/METHODOLOGY_ADDENDUM_OE.md`.

### 9.3 The data, verified
The file every reported Norman result used is pertpy's **`norman_2019_raw`** (re-verified
2026-10-18 by reading both pertpy files' schemas): 111,445 cells x 33,694 genes, integer counts,
Ensembl var_names plus a `gene_symbols` column, `guide_identity` plus 106 `guide_<GENE>` columns, no
`perturbation_name`. (pertpy's processed `norman_2019` has 111,255 cells, 19,018 genes, symbol
var_names, a `perturbation_name` column, and log-normalized X with raw counts only in
`layers['counts']`.) After `prep_norman.py`: **11,855 control cells**, 57,831 single-gene cells,
41,759 dual-gene cells, 105 distinct single genes, **102 testable** (>= 30 cells, symbol in the
matrix). KLF1 1,960 cells; ETS2 1,201; SPI1 613; CEBPA 460. **GATA1 and RUNX1 are not in the
library** (0 cells); RUNX1T1 is.

Three prep failures on vast.ai before it worked: (1) the schema guard expected `perturbation_name`;
(2) Ensembl var_names vs symbol labels gave 0 testable perturbations and a ZeroDivisionError inside
`multipletests([])` (fixed: empty guard `a356749` + symbol remap); (3) setting var_names from a
column of the same name broke the h5ad write and corrupted the file (X came back None) (fixed:
delete the column, assign a list). A harmless "variable names are not unique" warning remains
(duplicate symbols after the remap).

### 9.4 The first run, logged as a provisional positive
With no null (`--n-shuffle 0`): up-rate mean ~0.121; KLF1 grounded in 2/3 seeds, ETS2 3/3;
seed-stable set KLF1, SLC4A1, ETS2, RUNX1T1, CBFA2T3, hematopoietic rather than housekeeping.

### 9.5 The null run (2026-08-10): negative
Latent 2048, k 32, `--n-hvg 4000` (4,047 genes), 102 testable, direction up:

| Seed | Grounded | Shuffle null | Gates passed (match / floor / FDR / effect-size) | rho(effect, movement) |
|---|---|---|---|---|
| 0 | 15/102 = 0.147 | mean 0.175, **p = 0.902** | 102 / 87 / 96 / 15 | -0.80 |
| 1 | 13/102 = 0.127 | not run | 102 / 94 / 98 / 13 | -0.81 |
| 2 | 20/102 = 0.196 | not run | 102 / 94 / 98 / 20 | -0.78 |

Median matched AUC ~0.67-0.68. Grounded perturbations had smaller total effects than the rest
(2.13 vs 2.75, 2.03 vs 2.75, 2.38 vs 2.73). Seed-stable across all three seeds: BAK1, ETS2, KLF1,
PTPN12, SLC4A1, TBX2, UBASH3B, ZBTB25, MAP2K3, MAP2K6.

**Inferred mechanism of the null** (from the gate counts above; not separately tested): the match,
floor and FDR gates pass 85-100% of perturbations, so the grounding rate is set almost entirely by
the effect-size gate. That gate is relative (robust-z among effect-matched peers, alpha 0.10), so it
passes a roughly fixed upper tail whatever the labels are, and the shuffle reproduces ~0.175. On
Replogle the floor and FDR gates were binding (35 and 24 of 162), so grounding there tracked real
suppression and beat the null. Corollary for any future use: a grounding rate is only meaningful
when the absolute gates bind.

### 9.6 What the features actually are (the `top_program_genes` readout)
- **Erythroid/hemoglobin program** (features 537 / 630 / 1982 in seeds 0 / 1 / 2): HBZ, ALAS2, HBG2,
  SLC25A37, HBA1, FTH1, HBG1, GYPB, BLVRB (plus HIST1H1C in seed 2). The same program in every seed:
  real, hand-verifiable biology. Grounded by **CBL** (seeds 0, 1), **PTPN12**, **UBASH3B**, **ZBTB25**
  (all three seeds) and ZNF318 (seed 2), not by KLF1.
- **KLF1's feature:** TMSB10, PNMT, KLF1, S100A10 / EMP3, HBZ, ACP5, S100A6, GAL (self-dominated,
  one globin).
- **ETS2's feature:** AIF1, ETS2, SH3BGRL3, TMSB4X, CD52, ARHGDIB, AP1S2, GMFG (myeloid-like);
  MAPK1 shares it in seed 2.
- **Stress/proteostasis feature:** PRDX1, LDHB, PSMA4, UBB, CACYBP, LDHA, HSP90AA1 (MAP2K3, MAP2K6,
  CSRNP1, SLC38A2): the Replogle-style housekeeping grounding again.
- **IRF1** (seed 2): PSMB8 and PSMB9 (immunoproteasome, textbook IRF1 targets) among its top genes.
- TBX2/TBX3: AVP, GATA2, TBX3; BAK1: LHX1, KLRB1, CSF2 (no coherent biology); SLC4A1 (erythroid
  Band 3) grounds a non-erythroid feature (KRT79, IFITM2, IRX3).

Interpretation: the SAE isolates K562's dominant axes (erythroid, stress) cleanly, and the matching
step grounds almost any grouping of cells to one of them. **Correction (2026-10-18):** the August
note that the hemoglobin feature was "NOT [grounded] by GATA1 (doesn't ground)" was wrong as
stated. GATA1 was never perturbed in Norman; it could not ground.

### 9.7 Regulon recovery v1 (`2f5578a`, spec `bbfac449`): uninformative
Design: for each overexpressed TF with a TRRUST regulon (>= 5 targets in panel, self excluded),
select the feature whose top-50 control-cell gene associations are most enriched for those targets
(hypergeometric, Bonferroni over features), then test that feature's activation in the TF's cells
(Mann-Whitney, AUC >= 0.55), BH on both, label-shuffle null. Result: **0/10 recovered, null 0.00,
p = 1.000, all 3 seeds**. The diagnostic (`src/inspect_regulon_json.py`) showed the test never
fired: enrichment passed 0/10 in every seed; 21 TFs were low-power (too few targets in the
4,000-HVG panel); and promiscuous TRRUST hub targets (CDKN1A regulated by 115 TFs, MYC 74, CCND1 64,
VEGFA 58, TP53 51, CDH1 50, BCL2 46, AR 46, PTGS2 43) pulled every TF's selection onto the same
generic stress feature (all at AUC ~0.50). Power arithmetic: with 10 targets in a 4,047-gene panel,
a feature needs >= 4 of them in its top 50 to survive Bonferroni over 2,048 features.

### 9.8 Regulon recovery v2 (`8bc8734`, spec `c653f6d1`): a real negative, line closed
Fixes: force-keep regulon targets in the panel (`hvg_keep_mask`) and drop hub targets regulated by
more than 10 TFs (133 of 2,492). Panel 4,047 -> 5,791 genes; tested pool 10 -> 14 (AHR, CEBPA,
CEBPB, CEBPE, EGR1, ETS2, FOXA1, HNF4A, IRF1, JUN, POU3F2, SNAI1, SPI1, TP73); 17 still low-power
(including KLF1, CBFA2T3, LYL1, MEIS1, PRDM1). Result: **0/14, every enrichment q = 1.0, all seeds**.
Some selected features were strongly activated (CEBPB seed 1 activation q ~1e-127, CEBPE seed 2
~5e-64, SPI1 seed 2 ~2e-38) but did not carry the regulon. Honest notes: (a) the v2 amendment
expected the force-keep to make the hematopoietic TFs testable; it did not, because their TRRUST
regulons are tiny (KLF1 has 2 targets), so the binding limit was TRRUST coverage, not the panel;
(b) per the v2 spec's own pre-commitment, the line is closed (no v3).

**Chance comparison, corrected 2026-10-18.** Simulating the selection itself (the best of 2,048
random top-50 gene sets from a 5,791-gene panel, `src/regulon_chance_baseline.py`), per TF and
seed from the pasted v2 output: **8 of 42** TF x seed cells sit above the chance distribution's mean,
**0 of 42** above its 95th percentile (SPI1 seed 0 sits exactly at it; its printed enrichment q of 1
implies exactly 5 hits). The August wording, "at or below chance for every regulon size", was too
strong; the conclusion stands: SAE features do not concentrate curated single-TF regulons beyond
chance. The script has not yet been run on the actual v2 JSONs (they were on the destroyed vast.ai
instance). Scope: this concerns single-TF curated-regulon alignment only; it does not contradict
Gate 0's lineage-level TRRUST positive (union of a lineage's TF targets, 512 latents, marrow data).
Precise synthesis: regulon alignment is detectable at lineage resolution and absent at single-TF
resolution with per-feature multiple-testing control. Likely reason (a hypothesis): TRRUST edges
are curated across many cell types; SAE features are co-expression modules in one cell line.

### 9.9 Program modulation: built, not yet run (the decisive test)
`src/program_modulation.py` (`4a3728e`, spec v1 `b396834d` -> v1.1 `ddb29bcd`). The question moves
from single-perturbation to program resolution, where Gate 0 found signal: is one externally
defined program causally modulable by a specific subset of perturbations?
- **Non-circular by construction.** Programs come from `config/human_markers.json`, curated for
  Gate 0 before any CRISPRa work. The feature is selected over control-cell associations only. Which
  perturbations modulate it is an output; CBL, UBASH3B and PTPN12 are never inputs.
- **Test.** For each of the 9 programs, pick the eligible feature (fires in >= 50 controls) most
  enriched for its markers; for each of the 105 perturbations (>= 20 cells), one-sided Mann-Whitney
  vs controls; BH over perturbations (not 2,048 features, which is why it has power where regulon
  recovery did not); modulates iff q < 0.05 and AUC >= 0.55.
- **Two nulls, both required.** A label shuffle (default "full": permute all labels including
  controls) and a feature-background null (the identical test on 100 random eligible features; the
  program's percentile). A positive needs shuffle p < 0.05 AND background percentile >= 90 AND a
  biologically coherent modulator list. The Gran, Progenitor and Cycling rows are built-in
  specificity controls.
- **Oracle** (`tests/test_program_modulation.py`): a specific world is detected (all true movers
  found, shuffle p 0.032, background 97th percentile); a dominant-axis world (one feature moved by
  nearly everything) reports 18 modulators but sits at the 20th percentile and reads as a negative;
  false-positive rate 0.031 over 8 seeds (nominal 0.05); feature selection invariant to label
  shuffling.
- **Smoke-tested on the real data 2026-10-18** with stand-in activations: 111,445 cells x 4,081
  genes, all 59/59 markers in the panel, 11,855 controls, 105 perturbations; analysis completes
  in seconds; random activations correctly give 0 modulators and p = 1.00. That checks plumbing
  only. **No real result exists in this repo.**
- **v1.1 amendment** (before any recorded run): primary analysis fits the SAE on all cells, as the
  frozen code does and as the only validated mode; control-only is a pre-declared sensitivity
  analysis; input must be `norman_2019_raw`.

### 9.10 Two null-design and training-design corrections
- **The shuffle null (2026-08-20).** Permuting labels among perturbed cells with controls fixed
  destroys perturbation identity but keeps the perturbed-vs-control contrast. For a
  modulator-count statistic that inflates the null (oracle: null 14.4 vs real 5 on a world with 4
  true modulators). For the Norman grounding rate it means p = 0.902 answers a specific question:
  perturbation identity adds nothing beyond the generic perturbation response. It does not show
  grounding is at absolute chance. That is still the informative negative for an identity claim.
  `program_modulation` defaults to the full permutation and relies on the feature-background null
  for the identity/dominant-axis concern.
- **Training cells (2026-10-18).** `config/causal_grounding_spec.json` said "trained on control
  cells only; the SAE never sees a perturbation". Every caller passed the full matrix, so every
  real run fit on all cells (labels unseen, so "labels never touch training or selection" holds;
  "never sees a perturbed cell" does not). On the pipeline's synthetic world, all-cell training
  grounds 8/8 planted regulators; control-only grounds **0/8** (signal present, matched AUC ~0.40
  vs 0.31, match-confidence p 0.04-0.16, best q ~0.078, just over the FDR cut). So the pre-registered
  design was never validated, and in this synthetic world it is underpowered. Handling: dated
  erratum in the spec (SHA `ee5dacb7` -> `aae1a9d9`); `train_rows()` + `--train-on {all,control}`
  with the default `all` (results byte-identical to before); `tests/test_pipeline_training.py`;
  erratum notes in the method docs. No past result changes sign; their interpretation narrows.

---

## 10. Scorecard: what is established, refuted, open

| Claim | Status | Evidence |
|---|---|---|
| Asymmetric modularity (granulocyte unified, erythroid distributed) | **Refuted** | 0/10 mouse, 0/10 human; survives in 1 of 10 method x species cells |
| The reversed direction (granulocyte more distributed) is biology | **Not established** | SAE-consistent in both species, but NMF disagrees in both: method-contingent |
| SAE-program modularity is method-dependent | **Established, frozen** | 5 decomposition families x 2 species |
| The Gate-0 negative is robust | **Established** | R1-R4 in both species; leave-one-marker-out 0/34 and 0/36 |
| SAE features recover lineage-level curated regulons (human) | **Supported (structural)** | TRRUST erythroid, 10/10 seeds, median p 6.4e-5 |
| SAE features align with single-TF curated regulons | **Not supported (quantified)** | Regulon recovery v2: 0/14; 0/42 cells beyond the chance 95th percentile |
| Expression-SAE causal grounding beats a stringent null (Replogle K562) | **Supported, housekeeping only** | 22/162 vs null 0.001; 5 seeds mean 0.195. Caveat: SAE fit on all cells |
| Grounding recovers lineage regulators on CRISPRi | **Not supported** | GATA1 the lone hit, lost under v4, effect-size p 0.435; no other master TF |
| The SAE adds grounding over PCA/NMF | **Not supported** | Ladder flat in K562 and RPE1 |
| Grounding is an essential-gene artifact | **Refuted** (the 10x stat was wrong-background) | 11.8% vs 12.6% |
| Grounding is an effect-size artifact | **Refuted** | M1 dropped 0 of 22; rho 0.28 |
| Overexpressed-TF identity matters for grounding (Norman) | **Not supported** | 0.147 vs stringent null 0.175, p 0.902 |
| A seed-stable erythroid/hemoglobin SAE program exists in K562 | **Established (descriptive)** | Same top genes in 3/3 seeds |
| CBL / UBASH3B / PTPN12 overexpression activates that program | **Hypothesis (post hoc, untested)** | They ground to it in Norman; mechanistically plausible (KIT/EPOR brakes) |
| Expression vs foundation-model embedding head-to-head (Arm A) | **Not run** | Pipeline supports `--rep embedding`; no embedding generator exists |
| Program modulation (the erythroid test) | **Not run** | Built, oracle-tested, pre-registered v1.1, smoke-tested on real data |
| Causally-supervised SAE | **Failed** | Synthetic: 0/4 vs plain 4/4; shelved |
| Control-only SAE training (the original design) | **Unvalidated** | Synthetic self-test 0/8 vs all-cell 8/8 |

---

## 11. Every bug, self-correction and honesty save (chronological)

This list is itself evidence for the project's framework claim: each item is a place where a
plausible result turned out to be an artifact, and how it was caught.

1. Pseudotime rooted inside a committed MEP (`iroot=0`, W31105) -> Palantir rooted at 10GMP (V0a).
2. V0b v1: cross-seed feature-index matching, `prob_19Lymph` omitted (~273 cells mislabeled),
   markers absent from the HVG set, biased entropy ratio -> v2.
3. The notebook trap: the canonically named, latest Drive notebook is v1 -> the repo script is
   canonical; never run the Drive notebooks.
4. `torch.manual_seed` NameError shipped without running (`ca1a746`) -> CLAUDE.md run-before-send.
5. V0b v2's NOT SUPPORTED came from a structurally blind test -> recorded as an artifact.
6. V0b v3's SUPPORTED came from count bias, noise-floor submodules and a liberal null -> v3.1.
7. "Paul15 is pre-log-transformed" is false (raw max 168) -> `--log1p auto` (audit A11).
8. Wrong globin symbols and missing Alas2/Ermap (my curation error) -> corrected; result reversed.
9. Under-regularized SAEs (L0 ~213, then ~63) -> L1 retuned at the real epoch count.
10. Overstated novelty ("SAEs on single-cell is novel") retracted; later the opposite overstatement
    ("the project is done") retracted once all competitors were verified to be embedding-space.
11. Pre-registration hash mismatch (trailing newline) -> canonical bytes; content unchanged.
12. Non-deterministic figures -> byte-identical reruns.
13. ICA at rank 512 did not converge and falsely flipped human to "supported" -> rank 50.
14. The mouse TRRUST null was kept out of the evidence (human database, uppercased symbols).
15. The causally-supervised SAE did not validate -> code deleted, idea shelved.
16. Response-based matching was circular (caught by the oracle) -> program-based matching.
17. Column specificity alone could not reject a global stressor (caught by tests) -> match-confidence gate.
18. A GitHub token pasted in plaintext -> revoked; repo made public.
19. Out-of-memory on the dense genome-scale matrix -> `--n-hvg`.
20. SD-standardized suppression unreachable for TopK activations (0/1789) -> Mann-Whitney AUC.
21. Double dipping -> cross-fitting.
22. Scoring all 1,789 perturbations including housekeeping -> `--tf-list`.
23. "~10x essential enrichment, p ~4e-7" used the wrong background -> no essential preference.
24. Column specificity (v4) removed GATA1 and kept the confound -> rolled back (v5).
25. M1 described as a possible "rescue" (an AND filter can only remove) -> corrected.
26. "Grounding just tracks effect-detectability" -> refuted by M1.
27. Norman prep failure 1: the schema guard expected `perturbation_name` -> parse `guide_identity`.
28. Norman prep failure 2: Ensembl var_names vs symbol labels -> 0 testable -> ZeroDivisionError in
    `multipletests([])` -> empty guard (`a356749`) + symbol remap.
29. Norman prep failure 3: var index name equal to a column name -> corrupt h5ad, X None -> fixed.
30. "The GPU is not being used" -> the trainer does use CUDA; the pipeline is CPU-dominated.
31. Container resets dropped numpy/anndata mid-session (recurring) -> reinstall with
    `python3 -m pip`; anndata needs `--ignore-installed packaging` on this image.
32. Norman "first genuinely encouraging result" -> retracted when the null came back p = 0.902.
33. Regulon recovery v1 never fired -> diagnosed with `inspect_regulon_json.py` -> v2.
34. Program-modulation oracle demanded an exact set from a stochastic test -> calibrated (FP 0.031).
35. Perturbed-only label shuffle inflates count nulls; Norman p = 0.902 re-scoped precisely ->
    full-permutation default for program modulation.
36. **(2026-10-18)** The SAE was fit on ALL cells in every real run, contrary to the frozen spec;
    control-only fails the synthetic self-test (0/8 vs 8/8) -> spec erratum, `--train-on`, tests.
37. **(2026-10-18)** My 2026-08-20 re-download line `pt.dt.norman_2019()` was sent without being
    run and points at the log-normalized file -> correct call `pt.dt.norman_2019_raw()`; August
    used raw (verified); `prep_norman.py` now refuses or repairs non-raw X; `tests/test_prep_norman.py`.
38. **(2026-10-18)** "GATA1 doesn't ground" in Norman -> GATA1 is not in the library (0 cells).
39. **(2026-10-18)** Regulon chance wording too strong -> 8/42 above the mean, 0/42 above the 95th percentile.
40. **(2026-10-18)** The v2 amendment's expectation (force-keep makes hematopoietic TFs testable)
    was wrong: TRRUST coverage, not the panel, was the limit.
41. **(2026-10-18)** The pipeline self-test had never run in this container (no torch) -> run today
    with PyPI torch on CPU (`download.pytorch.org` is policy-blocked here).
42. Open, cosmetic: the Norman diagnostic line prints "floor(auc<=0.45)" and "suppression" even for
    `direction=up` (labels only; values correct).
43. Open, cosmetic: the mouse decision JSON's `plain_finding` text says "erythroid (2)" for a median of 1.5.
44. Open, minor: duplicate gene symbols after the Norman remap ("variable names are not unique").
45. Open: PR #3's description and several docs are stale (section 17).

---

## 12. Frozen artifacts registry

All verified with `sha256sum` on 2026-10-18.

| Artifact | SHA-256 | Freezes | State |
|---|---|---|---|
| `config/preregistration_spec.json` | `fc342829191da5a0bb5244975c3f61c8934943cb512818446e5b5f0984f6011a` | Gate-0 metric, markers, thresholds, decision rule, Phase 2 null | intact |
| Mouse Gate-0 bundle (`data_g0/v0b_outputs/v0b_v3_1_decision.json` + `v0b_v3_1_per_seed.csv`) | `172861417a923a83705812a96e5ad1e555055e968323289879447e53beb84396` | executed mouse verdict | intact |
| Human Gate-0 bundle (`data_g0_human/v0b_outputs/...`) | `cc8159c882cf326bca5bbbc7f74d110d4a81338ec0ba07e208eb935bc8f729d3` | executed human verdict | intact |
| `config/causal_grounding_spec.json` | `aae1a9d9c6caaca00130cc4a71b37d135c6b7a2d45d18285e18601ba80fe128c` | causal metric v7 + amendment chain + 2026-10-18 erratum | current |
| `config/regulon_recovery_spec.json` | `c653f6d19e4c464d14fb8a6869d03fc490191720b7a7c939ae13c012fc4daed7` | regulon recovery v2 (line closed) | intact |
| `config/program_modulation_spec.json` | `ddb29bcdba1d9470af8d6d16b115e0a5ba0e0bf63fb2438b1e2230df17420325` | program modulation v1.1 (not yet run) | current |

Superseded hashes, kept for the record: causal spec v1 `b7d28fae`, v3 `b766d21c`, v4 `029ab066`, v5
`98104b2f`, v6 `81b44e04`, v7 `ee5dacb7`; regulon v1 `bbfac449`; program modulation v1 `b396834d`;
first V0b v2 assignments (not frozen) `f18b088a...`; Phase 1 annotation table `bd531092...`; Phase 1
matrix MD5 `60183a17...` (superseded by the marker-aware mouse MD5 `184ea5dd...` and human
`62796189...`); source data Paul15 `bd75eb90...`, Setty `be4f8603...`; TRRUST human file `a701172d...`
(per the notebook).

**Deliberately not frozen:** every Component 2 result (Replogle K562 and RPE1, ladder, M1, Norman,
regulon recovery) and every Gate-0 intermediate run. Raw Component 2 JSONs lived on destroyed
vast.ai instances and are regenerable from pinned code and seeds.

---

## 13. Complete file index

**Root**
- `README.md`: entry point and repo map (status table stale; section 17).
- `CLAUDE.md`: operating manual for agents: rules, conventions, status (sections 3 and 5 stale).
- `LAB_NOTEBOOK.md`: dated research log, 2026-06-29 to 2026-10-18.
- `requirements.txt`: numpy 2.4.6, pandas 3.0.3, scipy 1.17.1, statsmodels 0.14.6, torch >= 2.0,
  scanpy >= 1.10. pertpy and anndata are installed separately for the Perturb-seq work.
- `.gitignore`: large binaries out, small result artifacts in; `causal_out/` ignored.

**`config/`**
- `preregistration_spec.json` + `.sha256`: the Gate-0 pre-registration.
- `causal_grounding_spec.json` + `.sha256`: the causal metric, its amendment chain, the erratum.
- `regulon_recovery_spec.json` + `.sha256`: regulon recovery v2 (closed).
- `program_modulation_spec.json` + `.sha256`: program modulation v1.1 (pending).
- `human_markers.json`: human ortholog marker panel (Gate-0 human arm; program modulation).
- `gene_sets/hart_cegv2_core_essential.txt` (684), `hart_negv1_nonessential.txt` (927),
  `gene_sets/README.md` (still carries the uncorrected ~10x claim).
- `tf_lists/hs_hgnc_tfs.txt` (1,839 human TFs), `mm_mgi_tfs.txt` (1,721 mouse), from aertslab/pySCENIC.
- `trrust/trrust_rawdata.human.tsv`: TRRUST v2 human, 9,396 edges, 795 TFs.

**`data/README.md`**: manifest for the Phase 1 Drive inputs (refers to the superseded 128-latent
matrix and its MD5; lists Drive folder IDs).

**`data_g0/` (mouse) and `data_g0_human/` (human): executed Gate-0 evidence** (large matrices and
checkpoints gitignored): `gene_names.csv`, `preprocess_report.json`, `training_report.json`,
`v0b_outputs/v0b_v3_1_decision.json` + `v0b_v3_1_per_seed.csv` (frozen bundles),
`baselines/baselines_results.json`, `sensitivity/marker_sensitivity.csv` + `_summary.json`,
`robustness/robustness_decision_params.csv` + `_summary.json`, `r2_l0_band_summary.json`,
`r3_20seed_summary.json`, `r4_hvg_summary.json`, `grn/grn_baseline.json`, `ica/ica_baseline.json`,
`trrust/trrust_crossref.json`.

**`docs/`** (22 files with this one)

| Doc | What it is | Current? |
|---|---|---|
| `PROJECT_MASTER_SUMMARY.md` | This document | current |
| `PROJECT_HANDOFF.md` | Phase 1 closeout snapshot (historical) | historical by design |
| `PROJECT_AUDIT.md` | Corrections and rigor log, sections A-K (Gate 0) | current through Gate 0 |
| `STRATEGY_AND_POSITIONING.md` | Field situating, novelty audit, the reframe | framing partly superseded |
| `PREREGISTRATION.md` | Human-readable Gate-0 pre-registration | current |
| `COMPONENT0_STATUS.md` | Gate-0 status + Colab runner | historical |
| `COMPONENT0_RESULTS.md` | Every executed Gate-0 number + frozen hashes | current |
| `COMPONENT0_HARDENING.md` | R1-R4, TRRUST, GRN, ICA | current |
| `VERSION_HISTORY.md` | Metric lineage v1 -> v3.1 | current for Gate 0 |
| `REPRODUCIBILITY.md` | Scripted verification + expected hashes | current for Gate 0 |
| `ABSTRACT.md` | One-page competition summary (Gate 0) | results current; "Next" stale |
| `web/result_brief.html` | Self-contained visual brief (Gate 0) | current for Gate 0 |
| `ELEVATION_PLAN.md` | The causal pivot | defense 1 erratum added |
| `COMPONENT2_PLAN.md` | Early causal design (cosine version) | superseded by the built metric |
| `RUNPOD_EXECUTION.md` | GPU steps for the head-to-head | partly stale |
| `COMPONENT2_RESULTS.md` | Replogle results, confound saga, M1, ladder | stops before RPE1/Norman; erratum added |
| `COMPONENT2_NEXT_DIRECTIONS.md` | Direction map after the first run | partly stale |
| `UPGRADE_BACKLOG.md` | Tiers 0-6 + reframes | partly stale |
| `DIFFERENTIATION.md` | Verified competitors + rehearsed novelty answer | table claims stale |
| `IMPACT_STRATEGY.md` | Making the negative influential | its pitch is refuted |
| `METHODOLOGY_ADDENDUM_OE.md` | The overexpression pivot + `direction=up` | current; erratum added |
| `COMPUTE.md` | CPU vs GPU reasoning | current |
| `phase2_research_plan_v6.md` | Research plan v6 | stale (section 17) |
| `figures/fig_method_dependence.*`, `fig_robustness.*` | Gate-0 figures | current |
| `figures/fig_component2_diagnostics.*` | Component 2 diagnostics | panel B shows the retracted 10x claim |

**`src/`** (22 scripts)

| Script | Purpose |
|---|---|
| `preprocess_paul15.py` | Mouse marker-aware preprocessing (`--log1p auto`, h5 fallback loader) |
| `preprocess_human_marrow.py` | Human Setty preprocessing (normalize + log1p, ortholog markers) |
| `train_sae.py` | Gate-0 L1 SAE trainer (ReLU, AdamW, cosine schedule, resampling, `--seed-start`) |
| `v0b_module_definitions.py` | V0b v2 + shared loaders and markers (`V0B_MARKERS_JSON` override) |
| `v0b_v3_loading.py` | V0b v3 continuous decoder loading |
| `v0b_v3_1_decision.py` | The pre-registered Gate-0 metric |
| `baselines_nmf_pca.py` | PCA/NMF participation ratio + identical modularity |
| `marker_sensitivity.py` | Leave-one-marker-out |
| `robustness_sweeps.py` | R1 decision-knob sweep |
| `grn_baseline.py` | Co-expression GRN as the 4th family |
| `ica_baseline.py` | FastICA as the 5th family |
| `trrust_crossref.py` | TRRUST lineage-level cross-reference (+ `load_trrust`) |
| `make_figures.py` | Deterministic Gate-0 figures |
| `causal_grounding.py` | The causal-grounding metric (v7) |
| `causal_pipeline.py` | TopK SAE trainer, Perturb-seq loader, representations, shuffle null, `train_rows`, CLI |
| `make_component2_figures.py` | Component 2 diagnostic figure (hardcoded numbers; panel B stale) |
| `inspect_perturb_h5ad.py` | Schema inspector that proposes `--pert-col` / `--control-value` |
| `prep_norman.py` | Norman guide parsing, symbol remap, raw-count guard |
| `regulon_recovery.py` | Regulon-recovery test (v2, closed) |
| `inspect_regulon_json.py` | Which gate zeroed a regulon run |
| `regulon_chance_baseline.py` | Best-of-N random-gene-set chance comparison |
| `program_modulation.py` | The pending decisive test |

**`tests/`** (9 suites, all passing 2026-10-18): `test_v0b_logic.py`, `test_v0b_v3_logic.py`,
`test_v0b_v3_1_logic.py`, `test_baselines_logic.py`, `test_causal_grounding.py`,
`test_regulon_recovery.py`, `test_program_modulation.py`, `test_prep_norman.py` (new),
`test_pipeline_training.py` (new; trains real SAEs when torch is present).

---

## 14. Datasets and external resources

| Resource | What | Size | How obtained | Used for | Status |
|---|---|---|---|---|---|
| Paul et al. 2015 (Cell) | mouse myeloid progenitors, 19 clusters | 2,730 cells | git-LFS scanpy mirror (direct download proxy-blocked), sha256 `bd75eb90...` | Phase 1, Gate 0 mouse | done |
| Setty et al. 2019 (Nat Biotech) | human CD34+ marrow (Palantir repo) | 4,142 cells after preprocessing | git-LFS, `marrow_sample_scseq_counts.h5ad`, sha256 `be4f8603...` | Gate 0 human | done |
| CELLxGENE Census bone marrow | human, 759k pool | 25,000-cell subsample | Census (release date never pinned) | Phase 1 human only | superseded |
| Replogle et al. 2022, K562 essential | CRISPRi Perturb-seq | 310,385 cells | `pertpy.data.replogle_2022_k562_essential()` | Component 2 | done (housekeeping-only) |
| Replogle et al. 2022, RPE1 | CRISPRi Perturb-seq | size not recorded | pertpy (`replogle_2022_rpe1`) | replication | done |
| Norman et al. 2019 (Science) | K562 CRISPRa, single + dual guides | 111,445 cells x 33,694 genes; 102 testable | `pertpy.data.norman_2019_raw()` (= scverse S3 `norman_2019_raw.h5ad`, 828,160,302 bytes) | overexpression arm | done (negatives); program modulation pending |
| Joung et al. 2023 (Cell), GSE216481 | hESC TF ORF overexpression | ~3,550 TFs, ~254,519 cells | GEO raw TAR, 16.7 GB | intended substrate | not used (too heavy) |
| TRRUST v2 human | curated TF -> target edges | 9,396 edges, 795 TFs | GitHub mirror (grnpedia.org blocked) | Gate-0 cross-ref; regulon recovery | committed |
| TRRUST mouse | | | all mirrors 404 | | unavailable |
| Hart CEGv2 / NEGv1 | core-essential / nonessential | 684 / 927 | hart-lab/bagel mirror | confound analysis | committed |
| aertslab TF lists | human / mouse TFs | 1,839 / 1,721 | GitHub | TF restriction, GRN | committed |
| Lambert et al. 2018 TF list | sequence-specific TFs | ~1,639 | humantfs.ccbr.utoronto.ca (egress-blocked) | backlog S.1 | not obtained |
| cisTarget databases | SCENIC motif rankings | GBs | resources.aertslab.org (blocked) | SCENIC baseline | not obtained |
| van Galen 2019 (GSE116256), Hegde 2025 (GSE270148), Park 2020 (E-MTAB-8581) | plan v6 Components 1, 3, 5 | | | | never run |
| scGPT, Geneformer checkpoints | foundation models for Arm A | | | | never run |

---

## 15. Compute, environment, operations

- **This container:** ephemeral cloud VM, 4 CPU cores, 15 GB RAM, no GPU, proxied egress. Reachable:
  PyPI, GitHub, the scverse S3 bucket. Blocked at various times: `download.pytorch.org`
  (policy 403), grnpedia.org, resources.aertslab.org, humantfs.ccbr.utoronto.ca, GEO/figshare.
  Installed packages vanish on resets (twice mid-session). torch installs from PyPI (2.14.1, CUDA
  build, runs on CPU); anndata needs `--ignore-installed packaging`; always use `python3 -m pip`.
- **Measured timings:** Gate-0 SAE ~50 s/seed (mouse), ~75 s/seed (human), so 10 seeds in 8-12 min;
  causal synthetic self-test 23 s; Norman prep 35 s; Norman loader 37 s; Norman TopK SAE step
  332 ms on this CPU, so an all-cell fit (21,800 steps) takes **2.0 h per seed** and a control-only
  fit (2,400 steps) 0.22 h per seed.
- **Rented GPU:** vast.ai RTX A5000 (24 GB box), July-August. Workflow: `nohup` scripts, ntfy.sh
  phone alerts, SHA-pinned checkouts with grep guards, outputs pasted back. The ntfy topic was
  left as the literal placeholder `YOUR_TOPIC_HERE`, so alerts went to a public topic nobody was
  subscribed to; ntfy topics are public, so use a long random string. Costs were not recorded.
  The instances are gone; their raw JSONs with them.
- **Other hardware:** the local AMD RX 5700 XT is not practical (ROCm does not support gfx1010);
  Colab T4 works for small jobs.
- **Security:** the repo is public precisely so no credential ever goes on rented compute. Keep it
  that way.

---

## 16. Field, competitors, and what the differentiation now honestly is

Verified competitor landscape (`docs/DIFFERENTIATION.md`, abstracts checked 2026-07-18; arXiv full
text returned 403, so re-read the PDFs before submission):
- Kendiukhov 2026, arXiv 2603.02952: TopK SAEs on Geneformer V2-316M and scGPT residual streams
  (82,525 and 24,527 features); 29-59% annotate to GO/KEGG/TRRUST; only 3/48 TFs (6.2%; 10.4%
  multi-tissue) show target-specific CRISPRi responses. The project's original null.
- Kendiukhov 2026, arXiv 2603.01752 (causal circuit tracing; 56.4% directional CRISPRi accuracy,
  "co-expression, not causal encoding") and 2603.11940 (exhaustive circuit mapping). Embedding space.
- arXiv 2603.10261 (hematopoietic manifold extracted from scGPT internals): adjacent, not a collision.
- CytoSAE, arXiv 2507.12464: SAEs on blood-cell microscopy images: not a collision.
- Also: bioRxiv 2025.10.22.681631; arXiv 2510.25807; Schuster, arXiv 2410.11468 (closest to
  expression/latent-variable-model space); Jain and Sharma, arXiv 2604.02511 (TF Atlas re-analysis).
- STS template cited in the strategy doc: PLI-Analyzer (an audit showing AI protein-complex
  predictors are ~50% accurate) reached the 2026 STS top 10.

**What held up and what did not, of the five differentiators `DIFFERENTIATION.md` claims:**
(1) expression space rather than embedding space: still true, and still unique in that cluster;
(2) the essential-gene confound and column-specificity metric: no longer a contribution (the 10x
statistic was wrong-background and the gate failed); (3) a within-study expression-vs-embedding
head-to-head: never run; (4) the Gate-0 artifact-vs-signal framework: true and frozen; (5) GATA1 as
a positive control: did not hold. The honest differentiator today is expression-space SAEs plus a
pre-registered framework that has now documented a concrete catalog of ways SAE metrics mislead
on single-cell data: method-dependence (Gate 0), dominant-axis promiscuity and gate saturation
(Norman), hub-target hijacking of regulon selection, wrong-background enrichment, the TopK
competition artifact, and training-cell leakage, each with the diagnostic that exposed it. The
program-modulation test is the one pending positive-shaped claim. The rehearsed novelty paragraph
must be rewritten to match.

---

## 17. Stale or wrong statements in other files (fix before writing anything)

1. `README.md` status table stops at the Arm B negative; it lists the causal spec as v7 `ee5dacb7`
   (now `aae1a9d9` with the erratum) and has nothing on RPE1, Norman, regulon recovery, or program
   modulation.
2. `CLAUDE.md` section 3 says "Paul15 is pre-log-transformed: skip normalize_total/log1p" (false
   for `sc.datasets.paul15()`, audit A11) and "MD5 `60183a17...` MUST be verified at load" (that is
   the superseded Phase 1 matrix; Gate 0 uses `184ea5dd...` and `62796189...` via the dual-source
   check). Section 5 says "Next: Component 2" (it ran).
3. `docs/IMPACT_STRATEGY.md`: its thesis and one-sentence pitch say SAE grounding tracks
   perturbation effect size; M1 refuted that the next day.
4. `docs/DIFFERENTIATION.md`: "essential-gene confound control: yes, quantified (~10x) and gated
   (column specificity)"; "GATA1-anchored"; the rehearsed paragraph cites the confound fix. All
   outdated (section 16).
5. `config/gene_sets/README.md` repeats the uncorrected "~10x enriched, p ~4e-7".
6. `docs/figures/fig_component2_diagnostics.*` panel B is titled "~10x essential-enriched (the
   confound)" from numbers hardcoded in `src/make_component2_figures.py`. Do not use panel B.
7. `docs/COMPONENT2_RESULTS.md` stops at the representation ladder (lists RPE1 as not done; no
   Norman, regulon or program-modulation sections). Its training line now carries an erratum.
8. `docs/COMPONENT2_NEXT_DIRECTIONS.md` and `docs/UPGRADE_BACKLOG.md`: "S.2 kills the confound"
   (it failed); the Part 2 re-run and the GATA1 positive control are superseded.
9. `docs/ABSTRACT.md` "Next (the centerpiece)" describes Component 2 as unrun and GPU-blocked.
10. `docs/phase2_research_plan_v6.md`: states the reversed direction as a finding (it is
    method-contingent); describes the human arm as the 25k CELLxGENE subsample (the executed arm is
    Setty, 4,142 cells); Components 1, 3, 4, 5 never ran; RQ1's 24-TF annotated-feature design was
    replaced by the matching metric.
11. `docs/STRATEGY_AND_POSITIONING.md`: the reframe promises "a corrected, causally-tested result
    (erythroid unified ...)". Nothing structural was causally tested, and the direction is
    method-contingent.
12. `docs/RUNPOD_EXECUTION.md`: cites spec `81b44e04` (v6); Arm A and its scGPT `embed_data` stub
    were never run or tested.
13. `docs/COMPONENT2_PLAN.md`: the cosine-based design was superseded by the built metric.
14. `data/README.md` describes the Phase 1 inputs (128-latent SAEs, MD5 `60183a17...`).
15. PR #3's title and body describe only Gate-0 hardening; the branch now carries all of Component 2.
16. `LAB_NOTEBOOK.md` entries of 2026-07-05 and 2026-08-10 contain the control-only claim, the
    GATA1-in-Norman implication and the over-strong chance wording. They are left as written
    (history) and corrected by the 2026-10-18 entry.

---

## 18. Future directions

### 18.0 The decision in front of you
STS closes in 18 days. Two things compete for those days: one decisive experiment (about an hour
of GPU time) and the research report (most of the time). Do both, in parallel, in that order: launch
the run on day 1, write from day 1. About 80% of the report is identical whichever way the run
comes out (Gate 0, the four causal tests, the artifact catalog); only the centerpiece section
depends on it. Do not start new engineering (Arm A, new datasets) before November 5.

### 18.1 The next action: run program modulation (verified recipe)
Verified in this container on 2026-10-18, except the full-size GPU SAE fits: the download, the
size check and `prep_norman.py` were run on the real file; the loader and analysis ran at full
scale with stand-in activations; all three Perturb-seq CLIs ran end to end on a real-data
subsample with a tiny SAE; the pipeline self-test trained real SAEs on CPU. Every result JSON now
records `sae.trained_on` and `sae.n_train_cells`. Run on a rented NVIDIA GPU (PyTorch image); no
token on the box (if pip complains about `packaging`, add `--ignore-installed packaging`):

```bash
git clone https://github.com/jbar1125/SAE-scRNA-seq-Research.git && cd SAE-scRNA-seq-Research
git checkout claude/code-handoff-doc-2v3jt3
grep -q "def ensure_raw_counts" src/prep_norman.py && grep -q "def train_rows" src/causal_pipeline.py \
  && echo "CODE OK" || { echo "STALE - abort"; exit 1; }
python3 -m pip install -q scanpy anndata statsmodels
python3 tests/test_program_modulation.py | tail -1    # ALL PROGRAM-MODULATION TESTS PASSED
python3 tests/test_prep_norman.py | tail -1           # ALL PREP-NORMAN TESTS PASSED
python3 tests/test_pipeline_training.py | tail -1     # trains real SAEs; confirms torch works

# raw counts (NOT norman_2019(), which is log-normalized)
curl -sS --fail -o norman_raw.h5ad \
  https://scverse-exampledata.s3.eu-west-1.amazonaws.com/pertpy/norman_2019_raw.h5ad
[ "$(stat -c%s norman_raw.h5ad)" = "828160302" ] && echo "DATA OK" || echo "DATA SIZE MISMATCH"
python3 src/prep_norman.py --in norman_raw.h5ad --out norman_pert.h5ad   # expect: 102 testable

mkdir -p causal_out
for S in 0 1 2; do   # PRIMARY (spec v1.1): SAE fit on all cells
  python3 src/program_modulation.py --adata norman_pert.h5ad --pert-col perturbation \
    --control-value control --latent 2048 --k 32 --n-hvg 4000 --n-shuffle 50 \
    --train-on all --seed $S --out causal_out/progmod_seed$S.json
done
for S in 0 1 2; do   # PRE-DECLARED SENSITIVITY: control-only fit (unvalidated mode)
  python3 src/program_modulation.py --adata norman_pert.h5ad --pert-col perturbation \
    --control-value control --latent 2048 --k 32 --n-hvg 4000 --n-shuffle 50 \
    --train-on control --seed $S --out causal_out/progmod_ctrl_seed$S.json
done
```
Needs ~4 GB disk for the data (the prepped file is 3.0 GB). `pertpy.data.norman_2019_raw()`
fetches the same file if you prefer pertpy. On this CPU container the primary fits take 2 h per
seed, too long for its 2-hour background limit; on a GPU, minutes.

### 18.2 How to read it (pre-registered; decide before looking)
- **Positive:** an erythroid row (Ery_Effector or Ery_Heme most likely) with shuffle p < 0.05 AND
  background percentile >= 90 AND a coherent modulator list, in >= 2/3 seeds, while Gran,
  Progenitor and Cycling do not show the same. Then the project has a positive centerpiece: an
  unsupervised SAE recovers a canonical erythroid program whose activity is causally modulated by
  a specific perturbation set. If that set is CBL / UBASH3B / PTPN12, verify against the primary
  literature (cite it; do not assert it) that these KIT/EPOR-signaling brakes push K562 toward
  erythroid differentiation, and check where KLF1 (the one erythroid master TF in the library)
  lands.
- **Dominant-axis negative:** many modulators at an unremarkable background percentile. That is the
  K562 erythroid axis behaving like the stress axis did, and it is a negative. Given the Norman
  gate-saturation mechanism (section 9.5), this outcome is a live possibility; the background
  null exists to catch it.
- **Clean negative:** few modulators and shuffle p not small.
- **Mixed across seeds:** not established.
- The control-only sensitivity run is reported either way and is not a gate.

### 18.3 What to write for STS (three shapes, ranked by honest ceiling)
1. **Positive-centered** (if 18.2 is positive): "SAE gene programs in hematopoiesis are mostly
   method artifacts, and single-perturbation causal grounding fails across CRISPRi and CRISPRa,
   but a pre-registered program-level test recovers a causally modulated erythroid program."
2. **The auditor** (if negative): "A pre-registered audit of sparse-autoencoder interpretability in
   single-cell biology: structure is method-dependent; causal grounding fails at single-perturbation
   resolution in two screens and two cell lines; here is the catalog of how each naive metric
   misleads and the diagnostics that separate 'no signal' from 'test never fired'." This is the
   PLI-Analyzer shape. You said you would not submit a null; this is a methods-and-audit result
   with real positive statements (method-dependence, a seed-stable erythroid program, the
   dominant-axis mechanism, the toolkit), but it is negative about the original hypotheses.
   Whether to submit it to STS or wait for ISEF with more work is your call.
3. **The floor:** Gate 0 alone (frozen, cross-species, five methods).
Check the official requirements page for the report's limits and the other components.

### 18.4 After STS, toward ISEF (in priority order)
- **P0: Arm A, the expression-vs-embedding head-to-head on Replogle K562.** The original
  centerpiece, never run. Needs an scGPT (or Geneformer) embedding generator (not written; the
  `RUNPOD_EXECUTION.md` stub is untested) and a GPU; then `--rep embedding` with the same latent,
  k, and training mode as Arm B. Expression grounds ~13-20% of 162 TFs (all-cell SAE). If embeddings
  ground ~6% under the identical metric, "expression space is more causally legible" is a real
  comparative result, even though the programs are housekeeping. Run it as a matched grid.
- **P1: settle the training-mode question.** Either make control-only training work (diagnose its
  power loss; validate it on the synthetic world) or adopt all-cell training as the design and
  re-argue the triviality defense on its remaining legs (labels unseen, programs defined on
  controls, perturbed gene excluded, cross-fitting). Document whichever is chosen.
- **P1: if program modulation is positive:** dual-guide combinations (41,759 dual-guide cells),
  guide-level dose-response, and an independent dataset for replication.
- **P1: re-run `regulon_chance_baseline.py` on real v2 JSONs** (re-run v2 first; cheap on GPU).
- **P2:** the Lambert 2018 TF list; the Joung TF Atlas from a processed source if one is reachable;
  the AML clinical arm only if there is a positive program to translate.
- **P3:** causally-supervised SAE take 2; packaging the metric as a benchmark; a bioRxiv or
  workshop preprint (MLCB, LMRL).

### 18.5 Cheap hygiene before any writing
Fix the section 17 list (above all IMPACT_STRATEGY, DIFFERENTIATION, CLAUDE.md, README, the
gene-set README, panel B, PR #3); verify the citations the audit still marks VERIFY (A4:
Beneyto-Calabuig 2023, Lasry 2023, Petti 2019); re-read the Kendiukhov PDFs for their exact 48-TF
panel and grounding definition; record compute hours and cost.

### 18.6 What not to do
- No regulon recovery v3 (pre-committed) and no more patches to the binary grounding rate.
- No new engineering before November 5 unless the program-modulation result is in hand and the
  report is on track.
- Never feed pertpy's processed `norman_2019()` downstream.
- Do not cite: the 0.37 crossover, the 49 bimodal features, any Phase 1 participation ratio, the
  "~10x essential enrichment", "grounding tracks effect size", "the SAE was trained on control
  cells", or "GATA1 failed to ground in Norman".

---

## 19. How to resume (for you or a future agent)

1. Read section 0, then 17, then 18 of this file; then `CLAUDE.md`.
2. Environment (fresh container):
   `python3 -m pip install numpy scipy statsmodels pandas scanpy torch` and
   `python3 -m pip install --ignore-installed packaging anndata`.
3. Verify the repo:
   ```bash
   python3 -m py_compile src/*.py tests/*.py
   for t in tests/test_*.py; do python3 "$t" | tail -1; done      # 9 x "ALL ... PASSED"
   sha256sum config/*_spec.json                                     # compare to section 12
   cat data_g0/v0b_outputs/v0b_v3_1_decision.json data_g0/v0b_outputs/v0b_v3_1_per_seed.csv | sha256sum
   cat data_g0_human/v0b_outputs/v0b_v3_1_decision.json data_g0_human/v0b_outputs/v0b_v3_1_per_seed.csv | sha256sum
   ```
4. Follow CLAUDE.md: run before presenting, a verification table, SHA-pinned pastes with grep
   guards, lead with bad news, never freeze inconclusive output.
5. New results go in `LAB_NOTEBOOK.md` (dated), the relevant results doc, and section 10 here.

---

## 20. Glossary

- **SAE (sparse autoencoder):** learns an overcomplete dictionary so each cell is reconstructed from
  a few active features. **L1 SAE** (Gate 0) penalizes activation magnitude; **TopK SAE**
  (Component 2) keeps only the k largest activations per cell.
- **L0:** mean number of active features per cell; Gate 0's quality band is 20-50.
- **Decoder atom / program:** a feature's decoder column (Gate 0 loadings), or its correlation with
  each gene over control cells (Component 2 program).
- **Participation ratio:** (sum of eigenvalues)^2 / sum of squared eigenvalues, an effective
  dimensionality.
- **Per-seed winner-take-all:** assign features to submodules within each seed separately, because
  SAE feature identities do not persist across seeds.
- **BH-FDR:** Benjamini-Hochberg false-discovery-rate control within a test family.
- **Abundance-matched null:** random gene sets matched on decoder mass/expression, so "expressed
  genes attract decoder mass" cannot masquerade as a program.
- **Excess over baseline:** a submodule's strength minus the Progenitor control's (floor 0.10).
- **Mann-Whitney AUC:** U / (n1 x n2), the probability a perturbed cell's activation exceeds a
  control cell's; scale-free, robust to sparse activations.
- **Cross-fitting:** select on one half of a perturbation's cells, test on the other.
- **Robust z (median/MAD):** outlier score used by the match-confidence, column-specificity and
  effect-size gates.
- **Column specificity:** a perturbation must move a feature more than other perturbations do.
- **Effect-size control (M1):** movement must exceed what effect-size-matched peers show.
- **Label-shuffle null:** relabel cells and re-run. "Perturbed-only" keeps controls fixed (asks
  whether identity matters); "full" permutes everything (asks whether any group structure exists).
- **Feature-background null:** the same per-perturbation test on random features, to show a
  program is special rather than a dominant axis.
- **Gate saturation:** when the absolute gates pass nearly everything, a relative gate sets the
  rate and a shuffle reproduces it.
- **Dominant axis:** a direction of variation (in K562: erythroid, stress) large enough that most
  groupings of cells differ from controls along it.
- **TopK competition artifact:** an activated program pushes competitors out of the top k, so they
  look suppressed.
- **CRISPRi / CRISPRa / Perturb-seq:** CRISPR knockdown / activation with single-cell RNA readout.
- **Essential library:** Replogle's screen of essential genes; few lineage TFs.
- **TRRUST, regulon, hub target:** curated TF -> target database; a TF's target set; a target shared
  by many TFs (CDKN1A: 115).
- **HVG:** highly variable gene. **K562:** erythroleukemia line. **RPE1:** retinal epithelial line.
- **Pre-registration / amendment / erratum:** design frozen by SHA before running; a dated change
  made before a real run; a dated correction of a factual misstatement.

---

*Sources for the STS dates: [Regeneron STS application requirements](https://www.societyforscience.org/regeneron-sts/application-requirements), [Regeneron STS FAQ](https://societyforscience.org/regeneron-sts/frequently-asked-questions).*
