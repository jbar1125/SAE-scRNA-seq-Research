# UPGRADE BACKLOG — exhaustive list of ways to elevate this project

Purpose: a single, prioritized, brutally-honest backlog of every upgrade and
redirection worth considering, big or small. Written 2026-07-18 against the real
state of the repo (Gate 0 frozen; the causal head-to-head metric just fixed and
about to produce its first real number). Cross-references the gap list in
`STRATEGY_AND_POSITIONING.md` (V-A..V-L) and the plan in `ELEVATION_PLAN.md`.

How to read this: tiers are ordered by LEVERAGE (payoff per unit effort), not by
size. Tier 0 is cheap and near-mandatory. Tiers 2-3 are the big swings the user
asked for. Each item says what it is, why it helps, and rough effort. "DONE"/
"PARTIAL" tags mark what already exists so nothing here is double-counted.

Honest frame: the project's one winnable differentiator is the CAUSAL head-to-head
(is regulatory logic legible in expression but compressed away in foundation-model
embeddings?). Almost everything below is in service of making that result
unassailable, more general, or more clearly a built thing rather than a number.

---

## SURFACED FROM THE FIRST REAL RUN (2026-07-18) — do these before any specificity claim

The first real Arm-B run grounded 22/162 = 13.6% (survives the shuffle null), BUT the
grounded set is dominated by essential / general-transcription-machinery / non-TF
knockdowns (TBP, TAF7, GTF2A2, POLD2, RFC2, SRP9, CSNK2B, PTPMT1, ...) with NO lineage
regulators (GATA1/SPI1/TAL1/RUNX1). That is essential-gene stress collapse, not regulatory
specificity (LAB_NOTEBOOK 2026-07-18). These three fixes are now REQUIRED, not optional:

S.1 **Stricter TF list.** Replace `config/tf_lists/hs_hgnc_tfs.txt` (contains basal
machinery) with the Lambert et al. 2018 sequence-specific-DBD human TF list (~1600). Commit
it to `config/tf_lists/lambert2018_hs_tfs.txt` and re-run. Expected to drop most artifacts.

S.2 **A specificity filter that actually bites.** In `causal_grounding.py`: add a
program-coherence requirement (matched feature's top decoder genes form a GO-enriched /
non-diffuse module) and tighten the match gate — right now it passes 162/162, so it does no
filtering. The essential-gene "global lowness" features should be filtered here.

S.3 **Essential-gene control arm.** Score grounding separately for essential vs
non-essential genes using the DepMap common-essential list (commit it to `config/`). If the
rate is high only among essentials, report THAT (a different, still-honest finding).

---

## TIER 0 — cheap, critical, do before the next writeup (mostly no GPU)

0.1 **Differentiation table vs Kendiukhov et al.** (V-A, flagged CRITICAL, still not
built). A one-screen table: rows = this project vs Kendiukhov 2026 vs the 2025-26
SAE-single-cell wave; columns = dataset, representation space (expression vs
foundation-model embedding), trajectory-resolved?, causal test?, rigor devices. A
judge who knows the paper WILL ask "how is this not their atlas?" You need this
rehearsed and on paper. Effort: hours.

0.2 **Verify the two adjacent hematology+SAE papers.** CytoSAE (arXiv 2507.12464)
and "Hematopoietic Manifold in scGPT" (arXiv 2603.10261). The strategy doc flags
both as UNVERIFIED potential collisions. If either already does expression-space +
causal grounding, you must know before you write. Effort: an afternoon of reading.

0.3 **Threshold sensitivity sweep on the causal metric.** Report grounding rate as a
function of `auc_floor` (0.45), `match_alpha` (0.10), `min_active_cells` (50),
`latent` (2048), `k` (32). Gate 0 already did this for its decision params (R1); the
causal arm MUST too, or "you tuned 0.45 to get the answer" is a fair hit. Pre-register
the primary values, report the surface. Effort: a sweep script + a heatmap.

0.4 **Positive + negative control perturbations on REAL data (not just synthetic).**
Positive: pick TFs with textbook K562 targets (GATA1 -> globin/heme; this ties K562
back to your erythroid story) and confirm the KD suppresses the hemoglobin-program
feature by hand. Negative: confirm non-targeting / safe-harbor / intergenic guides do
NOT get grounded. A metric that self-evidently passes known positives and rejects
known negatives ON REAL DATA is far more convincing than a synthetic oracle alone.
Effort: low once the run exists.

0.5 **GO / pathway enrichment on every grounded program (orthogonal, blinded).** V-D
asks for non-marker validation. A grounded erythroid program should enrich heme
biosynthesis; a granulocyte one should enrich azurophil-granule genes. Blinded GO on
decoder-top genes is independent evidence the programs are biologically real, not
metric artifacts. Effort: low (gseapy/enrichr on decoder genes).

---

## TIER 1 — make the causal centerpiece unassailable (the live differentiator)

1.1 **Run Arm A (embedding-space) YOURSELF; stop leaning on the cited 6.2%.** The
single highest-value technical upgrade. Embed the SAME Replogle control cells with
scGPT, train the identical TopK SAE on the embeddings, run the identical grounding
metric. Now the head-to-head is apples-to-apples and a judge cannot say "different
metric." COMPONENT2_PLAN calls this out as the honesty gap; close it. Effort: the
scGPT inference step (GPU, hours) is the one real lift.

1.2 **Two foundation models, not one.** Do scGPT AND Geneformer for Arm A (Kendiukhov
did both). If expression beats both, "it's the embedding compression" is airtight. If
it beats one and not the other, that split is itself an interesting finding. Effort:
+1 embedding run.

1.3 **SAE-vs-PCA-vs-NMF-vs-raw causal control (V-G, the first ML-reviewer question).**
Run the grounding test on PCA components and NMF factors of the SAME expression, and
on raw HVG expression. Answers: is the causal legibility about the SAE, or just about
staying in expression space? If expression-SAE >> expression-PCA, the SAE earns its
place. If PCA ties it, the honest story becomes "expression space, not the SAE" — still
publishable, but you must KNOW. Effort: moderate (reuse baselines_nmf_pca.py).

1.4 **RPE1 as a second cell line.** Replogle is genome-scale in K562 AND RPE1. Re-run
the whole head-to-head on RPE1. Cross-cell-line replication is the causal analogue of
the mouse+human replication you already did in Gate 0 — a robustness win that costs
only compute. Effort: repeat the pipeline on a second h5ad.

1.5 **Multi-seed with the >=3/5 rule + bootstrap CI on the DIFFERENCE.** Require a
perturbation to ground in >=3/5 seeds (your standing convention). Put a bootstrap CI on
the expression-minus-embedding grounding-rate difference; the claim holds only if that
CI excludes 0 (COMPONENT2_PLAN already says this — enforce it). Effort: seeds 0-4 loop
(already scripted) + a CI script.

1.6 **The per-TF causal-certificate table = the deliverable.** Don't report only a
rate. Emit a table: TF, matched feature, program genes, AUC, q, seed-support. This IS
the artifact (ELEVATION_PLAN §4) — a resource, not a number. Effort: format the
existing `per_perturbation` output.

1.7 **Cross-fit K-repeat stability.** You split A|B once (seeded). Repeat the split K
times and average, so no result rides on one lucky split. Closes a real objection
cheaply. Effort: a loop around the existing split.

1.8 **Dose-response.** Replogle carries per-guide knockdown efficiency. Do features
whose regulator is knocked down harder show stronger suppression? A monotone
dose-response is a strong causal signature that rules out coincidence. Effort: moderate.

---

## TIER 2 — the flagship redirection (bigger swing, highest upside)

2.1 **The representation-ladder causal-decay CURVE (the money figure).** Generalize the
binary expression-vs-embedding into a LADDER: raw expression -> HVG -> PCA -> NMF -> VAE
latent -> scGPT -> Geneformer. Plot causal-grounding rate along the ladder. If it decays
monotonically as representations abstract away from expression, the headline stops being
"we beat 6.2%" and becomes a LAW: "causal regulatory legibility decays as learned
representations abstract away from gene expression." That is a figure a judge remembers
and a genuinely novel scientific statement, and it makes ANY monotone outcome a result
(de-risks the "what if expression isn't much higher" fear). Highest-upside single idea
in this document. Effort: it is mostly Tier-1 arms plus 1-2 more representations, framed
as a curve.

2.2 **Held-out perturbation PREDICTION test.** The strongest possible causal claim.
Hold out a set of TFs entirely; ask whether the SAE's programs PREDICT their unseen
knockdown signatures (feature -> decoder direction -> predicted DE, scored vs the real
held-out DE). If expression-space programs predict held-out perturbations better than
embedding-space, that is prediction, not correlation. Effort: moderate; reuses the
matching machinery on a held-out split.

2.3 **Clinical payload — map grounded programs onto leukemia patient data (V-L human
relevance).** K562 is CML. Take AML/MDS patient scRNA-seq (van Galen 2019,
Beneyto-Calabuig 2023 — already in the audit's reference list) and ask: are the
CAUSALLY-GROUNDED programs the ones dysregulated in leukemia? A program that is both
causally grounded AND a leukemia axis = causality + biology + clinical, the trifecta
that turns "nice method" into "STS-finalist biology." Effort: moderate-high (a second
dataset + projection), but it is the biggest credibility jump for a multidisciplinary
panel.

---

## TIER 3 — method contribution (high risk / high reward)

3.1 **Causally-supervised SAE, take 2.** The first attempt failed and was honestly
shelved (ELEVATION_PLAN §4, LAB_NOTEBOOK 2026-07-05). But a working version is the
"I built a new method" headline. New formulations to try, each with a held-out-perturbation
eval so it cannot fool itself: (a) contrastive — pull features toward perturbation-response
directions using an auxiliary loss on one held-out perturbation set, evaluate on a
DIFFERENT held-out set; (b) a perturbation-consistency regularizer penalizing features
whose program never moves under any perturbation; (c) a decoder-Jacobian alignment
penalty. If ANY variant beats plain SAE on held-out grounding, that is a real method
result. If none do, report it (you already have the honesty muscle for this). Effort:
high; genuinely uncertain payoff.

3.2 **Package the grounding test as a reusable benchmark ("CausalGround" or similar).**
Turn `causal_grounding.py` + `causal_pipeline.py` into a standalone, documented benchmark
anyone can run on any representation of any Perturb-seq dataset. Benchmarks are high-impact,
citable, and judge-legible ("she built the causal-validation benchmark the field was
missing"). This is the "auditor tool" deliverable in package form. Effort: moderate
(API polish, README, examples).

---

## TIER 4 — artifact / resource (low-risk wins once the number lands)

4.1 **Interactive causally-grounded myeloid program atlas (web).** You already have
`docs/web/`. Build a browsable page: each grounded program -> its genes, GO enrichment,
causal certificate (the perturbation that validates it), and activity along pseudotime.
This is the "built resource" that wins (ELEVATION_PLAN §4), and it is low-risk formatting
once the head-to-head output exists. Effort: moderate, mostly front-end.

4.2 **Program cards (named, human-readable).** "Program 17 = heme biosynthesis, grounded
by GATA1 knockdown (AUC 0.31, q 1e-4), active at erythroid commitment." A judge can read
these; a rate table they cannot. Effort: low.

4.3 **Druggability / actionability overlay.** Map grounded programs' regulators to DGIdb.
"Causally-validated, druggable gene programs in myeloid cells" adds a translational hook.
Effort: low.

---

## TIER 5 — rigor, reproducibility, engineering

5.1 **Random-dictionary + SCENIC/TRRUST orthogonal baselines on the CAUSAL arm.** Random
decoder directions should ground near the FDR floor (confirms the test isn't structurally
inflated). And you already built `grn_baseline.py` + `trrust_crossref.py` for Gate 0 — ask
whether your grounded TF->program pairs overlap SCENIC/TRRUST regulons above chance
(orthogonal validation from an established tool). Effort: low-moderate (code exists).

5.2 **Power analysis.** How many perturbations are needed to detect a given expression-vs-
embedding gap at your FDR? Almost no high-schooler reports one; it reads as graduate-level.
Effort: low.

5.3 **FDR calibration on real data.** Show the label-shuffle grounding rate is <= your
nominal 0.05. You run a shuffle null; make the calibration explicit as a reported number.
Effort: trivial (already computed).

5.4 **Effect sizes with CIs everywhere (V-E).** AUC / Cohen's d with confidence intervals,
not bare p-values, throughout the causal arm. Effort: low.

5.5 **Pin the causal-arm environment + one-command reproduction.** Freeze pertpy version,
Replogle Figshare DOI, scGPT checkpoint hash, CUDA, seeds; add a `make causal` that runs
the whole head-to-head end to end (as REPRODUCIBILITY.md does for Gate 0). A pinned
RunPod/Docker template so a mentor or judge can re-run. Reproducibility is already a
project strength — extend it to the differentiator. Effort: moderate.

5.6 **Log total compute + dollar cost.** You are running on rented GPU. A methods line
"total compute: X GPU-hours, $Y, fully public data" is honest, rare, and judge-friendly.
Effort: trivial.

5.7 **CI: run the logic tests on every push (GitHub Actions).** A green "tests pass" badge
is cheap credibility and guards against the exact "shipped untested code" failure CLAUDE.md
was written around. Effort: low.

---

## TIER 6 — narrative, positioning, competition

6.1 **Lead with the "AI-interpretability auditor" framing (STRATEGY §3, the PLI-Analyzer
template).** A recent STS top-10 won by rigorously auditing whether a hyped AI method
actually works. Your one-line hook: "everyone is applying AI interpretability to biology
faster than they validate it; I built the causal test — and found where regulatory logic
actually lives." Rehearse it. Effort: writing.

6.2 **Kill the v-number tangle (communicability gap).** One clean narrative in the
manuscript; v1->v2->v3.1 history goes to a version-history appendix (VERSION_HISTORY.md
already exists — point the reader there and keep the body clean). Effort: writing.

6.3 **Preprint BEFORE STS/ISEF (Gate 3).** A citable bioRxiv / ML-for-genomics workshop
paper (MLCB, LMRL, a NeurIPS workshop) massively strengthens a competition entry. Effort:
the writeup, once the number lands.

6.4 **A named "Limitations / what does not survive" section.** K562-not-primary, one
matching operationalization, the shelved supervised-SAE. Judges trust projects that state
their own limits. Effort: writing (you already have the honest content in the docs).

6.5 **Rehearsed one-paragraph novelty answer** to "how is this not Kendiukhov's atlas?"
Memorized, crisp. Effort: trivial once 0.1 exists.

---

## REDIRECTIONS — the big reframes (choose one to organize around)

R1 **"The causal legibility of learned representations" as the unifying thesis.** Fold
Gate 0 (structure is method-dependent) and Component 2 (causal legibility decays with
abstraction, per the ladder 2.1) into ONE arc: naive structural claims are method-
dependent AND causal legibility is highest in expression and lost under abstraction —
therefore validate causally, in expression space. This is a stronger, more general thesis
than the current two-act framing and it makes Gate 0 setup rather than a detour.

R2 **Dual deliverable: atlas (resource) + CausalGround (tool).** Reframe the output from
"a number" to two built things — the causally-grounded myeloid program atlas (4.1) and the
open grounding benchmark (3.2). "I built the test AND the resource it produces" is a
maturity signal.

R3 **Go up a level to a general benchmark.** Position CausalGround as the standard any
representation-learning method for single cell should be scored on. Highest ceiling,
highest effort; only pursue if the core result is solid.

---

## Honest risk register

- **Expression may not beat embedding by much.** Mitigations already designed: the ladder
  (2.1) makes any monotone trend a result; the PCA/NMF control (1.3) localizes where the
  signal is regardless; internal head-to-head means you never depend on matching the cited
  6.2%. If expression ties embedding, that corroborates Kendiukhov from a new angle — still
  publishable, reported honestly.
- **The supervised SAE (3.1) may fail again.** It already did once. Treat it as optional
  upside, never as a load-bearing claim.
- **Clinical mapping (2.3) and the ladder (2.4/2.1) add scope.** Do them only after the
  core K562 head-to-head is locked; do not let scope sprawl delay the one result that
  matters (STRATEGY §7 warns about exactly this).
- **Everything here is subordinate to one thing:** get the first honest, multi-seed,
  Arm-A-included head-to-head number. Most of this backlog only matters if that lands.
