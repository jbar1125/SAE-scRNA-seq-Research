# STRATEGY & POSITIONING: publication / STS / ISEF deep review

> **CURRENT STATUS, 2026-10-18:** this is a historical strategy review. Its staged
> gates, competition predictions and unexecuted components are proposals, not a
> record of completion. Gate 0 is complete and frozen on Paul15 mouse and Setty 2019
> human CD34+ marrow (4,142 cells), with five method families. Its durable finding
> is method-dependence; the reversed lineage direction is not invariant or causally
> validated. Replogle K562/RPE1, Norman and regulon analyses ran later without
> establishing the intended lineage-regulatory claim. The causal ladder used
> PCA/NMF **as inputs to an SAE**, not direct factor baselines. Foundation-model
> embeddings, AML and tumor translation never ran.
>
> Every historical causal SAE fit used all cells, despite control-only design
> language; causal raw JSONs from rented compute are missing. Recover/regenerate
> these before publication claims depend on their exact reported numbers. The
> wrong-background essential enrichment and failed specificity gate are not solved
> confounds, and M1 did not support the effect-size explanation. Current evidence
> and priorities: `PROJECT_MASTER_SUMMARY.md` sections 8-10 and 18,
> `IMPACT_STRATEGY.md`, and `DIFFERENTIATION.md`. The remaining historical plan below
> does not establish pre-registration of later unexecuted analyses.

Brutally honest external-reviewer analysis of where this project stands in the
field, what will win and what will sink it, and the concrete reframe + revised
plan. No sugarcoating (per CLAUDE.md). Read with PROJECT_AUDIT.md.

Bottom line up front: the original novelty claim is largely gone, but the project
has *organically become* the kind of rigorous "does the hyped method actually
recover real biology, and how would you know?" project that has won STS top-10
slots. The path to competitive is to embrace that, lead with the
artifact-vs-signal framework + a causal test, and stop competing on "we applied
SAEs to cells first." It cannot be made "impossible to lose" — that is not a real
state for STS/ISEF. It can be made genuinely top-tier. Below is how.

---

## 1. Field review recorded in 2026 (historical)

SAEs on single-cell data is an active, published area, not white space:

- Kendiukhov 2026, "Sparse autoencoders reveal organized biological knowledge but
  minimal regulatory logic in single-cell foundation models: a comparative atlas
  of Geneformer and scGPT" (arXiv 2603.02952). This is the project's own 6.2%-null
  source. Its thesis (features encode real biology but weak causal/regulatory
  logic) is adjacent to this project's emerging conclusion. TREAT AS THE PRIMARY
  COMPETITOR, cite explicitly, and differentiate.
- "Sparse Autoencoders Reveal Interpretable Features in Single-Cell Foundation
  Models" (bioRxiv 2025.10.22.681631).
- "Discovering Interpretable Biological Concepts in Single-cell RNA-seq Foundation
  Models" (arXiv 2510.25807).
- "Can sparse autoencoders make sense of gene expression latent variable models?"
  (arXiv 2410.11468) — closest to expression/LVM space rather than embeddings.
- Foundations: Cunningham/Bricken 2023 (SAEs find interpretable LM features,
  arXiv 2309.08600); Sparsely-connected AE for scRNA-seq (npj Syst Biol Appl 2020).
- Field leadership: Fabian Theis (2025 ISCB Innovator, single-cell + ML), James
  Zou (2025 ISCB Overton, ML for biology). Judges in this niche may know this work.

Implication: "apply SAEs to single-cell data to find interpretable features" is
established (2024-2026). But — CORRECTION to an earlier overstatement in this doc:
the reviewed Kendiukhov cluster is on FOUNDATION-MODEL EMBEDDING/activation space
at the atlas level. Kendiukhov trains
SAEs on Geneformer/scGPT residual streams (~82k / ~24k features), not on raw
expression, and its 6.2% is the reported rate OF FOUNDATION-MODEL features under its
own metric and panel. This project has not executed the proposed integrated
expression-space, developmental, matched-causal and clinical program, nor established
that the entire literature lacks it. The base technique is not novel. The earlier
"the project is done" framing was too strong and is retracted.

Two adjacent papers to personally verify (I could not fully read them):
- CytoSAE: Interpretable Cell Embeddings for Hematology (arXiv 2507.12464).
- Discovery of a Hematopoietic Manifold in scGPT (arXiv 2603.10261).
Both are hematology + SAE/foundation-model; check whether either does
expression-space + causal + clinical before claiming that ground.

## 2. Novelty audit — what dies, what survives

DEAD as novelty (do not claim these):
- "SAEs recover interpretable features from single-cell data" as a general idea.
- "Embedding-space SAEs on Geneformer/scGPT" (Phase 2 Component 4). This IS
  Kendiukhov's atlas. Frame Component 4 as a controlled comparison/replication, not
  a novel contribution.

Historical candidate differentiators (current completion/scope is stated above):
1. **The artifact-vs-signal framework.** The v1->v2->v3->v3.1 saga is not
   embarrassing history; it is the contribution. You have a documented, executed,
   pre-registered method for deciding whether an SAE "program" is real or a metric
   artifact (control-referenced strength, count-fair concentration, calibrated
   null, marker-curation audit). No SAE-single-cell paper does this rigorously.
   This is the novel, defensible core. Most of the field reports features and
   asserts interpretability; you built the adversarial test.
2. **The head-to-head causal question (the sharpest novel claim).** Kendiukhov
   showed FOUNDATION-MODEL SAE features have minimal causal logic (6.2% of TFs).
   The open, unpublished question: do EXPRESSION-SPACE SAE features, annotated along
   a developmental trajectory, have BETTER causal grounding when tested against the
   same Replogle CRISPRi data? This was a proposed controlled head-to-head against
   a literature reference. A valid comparison requires a within-study embedding run
   under the same metric, panel and settings; exceeding a percentage from another
   paper would not by itself locate causal structure. This head-to-head remains a
   proposal.
3. **Expression-space + developmental trajectory.** Most of the field does
   embedding-space, atlas-level "what features exist." You do direct expression
   SAEs along a commitment trajectory (pseudotime + branch probabilities +
   submodule structure). More mechanistic, less crowded.
4. **A decomposition-dependent structural observation:** the SAE gives fewer
   erythroid than granulocyte programs in both executed species. The five-method
   comparison does not establish an invariant biological modularity direction.
5. **Expression<->embedding concordance** (does a program found in raw expression
   survive re-derivation in scGPT space?). This is a real cross-representation
   question the atlas papers do not answer.

## 3. Benchmark vs recent STS / ISEF winners

Recent top computational-bio work (Society for Science, 2024-2026):
- Frances Liang, "PLI-Analyzer" (2026 STS top 10): a tool to test whether
  AI-predicted protein complexes are biologically plausible; found AlphaFold3 and
  Boltz-2 ~50% accurate. THIS IS THE TEMPLATE. It won by rigorously auditing
  whether a hyped AI method actually works. Your project is the SAE analogue.
- Mythreya Dharani, "I2-CISSP" (2026 STS top 10): ML chemo-response prediction with
  interpretability of which genes drive the call. Wins on interpretability + clinical
  framing.
- Jerry Xu (2026 STS top 10): efficient protein-structure comparison via learned
  compression.
- Sophie Chen (2024 finalist): CNN for squamous-cell-carcinoma histology.
- Finan Gammell (2025 finalist): statistical model for cancer-progression genes.
- ISEF 2025 comp-bio: StAIR (atom-level protein reconstruction), deep-learning drug
  discovery on olfactory receptors, ALS biomarker + therapeutic.

What wins (pattern): addresses a REAL GAP in existing tools; combines ML with
biological insight; delivers INTERPRETABILITY or VALIDATION; report reads like a
graduate thesis; a crisp, defensible headline a non-specialist judge can repeat.

Where this project currently stands against that bar:
- Rigor/validation: potentially ABOVE the median winner (few high-schoolers
  pre-register, audit their own artifacts, and run causal tests). Strength.
- Novelty of method: BELOW (the method is now published). Must pivot the novelty
  to the framework + causal test, not the SAE application.
- Status at this review's drafting: mouse-only, n=5, causal arm unbuilt. Update:
  Gate 0 subsequently ran 10 seeds in both species and causal expression analyses
  ran. Their limitations are stated in the current-status notice above.
- Communicability: currently a tangle of v-numbers. Needs a single clean narrative.

## 4. Vulnerabilities (the full gap list, ranked by how badly a judge/reviewer hits it)

CRITICAL (fix or the project is beatable on contact):
- V-A. Novelty collision with Kendiukhov et al. Not addressed anywhere in current
  materials. A judge who knows it will ask "how is this not their atlas?" You need
  a one-paragraph, rehearsed answer and an explicit differentiation table.
- V-B. Metric-dependent verdict. The headline read NO/YES/reversed across metrics.
  Unmanaged, this looks like p-hacking. Managed (one pre-registered metric, the
  whole point being that naive metrics disagree and here is the principled one), it
  is a strength. You MUST pre-register v3.1 and freeze before Phase 2.
- V-C. The causal arm was unbuilt at drafting. It subsequently ran in expression
  space; the matched foundation-model comparison remains unexecuted.
- V-D. Marker-set dependence. The globin bug proved the entire conclusion can hinge
  on marker curation. Reviewers will suspect cherry-picking. Needs: pre-registered,
  literature-cited, non-overlapping marker sets; a leave-one-marker-out sensitivity
  analysis; and orthogonal (non-marker) validation (e.g., GO/pathway enrichment on
  decoder genes, blinded).

HIGH:
- V-E. Statistical power: n=5 seeds, n=5 bootstrap, mouse-only. Need >=10 seeds,
  permutation nulls, human replication, effect sizes with CIs everywhere.
- V-F. SAE not sparse in the upgraded run (L0~63 > 20-50). A non-sparse "sparse"
  autoencoder is a self-inflicted wound a reviewer will flag instantly. Tune L1 to
  the QC band and report an L0/L1 sweep.
- V-G. Baseline comparison (NMF/PCA) incomplete. "Does the SAE structure beat matrix
  factorization?" is the first question any ML reviewer asks. V11/V12 must be
  finished and the SAE must be shown to add something PCA/NMF do not.
- V-H. Compressive vs overcomplete framing. Phase 1 used 128<2000; executed Gate 0
  used 512<2012 (mouse) and 512<2032 (human). These are undercomplete relative to
  input genes. Code supporting larger dictionaries does not make those runs
  overcomplete. Justify the actual dictionary and report declared sweeps.

MEDIUM:
- V-I. Provenance/reproducibility gaps: the MD5-locked matrix's true preprocessing
  is uncertain (the "pre-log-transformed" claim was false); CELLxGENE census date
  unpinned; exact Phase-1 training config unpinned. Fix before any submission.
- V-J. Citation integrity: Kendiukhov year vs arXiv-id mismatch; unverified
  Beneyto-Calabuig/Lasry/Petti. A single wrong citation dents credibility with a
  thesis-level judge.
- V-K. Overclaiming risk: "asymmetric modularity" language implies mechanism the
  data does not yet support. Discipline the claims to what survived v3.1.
- V-L. ISEF/STS logistics: fully computational + public data is fine, but the
  Research Plan must nail the "why this matters" and human-relevance framing for a
  multidisciplinary panel (most judges are not ML people).

## 5. THE REFRAME (new thesis and narrative)

Old (dead) thesis: "SAEs recover interpretable, asymmetrically-modular gene
programs in hematopoiesis."

New thesis (defensible, novel, winning-shaped):

> "Interpretability methods are being applied to single-cell data faster than they
> are being validated. I build an adversarial framework that decides whether an
> SAE-derived gene program is a real biological module or a measurement artifact —
> using control-referenced effect sizes, calibrated nulls, and perturbation tests
> whose evolving design and execution boundaries are documented — and apply it to hematopoietic lineage
> commitment. The framework overturns my own initial finding (showing it was a
> marker-curation and metric artifact), shows that lineage modularity depends on
> the decomposition, and tests the separate hypothesis that annotated expression
> features identify perturbation-specific regulatory programs. The causal tests
> have not validated the structural reversal or established the intended
> lineage-regulatory claim."

Why this wins a multidisciplinary panel:
- The hook is understandable to any judge: "everyone's using AI interpretability on
  biology; does it actually work? I built the test." (Same shape as PLI-Analyzer.)
- It turns your biggest weakness (the flip-flopping finding) into the evidence for
  your contribution (naive metrics disagree; here is the principled one).
- It has a causal spine (CRISPRi) that most interpretability work lacks.
- It ends on bounded structural evidence and explicitly reported causal-test limits.
- The self-correction (you falsified your own headline) reads as scientific
  maturity, which top judges reward.

Narrative arc for the manuscript / presentation:
1. Motivation: interpretability-for-biology is outpacing validation (cite the
   2025-2026 SAE-single-cell wave). Pose the question: real modules or artifacts?
2. System: SAEs on hematopoiesis expression + trajectory (why this is a good
   testbed: known biology, sharp lineage structure, CRISPRi available).
3. The trap: naive per-feature enrichment gives a headline that flips with metric
   choice and marker curation (show v1->v2->v3). This is the problem.
4. The framework: control-referenced strength, count-fair concentration, calibrated
   permutation null, marker-curation audit, pre-registration. (The contribution.)
5. Structural result: modularity answers differ across five decomposition families
   in mouse and human. Do not claim the SAE beats PCA/NMF at biological truth.
6. Causal tests: report Replogle/Norman and regulon outcomes with the actual metric,
   training mode and archive limitations. The 6.2% literature rate is context;
   an internally matched embedding comparator remains future work.
7. Honest scope: what survives, what does not, what the method says about the
   broader SAE-single-cell literature.

## 6. Revised research plan (concrete, phased, to bulletproof)

GATE 0 — Lock the foundation (weeks 0-2). Non-negotiable before anything else.
- Pre-register ONE metric: v3.1 control-referenced, Progenitor-only baseline (drop
  Cycling — cell cycle is a real program, a bad null), abundance-matched
  permutation null, effect-size floor. SHA-freeze the decision spec BEFORE re-runs.
- Fix sparsity: L1 sweep to L0 in 20-50; report the sweep. Retrain >=10 seeds.
- Finish NMF/PCA baselines (V11/V12) on identical inputs; show what SAE adds.
- Marker sets: pre-register from cited literature, non-overlapping, with a
  leave-one-out sensitivity analysis and an orthogonal blinded GO/pathway check.
- Reproducibility: pin census date, exact training config, resolve the
  log-transform provenance, fix citations.

GATE 1 — Lock the finding (weeks 2-5).
- Human replication was subsequently executed on Setty 2019 CD34+ marrow, 4,142
  cells, not the older 25k CELLxGENE SAEs. The direction remains method-dependent.
- Robustness: dictionary-size sweep (256/512/1024/2048), TOP_K sweep, seed
  stability of the verdict. Report the verdict's sensitivity, not one number.
- Freeze the module assignment table (SHA) as the causal pre-registration.

GATE 2 — The causal centerpiece (weeks 5-9). This is the differentiator; do not cut.
- Replogle K562 CRISPRi: project SAEs onto NT + each TF knockdown, standardized
  effect, top-k suppression, one-sided binomial vs 6.2%. OOD gate active (fix V14).
- Report per-TF causal validation table; frame honestly against Kendiukhov's 6.2%.
- If it matches the null, that is STILL a publishable, honest result (and it is the
  same shape as PLI-Analyzer's "these tools are only 50% right").

GATE 3 — Situate + write (weeks 9-12).
- Explicit differentiation section vs Kendiukhov et al. and the 2025-2026 SAE
  wave (a table: dataset, space, trajectory, causal test, validation rigor).
- Manuscript to the arc in section 5. Target a real venue (an ML-for-genomics
  workshop or bioRxiv) BEFORE STS/ISEF — a preprint massively strengthens both.

## 7. What to cut, what to add

CUT / de-emphasize:
- Any "first to apply SAEs to single-cell" claim.
- The embedding-SAE atlas as a novelty (keep only as a concordance check).
- The Shannon-entropy-ratio and largest-fraction-alone metrics (superseded; keep
  only as documented negative controls showing why naive metrics fail).
- Scope sprawl (Park thymus Part 5, extra datasets) until the core is bulletproof.

ADD:
- The causal CRISPRi arm as the centerpiece (currently underweighted).
- Human replication as a hard gate.
- The differentiation-vs-prior-work section.
- A "limitations and what does not survive" section (judges trust projects that
  state their own limits).
- Orthogonal, non-marker validation of features (blinded GO/pathway).

## 8. Honest odds and non-negotiables

Odds, stated plainly: with GATE 0-2 completed (pre-registered metric, human
replication, causal test done, baselines beaten, differentiation clear), this is a
credible STS finalist / ISEF top-award-caliber project in Computational Biology and
Bioinformatics, and a plausible workshop preprint. Without the causal arm and human
replication, it is a late replication of published work and will not stand out.
Nothing here guarantees a win; anyone who tells you a project is "impossible to
lose" is selling something. What is real: this can be made genuinely excellent, and
excellence + a rehearsed, honest narrative is what maximizes the odds.

Non-negotiables (if only these get done):
1. Pre-register v3.1 (one metric), freeze it, then never move it.
2. Human replication of the corrected finding.
3. The CRISPRi causal test against the 6.2% null.
4. Explicit, unhidden differentiation from Kendiukhov et al.
5. Discipline every claim down to what survived the framework.
