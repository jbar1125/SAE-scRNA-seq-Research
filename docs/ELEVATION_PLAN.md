# ELEVATION PLAN — turning a rigorous null into a winning positive result

Written after the user's honest call: "It is not a strong result as of now and will
not win anything." That is correct. This document diagnoses why, then designs the
transformation, grounded in the current (2025-2026) literature. Brutal honesty
throughout, per CLAUDE.md.

## 1. Why the current project does not win (diagnosis)

- **It is a negative result about a niche claim.** "Asymmetric modularity is not
  supported" answers a question few judges have heard of, with a "no."
- **"Methods disagree" is hygiene, not discovery.** That PCA/NMF/ICA/SAE give
  different decompositions is expected; demonstrating it rigorously is good practice
  but not a finding a judge remembers.
- **Zero causal or biological payload.** Everything is structural (decoder-mass
  concentration). No new biology, no prediction, no built artifact.
- **The framing is defensive** ("prior claims may be artifacts") rather than
  generative ("we found / built X"). Grand-prize projects are generative.
- **The substrate is crowded.** "SAEs on single-cell" is now a wave (Kendiukhov 2026;
  the Geneformer/scGPT SAE atlas, arXiv 2603.02952; CytoSAE; biorxiv 2025.10.22).
  Being one more entrant is not a differentiator.

Verdict: the rigor is real and worth keeping, but it is a *foundation*, not a result.

## 2. The reframe — a positive, benchmark-beating central question

The most important fact in the recent literature: SAEs on single-cell **foundation
models** (Geneformer, scGPT) encode organized biological knowledge but **minimal
causal regulatory logic** — only **3 of 48 TFs (6.2%)**, rising to just **10.4%** with
a multi-tissue control, show regulatory-target-specific feature responses when tested
against genome-scale CRISPRi (Replogle) perturbations (arXiv 2603.02952).

Nobody has asked the obvious next question, and it is the winnable one:

> **Is that 6.2% causal ceiling a limitation of SAEs, or of foundation-model
> embeddings? Do SAEs trained directly on gene EXPRESSION recover the causal
> regulatory logic that foundation models discard?**

**Hypothesis (positive, falsifiable):** Expression-space SAE features are causally
grounded at a substantially higher rate than 6.2-10%, because raw expression retains
the direct co-perturbation structure that foundation-model compression throws away.

If true, the headline is: **"Expression-space sparse autoencoders recover causal gene
programs at Nx the rate of foundation-model SAEs, on the field's own CRISPRi
benchmark — locating where causal regulatory logic actually lives."** That is
positive, novel, quantitative, and directly extends a very recent paper.

## 3. The experiment — an apples-to-apples causal head-to-head

Data: **Replogle 2022 K562 genome-scale Perturb-seq** (CRISPRi; ~2k essential or ~7k
genome-wide gene knockdowns, single-cell). K562 is a myeloid leukemia line, so the
hematopoiesis framing survives. Access: `pertpy.data.replogle_2022_k562_*` / Figshare.

- **Arm A (embedding-space; reproduce the field, self-run).** Embed Replogle control
  cells with scGPT (and/or Geneformer) -> train a TopK SAE on the residual-stream
  activations -> run the causal-grounding test. Target: reproduce ~6-10% (validates
  our pipeline against the published number).
- **Arm B (expression-space; the contribution).** Train a TopK SAE directly on the
  log-normalized expression of the same control cells -> run the IDENTICAL test.
- **Compare** the causal-grounding rate. Prediction: Arm B >> Arm A.

**Causal-grounding metric (pre-registered before running; precise):** For each tested
perturbation p = knockdown of gene t:
1. From the perturbation data, compute p's differential-expression signature.
2. Find the SAE feature whose decoder direction best aligns with p's signature.
3. "Causally grounded for p" iff that feature's activation shifts **specifically** in
   p's knockdown cells (vs control AND vs an abundance/size-matched set of other
   knockdowns), at a pre-set effect size and FDR.
4. Causal-grounding rate = fraction of tested perturbations that are grounded.

This is clean, quantitative, and inherits every rigor device we already built
(abundance-matched null, BH-FDR, effect-size floor, multi-seed, pre-registration).

### 3b. The triviality trap (and how the design defeats it) — critical

A sharp judge's objection: *"Of course an expression-space SAE beats a compressed
embedding at tracking a perturbation — it literally sees the genes. Knock down gene t
and t's mRNA drops; any feature containing t moves. That is tautology, not causal
grounding."* If the test were "does a feature track the KD signature," Arm B would win
trivially and the result would be worthless. The design must make the win
NON-trivial, on three axes:

1. **Train on CONTROL cells only.** The SAE never sees a single perturbation. Its
   programs come purely from co-expression in unperturbed cells. The causal test then
   asks whether those *unsupervised* programs behave as *causal units* when their
   putative regulator is knocked down. That is a real test of "are the programs
   causally real," not curve-fitting to perturbations.
2. **Score at the PROGRAM level, not the target gene.** Exclude the knocked-down gene
   t itself (and guides for it) from the feature's scored genes; require a *coherent
   response of the downstream module* (t's regulatory targets), so we measure
   regulatory logic, not the trivial self-drop of t.
3. **Demand SPECIFICITY.** KD of t must move t's program feature more than it moves
   random features AND more than other KDs move that feature (a perturbation x feature
   specificity matrix; grounded = the diagonal beats the off-diagonal at FDR).

With these, the residual expression-space advantage is exactly the SCIENTIFIC CLAIM,
not a confound: **causal regulatory structure is legible in expression and illegible
in foundation-model embeddings.** That is a substantive, contrarian critique of the
field's rush toward foundation-model interpretability — recoverable causal logic is
being compressed away. Stated plainly, defensible, and it is the thing that makes the
comparison worth running rather than obvious.

## 4. The ambitious extension (the "million times better" version)

If Arm B beats the baseline, go beyond measurement to a **built artifact**:

- **A causally-grounded gene-program atlas** for myeloid cells: the expression-space
  SAE dictionary, where each retained feature carries a *causal certificate* — the
  perturbation(s) that validate it. A validated resource, not just a number.
- **A causally-supervised SAE (method contribution):** add a perturbation-consistency
  term to the SAE objective (held-out knockdowns), training features to be causally
  grounded by construction, and show the grounding rate rises further. "The first
  causally-supervised gene-program SAE."

Either extension turns a comparison into a method/resource — the kind of thing that
wins.

**Honesty update (2026-07-05): the causally-supervised SAE was PROTOTYPED and it did
NOT work.** A first design (add a decoder-side reconstruction of held-out perturbation
DE signatures to the SAE objective) failed its own synthetic self-test — it did not
improve, and sometimes hurt, held-out causal grounding, on both cleanly-separable and
entangled synthetic programs. It is shelved pending a better formulation and is NOT
claimed as a result. This does not touch the centerpiece: the expression-vs-embedding
HEAD-TO-HEAD (sections 2-3) stands entirely on its own and is the winnable result.
The causal ATLAS extension (emit the grounded programs with their perturbation
certificates) remains viable and is low-risk, since it is just a formatting of the
head-to-head output.

## 5. What we KEEP from the current work (it becomes the credibility spine)

Nothing is wasted; it is repositioned:
- Pre-registration + SHA freeze, abundance-matched nulls, BH-FDR, effect-size floor,
  multi-seed, leave-one-out, PCA/NMF/ICA/GRN baselines, cross-species replication.
- The artifact-vs-signal finding is now the MOTIVATION: "single-method, non-causal
  program claims are method-dependent (we proved it) -> therefore causal validation is
  necessary -> here it is." Gate 0 sets up the causal centerpiece instead of standing
  alone.

## 6. Honest risks and mitigations

- **Arm B might also be low.** If expression-space grounding is ~10%, the win weakens
  to "causal grounding is hard for SAEs generally." Mitigation: strong mechanistic
  prior it will be higher (expression carries the direct signal); even a 2-3x lift is
  a result; and we run both arms internally so the comparison holds regardless of
  matching the exact published 6.2%.
- **Reproducing scGPT/Geneformer + their SAE + their metric is real engineering.**
  Mitigation: internal head-to-head (do not depend on their number); start with scGPT
  only; keep the metric identical across arms.
- **K562 is a leukemia line, not primary hematopoiesis.** Mitigation: frame as myeloid
  / K562 explicitly; optionally add a primary-hematopoiesis Perturb-seq later.
- **Compute.** scGPT inference + SAE training on 1e5-1e6 cells needs a real GPU
  (RunPod). Budget below.

## 7. Execution plan (RunPod)

Phase 0 (in-container, CPU, now): finalize the pre-registered causal metric as code +
synthetic tests; script the full pipeline; verify on a tiny synthetic perturbation set;
secure Replogle access path. NO heavy compute.

Phase 1 (RunPod GPU): download Replogle; run scGPT embeddings (Arm A); train
expression-space + embedding-space TopK SAEs (multi-seed); compute causal-grounding
rates; the head-to-head number.

Phase 2 (RunPod GPU): the extension (causal atlas and/or causally-supervised SAE).

Phase 3: reformat the whole writeup around the positive causal result; regenerate
figures; update the abstract and brief.

Suggested RunPod: 1x A100 40-80GB (or L40S) + ~100GB disk; PyTorch image. Rough
compute: scGPT inference on ~300k cells (hours), SAE training (minutes-hours/seed).

## 8. Success criteria (what "winning-caliber" means here)

- Arm A reproduces the ~6-10% embedding-space ceiling (pipeline validated).
- Arm B (expression-space) shows a clearly higher, statistically robust
  causal-grounding rate (target: a multiple of the baseline), pre-registered.
- At least one built artifact (causal atlas or causally-supervised SAE).
- The whole thing reproducible, pre-registered, and honestly caveated.

If Arm B does NOT beat the baseline, we report that honestly — but the mechanistic
prior and the internal design make a positive result the likely outcome.
