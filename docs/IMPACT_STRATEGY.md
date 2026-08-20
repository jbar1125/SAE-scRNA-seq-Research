# IMPACT STRATEGY — making this project impressive and influential (honest)

Written 2026-07-18, after the causal arm came back a negative (expression-space SAE
grounding tracks perturbation effect-size, not lineage-regulatory logic; no master TF
robustly grounds on Replogle-K562-essential). This document answers one question: given
where we actually are, what are the highest-leverage moves to make this project genuinely
impressive and influential? No spin — the plan leans INTO the negative, because that is
where the real, defensible contribution now lives.

---

## 1. Why the old headline can't carry it (and why that's OK)

The dream was "expression-space SAEs recover causal regulatory logic that foundation
models discard, beating the 6.2% benchmark." Tonight that's on the ropes: on this dataset
the metric finds big-effect housekeeping knockdowns, not regulators. Chasing that headline
harder is the low-probability path.

But three things are now TRUE and STRONG:
- Gate 0: SAE modularity claims are method-dependent (rigorous, frozen, cross-species).
- Component 2: SAE causal grounding is confounded by perturbation effect-size, and once you
  see that, the naive "X% grounded" numbers in the field look shaky.
- Independent convergence: Kendiukhov (embedding space) concludes SAE features encode
  co-expression, not causal logic. We reach the same wall from EXPRESSION space.

That is a coherent, important, contrarian story: **interpretability methods are being
trusted to recover causal biology faster than anyone has checked whether they can — and a
rigorous check says they largely can't, in either representation.** That is more
influential than one more benchmark number.

## 2. The reframe (the impressive thesis)

> **"When does single-cell interpretability recover CAUSAL biology? A pre-registered,
> confound-aware audit — showing that SAE gene programs are method-dependent (structure)
> and effect-size-confounded (causality) in both expression and foundation-model space,
> and a benchmark that separates genuine causal grounding from artifact."**

The project's identity becomes **the auditor and the standard**, not one more method. This
is the PLI-Analyzer shape (audit a hyped AI method, rigorously, and report what it can't
do) — a proven STS-top-10 template — applied to the single-cell SAE wave.

## 3. Why a rigorous negative is influential here

- **Field gap:** there is NO standard causal-validation benchmark for single-cell
  interpretability. Papers assert "our features are interpretable/biological" and stop.
  A confound-aware causal test fills a real hole.
- **Timing:** the SAE-single-cell wave is 2025-2026 and accelerating (Kendiukhov cluster,
  CytoSAE, etc.). A rigorous skeptic with tools, RIGHT NOW, is maximally useful.
- **Saves effort:** a clean "don't trust SAE features as causal without this test" result
  redirects a whole subfield. Negatives that redirect fields get cited.
- **Judge-legible:** "everyone's using AI interpretability on cells; does it actually
  recover the biology? I built the test and found where it breaks" — repeatable by any
  non-specialist judge.

## 4. The highest-leverage moves (prioritized)

### M1 — The effect-size-residual grounding metric (THE methodological centerpiece)
The single most valuable move, and it does double duty. Tonight's finding is that grounding
tracks perturbation effect magnitude. Turn that from an embarrassment into the contribution:
- **Diagnose it quantitatively:** plot grounding (and matched-feature AUC) vs each
  perturbation's total transcriptional effect size. Show the correlation. This IS a finding:
  "naive causal-grounding rate is largely a readout of effect size."
- **Fix it at the metric level:** define grounding on the RESIDUAL — a perturbation is
  causally grounded only if its matched feature is suppressed MORE than its overall effect
  magnitude predicts (regress out / stratify by effect size, or compare to effect-matched
  perturbations). 
- **Why this is impressive:** it is a genuine methodological contribution nobody has made
  (the embedding-space papers don't control for it either), AND it is the last honest shot
  at a positive: a specific regulator like GATA1 may have a modest TOTAL effect but a strong
  SPECIFIC feature suppression → high residual → grounded even when raw grounding misses it.
  If the residual metric surfaces lineage regulators that the raw metric buried, that is the
  rescued, stronger result. If it doesn't, the effect-size-confound finding stands on its own.

### M2 — Comprehensive convergent evidence: the head-to-head + the representation ladder
Convert one negative into a general statement. Run (all already wired: `--rep`):
- Arm A embedding (scGPT/Geneformer) at matched settings.
- The ladder: expression → PCA → NMF → embedding.
Outcomes, all publishable: (a) flat-and-low everywhere → "SAE causal grounding fails
regardless of representation" (strong, convergent with Kendiukhov); (b) a decay curve →
"causal legibility lives in expression"; (c) expression clearly > embedding even if low →
the original comparative claim survives in weakened form. A curve/《comparison beats a
single number and de-risks the whole thing.

### M3 — "Interpretability illusions in single-cell": a confound taxonomy
Package the artifacts as a catalog, each demonstrated on real data:
- method-dependence (Gate 0),
- effect-size confounding (Component 2, M1),
- dataset-composition traps (the essential-library skew; the wrong-background enrichment
  error I made and caught),
- marker-curation artifacts (the globin-symbol bug from Gate 0).
A named taxonomy of "ways SAE interpretability misleads, and the control for each" is a
genuinely useful, citable contribution (cf. the "interpretability illusions" line in NLP
mech-interp). It also turns every mistake this project made into demonstrated pedagogy.

### M4 — The GATA1 boundary case (make the negative scientifically rich)
A negative with a boundary condition is far stronger than a blanket one. GATA1 is the one
lineage TF that flickered. Investigate: is its matched feature the erythroid/heme module
(HBB/HBA/ALAS2/GYPA/KLF1)? Under the residual metric (M1), does it ground cleanly? The
statement "SAE grounding recovers regulators only when their program is dominant, single-
regulator, and high-effect (e.g. GATA1 in erythroleukemia) and fails otherwise" is a
precise, mechanistic, memorable boundary — much better than "it doesn't work."

### M5 — Ship "CausalGround" as an open benchmark (the influence play)
Package `causal_grounding.py` + `causal_pipeline.py` + the reference sets + the effect-size
control into a clean, documented, pip-installable benchmark: given any representation of any
Perturb-seq dataset, report the confound-aware causal-grounding rate with calibrated nulls.
Benchmarks are the highest-leverage influence artifact — they get adopted and cited, and
they make the project infrastructure, not a one-off. Include a small leaderboard table
(expression/PCA/NMF/scGPT/Geneformer).

### M6 — Preprint + position as an independent extension of Kendiukhov
A citable bioRxiv / ML-for-genomics workshop paper (MLCB, LMRL) of the framework + the
negative + the benchmark. Frame explicitly as: independent, expression-space extension of
Kendiukhov's embedding-space finding, plus the effect-size-confound methodology they lack.
A high-schooler independently reaching and rigorously extending a March-2026 result is
itself a notable story. A preprint massively strengthens any competition entry.

### M7 — Communication as a first-class deliverable
- One clean narrative (kill the v-number tangle; the arc is: hyped method → does it recover
  causal biology? → rigorous test → here's where it breaks + the tool).
- The money figures: the effect-size-vs-grounding scatter (M1), the representation ladder
  (M2), the confound taxonomy (M3).
- An accessible thread/blog explaining the negative to a general audience. Science
  communication of a rigorous negative is rare and memorable.

## 5. Three shapes the final project can take (all honest, ranked by ceiling)

1. **"The Auditor + The Benchmark" (recommended; highest honest ceiling).** M1+M2+M3+M5+M6.
   A pre-registered, confound-aware framework + open benchmark showing SAE interpretability
   is method-dependent and effect-size-confounded across representations. Contrarian,
   tool-backed, field-relevant, judge-legible. Wins on rigor + utility, not on a fragile
   positive.
2. **"The Rescued Positive" (higher ceiling, lower probability).** If M1's residual metric
   or a lineage-TF dataset surfaces genuine regulator grounding, the head-to-head becomes a
   real "where causal logic lives" result. Pursue M1/M4 first; if the signal appears, pivot
   here. If not, fall back to shape 1 with no loss.
3. **"The Rigor Floor" (safety net, already banked).** Gate 0 alone — method-dependence,
   cross-species, 5 methods — is a complete honest project. This exists regardless.

The move: pursue shape 1 as the plan, run M1/M4 as the option on shape 2, with shape 3 as
the guaranteed floor. You cannot end below "complete and rigorous."

## 6. Honest risks / when to abandon a direction

- **M1 residual metric may still be flat** → then the effect-size-confound finding is the
  contribution; do not keep patching. (We have already burned two metric patches; M1 is the
  principled one, but it is the LAST metric iteration before we accept the negative.)
- **The lineage-TF dataset (M4/shape 2) may not exist / be reachable** → the boundary-case
  analysis on GATA1 within Replogle still stands.
- **Benchmark adoption is not guaranteed** → but shipping it is cheap and the artifact
  itself is the value, adoption or not.
- **Do not oversell convergence with Kendiukhov** → it is independent corroboration +
  extension, not a claim of priority.

## 7. The one-sentence pitch (the test of whether this is impressive)

"Everyone is using sparse autoencoders to interpret single-cell data; I built a
pre-registered, confound-aware test of whether those features actually recover CAUSAL
biology, and found they track perturbation effect-size — not regulatory logic — in both
gene-expression and foundation-model space, and I shipped the benchmark so the field can
check its own claims."

If that sentence lands, the negative is the win.
