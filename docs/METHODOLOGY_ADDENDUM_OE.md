# METHODOLOGY ADDENDUM — Overexpression causal grounding on the TF Atlas

Added 2026-07-18. This addendum documents the lens change from CRISPRi knockdown (Replogle)
to TF **overexpression** (the Joung Transcription-Factor Atlas), the metric generalization
that supports it, and the pre-registration for the run. It extends `PREREGISTRATION.md` and
the causal-grounding spec (`config/causal_grounding_spec.json`); nothing here supersedes the
frozen Gate-0 work.

---

## 1. Motivation (why this pivot, honestly)

The Replogle CRISPRi runs (K562-essential, RPE1) established, rigorously, that expression-
space SAEs recover genuine but **housekeeping/proliferation** causal programs, not lineage
regulatory logic, and no better than PCA/NMF (see `COMPONENT2_RESULTS.md`). The dominant
reason is the **substrate**: these are ESSENTIAL-gene screens. They knock down housekeeping
genes; the hematopoietic lineage transcription factors (GATA1, SPI1, TAL1, CEBPA, RUNX1,
KLF1) are not perturbed, so there is no lineage signal to ground. The null was a
wrong-substrate result, not proof the method cannot find regulatory logic.

The fix is to test on data where lineage/regulatory TFs ARE perturbed and the signal is
documented to exist.

## 2. The dataset

**Joung et al. 2023, "A Transcription Factor Atlas of Directed Differentiation," Cell**
(GEO **GSE216481**). A MORF pooled **overexpression** screen: ~3,550 TF open reading frames
introduced into human stem cells, ~254,519 single cells. It contains the hematopoietic
lineage TFs as named perturbations.

Independent validation that the signal is present: **Jain & Sharma 2026 (arXiv 2604.02511)**
re-analyzed GSE216481 and recovered **TF-specific signatures for 59 of 61 testable TFs**
(vs 27 by naive one-vs-rest), using embryoid-body cells as an external baseline plus
background subtraction to handle the deposit's missing intra-pool GFP/mCherry controls. So a
positive result is achievable here, unlike the essential screens.

## 3. Why overexpression is a better causal test

- **Cleaner causal direction.** TF ON -> its program ON is a direct, specific activation.
  It **sidesteps the essential-gene / global-collapse confound** entirely: overexpressing
  GATA1 drives the erythroid program; it does not cause the broad transcriptional collapse
  that made essential-gene knockdowns ground non-specifically.
- **The biology is present.** Lineage TFs are the perturbations, so "does the SAE recover
  GATA1 -> erythroid?" is finally a question the data can answer.
- **A specificity anchor.** Because each perturbation is a single named TF with a
  literature-known target program, known links (GATA1->globin/heme, SPI1->myeloid, MYOD->
  muscle) are hand-checkable positive controls.

## 4. The metric generalization (`direction` parameter)

`causal_grounding(..., direction=)` — `'down'` is the original CRISPRi metric (unchanged
default); `'up'` is overexpression. The sign flips, all cross-fit / FDR / gate machinery
identical:

| step | down (knockdown) | up (overexpression) |
|------|------------------|---------------------|
| match target | `-DE` (program SUPPRESSED) | `+DE` (program ACTIVATED) |
| Mann-Whitney (split B) | `alternative='less'` | `alternative='greater'` |
| AUC = P(perturbed > control) | expect LOW | expect HIGH |
| suppression/activation floor | `auc <= auc_floor` (0.45) | `auc >= 1 - auc_floor` (0.55) |
| directional magnitude `move` | `0.5 - auc` | `auc - 0.5` |
| column specificity tail | lower | upper |
| effect-size control (M1) | `move` outlier vs effect-matched peers | identical (uses `move`) |
| perturbed gene t excluded | yes (self-drop trivial) | yes (own ORF transcript trivial) |

Cross-fitting (match on split A, test on held-out split B), BH-FDR, the effect-size control
(default on), and the label-shuffle null all carry over unchanged.

**Validation (synthetic, `tests/test_causal_grounding.py::test_overexpression_direction`):**
in a world where each TF's overexpression ACTIVATES its program, `direction='up'` grounds
6/6 TFs with 0 background, while the wrong-sign `direction='down'` grounds 0/6 -- proving the
generalization is real and non-trivial. The pipeline was smoke-tested end-to-end in `up`
mode (SAE training + grounding).

## 4b. TopK feature-competition subtlety (important, and why `up` is the correct direction)

A TopK SAE keeps only the k largest features active per cell. So when an overexpressed TF
drives its program feature strongly, that feature wins a top-k slot and can push COMPETITOR
features OUT of the top-k -- making them appear "suppressed" in those cells. Consequence,
verified on synthetic (trained SAE): on overexpression data the WRONG-sign `direction='down'`
grounds the competition-suppressed competitors (6/6 in a toy), while `direction='up'` grounds
the genuinely activated program (3/6 in the same toy, limited only by the toy SAE's imperfect
recovery of 6 programs). This is NOT a metric error -- the metric-level test on PLANTED
activations (no competition) cleanly gives up=6/6, down=0/6.

Implications, pre-registered:
- **Use `direction='up'` on overexpression data.** It targets the activated program; the
  `down` direction would report the top-k competition artifact.
- **Guard with the match step + positive controls.** The match requires the grounded feature
  to align with the TF's OWN +DE signature, so `up` grounds the TF to its actual program, not
  a random activated feature. The named positive controls (GATA1->erythroid etc., Section 6)
  confirm this on real data.
- Report this competition property explicitly; it is itself a finding about TopK SAEs on
  perturbation data.

## 5. Differentiation from the prior re-analysis (novelty)

Jain & Sharma (2604.02511) recovered TF signatures with **one-vs-rest differential
expression** -- a gene list per TF. Our contribution is different and complementary:
1. **Interpretable PROGRAMS, not gene lists.** SAE features are (aim to be) monosemantic
   programs; grounding asks whether an unsupervised, control-trained program behaves as a
   causal unit under its TF's overexpression.
2. **Program-level structure.** Which TFs CONVERGE on shared SAE features -> a modular map
   of TF regulation that per-TF DE cannot produce.
3. **A confound-aware causal metric** (cross-fitting, effect-size control, calibrated
   nulls, pre-registration) rather than DE + background subtraction.
4. **The head-to-head frame** (expression vs foundation-model embedding) still applies.

## 6. Pre-registration for the run (freeze BEFORE running)

- **Metric:** the `direction='up'` causal-grounding spec (amend `config/causal_grounding_spec.json`
  to v7 and re-freeze the SHA) BEFORE the first atlas run.
- **Baseline/control:** use the deposit's negative controls if recoverable, else the
  embryoid-body external baseline (as in 2604.02511). Document which, and map it to
  `CONTROL_LABEL`. This choice is frozen before running.
- **Primary readouts (declared in advance):**
  1. Overall grounding rate (fraction of tested TFs whose program is causally grounded),
     vs the label-shuffle null.
  2. **Named positive controls:** do GATA1, SPI1, TAL1, KLF1, CEBPA ground, and to
     literature-consistent programs? (Hand-checkable.)
  3. Program convergence: how many distinct SAE features do the grounded TFs map to
     (module structure).
  4. SAE vs PCA/NMF at matched settings (does the SAE add over matrix factorization HERE,
     where signal exists -- the ladder that came back flat on the essential screens).
- **Honest success criterion:** a positive is "canonical lineage TFs ground to their known
  programs above the null, and the SAE adds interpretable module structure." A null (TFs do
  not ground even here) would be a strong, publishable statement that SAE program grounding
  fails even where the perturbation signal is documented -- but the mechanistic prior
  (overexpression drives specific programs; 59/61 recoverable by DE) favors a positive.

## 7. Data-access plan (the tractability step)

GSE216481 is a GEO deposit (count matrix + guide-barcode -> TF assignment), not a pertpy
one-liner. Steps (to be scripted next):
1. Obtain the processed atlas (GEO supplementary h5ad / a cellxgene or Figshare mirror if
   available; else raw matrix + barcode map).
2. Build `obs['TF']` = overexpressed TF per cell from the guide-barcode assignment; set the
   control label to the chosen baseline (Section 6).
3. Feed to `load_perturbseq` (normalize_total + log1p, HVG-union of TF-target genes, memory
   bound with `--n-hvg`) and run with `--direction up`.
This is the one hands-on step; it is deliberately separated from the metric work (which is
built and tested) so a data hiccup cannot corrupt the analysis.

## 8. Honest limitations

- The atlas is hESC-context overexpression, not primary hematopoiesis. Framing: "TF-driven
  programs, including the hematopoietic lineage TFs," not "in vivo hematopoiesis."
- Overexpression is not physiological dosage; a TF may drive programs it would not at native
  levels. Grounding a KNOWN link (GATA1->erythroid) is the guard against over-reading.
- The re-analysis competitor exists; our differentiation is the interpretable-program +
  confound-aware-metric angle (Section 5), which must be stated explicitly.
- No result is claimed yet. This addendum is the pre-registration; the number comes from the
  frozen run.
