# Component 2: perturbation analyses and evidence status

**Current audit, 2026-10-18:** no new lineage-regulatory discovery or expression-
versus-foundation-model advantage is established. The original causal result JSONs
were lost with ephemeral GPU instances. The dated tables and notebook reports
below preserve the research history, but are not a substitute for archived reruns.

## Current interpretation

| Analysis | Recorded result | What it establishes and what remains missing |
|---|---|---|
| Replogle K562, original five-seed run | 22/162, 31/162, 42/162, 31/162, 32/162 grounded perturbations | Historical v3 results, mainly housekeeping/machinery hits; original per-perturbation JSONs unavailable. |
| Essentiality check after correcting the background | v4 Hart-core-essential-tagged 6/51, other tested TFs 14/111 | Similar observed grounding fractions; the untagged group is not a validated nonessential reference set. The former 10-fold enrichment used the wrong background and is withdrawn. These v4 counts are not a reanalysis of the entire v3 five-seed set. |
| Column-specificity amendment | GATA1 decreased from at least 3/5 seeds to 2/5; housekeeping hits remained | The proposed fix did not establish the intended specificity and was disabled by default. |
| Effect-size control, K562 seed 0 | Raw and controlled sets both 22/162; Spearman rho 0.28 | This control removed no calls in the recorded comparison. The broad claim that grounding merely tracks total effect size is unsupported; universal freedom from confounding is also not established. |
| Representation ladder, K562 | Mean 0.121 expression, 0.130 PCA-input, 0.111 NMF-input | Every arm fits an SAE. The experiment compares SAE input representations, not SAEs against plain PCA/NMF. NMF-input has one recorded seed. |
| RPE1 ladder | Recorded rates 0.076 expression, 0.070 PCA-input, 0.076 NMF-input | Similar observed rates and housekeeping/chromatin hits in a second cell line. This is a historical report requiring archived regeneration. |
| Norman CRISPRa rate | Recorded perturbation-identity shuffle p = 0.902 | No evidence that the rate distinguishes perturbation identity beyond generic perturbed-versus-control structure under that shuffle. It is not proof that all biological signal is absent. |
| Norman regulon recovery v2 | 0/14 TFs pass; 0/42 TF-by-seed comparisons exceed the simulated chance 95th percentile | Unsuccessful single-TF curated-regulon alignment under this analysis. TRRUST coverage and the selection baseline limit the claim; exact real-result JSONs are missing. |
| Program-level modulation | Implemented; full real-data discovery not established | Existing plans and synthetic checks do not establish a positive result. Full-shuffle counts and a feature-background percentile are exploratory diagnostics; perturbation identity requires a separate test and does not establish pathway specificity. |
| Foundation-model embedding arm | Not run | The published 6.2% number is not a matched comparison. |

These causal numbers are historical reports traced to `LAB_NOTEBOOK.md` and
`PROJECT_MASTER_SUMMARY.md` sections 8-9. They were not recomputed from missing
original outputs during this documentation correction. In particular, similar
rates do not constitute a statistical equivalence result or a proof of no signal.

## Training and analysis scope

Every historical causal SAE was fit on **all cells**, without perturbation labels.
Control cells were used to compute feature-gene associations. The code did not
implement the original control-only training statement. The erratum is recorded in
`config/causal_grounding_spec.json` (current SHA `aae1a9d9...`).

Cross-fitting separates perturbation-specific matching from suppression testing,
but test cells still contributed to SAE training. The current control-only mode
is a sensitivity analysis: its actual-training synthetic test recovers 0/8 planted
regulators, compared with 8/8 for all-cell fitting. That result limits the current
validation; it does not prove all control-only approaches are unsuitable.

Causal metric amendments were made after inspecting real results. Describe these
analyses as exploratory and report their chronology. The single-perturbation
results neither establish nor causally validate a universal hematopoietic modularity
direction. K562 and RPE1 are cell lines, and transfer to primary hematopoiesis is
not established.

The Norman erythroid-like feature is a descriptive lead. GATA1 is absent from the
Norman perturbation library, so it did not fail a Norman grounding test. Previously
noticed CBL/UBASH3B/PTPN12 responses require appropriately calibrated testing,
primary-literature review, and independent validation before any new mechanism is
claimed. Program-level selection, multiplicity, perturbation identity, and generic
response axes must be checked explicitly.

## Publication and reproduction requirements

Archive reruns with code and data hashes, resolved environment, exact arguments,
training diagnostics, per-perturbation results, negative results, and null outputs.
Record seed variability and uncertainty; use comparable baselines. The original
causal diagnostic figure was based on hardcoded recorded numbers and contained a
wrong-background enrichment panel; it cannot serve as independent result evidence.
See [REPRODUCIBILITY.md](REPRODUCIBILITY.md).

## Historical record: 2026-07-18 K562 analyses

The sections below preserve the sequence of provisional interpretations and
self-corrections. Statements within this record about essentiality confounding,
effect detectability, validated GATA1 controls, unrun RPE1, or definitive regulatory
mechanisms are **superseded by the current interpretation above**. Retaining them
records how the hypotheses changed; it does not endorse them as current claims.
The historical phrase “non-essential 14/111” below means not in the Hart CEGv2
tag set, not validated nonessential genes. The representation-ladder paragraph and
exact-regeneration claim are corrected
explicitly because they misdescribed the implementation and evidence status.

---

> **STATUS: PRELIMINARY. NOT FROZEN.** Arm B (expression-space) only — the
> embedding-space Arm A that makes this a head-to-head is NOT yet run. Single cell
> line (K562). The grounded set is essential-gene-confounded (see below). Do NOT
> report these as a final finding or SHA-freeze them. Recorded 2026-07-18 as an
> honest interim checkpoint.

## Run configuration

- Data: Replogle 2022 K562 essential Perturb-seq (`pertpy.data.replogle_2022_k562_essential`), 310,385 cells.
- Genes: top-2000 HVG UNION perturbed-TF genes = 2,136.
- Scored perturbations: 162 TFs (`config/tf_lists/hs_hgnc_tfs.txt`, >= min_cells).
- SAE: TopK, GPU (cuda). **ERRATUM 2026-10-18:** this line originally said "trained on control
  (non-targeting) cells only". That was the design intent, not what ran: every caller fit the SAE
  on ALL cells (perturbation labels unseen). Every number in this file is from all-cell SAEs. See
  the erratum in `config/causal_grounding_spec.json` and LAB_NOTEBOOK 2026-10-18.
- Metric: cross-fit Mann-Whitney grounding, spec SHA `b766d21c` (`config/causal_grounding_spec.json`).

## Seed stability (latent 2048, k 32)

| seed | grounded | rate |
|-----:|:--------:|-----:|
| 0 | 22/162 | 0.136 |
| 1 | 31/162 | 0.191 |
| 2 | 42/162 | 0.259 |
| 3 | 31/162 | 0.191 |
| 4 | 32/162 | 0.198 |
| **mean** | | **0.195 (sd 0.039, range 0.136-0.259)** |

All five seeds sit far above the label-shuffle null (~0.00). The finding "expression-
space SAE features show real perturbation-specific grounding" is SEED-ROBUST. The exact
rate is noisy (~20% coefficient of variation); the effect is not.

## Seed-stable grounded set (grounded in >= 3/5 seeds): 26

CNOT3, CSNK2B, CXXC1, DMAP1, E4F1, **GATA1**, GTF2A2, HINFP, HSF1, ILF2, MAX, NCBP2,
NELFB, POLD2, PSMD12, RFC2, RFC3, SFPQ, SMUG1, SRP9, TBP, TFDP1, THOC2, TIMELESS, ZMAT2,
ZNF335.

Honest categorization:
- **~24/26 are essential / general-transcription machinery / non-TF:** TBP, GTF2A2
  (basal apparatus), POLD2/RFC2/RFC3 (replication), SRP9 (signal recognition particle),
  PSMD12 (proteasome), THOC2/NCBP2/SFPQ (RNA processing/export), CSNK2B (kinase), SMUG1
  (base-excision repair), CNOT3, DMAP1, TIMELESS, ZMAT2, NELFB, ILF2, HINFP.
- **Notable real regulators that DO ground robustly:** **GATA1** — the master erythroid
  TF — grounds in >=3/5 seeds in an erythroleukemia line. This is a positive control
  PASSING: the metric catches a genuine lineage regulator. Plus MAX (MYC network), HSF1
  (heat-shock), E4F1, CXXC1, ZNF335 as defensible sequence-specific TFs.
- **Interpretation:** regulatory signal is PRESENT but BURIED under essential-gene
  transcriptional collapse. The three specificity fixes (S.1-S.3 in `UPGRADE_BACKLOG.md`)
  are required before any regulatory-specificity claim.

## Hyperparameter sensitivity (seed 0)

| config | grounded | rate |
|:------|:--------:|-----:|
| latent 512, k 32 | 64/162 | 0.395 |
| latent 1024, k 32 | 46/162 | 0.284 |
| latent 2048, k 32 | 22/162 | 0.136 |
| latent 2048, k 16 | 31/162 | 0.191 |
| latent 2048, k 64 | 58/162 | 0.358 |

- **Dictionary size (at k=32) is inversely related to the rate:** 512 -> 0.395, 1024 ->
  0.284, 2048 -> 0.136. Smaller/broader dictionaries ground MORE, consistent with the
  confound: broad features capture global-state axes that essential knockdowns suppress.
  A higher rate is therefore NOT "better" — it can mean less specific features.
- **k (at latent 2048) trends up with sparsity budget** (k64 = 0.358 highest), noisy at
  single seed.
- **The rate spans 0.136-0.395 from hyperparameters ALONE** (~3x). The absolute grounding
  rate is not a standalone quotable quantity.

## Honest conclusions

1. Grounding is real and seed-robust (all seeds >> null), but the ABSOLUTE rate is
   hyperparameter-dominated. The expression-vs-embedding head-to-head MUST fix (latent, k)
   identical across arms and report the rate across a grid, not a single number. Comparing
   the raw rate to Kendiukhov's 6.2% is invalid until Arm A is run at matched settings.
2. The grounded set is essential-gene-dominated, so the current rate largely reflects
   essential-knockdown transcriptional collapse, not regulatory specificity — BUT GATA1
   (and a few other real TFs) ground robustly, so the metric does detect genuine
   regulators; they are swamped, not absent.
3. Smaller dictionaries grounding more reinforces conclusion 2: the rate is inflated by
   broad, global-state features.

Diagnostic figure: `docs/figures/fig_component2_diagnostics.png` (generated by
`src/make_component2_figures.py`) -- panel A the hyperparameter dominance, panel B the
essential-gene enrichment. Marked preliminary/pre-v4-rerun.

## Addendum: quantified confound analysis + falsifiable predictions (2026-07-18)

Cross-referencing the 26 seed-stable hits against Hart CEGv2 (core-essential) and NEGv1
(reference-nonessential):

- **Two-sided essentiality bias.** Grounded genes are ENRICHED for core-essential
  (8/26 = 31% vs a 3.0% pool base rate, ~10x, hypergeometric p~4e-7) AND DEPLETED of
  reference-nonessential (0/26 = 0% vs a 5.3% pool rate). Both directions point the same
  way: grounding tracks essentiality, not regulatory role.
- **CEGv2 undercounts the confound.** By curated function, 20/26 are housekeeping/machinery
  (replication POLD2/RFC2/RFC3/TIMELESS, basal TBP/GTF2A2, RNA-proc SFPQ/NCBP2/THOC2/ILF2/
  ZMAT2, proteasome PSMD12, kinase CSNK2B, SRP9, chromatin CXXC1/DMAP1, ...); only 6/26 are
  sequence-specific regulators (GATA1, MAX, HSF1, E4F1, TFDP1, ZNF335), and only GATA1 is
  lineage-specific. So the true housekeeping fraction (~77%) exceeds the strict CEGv2 31%.
- **Honest background caveat.** The correct enrichment denominator is the 162 TESTED TFs
  (which are HVGs in K562 and thus already essential-enriched), not all 1839 pySCENIC TFs;
  so the ~10x magnitude is an upper bound and the true enrichment is smaller. The
  nonessential DEPLETION (0/26) is background-robust evidence regardless.

**Two falsifiable predictions for the v4 (column-specificity) re-run** -- recorded now so
they cannot be retrofitted:

- **P1 (mechanism of the sweep).** The v3 rate's inverse dependence on dictionary size
  (512->0.395 ... 2048->0.136) is hypothesized to be DRIVEN by the confound: smaller
  dictionaries learn broader, more polysemantic features -> more global 'cell-health' axes
  -> more essential-KD collapse registers as grounding. PREDICTION: under v4, the
  dictionary-size dependence FLATTENS substantially. If it does, the confound explains the
  sweep. If it persists, there is an additional (capacity) cause worth investigating.
- **P2 (seed churn).** The v3 grounded set churns (seed rates 22/31/42/31/32; only 26 of a
  larger union are >=3/5-stable), consistent with threshold-borderline confound hits.
  PREDICTION: under v4 the surviving set is MORE seed-stable (a higher fraction of grounded
  perturbations are >=3/5), because real regulators are reproducible and confound hits are
  not.

If P1 and P2 both hold on the v4 re-run, the confound diagnosis is confirmed end to end.

## UPDATE 2026-07-18 (later): the column-specificity fix FAILED; two corrections

Ran the v4 (column-specificity) metric on real Replogle, single-seed then 5-seed. Honest
outcome, leading with the bad news:

- **The fix did not work and made things worse.** v4 grounded 20-34/162 (mean ~0.157),
  the SAME housekeeping/machinery set. The 5-seed-stable set (22) is entirely housekeeping;
  **GATA1 dropped from >=3/5 (v3) to 2/5 (v4)** -- column specificity removed the one
  lineage regulator -- and every other master TF (SPI1, TAL1, RUNX1, KLF1, CEBPA, GFI1,
  GATA2, MYB) is 0/5. Root cause: distinct essential KDs hit DISTINCT stress features (each
  is its own top suppressor, so the gate passes them), while GATA1's erythroid program in
  K562 has multiple regulators (so GATA1 looks non-specific and gets cut). Rolled back to
  default OFF (spec v5).
- **Correction to the "~10x essential enrichment" stat.** That used the wrong background
  (all 1839 TFs). The v4 tag-breakdown gives the right one directly: essential 6/51 =
  11.8%, non-essential 14/111 = 12.6% -- **about EQUAL, no essential preference.** The
  tested pool is ~31% essential because this is the Replogle ESSENTIAL library (essential
  knockdowns by construction); the grounded set is housekeeping-heavy because the POOL is.

**Honest bottom line for the causal centerpiece as currently designed:** on this dataset,
expression-space SAE grounding tracks perturbation **effect-detectability** (the biggest,
cleanest transcriptional collapses = housekeeping/essential genes), NOT lineage-regulatory
specificity. No master hematopoietic regulator is robustly grounded. This is a real,
honest negative for the "expression recovers regulatory logic" hypothesis on Replogle-K562-
essential; it is not the benchmark-beating headline that was hoped for.

**What is NOT yet ruled out:** (a) the head-to-head -- does the EMBEDDING arm do even worse?
The comparison, not the absolute rate, is the deliverable and it is unrun; (b) a dataset
that actually perturbs lineage TFs (genome-wide Replogle / RPE1 / a hematopoietic
Perturb-seq) -- the essential library barely contains regulators like SPI1/CEBPA/GATA2.

## UPDATE 2026-07-18 (M1 result): effect-size is NOT the confound; grounding is genuine but housekeeping

Ran the raw vs effect-size-controlled metric on real Replogle (seed 0). Result:

- **The effect-size control changed nothing:** raw 22/162 and controlled 22/162 are the
  IDENTICAL set (0 dropped, 0 "rescued"). Note: the gate is an AND-filter, so it can only
  remove, never add -- "rescue" was a mis-framing on our part.
- **The dropped-zero is the finding:** all 22 grounded perturbations PASS the effect-size
  control -- they suppress their matched feature MORE than effect-matched peers do. So they
  are NOT effect-size artifacts. The diagnostic confirms a WEAK effect-size relationship
  (Spearman rho(effect, suppression) = 0.28; mean effect grounded 3.81 vs not-grounded
  3.07). This REFUTES the earlier "grounding just tracks effect-detectability" worry -- that
  claim was overstated.
- **GATA1 still does not ground, and not because of a confound:** its effsize_p = 0.435 (not
  even an effect-relative suppression outlier). The lineage signal is simply not present/
  suppressed in this dataset; K562 is a fixed leukemic state and the essential library does
  not perturb lineage TFs. (MYC is an effect-relative outlier, effsize_p 0.018, but still
  fails the absolute suppression floor.)

**Honest reinterpretation (supersedes the "all confound" read above):** the expression-space
SAE genuinely recovers causally-grounded, effect-size-controlled gene programs (22/162,
above null) -- but they are PROLIFERATION/HOUSEKEEPING programs, the biology that is variable
and perturbed in K562-essential, not hematopoietic lineage programs. The method works; the
dataset does not contain the lineage biology the project set out to test. The two real next
steps are therefore (1) the head-to-head (do embeddings ground these same programs, or
fewer?) and (2) a dataset that perturbs lineage regulators.

M1 status: the effect-size control stays in the metric (default on) as a validated
robustness CHECK -- on this data it is a no-op, but it defends the grounded set against the
effect-size critique and would filter the confound in datasets where it is present.

## UPDATE 2026-07-18 (representation ladder): similar recorded rates for SAEs on different input representations

Ran the causal metric with the SAE trained on expression vs PCA-128 vs NMF-128 of the same
data (`--rep`). Result:

**Interpretation corrected 2026-10-18 after inspecting the implementation:** every
arm trains a TopK SAE. `build_representation()` constructs expression, PCA-128, or
NMF-128 inputs, and `run()` then trains an SAE on that input. The original wording
claiming equivalence to *plain* PCA/NMF was incorrect.

- Recorded mean grounding rates: SAE on expression 0.121 (seeds
  0.136/0.093/0.136), SAE on PCA-128 0.130 (0.130/0.111/0.148), SAE on
  NMF-128 0.111 (seed 0 only). These are similar observed rates, without a formal
  equivalence test; the NMF-input comparison is under-replicated.
- Recorded seed-0 grounded counts were expression 22, PCA-input 21, NMF-input 18.
  Expression/PCA-input overlap was 6 and expression/NMF-input overlap 2. These
  are historical reports, not currently inspectable per-perturbation JSON evidence.
- This experiment does not establish that an SAE offers no benefit over direct
  matrix-factorization features. Direct PCA/NMF feature baselines and a matched
  foundation-model embedding arm remain separate, uncompleted comparisons.

## Historical next-step list (superseded; see current status above)

- **Arm A** (scGPT / Geneformer embedding-space SAE) at MATCHED hyperparameters — the
  actual head-to-head. This is the deliverable; the number above is not.
- The three specificity fixes S.1-S.3 (Lambert 2018 sequence-specific TF list; a
  program-coherence specificity filter that actually bites; an essential-gene control
  split via DepMap) — see `UPGRADE_BACKLOG.md`.
- Second cell line (RPE1), effect-size CIs, cross-fit K-repeat stability.

## Historical reproduction command (not evidence of exact regeneration)

The original JSONs were not preserved. This command records the intended run shape;
current defaults include later amendments and must be matched to the specific historical
revision before a reproduction is claimed. Exact regenerated outputs have not been verified:

```bash
for s in 0 1 2 3 4; do
  python3 src/causal_pipeline.py --adata replogle.h5ad --pert-col gene \
    --control-value non-targeting --rep expression --latent 2048 --k 32 --seed $s \
    --n-hvg 2000 --tf-list config/tf_lists/hs_hgnc_tfs.txt --n-shuffle 50 \
    --out causal_out/expr_grounding_seed$s.json
done
```
