# COMPONENT 2 — NEXT DIRECTIONS (post-first-run analysis + execution map)

Written 2026-07-18 after the first real Arm-B run + the overnight multi-seed/sweep +
the metric-hardening session. This is the RESULTS-DRIVEN plan: what the numbers dictate,
in sequence, with decision branches. For the full unprioritized menu see
`UPGRADE_BACKLOG.md`; for the raw numbers see `COMPONENT2_RESULTS.md`.

---

## PART 1 — What the first run actually established (synthesis)

1. **Grounding is real and seed-robust, but hyperparameter-DOMINATED.** Mean 0.195 (sd
   0.039) across 5 seeds at latent2048/k32; every seed far above the ~0 shuffle null. But
   the sweep shows the rate spans **0.136-0.395 from hyperparameters alone** (latent 512
   -> 0.395, 2048 -> 0.136 at fixed k), and *smaller dictionaries ground more*. The
   absolute rate is therefore NOT a quotable quantity. Consequence: the head-to-head must
   be a matched (latent, k) GRID, never a single number.
2. **The grounded set is an essential-gene artifact, quantified.** 8/26 seed-stable hits
   are Hart CEGv2 core-essential (~10x enrichment over a 3.0% base rate, hypergeometric
   p~4e-7), and CEGv2 is the STRICT 684-gene core -- the true housekeeping fraction is
   higher. Mechanism: essential-KD collapses transcription globally, dropping a shared
   'cell-health' feature that many knockdowns drop; suppression-vs-control alone cannot
   tell that from real regulation.
3. **The metric CAN find a real regulator.** GATA1 -- the master erythroid TF -- grounds
   in >=3/5 seeds in an erythroleukemia line, and sits in the NON-essential bucket. So the
   regulatory signal is present but buried, not absent.
4. **Two fixes are now built + synthetic-validated (this session):**
   - S.2 column specificity (metric v4, spec `029ab066`): reject features suppressed by
     MANY perturbations; synthetic proves it kills 5/5 shared stressors while keeping true
     regulators.
   - S.3 essential handling (`--exclude-list`, `--tag-list`): score non-essential TFs
     only, or split the rate by essential membership in one run.
   Neither has been run on real data yet -- that is step 1 below.

**The one-sentence state:** we have a real but confounded signal, a metric now hardened
against the confound, and a clean positive control (GATA1) -- but no valid head-to-head
number yet.

---

## PART 2 — The immediate re-run (do FIRST next GPU session)

Goal: does grounding survive the v4 column-specificity gate + essential removal, and does
GATA1 survive? This tells us whether there is ANY real regulatory grounding to compare.

```bash
git pull origin claude/code-handoff-doc-2v3jt3   # get metric v4 + gene sets
# A. v4 metric, ALL TFs, TAG essentials -> see the essential vs non-essential split
python3 src/causal_pipeline.py --adata replogle.h5ad --pert-col gene \
  --control-value non-targeting --rep expression --latent 2048 --k 32 --seed 0 \
  --n-hvg 2000 --tf-list config/tf_lists/hs_hgnc_tfs.txt \
  --tag-list config/gene_sets/hart_cegv2_core_essential.txt \
  --n-shuffle 50 --out causal_out/v4_tagged_seed0.json
# B. v4 metric, NON-ESSENTIAL TFs only, multi-seed
for s in 0 1 2 3 4; do
  python3 src/causal_pipeline.py --adata replogle.h5ad --pert-col gene \
    --control-value non-targeting --rep expression --latent 2048 --k 32 --seed $s \
    --n-hvg 2000 --tf-list config/tf_lists/hs_hgnc_tfs.txt \
    --exclude-list config/gene_sets/hart_cegv2_core_essential.txt \
    --n-shuffle 50 --out causal_out/v4_noness_seed$s.json
done
```

Decision branches:
- **If the tag breakdown shows in-essential >> non-essential under v4** -> column
  specificity did NOT fully remove the confound; escalate the gate (lower colspec_alpha,
  or add a program-coherence filter) before trusting any rate.
- **If v4 collapses the rate toward the null AND GATA1/other lineage TFs disappear** ->
  the entire signal was the confound; the honest finding becomes "expression-space SAE
  grounding on this data is essentially all essential-gene collapse" -- still publishable,
  reframes the story toward the confound as the result.
- **If v4 keeps a meaningful non-essential, GATA1-containing set above the null** -> THAT
  is the real grounded set; proceed to the head-to-head (Part 3). Most likely and most
  desirable outcome.

Also run the negative-control check (should ground ~0):
```bash
# nonessential-reference genes are not TFs; expect near-zero grounding (sanity)
python3 src/causal_pipeline.py ... --tf-list config/gene_sets/hart_negv1_nonessential.txt ...
```

---

## PART 3 — The head-to-head (the deliverable), redesigned as a grid

The results forced a redesign: because the rate is hyperparameter-dominated, the
comparison is a PAIRED GRID, identical (latent, k) for both arms.

- **Arm B (expression):** already runnable; run the (latent, k) grid {512,1024,2048} x
  {16,32,64} at v4, non-essential TFs.
- **Arm A (embedding):** scGPT AND Geneformer embeddings of the SAME control cells; train
  the SAME TopK SAE at the SAME grid; run the SAME v4 metric. This is the hands-on setup
  (checkpoint download, `embed_data`), the one real engineering lift left.
- **Comparison statistic:** at each grid point, expression_rate - embedding_rate; report
  the paired difference across the grid with a bootstrap CI. Positive result iff the CI
  excludes 0 across the grid (not a single lucky point). Also report the CURVE shape.
- Reproduce Arm A ~6-10% only as a pipeline sanity check, NOT as the comparison baseline
  (their hyperparameters differ; only the internal matched comparison is valid).

---

## PART 4 — The exhaustive direction map (by research question)

Every direction worth pursuing, grouped by the QUESTION it answers, re-prioritized by
tonight's results. Priority: [P0]=critical path, [P1]=high, [P2]=valuable, [P3]=optional.

### Q-A. Is there real (non-confound) regulatory grounding?  [P0]
- A1 [P0] The Part-2 re-run (v4 + essential split). The gate to everything else.
- A2 [P1] Program-coherence filter (S.2b): require the matched feature's decoder-top genes
  to be a concentrated, GO-enriched module, not a diffuse global axis. A second, orthogonal
  defense against the confound if column specificity is insufficient.
- A3 [P1] Stricter TF panel (S.1): fetch the Lambert 2018 sequence-specific-DBD list
  (~1639) to replace the pySCENIC 1839 that includes basal machinery (TBP/TAFs/GTFs).
  Host was egress-blocked here; obtain and commit as
  `config/tf_lists/lambert2018_hs_seqspecific_tfs.txt`.
- A4 [P2] Regress out per-cell global transcriptional level (library complexity) before
  the suppression test, so 'the cell is globally sick' cannot register as feature-specific
  suppression. The most direct attack on the confound mechanism; design carefully.

### Q-B. Does expression beat foundation-model embeddings?  [P0 deliverable]
- B1 [P0] Arm A scGPT at the matched grid (Part 3).
- B2 [P1] Arm A Geneformer too (beat BOTH -> airtight).
- B3 [P1] Paired bootstrap CI on the grid-wise difference; >=3/5 seeds each.

### Q-C. Where does causal legibility LIVE? (the flagship reframe)  [P1]
- C1 [P1] The representation-ladder curve: raw -> HVG -> PCA -> NMF -> VAE -> scGPT ->
  Geneformer, grounding rate along it. If monotone-decaying, the headline becomes a LAW,
  and any monotone outcome is a result (de-risks a weak expression-vs-embedding gap).
- C2 [P2] SAE-vs-PCA/NMF at matched dimensionality within expression space: does the SAE
  add over plain matrix factorization? (The first ML-reviewer question.)

### Q-D. Is the grounding biologically real?  [P1]
- D1 [P1] GATA1 deep-dive: is its matched feature the erythroid/heme module (HBB, HBA,
  ALAS2, GYPA, KLF1)? A hand-verifiable positive control anchors the whole result.
- D2 [P1] Blinded GO/pathway enrichment on every grounded feature's program.
- D3 [P2] Cross-reference grounded TF->program pairs against TRRUST/SCENIC regulons
  (code already exists from Gate 0) -- orthogonal validation.
- D4 [P2] Positive-control panel: GATA1, TAL1, KLF1 (erythroid), SPI1, CEBPA (myeloid) --
  do the known master regulators ground, and to the expected modules?

### Q-E. Is it robust?  [P1]
- E1 [P1] The (latent, k) grid at v4 (also feeds Part 3).
- E2 [P1] RPE1 second cell line -- causal analogue of the mouse+human Gate-0 replication.
- E3 [P2] Matching-rule robustness: program-DE alignment vs decoder-top-gene cosine vs MI;
  result should not hinge on the operationalization.
- E4 [P2] Cross-fit K-repeat (average over K A|B splits) so nothing rides on one split.
- E5 [P2] Threshold sensitivity surface for auc_floor / colspec_alpha / match_alpha.

### Q-F. Can the programs PREDICT unseen perturbations?  [P2]
- F1 [P2] Hold out a TF set entirely; predict their KD signatures from the SAE programs;
  compare expression vs embedding. Prediction > correlation = the strongest claim.

### Q-G. Can we build a causally-supervised SAE (method)?  [P3]
- G1 [P3] Take 2 (v1 failed, honestly shelved). New formulations: contrastive to
  perturbation directions with a held-out eval; perturbation-consistency regularizer;
  decoder-Jacobian alignment. Optional upside; never load-bearing.

### Q-H. Can we ship a resource / benchmark?  [P2]
- H1 [P2] Causally-grounded myeloid program atlas (web): each grounded program + genes +
  GO + causal certificate + pseudotime activity. Low-risk once a real grounded set exists.
- H2 [P3] Package the metric as a reusable benchmark ("CausalGround") anyone can run on any
  representation of any Perturb-seq dataset.

### Q-I. Is it positioned to win?  [P1]
- I1 [P1] Kendiukhov differentiation table (still not built; a judge WILL ask).
- I2 [P1] Verify CytoSAE (2507.12464) + scGPT-hematopoietic-manifold (2603.10261) for
  expression+causal collision.
- I3 [P2] Preprint (bioRxiv / MLCB / LMRL) before STS/ISEF.
- I4 [P2] The 'AI-interpretability auditor' narrative (PLI-Analyzer template); kill the
  v-number tangle in the writeup.

---

## PART 5 — Recommended sequence (what to actually do, in order)

1. [P0] Part-2 re-run (v4 + essential split + GATA1 check). One GPU session. Decides the
   story.
2. [P0/P1] If real grounding survives: Arm A scGPT at the matched grid (Part 3) -> the
   head-to-head number. Add Geneformer.
3. [P1] Biology anchor (GATA1 deep-dive + blinded GO) + the differentiation table +
   competitor verification (can be done off-GPU, in parallel).
4. [P1] The representation-ladder curve (C1) -- turns the head-to-head into the flagship
   'where causal legibility lives' result.
5. [P1] RPE1 replication + the (latent,k) grid CI.
6. [P2+] Atlas, prediction test, preprint; [P3] supervised-SAE take 2 only as upside.

**Fastest path to a defensible headline:** steps 1-4. If step-1 grounding survives and
step-2 shows expression > embedding on the matched grid with GATA1-anchored biology, that
is the winnable result the whole project has been driving toward -- reported honestly,
including the essential-gene confound we found and fixed, which is itself a rigor story.
