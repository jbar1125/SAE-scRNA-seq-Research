# COMPONENT 2 — NEXT DIRECTIONS (post-first-run analysis + execution map)

> **CURRENT STATUS, 2026-10-18:** the July execution map below is historical and
> superseded. The v4 re-run happened: column specificity retained housekeeping hits,
> reduced GATA1 to 2/5 seeds and was rolled back. The claimed ~10x enrichment used
> all 1,839 TFs rather than the 162 tested; the reported v4 rates were 6/51
> core-essential-tagged versus 14/111 other tested TFs, not a validated nonessential
> panel. M1 retained all 22 seed-0 hits, so the effect-size explanation was not
> supported. Replogle RPE1 and the PCA/NMF input ladder were subsequently run;
> foundation-model embeddings never ran. All causal SAE fits used **all cells**.
>
> PCA/NMF ladder rungs feed the factorization into an SAE, so they cannot establish
> SAE superiority or equivalence to direct PCA/NMF causal baselines. GATA1 was a
> Replogle hit under earlier settings, not a validated general positive control;
> it is absent from Norman. Historical causal raw JSONs are missing. Current next
> steps are archive recovery/regeneration, the program-modulation test with declared
> nulls and sensitivities, and bounded publication writing. See
> `PROJECT_MASTER_SUMMARY.md` sections 8-10 and 18. Historical commands below do not
> reconstruct v4 simply by naming the output file v4; current defaults differ.

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
2. **Housekeeping-heavy hits motivated a confound hypothesis.** The later tested-pool
   comparison did not support the claimed essential enrichment. Global-collapse
   explanations were hypotheses, not conclusions established by the invalid ~10x
   calculation or by the hit list's functional annotations.
3. **GATA1 was a v3 Replogle hit in >=3/5 seeds.** It fell to 2/5 under v4 and did
   not establish broad lineage-regulatory recovery. Being outside CEGv2 does not
   identify a perturbation as a validated nonessential control.
4. **Two fixes are now built + synthetic-validated (this session):**
   - S.2 column specificity (metric v4, spec `029ab066`): reject features suppressed by
     MANY perturbations; synthetic proves it kills 5/5 shared stressors while keeping true
     regulators in that synthetic test. It later failed on real data and was disabled.
   - S.3 essential handling (`--exclude-list`, `--tag-list`): score non-essential TFs
     only, or split the rate by essential membership in one run.
   At drafting these were synthetic-only. Their real-data outcomes are corrected above.

**Current interpretation:** reported perturbation response did not establish the
intended lineage-regulatory claim or a matched expression-vs-embedding comparison.

---

## PART 2 — Archived v4 re-run proposal (completed and superseded)

Historical question: would v4 remove housekeeping responses while retaining GATA1?
It did not. The command block and decision branches below record the proposal, not
the current next action or a supported interpretation of its failure.

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

### M1 (built 2026-07-18; subsequently tested with no hits removed)
`causal_grounding.py` now has an effect-size-control gate (spec v6, SHA `81b44e04`, default
on, `--effsize-alpha 1.0` disables) + an effect-size-vs-suppression diagnostic. Synthetic-
validated (crux: identical suppression at low vs high effect gets opposite calls). The open
empirical question was whether effect-matched comparison removed the reported hits.
It removed 0/22; rho = 0.28. The AND-filter cannot add or rescue a hit. The archived
commands below record that comparison, not a pending rescue experiment:
```bash
for a in 1.0 0.10; do   # 1.0 = raw, 0.10 = effect-size-controlled
  python3 src/causal_pipeline.py --adata replogle.h5ad --pert-col gene \
    --control-value non-targeting --rep expression --latent 2048 --k 32 --seed 0 \
    --n-hvg 2000 --tf-list config/tf_lists/hs_hgnc_tfs.txt --effsize-alpha $a \
    --n-shuffle 50 --out causal_out/m1_effsize_$a.json
done
```
The result did not support the proposed effect-size explanation. Do not represent
the historical expectation as a finding.

### Q-A. Is there real (non-confound) regulatory grounding?  [P0]
- A1 [completed/superseded] Part-2 rerun and M1. Do not iterate gates until a desired regulator appears.
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
- C1 [P1] The representation-ladder curve: expression -> PCA -> NMF -> (VAE) -> scGPT ->
  Geneformer, grounding rate along it. A decay was a hypothesis; a monotone curve
  from a selected grid would not by itself establish a general law.
  **The PCA/NMF rungs are BUILT: `--rep pca` / `--rep nmf` (+ `--rep-dim`) in the pipeline,
  synthetic-tested.** Run the expression + pca + nmf rungs now; add scGPT/Geneformer with
  Arm A:
  ```bash
  for r in expression pca nmf; do
    python3 src/causal_pipeline.py --adata replogle.h5ad --pert-col gene \
      --control-value non-targeting --rep $r --rep-dim 128 --latent 2048 --k 32 --seed 0 \
      --n-hvg 2000 --tf-list config/tf_lists/hs_hgnc_tfs.txt \
      --exclude-list config/gene_sets/hart_cegv2_core_essential.txt \
      --n-shuffle 50 --out causal_out/ladder_$r.json
  done
  ```
- C2 [P2, unimplemented] Direct SAE-vs-PCA/NMF causal baseline: does SAE training
  add over scoring factors directly? Existing `--rep pca/nmf` trains another SAE
  and does not answer this question.

### Q-D. Is the grounding biologically real?  [P1]
- D1 [historical hypothesis] Inspect the Replogle GATA1 feature and genes. This was
  not established as a robust positive control; Norman contains no GATA1 perturbation.
- D2 [P1] Blinded GO/pathway enrichment on every grounded feature's program.
- D3 [P2] Cross-reference grounded TF->program pairs against TRRUST/SCENIC regulons
  (code already exists from Gate 0) -- orthogonal validation.
- D4 [P2] Positive-control panel: GATA1, TAL1, KLF1 (erythroid), SPI1, CEBPA (myeloid) --
  do the known master regulators ground, and to the expected modules?

### Q-E. Is it robust?  [P1]
- E1 [P1] The (latent, k) grid at v4 (also feeds Part 3).
- E2 [subsequently run] RPE1 second cell line: reported housekeeping/chromatin hits,
  not lineage-regulatory recovery. The complete historical raw outputs need archiving.
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
- I1 [built, corrected 2026-10-18] Kendiukhov differentiation table in `DIFFERENTIATION.md`.
- I2 [P1] Verify CytoSAE (2507.12464) + scGPT-hematopoietic-manifold (2603.10261) for
  expression+causal collision.
- I3 [P2] Preprint (bioRxiv / MLCB / LMRL) before STS/ISEF.
- I4 [P2] The 'AI-interpretability auditor' narrative (PLI-Analyzer template); kill the
  v-number tangle in the writeup.

---

## PART 5 — Current sequence (2026-10-18)

1. Recover or regenerate the missing causal raw outputs with data/code/spec hashes,
   training metadata, seeds and complete null results.
2. Run the declared program-modulation experiment and report every seed, background
   null, specificity control and training-mode sensitivity, regardless of outcome.
3. Write the structural audit and causal limitations in parallel. Verify citations
   from full texts before asserting novelty or repeating competitor percentages.
4. Treat matched embeddings, direct-factor causal baselines, independent replication
   and held-out training as separately designed future experiments. Do not relabel
   the existing input ladder as any of those results.
