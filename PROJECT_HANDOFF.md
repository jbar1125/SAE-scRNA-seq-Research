# PROJECT HANDOFF: SAE Features in Hematopoietic Lineage Commitment

> **AUDIT NOTE:** This is the Phase-1-closeout snapshot, left intact as a historical record. For verified corrections, the V0b notebook trap, and rigor upgrades, read [`PROJECT_AUDIT.md`](./PROJECT_AUDIT.md). Where the two conflict, the audit is authoritative.

**Handoff target:** Claude Code (repo-based agentic work)
**Prepared:** end of Phase 1, before summer Phase 2
**Researcher:** Jacob Barzideh (jbarz), high school, Jericho Science Research
**Competition targets:** ISEF 2027, Regeneron STS 2027
**Mentor:** Mrs. Han (faculty), remote affiliation with Cold Spring Harbor Laboratory (mentored, not employed)
**Project type:** Fully computational. No wet bench (mentor confirmed student cannot do independent wet-bench work).

---

## 0. HOW TO USE THIS DOCUMENT

This is the single source of truth for the project state at Phase 1 closeout. It is exhaustive on purpose. Sections:

1. One-paragraph thesis
2. Scientific narrative and the headline finding
3. Current status board (done / in progress / blocked)
4. Phase 1: complete technical specification
5. The V0-V15 validation plan, item by item
6. The two active bugs (V14, V15) in full detail
7. Key quantitative findings to date
8. Methodological decisions and hard-won nuances (READ THIS to avoid repeating mistakes)
9. Phase 2 research plan (full structure, 5 parts, 4 gates)
10. ISEF / STS paperwork status
11. Every document version produced and why
12. File inventory
13. Verified references with accessions
14. Immediate next actions
15. Researcher working style and communication preferences

---

## 1. THESIS

This project applies sparse autoencoders (SAEs) from mechanistic interpretability to single-cell RNA sequencing (scRNA-seq) data to discover interpretable, sparsely-activated gene programs in hematopoiesis, then tests whether those programs are causally grounded (via CRISPRi perturbation), clinically meaningful (via AML and tumor bone-marrow data), and representation-invariant (via foundation-model embedding-space SAEs). The central scientific claim is "asymmetric modularity": granulocyte lineage commitment is captured as a single unified program while erythroid commitment is distributed across multiple interacting sub-modules.

---

## 2. SCIENTIFIC NARRATIVE AND HEADLINE FINDING

### The gap
Hematopoiesis is a continuous, branching process. Standard clustering and dimensionality reduction (PCA, scVI, Leiden) group cells but do not cleanly recover the overlapping gene programs running inside each cell because those programs are in "superposition." SAEs, via L1-regularized sparse dictionary learning, are designed to recover monosemantic features (one feature = one concept). The question is whether SAEs trained on scRNA-seq recover biologically meaningful, causally grounded programs.

### The finding (asymmetric modularity)
When SAE features are mapped onto developmental trajectories, the two major myeloid lineages decompose differently:
- **Granulocyte commitment**: a single unified module (mean z = -5.53 in mouse, -10.68 in human). One dominant feature carries the program.
- **Erythroid commitment**: distributed hierarchically across multiple interacting sub-features, with sub-additive double-knockout effects.

This asymmetry replicates across species (mouse Paul15 -> human CELLxGENE bone marrow) with stronger effect sizes in human (3/3 human seeds).

### Why it might matter
If asymmetric modularity is real and not a representation artifact, it says something about how the regulatory architecture of these two lineages differs. The Phase 2 plan tests causal grounding (Replogle CRISPRi), disease relevance (van Galen AML, Hegde tumor BM), and representation invariance (scGPT/Geneformer embedding-space SAEs).

### The honest caveat baked into the whole project
Identifying a biologically plausible feature does NOT prove causal regulatory logic. Kendiukhov (2025) showed only 6.2% of TF-associated features in foundation-model SAEs respond specifically to CRISPRi knockdown of that TF. This 6.2% is the pre-registered null hypothesis for Phase 2 Part 2. The project is built to report null results honestly if the features fail causal validation.

---

## 3. CURRENT STATUS BOARD

### DONE
- Mouse SAE training (5 seeds, 128 latents). Reproducible 5/5.
- Human SAE training (3 seeds, 128 latents, 25k CELLxGENE BM subsample).
- NMF baseline (5 seeds, rank 128).
- Preprocessing MD5-verified and byte-reproducible.
- V0a (Palantir pseudotime root fix): Spearman 0.733, 1Ery self-prob 0.90. 4/4 pre-registered criteria passed.
- Participation ratio: SAE and PCA computed.
- Phase 2 research plan: at v5 (current), all structural directives applied.
- ISEF paperwork: synopsis, all-forms doc, addendums doc all drafted.
- ICML 2026 mech interp workshop paper submitted (team paper, "Semantic Drift Recovery (SDR) Pipeline"). Jacob wrote part of intro + results. Decision ~6/12.

### IN PROGRESS
- V0b (gene-content module definitions): notebook at v2, 30 cells, delivered. Must be run to completion on Colab and the module_assignments_v0b.csv SHA-256 frozen. **This is the immediate blocker for locking the asymmetric modularity finding.**
- V13 cross-species human concordance: partially done (human SAE trained; concordance test not run).

### BLOCKED / NOT STARTED
- V14: OOD threshold bug (computes as 0.000). NOT FIXED.
- V15: perturbation delta scale bug (deltas ~0.001-0.005). NOT FIXED.
- V0d: non-circular perturbation pipeline. NOT BUILT.
- V0c: statistical bimodality (Hartigan dip test). NOT BUILT.
- V11 NMF participation ratio. NOT COMPUTED.
- V12 NMF rigorous modularity comparison. NOT DONE.
- V9 trajectory shape (logistic fits to replace "0.37 crossover"). NOT DONE.

---

## 4. PHASE 1: COMPLETE TECHNICAL SPECIFICATION

### 4.1 Mouse dataset (Paul15)
- Source: `scanpy.datasets.paul15()`
- Shape after preprocessing: **2,730 cells x 2,000 HVGs**
- 19 myeloid progenitor clusters (paul15_clusters): 1Ery-6Ery, 7MEP, 8Mk, 9GMP, 10GMP, 11DC, 12Baso, 13Baso, 14Mo, 15Mo, 16Neu, 17Neu, 18Eos, 19Lymph
- **CRITICAL: Paul15 data is pre-log-transformed.** Do NOT run normalize_total or log1p. Doing so double-transforms.
- HVG selection: variance-based (NOT seurat_v3, which errors on this data). HVG selection happens BEFORE scaling.
- Canonical marker genes force-included in HVG set.
- Expression matrix is scaled (zero-mean, unit-variance) after HVG selection.
- **MD5 of preprocessed expression matrix: `60183a17983c8b977d036e0f3a58da61`** (verify this hash at the start of every notebook).

### 4.2 Human dataset (CELLxGENE Census)
- Source: CZ CELLxGENE Discover Census, stable release (CZI 2023).
- Pool: 759,000 healthy human bone marrow cells (tissue_general = bone marrow, is_primary_data = True, disease = normal, hematopoietic lineage only).
- Training subsample: **25,000 cells, stratified.**
- **CRITICAL: CELLxGENE pipeline REQUIRES normalize_total + log1p before HVG selection** (opposite of Paul15, which is already log-transformed). Do not conflate the two pipelines.
- **OUTSTANDING: pin the exact Census release date used for training.** It is recorded in the training notebook. This must go into the Phase 1 reproducibility doc before Phase 2.

### 4.3 SAE architecture and hyperparameters
```
class SparseAutoencoder(nn.Module):
    encoder = nn.Linear(input_dim, latent_dim, bias=True)   # 2000 -> 128
    decoder = nn.Linear(latent_dim, input_dim, bias=True)   # 128 -> 2000
    forward: z = ReLU(encoder(x)); return decoder(z), z
```
- Latent dim: **128**
- Optimizer: **AdamW, lr = 5e-4**
- L1 penalty: **lambda_L1 = 0.1**
- Epochs: **300**
- Decoder columns unit-normalized (normalize_decoder) at every step.
- Dead-neuron resampling every 500 steps.
- Batch size 256, weight_decay 1e-5, gradient clip norm 1.0, cosine annealing with 5-epoch warmup (these last specifics are the Phase 2 embedding-SAE config; mouse Phase 1 used the core AdamW/lr/lambda/epochs above).
- Mouse: 5 seeds (0-4). Human: 3 seeds (0-2).

### 4.4 Reconstruction quality (mouse, per seed)
| Seed | Recon MSE | Mean active features/cell |
|------|-----------|---------------------------|
| 0 | 0.6157 | 46.4 |
| 1 | 0.6076 | 45.0 |
| 2 | 0.6143 | 46.2 |
| 3 | 0.6224 | 48.0 |
| 4 | 0.6186 | 46.9 |
Mean recon MSE 0.6157, std 0.0055. L0 sparsity ~45-48 active features per cell. This 20-50 range is the QC gate for all Phase 2 projections.

### 4.5 NMF baseline
- Vanilla sklearn NMF, NNDSVD init, rank 128, 5 seeds.
- Files: nmf_seed0-4.npz, nmf_stability_results.csv.
- Participation ratio NOT yet computed (V11 task).

### 4.6 V0a Palantir pseudotime (DONE, replaces broken DPT)
- **Why V0a exists:** original DPT used iroot=0, which is cell W31105, a 7MEP (already committed erythroid), not a progenitor. Every downstream pseudotime number was rooted inside the erythroid lineage. This was a real bug.
- Fix: Palantir with root = centroid of 10GMP cluster, cell **W39126 (idx 2690)**.
- 7 manual terminal states:
  - W31534 -> 8Mk
  - W31694 -> 19Lymph
  - W37539 -> 1Ery
  - W37872 -> 14Mo
  - W38318 -> 16Neu
  - W38981 -> 13Baso
  - W39151 -> 11DC
- Pseudotime range [0.0, 1.0]; entropy range [0.2768, 1.7467].
- **Spearman vs canonical cluster ordering: 0.733.**
- 1Ery branch self-probability: 0.90 (0.90/0.86/0.84/0.83 across 1Ery-4Ery).
- Known anomalies (documented, contained): 7MEP -> 19Lymph at 0.47 (backwards); 6Ery has low 1Ery prob (0.27).
- Output file: `cell_metadata_palantir.csv` with columns: cell_id, paul15_clusters, palantir_pseudotime, palantir_entropy, and 7 branch-prob columns suffixed with terminal cell IDs (prob_1Ery_W37539, prob_14Mo_W37872, prob_16Neu_W38318, prob_13Baso_W38981, prob_11DC_W39151, prob_8Mk_W31534, prob_19Lymph_W31694).

### 4.7 V0b gene-content module definitions (v2 notebook delivered, MUST BE RUN)
- **Why V0b exists:** the original "module" definitions (SAE_rebuilt.ipynb Cell 36) defined erythroid vs granulocyte modules purely by pseudotime correlation SIGN, not by gene content. Combined with the V0a root bug, the "early" features were not progenitor features. The asymmetric modularity claim rested on this. V0b re-derives modules from canonical marker gene content.
- Notebook: `v0b_module_definitions.ipynb`, 30 cells, v2.
- Approach: per (seed, feature, submodule) hypergeometric enrichment on top-30 decoder genes, BH FDR-corrected. Per-seed winner-take-all assignment (NOT cross-seed feature-index matching, because SAEs do not preserve feature identity across seeds). Cells filtered into committed groups by Palantir branch probability.

#### V0b v2 marker sets (HVG-present only)
Paul15 variance HVG selection filtered out hemoglobins, late granulocyte markers (Ltf, Lcn2, S100a8/9, Mmp8/9, Camp, Ngp), MK granule proteins, and most monocyte/baso/DC/lymph terminal markers. Marker sets restricted to what is actually in the 2000-HVG set:
- Ery_TF (8): Gata1, Klf1, Tal1, Lmo2, Zfpm1, Stat5a, Bcl11a, Myb
- Ery_Heme (6): Fech, Hmbs, Ppox, Cpox, Urod, Alad
- Ery_Membrane (6): Gypc, Ank1, Rhag, Aqp1, Epor, Tspo
- Gran_TF (3): Cebpa, Cebpe, Runx1
- Gran_Primary (4): Mpo, Elane, Prtn3, Ctsg
- Progenitor (4, control): Kit, Cd34, Mllt3, Eif4ebp1
- Cycling (16, control): Top2a, Pcna, Mcm2-7, Cdk1/4/6, Cenpe, Cenpf, Aurkb, Birc5, Ccnb2

#### V0b v2 key constants
- TOP_K = 30 (top decoder genes per feature)
- ENRICHMENT_Q_THRESHOLD = 0.05
- WINNER_TAKE_ALL_MIN_OVERLAP = 2 (require >= 2 marker genes in top-30)
- BRANCH_PROB_COMMITTED = 0.7
- BRANCH_PROB_UNCOMMITTED = 0.5
- N_BOOTSTRAP = 1000
- SEED = 42

#### V0b v2 cell-group filtering (committed lineage)
- committed_erythroid: prob_1Ery > 0.7
- committed_granulocyte: (prob_14Mo + prob_16Neu + prob_13Baso) > 0.7
- uncommitted_progenitor: max over ALL 7 terminal probs < 0.5
- intermediate: everything else
- **CRITICAL FIX in v2:** max_branch_prob MUST include prob_19Lymph. v1 omitted it, which mislabeled ~273 high-Lymph cells as uncommitted_progenitor.

#### V0b v2 predicted cell counts (run on actual data)
- committed_erythroid: 751 (mostly 1Ery-5Ery)
- committed_granulocyte: 897 (mostly 13Baso, 14Mo, 15Mo, 16Neu)
- uncommitted_progenitor: ~700 (v2 with 19Lymph fix; v1 wrongly gave 973)
- intermediate: ~382 (v1 wrongly gave 109)
- committed_erythroid and committed_granulocyte counts are stable across v1/v2 (determined by prob_1Ery and prob_gran_combined, not max_branch_prob).

#### V0b v2 asymmetric modularity test (redesigned)
v1 used Shannon entropy ratio (ery entropy / gran entropy), but erythroid has 3 submodules and granulocyte has 2, biasing the ratio. v2 uses:
- **Normalized entropy** (entropy / log(n_submodules)) in [0,1]: 1.0 = fully distributed, 0.0 = unified.
- **Largest-submodule fraction**: max(counts)/sum(counts). 1.0 = unified.
Both computed per seed. Bootstrap CI across 5 seeds.
Asymmetric modularity is SUPPORTED if BOTH hold in >= 3 of 5 seeds AND both bootstrap CIs exclude 0 in the same direction:
- ery_norm_entropy > gran_norm_entropy (erythroid more distributed)
- ery_largest_frac < gran_largest_frac (erythroid less concentrated)
If neither holds: the correlation-sign-based finding does not survive the rigorous gene-content test, and the claim must be softened. **This is a real possibility and must be reported honestly.**

#### V0b output files
module_assignments_v0b.csv (the pre-registration table to SHA-256 freeze), enrichment_full_v0b.csv, modularity_metrics_per_seed.csv, submodule_dynamics.csv, cell_metadata_v0b.csv, v0b_provenance.json.

#### V0b run instructions (KNOWN GOTCHA)
The user repeatedly ran v1 instead of v2. To verify v2 before running:
- Cell 4 MARKER_SETS starts with `'Ery_TF': ['Gata1', 'Klf1', 'Tal1', 'Lmo2', 'Zfpm1', 'Stat5a', 'Bcl11a', 'Myb']` (NOT hemoglobins like Hba-a1).
- Cell 1 has `WINNER_TAKE_ALL_MIN_OVERLAP = 2` (NOT `CROSS_SEED_THRESHOLD = 3`).
- Cell 9 terminal_cols INCLUDES `prob_19Lymph`.
File goes in Colab Drive at the path matching the loader (BASE = `/content/drive/MyDrive/SAE_phase1`) alongside expression_matrix.npy, gene_names.csv, cell_metadata_palantir.csv, sae_seed0-4.pt. Run top to bottom from a clean runtime.

---

## 5. THE V0-V15 VALIDATION PLAN (ITEM BY ITEM)

### Foundation fixes (precede everything)
- **V0a** Pseudotime root fix. DONE.
- **V0b** Gene-content modules. v2 notebook delivered; MUST RUN + freeze SHA-256.
- **V0c** Statistical bimodality (Hartigan dip test + BH FDR, replaces sigma=3 peak-finding). NOT BUILT. Small notebook (~3 hrs).
- **V0d** Non-circular perturbation (HVG markers scored on perturbed reconstruction without re-encoding). NOT BUILT. Tied to V14/V15.

### Foundation validation
- **V1** Data integrity. MOSTLY DONE (MD5 verified; cluster cross-check vs Paul15 Fig 1 pending).
- **V2** Preprocessing reproducibility. DONE (byte-identical across fresh runs).
- **V3** SAE training reproducibility. PARTIAL (5/5 reproducibility passed; formal single-seed identical-loss-curve check not run).

### Feature recovery validation
- **V4** Marker recovery quantification (10 TFs x 5 seeds x 4 thresholds K=10/20/30/50 = 200-entry matrix). PARTIAL (V0b covers K=30 only; threshold sweep pending).
- **V5** Feature stability across seeds (Hungarian/cosine matching, per-feature matching-score distribution). NOT DONE. Distinct from V0b.
- **V6** Bimodality = V0c.

### Pseudotime validation
- **V7** Root cell sensitivity (10 roots from progenitor cluster). PARTIAL (9GMP vs 10GMP correlation 0.9931; 8 more roots pending).
- **V8** Pseudotime method comparison (DPT, Palantir, scFates, CellRank). PARTIAL (Palantir + CellRank done; CellRank GPCCA flagged 7MEP as terminal incorrectly and was uninformative; scFates not run).
- **V9** Trajectory shape with logistic fits (midpoint, slope, sharpness per lineage). NOT DONE. **This replaces the "0.37 crossover" headline with the real quantitative finding.** The crossover position is a derived quantity, not the headline.

### Structural finding validation
- **V10** Modularity ratio robustness. DONE inside V0b v2 (4 metrics + bootstrap).
- **V11** Participation ratio (SAE/PCA/NMF). PARTIAL. SAE = 24.2, PCA = 43.4 (memory) / 42.3 (one JSR doc - reconcile this discrepancy). NMF pending. Formal bootstrap pending.
- **V12** NMF rigorous comparison (same feature-interpretation pipeline on NMF components, modularity ratio for NMF). PARTIAL.

### Cross-species
- **V13** Human replication audit (re-run feature recovery on 25k subset, V0b-style asymmetry in human, cross-species concordance vs mouse). PARTIAL.

### Perturbation bugs (blocking, not validation)
- **V14** OOD threshold bug. NOT FIXED. See section 6.
- **V15** Perturbation delta scale bug. NOT FIXED. See section 6.

### Priority order for Phase 1 closeout
1. Run V0b v2 + SHA-256 freeze.
2. Pin CELLxGENE Census date in reproducibility doc.
3. V0c (bimodality).
4. V14 + V15 + V0d (the hard perturbation sprint).
5. V11 NMF + V12.
6. V9 trajectory shape.
7. V13 cross-species concordance.
Items V5, V7 (full 10-root), V8 scFates, V4 threshold sweep, V1 cluster cross-check can spill into Phase 2 buffer.

---

## 6. THE TWO ACTIVE BUGS (FULL DETAIL)

### V14: OOD threshold computing as 0.000
- **Symptom:** the PCA-based out-of-distribution gate flags ALL cells as OOD because the threshold computes to 0.000.
- **Root cause:** incorrect self-distance exclusion in the nearest-neighbor null distribution. When building the null distribution of nearest-neighbor distances among real cells, each cell's nearest neighbor is itself (distance 0) because self is not excluded. The 95th percentile of a distribution full of zeros is ~0.
- **Fix:** exclude self when computing each cell's nearest-neighbor distance (k-NN with k starting at the second neighbor, or mask the diagonal of the distance matrix).
- **Verification:** re-run; threshold should be a positive number. Test on synthetic data where OOD membership is known (inject clearly-out-of-distribution points and confirm only those flag).

### V15: Perturbation delta values extremely small (~0.001-0.005)
- **Symptom:** knocking out a gene and re-scoring produces near-zero probability shifts.
- **Root cause:** scale mismatch between the decoder output and the X_train used to fit the downstream classifier. The decoder reconstructs in one scale; the classifier was fit on X_train in another scale (likely the scaled expression matrix). Feeding decoder output to a classifier expecting scaled input produces tiny deltas.
- **Fix:** trace the exact preprocessing chain. Identify where decoder output and classifier-input scales diverge. Apply the same scaling transform to decoder output that was applied to X_train before classification.

### V0d: the deeper circularity problem (must accompany V14/V15)
- The current perturbation pipeline decodes a perturbed latent, RE-ENCODES through the SAME SAE, and scores with the SAME features. Output is bounded by what the SAE can represent. This is methodologically circular.
- **Non-circular fix:** score perturbations using independent marker genes that exist in the HVG set (Klf1, Gata1, Cebpa, Cebpe, Elane, Mpo), measured DIRECTLY on the perturbed reconstruction, without re-encoding.
- The classifier-based framework was ALSO identified as circular before implementation; marker-gene scoring was adopted instead.
- **Reporting nuance:** perturbation analysis shows FEATURE SPECIFICITY, not causal bistability. Do not overclaim. The Replogle CRISPRi work in Phase 2 is the real causal test.

### Three-step perturbation methodology (target design)
1. Independent marker-gene scoring (Klf1, Gata1, Cebpa, Cebpe, Elane, Mpo).
2. Logistic regression classifier on paul15_clusters (~0.847 cross-validated accuracy) — use with care given circularity caveat.
3. PCA-based OOD gate (flags reconstructions where >10% of cells fall outside the 95th percentile of real-cell nearest-neighbor distances) — fix V14 first.

---

## 7. KEY QUANTITATIVE FINDINGS TO DATE

- Pseudotime correlations up to r = -0.828 (Feature 73, which contains Id1).
- Unsupervised recovery of canonical TFs and effectors: GATA1, TAL1, FECH, Aqp1, erythroferrone (Erfe), GFI1B, Id1.
- Progenitor-to-erythroid commitment crossover consistently at pseudotime ~0.37 (BEING REPLACED by V9 logistic-fit characterization; treat 0.37 as provisional).
- 49 bimodal features suggesting transition-reactivated gene programs (BEING REPLACED by V0c proper statistics; 49 is a heuristic count from sigma=3/height=max*0.25/distance=10 peak-finding, NOT a statistical result).
- Participation ratio: SAE 24.2 vs PCA 43.4. NMF reported 1.8x higher than SAE (~34.94). (Reconcile PCA 43.4 vs 42.3 discrepancy across documents.)
- Asymmetric modularity z-scores: granulocyte unified module mean z = -5.53 (mouse), -10.68 (human).
- Human generalization replicates with stronger effect sizes (3/3 human seeds).
- NMF does NOT recover the sub-additive erythroid hierarchy (1/5 NMF seeds vs 4/5 SAE seeds).
- CEBPA/CEBPE recovery is run-dependent (not stable across all seeds).

### Annotation table (generate_annotation_table.py)
- Generates feature-to-TF annotation table from SAE checkpoints.
- Demo run on mouse checkpoints: annotation_table_mouse.csv, 640 entries (128 features x 5 seeds), SHA-256 `bd5310926e2b0d735856cb0f4317ad8f401917267fdae71dcd254f6b4bd06373`.
- **Known bug:** GATA1/TAL1 get absorbed into GFI1B due to marker-set overlap and short-marker-set normalization bias. Marker refinement is a Phase 1 closeout task. (V0b's curated non-overlapping marker sets are the fix direction.)

---

## 8. METHODOLOGICAL DECISIONS AND HARD-WON NUANCES

**Read this section before touching the pipeline. These are mistakes already made and resolved.**

### Data
- Paul15 is PRE-LOG-TRANSFORMED. Skip normalize_total and log1p. (CELLxGENE is NOT; it needs both before HVG.)
- Use variance-based HVG selection for mouse. seurat_v3 ERRORS on Paul15.
- HVG selection BEFORE scaling.
- Force-include canonical marker genes in the HVG set or they get filtered out (this is why late-granulocyte and hemoglobin markers are still missing; the force-include list was not comprehensive).
- Save checkpoints immediately after good runs.
- np.ptp is deprecated in NumPy 2.0.
- Colab session resets cause NameError failures; code must survive re-execution (re-import, re-define, re-load at the top of each notebook).
- Dataset provenance must be traceable before inclusion in application materials.

### Methodological
- SAEs do NOT preserve feature identity across random initializations. Cross-seed logic must use per-seed winner-take-all, NOT cross-seed feature-index matching. (This broke V0b v1.)
- Perturbation analysis shows feature specificity, not causal bistability. Do not overclaim.
- Re-encoding perturbed expression through the same SAE is circular. Use independent marker-gene scoring on the reconstruction.
- Classifier-based perturbation framework is circular; marker-gene scoring adopted.
- The "0.37 crossover" and "49 bimodal features" are provisional/heuristic, not yet statistically grounded. V9 and V0c fix these.
- The original DPT root was a committed MEP, not a progenitor. Always root pseudotime at a genuine progenitor (10GMP).
- Modularity must be defined by GENE CONTENT, not pseudotime-correlation sign.

### Statistical conventions
- BH-FDR per test family.
- Effect sizes (Cohen's d) reported alongside p-values.
- Cross-seed validation: >= 2 of 3 human seeds, or >= 3 of 5 mouse seeds, for any reported finding.
- Pre-specify thresholds before data inspection (pre-registration via SHA-256 freezes).

---

## 9. PHASE 2 RESEARCH PLAN (FULL STRUCTURE)

**Current document: `phase2_research_plan_v5.docx` (260 paragraphs).**

Title: "Causal Validation and Clinical Translation of Sparse Autoencoder Features in Hematopoietic Lineage Commitment"

### Structure
- 5 thematic Parts (renamed from "Component" to "Part" globally), 12-week timeline, 4 decision gates.
- Random seed fixed at 42. Software versions pinned.
- ROL is a 5-paragraph funnel, 36 references, APA. Final ROL paragraph funnels into 5 strands mapping 1:1 to the 5 RQs and 5 objectives.
- RQs and Objectives are strictly parallel (1:1 mirrored language).

### Part 1: van Galen 2019 AML application (weeks 1-3)
- Clinical translation arm. Project human SAEs onto AML cells.
- Week 1: download GSE116256, QC, normalize 1e4 + log1p, HVG seurat_v3 n=2000, scVI integration (sample_id batch, n_latent=30), celltypist annotation (Immune_All_Low.pkl), restrict to 7 in-distribution bins.
- Week 2: project SAE, QC gate (L0 in 20-50, MSE within +50%), Wilcoxon AML vs healthy per (feature, bin) = 896 tests/seed, BH-FDR, Cohen's d, robustness >= 2/3 seeds.
- Week 3: patient pseudobulk, MixedLM (feature ~ disease_status + cell_type + (1|patient_id)), mutation stratification (NPM1, FLT3-ITD, IDH1, IDH2, CEBPA, TP53, RUNX1), Mann-Whitney U on subtypes with >= 3 patients, Cohen 1988 power analysis.
- Gate 1: PASS if >= 10 features q<0.05 & |d|>0.3 in >= 2 seeds AND >= 1 mutation distinction q<0.10.

### Part 2: Replogle CRISPRi causal validation (weeks 4-6)
- The most rigorous test. Null = Kendiukhov 6.2%.
- Week 4: download K562_gwps_normalized_singlecell_01.h5ad, OOD calibration (healthy BM vs K562 NT, flag if MSE shifts >50%), subset 24 TFs + NT 1:3, target 300k-500k cells.
- Week 5: project onto NT baseline + each KD, standardized effect d = (mean_KD - mean_NT)/pooled_std, top-k most-suppressed (k=1,3,5), check pre-registered annotated feature against top-k.
- Week 6: off-target full effect vector, permutation 10,000 shuffles, BH-FDR across 3,072 tests (24x128)/seed, causal validation table, one-sided binomial vs 6.2%.
- Gate 2: STRONG PASS >= 8/24 TFs at top-3, MODEST 3-7, MATCH NULL 0-2 (reframe as descriptive).

### Part 3: Hegde 2025 tumor BM projection (week 7)
- Mouse-to-mouse, no cross-species mapping.
- Download GSE270148, subset KP tumor + naive C57BL/6, exclude lung TME, scVI integration with Paul15, scANVI label transfer, project mouse SAE, Wilcoxon tumor vs naive, BH-FDR.
- Asymmetry test (RQ4): chi-squared 2x2 (erythroid suppressed vs gran elevated), Yates corrected. Concordance: Pearson within each feature set.

### Part 4: scGPT (+ optional Geneformer) embedding-space SAE (week 8)
- Representation invariance test. Hyperparameters matched to Phase 1.
- scGPT inference on Paul15 (HCOP ortholog mapping mouse->human), extract embeddings (n_cells, 512) from best-silhouette layer, Leiden ARI validation, train SAE on embeddings (matched hyperparams, 5 seeds), captum IntegratedGradients attribution, annotate vs canonical markers.
- **CCA appears in BOTH data prep AND analysis (v5 directive):** prep step = subset both matrices to matched 2,730 cells + column-wise standardize; analysis = sklearn CCA, first 5 canonical correlations + bootstrap CIs.
- Hungarian matching (scipy linear_sum_assignment), modularity ratio per architecture + 1,000-iteration bootstrap CIs.
- Gate 3: PASS if embedding-space asymmetry same direction as expression-space.

### Part 5 (optional): cross-system replication (week 9)
- Park 2020 human thymus (E-MTAB-8581, ~200k cells). Test whether asymmetric modularity generalizes to T-cell development.

### Buffer/writing (weeks 9-12)
- Week 9: Part 5 or gap-fill. Week 10: Methods + figures. Week 11: full manuscript + mentor review. Week 12: polish + supplementaries.
- Gate 4: PASS if 3+/4 primary parts cleared their gates.

### Statistical test families (for reference)
- Cell-bin family: 896 tests (128 features x 7 in-distribution bins).
- Replogle FDR family: 3,072 tests (24 TFs x 128 features) per seed.
- Permutation: 10,000 shuffles per (TF, feature) pair.
- Bootstrap: 1,000 iterations for modularity CIs and CCA CIs.

### 24 canonical hematopoietic TFs (Replogle subset)
GATA1, KLF1, TAL1, FLI1, CEBPA, CEBPE, GFI1, GFI1B, SPI1, RUNX1, LMO2, NFE2, MYB, BCL11A, IKZF1, MEF2C, HOXA9, MEIS1, ZBTB7A, IRF8, ETV6, EVI1, FOXO1, ID1.

### Phase 2 plan formatting conventions (v5)
- "Component" renamed to "Part" globally.
- Dataset names bold.
- Software names italic.
- Each dataset entry includes species label (Human/Mouse), publication year, size with units.
- scGPT explicitly labeled open source, MIT license.
- Geneformer explicitly labeled open source, Apache 2.0 license.
- Risk and Safety reduced to N/A with Form 3 cross-reference.
- Hotelling 1936 added to bibliography for CCA provenance.
- No preface throat-clearing sentences under section headers.

---

## 10. ISEF / STS PAPERWORK STATUS

Project is computational, no wet bench. Form determinations:
- **NEEDED:** Form 1 (Adult Sponsor), Form 1A (Student Checklist + Research Plan), Form 1B (Approval), Form 3 (Risk Assessment, school-policy override for computer use), Official Abstract.
- **DEFENSIVE:** Form 2 (Qualified Scientist) if mentor signs.
- **NOT NEEDED:** Forms 1C, 4, 5A, 5B, 6A, 6B, 7. (Form 7 pending local fair coordinator confirmation since Phase 1 was internal-only.)
- **Rules Wizard:** answer "None of the Above" on the 6 checkboxes.

### Addendums (all assessed)
- AI Usage Policy: RELEVANT. Jacob documented AI usage per ISEF policy with prompt logbook + citation. Future AI use restricted to acceptable categories (lit summarization, grammar refinement of own text, code refinement with citation, statistical method ID, bibliography formatting).
- Field Work Safety Plan: N/A (no field work).
- BSL-1 checklist: N/A (no microorganisms/cell culture).
- BSL-2 checklist: N/A (no biological agents).

### Documents produced
- isef_paperwork_synopsis.docx (144 paragraphs)
- isef_all_forms.docx (57 paragraphs, all 13 forms with who-signs/when/justification)
- isef_addendums.docx (19 paragraphs, 4 addendums with summaries)

### Subject-specific guideline determinations (in research plan)
- Human Participants: N/A (de-identified public data, not human participants research per ISEF).
- Vertebrate Animals: N/A (mouse data generated by others under Weizmann and Mount Sinai IACUC).
- Hazardous Biological Agents: N/A.
- Hazardous Chemicals/Activities/Devices: N/A (Form 3 covers computer/cloud use).

---

## 11. EVERY DOCUMENT VERSION PRODUCED (AND WHY)

### Phase 2 research plan lineage
- phase2_research_plan.docx (v1, 247 para) — first full draft.
- phase2_research_plan_v2.docx (299 para) — expansion.
- phase2_research_plan_v3.docx (279 para) — revision.
- phase2_research_plan_v4.docx (285 para) — CELLxGENE source documented.
- phase2_research_plan_conversational.docx (263 para) — conversational voice attempt, better section linking, Risk trimmed. (Superseded: moved ROL too far toward first person.)
- phase2_research_plan_final.docx (257 para) — restored formal ROL, ISEF norms matched, preface sentences cut, "Component" still present.
- **phase2_research_plan_v5.docx (260 para) — CURRENT.** ROL funneled into 5 strands; RQs/objectives strict 1:1; CCA in prep + analysis; species/year/units on datasets; scGPT MIT + Geneformer Apache 2.0 explicit; "Component"->"Part"; datasets bold; software italic; Risk N/A; Hotelling 1936 added.

### Notebooks and scripts
- v0b_module_definitions.ipynb (v2, 30 cells) — CURRENT. Gene-content modules. MUST RUN.
- generate_annotation_table.py — feature-to-TF annotation table generator. Has the GATA1/TAL1->GFI1B absorption bug.

### ISEF docs
- isef_paperwork_synopsis.docx, isef_all_forms.docx, isef_addendums.docx (see section 10).

---

## 12. FILE INVENTORY

### Inputs (Phase 1 working set)
- expression_matrix.npy — 2730 x 2000 float32, scaled, MD5 60183a17983c8b977d036e0f3a58da61
- cell_metadata.csv — cell_id, paul15_clusters, dpt_pseudotime (old DPT, superseded by V0a)
- cell_metadata_palantir.csv — V0a output: cell_id, paul15_clusters, palantir_pseudotime, palantir_entropy, 7 suffixed branch-prob columns
- gene_names.csv — 2000 mouse HVGs, MGI symbol case

### Checkpoints
- sae_seed0.pt ... sae_seed4.pt + best_model.pt (mouse SAE, 5 seeds)
- sae_human_seed0.pt ... sae_human_seed2.pt + best_model.pt (human SAE, 3 seeds)
- nmf_seed0.npz ... nmf_seed4.npz (NMF baseline)
- Drive paths: Drive/MyDrive/SAE_stability/sae_seed{0-4}.pt; Drive/MyDrive/SAE_human/sae_human_seed{0-2}.pt; Phase 1 working dir Drive/MyDrive/SAE_phase1/

### Results files
- nmf_stability_results.csv, human_stability_results.csv
- annotation_table_mouse.csv (640 entries, SHA bd5310...)

### Notebooks
- SAE_rebuilt.ipynb — canonical Phase 1 notebook (Cell 6 = old DPT root bug; Cell 36 = correlation-sign modules; Cell 38 = heuristic bimodality; Cell 44 = circular perturbation). Also SAE_all_in_one variants.
- v0b_module_definitions.ipynb — current V0b.

### Template/reference PDFs
- JSR_Research_Plan_Template.pdf, Copy_of_JSR_Research_Plan_Template.pdf (current plan as PDF), Barzideh__Jacob_-_JSOTR18.pdf (cycle reflection).

---

## 13. VERIFIED REFERENCES WITH ACCESSIONS

- Replogle et al. 2022, Cell 185(14):2559-2575.e28. K562 data: figshare+ DOI 10.25452/figshare.plus.20029387. 9,866 genes.
- van Galen et al. 2019, Cell 176(6):1265-1281.e24. GSE116256. 38,410 cells, 40 aspirates (16 AML + 5 healthy). Seq-Well. Processed Seurat object: petervangalen/reanalyze-aml2019.
- Hegde et al. 2025, Nature 646(8087):1214-1222. GSE270148. Code: github.com/Merad-Lab/Hegde_Myelopoiesis_Epigenetics. Miriam Merad senior author.
- Zeng et al. 2022, Nature Medicine 28:1212-1223. **Primary AML companion paper (READ THIS for Part 1 framing).** John Dick lab. Hierarchy + drug-response deconvolution of 1000+ patients into Primitive/Mature/GMP/Intermediate classes. Methodologically linked to van Galen.
- Secondary AML reads: Beneyto-Calabuig et al. 2023 (Cell Stem Cell, clonal multi-omics); Lasry et al. 2023 (Nature Cancer, immune microenvironment + risk). Petti et al. 2019 (Nat Comms) only if calling mutations from expression.
- Kendiukhov 2025, arXiv 2603.02952 (atlas, 6.2% TF validation rate, the null). Kendiukhov 2026a arXiv 2603.01752 (causal circuit tracing). Kendiukhov 2026b arXiv 2603.11940 (exhaustive circuit mapping).
- Cui scGPT 2024, Nature Methods 21(8):1470-1480. MIT license. Zenodo DOI 10.5281/zenodo.10466117. 12 layers, dim 512, ~33M cells.
- Theodoris Geneformer 2023, Nature 618(7965):616-624. Chen et al. 2024 V2-316M quantized (bioRxiv 2024.08.16.608180). Apache 2.0. 18 layers, dim 1152, ~95M cells.
- Setty Palantir 2019, Nature Biotech 37(4):451-460.
- Lange CellRank 2022, Nature Methods 19(2):159-170. Weiler CellRank 2 2024, Nature Methods.
- Paul 2015, Cell 163(7):1663-1677.
- CZI CELLxGENE Census 2023, bioRxiv 10.1101/2023.10.30.563174.
- Park 2020, Science 367(6480):eaay3224. E-MTAB-8581.
- Adrover 2025, Nature 645:484-495 (neutrophils, vascular occlusion, metastasis — wet-bench extension hypothesis link).
- LaMarche 2024, Nature 625:166-174 (IL-4 axis, pro-tumorigenic myelopoiesis).
- Bricken 2023, Cunningham 2023, Templeton 2024 (SAE interpretability foundations). Makhzani & Frey 2014 (k-sparse AE). Rajamanoharan 2024 (JumpReLU). Bussmann 2024 (BatchTopK), 2025 (Matryoshka).
- Benjamini & Hochberg 1995 (FDR). Cohen 1988 (power). Hotelling 1936 (CCA). Haghverdi 2016 (DPT). Wolf 2018 (scanpy). Lopez 2018 (scVI). Orkin & Zon 2008, Velten 2017 (hematopoiesis). Schuster 2025, Claye 2025 (SAEs for gene expression). Eraslan 2019, Ma & Xu 2022 (DL for genomics).

---

## 14. IMMEDIATE NEXT ACTIONS (FOR CLAUDE CODE)

In priority order:

1. **Run V0b v2 to completion.** Verify v2 (not v1) per section 4.7 gotcha. Confirm cell counts (~751/897/700/382). Read Cell 7 (per-seed module counts), Cell 9 (cell groups), Cell 12 (per-seed asymmetric modularity metrics), Cell 13 (bootstrap CIs). Determine whether asymmetric modularity holds under the gene-content test. SHA-256 freeze module_assignments_v0b.csv.
2. **Pin CELLxGENE Census release date** in a Phase 1 reproducibility document (requirements.txt + MD5 + census version).
3. **Reconcile the PCA participation ratio discrepancy** (43.4 vs 42.3 across documents). Recompute with one documented formula.
4. **Build V0c** (Hartigan dip test + BH FDR for bimodality). Report how many of the 49 features survive proper statistics.
5. **Fix V14** (OOD self-distance exclusion), **fix V15** (decoder/classifier scale match), **build V0d** (non-circular marker-gene perturbation scoring). Then run Klf1/Gata1/Cebpa/Cebpe knockouts with OOD gate active, report per-marker shifts + bootstrap CIs.
6. **Compute V11 NMF participation ratio** (rank 128, 5 seeds) + **V12** (NMF modularity via same pipeline). Output one table with bootstrap 95% CIs across SAE/PCA/NMF.
7. **Build V9** (logistic-fit trajectory shape) to replace the provisional 0.37 crossover with midpoint/slope/sharpness per lineage.
8. **Run V13 cross-species concordance**: V0b on 3 human seeds (HCOP-mapped markers, celltypist annotations as lineage proxy), compare human asymmetry direction vs mouse.
9. **Fix the annotation table marker absorption bug** (GATA1/TAL1 -> GFI1B) via curated non-overlapping marker sets; regenerate on all 8 seeds; freeze new SHA-256.
10. **Download van Galen GSE116256** and **read Zeng 2022** to prep Phase 2 Part 1.

### Compute constraints
- Free Colab (T4, ~12 GB RAM, 90-min idle timeout, 12-hr max session) for Phase 1. Pro+ expected for summer Phase 2.
- Paul15 (2,730 cells) work is trivial. 25k human subset fine. Full 759k human inference needs batching. Full Replogle (~2.5M cells) and Geneformer V2-316M need A100 (RunPod).
- All code must survive Colab kernel resets (re-import/re-define/re-load at top of each notebook).

---

## 15. RESEARCHER WORKING STYLE AND COMMUNICATION PREFERENCES

- Extremely concise responses. No fluff.
- **No em-dashes.** No AI-sounding language.
- Per-bullet formatting for structured answers.
- First person comfortable in operational contexts.
- Code delivered ready-to-run without lengthy explanation.
- Values honest methodological self-assessment and brutal constructive criticism over validation. ("BE AS HONEST AS POSSIBLE, CONSTRUCTIVE CRITICISM TO ALL HEAVENS.")
- ROL voice: formal third-person/passive academic English.
- Methodology voice: imperative future tense.
- Concrete next steps phrased as "I must...".

### Voice calibration for documents
- ISEF research plans use direct factual statements without preface throat-clearing (confirmed via research on real ISEF plans).
- Do not add filler sentences under section headers (e.g., no "This section describes..." or "Each subsection is required...").

### Collaboration context
- Co-submitted ICML 2026 mech interp workshop paper (SDR pipeline) with a team. Jacob wrote part of intro + results.
- Splunk Agentic Ops hackathon with collaborator Aaron Yang.
- Lab notebook: JSR-compliant bound carbon-copy format + Benchling/GitHub electronic.

### Possible future directions (not committed)
- Wet-bench extension (would trigger ISEF Forms 5b and 6b). Hypothesis: SAE commitment crossover point connects to tumor-induced myeloid skewing (Adrover 2025).
- Outreach to Mikala Egeblad at Cold Spring Harbor Laboratory for collaboration.
- Phase 3 (GSE270148 Hegde KP tumor BM) and Phase 4 already sketched as extensions.

---

*End of handoff. This document reflects project state at Phase 1 closeout. The single most important immediate action is running V0b v2 and determining whether the asymmetric modularity finding survives the gene-content test, because every downstream Phase 2 claim depends on it.*
