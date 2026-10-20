# Research Plan (v6)

> **CURRENT STATUS, 2026-10-18:** this is the historical Phase 2 proposal, with
> factual corrections below. It is not a record that its five components ran, nor
> proof that later experiments were pre-registered under this design. Gate 0 was
> executed and frozen using Paul15 mouse and Setty 2019 CD34+ human marrow (4,142
> cells), not the older 25,000-cell CELLxGENE sample. The structural result is
> method-dependence; the reversed lineage direction is not invariant or causally
> validated. Gate-0 dictionaries have 512 latents for 2,012/2,032 genes and are
> undercomplete relative to the gene input.
>
> Component 2 was replaced by cross-fit matching/Mann-Whitney tests per perturbation,
> not the original 24-TF annotated-feature/SD-drop design. Its causal SAEs trained
> on all cells, without perturbation labels, despite the earlier control-only intent.
> Replogle K562/RPE1 and Norman analyses ran; historical raw JSONs are missing.
> PCA/NMF ladder rungs were inputs to an SAE, not direct factor causal baselines.
> Components 1, 3, 4 and 5, including the embedding head-to-head, never ran.
> The cited 6.2% is literature context from a different metric/panel, not a valid
> common null. Current results and priorities: `PROJECT_MASTER_SUMMARY.md`
> sections 8-10 and 18. All other objectives, thresholds and timelines below are
> historical proposals unless explicitly identified as executed.

**Causal Grounding of Expression-Space Sparse Autoencoder Features in Hematopoietic Lineage Commitment: A Head-to-Head Test Against Foundation-Model Interpretability, with Clinical Translation**

> v6 originally superseded v5 by (1) reversing the proposed Phase 1 lineage direction
> under the SAE analysis, now known to be method-dependent; (2) proposing a controlled head-to-head
> between expression-space and foundation-model SAE features against the same CRISPRi
> null; (3) an explicit differentiation from Kendiukhov (2026) is added; (4) Phase 1
> intended inputs are corrected/upgraded artifacts (512-latent SAEs, log-normalized
> marker-aware preprocessing, verified erythroid marker symbols, a pre-registered
> control-referenced modularity metric). See PROJECT_AUDIT.md (section J) and
> STRATEGY_AND_POSITIONING.md.

---

## Review of Literature / Rationale

Establishing the latent structure of gene expression in hematopoietic lineage
commitment requires decomposing high-dimensional single-cell RNA sequencing
(scRNA-seq) data into discrete, interpretable regulatory programs. Hematopoiesis is
governed by a continuous, branching continuum of cell states (Orkin & Zon, 2008;
Velten et al., 2017), which traditional clustering resolves poorly because
overlapping gene modules exist in "superposition" (Schuster, 2025). Standard deep
learning pipelines provide necessary dimensionality reduction (Wolf et al., 2018;
Lopez et al., 2018; Eraslan et al., 2019; Ma & Xu, 2022) but frequently lack
"monosemanticity," the association of one feature with one biological concept
required for mechanistic insight (Bricken et al., 2023; Cunningham et al., 2023).
Sparse autoencoders (SAEs) address this by using L1-regularized sparse dictionary
learning to recover interpretable features from dense representations (Makhzani &
Frey, 2014; Rajamanoharan et al., 2024; Templeton et al., 2024).

In the corrected Gate-0 analysis, 512-latent SAEs were trained directly on log-normalized scRNA-seq
*expression* data — distinct from the recent wave of work that trains SAEs on the
internal activations of single-cell foundation models (Cui et al., 2024; Theodoris
et al., 2023; Kendiukhov, 2026a) — to isolate latent programs along the myeloid
commitment trajectory. The executed Gate 0 analyzed the Paul et al. (2015) mouse
atlas and Setty et al. (2019) human CD34+ marrow, 4,142 cells. An earlier exploratory
Phase 1 used a 25,000-cell CELLxGENE sample; it is not the frozen human arm.
Palantir/CellRank trajectory analyses belong to that earlier exploration rather than
establishing the Gate-0 cross-method result.
Crucially, the finding was established through an adversarial, pre-registered
framework that distinguishes real programs from measurement artifacts:
control-referenced feature strengths, count-fair concentration metrics, calibrated
permutation nulls, and a documented marker-curation audit. Under the frozen SAE
analysis, fewer erythroid than granulocyte programs are detected. Across five
decomposition families and two species the original erythroid-more-distributed
claim survives only in mouse NMF (1/10 method-by-species cells). The defensible
finding is that measured modularity depends on the decomposition. Neither direction
has been established as a metric-independent biological organization or causally
validated by the later perturbation tests.

Recovering biologically plausible features does not, however, guarantee causal
regulatory logic, and this is where the field's current frontier and this project's
central question lie. Systematic evaluation of single-cell foundation models
(scGPT: Cui et al., 2024; Geneformer: Theodoris et al., 2023) shows severe
mechanistic limits: Kendiukhov (2026a) trained SAEs on the residual-stream
activations of Geneformer and scGPT and found that only 6.2 percent (3 of 48) of
transcription-factor-associated *foundation-model* features respond specifically to
CRISPRi knockdown of their target, with organized biological knowledge but minimal
causal regulatory logic (Kendiukhov, 2026b; Kendiukhov, 2026c). This establishes a
concrete, published null — but it is a null about features living in a foundation
model's internal representation. It remains untested whether features learned
directly in *expression space*, annotated along a developmental trajectory, carry
more causal grounding than foundation-model features against the same genome-scale
perturbational data (Replogle et al., 2022). That head-to-head was proposed here
but remains unexecuted. The published percentage does not substitute for a matched
within-study comparison under the implemented metric.

The importance of closing it is amplified by two unanswered questions about where
these programs matter clinically and whether they are representation artifacts.
Single-cell atlases of Acute Myeloid Leukemia (AML) show that malignant hierarchies
profoundly dysregulate healthy hematopoiesis (van Galen et al., 2019), and tumors
reprogram bone-marrow progenitors through specific signaling axes to drive
immunosuppression and metastasis (LaMarche et al., 2024; Adrover et al., 2025). It
is unknown whether tumor-driven signals asymmetrically hijack the specific lineage
programs identified here (Hegde et al., 2025), or whether the observed modularity is
merely a representation-dependent artifact of expression-space processing that would
not survive re-derivation in a foundation-model embedding. Because the direction of
the asymmetry was itself corrected during Phase 1 (an earlier draft had it reversed,
until a marker-curation and metric audit overturned it), the project treats every
directional claim as a hypothesis to be re-tested, not assumed.

The original Phase 2 proposal had three pillars: causal grounding, clinical
translation and representation invariance. It proposed projecting expression SAEs
onto Replogle, comparing matched foundation-model SAEs, translating programs to
van Galen AML/Hegde tumor data and comparing modularity across embedding spaces.
Only the causal pillar was implemented, with a replacement metric and newly trained
dataset-specific TopK SAEs. The other pillars remain unexecuted. The proposal's
references to a "6.2 percent null" are superseded: that literature result used a
different panel and analysis. A future controlled comparison must be designed and
declared before its own execution, with appropriate uncertainty, multiplicity
controls and explicit negative/mixed-result reporting.

### Differentiation from prior SAE / single-cell interpretability work

The base technique (SAEs recover interpretable features from single-cell data) is
established (Schuster, 2025; Kendiukhov, 2026a; Cunningham et al., 2023). This
project does not claim it. Executed differentiation is direct expression-space
analysis plus an artifact-vs-signal framework that overturned its initial finding
and demonstrated structural method-dependence. A matched embedding head-to-head
and clinical/tumor translation were proposed differentiators, not completed
contributions. The reviewed abstracts do not establish that no prior work addresses
these questions. Kendiukhov (2026a) supplies literature context, not a directly
comparable statistical null; read its full methods before any numerical comparison.

---

## Research Question(s)

1. Implemented replacement for the original 24-TF question: do testable
   perturbations produce a held-out, correctly signed response in an SAE feature
   matched to their differential-expression signature, beyond declared controls?
   The implementation scores perturbations, not a fixed table of annotated feature
   indices. A matched embedding arm under the same metric remains a future question;
   the literature's 6.2 percent is not its statistical null.
2. Within matched cell-type bins (HSC, MPP, GMP, MEP, granulocytic, erythroid,
   monocyte), do SAE activations differ between AML and healthy cells with BH-FDR
   q < 0.05 in at least 10 features at |d| > 0.3?
3. Do feature activations stratify AML patients by NPM1, FLT3-ITD, TP53, RUNX1
   status at pseudobulk patient-level resolution, given typical sample sizes of 3 to
   6 patients per subtype?
4. In Hegde (2025) tumor-bearing bone marrow, does the erythroid unified module
   show coordinated suppression (its single dominant feature decreased) while
   granulocyte sub-module features show coordinated elevation (multiple distributed
   features increased), relative to naive controls?
5. Does an SAE trained on scGPT (and optionally Geneformer V2-316M) embeddings of
   Paul15 recover a modularity asymmetry in the same direction as the corrected
   expression-space SAE (erythroid unified, granulocyte distributed)?

---

## Objectives

- Revised implemented objective: test perturbation-specific responses with
  dataset-specific TopK SAEs and cross-fit matching/Mann-Whitney analysis. The
  Phase 1 projection, fixed 24-TF feature table and binomial-vs-6.2% design were
  replaced. The within-study foundation-model comparison was not completed.
- To apply Phase 1 SAEs to the van Galen (2019) AML atlas and quantify
  AML-versus-healthy differential activation within matched cell-type bins, with
  Benjamini-Hochberg FDR control across all (feature, cell-type) tests.
- To test patient-level stratification of AML cases by mutational subtype using
  pseudobulk SAE feature activation and mixed-effects modeling accounting for
  within-patient correlation.
- To test asymmetric reprogramming of erythroid versus granulocyte programs in Hegde
  (2025) tumor-bearing bone marrow versus naive controls, in the corrected
  direction (erythroid unified module suppressed, granulocyte distributed modules
  elevated).
- To test representation invariance of the modularity asymmetry by training parallel
  SAEs on scGPT (and optionally Geneformer V2-316M) embeddings of Paul15 and
  comparing recovered modularity direction and ratios.
- To produce a publication-ready manuscript integrating Phase 1 and Phase 2 with
  pre-specified statistical thresholds and explicit null-result reporting for any
  failed gate.

---

## Engineering Goals

N/A. This is a scientific investigation rather than an engineering project.

---

## Expected Outcomes

- Causal validation table covering 24 canonical hematopoietic TFs: per-TF effect
  size, permutation p-value, BH q-value, top-rank of the annotated feature, and
  pass/fail against pre-registered thresholds, reported in parallel for
  expression-space and foundation-model-embedding SAE features. Direct comparison
  plot of both validation rates against the Kendiukhov 6.2 percent benchmark.
- AML-versus-healthy feature activation matrix: features x cell-type bins tests per
  seed, with q-values and Cohen's d.
- Patient-level pseudobulk feature activation table with mixed-effects model output
  and at minimum one mutation-subtype distinction with power analysis.
- Hegde (2025) tumor-versus-naive differential activation per Phase 1 module
  assignment, with a chi-squared test of the corrected asymmetric prediction
  (erythroid unified module down, granulocyte distributed modules up).
- Cross-architecture comparison: CCA output, Hungarian-assignment feature-pair
  matching, and modularity-ratio bootstrap CIs for expression-space,
  scGPT-embedding-space, and (if compute permits) Geneformer-V2-embedding-space SAEs.

---

## Materials

### Equipment (make/model/quantity)
- Personal laptop/desktop, 1 unit: Windows 11 Home.
- Google Colab Pro+ (or access to a T4 / A100 / H200 GPU; equivalent compute
  acceptable), 51 GB CPU-memory tier.
- Google Drive 200 GB tier (or USB drive) for data storage.
- RunPod cloud GPU, NVIDIA A100 80 GB, on-demand: for full Replogle inference and
  Geneformer embedding-SAE training.

### Software (program and version)
- Python 3.11.
- PyTorch 2.4+ (validated with 2.12 in the current environment), CUDA 12.1+.
- scanpy 1.10+, anndata 0.10+, scvi-tools 1.1+, scikit-learn 1.4+.
- celltypist 1.6+ with Immune_All_Low.pkl reference model.
- statsmodels 0.14+ (MixedLM); pymer4 0.8+ as alternative.
- numpy 1.26+, pandas 2.1+, scipy 1.11+, matplotlib 3.8+, seaborn 0.13+.
- scGPT whole-human pretrained checkpoint (Cui et al., 2024): github.com/bowang-lab/scGPT;
  weights Zenodo DOI 10.5281/zenodo.10466117.
- Geneformer V2-316M pretrained checkpoint (Theodoris et al., 2023; Chen et al.,
  2024 quantized variant): huggingface.co/ctheodoris/Geneformer.
- Palantir 1.4+, CellRank 2.0+ (carried from Phase 1).
- captum 0.7+ for gradient-based attribution on embeddings.
- Git 2.40+, GitHub (private repository during development).
- GitHub electronic lab notebook + carbon-copy bound lab notebook (JSR-required).

### Datasets
- Replogle et al. (2022) K562 genome-wide Perturb-seq. figshare+ DOI
  10.25452/figshare.plus.20029387. Subset to 24 canonical hematopoietic TFs (GATA1,
  KLF1, TAL1, FLI1, CEBPA, CEBPE, GFI1, GFI1B, SPI1, RUNX1, LMO2, NFE2, MYB, BCL11A,
  IKZF1, MEF2C, HOXA9, MEIS1, ZBTB7A, IRF8, ETV6, EVI1, FOXO1, ID1) + non-targeting
  controls 1:3. Target 300,000-500,000 cells.
- van Galen et al. (2019) AML atlas. GEO GSE116256. 38,410 cells, 40 aspirates (16
  AML incl. timepoints; 5 healthy). Seq-Well. Processed Seurat object:
  petervangalen/reanalyze-aml2019.
- Hegde et al. (2025) tumor bone marrow. GEO GSE270148. Mouse scRNA-seq subset only
  (KP tumor-bearing + naive C57BL/6). Code: github.com/Merad-Lab/Hegde_Myelopoiesis_Epigenetics.
- scGPT whole-human pretrained model (12 layers, dim 512, ~33M cells; Cui et al., 2024).
- Geneformer V2-316M (18 layers, dim 1,152, ~95M cells; Theodoris et al., 2023;
  Chen et al., 2024). Stretch goal.
- **Executed Gate-0 checkpoints (not the causal pipeline's training inputs).**
  Mouse SAEs trained on marker-aware, log-normalized Paul15 (globins Hba-a2/Hbb-b1,
  Alas2, Ermap force-included; expression matrix MD5 recorded in
  preprocess_report.json), 512 latents (undercomplete relative to 2,012 genes,
  sparsity tuned to L0 in 20-50), 10 seeds. Human SAEs trained on Setty 2019 CD34+
  marrow, 4,142 cells x 2,032 genes, 512 latents and 10 seeds. The earlier
  128-latent, variance-only-HVG SAEs (MD5 60183a17983c8b977d036e0f3a58da61) are
  SUPERSEDED and retained only as a documented baseline.
- Legacy CZ CELLxGENE sample (CZI, 2023):
  chanzuckerberg.github.io/cellxgene-census/. Its Census release was not pinned;
  it was superseded for Gate 0 by the sourced Setty dataset. Do not invent a pin.
- Paul et al. (2015) mouse hematopoiesis atlas via scanpy.datasets.paul15(). 2,730
  cells, 19 myeloid progenitor clusters. NOTE: this source is RAW-scale (max ~168);
  it is log1p-normalized in preprocessing (correcting an earlier "pre-log" error).
- Optional Component 5: Park et al. (2020) human thymus atlas, ArrayExpress
  E-MTAB-8581 (~200,000 cells).
- HCOP ortholog mapping (Human-Mouse): genenames.org/tools/hcop/. Ensembl Compara
  orthology table as backup.
- Historical proposal: a frozen Phase 1 feature-to-TF annotation table. The
  implemented causal analysis instead matches features separately in each fitted
  dataset/seed; this document does not establish that the proposed table was frozen.

### Chemicals — N/A. Purely computational, public de-identified scRNA-seq.
### Kits — N/A. Purely computational, public de-identified scRNA-seq.
### Tissues/Cells — N/A. Public de-identified datasets only.
### Organisms — N/A. Public de-identified datasets only.

---

## Methodology

Five thematic components across a compressible ~12-week timeline. Each has
pre-specified pass/fail decision gates with FDR control and effect-size minima.
Random seed fixed at 42. Software versions pinned via Materials. Phase 1 corrected
deliverables are the inputs.

### Component 0: Phase 1 finalization (prerequisite; must complete before gates)
- Pre-register ONE modularity metric: the control-referenced decision
  (Progenitor-only baseline, abundance-matched permutation null, effect-size floor).
  SHA-freeze the decision spec before any Phase 2 re-run.
- Fix sparsity: L1 swept so mean L0 lands in 20-50; retrain >=10 seeds; report the
  L0/L1 sweep and dictionary-size sweep (256/512/1024/2048).
- Finish PCA and NMF baselines on identical inputs; report what the SAE adds.
- Marker sets: pre-registered, literature-cited, non-overlapping, with a
  leave-one-marker-out sensitivity analysis and a blinded GO/pathway cross-check.
- Human replication of the corrected asymmetry direction (mouse AND human).
- Freeze the corrected feature-to-TF annotation table (SHA-256).

### Component 1: van Galen 2019 AML application (weeks 1-3)
- Week 1 — acquisition/preprocessing: download GSE116256 (GEOparse or NCBI FTP);
  build unified AnnData with sample_id, donor_id, disease_status, timepoint. QC:
  gene_count 200-7,000, mito < 10 percent, UMI > 500 (document removals). Normalize
  per cell to 1e4, log1p. HVG via sc.pp.highly_variable_genes(flavor='seurat_v3',
  n_top_genes=2000) on raw counts. Integrate with scVI (batch=sample_id, n_latent=30,
  n_layers=2, dropout=0.1, max_epochs=400, patience=45). Annotate with celltypist
  (Immune_All_Low.pkl); validate by markers per bin (HSC: CD34/KIT/AVP; MPP:
  CD34/FLT3/SPINK2; GMP: ELANE/MPO/AZU1; MEP: GATA1/KLF1/GYPA; erythroid:
  HBB/HBA1/ALAS2; granulocytic: MPO/ELANE/S100A8/S100A9/CTSG; monocyte:
  CD14/LYZ/FCN1/VCAN; DC: FCER1A/CLEC4C). Restrict downstream analysis to the 7
  in-distribution bins matching Phase 1 human training (HSC, MPP, GMP, MEP,
  erythroid, granulocytic, monocyte); annotate but do not analyze OOD populations.
  Save AnnData + provenance JSON.
- Week 2 — projection/statistics: match Phase 1 SAE input gene set to van Galen;
  zero-impute missing genes only if missing fraction < 5 percent. Project features
  via the trained encoder per seed. QC gate: per-cell L0 in 20-50; reconstruction
  MSE within +50 percent of Phase 1 mean (debug before continuing if either fails).
  Per (feature, bin), Wilcoxon rank-sum AML vs healthy (128-or-more features x 7
  bins per seed). BH-FDR per seed. Cohen's d per (feature, bin, seed). Robustness:
  report only features significant (q < 0.05, |d| > 0.3) in >= 2 of 3 human seeds
  (or 3 of 5 mouse, cross-species).
- Week 3 — patient-level + mutations: aggregate per-patient per-bin per-feature mean
  activation. Fit statsmodels.MixedLM: feature_activation ~ disease_status +
  cell_type + (1 | patient_id); extract disease_status coefficient, Wald p, 95% CI.
  BH-FDR across features. Mutation stratification (NPM1, FLT3-ITD, IDH1, IDH2,
  CEBPA, TP53, RUNX1 from van Galen Supp. Table 1); Mann-Whitney U on pseudobulk for
  each pair with >= 3 patients/group. Power analysis per Cohen (1988). Heatmap
  (patients x features, annotated by status/mutation/timepoint).
- Decision Gate 1: PASS if >= 10 features q < 0.05, |d| > 0.3 in >= 2 seeds AND >= 1
  mutation distinction q < 0.10. PARTIAL 3-9 features or no mutation distinction
  (continue with documented limitation). FAIL < 3 features (stop; debug transfer/
  scaling).

### Component 2: implemented replacement and execution boundary

- The original 24-TF projection, SD-drop/top-rank table and 8/24 gate were replaced,
  not executed as this plan proposed. The historical text remains in Git history.
- Implemented code: `src/causal_pipeline.py` trains a dataset-specific TopK SAE;
  `src/causal_grounding.py` defines control-cell gene associations, cross-fit feature
  matching and a held-out Mann-Whitney response test with direction/floor, FDR,
  match-confidence and effect-size gates plus label-shuffle controls. Consult
  `config/causal_grounding_spec.json` (v7 plus erratum, SHA prefix `aae1a9d9`) for
  the complete rule and amendment chronology.
- As executed, SAE training used all cells, with perturbation labels unseen.
  Control-only training was intended but not implemented in the reported runs;
  the subsequently exposed mode fails the synthetic pipeline self-test (0/8).
- Reported data: Replogle K562 essential (162 tested TFs), RPE1, and later Norman
  CRISPRa. These are not a completed matched embedding comparison. Historical raw
  JSONs need recovery or regeneration before publication uses exact causal values.
- PCA/NMF runs transform the SAE input, then train the same class of SAE. They do
  not measure direct PCA/NMF causal-factor performance.
- No superiority claim follows from comparison with the literature's 6.2 percent.
  A future embedding arm needs matched cells, panel, training mode and settings,
  independent controls and a declared paired comparison with uncertainty.

### Component 3: Hegde 2025 tumor bone marrow projection (week 7)
- Download Hegde scRNA-seq (GSE270148); subset KP tumor-bearing + naive C57BL/6 BM
  myeloid progenitors and circulating monocytes; exclude lung TME. QC as Component 1.
  Both Paul15 and Hegde are mouse C57BL/6 (no cross-species mapping); verify MGI
  symbols. Integrate with Paul15 via scVI (batch=study_id). Annotate by scANVI label
  transfer to Paul15 clusters. Project the mouse SAE per seed. QC gate (L0 20-50; MSE
  within +50 percent). Per (bin, feature) Wilcoxon tumor vs naive; BH-FDR. Identify
  the erythroid unified module feature(s) and granulocyte distributed sub-module
  features from the annotation table.
- Asymmetry test (RQ4, corrected direction): chi-squared 2x2 on (erythroid unified
  feature suppressed at q < 0.05, unchanged) vs (granulocyte sub-module features
  elevated at q < 0.05, unchanged), Yates corrected. Concordance: Pearson of effect
  sizes within the granulocyte sub-module feature set (coordinated distributed
  elevation) and the erythroid module (coordinated suppression). Visualize per-feature
  tumor-vs-naive effect annotated by module.

### Component 4: Embedding-space SAE — controlled comparison (week 8)
- Framing: this component is a CONTROLLED COMPARISON, not a novel method (SAEs on
  scGPT/Geneformer embeddings are established; Kendiukhov, 2026a). Its outputs feed
  Component 2's head-to-head.
- scGPT (required): download weights (Zenodo). HCOP one-to-one ortholog map Paul15
  mouse -> human. Run scGPT inference; extract per-cell embeddings (n_cells, 512)
  from the best-silhouette layer on paul15_clusters. Validate embedding (Leiden
  res=1.0; report ARI vs paul15_clusters). Train an SAE on embeddings with matched
  proposed hyperparameters (a declared dictionary size, lambda tuned to L0 20-50, AdamW
  lr=5e-4, cosine annealing 5-epoch warmup, weight_decay=1e-5, grad-clip 1.0,
  batch 256, 300 epochs, dead-neuron resampling every 500 steps, decoder
  unit-normalized, >=5 seeds). Verify training (MSE, L0 20-50, zero dead features).
  Attribute features to top input genes via captum.IntegratedGradients; annotate
  against the same corrected canonical marker sets. Compute CCA on expression-space
  vs embedding-space activation matrices (first 5 canonical correlations). Hungarian
  matching on cosine similarity of annotation vectors. Modularity ratio per
  architecture with 1,000-iteration bootstrap 95% CIs. **Feed the annotated
  embedding-space features into Component 2 for the causal head-to-head.**
- Geneformer V2-316M (stretch, if RunPod A100 available): repeat on Geneformer
  embeddings (dim 1,152). Provides a direct methodological parallel to Kendiukhov
  (2026a), who used both models.
- Decision Gate 3: PASS if embedding-space modularity is in the same direction as
  the corrected expression-space finding (erythroid unified, granulocyte
  distributed) with same-sign bootstrap CIs. FAIL if different (publishable as a
  limitation: expression- and embedding-space SAEs capture different structure).

### Component 5 (optional): cross-system replication (week 9)
- Only if Components 1-4 are on schedule. Download Park et al. (2020) thymus atlas
  (E-MTAB-8581, ~200k cells). Identical preprocessing to Phase 1 (QC, log1p, scale
  max_value=10, HVG top 2,000 by variance). Train SAE (matched hyperparameters, 5
  seeds). Annotate against T-cell development markers (CD3D, CD3E, CD4, CD8A, CD8B,
  RAG1, RAG2, PTCRA, BCL11B, TCF7, GATA3, RORC, FOXP3). Test whether an analogous
  modularity asymmetry replicates; modularity ratio with bootstrap CIs.

### Buffer / post-experimental (weeks 9-12)
- Week 9: Component 5 if on track, else gap-fill.
- Week 10: write Methods (Phase 1 brief, Phase 2 full), IMRaD. Figures (matplotlib/
  seaborn): (1) corrected asymmetric modularity, (2) Replogle causal head-to-head
  vs 6.2 percent, (3) AML cross-disease + mutation stratification, (4) Hegde
  reprogramming, (5) expression-vs-embedding comparison.
- Week 11: assemble Intro/Methods/Results/Discussion/Limitations/Conclusion;
  mentor review; revise.
- Week 12: final formatting, supplementary tables, code archive, provenance.
- Decision Gate 4 (before writing): PASS if 3+/4 primary components met gates.
  PARTIAL 2/4 (compress scope). FAIL < 2 (reframe around the Phase 1 methodology +
  the causal head-to-head + partial Phase 2).

---

## Risk and Safety

Computational project on cloud infrastructure using publicly available,
de-identified data.

- Risk 1 — Extended computer use: digital eye strain, musculoskeletal discomfort.
  Precautions: 20-20-20 rule; ergonomic seating; breaks every 45-60 min; <= 8 h/day.
  Storage/disposal: N/A.
- Risk 2 — Data security: unauthorized access to work/credentials. Precautions: 2FA
  on Google/GitHub; private repo during development; no hardcoded credentials; tokens
  in Colab secrets; password manager. Storage: Google Drive with access logs +
  version history. Disposal: sensitive intermediates deleted post-submission;
  records retained 5+ years.
- Risk 3 — Reproducibility drift: Colab package updates break pipelines.
  Precautions: pin versions in requirements.txt; weekly MD5 verification of the
  matrix hash recorded in preprocess_report.json; reproducibility notebook in repo.
  Storage: requirements.txt in GitHub. Disposal: N/A.
- Risk 4 — Multiple-testing error: thousands of tests across components.
  Precautions: BH-FDR per component family; pre-specified thresholds before data
  inspection; cross-seed validation (2/3 human or 3/5 mouse). Storage: all p/q/effect
  matrices archived. Disposal: N/A.
- Risk 5 — Work overload / mental health: 12-week intensive schedule.
  Precautions: >= 1 day off/week; 7+ h sleep; weekly mentor check-ins; pause if
  warning signs appear. Storage/disposal: N/A.

---

## Data Analysis

Shared across components: BH-FDR control, Cohen's d effect sizes, mixed-effects
modeling for hierarchical data, cross-seed validation.

- Component 1 (van Galen AML): per-cell Wilcoxon within bins, BH-FDR across the test
  family; Cohen's d = (mean_AML - mean_healthy)/pooled_std; patient-level
  feature_activation ~ disease_status + cell_type + (1 | patient_id) via
  statsmodels.MixedLM; mutation stratification via Mann-Whitney U / Kruskal-Wallis
  on pseudobulk with Cohen (1988) power analysis; robustness >= 2/3 (or 3/5) seeds.
- Component 2 (implemented replacement): cross-fit matching and held-out
  Mann-Whitney AUC with the gates/nulls in the current causal spec. The 24x128
  annotated-feature family and binomial test against 6.2 percent were not the
  analysis that ran. See the corrected Component 2 section above.
- Component 3 (Hegde): per-bin Wilcoxon tumor vs naive, BH-FDR; corrected-direction
  asymmetry chi-squared (Yates); Pearson concordance within module sets; per-module
  effect with 95% bootstrap CIs.
- Component 4 (scGPT/Geneformer): CCA (first 5 canonical correlations); Hungarian
  assignment (scipy.optimize.linear_sum_assignment) on cosine similarity of
  annotation vectors; modularity ratio per architecture with 1,000-iteration
  bootstrap 95% CIs.
- Component 5 (optional): modularity ratio in thymus; one-sample t-test on the
  bootstrap distribution vs null ratio 1.0.

---

## Subject Specific Guidelines

- Human Participants Research: N/A. The project uses de-identified, publicly
  available scRNA-seq datasets; per ISEF this is not human-participants research.
- Vertebrate Animal Research: N/A. No animals used. Publicly available mouse data
  (Paul 2015; Hegde 2025) were generated by others under the Weizmann Institute and
  Icahn School of Medicine at Mount Sinai IACUC approvals.
- Potentially Hazardous Biological Agents: N/A. Fully computational; no biological
  agents.
- Hazardous Chemicals, Activities, and Devices: N/A. No hazardous chemicals,
  activities, or devices beyond standard computer/cloud use addressed under Risk and
  Safety.

---

## Bibliography

Adrover, J. M., Han, X., Sun, L., Fujii, T., Sivetz, N., Dassler-Plenker, J., Evans, C., Peters, J., He, X.-Y., Cannon, C. D., Ho, W. J., Raptis, G., Powers, R. S., & Egeblad, M. (2025). Neutrophils drive vascular occlusion, tumour necrosis and metastasis. Nature, 645, 484-495. https://doi.org/10.1038/s41586-025-09278-3

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: A practical and powerful approach to multiple testing. Journal of the Royal Statistical Society: Series B, 57(1), 289-300.

Bricken, T., Templeton, A., Batson, J., Chen, B., Jermyn, A., Conerly, T., Turner, N. L., Anil, C., Denison, C., Askell, A., Lasenby, R., Wu, Y., Kravec, S., Schiefer, N., Maxwell, T., Joseph, N., Tamkin, A., Nguyen, K., McLean, B., ... Olah, C. (2023). Towards monosemanticity: Decomposing language models with dictionary learning. Transformer Circuits Thread.

Bussmann, B., Leask, P., & Nanda, N. (2024). BatchTopK sparse autoencoders. arXiv preprint arXiv:2412.06410.

Bussmann, B., Nabeshima, N., Karvonen, A., & Nanda, N. (2025). Learning multi-level features with Matryoshka sparse autoencoders. arXiv preprint arXiv:2503.17547.

Chen, H., Venkatesh, M. S., Gomez Ortega, J., Mahesh, S. V., Nandi, T., Madduri, R., Pelka, K., & Theodoris, C. V. (2024). Quantized multi-task learning for context-specific representations of gene network dynamics. bioRxiv 2024.08.16.608180. https://doi.org/10.1101/2024.08.16.608180

Claye, C., Marschall, P., Ouerdane, W., Hudelot, C., & Duquesne, J. (2025). A framework to extract and interpret biological concepts from scRNAseq generative foundation models. ICML 2025 Generative AI and Biology (GenBio) Workshop.

Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). Lawrence Erlbaum Associates.

Cui, H., Wang, C., Maan, H., Pang, K., Luo, F., Duan, N., & Wang, B. (2024). scGPT: Toward building a foundation model for single-cell multi-omics using generative AI. Nature Methods, 21(8), 1470-1480. https://doi.org/10.1038/s41592-024-02201-0

Cunningham, H., Ewart, A., Riggs, L., Huben, R., & Sharkey, L. (2023). Sparse autoencoders find highly interpretable features in language models. arXiv preprint arXiv:2309.08600.

CZI Cell Science Program. (2023). CZ CELLxGENE Discover: A single-cell data platform for scalable exploration, analysis and modeling of aggregated data. bioRxiv. https://doi.org/10.1101/2023.10.30.563174

Eraslan, G., Avsec, Z., Gagneur, J., & Theis, F. J. (2019). Deep learning: new computational modelling techniques for genomics. Nature Reviews Genetics, 20(7), 389-403.

Haghverdi, L., Buttner, M., Wolf, F. A., Buettner, F., & Theis, F. J. (2016). Diffusion pseudotime robustly reconstructs lineage branching. Nature Methods, 13(10), 845-848.

Hegde, S., Giotti, B., Soong, B. Y., Halasz, L., Le Berichel, J., Schaefer, M. M., Kloeckner, B., Mattiuz, R., Park, M. D., Magen, A., Marks, A., Belabed, M., Hamon, P., Chin, T., Troncoso, L., Lee, J. J., Fan, K., Ahimovic, D., Bale, M. J., ... Merad, M. (2025). Myeloid progenitor dysregulation fuels immunosuppressive macrophages in tumours. Nature, 646(8087), 1214-1222. https://doi.org/10.1038/s41586-025-09493-y

Kendiukhov, I. (2026a). Sparse autoencoders reveal organized biological knowledge but minimal regulatory logic in single-cell foundation models: A comparative atlas of Geneformer and scGPT. arXiv preprint arXiv:2603.02952.

Kendiukhov, I. (2026b). Causal circuit tracing reveals distinct computational architectures in single-cell foundation models. arXiv preprint arXiv:2603.01752.

Kendiukhov, I. (2026c). Exhaustive circuit mapping of a single-cell foundation model reveals massive redundancy, heavy-tailed hub architecture, and layer-dependent differentiation control. arXiv preprint arXiv:2603.11940.

LaMarche, N. M., Hegde, S., Park, M. D., Maier, B. B., Troncoso, L., Le Berichel, J., Hamon, P., Belabed, M., Mattiuz, R., Hennequin, C., Chin, T., Reid, A. M., Reyes-Torres, I., Nemeth, E., Zhang, R., Olson, O. C., Doroshow, D. B., Rohs, N. C., Gomez, J. E., ... Merad, M. (2024). An IL-4 signalling axis in bone marrow drives pro-tumorigenic myelopoiesis. Nature, 625, 166-174. https://doi.org/10.1038/s41586-023-06797-9

Lange, M., Bergen, V., Klein, M., Setty, M., Reuter, B., Bakhti, M., Lickert, H., Ansari, M., Schniering, J., Schiller, H. B., Pe'er, D., & Theis, F. J. (2022). CellRank for directed single-cell fate mapping. Nature Methods, 19(2), 159-170. https://doi.org/10.1038/s41592-021-01346-6

Lopez, R., Regier, J., Cole, M. B., Jordan, M. I., & Yosef, N. (2018). Deep generative modeling for single-cell transcriptomics. Nature Methods, 15(12), 1053-1058.

Ma, A., & Xu, D. (2022). Deep learning for single-cell data analysis. Nature Reviews Molecular Cell Biology, 23(6), 405-426.

Makhzani, A., & Frey, B. J. (2014). k-sparse autoencoders. arXiv preprint arXiv:1312.5663.

Orkin, S. H., & Zon, L. I. (2008). Hematopoiesis: An evolving paradigm for stem cell biology. Cell, 132(4), 631-644. https://doi.org/10.1016/j.cell.2008.01.025

Park, J. E., Botting, R. A., Dominguez Conde, C., Popescu, D. M., Lavaert, M., Kunz, D. J., Goh, I., Stephenson, E., Ragazzini, R., Tuck, E., Wilbrey-Clark, A., Roberts, K., Kedlian, V. R., Ferdinand, J. R., He, X., Webb, S., Maunder, D., Vandamme, N., Mahbubani, K. T., ... Teichmann, S. A. (2020). A cell atlas of human thymic development defines T cell repertoire formation. Science, 367(6480), eaay3224. https://doi.org/10.1126/science.aay3224

Paul, F., Arkin, Y., Giladi, A., Jaitin, D. A., Kenigsberg, E., Keren-Shaul, H., Winter, D., Lara-Astiaso, D., Gury, M., Weiner, A., David, E., Cohen, N., Lauridsen, F. K. B., Haas, S., Schlitzer, A., Mildner, A., Ginhoux, F., Jung, S., Trumpp, A., Porse, B. T., Tanay, A., & Amit, I. (2015). Transcriptional heterogeneity and lineage commitment in myeloid progenitors. Cell, 163(7), 1663-1677. https://doi.org/10.1016/j.cell.2015.11.013

Rajamanoharan, S., Lieberum, T., Sonnerat, N., Conmy, A., Varma, V., Kramar, J., & Nanda, N. (2024). Jumping ahead: Improving reconstruction fidelity with JumpReLU sparse autoencoders. arXiv preprint arXiv:2407.14435.

Replogle, J. M., Saunders, R. A., Pogson, A. N., Hussmann, J. A., Lenail, A., Guna, A., Mascibroda, L., Wagner, E. J., Adelman, K., Lithwick-Yanai, G., Iremadze, N., Oberstrass, F., Lipson, D., Bonnar, J. L., Jost, M., Norman, T. M., & Weissman, J. S. (2022). Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq. Cell, 185(14), 2559-2575.e28. https://doi.org/10.1016/j.cell.2022.05.013

Schuster, V. (2025). Can sparse autoencoders make sense of gene expression latent variable models? arXiv preprint arXiv:2410.11468v3.

Setty, M., Kiseliovas, V., Levine, J., Gayoso, A., Mazutis, L., & Pe'er, D. (2019). Characterization of cell fate probabilities in single-cell data with Palantir. Nature Biotechnology, 37(4), 451-460. https://doi.org/10.1038/s41587-019-0068-4

Templeton, A., Conerly, T., Marcus, J., Lindsey, J., Bricken, T., Chen, B., Pearce, A., Citro, C., Ameisen, E., Jones, A., Cunningham, H., Turner, N. L., McDougall, C., MacDiarmid, M., Freeman, C. D., Sumers, T. R., Rees, E., Batson, J., Jermyn, A., ... Olah, C. (2024). Scaling monosemanticity: Extracting interpretable features from Claude 3 Sonnet. Transformer Circuits Thread.

Theodoris, C. V., Xiao, L., Chopra, A., Chaffin, M. D., Al Sayed, Z. R., Hill, M. C., Mantineo, H., Brydon, E. M., Zeng, Z., Liu, X. S., & Ellinor, P. T. (2023). Transfer learning enables predictions in network biology. Nature, 618(7965), 616-624. https://doi.org/10.1038/s41586-023-06139-9

van Galen, P., Hovestadt, V., Wadsworth II, M. H., Hughes, T. K., Griffin, G. K., Battaglia, S., Verga, J. A., Stephansky, J., Pastika, T. J., Lombardi Story, J., Pinkus, G. S., Pozdnyakova, O., Galinsky, I., Stone, R. M., Graubert, T. A., Shalek, A. K., Aster, J. C., Lane, A. A., & Bernstein, B. E. (2019). Single-cell RNA-seq reveals AML hierarchies relevant to disease progression and immunity. Cell, 176(6), 1265-1281.e24. https://doi.org/10.1016/j.cell.2019.01.031

Velten, L., Haas, S. F., Raffel, S., Blaszkiewicz, S., Islam, S., Hennig, B. P., Hirche, C., Lutz, C., Buss, E. C., Nowak, D., Boch, T., Hofmann, W.-K., Ho, A. D., Huber, W., Trumpp, A., Essers, M. A. G., & Steinmetz, L. M. (2017). Human haematopoietic stem cell lineage commitment is a continuous process. Nature Cell Biology, 19(4), 271-281.

Weiler, P., Lange, M., Klein, M., Pe'er, D., & Theis, F. (2024). CellRank 2: Unified fate mapping in multiview single-cell data. Nature Methods. https://doi.org/10.1038/s41592-024-02303-9

Wolf, F. A., Angerer, P., & Theis, F. J. (2018). SCANPY: Large-scale single-cell gene expression data analysis. Genome Biology, 19(1), 15. https://doi.org/10.1186/s13059-017-1382-0

---

## Addendum

Repository-status correction dated 2026-10-18: the notices and implemented-method
replacement above distinguish historical proposals from executed work. They do not
retroactively pre-register results or assert the status of an official ISEF filing.
The original version of this plan is preserved in Git history.
