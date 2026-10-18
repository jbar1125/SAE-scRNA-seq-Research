# DIFFERENTIATION: this project and the SAE / single-cell field

**CURRENT STATUS, 2026-10-18.** The July claims of a solved essential-gene confound,
a completed expression-vs-embedding comparison and a validated GATA1 anchor are
withdrawn. The supported contribution is an expression-space interpretability audit.
See [PROJECT_MASTER_SUMMARY.md](PROJECT_MASTER_SUMMARY.md), sections 8-10 and 16.

## Competitor landscape recorded on 2026-07-18

These descriptions were checked against search abstracts, not full texts; direct
arXiv full-text retrieval returned 403 in that session. Re-read the papers and update
the search before submission. This table cannot establish field-wide priority.

| paper | arXiv | reported scope | representation |
|---|---|---|---|
| Kendiukhov SAE atlas of Geneformer/scGPT | 2603.02952 | TopK SAEs; reported 3/48 TFs (6.2%, 10.4% multi-tissue) with regulatory-target-specific CRISPRi response | foundation-model residual streams |
| Kendiukhov Causal Circuit Tracing | 2603.01752 | Feature ablation; reported 56.4% directional CRISPRi accuracy | foundation-model embeddings |
| Kendiukhov Exhaustive Circuit Mapping | 2603.11940 | Redundancy and hub structure of feature circuits | foundation-model embeddings |
| Hematopoietic manifold in scGPT | 2603.10261 | Extraction of a compact hematopoietic manifold | scGPT internals |
| CytoSAE | 2507.12464 | Blood-cell microscopy image features | images |

The Kendiukhov cluster tests a different input representation from this project's
direct expression SAEs. CytoSAE uses a different modality. Other relevant literature
includes arXiv 2410.11468, arXiv 2510.25807 and bioRxiv 2025.10.22.681631; their
methods must be reviewed before saying expression-space SAEs or the audit are unique.

## Differentiation and its evidence boundary

| axis | This project as executed | Boundary |
|---|---|---|
| SAE input | Gene expression; PCA/NMF transformations also fed to SAEs | scGPT/Geneformer arm never run |
| Structural validation | Frozen Gate 0, two species, five method families, marker and threshold checks | The reversed erythroid/granulocyte direction is method-dependent |
| Perturbation validation | Reported Replogle K562/RPE1, Norman rate and regulon tests | Single-perturbation lineage-regulatory claim not established; historical raw JSONs missing |
| Essential-gene analysis | Wrong-background enrichment detected and corrected | No supported "~10x enrichment"; the column-specificity fix failed |
| Training design | Prior causal SAEs trained on all cells without perturbation labels | Control-only training was intended but not executed; its synthetic self-test fails |
| Hematopoietic program | Stable erythroid expression program observed in K562 | Program-specific null test still needed; GATA1 absent in Norman |
| Causal baselines | Expression -> SAE, PCA -> SAE and NMF -> SAE ladder | Does not compare SAE activations against direct PCA/NMF factors |
| Expression vs embedding | Proposed matched experiment | No within-study superiority claim; 6.2% is literature context, not a shared null |

The earlier enrichment compared 8/26 seed-stable hits with all 1,839 listed TFs,
although only 162 were tested. The reported v4 split is 6/51 core-essential-tagged
versus 14/111 other tested TFs; the second group is not a validated nonessential
reference panel. Column specificity retained housekeeping hits and reduced GATA1's
support to 2/5 seeds, so it was disabled by default. The M1 effect-size filter
retained the same 22/162 seed-0 hits; it does not support an effect-size-confound pitch.

## The contribution that can currently be defended

1. **A bounded structural result:** the chosen decomposition changes measured
   lineage modularity on identical matrices, even when seeds and analysis thresholds
   are checked. This is an empirical interpretability limitation.
2. **An expression-space perturbation audit:** reported tests expose the gap between
   annotated features and the particular regulatory claim tested. This complements
   embedding-space work but is not an identical reproduction of its metric.
3. **A documented diagnostic record:** wrong-background enrichment, failed column
   specificity, gate saturation, hub-target selection, TopK competition and the
   training-mode mismatch each have explicit checks. Do not describe all candidate
   mechanisms as experimentally established causes.

A general benchmark, causal program atlas, clinical translation and matched
foundation-model comparison are future deliverables, not completed differentiators.

## A usable one-paragraph answer

"Kendiukhov studies sparse autoencoders on the internal representations of single-cell
foundation models. I study expression-space features and ask how their biological
interpretation changes under adversarial structural and perturbation tests. The
frozen structural analysis shows that inferred lineage modularity depends on the
decomposition. The later perturbation analyses did not establish the intended
lineage-regulatory claim and exposed several metric and provenance failures. My
contribution is that audit and its diagnostics. I have not yet run a matched
expression-vs-embedding comparison or shown that expression SAEs outperform it."

## Before submission

- Read the competitor full texts and compare the exact TF panels, training cells,
  feature definitions and nulls. Do not compare percentages as though the metrics match.
- Recover or regenerate the causal raw-output archive and archive the complete
  program-modulation test. Preserve failed and mixed results.
- Distinguish K562 leukemia biology from primary hematopoiesis and document each
  dataset's perturbation coverage. GATA1 in Replogle does not imply GATA1 in Norman.
- Refresh citations and novelty claims without promising publication, priority or
  competition outcomes.
