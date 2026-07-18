# DIFFERENTIATION — this project vs the 2025-2026 SAE / single-cell field

The explicit differentiation the strategy doc marks CRITICAL (V-A) and Gate 3 requires.
Built 2026-07-18 from VERIFIED sources (arXiv abstracts confirmed via search; see the
verification log at the bottom). A judge who knows this field WILL ask "how is this not
Kendiukhov's atlas?" -- this is the rehearsed, sourced answer.

## The verified competitor landscape

There is now an active cluster of SAE + single-cell-foundation-model papers, most from one
prolific group (Ihor Kendiukhov, Tubingen / Biodyn-AI):

| paper | arXiv | what it does | space |
|-------|-------|--------------|-------|
| SAE atlas of Geneformer/scGPT (the 6.2% source) | 2603.02952 | TopK SAEs on Geneformer V2-316M + scGPT residual streams (82,525 / 24,527 features); 29-59% annotate to GO/KEGG/TRRUST; **only 3/48 TFs (6.2%, 10.4% multi-tissue) show regulatory-target-specific responses vs Replogle CRISPRi** | embedding |
| Causal Circuit Tracing | 2603.01752 | Ablates SAE features in the same models, traces downstream; **gene-level CRISPRi validation 56.4% directional accuracy -> "co-expression, not causal encoding"** | embedding |
| Exhaustive Circuit Mapping | 2603.11940 | Redundancy / hub architecture of the same SAE feature circuits | embedding |
| Hematopoietic manifold in scGPT | 2603.10261 | Extracts a compact ~8-10D hematopoietic manifold from scGPT internals; a performant algorithm from FM internals | embedding |
| CytoSAE | 2507.12464 | SAE on 40k+ blood-cell **microscopy IMAGES** (MICCAI 2025 cytology) | image (not transcriptomics) |

Two would-be collisions, resolved:
- **CytoSAE is NOT a collision** -- it is computer vision on blood-cell images, a different
  modality entirely. Hematology + SAE, but nothing to do with expression-space grounding.
- **The scGPT hematopoietic-manifold paper is ADJACENT, not a collision** -- it extracts
  structure FROM scGPT's internals (embedding space); it does not train SAEs on raw
  expression or run an expression-vs-embedding causal head-to-head.

**The load-bearing fact: every one of these works in foundation-model EMBEDDING space.
None trains SAEs directly on gene EXPRESSION. None runs the expression-vs-embedding causal
head-to-head. That question is open, and their own conclusion ("embeddings encode
co-expression, not causal logic") is exactly what motivates it.**

## The differentiation table

| axis | This project | Kendiukhov atlas (2603.02952) | Kendiukhov circuit tracing (2603.01752) | scGPT manifold (2603.10261) |
|------|-------------|------------------------------|----------------------------------------|----------------------------|
| SAE input space | **gene EXPRESSION** (+ embedding as a within-study control arm) | Geneformer/scGPT residual stream | Geneformer/scGPT residual stream | scGPT internals |
| Core question | does causal logic live in expression that embeddings discard? | do FM features encode causal logic? (no) | same, via ablation circuits | can we extract an algorithm from scGPT? |
| Causal test | cross-fit + Mann-Whitney + **column-specificity** grounding on Replogle | target-specific response, 48 TFs | feature-ablation directional accuracy | none (structure extraction) |
| Essential-gene confound control | **yes -- quantified (~10x) and gated (column specificity, CEGv2 split)** | not reported | not reported | n/a |
| Head-to-head expression vs embedding | **yes (the deliverable)** | no (embedding only) | no (embedding only) | no |
| Hematopoiesis lineage biology | **yes (GATA1-anchored, K562 myeloid + Gate-0 mouse/human)** | pan-tissue atlas | pan-tissue | hematopoiesis (manifold) |
| Artifact-vs-signal rigor framework | **yes (Gate 0: pre-registration, method-dependence, calibrated nulls)** | no | no | partial (external panels) |

## Five things that are genuinely ours (rank-ordered)

1. **Expression-space, not embedding-space.** The entire competitor cluster is embedding-
   space. We ask the inverted, unasked question their results set up.
2. **The essential-gene confound + column-specificity metric.** We found (quantified ~10x)
   and fixed a confound -- essential-knockdown global collapse masquerading as grounding --
   that the embedding papers do not appear to control for. A confound-aware causal metric
   is a methodological contribution in its own right.
3. **A within-study head-to-head** (same cells, same SAE, same metric, matched grid) rather
   than comparing across papers -- internally valid regardless of reproducing their 6.2%.
4. **The artifact-vs-signal framework from Gate 0** (pre-registration, cross-species/method
   modularity method-dependence, calibrated nulls) as the credibility spine.
5. **Hematopoietic lineage biology with a real positive control (GATA1).**

## The rehearsed one-paragraph answer

"Kendiukhov's atlas and circuit-tracing papers train sparse autoencoders on the internal
activations of single-cell foundation models -- Geneformer and scGPT -- and find that those
features encode organized biological knowledge but almost no causal regulatory logic: only
about six percent of transcription factors show target-specific responses to CRISPRi. My
project asks the question their result makes obvious but nobody has run: is that a limit of
sparse autoencoders, or of the foundation-model embedding? I train the same kind of sparse
autoencoder directly on gene expression and run the identical causal test, head-to-head
against the embedding version on the same cells. Along the way I found and fixed a confound
their metric doesn't control for -- knocking down an essential gene collapses transcription
globally and fakes a grounding signal -- so my metric adds a column-specificity test that
demands a perturbation suppress its own program more than others do. The contribution is not
'SAEs on single cells' -- that's published -- it's locating WHERE causal regulatory logic is
legible, with a confound-aware metric and a hematopoietic positive control."

## Honest risks

- **The competitor is prolific and fast** (multiple papers in early 2026, public code at
  Biodyn-AI). They could run the expression-space comparison themselves. Mitigation: move
  fast; lean on the differentiators they are unlikely to replicate quickly -- the
  confound-aware metric, the hematopoiesis biology, the Gate-0 framework, and the
  representation-ladder framing (COMPONENT2_NEXT_DIRECTIONS C1).
- **K562 is a leukemia line**, not primary hematopoiesis (stated everywhere; RPE1 + a
  primary set are on the roadmap).
- **Metric comparability.** Our grounding metric differs from theirs by design (more
  conservative). We must run Arm A ourselves with OUR metric so the head-to-head is
  apples-to-apples, and NOT claim to beat their 6.2% until we do.

## Verification log (2026-07-18)

- 2603.02952 (Kendiukhov atlas): CONFIRMED real; 6.2% (3/48 TFs), 10.4% multi-tissue;
  TopK SAEs on Geneformer/scGPT residual streams. Matches the project's cited null.
- 2603.01752 (Causal Circuit Tracing, Kendiukhov): CONFIRMED; feature-ablation CRISPRi
  validation 56.4% directional -> "co-expression not causal." A second embedding-space
  causal paper; still not expression-space.
- 2603.11940 (Exhaustive Circuit Mapping): CONFIRMED same space/group.
- 2507.12464 (CytoSAE): CONFIRMED = microscopy images, MICCAI 2025. NOT a collision.
- 2603.10261 (scGPT hematopoietic manifold): CONFIRMED = extraction from scGPT internals.
  Adjacent, not a collision.
- arXiv full text was 403 via direct fetch here; claims verified via search abstracts.
  Re-read the PDFs directly before any submission (cross-check the exact 48-TF panel and
  their grounding definition against ours).
