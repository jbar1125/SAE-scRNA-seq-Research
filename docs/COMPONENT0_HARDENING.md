# COMPONENT 0 HARDENING — stress-testing the frozen Gate-0 verdict

**CURRENT STATUS, 2026-10-18.** The structural checks below concern Gate 0 and its
SAE analysis, not a causal validation of lineage organization. The executed human
arm is Setty 2019 CD34+ marrow (4,142 cells). Replogle/Norman causal work ran later;
it did not establish the intended lineage-regulatory claim, and the embedding arm
remains unrun. See `PROJECT_MASTER_SUMMARY.md` sections 5 and 8-10.

Follow-up to the frozen Gate 0 (`COMPONENT0_RESULTS.md`). Gate 0 concluded: the
original "asymmetric modularity" claim is NOT SUPPORTED in mouse and human; the data
show granulocyte MORE distributed than erythroid under the SAE analysis; and the modularity verdict is
method-dependent. This document collects the robustness battery that stress-tests
that conclusion. Every unit was executed in-container on CPU; artifacts are the small
JSON summaries under `data_g0*/robustness/` and `data_g0*/trrust/`.

**The SAE rejection survives the tested robustness settings.** No tested seed count,
sparsity level, HVG count or decision threshold resurrects the original claim;
granulocyte counts are at least erythroid counts, with a tie in some settings.
This does not make the reversed direction invariant across decomposition methods.

## Figures

Regenerate with `python3 src/make_figures.py` (reads only the committed JSON artifacts).

![Method-dependence across 5 families x 2 species](figures/fig_method_dependence.png)

*Fig 1. Median real-program counts (erythroid/granulocyte) per decomposition. The
original claim needs erythroid > granulocyte (a bar pointing right, warm); it appears
in exactly one of ten method x species cells (mouse NMF).*

![Mouse verdict stable across every stress test](figures/fig_robustness.png)

*Fig 2. Granulocyte minus erythroid real-program count across every robustness
condition; nonnegative values mean these settings do not support the original claim.*

## Battery

| Unit | What it varies | Result |
|------|----------------|--------|
| **R1** decision-parameter sensitivity | effect-size floor x significance threshold (12 configs) | **0/12 support the original claim; granulocyte >= erythroid at every point, BOTH species.** Not a knob-tuning artifact. |
| **R2** L0-band sensitivity | L1 -> mean L0 across (and just past) the 20-50 band | BOTH species. Mouse: NOT SUPPORTED at L0 27.7 / 42.2 / 52.0 (gran >= ery throughout). Human: NOT SUPPORTED at L0 21.4 / 32.5 / 41.5 (gran 2 / ery 1 throughout). Not specific to the frozen L0. |
| **R3** seed stability | 10 seeds -> 20 seeds (via `--seed-start`) | BOTH species. Mouse: gran_n_real = 3 in **20/20** seeds, 0/20 support the claim. Human: gran 2 / ery 1 stable over 20 seeds, 0/20 support. n=10 was sufficient in both. |
| **R4** HVG-count sensitivity | 1000 / 2000 / 3000 highly-variable genes | BOTH species NOT SUPPORTED at all three counts; granulocyte >= erythroid throughout. (Caveat: at fixed l1=0.8, L0 leaves the band at the gene-count extremes since L0 depends on gene count; at human 1000-genes both lineages collapse to 1 (a tie, still not the claim). Verdict/direction stable regardless.) |
| **TRRUST** curated-reference cross-reference | curated TF->target regulon (TRRUST v2) | Human: erythroid TF targets enriched in the erythroid-TF feature in **100% of seeds** (median p=6.4e-5); granulocyte 50%. This is lineage-level annotation agreement, not perturbation validation or proof of causal regulation. Mouse cross-ref underpowered (human DB, small panel) and not used as evidence. |
| **GRN** 4th decomposition family | co-expression regulons (|Pearson|; NOT motif-pruned SCENIC) run through the same v3.1 modularity | Adds a 4th answer: mouse ery=2/gran=2 (tie); human ery=1/gran=2. Still NOT SUPPORTED. |
| **ICA** 5th decomposition family | FastICA independent components (rank 50, converged 3/3; rank 512 is non-convergent) through the same v3.1 modularity | mouse ery=1/gran=3, human ery=2/gran=3. NOT SUPPORTED; agrees with the SAE granulocyte direction in both species. |

### Method-dependence, across 5 decomposition families (both species)

ery_n_real / gran_n_real per method; the original claim needs erythroid > granulocyte.

| method | mouse (ery/gran) | human (ery/gran) | supports original claim? |
|--------|:---------------:|:----------------:|--------------------------|
| PCA (linear) | 1 / 1 | 1 / 3 | no (mouse resolves nothing; human granulocyte) |
| ICA (independent, rank 50*) | 1 / 3 | 2 / 3 | no (granulocyte, both) |
| NMF (parts-based) | 2.5 / 2 | 0 / 0 | mouse **yes** (the lone outlier); human nothing |
| SAE (sparse) | 1.5 / 3 | 1 / 2 | no (granulocyte, both) |
| GRN (co-expression) | 2 / 2 | 1 / 2 | no (mouse tie; human granulocyte) |

*ICA is run at rank 50 (its convergent regime): FastICA at rank 512 does not converge
(0-1/3 seeds), and non-converged components are unreliable; at rank 50 it converges
3/3 in both species. The rank differs from the other families (512), stated openly.
The non-converged rank-512 run had falsely flipped human ICA to "supported" - a
convergence artifact, corrected here (see LAB_NOTEBOOK.md).

Five decomposition families give divergent modularity answers on the identical matrix.
The original claim (erythroid more distributed) survives in exactly **1 of 10
method x species combinations** (mouse NMF); every
other combination rejects it, and the SAE's granulocyte-more-distributed direction is
echoed by ICA (both species) and human PCA/GRN. This is the artifact-vs-signal thesis,
now demonstrated across five method families and two species: "gene-program modularity"
is co-determined by the chosen decomposition, so no single-method modularity claim is
biological ground truth.

## Interpretation

- The original-claim rejection is stable over the tested SAE thresholds, sparsity,
  seeds and gene counts. The cross-method table limits that statement's scope.
- TRRUST supplies lineage-level annotation agreement under a curated reference.
  It neither tests perturbation response nor establishes decomposition-independent
  biological modules. The later single-TF regulon test is a different experiment.
- The method-dependence caveat from Gate 0 still stands (PCA/NMF/SAE disagree); the
  hardening does not remove it. Later Component 2 tests did not causally validate
  the structural reversal.

## Deferred / not run here (documented, not hidden)

- **Motif-pruned SCENIC baseline**: cisTarget motif databases were proxy-blocked; full SCENIC
  is infeasible in-container. TRRUST covers the curated-regulon validation. Runner in
  `LAB_NOTEBOOK.md` / to run on an unblocked box.
- **Second human dataset** (non-CD34 marrow): would test whether the human globin
  non-significance is a CD34-selection artifact. Needs reachable data + CPU; deferred.
- **Component 2 update:** expression-space Replogle and Norman analyses ran later;
  the matched embedding comparison remains unrun. See `COMPONENT2_RESULTS.md` and
  `PROJECT_MASTER_SUMMARY.md`; `COMPONENT2_PLAN.md` preserves a superseded design.

## Reproduce

```bash
python3 src/robustness_sweeps.py --data-dir ./data_g0        --out ./data_g0/robustness        # R1 (mouse)
V0B_MARKERS_JSON=config/human_markers.json \
  python3 src/robustness_sweeps.py --data-dir ./data_g0_human --out ./data_g0_human/robustness  # R1 (human)
# R2/R3/R4 are retraining sweeps (see LAB_NOTEBOOK.md for the exact drivers)
V0B_MARKERS_JSON=config/human_markers.json \
  python3 src/trrust_crossref.py --data-dir ./data_g0_human --species human --out ./data_g0_human/trrust
```
