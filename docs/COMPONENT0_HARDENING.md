# COMPONENT 0 HARDENING — stress-testing the frozen Gate-0 verdict

Follow-up to the frozen Gate 0 (`COMPONENT0_RESULTS.md`). Gate 0 concluded: the
original "asymmetric modularity" claim is NOT SUPPORTED in mouse and human; the data
show granulocyte MORE distributed than erythroid; and the modularity verdict is
method-dependent. This document collects the robustness battery that stress-tests
that conclusion. Every unit was executed in-container on CPU; artifacts are the small
JSON summaries under `data_g0*/robustness/` and `data_g0*/trrust/`.

**Bottom line: the verdict survives every stress test.** No seed count, sparsity
level, HVG count, or decision-threshold choice resurrects the original claim, and the
granulocyte-more-distributed direction holds throughout.

## Battery

| Unit | What it varies | Result |
|------|----------------|--------|
| **R1** decision-parameter sensitivity | effect-size floor x significance threshold (12 configs) | **0/12 support the original claim; granulocyte >= erythroid at every point, BOTH species.** Not a knob-tuning artifact. |
| **R2** L0-band sensitivity | L1 -> mean L0 across (and just past) the 20-50 band | NOT SUPPORTED at L0 27.7 / 42.2 / 52.0; granulocyte >= erythroid throughout. Not specific to L0=37. |
| **R3** seed stability | 10 seeds -> 20 seeds (via `--seed-start`) | gran_n_real = 3 in **20/20** seeds; granulocyte > erythroid in 20/20; 0/20 support the claim. n=10 was sufficient. |
| **R4** HVG-count sensitivity | 1000 / 2000 / 3000 highly-variable genes | NOT SUPPORTED at all three; granulocyte 3 / erythroid 1-2. (Caveat: at fixed l1=0.8, L0 leaves the band at 1000 and 3000 genes since L0 depends on gene count; verdict stable regardless.) |
| **TRRUST** ground-truth cross-reference | curated TF->target regulon (TRRUST v2) | Human: erythroid TF targets enriched in the erythroid-TF feature in **100% of seeds** (median p=6.4e-5); granulocyte 50%. SAE features recover real regulatory structure, and the unified erythroid program captures its regulon more tightly than the distributed granulocyte one. Mouse cross-ref underpowered (human DB, small panel) and not used as evidence. |

## Interpretation

- The core Gate-0 negative is not fragile: it is invariant to the analysis choices a
  skeptic would probe (thresholds, sparsity, seeds, gene count).
- TRRUST adds an independent, orthogonal line of support: the SAE is not hallucinating
  gene sets; the feature that captures a lineage's TFs also captures their known
  targets, and it does so more consistently for the unified (erythroid) program than
  the distributed (granulocyte) one, which is exactly what the modularity result
  predicts.
- The method-dependence caveat from Gate 0 still stands (PCA/NMF/SAE disagree); the
  hardening does not remove it, and the causal test (Component 2) remains the arbiter.

## Deferred / not run here (documented, not hidden)

- **SCENIC 5th baseline**: cisTarget motif databases are proxy-blocked; full SCENIC
  is infeasible in-container. TRRUST covers the curated-regulon validation. Runner in
  `LAB_NOTEBOOK.md` / to run on an unblocked box.
- **Second human dataset** (non-CD34 marrow): would test whether the human globin
  non-significance is a CD34-selection artifact. Needs reachable data + CPU; deferred.
- **Component 2** (causal CRISPRi head-to-head vs the 6.2% null): GPU-blocked; fully
  specified in `docs/COMPONENT2_PLAN.md`.

## Reproduce

```bash
python3 src/robustness_sweeps.py --data-dir ./data_g0        --out ./data_g0/robustness        # R1 (mouse)
V0B_MARKERS_JSON=config/human_markers.json \
  python3 src/robustness_sweeps.py --data-dir ./data_g0_human --out ./data_g0_human/robustness  # R1 (human)
# R2/R3/R4 are retraining sweeps (see LAB_NOTEBOOK.md for the exact drivers)
V0B_MARKERS_JSON=config/human_markers.json \
  python3 src/trrust_crossref.py --data-dir ./data_g0_human --species human --out ./data_g0_human/trrust
```
