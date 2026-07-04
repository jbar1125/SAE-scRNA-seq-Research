# VERSION HISTORY

The version lineage of the analysis, so no prior version is "lost" and anyone can
see why each was replaced. The authoritative byte-level record is `git log`; this is
the human-readable map. Nothing here is deleted from the repo: every version is
recoverable from git, and the metric evolution is documented rather than overwritten.

## Metric (the decision that says SUPPORTED / NOT SUPPORTED)

| Version | File | Method | Why replaced |
|---------|------|--------|--------------|
| v1 (Drive) | `v0b_module_definitions.ipynb` (NOT in repo; Drive only) | cross-seed index matching, `H_ery/H_gran` ratio, markers absent from HVG | Multiple bugs; is the "notebook trap". Superseded, never run for a real verdict. |
| v2 | `src/v0b_module_definitions.py` | hypergeometric + BH-FDR + per-seed winner-take-all (gated) | The q<0.05 gate is structurally blind to distributed signal; its NOT-SUPPORTED is an artifact, not a clean negative. |
| v3 | `src/v0b_v3_loading.py` | continuous decoder-loading, no gate | Removes the blindness, but has a count/abundance bias and no baseline; its SUPPORTED verdict is an artifact. |
| **v3.1 (current, pre-registered)** | `src/v0b_v3_1_decision.py` | Progenitor-referenced + abundance-matched permutation null + BH-FDR + effect-size floor + 60% rule | THE metric. Frozen spec SHA `fc342829` (`config/preregistration_spec.*`). |

v2 and v3 are retained in `src/` on purpose: v3.1 imports v3's loading functions,
and both are part of the documented "the metric is the experiment" story.

## Data pipeline

| Version | Change |
|---------|--------|
| Original | generic variance HVG; wrong globin symbols (`Hba-a1`/`Hbb-bs`); no `Alas2`. Produced the false "erythroid undetected" result. |
| Marker-aware (current, mouse) | `src/preprocess_paul15.py`: force-include the full marker union; corrected Paul15 symbols (`Hba-a2`/`Hbb-b1`) + `Alas2` + `Ermap`; `--log1p auto` (Paul15 loads raw, max ~168). |
| Human arm | `src/preprocess_human_marrow.py`: CELLxGENE pipeline (normalize_total + log1p) + 1:1 human orthologs (`config/human_markers.json`) via the `V0B_MARKERS_JSON` override. |

## Executed runs (frozen)

| Run | Data | SAE | Verdict | Frozen bundle sha256 |
|-----|------|-----|---------|----------------------|
| Mouse Gate-0 | Paul15 2730x2012 | 512 latents, l1=0.8, 10 seeds, L0 35.75 | NOT SUPPORTED (reverse; method-dependent) | `172861...` |
| Human Gate-0 | Setty 4142x2032 | 512 latents, l1=0.8, 10 seeds, L0 32.5 | NOT SUPPORTED (reverse replicates) | `cc8159c8...` |

Only pre-registered, sparsity-corrected, >=10-seed runs are ever frozen. Earlier
under-/over-regularized runs (n=5, out-of-band L0) are recorded in `PROJECT_AUDIT.md`
and the `LAB_NOTEBOOK.md` but are explicitly NOT frozen and NOT cited as findings.

## How to recover any prior version

```bash
git log --oneline            # find the commit
git show <sha>:<path>        # view a file at that version
git checkout <sha> -- <path> # restore it (then commit)
```
