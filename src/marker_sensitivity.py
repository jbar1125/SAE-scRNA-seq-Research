#!/usr/bin/env python3
"""Component 0: leave-one-marker-out sensitivity of the v3.1 verdict.

Addresses the "did you cherry-pick markers?" reviewer concern (made concrete by the
globin-symbol bug). For every present marker gene, drop it from its submodule,
re-run the pre-registered v3.1 decision, and record whether the program counts and
the SUPPORTED/NOT verdict change. A robust conclusion should not flip on any single
marker.

Run:
  python marker_sensitivity.py --data-dir <dir> --out <dir>/sensitivity --n-perm 200
Needs gene_names.csv + sae_seed*.pt (expression_matrix.npy optional for MD5).
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

import numpy as np
import pandas as pd

import v0b_module_definitions as v0b
import v0b_v3_loading as v3
import v0b_v3_1_decision as v31


def _decision_for_coverage(decoder_weights, coverage, n_genes, n_perm, rng):
    usable = [m for m, c in coverage.items() if c["present"] >= v0b.MIN_MARKERS_PRESENT]
    ery = [m for m in v0b.ERY_SUBMODULES if m in usable]
    gran = [m for m in v0b.GRAN_SUBMODULES if m in usable]
    enrich = v3.loading_enrichment(decoder_weights, coverage, usable, n_genes)
    strength = v3.submodule_strength(enrich, usable)
    perm = v31.abundance_matched_null(decoder_weights, coverage, usable, n_genes, n_perm, rng)
    cr = v31.per_seed_control_referenced(strength, ery, gran, perm, usable)
    return v31.decide(cr)


def _drop_marker(coverage, module, gene_idx):
    cov = copy.deepcopy(coverage)
    present_idx = cov[module]["present_indices"]
    keep = [i for i in present_idx if i != gene_idx]
    cov[module]["present_indices"] = keep
    cov[module]["present"] = len(keep)
    return cov


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default="./data")
    ap.add_argument("--out", default="./sensitivity")
    ap.add_argument("--n-perm", type=int, default=200)
    args = ap.parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)

    data_dir = Path(args.data_dir)
    gene_names = v0b.load_gene_names(data_dir)
    n_genes = len(gene_names)
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    v0b.load_expression(data_dir)
    decoder_weights, _ = v0b.load_checkpoints(data_dir, input_dim=n_genes)
    coverage, _ = v0b.compute_coverage(gene_names)

    rng = np.random.default_rng(v0b.SEED)
    base = _decision_for_coverage(decoder_weights, coverage, n_genes, args.n_perm, rng)

    rows = []
    lineage_mods = [m for m in v0b.ERY_SUBMODULES + v0b.GRAN_SUBMODULES
                    if coverage[m]["present"] >= v0b.MIN_MARKERS_PRESENT]
    for module in lineage_mods:
        for gene in coverage[module]["present_genes"]:
            cov = _drop_marker(coverage, module, gene_to_idx[gene])
            if cov[module]["present"] < v0b.MIN_MARKERS_PRESENT:
                rows.append({"module": module, "dropped": gene, "note": "submodule dropped",
                             "supported": None, "ery_n_real": None, "gran_n_real": None})
                continue
            d = _decision_for_coverage(decoder_weights, cov, n_genes, args.n_perm,
                                       np.random.default_rng(v0b.SEED))
            rows.append({"module": module, "dropped": gene,
                         "supported": d["asymmetric_modularity_supported"],
                         "ery_n_real": d["median_ery_n_real"],
                         "gran_n_real": d["median_gran_n_real"],
                         "flips_verdict": d["asymmetric_modularity_supported"]
                                          != base["asymmetric_modularity_supported"]})
    df = pd.DataFrame(rows)
    df.to_csv(out / "marker_sensitivity.csv", index=False)
    n_flip = int(df.get("flips_verdict", pd.Series(dtype=bool)).sum())
    summary = {"baseline_supported": base["asymmetric_modularity_supported"],
               "baseline_finding": base["plain_finding"],
               "n_markers_tested": int(df["dropped"].notna().sum()),
               "n_markers_that_flip_verdict": n_flip,
               "verdict_robust_to_single_marker": n_flip == 0}
    with open(out / "marker_sensitivity_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("Baseline:", base["plain_finding"])
    print(f"Markers tested: {summary['n_markers_tested']}; "
          f"verdict flips on {n_flip}; robust={summary['verdict_robust_to_single_marker']}")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
