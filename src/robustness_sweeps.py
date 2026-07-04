#!/usr/bin/env python3
"""Component 0 hardening: robustness of the frozen v3.1 verdict to its free knobs.

R1 (this script): decision-parameter sensitivity. The pre-registered v3.1 verdict
uses two thresholds that a skeptic could call arbitrary: the effect-size floor
(EXCESS_FLOOR, default 0.10) and the significance threshold (PERM_Q, default 0.05).
If the NOT-SUPPORTED verdict and the granulocyte-more-distributed direction only
hold at one lucky choice of these, the result is a knob-tuning artifact. This sweeps
a grid of both and reports the verdict at every point.

Efficiency: the expensive abundance-matched permutation null depends only on the
data + N_PERM + N_ABUND_BINS, NOT on the two thresholds. So it is computed ONCE per
dataset and the whole threshold grid is evaluated cheaply on top.

Run (from repo root):
  python3 src/robustness_sweeps.py --data-dir ./data_g0        --out ./data_g0/robustness
  V0B_MARKERS_JSON=config/human_markers.json \
    python3 src/robustness_sweeps.py --data-dir ./data_g0_human --out ./data_g0_human/robustness
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import v0b_module_definitions as v0b
import v0b_v3_loading as v3
import v0b_v3_1_decision as v31

EXCESS_FLOORS = [0.05, 0.10, 0.15, 0.20]
PERM_QS = [0.01, 0.05, 0.10]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default="./data_g0")
    ap.add_argument("--out", default="./data_g0/robustness")
    ap.add_argument("--n-perm", type=int, default=500)
    args = ap.parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(v0b.SEED)

    data_dir = Path(args.data_dir)
    gene_names = v0b.load_gene_names(data_dir)
    n_genes = len(gene_names)
    v0b.load_expression(data_dir)
    dw, _ = v0b.load_checkpoints(data_dir, input_dim=n_genes)
    coverage, usable = v0b.compute_coverage(gene_names)
    ery = [m for m in v0b.ERY_SUBMODULES if m in usable]
    gran = [m for m in v0b.GRAN_SUBMODULES if m in usable]

    # expensive parts computed ONCE
    enrich = v3.loading_enrichment(dw, coverage, usable, n_genes)
    strength = v3.submodule_strength(enrich, usable)
    perm = v31.abundance_matched_null(dw, coverage, usable, n_genes, args.n_perm, rng)

    orig_floor, orig_q = v31.EXCESS_FLOOR, v31.PERM_Q
    rows = []
    for floor in EXCESS_FLOORS:
        for q in PERM_QS:
            v31.EXCESS_FLOOR, v31.PERM_Q = floor, q          # read at call time
            cr = v31.per_seed_control_referenced(strength, ery, gran, perm, usable)
            d = v31.decide(cr)
            rows.append({
                "excess_floor": floor, "perm_q": q,
                "median_ery_n_real": d["median_ery_n_real"],
                "median_gran_n_real": d["median_gran_n_real"],
                "gran_ge_ery": d["median_gran_n_real"] >= d["median_ery_n_real"],
                "n_seeds_supported_pattern": d["n_seeds_supported_pattern"],
                "asymmetric_modularity_supported": d["asymmetric_modularity_supported"],
            })
    v31.EXCESS_FLOOR, v31.PERM_Q = orig_floor, orig_q         # restore

    df = pd.DataFrame(rows)
    n = len(df)
    any_supported = bool(df["asymmetric_modularity_supported"].any())
    gran_ge_ery_all = bool(df["gran_ge_ery"].all())
    summary = {
        "dataset": str(data_dir),
        "grid_points": n,
        "excess_floors": EXCESS_FLOORS, "perm_qs": PERM_QS, "n_perm": args.n_perm,
        "any_config_supports_original_claim": any_supported,
        "granulocyte_ge_erythroid_in_all_configs": gran_ge_ery_all,
        "verdict_robust_to_thresholds": (not any_supported) and gran_ge_ery_all,
    }
    df.to_csv(out / "robustness_decision_params.csv", index=False)
    with open(out / "robustness_decision_params_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("\n=== R1 decision-parameter sensitivity ===")
    print(df.to_string(index=False))
    print("\nsummary:", json.dumps(summary, indent=2))
    print(f"\nVERDICT ROBUST TO THRESHOLDS: {summary['verdict_robust_to_thresholds']} "
          f"(no config supports the original claim; granulocyte >= erythroid "
          f"real-program count at every one of {n} grid points)")


if __name__ == "__main__":
    main()
