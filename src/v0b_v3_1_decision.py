#!/usr/bin/env python3
"""V0b v3.1: control-referenced, count-fair, pre-registered asymmetry decision.

This is the SINGLE pre-registered modularity metric (Component 0, Gate 0). It
fixes the v3 verdict artifact (PROJECT_AUDIT.md H) and locks the decision rule
before the final runs. Changes vs the earlier v3.1 draft:

1. Progenitor-ONLY control baseline. Cycling is dropped from the baseline because
   cell cycle is a real biological program (a poor null); it is still reported as a
   secondary control that should sit near the baseline.
2. Abundance-matched permutation null. Random gene sets are matched to each marker
   set's per-gene decoder mass (binned), so "significant" means a coherent program,
   not merely "high-decoder-mass genes." This addresses the mis-calibrated null in
   which the Cycling control tested significant.
3. Effect-size floor. A submodule is a REAL program only if its excess strength
   over the Progenitor baseline exceeds EXCESS_FLOOR AND its abundance-matched
   permutation q < 0.05.

Decision rule (pre-registered; do not change after freezing):
- Per seed, per lineage: excess[m] = max(0, strength[m] - Progenitor_baseline).
- n_real[lineage] = submodules with excess > EXCESS_FLOOR AND perm_q < 0.05.
- The original claim (erythroid MORE distributed than granulocyte) is SUPPORTED
  only if, in >= 3 of 5 seeds: ery_n_real >= 2 AND ery_n_real >= gran_n_real AND
  ery has higher excess-entropy than gran. A lineage with < 2 real programs cannot
  be "distributed"; it is unified or undetected. The plain finding (which lineage
  has how many real programs) is reported regardless of the SUPPORTED/NOT verdict.

Run:
  python v0b_v3_1_decision.py --data-dir <dir> --output-dir <dir>/v0b_outputs
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

import v0b_module_definitions as v0b
import v0b_v3_loading as v3

CONTROL_SUBMODULE = "Progenitor"      # primary baseline (Cycling dropped: real program)
SECONDARY_CONTROLS = ["Cycling"]      # reported, should sit near baseline
EXCESS_FLOOR = 0.10
PERM_Q = 0.05
N_PERM = 1000
N_ABUND_BINS = 20
SEED = v0b.SEED


def _norm_largest_frac(vals):
    n = len(vals); total = float(sum(vals))
    if n <= 1 or total == 0:
        return np.nan
    frac = max(vals) / total
    return (frac - 1.0 / n) / (1.0 - 1.0 / n)


def abundance_matched_null(decoder_weights, coverage, usable, n_genes, n_perm, rng,
                           n_bins=N_ABUND_BINS):
    """Per-submodule empirical p: is the submodule's mean-across-seeds strength above
    what abundance-matched random gene sets achieve? Matching is on per-gene total
    |decoder weight| (binned), so significance reflects coherent concentration, not
    gene abundance."""
    abund = sum(np.abs(decoder_weights[s]).sum(axis=1) for s in decoder_weights)  # (n_genes,)
    edges = np.quantile(abund, np.linspace(0, 1, n_bins + 1)[1:-1])
    gene_bin = np.digitize(abund, edges)                    # bin index per gene
    genes_by_bin = {b: np.where(gene_bin == b)[0] for b in np.unique(gene_bin)}
    Wabs = {s: np.abs(decoder_weights[s]) for s in decoder_weights}
    total = {s: np.where(Wabs[s].sum(0) == 0, np.nan, Wabs[s].sum(0)) for s in decoder_weights}

    def strength_for(idx, K):
        per_seed = []
        for s in decoder_weights:
            enr = (Wabs[s][idx, :].sum(0) / total[s]) * (n_genes / K)
            per_seed.append(v3._top_mean(enr))
        return float(np.mean(per_seed))

    modules, pvals = [], []
    for m in usable:
        idx = np.asarray(coverage[m]["present_indices"])
        K = len(idx)
        obs = strength_for(idx, K)
        null = np.array([strength_for(
            np.array([rng.choice(genes_by_bin[gene_bin[g]]) for g in idx]), K)
            for _ in range(n_perm)])
        modules.append(m)
        pvals.append((1 + int((null >= obs).sum())) / (1 + n_perm))
    q = multipletests(pvals, method="fdr_bh")[1]
    return {m: {"perm_p": float(p), "perm_q": float(qq)}
            for m, p, qq in zip(modules, pvals, q)}


def per_seed_control_referenced(strength_df, ery_subs, gran_subs, perm, usable):
    real = {m for m in usable if perm.get(m, {}).get("perm_q", 1.0) < PERM_Q}
    rows = []
    for seed in range(v0b.N_SEEDS):
        s = strength_df[strength_df.seed == seed].set_index("module")["strength"]
        baseline = float(s.get(CONTROL_SUBMODULE, 0.0))

        def excess(subs):
            return [max(0.0, float(s.get(m, 0.0)) - baseline) for m in subs]

        def n_real(subs, x):
            return int(sum((xi > EXCESS_FLOOR) and (m in real)
                           for m, xi in zip(subs, x)))

        ery_x, gran_x = excess(ery_subs), excess(gran_subs)
        rows.append({
            "seed": seed, "progenitor_baseline": baseline,
            "cycling_strength": float(s.get("Cycling", np.nan)),
            "ery_excess": ery_x, "gran_excess": gran_x,
            "ery_n_real": n_real(ery_subs, ery_x),
            "gran_n_real": n_real(gran_subs, gran_x),
            "ery_excess_entropy": v0b._norm_entropy(ery_x),
            "gran_excess_entropy": v0b._norm_entropy(gran_x),
            "ery_conc": _norm_largest_frac(ery_x),
            "gran_conc": _norm_largest_frac(gran_x),
        })
    return pd.DataFrame(rows)


def decide(cr_df):
    both = ((cr_df["ery_excess_entropy"] > cr_df["gran_excess_entropy"]) &
            (cr_df["ery_n_real"] >= 2) &
            (cr_df["ery_n_real"] >= cr_df["gran_n_real"]))
    # Pre-registered >=3/5 rule as a >=60% majority (scales to >=10 seeds).
    support_threshold = int(np.ceil(0.6 * len(cr_df)))
    e, g = cr_df["ery_n_real"].median(), cr_df["gran_n_real"].median()
    if e == 0 and g == 0:
        finding = "Neither lineage has a real above-control program (undetected)."
    elif e >= g and e >= 2:
        finding = f"Erythroid distributed across {e:.0f} real programs vs granulocyte {g:.0f}."
    elif g > e:
        finding = (f"Granulocyte has more real programs (median {g:.0f}) than erythroid "
                   f"({e:.0f}); erythroid is unified/undetected. NOT the original claim.")
    else:
        finding = f"Erythroid {e:.0f} real programs vs granulocyte {g:.0f}."
    return {
        "control_baseline": CONTROL_SUBMODULE,
        "median_ery_n_real": float(e), "median_gran_n_real": float(g),
        "n_seeds_supported_pattern": int(both.sum()),
        "support_threshold_60pct": support_threshold,
        "asymmetric_modularity_supported": bool(int(both.sum()) >= support_threshold),
        "plain_finding": finding,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default=os.environ.get("V0B_DATA_DIR", "./data"))
    ap.add_argument("--output-dir", default=os.environ.get("V0B_OUTPUT_DIR", "./v0b_outputs"))
    ap.add_argument("--n-perm", type=int, default=N_PERM)
    args = ap.parse_args()
    out = Path(args.output_dir); out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    data_dir = Path(args.data_dir)
    gene_names = v0b.load_gene_names(data_dir)
    n_genes = len(gene_names)
    v0b.load_expression(data_dir)
    decoder_weights, _ = v0b.load_checkpoints(data_dir, input_dim=n_genes)
    coverage, usable = v0b.compute_coverage(gene_names)
    ery_subs = [m for m in v0b.ERY_SUBMODULES if m in usable]
    gran_subs = [m for m in v0b.GRAN_SUBMODULES if m in usable]

    enrich = v3.loading_enrichment(decoder_weights, coverage, usable, n_genes)
    strength = v3.submodule_strength(enrich, usable)
    perm = abundance_matched_null(decoder_weights, coverage, usable, n_genes, args.n_perm, rng)
    cr = per_seed_control_referenced(strength, ery_subs, gran_subs, perm, usable)
    decision = decide(cr)

    strength_summary = (strength.groupby("module")["strength"].mean()
                        .reset_index().sort_values("strength", ascending=False))
    strength_summary["perm_q"] = strength_summary["module"].map(
        lambda m: perm.get(m, {}).get("perm_q", np.nan))
    cr.to_csv(out / "v0b_v3_1_per_seed.csv", index=False)
    with open(out / "v0b_v3_1_decision.json", "w") as f:
        json.dump({"excess_floor": EXCESS_FLOOR, "perm_q_threshold": PERM_Q,
                   "control_baseline": CONTROL_SUBMODULE, "decision": decision,
                   "submodule": strength_summary.set_index("module").to_dict("index")},
                  f, indent=2)

    print("\n" + "=" * 72)
    print("V0b v3.1 (pre-registered) SUBMODULE STRENGTH + abundance-matched perm q")
    print("=" * 72)
    print(strength_summary.to_string(index=False))
    print("\n" + "=" * 72)
    print("CONTROL-REFERENCED DECISION (Progenitor baseline)")
    print("=" * 72)
    print(cr[["seed", "progenitor_baseline", "ery_n_real", "gran_n_real",
              "ery_excess_entropy", "gran_excess_entropy"]].to_string(index=False))
    for k, v in decision.items():
        print(f"  {k}: {v}")
    print(f"\n  ASYMMETRIC MODULARITY (pre-registered): "
          f"{'SUPPORTED' if decision['asymmetric_modularity_supported'] else 'NOT SUPPORTED'}")
    print(f"  FINDING: {decision['plain_finding']}")


if __name__ == "__main__":
    main()
