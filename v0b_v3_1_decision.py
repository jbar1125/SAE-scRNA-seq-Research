#!/usr/bin/env python3
"""V0b v3.1: control-referenced, count-fair asymmetry decision.

v3 removed the significance gate but its verdict was still an artifact
(PROJECT_AUDIT.md H): the `largest_fraction` criterion has a 1/n count bias, and
"distribution" was measured on submodules sitting at the noise floor. v3.1 fixes
the decision, reusing v3's continuous loading strengths:

1. Control baseline. Per seed, baseline = mean strength of the control submodules
   (Progenitor, Cycling). Real programs must clear this noise floor.
2. Excess strength = max(0, strength - baseline). Submodules at/below control
   contribute 0, so noise-floor submodules can no longer masquerade as
   "distribution".
3. n_real = submodules per lineage with excess > EXCESS_FLOOR. This is the honest
   count of detected programs.
4. Count-fair concentration = (largest_frac(excess) - 1/n)/(1 - 1/n), so uniform
   maps to 0 and the metric is comparable across lineages with different submodule
   counts. Normalized entropy is also on excess.
5. Decision. The original claim (erythroid MORE distributed than granulocyte) is
   SUPPORTED only if erythroid actually has real programs to distribute:
   ery_n_real >= 2 AND ery_n_real >= gran_n_real AND ery more distributed
   (higher excess entropy) in >= 3/5 seeds. A lineage with 0 real programs cannot
   be "distributed"; it is undetected.

Reports the plain finding regardless of verdict: which lineage has real
above-control programs, and how many.

Run (after training on marker-aware data):
  python v0b_v3_1_decision.py --data-dir /content/drive/MyDrive/data_ma \
      --output-dir /content/drive/MyDrive/data_ma/v0b_outputs
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

import v0b_module_definitions as v0b
import v0b_v3_loading as v3

CONTROL_SUBMODULES = ["Progenitor", "Cycling"]
EXCESS_FLOOR = 0.10          # excess strength above control to count as a real program
SEED = v0b.SEED


def _norm_largest_frac(vals):
    """Count-normalized concentration: uniform -> 0, fully concentrated -> 1."""
    vals = [v for v in vals]
    n = len(vals)
    total = float(sum(vals))
    if n <= 1 or total == 0:
        return np.nan
    frac = max(vals) / total
    return (frac - 1.0 / n) / (1.0 - 1.0 / n)


def per_seed_control_referenced(strength_df, ery_subs, gran_subs, control_subs):
    rows = []
    for seed in range(v0b.N_SEEDS):
        s = strength_df[strength_df.seed == seed].set_index("module")["strength"]
        baseline = float(np.mean([s.get(m, np.nan) for m in control_subs])) if control_subs else 0.0

        def excess(subs):
            return [max(0.0, float(s.get(m, 0.0)) - baseline) for m in subs]

        ery_x, gran_x = excess(ery_subs), excess(gran_subs)
        rows.append({
            "seed": seed, "control_baseline": baseline,
            "ery_excess": ery_x, "gran_excess": gran_x,
            "ery_n_real": int(sum(x > EXCESS_FLOOR for x in ery_x)),
            "gran_n_real": int(sum(x > EXCESS_FLOOR for x in gran_x)),
            "ery_norm_entropy": v0b._norm_entropy(ery_x),
            "gran_norm_entropy": v0b._norm_entropy(gran_x),
            "ery_conc": _norm_largest_frac(ery_x),
            "gran_conc": _norm_largest_frac(gran_x),
        })
    return pd.DataFrame(rows)


def decide(cr_df):
    ery_more_distributed = (cr_df["ery_norm_entropy"] > cr_df["gran_norm_entropy"])
    ery_has_programs = (cr_df["ery_n_real"] >= 2)
    ery_ge_gran = (cr_df["ery_n_real"] >= cr_df["gran_n_real"])
    both = ery_more_distributed & ery_has_programs & ery_ge_gran
    supported = int(both.sum()) >= 3
    return {
        "median_ery_n_real": float(cr_df["ery_n_real"].median()),
        "median_gran_n_real": float(cr_df["gran_n_real"].median()),
        "n_seeds_ery_has>=2_programs": int(ery_has_programs.sum()),
        "n_seeds_supported_pattern": int(both.sum()),
        "asymmetric_modularity_supported": bool(supported),
        "plain_finding": None,   # filled below
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default=os.environ.get("V0B_DATA_DIR", "./data"))
    ap.add_argument("--output-dir", default=os.environ.get("V0B_OUTPUT_DIR", "./v0b_outputs"))
    args = ap.parse_args()
    out = Path(args.output_dir); out.mkdir(parents=True, exist_ok=True)

    data_dir = Path(args.data_dir)
    gene_names = v0b.load_gene_names(data_dir)
    n_genes = len(gene_names)
    v0b.load_expression(data_dir)                       # MD5 check if present
    decoder_weights, _ = v0b.load_checkpoints(data_dir, input_dim=n_genes)
    coverage, usable = v0b.compute_coverage(gene_names)

    ery_subs = [m for m in v0b.ERY_SUBMODULES if m in usable]
    gran_subs = [m for m in v0b.GRAN_SUBMODULES if m in usable]
    control_subs = [m for m in CONTROL_SUBMODULES if m in usable]

    enrich = v3.loading_enrichment(decoder_weights, coverage, usable, n_genes)
    strength = v3.submodule_strength(enrich, usable)
    cr = per_seed_control_referenced(strength, ery_subs, gran_subs, control_subs)
    decision = decide(cr)

    e_real, g_real = decision["median_ery_n_real"], decision["median_gran_n_real"]
    if e_real == 0 and g_real == 0:
        finding = "Neither lineage has an above-control program (undetected)."
    elif e_real == 0:
        finding = (f"Only granulocyte has real above-control programs "
                   f"(median {g_real:.0f}); erythroid is undetected. This is NOT "
                   f"the 'erythroid distributed' claim.")
    elif g_real == 0:
        finding = (f"Only erythroid has real programs (median {e_real:.0f}); "
                   f"granulocyte undetected.")
    else:
        finding = (f"Both lineages detected: erythroid median {e_real:.0f} vs "
                   f"granulocyte {g_real:.0f} above-control programs.")
    decision["plain_finding"] = finding

    strength_summary = (strength.groupby("module")["strength"].mean()
                        .reset_index().sort_values("strength", ascending=False))
    cr.to_csv(out / "v0b_v3_1_per_seed.csv", index=False)
    with open(out / "v0b_v3_1_decision.json", "w") as f:
        json.dump({"control_submodules": control_subs, "excess_floor": EXCESS_FLOOR,
                   "decision": decision,
                   "mean_strength": strength_summary.set_index("module")["strength"].to_dict()},
                  f, indent=2)

    print("\n" + "=" * 72)
    print("V0b v3.1 CONTROL-REFERENCED DECISION")
    print("=" * 72)
    print(cr[["seed", "control_baseline", "ery_n_real", "gran_n_real",
              "ery_norm_entropy", "gran_norm_entropy"]].to_string(index=False))
    print("\n  above-control programs (median): "
          f"erythroid {e_real:.0f}, granulocyte {g_real:.0f}")
    for k, v in decision.items():
        print(f"  {k}: {v}")
    print(f"\n  ASYMMETRIC MODULARITY (control-referenced): "
          f"{'SUPPORTED' if decision['asymmetric_modularity_supported'] else 'NOT SUPPORTED'}")
    print(f"  FINDING: {finding}")


if __name__ == "__main__":
    main()
