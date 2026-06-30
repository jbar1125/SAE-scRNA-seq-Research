#!/usr/bin/env python3
"""V0b v3: continuous decoder-loading test of submodule structure.

Why v3
------
v0b (v2) assigns each feature to one submodule via hypergeometric enrichment +
BH-FDR + winner-take-all. The 2026-06-30 run and diagnostic (PROJECT_AUDIT.md G)
showed that gate discards diluted-but-real signal: erythroid TF/membrane and
granulocyte TF markers co-occur at overlap 2 but never clear q<0.05 across ~896
mostly-null tests, so erythroid scored as absent (0-1 features/seed) and the
asymmetry metric was uninformative. The only thing the gated test detects is
tightly co-concentrated effector programs (Gran_Primary, Ery_Heme).

v3 removes the gate. For each (seed, feature, submodule) it computes a continuous
loading enrichment from the FULL decoder column (no TOP_K, no significance
threshold): the fraction of the feature's total |decoder| mass that lands on the
submodule's markers, divided by the fraction expected from marker-set size
(=1 under the null). Submodule strength = mean of the top-FEAT_TOP features. The
within-lineage distribution is measured on these continuous strengths (normalized
entropy + largest fraction), so 'distributed' signal can register instead of being
thresholded to zero. A gene-label permutation null calibrates each submodule.

Two framings are reported: the original lineage axis (granulocyte vs erythroid)
and the axis the data actually shows (effector vs TF).

Honesty note
------------
This does NOT reinstate the headline. Per the diagnostic, the expectation is that
the lineage asymmetry still will not hold and the effector-vs-TF concentration
result will. The script reports whichever way it comes out.

Run
  python v0b_v3_loading.py --data-dir ./data --output-dir ./v0b_outputs
Needs only gene_names.csv + sae_seed0-4.pt. expression_matrix.npy is optional
(used only for the MD5 non-negotiable).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

import v0b_module_definitions as v0b

FEAT_TOP = 3          # submodule strength = mean of this many top features
N_PERM = 1000         # gene-label permutation null
SEED = v0b.SEED
TF_SUBMODULES = ["Ery_TF", "Gran_TF"]
EFFECTOR_SUBMODULES = ["Ery_Heme", "Ery_Membrane", "Gran_Primary"]


def loading_enrichment(decoder_weights, coverage, usable, n_genes) -> pd.DataFrame:
    """Continuous per (seed, feature, submodule) enrichment of decoder mass.

    enrichment = (mass on submodule markers / total feature mass) / (K / N_genes).
    Uses the full decoder column (no TOP_K, no threshold). 1.0 = null.
    """
    rows = []
    for seed in range(v0b.N_SEEDS):
        Wabs = np.abs(decoder_weights[seed])            # (N_genes, 128)
        total = Wabs.sum(axis=0)                        # (128,)
        total = np.where(total == 0, np.nan, total)
        for module in usable:
            idx = coverage[module]["present_indices"]
            K = len(idx)
            mass = Wabs[idx, :].sum(axis=0)             # (128,)
            enr = (mass / total) * (n_genes / K)        # (128,)
            for feat in range(v0b.LATENT_DIM):
                rows.append({"seed": seed, "feature_idx": feat,
                             "module": module, "enrichment": float(enr[feat])})
    return pd.DataFrame(rows)


def _top_mean(values, k=FEAT_TOP):
    v = np.sort(np.asarray(values, dtype=float))[::-1][:k]
    return float(np.mean(v)) if len(v) else 0.0


def submodule_strength(enrich_df, usable) -> pd.DataFrame:
    rows = []
    for seed in range(v0b.N_SEEDS):
        for module in usable:
            vals = enrich_df[(enrich_df.seed == seed) &
                             (enrich_df.module == module)]["enrichment"].values
            rows.append({"seed": seed, "module": module,
                         "strength": _top_mean(vals)})
    return pd.DataFrame(rows)


def per_seed_distribution(strength_df, ery_subs, gran_subs) -> pd.DataFrame:
    rows = []
    for seed in range(v0b.N_SEEDS):
        s = strength_df[strength_df.seed == seed].set_index("module")["strength"]
        ery = [float(s.get(m, 0.0)) for m in ery_subs]
        gran = [float(s.get(m, 0.0)) for m in gran_subs]
        rows.append({
            "seed": seed,
            "ery_total": sum(ery), "gran_total": sum(gran),
            "ery_norm_entropy": v0b._norm_entropy(ery),
            "gran_norm_entropy": v0b._norm_entropy(gran),
            "ery_largest_frac": v0b._largest_fraction(ery),
            "gran_largest_frac": v0b._largest_fraction(gran),
        })
    return pd.DataFrame(rows)


def permutation_null(decoder_weights, coverage, usable, n_genes, n_perm, rng):
    """Empirical p per submodule: is its mean-across-seeds strength above what a
    random gene set of the same size achieves? One-sided, (1+hits)/(1+n_perm)."""
    Wabs = {s: np.abs(decoder_weights[s]) for s in range(v0b.N_SEEDS)}
    total = {s: np.where(Wabs[s].sum(0) == 0, np.nan, Wabs[s].sum(0))
             for s in range(v0b.N_SEEDS)}

    def strength_for(idx, K):
        per_seed = []
        for s in range(v0b.N_SEEDS):
            enr = (Wabs[s][idx, :].sum(0) / total[s]) * (n_genes / K)
            per_seed.append(_top_mean(enr))
        return float(np.mean(per_seed))

    pvals = {}
    for module in usable:
        idx = coverage[module]["present_indices"]
        K = len(idx)
        obs = strength_for(np.asarray(idx), K)
        null = np.array([strength_for(rng.choice(n_genes, K, replace=False), K)
                         for _ in range(n_perm)])
        pvals[module] = {"observed_strength": obs,
                         "perm_p": float((1 + int((null >= obs).sum())) / (1 + n_perm))}
    return pvals


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default=os.environ.get("V0B_DATA_DIR", "./data"))
    ap.add_argument("--output-dir", default=os.environ.get("V0B_OUTPUT_DIR", "./v0b_outputs"))
    ap.add_argument("--n-perm", type=int, default=N_PERM)
    args = ap.parse_args()

    data_dir = Path(args.data_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    gene_names = v0b.load_gene_names(data_dir)
    n_genes = len(gene_names)
    X = v0b.load_expression(data_dir)              # MD5 check if present
    if X is not None and X.shape[1] != n_genes:
        raise ValueError(f"gene_names ({n_genes}) != expression cols ({X.shape[1]})")
    decoder_weights, _ = v0b.load_checkpoints(data_dir, input_dim=n_genes)

    coverage, usable = v0b.compute_coverage(gene_names)
    ery_subs = [m for m in v0b.ERY_SUBMODULES if m in usable]
    gran_subs = [m for m in v0b.GRAN_SUBMODULES if m in usable]

    enrich = loading_enrichment(decoder_weights, coverage, usable, n_genes)
    strength = submodule_strength(enrich, usable)
    mod_df = per_seed_distribution(strength, ery_subs, gran_subs)
    decision = v0b.decide(mod_df, rng)
    perm = permutation_null(decoder_weights, coverage, usable, n_genes, args.n_perm, rng)

    # attach perm p to strength summary
    strength_summary = (strength.groupby("module")["strength"].mean()
                        .reset_index().rename(columns={"strength": "mean_strength"}))
    strength_summary["perm_p"] = strength_summary["module"].map(
        lambda m: perm[m]["perm_p"])
    strength_summary = strength_summary.sort_values("mean_strength", ascending=False)

    # effector vs TF framing
    def grp_mean(mods):
        mm = [m for m in mods if m in usable]
        return float(strength_summary.set_index("module").loc[mm, "mean_strength"].mean()) if mm else float("nan")
    effector_vs_tf = {
        "effector_mean_strength": grp_mean(EFFECTOR_SUBMODULES),
        "tf_mean_strength": grp_mean(TF_SUBMODULES),
    }

    enrich.to_csv(out_dir / "loading_enrichment_v0b_v3.csv", index=False)
    strength_summary.to_csv(out_dir / "submodule_strength_v0b_v3.csv", index=False)
    mod_df.to_csv(out_dir / "modularity_metrics_v3_per_seed.csv", index=False)

    provenance = {
        "version": "v3_continuous_loading",
        "constants": {"FEAT_TOP": FEAT_TOP, "N_PERM": args.n_perm, "SEED": SEED},
        "expression_matrix_md5_verified": X is not None,
        "n_genes": n_genes, "usable_submodules": usable,
        "submodule_strength": {r["module"]: {"mean_strength": float(r["mean_strength"]),
                                             "perm_p": float(r["perm_p"])}
                               for _, r in strength_summary.iterrows()},
        "effector_vs_tf": effector_vs_tf,
        "lineage_decision": decision,
    }
    with open(out_dir / "v0b_v3_provenance.json", "w") as f:
        json.dump(provenance, f, indent=2)

    print("\n" + "=" * 72)
    print("V0b v3 SUBMODULE STRENGTH (continuous decoder loading; 1.0 = null)")
    print("=" * 72)
    print(strength_summary.to_string(index=False))
    print(f"\nEffector mean strength: {effector_vs_tf['effector_mean_strength']:.3f}"
          f"   TF mean strength: {effector_vs_tf['tf_mean_strength']:.3f}")
    print("\n" + "=" * 72)
    print("LINEAGE AXIS (granulocyte vs erythroid) DECISION")
    print("=" * 72)
    for k, v in decision.items():
        print(f"  {k}: {v}")
    verdict = ("SUPPORTED" if decision["asymmetric_modularity_supported"]
               else "NOT SUPPORTED")
    print(f"\n  LINEAGE ASYMMETRY: {verdict}")
    print("  (Read alongside the effector-vs-TF strengths above; the real "
          "structure may be on that axis, not the lineage axis.)")


if __name__ == "__main__":
    main()
