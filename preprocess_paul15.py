#!/usr/bin/env python3
"""Marker-aware preprocessing of Paul15 (dataset upgrade).

The original pipeline let generic variance-based HVG selection drop hemoglobins
and late-granulocyte genes. That handicapped the erythroid (and granulocyte
secondary) lineages: the SAE never saw their strongest effector programs, so those
programs could not form features and the lineage looked weak/distributed for a
data reason, not a biological one (PROJECT_AUDIT.md, "code/compute" confounds).

This script force-includes the full canonical marker union (every gene in
v0b_module_definitions.MARKER_SETS, including the globin and late-granulocyte
submodules) on top of the top-variance HVGs, so no lineage is analyzed with its
effector program removed. Everything else matches the project pipeline.

CRITICAL: Paul15 is already log-transformed. Do NOT normalize_total/log1p.

Output (a NEW matrix, deliberately different from the MD5-locked original; it
carries its own recorded MD5):
  expression_matrix.npy   (n_cells x n_genes float32, scaled)
  gene_names.csv          selected genes, MGI symbols
  cell_metadata.csv       cell_id, paul15_clusters
  preprocess_report.json  n_hvg, rescued markers, per-submodule coverage, MD5

Run (Colab, where Paul15 downloads):
  python preprocess_paul15.py --n-hvg 2000 --out-dir ./data_markeraware
Then train (overcomplete): python train_sae.py --matrix .../expression_matrix.npy
  --latent-dim 512 --out-dir .../checkpoints
Then analyze:              python v0b_v3_loading.py --data-dir ...
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

import v0b_module_definitions as v0b


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-hvg", type=int, default=2000)
    ap.add_argument("--out-dir", default="./data_markeraware")
    args = ap.parse_args()

    import scanpy as sc  # lazy: only needed to run, not to import/syntax-check
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    ad = sc.datasets.paul15()                 # pre-log-transformed; do NOT log again
    genes = list(map(str, ad.var_names))
    X = np.asarray(ad.X, dtype=np.float32)    # (n_cells, n_genes_all)
    if float(X.max()) > 50:
        raise ValueError(f"Paul15 max {X.max():.1f} looks un-logged; expected "
                         f"pre-log values. Do not normalize/log this data.")

    gene_to_idx = {g: i for i, g in enumerate(genes)}
    variance = X.var(axis=0)
    top = list(np.argsort(variance)[::-1][:args.n_hvg])

    marker_union = []
    for gs in v0b.MARKER_SETS.values():
        for g in gs:
            if g in gene_to_idx and g not in marker_union:
                marker_union.append(g)
    forced = [gene_to_idx[g] for g in marker_union]
    rescued = [g for g in marker_union if gene_to_idx[g] not in set(top)]

    keep = sorted(set(top) | set(forced))
    kept_genes = [genes[i] for i in keep]
    Xk = X[:, keep]

    mu = Xk.mean(0); sd = Xk.std(0); sd[sd == 0] = 1.0
    Xs = ((Xk - mu) / sd).astype(np.float32)

    np.save(out / "expression_matrix.npy", Xs)
    pd.DataFrame({"gene": kept_genes}).to_csv(out / "gene_names.csv", index=False)
    pd.DataFrame({"cell_id": list(map(str, ad.obs_names)),
                  "paul15_clusters": list(map(str, ad.obs["paul15_clusters"]))}
                 ).to_csv(out / "cell_metadata.csv", index=False)

    md5 = hashlib.md5(Xs.tobytes()).hexdigest()
    coverage = {}
    kept_set = set(kept_genes)
    for m, gs in v0b.MARKER_SETS.items():
        present = [g for g in gs if g in kept_set]
        coverage[m] = {"present": len(present), "requested": len(gs),
                       "present_genes": present}
    report = {"n_cells": int(Xs.shape[0]), "n_genes": int(Xs.shape[1]),
              "n_hvg_requested": args.n_hvg, "expression_matrix_md5": md5,
              "n_markers_forced": len(forced),
              "rescued_markers_not_in_top_variance": rescued,
              "submodule_coverage": coverage}
    with open(out / "preprocess_report.json", "w") as f:
        json.dump(report, f, indent=2)

    print(f"Saved {Xs.shape} matrix. MD5 {md5}")
    print(f"Rescued {len(rescued)} markers that variance-HVG would have dropped: {rescued}")
    for m, c in coverage.items():
        print(f"  {m:14s} {c['present']:2d}/{c['requested']:2d}")


if __name__ == "__main__":
    main()
