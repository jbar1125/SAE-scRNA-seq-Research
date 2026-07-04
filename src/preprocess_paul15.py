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

NOTE ON LOG-TRANSFORM: the handoff calls Paul15 "pre-log-transformed", but
sc.datasets.paul15() actually returns RAW-scale data (max ~168). That claim is
also inconsistent with hemoglobins being HVG-filtered (only happens on log data).
So this script log1p's raw data by default (--log1p auto). See PROJECT_AUDIT.md.

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


def _load_paul15(sc):
    """Load Paul15 via scanpy; fall back to reading the local h5 directly if the
    installed scanpy parser is incompatible with the on-disk file version (the
    file layout has changed over scanpy releases). The local file is scanpy's
    datasetdir/paul15/paul15.h5."""
    try:
        return sc.datasets.paul15()
    except Exception as e:
        import h5py
        import anndata as ad_mod
        path = Path(sc.settings.datasetdir) / "paul15" / "paul15.h5"
        if not path.exists():
            raise
        print(f"scanpy paul15 parser failed ({type(e).__name__}); reading {path} directly.")
        with h5py.File(path, "r") as f:
            X = np.asarray(f["data.debatched"][()], dtype=np.float32)      # (cells, genes)
            rown = f["data.debatched_rownames"][()].astype(str)
            genes = [r.split(";")[0] for r in rown]                         # primary symbol
            clusters = np.asarray(f["cluster.id"][()]).reshape(-1).astype(int)
        # drop duplicate gene symbols (keep first occurrence)
        seen, keep, kept_genes = set(), [], []
        for i, g in enumerate(genes):
            if g not in seen:
                seen.add(g); keep.append(i); kept_genes.append(g)
        A = ad_mod.AnnData(X[:, keep])
        A.var_names = kept_genes
        A.obs["paul15_clusters"] = [str(c) for c in clusters]
        return A


def apply_log_policy(X, mode, raw_threshold=30.0):
    """Decide whether to log1p. sc.datasets.paul15() returns RAW-scale data
    (max ~168), contradicting the handoff's 'pre-log-transformed' claim. That
    claim is also inconsistent with hemoglobins being HVG-filtered, which only
    happens on log-scale data. So 'auto' log1p's when the data looks raw.
    Returns (X_out, did_log, raw_max)."""
    raw_max = float(X.max())
    do_log = (mode == "yes") or (mode == "auto" and raw_max > raw_threshold)
    if do_log:
        X = np.log1p(X)
    return X, do_log, raw_max


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-hvg", type=int, default=2000)
    ap.add_argument("--out-dir", default="./data_markeraware")
    ap.add_argument("--log1p", choices=["auto", "yes", "no"], default="auto",
                    help="auto: log1p when data looks raw (max>30). See PROJECT_AUDIT.md.")
    args = ap.parse_args()

    import scanpy as sc  # lazy: only needed to run, not to import/syntax-check
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    ad = _load_paul15(sc)
    genes = list(map(str, ad.var_names))
    X = np.asarray(ad.X, dtype=np.float32)    # (n_cells, n_genes_all)
    X, did_log, raw_max = apply_log_policy(X, args.log1p)
    if did_log:
        print(f"Applied log1p (raw max was {raw_max:.1f}). NOTE: handoff claims "
              f"Paul15 is pre-log-transformed, but this source is raw-scale; "
              f"see PROJECT_AUDIT.md.")
    else:
        print(f"Skipped log1p (raw max {raw_max:.1f}, --log1p={args.log1p}).")

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
              "log1p_applied": bool(did_log), "raw_max_before_log": raw_max,
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
