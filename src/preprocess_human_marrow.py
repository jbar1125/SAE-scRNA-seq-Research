#!/usr/bin/env python3
"""Marker-aware preprocessing of the human CD34+ bone-marrow hematopoiesis atlas
(Setty et al. 2019, distributed with Palantir as marrow_sample_scseq_counts.h5ad).

This is the HUMAN REPLICATION arm of Component 0 (Gate 0). It mirrors
preprocess_paul15.py structurally so the SAE sees the same kind of input (z-scaled
log-expression, marker-aware HVG), with ONE deliberate difference:

  Human data are raw UMI counts (max ~632), so they require library-size
  normalization BEFORE log1p (the standard CELLxGENE pipeline; CLAUDE.md non-
  negotiable "CELLxGENE requires both normalize_total and log1p before HVG").
  Paul15's debatched matrix is already depth-adjusted MARS-seq, so it got log1p
  only. That difference is expected and is recorded here, not hidden.

Marker set: HUMAN orthologs of the exact mouse submodules (human_markers.json),
so the pre-registered v3.1 metric and submodule structure are reused unchanged;
only the gene symbols differ. Run the downstream analysis with
  V0B_MARKERS_JSON=human_markers.json
so v0b_module_definitions overrides its mouse panel with the human orthologs.

Output (a NEW matrix carrying its own recorded MD5, accepted by
v0b.load_expression via preprocess_report.json):
  expression_matrix.npy   (n_cells x n_genes float32, z-scaled log-normalized)
  gene_names.csv          selected genes, HGNC symbols
  cell_metadata.csv       cell_id
  preprocess_report.json  n_hvg, rescued markers, per-submodule coverage, MD5

Run:
  python preprocess_human_marrow.py \
      --h5ad data_g0_human/marrow_sample_scseq_counts.h5ad \
      --markers human_markers.json --n-hvg 2000 --out-dir ./data_g0_human
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5ad", default="data_g0_human/marrow_sample_scseq_counts.h5ad")
    ap.add_argument("--markers", default="config/human_markers.json")
    ap.add_argument("--n-hvg", type=int, default=2000)
    ap.add_argument("--target-sum", type=float, default=1e4)
    ap.add_argument("--out-dir", default="./data_g0_human")
    args = ap.parse_args()

    import anndata as ad
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    A = ad.read_h5ad(args.h5ad)
    genes = list(map(str, A.var_names))
    X = np.asarray(A.X.todense() if hasattr(A.X, "todense") else A.X, dtype=np.float32)
    raw_max = float(X.max())

    # --- CELLxGENE pipeline: library-size normalize THEN log1p ---
    lib = X.sum(axis=1, keepdims=True)
    lib[lib == 0] = 1.0
    X = (X / lib) * args.target_sum
    X = np.log1p(X).astype(np.float32)
    print(f"Normalized to target_sum {args.target_sum:g} + log1p (raw max was {raw_max:.1f}).")

    marker_spec = json.load(open(args.markers))
    marker_sets = marker_spec["marker_sets"]

    gene_to_idx = {g: i for i, g in enumerate(genes)}
    variance = X.var(axis=0)
    top = list(np.argsort(variance)[::-1][:args.n_hvg])

    marker_union = []
    for gs in marker_sets.values():
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
    pd.DataFrame({"cell_id": list(map(str, A.obs_names))}
                 ).to_csv(out / "cell_metadata.csv", index=False)

    md5 = hashlib.md5(Xs.tobytes()).hexdigest()
    kept_set = set(kept_genes)
    coverage = {}
    for m, gs in marker_sets.items():
        present = [g for g in gs if g in kept_set]
        coverage[m] = {"present": len(present), "requested": len(gs),
                       "present_genes": present}
    report = {"species": "human", "source": "Setty2019 CD34+ marrow (Palantir)",
              "n_cells": int(Xs.shape[0]), "n_genes": int(Xs.shape[1]),
              "n_hvg_requested": args.n_hvg, "expression_matrix_md5": md5,
              "normalize_total_target_sum": args.target_sum,
              "log1p_applied": True, "raw_max_before_norm": raw_max,
              "n_markers_forced": len(forced),
              "rescued_markers_not_in_top_variance": rescued,
              "submodule_coverage": coverage}
    with open(out / "preprocess_report.json", "w") as f:
        json.dump(report, f, indent=2)

    print(f"Saved {Xs.shape} matrix. MD5 {md5}")
    print(f"Rescued {len(rescued)} markers variance-HVG would have dropped: {rescued}")
    for m, c in coverage.items():
        print(f"  {m:16s} {c['present']:2d}/{c['requested']:2d}")


if __name__ == "__main__":
    main()
