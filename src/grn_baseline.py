#!/usr/bin/env python3
"""Component 0 hardening: a co-expression GRN as a 4th decomposition family in the
method-dependence comparison (PCA, NMF, SAE, and now GRN).

This is the executable stand-in for the SCENIC baseline (full SCENIC is infeasible
in-container: the cisTarget motif-ranking databases are proxy-blocked). It is a
FAST, deterministic, CPU-friendly co-expression regulon network, and it is labeled
as exactly that, NOT as motif-pruned SCENIC:

  For each transcription factor T (aertslab TF list) present in the gene panel, its
  regulon-loading over genes = |Pearson correlation| of each gene with T across cells.
  This yields a genes x n_TF loading matrix, treated as a decomposition (components =
  TFs). The IDENTICAL v3.1 modularity pipeline used for PCA/NMF/SAE
  (baselines_nmf_pca._modularity_on_components) is then run on it, plus its
  participation ratio.

Purpose: if "gene-program modularity" were a property of the data, a GRN-style
decomposition should agree with the others. It does not have to, and the point of
the artifact-vs-signal thesis is precisely to show the answer is decomposition-
dependent. A correlation regulon is a legitimate (if lightweight) member of the
decomposition zoo; its motif-free nature is a stated limitation.

Run:
  python3 src/grn_baseline.py --data-dir ./data_g0 --tfs config/tf_lists/mm_mgi_tfs.txt \
      --out ./data_g0/grn
  V0B_MARKERS_JSON=config/human_markers.json python3 src/grn_baseline.py \
      --data-dir ./data_g0_human --tfs config/tf_lists/hs_hgnc_tfs.txt --out ./data_g0_human/grn
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import v0b_module_definitions as v0b
import baselines_nmf_pca as base


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default="./data_g0")
    ap.add_argument("--tfs", default="config/tf_lists/mm_mgi_tfs.txt")
    ap.add_argument("--out", default="./data_g0/grn")
    ap.add_argument("--n-perm", type=int, default=1000)
    args = ap.parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(v0b.SEED)

    data_dir = Path(args.data_dir)
    gene_names = v0b.load_gene_names(data_dir)
    n_genes = len(gene_names)
    X = v0b.load_expression(data_dir)
    if X is None:
        raise SystemExit("expression_matrix.npy required for the GRN baseline.")
    coverage, usable = v0b.compute_coverage(gene_names)     # respects V0B_MARKERS_JSON
    ery = [m for m in v0b.ERY_SUBMODULES if m in usable]
    gran = [m for m in v0b.GRAN_SUBMODULES if m in usable]

    tf_set = {t.strip() for t in open(args.tfs) if t.strip()}
    g2i = {g: i for i, g in enumerate(gene_names)}
    tf_idx = np.array([g2i[t] for t in gene_names if t in tf_set and t in g2i])
    tf_names = [gene_names[i] for i in tf_idx]
    if len(tf_idx) < 10:
        raise SystemExit(f"only {len(tf_idx)} panel TFs matched {args.tfs}; check symbols.")

    # correlation loading: W[gene, tf] = |corr(gene, tf)| across cells.
    Xc = X - X.mean(axis=0, keepdims=True)
    sd = X.std(axis=0, keepdims=True); sd[sd == 0] = 1.0
    Xn = Xc / sd
    corr = (Xn.T @ Xn[:, tf_idx]) / (X.shape[0] - 1)       # (n_genes, n_tf)
    W = np.abs(corr).astype(float)
    components = {0: W}                                     # one deterministic "seed"

    modularity = base._modularity_on_components(
        components, coverage, usable, n_genes, ery, gran, rng, n_perm=args.n_perm)
    pr = {"data": base.pr_of_representation(X),
          "grn_tf_activity": base.pr_of_representation(X[:, tf_idx])}

    result = {"method": "co-expression GRN (|Pearson| regulons; NOT motif-pruned SCENIC)",
              "n_tfs_used": int(len(tf_idx)), "n_genes": n_genes,
              "participation_ratio": pr,
              "modularity": modularity}
    with open(out / "grn_baseline.json", "w") as f:
        json.dump(result, f, indent=2)

    d = modularity["decision"]
    print(f"GRN baseline ({len(tf_idx)} TFs): PR(data)={pr['data']:.1f} "
          f"PR(TF-activity)={pr['grn_tf_activity']:.1f}")
    print(f"  modularity: ery_n_real={d['median_ery_n_real']} "
          f"gran_n_real={d['median_gran_n_real']} supported={d['asymmetric_modularity_supported']}")
    print(f"  finding: {d['plain_finding']}")


if __name__ == "__main__":
    main()
