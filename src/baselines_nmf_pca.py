#!/usr/bin/env python3
"""Component 0: PCA / NMF baselines (V11 participation ratio + V12 modularity).

Answers two Gate-0 questions:
1. Participation ratio (effective dimensionality) of the data/PCA, NMF, and SAE
   representations, with ONE documented formula, resolving the 43.4-vs-42.3 (PCA)
   and 1.8x-vs-34.94 (NMF) conflicts recorded in PROJECT_AUDIT.md (A1, A2).
2. Whether the v3.1 modularity structure is specific to the SAE, or also appears in
   PCA / NMF components (the "does the SAE beat matrix factorization" question every
   ML reviewer asks). Runs the identical v3.1 loading pipeline on PCA and NMF
   components.

Participation ratio (single formula, document this):
   PR = (sum_i lambda_i)^2 / sum_i lambda_i^2
where lambda_i are the eigenvalues of the covariance of the representation
(equivalently s_i^2 for singular values s_i of the centered matrix). PR is the
effective number of dimensions: PR = D for isotropic data, PR = 1 for rank-1.

Run:
  python baselines_nmf_pca.py --data-dir <dir> --rank 512 --out <dir>/baselines
Needs expression_matrix.npy + gene_names.csv. If SAE checkpoints are present, the
SAE participation ratio and modularity are computed for direct comparison.
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


def participation_ratio(eigs) -> float:
    """PR = (sum eig)^2 / sum(eig^2). eigs are non-negative variances."""
    eigs = np.asarray(eigs, dtype=float)
    eigs = eigs[eigs > 0]
    if eigs.size == 0:
        return float("nan")
    return float((eigs.sum() ** 2) / (eigs ** 2).sum())


def pr_of_representation(Z) -> float:
    """PR of the covariance eigenvalues of a (samples x features) representation."""
    Zc = Z - Z.mean(axis=0, keepdims=True)
    # eigenvalues of covariance = singular values^2 of centered Z
    s = np.linalg.svd(Zc, compute_uv=False)
    return participation_ratio(s ** 2)


def _modularity_on_components(components, coverage, usable, n_genes, ery, gran, rng,
                              n_perm=1000):
    """Run the v3.1 loading pipeline treating `components` (dict seed -> (n_genes, R))
    as a decoder. Temporarily sets N_SEEDS to the number of provided seeds."""
    old = v0b.N_SEEDS
    v0b.N_SEEDS = len(components)
    try:
        enrich = v3.loading_enrichment(components, coverage, usable, n_genes)
        strength = v3.submodule_strength(enrich, usable)
        perm = v31.abundance_matched_null(components, coverage, usable, n_genes, n_perm, rng)
        cr = v31.per_seed_control_referenced(strength, ery, gran, perm, usable)
        decision = v31.decide(cr)
        smry = strength.groupby("module")["strength"].mean().to_dict()
        return {"decision": decision, "mean_strength": smry,
                "perm_q": {m: perm[m]["perm_q"] for m in usable}}
    finally:
        v0b.N_SEEDS = old


def main():
    from sklearn.decomposition import PCA, NMF

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default="./data")
    ap.add_argument("--rank", type=int, default=512)
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--out", default="./baselines")
    args = ap.parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(v0b.SEED)

    data_dir = Path(args.data_dir)
    gene_names = v0b.load_gene_names(data_dir)
    n_genes = len(gene_names)
    X = v0b.load_expression(data_dir)
    if X is None:
        raise SystemExit("expression_matrix.npy required for baselines.")
    coverage, usable = v0b.compute_coverage(gene_names)
    ery = [m for m in v0b.ERY_SUBMODULES if m in usable]
    gran = [m for m in v0b.GRAN_SUBMODULES if m in usable]

    results = {"rank": args.rank, "participation_ratio": {}, "modularity": {}}

    # --- Participation ratios (one formula) ---
    results["participation_ratio"]["data_pca"] = pr_of_representation(X)  # PCA/data PR
    # NMF needs non-negative input; use the non-negative part of the (scaled) matrix
    Xnn = np.maximum(X, 0.0)
    nmf = NMF(n_components=args.rank, init="nndsvda", random_state=v0b.SEED, max_iter=400)
    W_nmf = nmf.fit_transform(Xnn)                    # (cells, rank) activations
    results["participation_ratio"]["nmf"] = pr_of_representation(W_nmf)
    results["participation_ratio"]["nmf_note"] = ("NMF fit on max(X,0); ideally the "
                                                  "log-normalized pre-scaling matrix.")

    # SAE PR + SAE modularity if checkpoints present
    sae_present = (data_dir / "sae_seed0.pt").exists()
    if sae_present:
        import torch
        dw, models = v0b.load_checkpoints(data_dir, input_dim=n_genes)
        with torch.no_grad():
            Z = models[0](torch.tensor(X, dtype=torch.float32))[1].numpy()
        results["participation_ratio"]["sae"] = pr_of_representation(Z)

    # --- Modularity on PCA / NMF components (V12) ---
    pca = PCA(n_components=args.rank, random_state=v0b.SEED).fit(X)
    pca_dec = {s: pca.components_.T.copy() for s in range(args.seeds)}    # (n_genes, rank)
    results["modularity"]["pca"] = _modularity_on_components(
        pca_dec, coverage, usable, n_genes, ery, gran, rng)

    nmf_dec = {}
    for s in range(args.seeds):
        m = NMF(n_components=args.rank, init="nndsvda", random_state=s, max_iter=400).fit(Xnn)
        nmf_dec[s] = m.components_.T.copy()
    results["modularity"]["nmf"] = _modularity_on_components(
        nmf_dec, coverage, usable, n_genes, ery, gran, rng)

    if sae_present:
        results["modularity"]["sae"] = _modularity_on_components(
            dw, coverage, usable, n_genes, ery, gran, rng)

    with open(out / "baselines_results.json", "w") as f:
        json.dump(results, f, indent=2, default=float)

    print("\n" + "=" * 72)
    print("PARTICIPATION RATIO (PR = (sum eig)^2 / sum eig^2)")
    print("=" * 72)
    for k, v in results["participation_ratio"].items():
        if isinstance(v, float):
            print(f"  {k:10s}: {v:.2f}")
    print("\n" + "=" * 72)
    print("MODULARITY BY REPRESENTATION (pre-registered v3.1 decision)")
    print("=" * 72)
    for rep, r in results["modularity"].items():
        d = r["decision"]
        print(f"  {rep:5s}: ery_n_real={d['median_ery_n_real']:.0f} "
              f"gran_n_real={d['median_gran_n_real']:.0f} "
              f"supported={d['asymmetric_modularity_supported']}")
    print(f"\nSaved: {out / 'baselines_results.json'}")


if __name__ == "__main__":
    main()
