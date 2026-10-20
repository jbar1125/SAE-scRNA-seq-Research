#!/usr/bin/env python3
"""Component 0 hardening: Independent Component Analysis (FastICA) as a 5th
decomposition family in the method-dependence comparison (PCA, NMF, SAE, GRN, ICA).

ICA is a standard, widely-used single-cell program-extraction method (it finds
statistically independent gene loadings rather than variance-maximizing (PCA),
parts-based (NMF), sparse (SAE), or correlation (GRN) ones). If "gene-program
modularity" were a property of the data, ICA should agree with the others. It does
not have to; the artifact-vs-signal thesis is precisely that the answer is
decomposition-dependent. Adding ICA broadens that demonstration from four families
to five.

Each FastICA component's gene loadings (|components_.T|) form the decoder analog fed
to the IDENTICAL v3.1 modularity pipeline used for PCA/NMF/SAE/GRN
(baselines_nmf_pca._modularity_on_components). A few random_states provide the
metric's seeds.

RANK NOTE (honesty): FastICA at rank 512 (matching the SAE dictionary) does NOT
converge on this data (0-1 of 3 seeds reach a fixed point in max_iter), and
non-converged ICA components are unreliable. 512 independent components is an
ill-posed ICA problem; ICA is standardly run with tens of components. So ICA is run
at rank 50, where it converges 3/3 robustly. This is a different rank from the other
families (512), which is stated, not hidden; the modularity metric counts real
programs among the marker submodules and does not require a matched rank. The
per-seed convergence flags are recorded in the output.

Run:
  python3 src/ica_baseline.py --data-dir ./data_g0 --out ./data_g0/ica
  V0B_MARKERS_JSON=config/human_markers.json python3 src/ica_baseline.py \
      --data-dir ./data_g0_human --out ./data_g0_human/ica
"""
from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np

import v0b_module_definitions as v0b
import baselines_nmf_pca as base


def main():
    from sklearn.decomposition import FastICA

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default="./data_g0")
    ap.add_argument("--rank", type=int, default=50,
                    help="ICA converges at ~50 components; 512 (SAE-matched) does not. See docstring.")
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--out", default="./data_g0/ica")
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--max-iter", type=int, default=1000)
    args = ap.parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(v0b.SEED)

    data_dir = Path(args.data_dir)
    gene_names = v0b.load_gene_names(data_dir)
    n_genes = len(gene_names)
    X = v0b.load_expression(data_dir)
    if X is None:
        raise SystemExit("expression_matrix.npy required for the ICA baseline.")
    coverage, usable = v0b.compute_coverage(gene_names)       # respects V0B_MARKERS_JSON
    ery = [m for m in v0b.ERY_SUBMODULES if m in usable]
    gran = [m for m in v0b.GRAN_SUBMODULES if m in usable]

    # FastICA at a few random_states -> per-seed gene-loading decoders (n_genes x rank).
    components, converged = {}, {}
    for s in range(args.seeds):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")                  # non-convergence is tolerable
            ica = FastICA(n_components=args.rank, random_state=s, max_iter=args.max_iter,
                          whiten="unit-variance")
            ica.fit(X)
        components[s] = np.abs(ica.components_.T).astype(float)   # (n_genes, rank)
        converged[s] = int(ica.n_iter_ < args.max_iter)

    modularity = base._modularity_on_components(
        components, coverage, usable, n_genes, ery, gran, rng, n_perm=args.n_perm)
    pr = {"data": base.pr_of_representation(X)}

    result = {"method": "FastICA (independent components)",
              "rank": args.rank, "seeds": args.seeds, "n_genes": n_genes,
              "converged_per_seed": converged,
              "participation_ratio": pr, "modularity": modularity}
    with open(out / "ica_baseline.json", "w") as f:
        json.dump(result, f, indent=2)

    d = modularity["decision"]
    print(f"ICA baseline (rank {args.rank}, {args.seeds} seeds, "
          f"converged {sum(converged.values())}/{args.seeds}):")
    print(f"  modularity: ery_n_real={d['median_ery_n_real']} "
          f"gran_n_real={d['median_gran_n_real']} supported={d['asymmetric_modularity_supported']}")
    print(f"  finding: {d['plain_finding']}")


if __name__ == "__main__":
    main()
