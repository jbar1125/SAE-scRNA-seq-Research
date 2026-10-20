#!/usr/bin/env python3
"""Prep the Norman 2019 CRISPRa Perturb-seq h5ad for causal_pipeline.py (`--direction up`).

Why this exists: the raw Norman deposit encodes the perturbation in a dual-guide
`guide_identity` string (e.g. "KLF1_NegCtrl0__...", "CEBPE_RUNX1T1_...") and its var_names
are ENSEMBL IDs, while perturbations are gene SYMBOLS. `load_perturbseq` tests a perturbation
only if its label is a var_name (a single gene present in the matrix), so without this prep 0
perturbations are testable (the ZeroDivisionError we hit). This script:

  1. Parses `guide_identity` -> a clean perturbation label: "control" for all-NegCtrl guides,
     the single gene for single-gene guides, "A+B" (sorted) for dual-gene guides. Dual labels
     are intentionally NOT single genes, so load_perturbseq drops them from the testable set --
     only single-gene overexpressions are grounded, which is what we want.
  2. Detects the var column of gene SYMBOLS by overlap with the perturbation gene tokens, maps
     var_names to symbols, and stashes the original Ensembl IDs in var["ensembl_id"].

It writes RAW counts (load_perturbseq applies normalize_total+log1p itself -- do NOT pre-normalize
here). Then run:
  python3 src/causal_pipeline.py --adata norman_pert.h5ad \
      --pert-col perturbation --control-value control --direction up ...

INPUT FILE (verified 2026-10-18 against both pertpy files): download with
`pertpy.data.norman_2019_raw()` (111,445 cells x 33,694 genes, integer counts, Ensembl
var_names + `gene_symbols`, no `perturbation_name`): this is the file every reported Norman
result used. `pertpy.data.norman_2019()` is the PROCESSED file (111,255 cells, log-normalized X,
raw counts only in layers['counts']); feeding its X downstream would double-normalize.
`ensure_raw_counts` below therefore checks X and swaps in layers['counts'] or refuses.

anndata + numpy only. Usage:
  python3 src/prep_norman.py --in norman_raw.h5ad --out norman_pert.h5ad [--guide-col guide_identity]
"""
import argparse

import numpy as np


def parse_guide(gi):
    """guide_identity string -> perturbation label ('control' | 'GENE' | 'A+B')."""
    head = str(gi).split("__")[0]
    genes = [t for t in head.split("_") if t and not t.startswith("NegCtrl")]
    if not genes:
        return "control"
    if len(genes) == 1:
        return genes[0]
    return "+".join(sorted(set(genes)))


def perturbation_tokens(labels):
    """The set of single gene symbols named by any (single or dual) perturbation label."""
    toks = set()
    for lab in labels:
        if lab == "control":
            continue
        toks.update(lab.split("+"))
    return toks


def detect_symbol_column(A, tokens):
    """Return the var column whose values best overlap the perturbation gene symbols, or None
    if the var_names themselves already overlap best (already symbols)."""
    best_col, best_hits = None, 0
    # baseline: do the current var_names already look like symbols?
    vn = set(map(str, A.var_names))
    base_hits = len(tokens & vn)
    for c in A.var.columns:
        vals = set(A.var[c].astype(str))
        hits = len(tokens & vals)
        if hits > best_hits:
            best_col, best_hits = c, hits
    # only remap if a column beats the current index materially
    if best_col is not None and best_hits > base_hits:
        return best_col, best_hits, base_hits
    return None, base_hits, base_hits


def _nonzero_sample(M, n=200000):
    """Up to n nonzero values of a sparse or dense matrix (for a cheap integer check)."""
    d = M.data if hasattr(M, "data") else np.asarray(M).ravel()
    d = np.asarray(d[:n], dtype=float)
    return d[d != 0]


def ensure_raw_counts(A):
    """Return (A, note) with A.X holding raw integer counts, or raise SystemExit.

    load_perturbseq applies normalize_total + log1p itself, so X written here must be raw. If X
    is not integer-valued but a 'counts' layer is (pertpy's processed norman_2019), that layer is
    swapped in; otherwise refuse, because writing normalized X would double-normalize every
    downstream result (the same error class as the Paul15 log1p mistake, PROJECT_AUDIT A11)."""
    d = _nonzero_sample(A.X)
    if d.size == 0 or np.allclose(d, np.round(d)):
        return A, "X is integer-valued (raw counts)"
    if "counts" in A.layers:
        dc = _nonzero_sample(A.layers["counts"])
        if dc.size and np.allclose(dc, np.round(dc)):
            A.X = A.layers["counts"].copy()
            return A, "X was normalized (max %.3g); replaced it with raw layers['counts']" % d.max()
    raise SystemExit("X is not raw integer counts (max %.3g, non-integer) and there is no raw "
                     "'counts' layer. Use pertpy.data.norman_2019_raw()." % d.max())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="inp", required=True, help="raw Norman h5ad (from pertpy)")
    ap.add_argument("--out", default="norman_pert.h5ad")
    ap.add_argument("--guide-col", default="guide_identity",
                    help="obs column holding the dual-guide identity string")
    ap.add_argument("--pert-col", default="perturbation",
                    help="name of the parsed perturbation column to create")
    args = ap.parse_args()

    import anndata as ad
    A = ad.read_h5ad(args.inp)
    print(f"loaded {A.shape[0]} cells x {A.shape[1]} genes")
    A, note = ensure_raw_counts(A)
    print(f"raw-count check: {note}")
    if args.guide_col not in A.obs.columns:
        raise SystemExit(f"--guide-col '{args.guide_col}' not in obs; columns: {list(A.obs.columns)}")

    labels = [parse_guide(g) for g in A.obs[args.guide_col].astype(str)]
    A.obs[args.pert_col] = labels
    from collections import Counter
    vc = Counter(labels)
    n_ctrl = vc.get("control", 0)
    n_single = sum(v for k, v in vc.items() if k != "control" and "+" not in k)
    n_dual = sum(v for k, v in vc.items() if "+" in k)
    print(f"parsed perturbations: control={n_ctrl} single-gene cells={n_single} "
          f"dual-gene cells={n_dual}; distinct single genes="
          f"{len({k for k in vc if k != 'control' and '+' not in k})}")

    tokens = perturbation_tokens(labels)
    col, hits, base = detect_symbol_column(A, tokens)
    if col is None:
        print(f"var_names already overlap perturbation symbols ({base}/{len(tokens)}); no remap")
    else:
        print(f"remapping var_names to symbols via var['{col}'] "
              f"(overlap {hits}/{len(tokens)} vs {base} for current index)")
        A.var["ensembl_id"] = list(map(str, A.var_names))
        syms = A.var[col].astype(str).tolist()
        del A.var[col]                      # avoid index-name == column-name write clash
        A.var_names = syms

    # sanity: how many single-gene perturbations are now testable (label is a var_name)?
    present = set(map(str, A.var_names))
    testable = sorted({k for k in vc if k != "control" and "+" not in k and k in present
                       and vc[k] >= 30})
    print(f"testable single-gene perturbations (>=30 cells, symbol in var): {len(testable)}")
    print(f"  examples: {testable[:15]}")
    for hema in ["KLF1", "CEBPA", "SPI1", "ETS2", "GATA1", "RUNX1"]:
        print(f"  {hema}: {'TESTABLE' if hema in testable else 'not testable'} "
              f"({vc.get(hema, 0)} cells)")

    A.write_h5ad(args.out)
    print(f"wrote {args.out} (raw counts; var_names=symbols). "
          f"Run causal_pipeline.py --pert-col {args.pert_col} --control-value control --direction up")


if __name__ == "__main__":
    main()
