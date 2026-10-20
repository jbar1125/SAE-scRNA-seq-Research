#!/usr/bin/env python3
"""Inspect a Perturb-seq / overexpression h5ad to CONFIGURE a grounding run.

Run this on a NEW dataset (e.g. the Joung TF Atlas GSE216481) BEFORE grounding, so
--pert-col / --control-value / --direction are set from the real obs schema rather than
guessed. It prints the shape, every obs column with its cardinality, and the top values of
the low-cardinality (perturbation-like) columns -- which reveals (a) which column names the
perturbed gene/TF and (b) what the control/baseline label is. It also flags candidate
control labels (non-targeting / GFP / mCherry / control / EB / baseline).

See docs/METHODOLOGY_ADDENDUM_OE.md (the deposit's missing-control caveat). numpy/anndata
only; reads backed (memory-light) so it works on the full 254k-cell atlas.

Usage: python3 src/inspect_perturb_h5ad.py --adata atlas.h5ad
"""
import argparse
import re

CONTROL_HINTS = re.compile(r"non.?targeting|control|ctrl|gfp|mcherry|baseline|\beb\b|"
                           r"unperturb|wildtype|wt|dmso|safe.?harbor|intergenic", re.I)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--adata", required=True)
    ap.add_argument("--top", type=int, default=15, help="top values to show per column")
    ap.add_argument("--max-card", type=int, default=6000,
                    help="only detail columns with <= this many unique values (perturbation-like)")
    args = ap.parse_args()

    import anndata as ad
    A = ad.read_h5ad(args.adata, backed="r")
    print(f"shape: {A.shape[0]} cells x {A.shape[1]} genes")
    print(f"var names (first 8): {list(map(str, A.var_names[:8]))}")
    print(f"obs columns ({len(A.obs.columns)}): {list(A.obs.columns)}")

    candidates = []
    for c in A.obs.columns:
        col = A.obs[c]
        if not (col.dtype == object or str(col.dtype).startswith("categ")):
            continue
        s = col.astype(str)
        nun = s.nunique()
        if nun > args.max_card:
            print(f"\n[{c}] categorical, nunique={nun} (too many to be the control axis; skipped)")
            continue
        vc = s.value_counts()
        hits = [v for v in vc.index if CONTROL_HINTS.search(str(v))]
        print(f"\n[{c}] nunique={nun}"
              + (f"  <-- candidate perturbation column" if 2 <= nun <= args.max_card else ""))
        print(vc.head(args.top).to_string())
        if hits:
            print(f"  CANDIDATE CONTROL LABEL(S) in this column: {hits[:8]}")
            candidates.append((c, hits, nun))

    print("\n" + "=" * 60)
    if candidates:
        print("LIKELY CONFIG:")
        for c, hits, nun in candidates:
            print(f"  --pert-col {c}  --control-value '{hits[0]}'   (col has {nun} values)")
        print("  --direction up   (overexpression atlas)  |  down for CRISPRi knockdown")
    else:
        print("No obvious control label found by keyword. Inspect the columns above by hand;")
        print("the control is whatever label marks unperturbed / baseline / GFP cells.")
    print("Then run src/causal_pipeline.py with those --pert-col/--control-value/--direction.")


if __name__ == "__main__":
    main()
