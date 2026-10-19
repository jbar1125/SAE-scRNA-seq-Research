#!/usr/bin/env python3
"""Sum duplicate symbols in prepared raw Norman counts into a separate SAE input.

This declared metadata correction preserves cell order, obs and exact library totals.
The original prepared file remains unchanged. X is the sole expression matrix in the
new file; source Ensembl identifiers are retained as semicolon-delimited provenance.
"""
import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from marker_program_baseline import _file_sha256, sum_duplicate_gene_counts


def prepare(inp, out):
    inp, out = Path(inp), Path(out)
    if out.exists() or out.with_suffix(".manifest.json").exists():
        raise FileExistsError(out)
    source_hash = _file_sha256(inp)
    a = ad.read_h5ad(inp)
    if "perturbation" not in a.obs:
        raise ValueError("input must first be prepared with prep_norman.py")
    original_genes = np.asarray(a.var_names).astype(str)
    matrix, genes = sum_duplicate_gene_counts(a.X, original_genes)
    # Keep the loader's float32 raw-count dtype only when the conversion is exact.
    for start in range(0, len(matrix.data), 1_000_000):
        part = matrix.data[start:start + 1_000_000]
        if not np.array_equal(part, part.astype(np.float32).astype(np.float64)):
            raise ValueError("raw aggregate counts cannot be represented exactly in float32")
    matrix = matrix.astype(np.float32)
    var = pd.DataFrame(index=pd.Index(genes, name="gene_symbol"))
    if "ensembl_id" in a.var:
        identifiers = {}
        for symbol, identifier in zip(original_genes, a.var["ensembl_id"].astype(str)):
            identifiers.setdefault(symbol, []).append(identifier)
        var["source_ensembl_ids"] = [";".join(identifiers[g]) for g in genes]
    b = ad.AnnData(X=matrix, obs=a.obs.copy(), var=var)
    note = {"source": str(inp.resolve()), "source_sha256": source_hash,
            "duplicate_policy": "sum raw counts by symbol before normalization",
            "input_genes": len(original_genes), "unique_genes": len(genes),
            "duplicate_columns_aggregated": len(original_genes) - len(genes),
            "cells": a.n_obs, "per_cell_library_totals_exactly_preserved": True,
            "dtype": str(matrix.dtype)}
    b.uns["symbol_aggregation"] = note
    out.parent.mkdir(parents=True, exist_ok=True)
    b.write_h5ad(out)
    note["output_sha256"] = _file_sha256(out)
    with out.with_suffix(".manifest.json").open("x") as handle:
        json.dump(note, handle, indent=2)
        handle.write("\n")
    return note


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--in", dest="inp", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(prepare(args.inp, args.out)), flush=True)
