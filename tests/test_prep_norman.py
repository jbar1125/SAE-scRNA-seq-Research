#!/usr/bin/env python3
"""Tests for src/prep_norman.py: guide parsing, symbol remap, and the raw-count guard.

The raw-count guard exists because pertpy ships two Norman files: norman_2019_raw (integer
counts; what every reported Norman result used) and norman_2019 (log-normalized X, raw counts in
layers['counts']). load_perturbseq normalizes on its own, so normalized X must never be written.
Needs numpy + anndata; no torch, no network. Run: python3 tests/test_prep_norman.py
"""
import os
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import prep_norman as pn  # noqa: E402

SCRIPT = os.path.join(os.path.dirname(__file__), "..", "src", "prep_norman.py")
GENES = ["KLF1", "CEBPA", "SPI1", "ETS2", "GATA1", "RUNX1", "RUNX1T1", "CEBPE", "HBB", "HBA1"]


def _raw_adata(seed=0, n=900):
    """Synthetic Norman-like AnnData: Ensembl var_names + a symbol column, integer counts."""
    import anndata as ad
    import pandas as pd
    rng = np.random.default_rng(seed)
    guides = []
    for _ in range(n):
        r = rng.random()
        if r < 0.4:
            guides.append(f"NegCtrl{rng.integers(0, 3)}_NegCtrl{rng.integers(0, 3)}__b")
        elif r < 0.85:
            guides.append(f"{rng.choice(GENES[:6])}_NegCtrl{rng.integers(0, 3)}__b")
        else:
            guides.append("CEBPE_RUNX1T1_")
    X = rng.poisson(1.5, size=(n, len(GENES))).astype("float32")
    var = pd.DataFrame({"gene_symbols": GENES}, index=[f"ENSG{100000 + i:06d}" for i in range(len(GENES))])
    obs = pd.DataFrame({"guide_identity": guides}, index=[f"c{i}" for i in range(n)])
    return ad.AnnData(X=X, obs=obs, var=var)


def _normalized(A):
    """Log-normalized copy of X, like pertpy's processed norman_2019."""
    X = A.X / A.X.sum(1, keepdims=True).clip(min=1) * 1e4
    return np.log1p(X).astype("float32")


def test_parse_guide():
    cases = {"NegCtrl0_NegCtrl1__abc": "control", "KLF1_NegCtrl0__xyz": "KLF1",
             "CEBPE_RUNX1T1_": "CEBPE+RUNX1T1", "RUNX1T1_CEBPE": "CEBPE+RUNX1T1",
             "SPI1": "SPI1", "NegCtrl9": "control"}
    for k, v in cases.items():
        assert pn.parse_guide(k) == v, (k, pn.parse_guide(k), v)
    print("parse_guide: control / single / dual (order-invariant) all correct")


def test_ensure_raw_counts_cases():
    A = _raw_adata()
    _, note = pn.ensure_raw_counts(A)
    assert "integer" in note, note

    B = _raw_adata()
    raw = B.X.copy()
    B.layers["counts"] = raw
    B.X = _normalized(B)
    B, note = pn.ensure_raw_counts(B)
    assert "replaced" in note and np.array_equal(B.X, raw), note

    C = _raw_adata()
    C.X = _normalized(C)
    try:
        pn.ensure_raw_counts(C)
        raise AssertionError("normalized X without a counts layer must be refused")
    except SystemExit as e:
        assert "norman_2019_raw" in str(e), e
    print("ensure_raw_counts: raw passes; processed is repaired from layers['counts']; "
          "normalized-without-counts is refused")


def _run_cli(A, tmp, tag):
    src, out = os.path.join(tmp, f"{tag}_in.h5ad"), os.path.join(tmp, f"{tag}_out.h5ad")
    A.write_h5ad(src)
    r = subprocess.run([sys.executable, SCRIPT, "--in", src, "--out", out],
                       capture_output=True, text=True)
    return r, out


def test_cli_end_to_end():
    import anndata as ad
    with tempfile.TemporaryDirectory() as tmp:
        A = _raw_adata()
        r, out = _run_cli(A, tmp, "raw")
        assert r.returncode == 0, r.stderr
        B = ad.read_h5ad(out)
        assert "KLF1" in list(B.var_names) and "ensembl_id" in B.var.columns
        assert "control" in set(B.obs["perturbation"])
        assert np.array_equal(B.X, A.X), "raw counts must be written unchanged"

        P = _raw_adata()
        raw = P.X.copy()
        P.layers["counts"] = raw
        P.X = _normalized(P)
        r, out = _run_cli(P, tmp, "processed")
        assert r.returncode == 0, r.stderr
        assert "replaced it with raw layers['counts']" in r.stdout, r.stdout
        assert np.array_equal(ad.read_h5ad(out).X, raw), "processed input must be written as raw"

        N = _raw_adata()
        N.X = _normalized(N)
        r, out = _run_cli(N, tmp, "normalized")
        assert r.returncode != 0 and not os.path.exists(out), "must refuse and write nothing"
    print("CLI end to end: raw written unchanged; processed written as raw; normalized refused "
          "with no output file")


if __name__ == "__main__":
    try:
        import anndata  # noqa: F401
    except ImportError:
        sys.exit("SKIP: anndata not installed (pip install anndata); prep_norman tests NOT run")
    test_parse_guide()
    test_ensure_raw_counts_cases()
    test_cli_end_to_end()
    print("ALL PREP-NORMAN TESTS PASSED")
