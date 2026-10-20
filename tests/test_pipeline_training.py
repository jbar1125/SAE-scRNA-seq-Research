#!/usr/bin/env python3
"""Tests for WHICH cells the causal pipelines fit the SAE on (train_rows / --train-on).

Finding (2026-10-18): config/causal_grounding_spec.json states the SAE is trained on control cells
only, but every caller passed the full matrix to train_topk_sae, so every real Component-2 run fit
the SAE on all cells (perturbation labels unseen). The synthetic self-test only validates that
all-cell mode: it grounds 8/8 planted regulators, while control-only grounds 0/8 (signal present,
AUC ~0.40, but below FDR). So 'all' stays the default and 'control' is exposed as UNVALIDATED.

These tests pin that behaviour: train_rows() masks; run() fits on exactly the selected rows and
still encodes every cell; the default reproduces the validated all-cell pipeline on a real TopK SAE;
the control-only result is printed for the record but not asserted (it is a known failure).
The row checks stub the trainer (no torch). Run: python3 tests/test_pipeline_training.py
"""
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import causal_pipeline as cp  # noqa: E402


def test_train_rows():
    labels = np.array(["control", "A", "control", "B", "control"])
    assert cp.train_rows(labels).all(), "default must be 'all' (the validated mode)"
    assert cp.train_rows(labels, "control").tolist() == [True, False, True, False, True]
    for bad in ("perturbed", None):
        try:
            cp.train_rows(labels, bad)
            raise AssertionError(f"train_on={bad!r} must be rejected")
        except ValueError:
            pass
    try:
        cp.train_rows(np.array(["A", "B"]), "control")
        raise AssertionError("control mode with no control cells must raise, not fit on nothing")
    except ValueError:
        pass
    print("train_rows: default 'all'; 'control' masks control rows; bad mode and zero-control "
          "inputs raise")


def test_run_fits_on_selected_rows_and_encodes_all():
    X, genes, labels, tested = cp._synthetic(n_prog=3, prog_size=6, n_bg=30)
    seen = {}
    real = cp.train_topk_sae

    def stub(rep, latent, k, seed=0, **kw):
        seen["n_fit"] = rep.shape[0]

        def encode(mat):
            seen["n_encoded"] = mat.shape[0]
            rng = np.random.default_rng(seed)
            return np.abs(rng.normal(0, 1, (mat.shape[0], latent))).astype(np.float32)
        return encode, {"latent_dim": latent, "k": k, "d": rep.shape[1], "device": "stub"}

    cp.train_topk_sae = stub
    try:
        with tempfile.TemporaryDirectory() as tmp:
            n_ctrl = int((labels == "control").sum())
            res = cp.run(X, X, genes, labels, tested, latent=16, k=4, seed=0,
                         out=os.path.join(tmp, "a.json"))
            assert seen["n_fit"] == X.shape[0] and seen["n_encoded"] == X.shape[0], seen
            assert res["sae"]["trained_on"] == "all" and res["sae"]["n_train_cells"] == X.shape[0]
            res = cp.run(X, X, genes, labels, tested, latent=16, k=4, seed=0,
                         out=os.path.join(tmp, "c.json"), train_on="control")
            assert seen["n_fit"] == n_ctrl and seen["n_encoded"] == X.shape[0], (seen, n_ctrl)
            assert res["sae"]["trained_on"] == "control" and res["sae"]["n_train_cells"] == n_ctrl
    finally:
        cp.train_topk_sae = real
    print(f"run(): default fits on all {X.shape[0]} cells; train_on='control' fits on the "
          f"{n_ctrl} control rows; every cell is encoded either way; trained_on is recorded")


def _self_test(train_on):
    X, genes, labels, tested = cp._synthetic()
    n_prog = sum(1 for t in tested if not t.startswith("bg"))
    with tempfile.TemporaryDirectory() as tmp:
        res = cp.run(X, X, genes, labels, tested, latent=128, k=8, seed=0,
                     out=os.path.join(tmp, f"{train_on}.json"), n_shuffle=20, train_on=train_on)
    pp = res["per_perturbation"]
    true_g = sum(pp[i]["grounded"] for i in range(n_prog))
    return res, n_prog, true_g, res["n_grounded"] - true_g, len(tested)


def test_real_sae_default_is_validated_and_control_is_recorded():
    try:
        import torch  # noqa: F401
    except ImportError:
        print("SKIP (torch not installed): real-SAE self-tests NOT run")
        return
    res, n_prog, true_g, false_g, n_tested = _self_test("all")
    assert res["sae"]["trained_on"] == "all"
    assert true_g >= n_prog - 2, (true_g, n_prog)
    assert false_g <= max(3, int(0.05 * (n_tested - n_prog))), false_g
    assert res["shuffle_control"]["null_rates_max"] < res["causal_grounding_rate"]
    print(f"real TopK SAE, default all-cell fit: {true_g}/{n_prog} true regulators, "
          f"{false_g} false positives (validated)")
    res, n_prog, true_g, false_g, _ = _self_test("control")
    assert res["sae"]["trained_on"] == "control"
    print(f"real TopK SAE, control-only fit (pre-registered intent, UNVALIDATED): "
          f"{true_g}/{n_prog} true regulators, {false_g} false positives "
          f"[known limitation, recorded, not asserted]")


if __name__ == "__main__":
    test_train_rows()
    test_run_fits_on_selected_rows_and_encodes_all()
    test_real_sae_default_is_validated_and_control_is_recorded()
    print("ALL PIPELINE-TRAINING TESTS PASSED")
