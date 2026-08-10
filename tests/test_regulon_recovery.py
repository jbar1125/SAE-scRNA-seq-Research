#!/usr/bin/env python3
"""Synthetic-oracle test for the regulon-recovery causal metric.

Plants a world where the ground truth is known, then checks the metric recovers exactly the
causal TFs and that the label-shuffle null collapses. Pure numpy; no torch/data needed.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import causal_grounding as cg  # noqa: E402
import regulon_recovery as rr  # noqa: E402


def _build(seed=0):
    """Construct activations/expression/labels + regulons with a known answer.

    - Genes 0..59. Feature f_L loads genes [10*L .. 10*L+9] (block L). 6 feature blocks.
    - TF "CAUSAL_L" has TRRUST regulon = block-L genes AND its OE cells strongly activate
      feature L (so its regulon-feature is up in OE) -> MUST recover.
    - TF "WRONGFEAT" has a real regulon (block 5) but its OE cells do NOT activate feature 5
      (activate nothing) -> must NOT recover (regulon feature not causally activated).
    - TF "NOREGULON" activates feature 0 but has no TRRUST regulon -> not tested.
    """
    rng = np.random.default_rng(seed)
    G, F = 60, 6
    gene_names = [f"g{i}" for i in range(G)]
    # feature -> block of genes it "programs"
    blocks = {f: list(range(10 * f, 10 * f + 10)) for f in range(F)}

    causal_tfs = [f"CAUSAL_{L}" for L in range(4)]         # 4 clean causal TFs (features 0..3)
    per_tf = 80
    labels = []
    rows_act = []   # activation rows (n_cells, F)
    rows_expr = []  # expression rows (n_cells, G)

    def emit(label, active_feats):
        labels.append(label)
        a = np.abs(rng.normal(0, 0.05, F))
        for f, amp in active_feats.items():
            a[f] = amp + rng.normal(0, 0.1)
        rows_act.append(a)
        # expression: control-cell co-expression must make gene_association put block genes on
        # their feature, so drive expression from the FULL activation vector every cell.
        x = rng.normal(0, 0.1, G)
        for f in range(F):
            for gi in blocks[f]:
                x[gi] += a[f] * (2.0 if gi != 10 * f else 2.0)  # feature f raises its block
        rows_expr.append(x)

    # controls: random low activation across all features (defines gene_association structure)
    for _ in range(400):
        emit("control", {f: abs(rng.normal(0.4, 0.15)) for f in range(F)})
    # causal TFs: OE strongly activates their own feature
    for L in range(4):
        for _ in range(per_tf):
            emit(f"CAUSAL_{L}", {L: 2.5})
    # wrong-feature TF: regulon is block 5, but OE activates nothing special
    for _ in range(per_tf):
        emit("WRONGFEAT", {})
    # no-regulon TF: activates feature 0 but absent from TRRUST
    for _ in range(per_tf):
        emit("NOREGULON", {0: 2.5})

    acts = np.array(rows_act)
    expr = np.array(rows_expr)
    labels = np.array(labels)
    regulons = {f"CAUSAL_{L}": {f"g{i}" for i in blocks[L]} for L in range(4)}
    regulons["WRONGFEAT"] = {f"g{i}" for i in blocks[5]}
    # NOREGULON deliberately absent
    return acts, expr, labels, gene_names, regulons, causal_tfs


def test_oracle_recovers_causal_only():
    acts, expr, labels, gnames, regulons, causal = _build()
    out = rr.regulon_recovery(acts, expr, labels, gnames, regulons,
                              top_k=10, min_targets=5, n_shuffle=0)
    rec = set(out["recovered_tfs"])
    assert rec == set(causal), f"expected exactly {causal}, got {sorted(rec)}"
    # WRONGFEAT was tested (has a regulon) but must not recover
    wf = next(r for r in out["per_tf"] if r["tf"] == "WRONGFEAT")
    assert not wf["recovers"], "WRONGFEAT activates no feature -> must not recover"
    # NOREGULON has no regulon -> not in the tested pool at all
    assert all(r["tf"] != "NOREGULON" for r in out["per_tf"]), "NOREGULON must be untested"
    print(f"oracle: recovered exactly the causal set {sorted(rec)}; WRONGFEAT/NOREGULON rejected")


def test_shuffle_null_collapses():
    acts, expr, labels, gnames, regulons, causal = _build()
    out = rr.regulon_recovery(acts, expr, labels, gnames, regulons,
                              top_k=10, min_targets=5, n_shuffle=40, seed=1)
    sn = out["shuffle_null"]
    assert out["n_recovered"] == len(causal), out["n_recovered"]
    assert sn["null_mean"] < out["n_recovered"], (sn["null_mean"], out["n_recovered"])
    assert sn["empirical_p"] < 0.05, sn["empirical_p"]
    print(f"shuffle null: real {out['n_recovered']} vs null mean {sn['null_mean']:.2f} "
          f"(max {sn['null_max']}), p={sn['empirical_p']:.3f}")


def test_low_power_flagged():
    acts, expr, labels, gnames, regulons, causal = _build()
    regulons["TINY"] = {"g0", "g1"}          # 2 targets -> below min_targets
    # give TINY some cells so it is present as a label
    extra_lab = np.array(["TINY"] * 40)
    acts2 = np.vstack([acts, acts[:40]])
    expr2 = np.vstack([expr, expr[:40]])
    labels2 = np.concatenate([labels, extra_lab])
    out = rr.regulon_recovery(acts2, expr2, labels2, gnames, regulons,
                              top_k=10, min_targets=5, n_shuffle=0)
    tiny = next(r for r in out["per_tf"] if r["tf"] == "TINY")
    assert tiny["low_power"] and not tiny["recovers"], "TINY must be low-power, not recovered"
    print("low-power: TINY (2 targets) flagged low_power, excluded from tested pool")


if __name__ == "__main__":
    test_oracle_recovers_causal_only()
    test_shuffle_null_collapses()
    test_low_power_flagged()
    print("ALL REGULON-RECOVERY TESTS PASSED")
