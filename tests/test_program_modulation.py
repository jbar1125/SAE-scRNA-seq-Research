#!/usr/bin/env python3
"""Synthetic-oracle test for the program-modulation metric.

Two worlds must be told apart, because confusing them is exactly how this project's earlier
metrics failed:
  A. SPECIFIC world -- one program is activated by a small, specific perturbation set. The
     metric must find those modulators, beat the shuffle null, and rank HIGH against the
     random-feature background.
  B. DOMINANT-AXIS world -- one feature moves for nearly EVERY perturbation (K562's real
     failure mode). The metric may report many modulators, but the feature-background
     percentile must NOT single it out, so the run reads as a negative rather than a discovery.
Pure numpy; no torch or data needed.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import program_modulation as pm  # noqa: E402


def _world(specific=True, seed=0):
    """Build activations/expression/labels with a known answer.

    Feature 0 carries the ERY program (genes 0..9, the marker set). In the SPECIFIC world only
    perturbations P0..P3 activate it; in the DOMINANT world nearly all perturbations do, and so
    do many OTHER features (so the background is high too)."""
    rng = np.random.default_rng(seed)
    G, F = 80, 40
    gene_names = [f"g{i}" for i in range(G)]
    blocks = {f: list(range(2 * f, 2 * f + 2)) for f in range(F)}
    blocks[0] = list(range(10))                     # feature 0 = the ERY program genes
    n_perts = 20
    perts = [f"P{i}" for i in range(n_perts)]
    movers = set(perts[:4]) if specific else set(perts[:18])

    labels, rows_a, rows_x = [], [], []

    def emit(lab, boost):
        labels.append(lab)
        a = np.abs(rng.normal(0.5, 0.2, F))         # every feature fires in control cells
        for f, amp in boost.items():
            a[f] += amp
        rows_a.append(a)
        x = rng.normal(0, 0.1, G)
        for f in range(F):
            for gi in blocks[f]:
                x[gi] += a[f] * 2.0
        rows_x.append(x)

    for _ in range(500):
        emit("control", {})
    for p in perts:
        for _ in range(60):
            if p in movers:
                if specific:
                    emit(p, {0: 2.5})               # only feature 0 moves
                else:
                    # dominant axis: feature 0 AND many others move for almost everything
                    emit(p, {0: 2.5, **{f: 2.0 for f in rng.choice(range(1, F), 12,
                                                                   replace=False)}})
            else:
                emit(p, {})

    return (np.array(rows_a), np.array(rows_x), np.array(labels), gene_names,
            {"ERY": [f"g{i}" for i in range(10)]}, sorted(movers))


def test_specific_world_detected():
    acts, expr, labels, gn, mk, movers = _world(specific=True)
    out = pm.program_modulation(acts, expr, labels, gn, mk, top_k=10, min_active_cells=10,
                                n_shuffle=30, n_random_features=30, seed=1)
    r = out["programs"]["ERY"]
    assert not r["low_power"]
    found = set(map(str, r["modulators"]))
    # every true mover must be found; a few BH false positives are expected and allowed
    # (measured FP rate over 8 seeds = 0.031 against a nominal alpha of 0.05)
    assert set(movers) <= found, f"missed true movers: {set(movers) - found}"
    assert len(found - set(movers)) <= 2, f"too many false positives: {found - set(movers)}"
    assert r["shuffle_null"]["empirical_p"] < 0.05, r["shuffle_null"]
    assert r["feature_background_percentile"] >= 90, r["feature_background_percentile"]
    print(f"specific world: found all movers {sorted(found)}, shuffle p="
          f"{r['shuffle_null']['empirical_p']:.3f}, background pct "
          f"{r['feature_background_percentile']:.0f} -> reads as POSITIVE")


def test_dominant_axis_not_flagged_as_discovery():
    acts, expr, labels, gn, mk, movers = _world(specific=False)
    out = pm.program_modulation(acts, expr, labels, gn, mk, top_k=10, min_active_cells=10,
                                n_shuffle=0, n_random_features=30, seed=1)
    r = out["programs"]["ERY"]
    bg = out["feature_background"]
    # many modulators, but the background is comparably high -> percentile must not single it out
    assert r["n_modulators"] >= 10, r["n_modulators"]
    assert bg["mean"] >= 5, bg
    assert r["feature_background_percentile"] < 90, (
        f"dominant axis must NOT rank as a discovery, got pct "
        f"{r['feature_background_percentile']:.0f} with bg mean {bg['mean']:.1f}")
    print(f"dominant-axis world: {r['n_modulators']} modulators BUT background mean "
          f"{bg['mean']:.1f} -> percentile {r['feature_background_percentile']:.0f} "
          f"-> correctly reads as NEGATIVE, not a discovery")


def test_low_power_program_flagged():
    acts, expr, labels, gn, mk, _ = _world(specific=True)
    mk = dict(mk); mk["TINY"] = ["g0", "NOT_IN_PANEL"]      # only 1 marker present
    out = pm.program_modulation(acts, expr, labels, gn, mk, top_k=10, min_active_cells=10,
                                n_shuffle=0, n_random_features=10, seed=1)
    assert out["programs"]["TINY"]["low_power"] is True
    print("low-power: program with <3 in-panel markers flagged, not scored")


def test_label_free_selection():
    """Feature selection must not depend on perturbation labels (shuffle-invariance)."""
    acts, expr, labels, gn, mk, _ = _world(specific=True)
    out1 = pm.program_modulation(acts, expr, labels, gn, mk, top_k=10, min_active_cells=10,
                                 n_shuffle=0, n_random_features=10, seed=1)
    rng = np.random.default_rng(7)
    shuffled = labels.copy()
    nz = np.where(labels != "control")[0]
    shuffled[nz] = rng.permutation(labels[nz])
    out2 = pm.program_modulation(acts, expr, shuffled, gn, mk, top_k=10, min_active_cells=10,
                                 n_shuffle=0, n_random_features=10, seed=1)
    assert out1["programs"]["ERY"]["feature"] == out2["programs"]["ERY"]["feature"]
    print("label-free: selected feature identical under label shuffle (non-circular)")


if __name__ == "__main__":
    test_specific_world_detected()
    test_dominant_axis_not_flagged_as_discovery()
    test_low_power_program_flagged()
    test_label_free_selection()
    print("ALL PROGRAM-MODULATION TESTS PASSED")
