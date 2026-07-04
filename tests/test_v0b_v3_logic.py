"""Synthetic smoke tests for the V0b v3 continuous-loading logic.

No real data, no torch. Engineers decoder weights with known concentration so the
loading-enrichment -> strength -> distribution -> decision path is verified.
Run: python tests/test_v0b_v3_logic.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import v0b_module_definitions as v0b
import v0b_v3_loading as v3


# Core v2 submodules only (exclude marker-aware-only Effector/Secondary).
_CORE = [m for m in v0b.MARKER_SETS if m not in ("Ery_Effector", "Gran_Secondary")]


def _build_gene_index():
    markers = []
    for m in _CORE:
        for g in v0b.MARKER_SETS[m]:
            if g not in markers:
                markers.append(g)
    filler = [f"FILLER_{i}" for i in range(300 - len(markers))]
    names = markers + filler
    return names, {g: i for i, g in enumerate(names)}


def _weights(plan, gene_to_idx, n_genes):
    """plan: {feature_idx: submodule}. Targeted features put large decoder mass on
    that submodule's markers; other genes get a tiny positive ramp so untargeted
    submodules sit near the null (enrichment ~1)."""
    n_filler = sum(1 for g in gene_to_idx if g.startswith("FILLER_"))
    W = np.zeros((n_genes, v0b.LATENT_DIM))
    for g, i in gene_to_idx.items():
        if g.startswith("FILLER_"):
            W[i, :] = 0.001 * (n_filler - int(g.split("_")[1]))
        else:
            W[i, :] = 0.001  # markers get small baseline in all features
    for feat, module in plan.items():
        for g in v0b.MARKER_SETS[module]:
            W[gene_to_idx[g], feat] = 10.0
    return {s: W.copy() for s in range(v0b.N_SEEDS)}


def _run(plan, n_perm=200):
    gene_names, gene_to_idx = _build_gene_index()
    n_genes = len(gene_names)
    dw = _weights(plan, gene_to_idx, n_genes)
    coverage, usable = v0b.compute_coverage(gene_names)
    enrich = v3.loading_enrichment(dw, coverage, usable, n_genes)
    strength = v3.submodule_strength(enrich, usable)
    ery = [m for m in v0b.ERY_SUBMODULES if m in usable]
    gran = [m for m in v0b.GRAN_SUBMODULES if m in usable]
    mod = v3.per_seed_distribution(strength, ery, gran)
    decision = v0b.decide(mod, np.random.default_rng(0))
    perm = v3.permutation_null(dw, coverage, usable, n_genes, n_perm, np.random.default_rng(0))
    return strength, mod, decision, perm


def _plan(per_module, start=0):
    plan, f = {}, start
    for module, n in per_module.items():
        for _ in range(n):
            plan[f] = module
            f += 1
    return plan


def test_strength_high_for_targeted_low_for_untargeted():
    strength, _, _, perm = _run(_plan({"Ery_TF": 3, "Gran_Primary": 3}))
    mean_s = strength.groupby("module")["strength"].mean()
    assert mean_s["Ery_TF"] > 5, mean_s.to_string()
    assert mean_s["Gran_Primary"] > 5, mean_s.to_string()
    # an untargeted submodule sits near the null (enrichment ~1)
    assert mean_s["Ery_Membrane"] < 2, mean_s.to_string()
    # targeted submodules are significant under the permutation null
    assert perm["Ery_TF"]["perm_p"] < 0.05, perm["Ery_TF"]
    print("PASS test_strength_high_for_targeted_low_for_untargeted")


def test_supported_when_ery_distributed_gran_unified():
    plan = _plan({"Ery_TF": 3, "Ery_Heme": 3, "Ery_Membrane": 3, "Gran_TF": 3})
    _, mod, decision, _ = _run(plan)
    assert decision["n_seeds_both_directions_hold"] == 5, mod.to_string()
    assert decision["asymmetric_modularity_supported"] is True, decision
    print("PASS test_supported_when_ery_distributed_gran_unified")


def test_not_supported_when_ery_unified():
    plan = _plan({"Ery_TF": 3, "Gran_TF": 3, "Gran_Primary": 3})
    _, _, decision, _ = _run(plan)
    assert decision["asymmetric_modularity_supported"] is False, decision
    print("PASS test_not_supported_when_ery_unified")


def test_metrics_bounds():
    _, mod, _, _ = _run(_plan({"Ery_TF": 3, "Ery_Heme": 3, "Gran_TF": 3}))
    for col in ["ery_norm_entropy", "gran_norm_entropy"]:
        assert (mod[col] >= -1e-9).all() and (mod[col] <= 1 + 1e-9).all(), col
    print("PASS test_metrics_bounds")


if __name__ == "__main__":
    test_strength_high_for_targeted_low_for_untargeted()
    test_supported_when_ery_distributed_gran_unified()
    test_not_supported_when_ery_unified()
    test_metrics_bounds()
    print("\nALL V0b v3 LOGIC TESTS PASSED")
