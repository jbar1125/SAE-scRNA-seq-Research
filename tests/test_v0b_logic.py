"""Synthetic smoke tests for the V0b decision logic.

These do not use real data or torch. They engineer decoder weights with known
enrichment so the enrichment -> winner-take-all -> modularity -> decision path
can be verified end to end. Run: python tests/test_v0b_logic.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import v0b_module_definitions as v0b


# Core v2 submodules (exclude the marker-aware-only Effector/Secondary so tests
# keep the original 3-ery / 2-gran structure they assert on).
_CORE = [m for m in v0b.MARKER_SETS if m not in ("Ery_Effector", "Gran_Secondary")]


def _build(gene_names):
    """gene_names with every core marker present + filler to 300 genes, plus an index."""
    markers = []
    for m in _CORE:
        for g in v0b.MARKER_SETS[m]:
            if g not in markers:
                markers.append(g)
    filler = [f"FILLER_{i}" for i in range(300 - len(markers))]
    names = markers + filler
    gene_names.extend(names)
    return {g: i for i, g in enumerate(names)}


def _weights_for(assignment_plan, gene_to_idx, n_genes):
    """assignment_plan: dict feature_idx -> submodule (or None). Build a (n_genes,128)
    decoder matrix where each planned feature's top-K are exactly that submodule's
    markers; non-target markers are forced to weight 0 so they never enter top-K."""
    n_filler = sum(1 for g in gene_to_idx if g.startswith("FILLER_"))
    W = np.zeros((n_genes, v0b.LATENT_DIM), dtype=np.float64)
    # Deterministic small positive ramp on filler so filler always outranks the
    # zero-weight non-target markers (keeps top-K free of stray markers).
    for g, i in gene_to_idx.items():
        if g.startswith("FILLER_"):
            W[i, :] = 0.001 * (n_filler - int(g.split("_")[1]))
    for feat, module in assignment_plan.items():
        if module is None:
            continue
        for g in v0b.MARKER_SETS[module]:
            W[gene_to_idx[g], feat] = 10.0  # dominates -> guaranteed in top-K
    return {s: W.copy() for s in range(v0b.N_SEEDS)}


def _run(assignment_plan):
    gene_names = []
    gene_to_idx = _build(gene_names)
    n_genes = len(gene_names)
    dw = _weights_for(assignment_plan, gene_to_idx, n_genes)
    coverage, usable = v0b.compute_coverage(gene_names)
    enrich = v0b.run_enrichment(dw, coverage, usable, n_genes)
    assignments = v0b.winner_take_all(enrich)
    ery = [m for m in v0b.ERY_SUBMODULES if m in usable]
    gran = [m for m in v0b.GRAN_SUBMODULES if m in usable]
    mod = v0b.per_seed_modularity(assignments, ery, gran)
    decision = v0b.decide(mod, np.random.default_rng(0))
    return assignments, mod, decision


def _plan(per_module: dict[str, int], start=0):
    """Build {feature_idx: submodule} planting `n` features per submodule.
    Enough co-significant tests so BH-FDR is not pathological (see audit note on
    small-marker-set conservativeness)."""
    plan, f = {}, start
    for module, n in per_module.items():
        for _ in range(n):
            plan[f] = module
            f += 1
    return plan


def test_assignment_recovers_engineered_modules():
    plan = _plan({"Ery_TF": 12, "Ery_Heme": 12, "Ery_Membrane": 12,
                  "Gran_TF": 12, "Gran_Primary": 12})
    assignments, _, _ = _run(plan)
    s0 = assignments[assignments["seed"] == 0].set_index("feature_idx")["module"]
    for feat, module in plan.items():
        assert s0.loc[feat] == module, f"feature {feat}: {s0.loc[feat]} != {module}"
    # A feature with no planted signal -> UNASSIGNED
    assert s0.loc[127] == "UNASSIGNED"
    print("PASS test_assignment_recovers_engineered_modules")


def test_supported_when_ery_distributed_gran_unified():
    # Erythroid spread across all 3 submodules; granulocyte only in Gran_TF.
    plan = _plan({"Ery_TF": 12, "Ery_Heme": 12, "Ery_Membrane": 12, "Gran_TF": 12})
    _, mod, decision = _run(plan)
    assert decision["n_seeds_both_directions_hold"] == 5, mod.to_string()
    assert decision["asymmetric_modularity_supported"] is True, decision
    print("PASS test_supported_when_ery_distributed_gran_unified")


def test_not_supported_when_ery_unified():
    # Erythroid concentrated in one submodule; granulocyte spread across its two.
    plan = _plan({"Ery_TF": 18, "Gran_TF": 12, "Gran_Primary": 12})
    _, _, decision = _run(plan)
    assert decision["asymmetric_modularity_supported"] is False, decision
    print("PASS test_not_supported_when_ery_unified")


def test_metrics_bounds():
    plan = _plan({"Ery_TF": 12, "Ery_Heme": 12, "Ery_Membrane": 12, "Gran_TF": 12})
    _, mod, _ = _run(plan)
    for col in ["ery_norm_entropy", "gran_norm_entropy"]:
        assert (mod[col] >= -1e-9).all() and (mod[col] <= 1 + 1e-9).all(), col
    print("PASS test_metrics_bounds")


if __name__ == "__main__":
    test_assignment_recovers_engineered_modules()
    test_supported_when_ery_distributed_gran_unified()
    test_not_supported_when_ery_unified()
    test_metrics_bounds()
    print("\nALL V0b LOGIC TESTS PASSED")
