"""Synthetic tests for the pre-registered v3.1 control-referenced decision.
Run: python tests/test_v0b_v3_1_logic.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # repo root for `tests.` cross-import
import v0b_module_definitions as v0b
import v0b_v3_loading as v3
import v0b_v3_1_decision as v31
from tests.test_v0b_v3_logic import _build_gene_index, _weights, _plan


def _run(plan, n_perm=300):
    names, idx = _build_gene_index()
    n = len(names)
    dw = _weights(plan, idx, n)
    cov, usable = v0b.compute_coverage(names)
    enrich = v3.loading_enrichment(dw, cov, usable, n)
    strength = v3.submodule_strength(enrich, usable)
    perm = v31.abundance_matched_null(dw, cov, usable, n, n_perm, np.random.default_rng(0))
    ery = [m for m in v0b.ERY_SUBMODULES if m in usable]
    gran = [m for m in v0b.GRAN_SUBMODULES if m in usable]
    cr = v31.per_seed_control_referenced(strength, ery, gran, perm, usable)
    return cr, v31.decide(cr), perm


def test_supported_when_ery_has_multiple_above_control():
    plan = _plan({"Ery_TF": 3, "Ery_Heme": 3, "Ery_Membrane": 3, "Gran_TF": 3})
    cr, d, perm = _run(plan)
    assert cr["ery_n_real"].median() >= 2, cr.to_string()
    assert d["asymmetric_modularity_supported"] is True, d
    # planted submodules are significant under the abundance-matched null
    assert perm["Ery_TF"]["perm_q"] < 0.05, perm["Ery_TF"]
    # Progenitor control is not a real program
    assert perm["Progenitor"]["perm_q"] >= 0.05, perm["Progenitor"]
    print("PASS test_supported_when_ery_has_multiple_above_control")


def test_not_supported_when_ery_undetected():
    plan = _plan({"Gran_Primary": 4, "Gran_TF": 3})
    cr, d, _ = _run(plan)
    assert cr["ery_n_real"].median() == 0, cr.to_string()
    assert d["asymmetric_modularity_supported"] is False, d
    print("PASS test_not_supported_when_ery_undetected")


def test_progenitor_only_baseline():
    # Cycling must NOT be used as the baseline; the control constant is Progenitor.
    assert v31.CONTROL_SUBMODULE == "Progenitor"
    assert "Cycling" in v31.SECONDARY_CONTROLS
    print("PASS test_progenitor_only_baseline")


def test_norm_largest_frac_bounds():
    assert abs(v31._norm_largest_frac([1.0, 1.0, 1.0]) - 0.0) < 1e-9
    assert abs(v31._norm_largest_frac([1.0, 0.0]) - 1.0) < 1e-9
    print("PASS test_norm_largest_frac_bounds")


if __name__ == "__main__":
    test_supported_when_ery_has_multiple_above_control()
    test_not_supported_when_ery_undetected()
    test_progenitor_only_baseline()
    test_norm_largest_frac_bounds()
    print("\nALL V0b v3.1 LOGIC TESTS PASSED")
