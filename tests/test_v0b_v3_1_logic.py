"""Synthetic tests for the v3.1 control-referenced decision.
Run: python tests/test_v0b_v3_1_logic.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import v0b_module_definitions as v0b
import v0b_v3_loading as v3
import v0b_v3_1_decision as v31
from tests.test_v0b_v3_logic import _build_gene_index, _weights, _plan


def _run(plan):
    names, idx = _build_gene_index()
    n = len(names)
    dw = _weights(plan, idx, n)
    cov, usable = v0b.compute_coverage(names)
    enrich = v3.loading_enrichment(dw, cov, usable, n)
    strength = v3.submodule_strength(enrich, usable)
    ery = [m for m in v0b.ERY_SUBMODULES if m in usable]
    gran = [m for m in v0b.GRAN_SUBMODULES if m in usable]
    ctrl = [m for m in v31.CONTROL_SUBMODULES if m in usable]
    cr = v31.per_seed_control_referenced(strength, ery, gran, ctrl)
    return cr, v31.decide(cr)


def test_supported_when_ery_has_multiple_above_control():
    plan = _plan({"Ery_TF": 3, "Ery_Heme": 3, "Ery_Membrane": 3, "Gran_TF": 3})
    cr, d = _run(plan)
    assert cr["ery_n_real"].median() >= 2, cr.to_string()
    assert d["asymmetric_modularity_supported"] is True, d
    print("PASS test_supported_when_ery_has_multiple_above_control")


def test_not_supported_when_ery_undetected():
    # mirrors the real data: only granulocyte has above-control programs
    plan = _plan({"Gran_Primary": 4, "Gran_TF": 3})
    cr, d = _run(plan)
    assert cr["ery_n_real"].median() == 0, cr.to_string()
    assert d["asymmetric_modularity_supported"] is False, d
    print("PASS test_not_supported_when_ery_undetected")


def test_norm_largest_frac_bounds():
    # uniform -> 0, fully concentrated -> 1
    assert abs(v31._norm_largest_frac([1.0, 1.0, 1.0]) - 0.0) < 1e-9
    assert abs(v31._norm_largest_frac([1.0, 0.0]) - 1.0) < 1e-9
    print("PASS test_norm_largest_frac_bounds")


if __name__ == "__main__":
    test_supported_when_ery_has_multiple_above_control()
    test_not_supported_when_ery_undetected()
    test_norm_largest_frac_bounds()
    print("\nALL V0b v3.1 LOGIC TESTS PASSED")
