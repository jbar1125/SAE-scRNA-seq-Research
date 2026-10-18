"""Counterexamples and calibration checks for the separate conditional identity test."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from program_modulation_specificity import program_modulation_specificity  # noqa: E402


def _world(kind="specific", seed=0):
    rng = np.random.default_rng(seed)
    labels = np.array(["control"] * 100 + [f"P{p}" for p in range(8) for _ in range(60)])
    # Every identity has exactly the same empirical response in the generic world.
    common = np.linspace(-0.8, 0.8, 60)
    scores = np.r_[rng.normal(0, 0.2, 100), np.tile(common, 8)]
    if kind == "specific":
        scores[np.isin(labels, ["P0", "P1"])] += 3
    elif kind == "generic":
        scores[labels != "control"] += 3
    a = np.c_[scores, rng.normal(size=len(labels))]
    return a, labels, {"ERY": {"feature": 0, "enrich_p_bonf": 1e-6}}


def test_specific_subset_and_equal_generic_response():
    a, labels, records = _world()
    specific = program_modulation_specificity(a, labels, records, n_shuffle=199, seed=3)["programs"]["ERY"]
    assert specific["conditional_heterogeneity"] and specific["empirical_p"] == 0.005
    a, labels, records = _world("generic")
    generic = program_modulation_specificity(a, labels, records, n_shuffle=199, seed=3)["programs"]["ERY"]
    assert generic["empirical_p"] == 1 and not generic["conditional_heterogeneity"]


def test_no_match_retained_and_missing_match_not_imputed():
    a, labels, records = _world()
    records.update({"NO_HITS": {"feature": -1, "enrich_p_bonf": 1.0},
                    "WEAK": {"feature": 1, "enrich_p_bonf": 0.8},
                    "UNKNOWN": {"feature": 1}, "LOW_POWER": {"low_power": True}})
    out = program_modulation_specificity(a, labels, records, n_shuffle=199, seed=3)
    assert out["n_planned_programs"] == 5
    for name in records.keys() - {"ERY"}:
        rec = out["programs"][name]
        assert rec["empirical_p"] == 1 and not rec["matched"] and not rec["conditional_heterogeneity"]
    assert out["programs"]["ERY"]["holm_p"] == 0.025


def test_duplicate_feature_is_not_duplicate_evidence_and_holm_family():
    a, labels, records = _world()
    records["ERY_ALIAS"] = dict(records["ERY"])
    out = program_modulation_specificity(a, labels, records, n_shuffle=19, seed=3)
    first, second = out["programs"].values()
    assert first["empirical_p"] == second["empirical_p"] == 0.05
    assert first["holm_p"] == second["holm_p"] == 0.1
    assert not first["conditional_heterogeneity"] and not second["conditional_heterogeneity"]
    assert first["between_group_rank_ss"] == second["between_group_rank_ss"]


def test_determinism_and_controls_combined_low_count_excluded():
    a, labels, records = _world()
    first = program_modulation_specificity(a, labels, records, n_shuffle=99, seed=4)
    assert first == program_modulation_specificity(a, labels, records, n_shuffle=99, seed=4)
    # Controls, combined perturbations, and groups below min_cells do not enter identity shuffles.
    extra = np.tile([1e8, -1e8], (60, 1))
    a2 = np.r_[a, extra]
    l2 = np.r_[labels, ["P0+P1"] * 50 + ["TOO_SMALL"] * 10]
    a2[labels.tolist().count("control") - 1, 0] = 1e7
    second = program_modulation_specificity(a2, l2, records, n_shuffle=99, seed=4)
    assert first["programs"] == second["programs"]
    assert first["perturbation_cell_counts"] == second["perturbation_cell_counts"]


def test_strata_block_confounding_and_shared_batch_effect():
    a, labels, records = _world()
    # Each identity occurs in its own stratum: no conditional identity contrasts exist.
    out = program_modulation_specificity(a, labels, records, strata=labels, n_shuffle=99)
    assert out["n_exchangeable_cells"] == 0
    rec = out["programs"]["ERY"]
    assert rec["status"] == "identity_confounded_with_strata" and rec["empirical_p"] == 1

    # Strong batch effect plus unequal batch composition creates an unstratified
    # association. Within-batch scores are identical, so the conditional test is null.
    rng = np.random.default_rng(4)
    labs, block, values = ["control"] * 50, ["ctrl"] * 50, rng.normal(size=50).tolist()
    for p, n_hi in (("P0", 70), ("P1", 10)):
        labs += [p] * 80
        block += ["high"] * n_hi + ["low"] * (80 - n_hi)
        values += [10.] * n_hi + [0.] * (80 - n_hi)
    mat = np.array(values)[:, None]
    plain = program_modulation_specificity(mat, labs, records, n_shuffle=199, seed=7)["programs"]["ERY"]
    strat = program_modulation_specificity(mat, labs, records, strata=block, n_shuffle=199, seed=7)["programs"]["ERY"]
    assert plain["conditional_heterogeneity"]
    assert strat["empirical_p"] == 1 and not strat["conditional_heterogeneity"]


def test_invalid_inputs_and_constant_feature():
    a, labels, records = _world()
    invalid = [
        dict(activations=np.full_like(a, np.nan), labels=labels, program_records=records),
        dict(activations=a, labels=labels[:-1], program_records=records),
        dict(activations=a, labels=[None] * len(labels), program_records=records),
        dict(activations=a, labels=np.full(len(labels), "P0"), program_records=records),
        dict(activations=a, labels=labels, program_records={"X": {"feature": 99, "enrich_p_bonf": .01}}),
        dict(activations=a, labels=labels, program_records={"X": {"feature": 0, "enrich_p_bonf": np.nan}}),
        dict(activations=a, labels=labels, program_records=records, strata=["batch"]),
        dict(activations=a, labels=labels, program_records=records, n_shuffle=0),
    ]
    for kwargs in invalid:
        try:
            program_modulation_specificity(**kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError(f"invalid input was accepted: {list(kwargs)}")
    a[:] = 1
    rec = program_modulation_specificity(a, labels, records, n_shuffle=19)["programs"]["ERY"]
    assert rec["status"] == "constant_feature" and rec["empirical_p"] == 1


if __name__ == "__main__":
    test_specific_subset_and_equal_generic_response()
    test_no_match_retained_and_missing_match_not_imputed()
    test_duplicate_feature_is_not_duplicate_evidence_and_holm_family()
    test_determinism_and_controls_combined_low_count_excluded()
    test_strata_block_confounding_and_shared_batch_effect()
    test_invalid_inputs_and_constant_feature()
    print("ALL PROGRAM-MODULATION SPECIFICITY TESTS PASSED (6 checks)")
