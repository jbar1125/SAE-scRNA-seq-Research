#!/usr/bin/env python3
"""Design validation of the frozen historical program-modulation metric.

The generic one-axis world is a counterexample to the claimed interpretation,
not evidence that biology behaves this way. Historical sources are loaded from
the immutable Git commit by the benchmark; current fixes do not erase history.
Run: python3 tests/test_program_modulation_counterexamples.py
"""
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import benchmark_program_modulation as bench
import program_modulation as current_pm
from program_modulation_specificity import program_modulation_specificity


def _run(scenario, n_shuffle=30):
    repo = Path(__file__).resolve().parents[1]
    historical, _ = bench.historical_module(repo)
    acts, expression, labels, genes, markers = bench.make_world(scenario, seed=0)
    result = historical.program_modulation(
        acts, expression, labels, genes, markers,
        top_k=bench.PARAMETERS["top_k"], n_random_features=30,
        n_shuffle=n_shuffle, shuffle_mode="full", seed=0,
    )
    record = result["programs"]["EXTERNAL_PROGRAM"]
    identity = bench.rank_identity_test(acts[:, record["feature"]], labels, 99, 10000)
    diagnostic = program_modulation_specificity(
        acts, labels, result["programs"], n_shuffle=99, seed=10000,
    )
    identity["production_diagnostic"] = diagnostic["programs"]["EXTERNAL_PROGRAM"]
    return record, identity


def test_generic_one_axis_is_legacy_counterexample():
    record, identity = _run("generic_one_axis")
    assert record["feature"] == 0
    assert record["enrich_p_bonf"] < 0.05
    assert record["n_modulators"] == 20
    assert record["shuffle_null"]["empirical_p"] < 0.05
    assert record["feature_background_percentile"] >= 90
    assert identity["permutation_p"] >= 0.05
    assert identity["production_diagnostic"]["matched"] is True
    assert identity["production_diagnostic"]["conditional_heterogeneity"] is False
    print("generic one-axis: legacy quantitative gates pass despite no identity effect")


def test_specific_subset_has_identity_signal():
    record, identity = _run("specific_subset")
    assert record["feature"] == 0
    assert set(record["modulators"]) == {f"P{i:02d}" for i in range(4)}
    assert identity["permutation_p"] < 0.05
    assert identity["production_diagnostic"]["conditional_heterogeneity"] is True
    print("specific subset: four planted modulators and identity heterogeneity recovered")


def test_null_world_does_not_pass_legacy_gate():
    record, identity = _run("null")
    assert record["n_modulators"] == 0
    assert record["shuffle_null"]["empirical_p"] == 1.0
    assert identity["permutation_p"] >= 0.05
    print("null world: no modulators and no identity signal at the fixed oracle seed")


def test_historical_no_marker_fallback_is_exposed():
    repo = Path(__file__).resolve().parents[1]
    historical, _ = bench.historical_module(repo)
    prog = np.array([[0, 0, 0, 3, 2, 1], [0, 0, 0, 2, 3, 1]], dtype=float)
    assert historical.select_program_feature(prog, [0, 1, 2], 6, 2, [0, 1]) == (-1, 1.0, 1.0)
    record, identity = _run("absent_markers_last_feature")
    assert record["feature"] == -1
    assert record["enrich_p_bonf"] == 1.0
    assert record["n_modulators"] == 20
    assert record["shuffle_null"]["empirical_p"] < 0.05
    assert record["feature_background_percentile"] >= 90
    assert identity["production_diagnostic"]["matched"] is False
    print("historical absent markers: index -1 silently scores an unrelated last feature")


def test_current_module_rejects_absent_and_weak_marker_matches():
    for scenario in ("absent_markers_last_feature", "generic_one_axis"):
        acts, expression, labels, genes, markers = bench.make_world(scenario, seed=0)
        if scenario == "generic_one_axis":
            # One of five markers is correlated; the other four are not expressed.
            # The historical selector has a real index but Bonferroni p=1.
            expression[:, 1:5] = 0.0
        output = current_pm.program_modulation(
            acts, expression, labels, genes, markers, top_k=5,
            n_shuffle=0, n_random_features=5, seed=0,
        )
        record = output["programs"]["EXTERNAL_PROGRAM"]
        assert record["low_power"] is True
        assert record["feature"] is None
        assert record["reason"] == "no_significant_program_match"
        assert not record.get("modulators", [])
        assert record.get("n_modulators", 0) == 0
        assert record["enrich_p_bonf"] > 0.05
    print("current module: absent and weak marker matches are explicitly untested")


def test_current_empty_eligible_set_preserves_planned_programs():
    acts, expression, labels, genes, markers = bench.make_world("null", seed=0)
    acts[:] = 0.0
    markers["SECOND_PLANNED_PROGRAM"] = genes[5:10]
    output = current_pm.program_modulation(
        acts, expression, labels, genes, markers, top_k=5, n_shuffle=0,
        n_random_features=5, seed=0,
    )
    assert set(output["programs"]) == set(markers)
    for record in output["programs"].values():
        assert record["low_power"] is True
        assert record["reason"] == "no_eligible_features"
    print("current empty feature set: every planned program is retained")


def test_benchmark_refuses_metadata_overwrite():
    script = Path(bench.__file__).resolve()
    with tempfile.TemporaryDirectory() as temp:
        directory = Path(temp)
        metadata = directory / "metadata.json"
        original = b'{"original_artifact": true}\n'
        metadata.write_bytes(original)
        result = subprocess.run(
            [sys.executable, str(script), "--out-dir", str(directory),
             "--seeds", "1", "--n-shuffle", "1"],
            text=True, capture_output=True,
        )
        assert result.returncode != 0
        assert metadata.read_bytes() == original
        assert not (directory / ".benchmark_claim").exists()
        assert not (directory / "results.jsonl").exists()
        assert "choose a fresh --out-dir" in result.stderr
    print("benchmark repeat: existing metadata is preserved before any output write")


if __name__ == "__main__":
    test_generic_one_axis_is_legacy_counterexample()
    test_specific_subset_has_identity_signal()
    test_null_world_does_not_pass_legacy_gate()
    test_historical_no_marker_fallback_is_exposed()
    test_current_module_rejects_absent_and_weak_marker_matches()
    test_current_empty_eligible_set_preserves_planned_programs()
    test_benchmark_refuses_metadata_overwrite()
    print("ALL PROGRAM-MODULATION COUNTEREXAMPLE TESTS PASSED")
