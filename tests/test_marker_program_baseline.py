"""Meaningful counterexamples for the independently declared marker-score baseline."""
import os
import sys

import numpy as np
from scipy import sparse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from marker_program_baseline import (marker_program_baseline, normalize_counts_to_markers,
                                     split_controls, sum_duplicate_gene_counts)  # noqa: E402


def _balanced_null():
    rng = np.random.default_rng(4)
    labels = np.array(["control"] * 200 + ["T"] * 80 + ["P1"] * 80)
    genes = np.array(["T", "g1", "g2", "g3"])
    discovery, inference = split_controls(labels)
    x = np.empty((len(labels), 4))
    x[discovery] = rng.uniform(1, 3, size=(len(discovery), 4))
    template = rng.uniform(1, 3, size=(20, 4))
    x[inference] = np.tile(template, (5, 1))
    x[labels == "T"] = x[labels == "P1"] = np.tile(template, (4, 1))
    markers = {f"program{j}": genes.tolist() for j in range(9)}
    return x, genes, labels, markers


def test_null_and_full_planned_family():
    x, genes, labels, panels = _balanced_null()
    out = marker_program_baseline(x, genes, labels, panels, n_identity_shuffle=0)
    assert out["family_size"] == 18 and out["n_perturbations"] == 2
    assert not any(r["modulates"] for p in out["programs"].values() for r in p["per_perturbation"])
    assert all(r["auc"] == .5 for p in out["programs"].values() for r in p["per_perturbation"])
    for j in range(1, 9):
        panels[f"program{j}"] = ["missing1", "missing2", "missing3"]
    missing = marker_program_baseline(x, genes, labels, panels, n_identity_shuffle=0)
    assert missing["family_size"] == 18
    assert sum(len(p["per_perturbation"]) for p in missing["programs"].values()) == 18
    assert all(r["mw_p"] == r["mw_q"] == 1 for name, p in missing["programs"].items()
               if name != "program0" for r in p["per_perturbation"])


def test_self_target_artifact_removed_and_omnibus_is_conditional():
    x, genes, labels, panels = _balanced_null()
    x[labels == "T", 0] += 20
    out = marker_program_baseline(x, genes, labels, panels, n_identity_shuffle=199)
    row = next(r for r in out["programs"]["program0"]["per_perturbation"] if r["pert"] == "T")
    assert row["modulates"] and row["auc"] == 1
    assert row["target_omission"]["target_removed"]
    assert row["target_omission"]["n_markers_remaining"] == 3
    assert row["target_omission"]["auc"] == .5
    assert not row["target_omission"]["modulates"] and not row["survives_target_omission"]
    identity = out["identity_heterogeneity"]
    assert identity["n_planned_programs"] == 9
    assert identity["programs"]["program0"]["conditional_heterogeneity"]
    assert "not a measured hypergeometric" in identity["external_marker_identity_note"]


def test_discovery_fit_does_not_use_inference_controls_or_perturbations():
    x, genes, labels, panels = _balanced_null()
    first = marker_program_baseline(x, genes, labels, panels, n_identity_shuffle=0)
    discovery, inference = split_controls(labels)
    changed = x.copy()
    changed[inference] = 100
    changed[labels != "control"] += 50
    second = marker_program_baseline(changed, genes, labels, panels, n_identity_shuffle=0)
    assert first["discovery_gene_parameters"] == second["discovery_gene_parameters"]
    assert first["control_split"] == second["control_split"]
    assert set(first["control_split"]["discovery_cell_ids"]).isdisjoint(first["control_split"]["inference_cell_ids"])
    assert len(discovery) == len(inference) == 100
    assert all(r["n_controls"] == 100 for p in first["programs"].values() for r in p["per_perturbation"])


def test_strata_remove_compositional_batch_association_and_no_overlap_is_untestable():
    labels = np.array(["control"] * 200 + ["T"] * 80)
    blocks = np.array(["high"] * 100 + ["low"] * 100 + ["high"] * 70 + ["low"] * 10)
    x = np.repeat(np.where(blocks == "high", 10., 1.)[:, None], 4, axis=1)
    genes = np.array(["T", "g1", "g2", "g3"])
    panels = {f"program{j}": genes.tolist() for j in range(9)}
    pooled = marker_program_baseline(x, genes, labels, panels, n_identity_shuffle=0)
    stratified = marker_program_baseline(x, genes, labels, panels, strata=blocks, n_identity_shuffle=0)
    assert pooled["programs"]["program0"]["per_perturbation"][0]["modulates"]
    row = stratified["programs"]["program0"]["per_perturbation"][0]
    assert row["mw_p"] == 1 and row["auc"] == .5 and not row["modulates"]
    blocks[labels == "T"] = "unseen"
    no_overlap = marker_program_baseline(x, genes, labels, panels, strata=blocks, n_identity_shuffle=0)
    row = no_overlap["programs"]["program0"]["per_perturbation"][0]
    assert row["status"] == "insufficient_comparison_cells" and row["mw_p"] == 1 and not row["modulates"]


def test_raw_duplicate_sum_and_all_gene_normalization():
    counts = sparse.csr_matrix([[1, 2, 3, 4, 90], [2, 3, 4, 5, 86]], dtype=float)
    genes = np.array(["A", "A", "B", "C", "background"])
    summed, unique = sum_duplicate_gene_counts(counts, genes)
    assert unique.tolist() == ["A", "B", "C", "background"]
    assert np.array_equal(summed.toarray(), [[3, 3, 4, 90], [5, 4, 5, 86]])
    assert np.array_equal(np.asarray(counts.sum(axis=1)), np.asarray(summed.sum(axis=1)))
    panels = {"program": ["A", "B", "C"]}
    try:
        normalize_counts_to_markers(counts, genes, panels)
    except ValueError as exc:
        assert "duplicate" in str(exc)
    else:
        raise AssertionError("duplicate symbols were silently accepted")
    expression, kept_genes, note = normalize_counts_to_markers(counts, genes, panels, duplicate_genes="sum")
    assert kept_genes.tolist() == ["A", "B", "C"]
    assert np.allclose(expression, np.log1p(np.array([[3, 3, 4], [5, 4, 5]]) * 100))
    assert note["duplicate_columns_aggregated"] == 1 and note["n_hvg"] == 0
    # Normalizing the marker subset first would instead use library totals 10/14.
    assert not np.isclose(expression[0, 0], np.log1p(3 / 10 * 10000))
    try:
        normalize_counts_to_markers(sparse.csr_matrix([[.1, .2, .3]]), ["A", "B", "C"], panels)
    except ValueError as exc:
        assert "raw integer counts" in str(exc)
    else:
        raise AssertionError("normalized counts were accepted")


def test_invalid_inputs_and_target_omission_low_power():
    x, genes, labels, panels = _balanced_null()
    invalid = [
        dict(expression=x, gene_names=["T", "T", "g2", "g3"], labels=labels, marker_sets=panels),
        dict(expression=x, gene_names=genes, labels=labels[:-1], marker_sets=panels),
        dict(expression=np.full_like(x, np.nan), gene_names=genes, labels=labels, marker_sets=panels),
        dict(expression=x, gene_names=genes, labels=np.full(len(labels), "T"), marker_sets=panels),
        dict(expression=x, gene_names=genes, labels=labels, marker_sets=panels, cell_ids=["same"] * len(labels)),
    ]
    for kwargs in invalid:
        try:
            marker_program_baseline(**kwargs, n_identity_shuffle=0)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid inputs accepted")
    panels = {"three_markers": ["T", "g1", "g2"]}
    x[labels == "T", 0] += 20
    out = marker_program_baseline(x, genes, labels, panels, n_identity_shuffle=0)
    row = next(r for r in out["programs"]["three_markers"]["per_perturbation"] if r["pert"] == "T")
    assert row["modulates"]
    assert row["target_omission"]["status"] == "low_power_after_target_omission"
    assert not row["survives_target_omission"]


if __name__ == "__main__":
    test_null_and_full_planned_family()
    test_self_target_artifact_removed_and_omnibus_is_conditional()
    test_discovery_fit_does_not_use_inference_controls_or_perturbations()
    test_strata_remove_compositional_batch_association_and_no_overlap_is_untestable()
    test_raw_duplicate_sum_and_all_gene_normalization()
    test_invalid_inputs_and_target_omission_low_power()
    print("ALL MARKER-PROGRAM BASELINE TESTS PASSED (6 checks)")
