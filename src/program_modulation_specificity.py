"""Conditional perturbation-identity heterogeneity for a control-selected feature.

This is a separate validation analysis, not a replacement for the frozen exploratory
program-modulation metric. A common response shared by every perturbation is not
identity heterogeneity. The null here therefore keeps controls fixed and permutes
identities only among eligible single-perturbation cells, optionally within strata.

The statistic is the weighted between-group sum of squares of pooled midranks.
Feature selection must use controls only and the fitted representation must be
invariant to perturbation-identity permutations. A label-blind all-cell fit satisfies
that invariance. Selection on these perturbation identities does not.

Inference is conditional on the observed cells and the specified exchangeability
strata. A rejection establishes rank-location heterogeneity among identities; it
does not establish pathway specificity, causal mechanism, biological replication,
SAE superiority, or a novel discovery. Unequal generic response severity can reject.
"""
from __future__ import annotations

from collections.abc import Mapping

import numpy as np
from scipy.stats import rankdata


def _positive_integer(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def _probability(value, name, *, allow_zero=False):
    if isinstance(value, (bool, np.bool_)) or not np.isscalar(value):
        raise ValueError(f"{name} must be a finite probability")
    try:
        value = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite probability") from exc
    lower_ok = value >= 0 if allow_zero else value > 0
    upper_ok = value <= 1 if allow_zero else value < 1
    if not np.isfinite(value) or not lower_ok or not upper_ok:
        raise ValueError(f"{name} must be a finite probability")
    return value


def _strings(values, name, n_rows):
    raw = np.asarray(values, dtype=object)
    if raw.ndim != 1 or len(raw) != n_rows:
        raise ValueError(f"{name} must be one-dimensional with one value per cell")
    if any(not isinstance(v, (str, np.str_)) or not str(v).strip()
           or str(v).strip().lower() in {"nan", "none", "null", "<na>", "nat"} for v in raw):
        raise ValueError(f"{name} must contain nonempty, nonmissing strings")
    return raw.astype(str)


def _holm(pvalues):
    """Holm family-wise adjusted p values, retaining every planned hypothesis."""
    p = np.asarray(pvalues, float)
    if not len(p):
        return p
    order = np.argsort(p, kind="stable")
    adjusted = np.minimum(1.0, np.maximum.accumulate(p[order] * np.arange(len(p), 0, -1)))
    out = np.empty_like(p)
    out[order] = adjusted
    return out


def program_modulation_specificity(activations, labels, program_records, *, strata=None,
                                  control_label="control", min_cells=20, n_shuffle=999,
                                  alpha=0.05, match_alpha=0.05, seed=0):
    """Test conditional identity heterogeneity for all planned named programs.

    ``program_records`` is the complete ``legacy_result['programs']`` mapping, not
    only the promising programs. Each matched record needs ``feature`` and
    ``enrich_p_bonf`` from control-only selection. ``low_power`` records, absent
    matches, feature -1, and marker-match p > ``match_alpha`` receive p=1 and stay
    in the Holm family. That marker gate is a diagnostic identity criterion; its
    hypergeometric p value is not itself a calibrated biological discovery test.

    ``strata``, when supplied, contains a nonempty string per original cell, such
    as batch or biological replicate. Identity assignments are shuffled within
    strata while preserving every group-by-stratum cell count. No permutation can
    resolve complete confounding of identity with strata; such analyses return
    p=1 and ``status='identity_confounded_with_strata'``. This function cannot
    infer which experimental units are exchangeable: callers must supply valid
    strata and distinguish conditional cell-level inference from replication.

    All programs use the same permutations. Duplicate selected features reuse the
    same computed result but remain separate planned hypotheses for correction.
    The returned ``conditional_heterogeneity`` is true only after Holm correction.
    """
    a = np.asarray(activations)
    if (a.ndim != 2 or not a.shape[0] or not a.shape[1]
            or not np.issubdtype(a.dtype, np.number) or np.iscomplexobj(a)):
        raise ValueError("activations must be a nonempty finite two-dimensional array")
    # Do not duplicate a full, potentially multi-gigabyte float32 activation matrix.
    if any(not np.isfinite(a[start:start + 8192]).all() for start in range(0, len(a), 8192)):
        raise ValueError("activations must be a nonempty finite two-dimensional array")
    lab = _strings(labels, "labels", a.shape[0])
    if not isinstance(control_label, str) or not control_label or not np.any(lab == control_label):
        raise ValueError("control_label must identify an observed control group")
    min_cells = _positive_integer(min_cells, "min_cells")
    n_shuffle = _positive_integer(n_shuffle, "n_shuffle")
    alpha = _probability(alpha, "alpha")
    match_alpha = _probability(match_alpha, "match_alpha")
    if isinstance(seed, (bool, np.bool_)) or not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("seed must be a nonnegative integer")
    seed = int(seed)
    if not isinstance(program_records, Mapping) or not program_records:
        raise ValueError("program_records must contain all planned named programs")
    if any(not isinstance(name, str) or not name for name in program_records):
        raise ValueError("program names must be nonempty strings")
    if strata is None:
        block_labels = np.full(a.shape[0], "all", dtype=str)
    else:
        block_labels = _strings(strata, "strata", a.shape[0])

    # Match validation occurs even when there are too few usable perturbations.
    records = {}
    for name, source in program_records.items():
        if not isinstance(source, Mapping):
            raise ValueError(f"program record {name!r} must be a mapping")
        rec = {"matched": False, "status": "no_marker_match", "empirical_p": 1.0,
               "holm_p": 1.0, "conditional_heterogeneity": False}
        feature = source.get("feature")
        p_match = source.get("enrich_p_bonf")
        if feature is not None:
            if isinstance(feature, (bool, np.bool_)) or not isinstance(feature, (int, np.integer)):
                raise ValueError(f"feature for {name!r} must be an integer or absent")
            if feature < -1 or feature >= a.shape[1]:
                raise ValueError(f"feature for {name!r} is outside the activation array")
        if p_match is not None:
            p_match = _probability(p_match, f"enrich_p_bonf for {name!r}", allow_zero=True)
        if source.get("low_power", False):
            rec["status"] = "low_power_marker_panel"
        elif feature is not None and feature >= 0 and p_match is not None and p_match <= match_alpha:
            rec.update(matched=True, status="pending", feature=int(feature),
                       marker_match_p_bonf=float(p_match))
        records[name] = rec

    eligible_names = [name for name in np.unique(lab)
                      if name != control_label and "+" not in name and np.sum(lab == name) >= min_cells]
    eligible_mask = np.isin(lab, eligible_names)
    eligible_rows = np.flatnonzero(eligible_mask)
    groups = lab[eligible_mask]
    _, group_ids = np.unique(groups, return_inverse=True)
    counts = np.bincount(group_ids, minlength=len(eligible_names))
    blocks = block_labels[eligible_mask]
    movable_blocks = []
    for block in np.unique(blocks):
        rows = np.flatnonzero(blocks == block)
        if len(np.unique(group_ids[rows])) > 1:
            movable_blocks.append(rows)
    movable_rows = (np.concatenate(movable_blocks) if movable_blocks else np.array([], dtype=int))
    movable_groups = set(groups[movable_rows])
    out = {
        "analysis": "conditional_cell_level_perturbation_identity_rank_heterogeneity",
        "planned_programs": list(program_records), "n_planned_programs": len(program_records),
        "n_eligible_perturbations": len(eligible_names), "n_eligible_cells": int(len(eligible_rows)),
        "perturbation_cell_counts": {name: int(n) for name, n in zip(eligible_names, counts)},
        "n_permutations_requested": n_shuffle, "n_permutations": 0, "seed": seed, "alpha": float(alpha),
        "marker_match_alpha": float(match_alpha), "multiplicity": "Holm over all planned programs",
        "conditioning": "Original control membership, selected features, eligible groups, and group-by-stratum counts fixed",
        "stratified": strata is not None, "n_strata": int(len(np.unique(blocks))),
        "n_movable_strata": len(movable_blocks), "n_exchangeable_cells": int(len(movable_rows)),
        "perturbations_without_exchangeable_cells": [name for name in eligible_names if name not in movable_groups],
        "scope": "Observed-cell rank-location heterogeneity; does not establish pathway specificity, mechanism, novelty, or biological replication",
        "programs": records,
    }
    reason = None
    if len(eligible_names) < 2:
        reason = "fewer_than_two_eligible_perturbations"
    elif not movable_blocks:
        reason = "identity_confounded_with_strata"
    if reason:
        for rec in records.values():
            if rec["matched"]:
                rec["status"] = reason
        return out

    features = sorted({r["feature"] for r in records.values() if r["matched"]})
    if not features:
        return out
    scores = rankdata(a[np.ix_(eligible_rows, features)], axis=0, method="average")
    scores -= scores.mean(axis=0, keepdims=True)
    total_ss = np.sum(scores * scores, axis=0)

    def statistic(ranks):
        return np.array([np.sum(np.bincount(group_ids, weights=ranks[:, j], minlength=len(counts)) ** 2 / counts)
                         for j in range(ranks.shape[1])])

    observed = statistic(scores)
    greater_equal = np.zeros(len(features), dtype=int)
    null_sum = np.zeros(len(features))
    null_max = np.zeros(len(features))
    rng = np.random.default_rng(seed)
    for _ in range(n_shuffle):
        permutation = np.arange(len(eligible_rows))
        for rows in movable_blocks:
            permutation[rows] = rng.permutation(rows)
        null = statistic(scores[permutation])
        # Summation order is constant. Tiny roundoff must not turn exact ties into evidence.
        greater_equal += null >= observed - np.maximum(1.0, np.abs(observed)) * 1e-12
        null_sum += null
        null_max = np.maximum(null_max, null)
    out["n_permutations"] = n_shuffle
    cache = {}
    for j, feature in enumerate(features):
        cache[feature] = {
            "status": "constant_feature" if total_ss[j] == 0 else "tested",
            "between_group_rank_ss": float(observed[j]),
            "rank_eta_squared": float(observed[j] / total_ss[j]) if total_ss[j] > 0 else 0.0,
            "empirical_p": float((greater_equal[j] + 1) / (n_shuffle + 1)),
            "null_mean": float(null_sum[j] / n_shuffle), "null_max": float(null_max[j]),
        }
    for rec in records.values():
        if rec["matched"]:
            rec.update(cache[rec["feature"]])
    adjusted = _holm([rec["empirical_p"] for rec in records.values()])
    for rec, adjusted_p in zip(records.values(), adjusted):
        rec["holm_p"] = float(adjusted_p)
        rec["conditional_heterogeneity"] = bool(rec["status"] == "tested" and adjusted_p <= alpha)
    return out
