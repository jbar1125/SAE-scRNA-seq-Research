#!/usr/bin/env python3
"""Heldout-control, external-marker baseline for Norman perturbation associations.

The plan is config/marker_program_baseline_plan.json. Scores use discovery controls
only; comparisons use the other controls. All planned programs and eligible
perturbations enter BH correction, including untestable hypotheses with p=1.
Results are conditional cell-level associations, not mechanisms or discoveries.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import sys

import numpy as np
from scipy import sparse
from scipy.stats import mannwhitneyu, norm, rankdata
from statsmodels.stats.multitest import multipletests


def _strings(values, name, n=None):
    values = np.asarray(values, dtype=object)
    if values.ndim != 1 or (n is not None and len(values) != n):
        raise ValueError(f"{name} must be a one-dimensional aligned array")
    if any(not isinstance(v, str) or not v.strip() or v.lower() in {"nan", "none", "<na>"} for v in values):
        raise ValueError(f"{name} must contain nonmissing strings")
    return values.astype(str)


def split_controls(labels, *, strata=None, seed=0, control_label="control"):
    """Outcome-blind, within-stratum 50/50 control split; discovery never tested."""
    labels = _strings(labels, "labels")
    blocks = np.full(len(labels), "all") if strata is None else _strings(strata, "strata", len(labels))
    if isinstance(seed, (bool, np.bool_)) or not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("seed must be a nonnegative integer")
    rng = np.random.default_rng(seed)
    discovery, inference = [], []
    for block in np.unique(blocks):
        controls = np.flatnonzero((labels == control_label) & (blocks == block))
        shuffled = rng.permutation(controls)
        cut = len(shuffled) // 2
        discovery.extend(shuffled[:cut].tolist())
        inference.extend(shuffled[cut:].tolist())
    return np.asarray(sorted(discovery), int), np.asarray(sorted(inference), int)


def sum_duplicate_gene_counts(counts, gene_names):
    """Sparse sum of raw columns by symbol, preserving exact library totals."""
    genes = _strings(gene_names, "gene_names")
    if counts.shape[1] != len(genes):
        raise ValueError("counts and genes are not aligned")
    matrix = sparse.csr_matrix(counts)
    _validate_counts(matrix)
    lookup, unique, mapping = {}, [], []
    for gene in genes:
        if gene not in lookup:
            lookup[gene] = len(unique)
            unique.append(gene)
        mapping.append(lookup[gene])
    grouping = sparse.csr_matrix((np.ones(len(genes), dtype=np.int64),
                                  (np.arange(len(genes)), mapping)), shape=(len(genes), len(unique)))
    # Float64 prevents float32 addition-order loss while retaining integer counts.
    result = matrix.astype(np.float64) @ grouping
    before = np.asarray(matrix.sum(axis=1, dtype=np.float64)).ravel()
    after = np.asarray(result.sum(axis=1, dtype=np.float64)).ravel()
    if not np.array_equal(before, after):
        raise ValueError("duplicate aggregation did not preserve exact library totals")
    return result.tocsr(), np.asarray(unique)


def _validate_counts(matrix):
    if not sparse.issparse(matrix):
        matrix = np.asarray(matrix)
    values = matrix.data if sparse.issparse(matrix) else matrix.ravel()
    for start in range(0, len(values), 1_000_000):
        part = values[start:start + 1_000_000]
        if not np.isfinite(part).all() or (part < 0).any() or not np.allclose(part, np.rint(part), rtol=0, atol=1e-6):
            raise ValueError("X must contain finite nonnegative raw integer counts; normalized X is not accepted")
    totals = np.asarray(matrix.sum(axis=1, dtype=np.float64)).ravel()
    if not np.isfinite(totals).all() or (totals <= 0).any():
        raise ValueError("every cell must have a finite positive raw-count library total")


def normalize_counts_to_markers(counts, gene_names, marker_sets, *, duplicate_genes="reject"):
    """Normalize/log all raw genes sparsely, then densify external marker columns."""
    import anndata as ad
    import scanpy as sc

    genes = _strings(gene_names, "gene_names")
    if counts.ndim != 2 or counts.shape[1] != len(genes):
        raise ValueError("counts and genes are not aligned")
    if duplicate_genes not in {"reject", "sum"}:
        raise ValueError("duplicate_genes must be reject or sum")
    matrix = sparse.csr_matrix(counts)
    _validate_counts(matrix)
    n_duplicate_columns = len(genes) - len(set(genes))
    if n_duplicate_columns and duplicate_genes == "reject":
        raise ValueError("duplicate gene symbols require explicit --duplicate-genes sum")
    if n_duplicate_columns:
        matrix, genes = sum_duplicate_gene_counts(matrix, genes)
    else:
        matrix = matrix.astype(np.float64, copy=True)
    annotated = ad.AnnData(X=matrix)
    sc.pp.normalize_total(annotated, target_sum=1e4)
    sc.pp.log1p(annotated)
    wanted = {gene for markers in marker_sets.values() for gene in markers}
    keep = np.array([gene in wanted for gene in genes])
    expression = np.asarray(annotated.X[:, keep].toarray(), dtype=np.float64)
    return expression, genes[keep], {
        "input_genes": int(counts.shape[1]), "unique_genes": len(genes),
        "marker_columns": int(keep.sum()), "duplicate_columns_aggregated": int(n_duplicate_columns),
        "duplicate_policy": duplicate_genes, "normalization": "all-gene normalize_total(10000), then log1p, then marker subset",
        "n_hvg": 0,
    }


def _compare_up(scores, pert_rows, control_rows, strata, min_cells):
    """One-sided pooled or stratified tie-corrected Mann-Whitney comparison."""
    result = {"auc": 0.5, "mw_p": 1.0, "n_cells": 0, "n_controls": 0,
              "status": "insufficient_comparison_cells"}
    if strata is None:
        if len(pert_rows) < min_cells or len(control_rows) < 20:
            return result
        a, b = scores[pert_rows], scores[control_rows]
        u, p = mannwhitneyu(a, b, alternative="greater", method="asymptotic", use_continuity=True)
        return dict(auc=float(u / (len(a) * len(b))), mw_p=float(p), n_cells=len(a),
                    n_controls=len(b), status="tested")
    u_total = expectation = variance = pair_count = 0.0
    used_p = used_c = 0
    for block in np.intersect1d(strata[pert_rows], strata[control_rows]):
        a = scores[pert_rows[strata[pert_rows] == block]]
        b = scores[control_rows[strata[control_rows] == block]]
        n1, n0 = len(a), len(b)
        ranks = rankdata(np.r_[a, b], method="average")
        u = ranks[:n1].sum() - n1 * (n1 + 1) / 2
        _, ties = np.unique(np.r_[a, b], return_counts=True)
        n = n1 + n0
        tie_term = np.sum(ties.astype(float) ** 3 - ties) / (n * (n - 1))
        u_total += u
        expectation += n1 * n0 / 2
        variance += n1 * n0 * ((n + 1) - tie_term) / 12
        pair_count += n1 * n0
        used_p += n1
        used_c += n0
    result.update(n_cells=used_p, n_controls=used_c)
    if used_p < min_cells or used_c < 20 or pair_count == 0:
        return result
    p = float(norm.sf((u_total - expectation - 0.5) / np.sqrt(variance))) if variance > 0 else 1.0
    result.update(auc=float(u_total / pair_count), mw_p=p, status="tested")
    return result


def marker_program_baseline(expression, gene_names, labels, marker_sets, *, strata=None,
                            cell_ids=None, control_label="control", seed=0, min_cells=20,
                            min_markers=3, min_control_nonzero=10, scale_floor=0.1,
                            alpha=0.05, auc_floor=0.55, n_identity_shuffle=999):
    """Analyze already log-normalized external-marker columns with heldout controls."""
    x = np.asarray(expression, dtype=float)
    if x.ndim != 2 or not len(x) or not np.isfinite(x).all() or (x < 0).any():
        raise ValueError("expression must be a finite nonnegative log-expression matrix")
    genes = _strings(gene_names, "gene_names", x.shape[1])
    if len(set(genes)) != len(genes):
        raise ValueError("gene_names must be unique; aggregate raw counts before normalization")
    labels = _strings(labels, "labels", len(x))
    blocks = None if strata is None else _strings(strata, "strata", len(x))
    ids = np.array([str(i) for i in range(len(x))]) if cell_ids is None else _strings(cell_ids, "cell_ids", len(x))
    if len(set(ids)) != len(ids):
        raise ValueError("cell_ids must be unique")
    if not isinstance(marker_sets, dict) or not marker_sets:
        raise ValueError("marker_sets must contain all planned program names")
    for name, panel in marker_sets.items():
        if not isinstance(name, str) or not name or not isinstance(panel, (list, tuple)) or not panel:
            raise ValueError("each named marker panel must be a nonempty list")
        _strings(panel, f"markers for {name}")
        if len(set(panel)) != len(panel):
            raise ValueError("marker panels must not contain duplicate genes")
    for name, value in (("min_cells", min_cells), ("min_markers", min_markers), ("min_control_nonzero", min_control_nonzero)):
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < 1:
            raise ValueError(f"{name} must be a positive integer")
    if not np.isfinite(scale_floor) or scale_floor <= 0 or not 0 < alpha < 1 or not 0.5 <= auc_floor <= 1:
        raise ValueError("invalid score scale, alpha, or AUC floor")
    if isinstance(n_identity_shuffle, bool) or not isinstance(n_identity_shuffle, (int, np.integer)) or n_identity_shuffle < 0:
        raise ValueError("n_identity_shuffle must be a nonnegative integer")
    discovery, inference = split_controls(labels, strata=blocks, seed=seed, control_label=control_label)
    if len(discovery) < 20 or len(inference) < 20:
        raise ValueError("at least 20 discovery and 20 inference controls are required")
    means = x[discovery].mean(axis=0)
    sd = x[discovery].std(axis=0, ddof=0)
    nonzero = (x[discovery] > 0).sum(axis=0)
    eligible_gene = (nonzero >= min_control_nonzero) & (sd > 1e-8)
    scales = np.maximum(sd, scale_floor)
    standardized = (x - means) / scales
    gidx = {gene: i for i, gene in enumerate(genes)}
    perturbations = [p for p in np.unique(labels) if p != control_label and "+" not in p and np.sum(labels == p) >= min_cells]
    pidx = {p: np.flatnonzero(labels == p) for p in perturbations}
    out = {
        "analysis": "heldout_control_external_marker_baseline", "seed": int(seed),
        "planned_programs": list(marker_sets), "n_perturbations": len(perturbations),
        "perturbations": perturbations, "family_size": len(marker_sets) * len(perturbations),
        "stratified": blocks is not None, "direction": "up", "alpha": alpha, "auc_floor": auc_floor,
        "control_split": {"n_discovery": len(discovery), "n_inference": len(inference),
                          "discovery_cell_ids": ids[discovery].tolist(), "inference_cell_ids": ids[inference].tolist()},
        "discovery_gene_parameters": {gene: {"mean": float(means[i]), "population_sd": float(sd[i]),
             "scale": float(scales[i]), "n_nonzero": int(nonzero[i]), "eligible": bool(eligible_gene[i])}
             for i, gene in enumerate(genes)},
        "multiplicity": "BH over all planned programs x eligible perturbations, separately for primary and target-omission sensitivity",
        "scope": "Conditional cell-level marker-score association; generic response, guide effects, and compositional effects remain possible. Not mechanism, novelty, independent replication, or SAE superiority.",
        "programs": {},
    }
    score_matrix = np.zeros((len(x), len(marker_sets)))
    omnibus_records = {}
    flat_rows = []
    for j, (name, markers) in enumerate(marker_sets.items()):
        indices = [gidx[g] for g in markers if g in gidx and eligible_gene[gidx[g]]]
        eligible_markers = [genes[i] for i in indices]
        sufficient = len(indices) >= min_markers
        score = standardized[:, indices].mean(axis=1) if sufficient else np.zeros(len(x))
        score_matrix[:, j] = score
        omnibus_records[name] = {"feature": j, "enrich_p_bonf": 0.0, "low_power": not sufficient}
        record = {"markers_requested": list(markers), "eligible_markers": eligible_markers,
                  "missing_markers": [g for g in markers if g not in gidx],
                  "ineligible_discovery_markers": [g for g in markers if g in gidx and not eligible_gene[gidx[g]]],
                  "low_power": not sufficient, "per_perturbation": []}
        for pert, rows in pidx.items():
            primary = (_compare_up(score, rows, inference, blocks, min_cells) if sufficient else
                       dict(auc=0.5, mw_p=1.0, n_cells=0, n_controls=0, status="low_power_marker_panel"))
            omitted = [i for i in indices if genes[i] != pert]
            target_removed = len(omitted) != len(indices)
            if len(omitted) < min_markers:
                sensitivity = dict(auc=0.5, mw_p=1.0, n_cells=0, n_controls=0, status="low_power_after_target_omission")
            elif not target_removed:
                sensitivity = dict(primary)
            else:
                omission_score = standardized[:, omitted].mean(axis=1)
                sensitivity = _compare_up(omission_score, rows, inference, blocks, min_cells)
            sensitivity.update(target_removed=target_removed, n_markers_remaining=len(omitted))
            row = {"pert": pert, **primary, "target_omission": sensitivity}
            record["per_perturbation"].append(row)
            flat_rows.append(row)
        out["programs"][name] = record
    if flat_rows:
        primary_q = multipletests([r["mw_p"] for r in flat_rows], method="fdr_bh")[1]
        sensitivity_q = multipletests([r["target_omission"]["mw_p"] for r in flat_rows], method="fdr_bh")[1]
        for row, q, oq in zip(flat_rows, primary_q, sensitivity_q):
            row.update(mw_q=float(q), modulates=bool(q < alpha and row["auc"] >= auc_floor and row["status"] == "tested"))
            omission = row["target_omission"]
            omission.update(mw_q=float(oq), modulates=bool(oq < alpha and omission["auc"] >= auc_floor and omission["status"] == "tested"))
            row["survives_target_omission"] = bool(row["modulates"] and omission["modulates"])
    for record in out["programs"].values():
        record["modulators"] = [r["pert"] for r in record["per_perturbation"] if r["modulates"]]
        record["target_omission_robust_modulators"] = [r["pert"] for r in record["per_perturbation"] if r["survives_target_omission"]]
    if n_identity_shuffle:
        from program_modulation_specificity import program_modulation_specificity
        keep = labels != control_label
        keep[inference] = True
        identity = program_modulation_specificity(score_matrix[keep], labels[keep], omnibus_records,
                 strata=None if blocks is None else blocks[keep], control_label=control_label,
                 seed=seed, min_cells=min_cells, n_shuffle=n_identity_shuffle, alpha=alpha)
        identity["external_marker_identity_note"] = "Match p=0 is stipulated for fixed external marker-score columns, not a measured hypergeometric enrichment p; low-power panels remain p=1."
        out["identity_heterogeneity"] = identity
    return out


def _file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--adata", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--markers", type=Path, default=Path("config/human_markers.json"))
    parser.add_argument("--plan", type=Path, default=Path("config/marker_program_baseline_plan.json"))
    parser.add_argument("--pert-col", default="perturbation")
    parser.add_argument("--control-value", default="control")
    parser.add_argument("--strata-col", default="gemgroup", help="technical exchangeability column; use 'none' for labeled pooled sensitivity")
    parser.add_argument("--duplicate-genes", choices=["reject", "sum"], default="reject")
    parser.add_argument("--n-hvg", type=int, default=0)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--n-identity-shuffle", type=int, default=999)
    args = parser.parse_args()
    if args.out.exists():
        parser.error(f"refusing to overwrite existing output: {args.out}")
    if args.n_hvg != 0:
        parser.error("the declared external-marker baseline uses all genes for normalization and n_hvg=0")
    marker_sets = json.loads(args.markers.read_text())["marker_sets"]
    if not args.plan.is_file():
        parser.error("the declared baseline plan must exist before running")
    source_paths = {"source": Path(__file__).resolve(), "plan": args.plan.resolve(),
                    "markers": args.markers.resolve(),
                    "identity_module": Path(__file__).resolve().with_name("program_modulation_specificity.py")}
    source_provenance = {name: {"path": str(path), "sha256": _file_sha256(path)}
                         for name, path in source_paths.items()}
    started_utc = datetime.now(timezone.utc).isoformat()
    import anndata as ad
    annotated = ad.read_h5ad(args.adata)
    labels = annotated.obs[args.pert_col].astype(str).to_numpy()
    labels[labels == args.control_value] = "control"
    blocks = None if args.strata_col.lower() == "none" else annotated.obs[args.strata_col].astype(str).to_numpy()
    print(json.dumps({"stage": "loaded_raw", "cells": annotated.n_obs, "genes": annotated.n_vars,
                      "strata_col": args.strata_col, "duplicate_policy": args.duplicate_genes}), flush=True)
    expression, genes, normalization = normalize_counts_to_markers(annotated.X, annotated.var_names.to_numpy(), marker_sets,
                                                                  duplicate_genes=args.duplicate_genes)
    ids = annotated.obs_names.astype(str).to_numpy()
    del annotated
    print(json.dumps({"stage": "normalized_all_genes_then_markers", **normalization}), flush=True)
    result = marker_program_baseline(expression, genes, labels, marker_sets, strata=blocks,
                                    cell_ids=ids, seed=args.seed, n_identity_shuffle=args.n_identity_shuffle)
    result["normalization"] = normalization
    for name, path in source_paths.items():
        if _file_sha256(path) != source_provenance[name]["sha256"]:
            raise RuntimeError(f"{name} changed during analysis; refusing to archive inconsistent provenance")
    result["provenance"] = {
        "input": {"path": str(args.adata.resolve()), "sha256": _file_sha256(args.adata)},
        **source_provenance,
        "started_utc": started_utc, "completed_utc": datetime.now(timezone.utc).isoformat(),
        "argv": sys.argv, "arguments": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        "versions": {"python": platform.python_version(), **{name: importlib.metadata.version(name)
                     for name in ("numpy", "scipy", "scanpy", "anndata", "statsmodels")}},
        "strata_col": args.strata_col, "duplicate_policy": args.duplicate_genes,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x") as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps({"stage": "complete", "out": str(args.out), "family_size": result["family_size"],
                      "n_perturbations": result["n_perturbations"],
                      "primary_modulators": {name: rec["modulators"] for name, rec in result["programs"].items()},
                      "target_omission_robust_modulators": {name: rec["target_omission_robust_modulators"] for name, rec in result["programs"].items()}}), flush=True)


if __name__ == "__main__":
    main()
