#!/usr/bin/env python3
"""Causal-grounding metric for gene-program SAEs (the elevated project's centerpiece).

Measures what fraction of genetic perturbations an SAE dictionary is CAUSALLY GROUNDED
for, on Perturb-seq (CRISPRi) data, defeating the triviality trap
(docs/ELEVATION_PLAN.md 3b):

- The SAE is trained on CONTROL cells only; it never sees a perturbation.
- Each feature's PROGRAM is its gene-association profile = correlation of the feature's
  activation with each gene across control cells. Defined identically for expression-
  and embedding-space SAEs, so the two arms are apples-to-apples.
- A perturbation p (knockdown of gene t) is MATCHED to the feature whose program aligns
  with p's downstream signature (the genes p moves), with t EXCLUDED so we match on the
  regulated module, not the target's self-drop. Matching is by program, not by response.
- p is GROUNDED iff its matched feature is (a) a CONFIDENT single-target match and (b)
  SPECIFICALLY suppressed in p's knockdown cells vs control, at BH-FDR.

CROSS-FIT + RANK TEST (v3, 2026-07-18 amendment; supersedes the v2 relative-change note).
Two changes fix two real bugs caught before any valid Replogle result:
  1. SD-standardization was invalid for sparse TopK activations. A feature fires in only
     ~k/L of cells, so its activation SD is large relative to its mean and even a full
     shutoff is only ~0.1 SD -- unreachable by any effect floor, forcing 0 grounded. The
     suppression test is now the one-sided Mann-Whitney U of the matched feature's
     activation (KD cells vs control cells), reported as AUC = U / (n_kd * n_ctrl); AUC
     < 0.5 means suppressed, and the floor is on AUC (<= auc_floor), a rank statistic
     with no scale assumption -- robust to sparse/zero-inflated activations.
  2. CROSS-FITTING removes double-dipping. Each perturbation's cells are split A|B. The
     feature is MATCHED on split A (which genes p moves) and the suppression Mann-Whitney
     is tested on the held-out split B, so the same cells never both select and validate
     the feature. Without this, a diffuse/background perturbation could be "grounded" by
     selecting whichever feature happened to dip in the very cells being tested.
Only features active in >= min_active_cells control cells are matchable (near-dead
features have unreliable programs). Global "is the rate above chance" is established by
the pipeline's label-shuffle null, not per-perturbation. See PREREGISTRATION amendment.

numpy/scipy only (no torch / GPU); unit-tested off the critical path.

Inputs
  activations  (n_cells, n_latents)   feature activations (control + perturbation cells)
  expression   (n_cells, n_genes)     log-normalized expression for the SAME cells
  pert_labels  (n_cells,)             gene knocked down per cell; CONTROL_LABEL for controls
  gene_names   list[str] length n_genes
  tested_perts list[str]              perturbations to score (should be in gene_names)
"""
from __future__ import annotations

import numpy as np
from scipy.stats import norm, mannwhitneyu
from statsmodels.stats.multitest import multipletests

CONTROL_LABEL = "control"
AUC_FLOOR = 0.45            # matched feature's KD-vs-ctrl AUC must be <= this (suppressed)
FDR = 0.05
MATCH_ALPHA = 0.10          # match-confidence gate: matched alignment is a clear outlier
MIN_CELLS = 20             # min cells per perturbation group to score at all
MIN_ACTIVE_CELLS = 50       # a feature must fire in >= this many control cells to be matchable


def _zscore_cols(M):
    mu = M.mean(0, keepdims=True)
    sd = M.std(0, keepdims=True); sd = np.where(sd == 0, np.nan, sd)
    return (M - mu) / sd


def gene_association(activations, expression, is_ctrl):
    """program[f, gene] = Pearson corr over CONTROL cells between feature f activation
    and gene expression. Identical for expression- and embedding-space SAEs."""
    A = _zscore_cols(activations[is_ctrl])
    X = _zscore_cols(expression[is_ctrl])
    n = A.shape[0]
    return (np.nan_to_num(A).T @ np.nan_to_num(X)) / max(n - 1, 1)   # (n_lat, n_genes)


def _robust_tail_p(x, pool, side="lower"):
    """One-sided robust-z tail p of x vs pool (median/MAD). side='lower': small = x
    unusually low; 'upper': small = x unusually high."""
    pool = pool[np.isfinite(pool)]
    if pool.size < 2:
        return 1.0
    med = np.median(pool)
    mad = np.median(np.abs(pool - med))
    scale = 1.4826 * mad if mad > 0 else (pool.std() if pool.std() > 0 else np.nan)
    if not np.isfinite(scale) or scale == 0:
        return 1.0
    z = (x - med) / scale
    return float(norm.cdf(z) if side == "lower" else norm.sf(z))


def causal_grounding(activations, expression, pert_labels, gene_names, tested_perts,
                     auc_floor=AUC_FLOOR, fdr=FDR, match_alpha=MATCH_ALPHA,
                     min_active_cells=MIN_ACTIVE_CELLS):
    """A perturbation p (KD of gene t) is GROUNDED iff:
      (1) MATCH: the feature f* whose program best aligns with p's downstream signature
          (t excluded) is a confident single-target match (align outlier, p_match <
          match_alpha) -- excludes diffuse/global perturbations; AND
      (2) SUPPRESSION: f*'s activation is stochastically LOWER in p's KD cells than in
          control (one-sided Mann-Whitney; rank-based, so robust to sparse/zero-inflated
          TopK activations), with AUC <= auc_floor (a real-sized effect); AND
      (3) the Mann-Whitney p survives BH-FDR across tested perturbations.
    Global specificity ("is the rate above chance") is established by the pipeline's
    label-shuffle control, not per-perturbation, which keeps this test simple and robust."""
    pert_labels = np.asarray(pert_labels)
    is_ctrl = pert_labels == CONTROL_LABEL
    g2i = {g: i for i, g in enumerate(gene_names)}

    prog = gene_association(activations, expression, is_ctrl)
    prog_n = prog / (np.linalg.norm(prog, axis=1, keepdims=True) + 1e-12)
    ctrl_act = activations[is_ctrl]
    reliable = (ctrl_act > 0).sum(0) >= min_active_cells
    mu_x_c = expression[is_ctrl].mean(0)

    rng = np.random.default_rng(0)
    rows, mw_p = [], np.ones(len(tested_perts))
    for pi, p in enumerate(tested_perts):
        idx = np.where(pert_labels == p)[0]
        if idx.size < MIN_CELLS or is_ctrl.sum() < MIN_CELLS or not reliable.any():
            rows.append({"pert": p, "matched_feature": -1, "grounded": False}); continue
        # CROSS-FIT: split A (match the feature) | B (test suppression) -> no double-dipping
        perm = rng.permutation(idx); half = idx.size // 2
        A, B = perm[:half], perm[half:]
        de = expression[A].mean(0) - mu_x_c
        t = g2i.get(p)
        if t is not None:
            de[t] = 0.0
        target = -de / (np.linalg.norm(de) + 1e-12)
        pm = prog_n.copy()
        if t is not None:
            pm[:, t] = 0.0
        align = pm @ target
        align[~reliable] = -np.inf
        f = int(np.argmax(align))                            # matched on split A
        al_rel = align[reliable]
        p_match = _robust_tail_p(align[f], al_rel[al_rel != align[f]] if (al_rel != align[f]).any()
                                 else al_rel, side="upper")
        kd = activations[B, f]; ct = ctrl_act[:, f]          # tested on held-out split B
        try:
            u, mwp = mannwhitneyu(kd, ct, alternative="less")
            auc = float(u) / (len(kd) * len(ct))             # <0.5 = suppressed
        except ValueError:
            mwp, auc = 1.0, 0.5
        mw_p[pi] = mwp
        rows.append({"pert": p, "matched_feature": f, "auc": float(auc),
                     "match_confidence_p": float(p_match), "mw_p": float(mwp),
                     "passes_match": bool(p_match < match_alpha),
                     "passes_floor": bool(auc <= auc_floor)})
    q = multipletests(mw_p, method="fdr_bh")[1]
    grounded = np.zeros(len(tested_perts), dtype=bool)
    for pi, r in enumerate(rows):
        r["mw_q"] = float(q[pi])
        r["grounded"] = bool(r.get("passes_match", False) and r.get("passes_floor", False)
                             and q[pi] < fdr)
        grounded[pi] = r["grounded"]
    return {
        "n_tested": int(len(tested_perts)),
        "n_grounded": int(grounded.sum()),
        "causal_grounding_rate": float(grounded.mean()) if len(tested_perts) else float("nan"),
        "auc_floor": auc_floor, "fdr": fdr, "match_alpha": match_alpha,
        "min_active_cells": min_active_cells, "n_reliable_features": int(reliable.sum()),
        "per_perturbation": rows,
    }
