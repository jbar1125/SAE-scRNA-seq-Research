#!/usr/bin/env python3
"""Causal-grounding metric for gene-program SAEs (the elevated project's centerpiece).

Measures what fraction of genetic perturbations an SAE dictionary is CAUSALLY GROUNDED
for, on Perturb-seq (CRISPRi) data, defeating the triviality trap
(docs/ELEVATION_PLAN.md 3b):

- The SAE is trained on CONTROL cells only; it never sees a perturbation.
- Each feature's PROGRAM is its gene-association profile = correlation of the feature's
  activation with each gene across control cells. This is defined identically for an
  expression-space SAE and an embedding-space (foundation-model) SAE, so the two arms
  are compared apples-to-apples.
- A perturbation p (knockdown of gene t) is MATCHED to the feature whose program aligns
  with p's downstream signature (the genes p moves), with t itself EXCLUDED so we match
  on the regulated module, not the target's self-drop. Matching is by program, NOT by
  response, so it is not circular with the grounding test.
- p is GROUNDED iff its matched feature is SPECIFICALLY suppressed by p's knockdown:
  the feature's activation drop under p exceeds an effect floor AND is a robust outlier
  versus how every other perturbation moves that same feature (column specificity), at
  BH-FDR. "More than other knockdowns move it" is the causal-specificity claim.

numpy/scipy only (no torch / GPU) so the logic is unit-tested off the critical path; on
RunPod the caller passes `activations` from the trained SAE and `expression` (log-norm
counts) for the same cells.

Inputs
  activations  (n_cells, n_latents)   feature activations (control + perturbation cells)
  expression   (n_cells, n_genes)     log-normalized expression for the SAME cells
  pert_labels  (n_cells,)             gene knocked down per cell; CONTROL_LABEL for controls
  gene_names   list[str] length n_genes
  tested_perts list[str]              perturbations to score (should be in gene_names)
"""
from __future__ import annotations

import numpy as np
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

CONTROL_LABEL = "control"
EFFECT_FLOOR = 0.25         # min standardized activation drop to count (pre-registered)
FDR = 0.05
MIN_CELLS = 3


def _zscore_cols(M):
    mu = M.mean(0, keepdims=True)
    sd = M.std(0, keepdims=True); sd = np.where(sd == 0, np.nan, sd)
    return (M - mu) / sd


def gene_association(activations, expression, is_ctrl):
    """program[f, gene] = Pearson corr over CONTROL cells between feature f activation
    and gene expression. Identical definition for expression- and embedding-space SAEs."""
    A = _zscore_cols(activations[is_ctrl])          # (n_ctrl, n_lat)
    X = _zscore_cols(expression[is_ctrl])           # (n_ctrl, n_genes)
    n = A.shape[0]
    prog = (np.nan_to_num(A).T @ np.nan_to_num(X)) / max(n - 1, 1)   # (n_lat, n_genes)
    return prog


def _robust_z(x, pool, side="lower"):
    """One-sided robust-z tail p of x vs a pool, using median/MAD. side='lower' returns
    P(below) (small = x unusually low); side='upper' returns P(above) (small = x
    unusually high)."""
    med = np.median(pool)
    mad = np.median(np.abs(pool - med))
    scale = 1.4826 * mad if mad > 0 else (pool.std() if pool.std() > 0 else np.nan)
    if not np.isfinite(scale) or scale == 0:
        return 1.0
    z = (x - med) / scale
    return float(norm.cdf(z) if side == "lower" else norm.sf(z))


def causal_grounding(activations, expression, pert_labels, gene_names, tested_perts,
                     effect_floor=EFFECT_FLOOR, fdr=FDR):
    pert_labels = np.asarray(pert_labels)
    is_ctrl = pert_labels == CONTROL_LABEL
    g2i = {g: i for i, g in enumerate(gene_names)}
    n_lat = activations.shape[1]

    prog = gene_association(activations, expression, is_ctrl)       # (n_lat, n_genes)
    prog_n = prog / (np.linalg.norm(prog, axis=1, keepdims=True) + 1e-12)

    ctrl_act = activations[is_ctrl]
    mu_c = ctrl_act.mean(0)
    sd_c = ctrl_act.std(0); sd_c = np.where(sd_c == 0, np.nan, sd_c)
    mu_x_c = expression[is_ctrl].mean(0)

    # standardized activation response R[p, f], matched feature, and alignment per pert
    R = np.full((len(tested_perts), n_lat), np.nan)
    matched = np.full(len(tested_perts), -1, dtype=int)
    aligns = [None] * len(tested_perts)
    for pi, p in enumerate(tested_perts):
        m = pert_labels == p
        if m.sum() < MIN_CELLS or is_ctrl.sum() < MIN_CELLS:
            continue
        R[pi] = (activations[m].mean(0) - mu_c) / sd_c
        de = expression[m].mean(0) - mu_x_c                        # DE signature
        t = g2i.get(p)
        if t is not None:
            de[t] = 0.0                                            # exclude the KD gene
        # match: feature whose program aligns with the SUPPRESSED genes (-de)
        target = -de / (np.linalg.norm(de) + 1e-12)
        pm = prog_n.copy()
        if t is not None:
            pm[:, t] = 0.0
        align = pm @ target                                       # (n_lat,)
        aligns[pi] = align
        matched[pi] = int(np.nanargmax(align))

    rows, raw_p = [], np.ones(len(tested_perts))
    for pi, p in enumerate(tested_perts):
        f = matched[pi]
        if f < 0 or np.isnan(R[pi, f]):
            rows.append({"pert": p, "matched_feature": int(f), "response": None,
                         "passes_floor": False}); continue
        resp = R[pi, f]
        col = R[:, f]; col = col[~np.isnan(col)]
        col_other = np.delete(col, np.where(np.isclose(col, resp))[0][:1]) if len(col) > 1 else col
        p_spec = _robust_z(resp, col_other, side="lower")         # suppressed vs other KDs
        # match confidence: the matched program must be a CLEAR target (excludes diffuse
        # / global perturbations whose alignment is spread across many features).
        al = aligns[pi]
        al_other = np.delete(al, f)
        p_match = _robust_z(al[f], al_other, side="upper")
        raw_p[pi] = max(p_spec, p_match)                          # need BOTH
        rows.append({"pert": p, "matched_feature": int(f), "response": float(resp),
                     "specificity_p": float(p_spec), "match_confidence_p": float(p_match),
                     "passes_floor": bool(resp <= -effect_floor)})
    q = multipletests(raw_p, method="fdr_bh")[1]
    grounded = np.zeros(len(tested_perts), dtype=bool)
    for pi in range(len(tested_perts)):
        rows[pi]["specificity_q"] = float(q[pi])
        rows[pi]["grounded"] = bool(rows[pi].get("passes_floor", False) and q[pi] < fdr)
        grounded[pi] = rows[pi]["grounded"]
    return {
        "n_tested": int(len(tested_perts)),
        "n_grounded": int(grounded.sum()),
        "causal_grounding_rate": float(grounded.mean()) if len(tested_perts) else float("nan"),
        "effect_floor": effect_floor, "fdr": fdr,
        "per_perturbation": rows,
    }
