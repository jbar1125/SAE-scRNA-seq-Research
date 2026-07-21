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
- p is GROUNDED iff its matched feature is (a) a CONFIDENT single-target match, (b)
  SPECIFICALLY suppressed in p's knockdown cells vs control, and (c) suppressed by p MORE
  than by other perturbations (column specificity), at BH-FDR.

COLUMN SPECIFICITY (v4, 2026-07-18 amendment). Gate (3) below is added because the v3
metric tested suppression only vs CONTROL, so an essential-gene knockdown that collapses
transcription globally -- dropping a shared 'cell-health' feature that MANY knockdowns
also drop -- passed independently for every such knockdown. The first real Replogle run
was ~10x enriched for core-essential genes as a result (docs/COMPONENT2_RESULTS.md). The
gate requires p to be a robust-z LOWER-tail outlier in its activation of the matched
feature vs how every OTHER tested perturbation moves that feature -- so a shared,
non-specific stressor (suppressed by many perturbations) is rejected, while a genuine
regulator (the sole/dominant suppressor of its own program) is kept. Validated on
synthetic (tests/test_causal_grounding.py::test_column_specificity): with the gate off, 5
shared stressors falsely ground; with it on, all 5 are rejected and the true regulators
survive.

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
COLSPEC_ALPHA = 1.0         # column-specificity gate DEFAULT OFF (v5 rollback): on real
                            # Replogle it removed GATA1 (the one lineage TF) without removing
                            # the housekeeping hits, because distinct essential KDs hit
                            # distinct features. Set to e.g. 0.10 to enable. See spec v5.
EFFSIZE_ALPHA = 0.10        # effect-size-control gate (v6, M1): matched-feature suppression
                            # must exceed what perturbations of comparable TOTAL effect size
                            # achieve. Set to 1.0 to disable (raw metric). See spec v6.
MIN_PERTS_EFFSIZE = 30      # effect-matching needs a population; below this the gate is a no-op
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


def _effect_matched_pvalue(supp, eff, k):
    """Effect-size CONTROL (M1). For each perturbation, robust upper-tail p of its matched-
    feature suppression vs its k nearest neighbours in TOTAL EFFECT SIZE. Small p = it
    suppresses its feature MORE than perturbations of comparable overall transcriptional
    effect -- i.e. the suppression is NOT merely explained by having a big effect. This
    controls the confound found on real Replogle (grounding tracks perturbation effect
    magnitude: a big-effect essential-gene knockdown always drops SOME feature). A modest-
    effect specific regulator stands out against its effect-matched peers; a big-effect
    knockdown whose suppression is typical-for-its-size does not. Returns an array of p."""
    supp = np.asarray(supp, float); eff = np.asarray(eff, float)
    n = supp.size
    out = np.ones(n)
    if n < 3:
        return out
    k = int(min(max(k, 3), n - 1))
    for i in range(n):
        d = np.abs(eff - eff[i]); d[i] = np.inf
        nb = np.argsort(d)[:k]
        out[i] = _robust_tail_p(supp[i], supp[nb], side="upper")
    return out


def causal_grounding(activations, expression, pert_labels, gene_names, tested_perts,
                     auc_floor=AUC_FLOOR, fdr=FDR, match_alpha=MATCH_ALPHA,
                     colspec_alpha=COLSPEC_ALPHA, effsize_alpha=EFFSIZE_ALPHA,
                     min_active_cells=MIN_ACTIVE_CELLS, direction="down"):
    """DIRECTION (2026-07-18, TF-atlas pivot). direction='down' = knockdown/CRISPRi: a
    regulator's program is SUPPRESSED (feature activation drops); this is the original
    metric, unchanged. direction='up' = OVEREXPRESSION (e.g. the Joung TF Atlas, GSE216481):
    a TF's program is ACTIVATED (feature activation rises). For 'up' the matching target is
    the +DE signature (genes the perturbation drives UP), the Mann-Whitney is one-sided
    'greater' (feature higher in perturbed cells), the floor is auc >= 1 - auc_floor, column
    specificity uses the UPPER tail, and the effect-size control uses the activation
    magnitude. The perturbed gene t is excluded either way (its own transcript moving is
    trivial). Overexpression is a cleaner, more specific causal test than knockdown and
    sidesteps the essential-gene / global-collapse confound entirely (see the methodology
    addendum docs/METHODOLOGY_ADDENDUM_OE.md).

    A perturbation p (KD of gene t) is GROUNDED iff:
      (1) MATCH: the feature f* whose program best aligns with p's downstream signature
          (t excluded) is a confident single-target match (align outlier, p_match <
          match_alpha) -- excludes diffuse/global perturbations; AND
      (2) SUPPRESSION: f*'s activation is stochastically LOWER in p's KD cells than in
          control (one-sided Mann-Whitney; rank-based, so robust to sparse/zero-inflated
          TopK activations), with AUC <= auc_floor (a real-sized effect); AND
      (3) COLUMN SPECIFICITY: p suppresses f* MORE than other perturbations do -- p's mean
          activation of f* is a robust-z LOWER-tail outlier vs how every other tested
          perturbation moves f* (p_colspec < colspec_alpha). This is what distinguishes a
          genuine regulator (suppresses its OWN program) from a non-specific stressor (e.g.
          an essential-gene knockdown that collapses transcription globally and drops a
          shared 'cell-health' feature that MANY knockdowns also drop). Suppression vs
          control alone (2) cannot tell these apart; column specificity can. See
          docs/COMPONENT2_RESULTS.md (the essential-gene confound); AND
      (4) EFFECT-SIZE CONTROL (v6, M1, DEFAULT ON): p's matched-feature suppression is an
          upper-tail outlier vs perturbations of comparable TOTAL transcriptional effect
          size (effsize_p < effsize_alpha). On real Replogle, grounding tracked effect
          magnitude -- a big-effect knockdown always drops SOME feature -- so raw grounding
          rewarded high-effect (housekeeping) knockdowns. This gate keeps only perturbations
          that suppress their feature MORE than their effect size predicts. No-op below
          MIN_PERTS_EFFSIZE perturbations (effect-matching needs a population); AND
      (5) the Mann-Whitney p survives BH-FDR across tested perturbations.
    Column specificity (3) is DEFAULT OFF (v5); effect-size control (4) is DEFAULT ON (v6).
    Global rate-vs-chance is still established by the pipeline's label-shuffle control."""
    pert_labels = np.asarray(pert_labels)
    is_ctrl = pert_labels == CONTROL_LABEL
    g2i = {g: i for i, g in enumerate(gene_names)}

    prog = gene_association(activations, expression, is_ctrl)
    prog_n = prog / (np.linalg.norm(prog, axis=1, keepdims=True) + 1e-12)
    ctrl_act = activations[is_ctrl]
    reliable = (ctrl_act > 0).sum(0) >= min_active_cells
    mu_x_c = expression[is_ctrl].mean(0)

    # per-perturbation mean activation of every feature (background for column specificity)
    pmean = np.full((len(tested_perts), activations.shape[1]), np.nan)
    for qi, q in enumerate(tested_perts):
        qidx = np.where(pert_labels == q)[0]
        if qidx.size:
            pmean[qi] = activations[qidx].mean(0)

    up = (direction == "up")                                 # overexpression/activation mode
    alt = "greater" if up else "less"
    col_side = "upper" if up else "lower"
    rng = np.random.default_rng(0)
    rows, mw_p = [], np.ones(len(tested_perts))
    for pi, p in enumerate(tested_perts):
        idx = np.where(pert_labels == p)[0]
        if idx.size < MIN_CELLS or is_ctrl.sum() < MIN_CELLS or not reliable.any():
            rows.append({"pert": p, "matched_feature": -1, "grounded": False}); continue
        # CROSS-FIT: split A (match the feature) | B (test the effect) -> no double-dipping
        perm = rng.permutation(idx); half = idx.size // 2
        A, B = perm[:half], perm[half:]
        de = expression[A].mean(0) - mu_x_c
        t = g2i.get(p)
        if t is not None:
            de[t] = 0.0
        # match the SUPPRESSED program (down) or the ACTIVATED program (up)
        target = (de if up else -de) / (np.linalg.norm(de) + 1e-12)
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
            u, mwp = mannwhitneyu(kd, ct, alternative=alt)
            auc = float(u) / (len(kd) * len(ct))             # = P(perturbed > control)
        except ValueError:
            mwp, auc = 1.0, 0.5
        # directional effect magnitude + floor: down wants auc low (suppressed), up wants auc high
        move = (auc - 0.5) if up else (0.5 - auc)
        passes_floor = (auc >= 1.0 - auc_floor) if up else (auc <= auc_floor)
        # COLUMN SPECIFICITY: is p a directional outlier for f vs how OTHER perts move f?
        col_bg = np.delete(pmean[:, f], pi)
        p_colspec = _robust_tail_p(float(np.mean(kd)), col_bg, side=col_side)
        # EFFECT SIZE: magnitude of p's downstream transcriptional shift (perturbed gene excluded)
        full_de = expression[idx].mean(0) - mu_x_c
        if t is not None:
            full_de[t] = 0.0
        effect_size = float(np.linalg.norm(full_de))
        mw_p[pi] = mwp
        rows.append({"pert": p, "matched_feature": f, "auc": float(auc), "move": float(move),
                     "match_confidence_p": float(p_match), "mw_p": float(mwp),
                     "colspec_p": float(p_colspec), "effect_size": effect_size,
                     "passes_match": bool(p_match < match_alpha),
                     "passes_floor": bool(passes_floor),
                     "passes_colspec": bool(p_colspec < colspec_alpha)})
    # EFFECT-SIZE CONTROL (M1): among scored perturbations, is each one's matched-feature
    # suppression an outlier vs perturbations of comparable total effect size?
    scored = [pi for pi, r in enumerate(rows) if r.get("matched_feature", -1) >= 0]
    if effsize_alpha < 1.0 and len(scored) >= MIN_PERTS_EFFSIZE:
        supp = np.array([max(0.0, rows[pi]["move"]) for pi in scored])   # directional magnitude
        eff = np.array([rows[pi]["effect_size"] for pi in scored])
        ep = _effect_matched_pvalue(supp, eff, k=max(8, len(scored) // 4))
        for j, pi in enumerate(scored):
            rows[pi]["effsize_p"] = float(ep[j])
            rows[pi]["passes_effsize"] = bool(ep[j] < effsize_alpha)
    else:                                                    # gate off / too few perts to match
        for pi in scored:
            rows[pi]["effsize_p"] = float("nan")
            rows[pi]["passes_effsize"] = True
    q = multipletests(mw_p, method="fdr_bh")[1]
    grounded = np.zeros(len(tested_perts), dtype=bool)
    for pi, r in enumerate(rows):
        r["mw_q"] = float(q[pi])
        r["grounded"] = bool(r.get("passes_match", False) and r.get("passes_floor", False)
                             and r.get("passes_colspec", False) and r.get("passes_effsize", False)
                             and q[pi] < fdr)
        grounded[pi] = r["grounded"]
    return {
        "n_tested": int(len(tested_perts)),
        "n_grounded": int(grounded.sum()),
        "causal_grounding_rate": float(grounded.mean()) if len(tested_perts) else float("nan"),
        "auc_floor": auc_floor, "fdr": fdr, "match_alpha": match_alpha,
        "colspec_alpha": colspec_alpha, "effsize_alpha": effsize_alpha,
        "direction": direction,
        "min_active_cells": min_active_cells, "n_reliable_features": int(reliable.sum()),
        "per_perturbation": rows,
    }
