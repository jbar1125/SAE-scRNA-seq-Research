#!/usr/bin/env python3
"""Program-modulation test: is a NAMED program causally modulable by a specific perturbation set?

WHY THIS QUESTION, AFTER FOUR NULLS. The Component-2 causal program established that
single-perturbation resolution is dead (LAB_NOTEBOOK 2026-08-10): the grounding rate does not
beat a shuffle null (p=0.902), and SAE features do not concentrate curated single-TF regulons
above chance. But Gate 0 DID find TRRUST regulon enrichment at LINEAGE granularity (p=6.4e-5).
That brackets the method's resolution limit and predicts the next test: drop to program level.

Instead of asking "does TF X ground its own feature" (2048-way search per TF, dead), this asks:
  For ONE pre-specified, externally-defined program (e.g. erythroid), which perturbations
  causally activate it, is that set larger than chance, and is it SPECIFIC?

NON-CIRCULARITY (the trap this design exists to avoid). The erythroid feature was previously
NOTICED via which perturbations grounded it (CBL/UBASH3B/PTPN12); re-testing those on the same
data would be circular and worthless. Here:
  - the program is defined by config/human_markers.json, curated for Gate 0 BEFORE any Norman
    work, so the gene set is external to everything measured here;
  - the feature is selected over CONTROL cells only (label-free), so perturbation labels never
    touch selection;
  - WHICH perturbations pass is a RESULT, never an input. No perturbation is named in advance.

THE THREE READOUTS (all pre-registered in config/program_modulation_spec.json):
  1. n_modulators: perturbations with BH-q < ALPHA and AUC >= AUC_FLOOR on the program feature.
     Multiple testing is over ~100 PERTURBATIONS, not 2048 features -- that is why this is
     powered where regulon recovery was not.
  2. LABEL-SHUFFLE NULL (mode "full", the default): permute ALL labels including control, so
     group membership is independent of biology, and recount. Answers "is the modulator set
     larger than chance?" NOTE: the legacy mode "perturbed" (permute among perturbed cells only,
     controls fixed) INFLATES this metric -- spreading strongly activating cells across every
     group lifts every group above the fixed control, so the null exceeds the real count (the
     oracle measured null 14.4 vs real 5). It is retained only to ask the different, stricter
     question "does perturbation identity matter beyond generic perturbation response?".
  3. FEATURE-BACKGROUND NULL: run the identical per-perturbation test on n_random_features other
     reliable features. Answers the promiscuity objection head-on -- K562 is dominated by ~2
     axes, so if EVERY feature has ~as many modulators, a large count means nothing. The program
     must rank high against this background to count.

HONEST PRE-STATEMENT OF THE LIKELY NEGATIVE: if the erythroid program simply IS K562's dominant
axis, many perturbations will nudge it and the feature-background percentile will be unremarkable.
That outcome is a real negative and must be reported as one, not spun as "many modulators found".
A positive requires a modulator set that is both above the shuffle null AND high against the
feature background AND biologically coherent on inspection.

Pure numpy/scipy; unit-tested against a synthetic oracle (tests/test_program_modulation.py).
"""
from __future__ import annotations

import numpy as np
from scipy.stats import hypergeom, mannwhitneyu
from statsmodels.stats.multitest import multipletests

import causal_grounding as cg

TOP_K = 50               # genes defining a feature's program, for marker enrichment
AUC_FLOOR = 0.55         # activation floor (up-direction)
ALPHA = 0.05             # BH-FDR over perturbations
MIN_CELLS = 20           # min cells for a perturbation to be tested
MIN_ACTIVE_CELLS = 50    # min control cells a feature must fire in to be eligible
N_RANDOM_FEATURES = 100  # feature-background null size


def _mw_auc_greater(a, b):
    """AUC = P(a > b) with one-sided Mann-Whitney p (H1: a greater)."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    if len(a) == 0 or len(b) == 0:
        return 0.5, 1.0
    if np.ptp(np.concatenate([a, b])) == 0:
        return 0.5, 1.0
    try:
        u, p = mannwhitneyu(a, b, alternative="greater")
    except ValueError:
        return 0.5, 1.0
    return float(u) / (len(a) * len(b)), float(p)


def _count_modulators(acts_f, ctrl_idx, pert_idx, auc_floor, alpha):
    """# perturbations that significantly ACTIVATE feature column acts_f. Returns (n, rows)."""
    ctrl = acts_f[ctrl_idx]
    rows, ps = [], []
    for p, idx in pert_idx.items():
        auc, mw_p = _mw_auc_greater(acts_f[idx], ctrl)
        rows.append({"pert": p, "auc": auc, "mw_p": mw_p, "n_cells": int(len(idx))})
        ps.append(mw_p)
    if not rows:
        return 0, []
    q = multipletests(ps, method="fdr_bh")[1]
    n = 0
    for r, qq in zip(rows, q):
        r["mw_q"] = float(qq)
        r["modulates"] = bool(qq < alpha and r["auc"] >= auc_floor)
        n += r["modulates"]
    return n, rows


def select_program_feature(prog, marker_idx, n_genes, top_k, eligible):
    """Feature whose top-`top_k` control-cell gene-associations are most marker-enriched.

    Label-free: `prog` is computed over control cells only. Returns (feature, raw p, bonf p)."""
    mset = set(int(i) for i in marker_idx)
    best_f, best_p = -1, 1.0
    for f in np.asarray(eligible):
        f = int(f)
        top = np.argsort(prog[f])[::-1][:top_k]
        hits = sum(1 for i in top if int(i) in mset)
        p = float(hypergeom.sf(hits - 1, n_genes, len(mset), top_k)) if hits > 0 else 1.0
        if p < best_p:
            best_p, best_f = p, f
    return best_f, best_p, min(1.0, best_p * max(1, len(eligible)))


def program_modulation(activations, expression, labels, gene_names, marker_sets,
                       top_k=TOP_K, auc_floor=AUC_FLOOR, alpha=ALPHA, min_cells=MIN_CELLS,
                       min_active_cells=MIN_ACTIVE_CELLS, n_shuffle=0,
                       n_random_features=N_RANDOM_FEATURES, shuffle_mode="full", seed=0):
    """See module docstring. marker_sets: {program_name: [gene symbols]}."""
    activations = np.asarray(activations, float)
    labels = np.asarray(labels).astype(str)
    gname = list(map(str, gene_names))
    gidx = {g: i for i, g in enumerate(gname)}
    n_genes = len(gname)
    rng = np.random.default_rng(seed)

    is_ctrl = labels == cg.CONTROL_LABEL
    ctrl_idx = np.where(is_ctrl)[0]
    prog = cg.gene_association(activations, expression, is_ctrl)

    # eligible features: fire in enough CONTROL cells (unreliable features have noisy programs)
    active_ctrl = (activations[ctrl_idx] > 0).sum(axis=0)
    eligible = np.where(active_ctrl >= min_active_cells)[0]

    pert_idx = {}
    for p in sorted(set(labels)):
        if p == cg.CONTROL_LABEL or "+" in p:
            continue
        idx = np.where(labels == p)[0]
        if len(idx) >= min_cells:
            pert_idx[p] = idx

    out = {"n_perturbations": len(pert_idx), "n_eligible_features": int(len(eligible)),
           "top_k": top_k, "auc_floor": auc_floor, "alpha": alpha, "programs": {}}
    if len(eligible) == 0 or not pert_idx:
        return out

    # ---- feature-background null: identical test on random eligible features ----
    n_rand = min(n_random_features, len(eligible))
    rand_feats = rng.choice(eligible, size=n_rand, replace=False)
    bg_counts = np.array([_count_modulators(activations[:, int(f)], ctrl_idx, pert_idx,
                                            auc_floor, alpha)[0] for f in rand_feats])
    out["feature_background"] = {
        "n_features": int(n_rand), "mean": float(bg_counts.mean()),
        "median": float(np.median(bg_counts)), "max": int(bg_counts.max())}

    for name, genes in marker_sets.items():
        midx = sorted({gidx[g] for g in genes if g in gidx})
        rec = {"n_markers_in_panel": len(midx), "markers_found": [gname[i] for i in midx]}
        if len(midx) < 3:
            rec["low_power"] = True
            out["programs"][name] = rec
            continue
        f, ep, ebonf = select_program_feature(prog, midx, n_genes, top_k, eligible)
        n_mod, rows = _count_modulators(activations[:, f], ctrl_idx, pert_idx, auc_floor, alpha)
        rec.update(low_power=False, feature=int(f), enrich_p=ep, enrich_p_bonf=ebonf,
                   n_modulators=int(n_mod),
                   modulators=[r["pert"] for r in sorted(rows, key=lambda r: r["mw_q"])
                               if r["modulates"]],
                   top_program_genes=[gname[int(i)] for i in np.argsort(prog[f])[::-1][:12]],
                   feature_background_percentile=float((bg_counts < n_mod).mean() * 100),
                   per_perturbation=rows)
        if n_shuffle:
            null = []
            nonctrl = np.where(~is_ctrl)[0]
            for _ in range(n_shuffle):
                perm = labels.copy()
                if shuffle_mode == "perturbed":
                    # legacy/stringent: identity destroyed, perturbed-vs-control contrast KEPT.
                    # Asks "does perturbation identity matter beyond generic perturbation
                    # response?" -- but it INFLATES this metric, because spreading strongly
                    # activating cells across every group lifts every group above control.
                    perm[nonctrl] = rng.permutation(labels[nonctrl])
                    cidx = ctrl_idx
                else:
                    # default "full": permute ALL labels incl. control, so group membership is
                    # independent of biology. The standard permutation null for "is there any
                    # group structure at all", and the one this metric is calibrated against.
                    perm = rng.permutation(labels)
                    cidx = np.where(perm == cg.CONTROL_LABEL)[0]
                pidx = {p: np.where(perm == p)[0] for p in pert_idx}
                pidx = {p: i for p, i in pidx.items() if len(i) >= min_cells}
                null.append(_count_modulators(activations[:, f], cidx, pidx,
                                              auc_floor, alpha)[0])
            null = np.array(null)
            rec["shuffle_null"] = {
                "mode": shuffle_mode, "n_shuffle": n_shuffle, "mean": float(null.mean()), "max": int(null.max()),
                "empirical_p": float((1 + int((null >= n_mod).sum())) / (1 + n_shuffle))}
        out["programs"][name] = rec
    return out


# ------------------------------------------------------------------ CLI ------
def main():
    import argparse
    import json
    from pathlib import Path

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--adata", required=True)
    ap.add_argument("--pert-col", default="perturbation")
    ap.add_argument("--control-value", default="control")
    ap.add_argument("--markers", default="config/human_markers.json")
    ap.add_argument("--latent", type=int, default=2048)
    ap.add_argument("--k", type=int, default=32)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--n-hvg", type=int, default=4000)
    ap.add_argument("--n-shuffle", type=int, default=50)
    ap.add_argument("--n-random-features", type=int, default=N_RANDOM_FEATURES)
    ap.add_argument("--shuffle-mode", choices=["full", "perturbed"], default="full",
                    help="full=permute all labels (standard); perturbed=legacy, keeps the "
                         "perturbed-vs-control contrast and INFLATES this metric")
    ap.add_argument("--out", default="causal_out/program_modulation.json")
    args = ap.parse_args()

    from causal_pipeline import load_perturbseq, train_topk_sae
    marker_sets = json.load(open(args.markers))["marker_sets"]
    keep = {g for gs in marker_sets.values() for g in gs}
    X, genes, labels, tested = load_perturbseq(args.adata, args.pert_col, args.control_value,
                                               n_hvg=args.n_hvg, keep_genes=keep)
    print(f"loaded {X.shape[0]} cells x {X.shape[1]} genes; {len(tested)} single-gene perts; "
          f"{len(marker_sets)} named programs")
    encode, info = train_topk_sae(X, args.latent, args.k, seed=args.seed)
    acts = encode(X)
    print(f"trained SAE {info}")
    res = program_modulation(acts, X, labels, genes, marker_sets, n_shuffle=args.n_shuffle,
                             n_random_features=args.n_random_features,
                             shuffle_mode=args.shuffle_mode, seed=args.seed)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    json.dump(res, open(args.out, "w"), indent=2)

    bg = res.get("feature_background", {})
    print(f"\n[{args.out}] {res['n_perturbations']} perturbations x "
          f"{res['n_eligible_features']} eligible features")
    print(f"  FEATURE-BACKGROUND (random features): mean {bg.get('mean', float('nan')):.1f} "
          f"median {bg.get('median', float('nan')):.1f} max {bg.get('max', 0)} modulators")
    print(f"\n  {'program':<16}{'mkrs':>5}{'feat':>7}{'modul':>7}{'bg pct':>8}{'shufP':>8}  top genes")
    for name, r in res["programs"].items():
        if r.get("low_power"):
            print(f"  {name:<16}{r['n_markers_in_panel']:>5}   LOW POWER (<3 markers in panel)")
            continue
        sp = r.get("shuffle_null", {}).get("empirical_p", float("nan"))
        print(f"  {name:<16}{r['n_markers_in_panel']:>5}{r['feature']:>7}{r['n_modulators']:>7}"
              f"{r['feature_background_percentile']:>7.0f}%{sp:>8.3f}  "
              f"{', '.join(r['top_program_genes'][:6])}")
    print("\n  READ: a POSITIVE needs shufP < 0.05 AND a high background percentile AND a "
          "coherent modulator list.\n  Many modulators at an unremarkable percentile = the "
          "dominant-axis artifact, i.e. a negative.")
    for name, r in res["programs"].items():
        if not r.get("low_power") and r.get("modulators"):
            print(f"\n  {name} modulators ({r['n_modulators']}): "
                  f"{', '.join(r['modulators'][:20])}")


if __name__ == "__main__":
    main()
