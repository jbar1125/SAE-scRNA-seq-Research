#!/usr/bin/env python3
"""Regulon-recovery causal test (Component 2, the promiscuity-robust redesign).

WHY THIS EXISTS. The binary causal-grounding RATE does not beat a label-shuffle null on
Norman CRISPRa (LAB_NOTEBOOK 2026-08-10): K562's variance is dominated by ~2 axes (erythroid,
stress), so the promiscuous "match a feature to each perturbation" step grounds ~any labeling.
This metric removes that failure mode by selecting each TF's feature from EXTERNAL ground truth
(TRRUST curated TF->target regulons) instead of from the data-driven matcher, then testing a
single causal question per TF.

PRE-REGISTERED TEST (frozen in config/regulon_recovery_spec.json before any real run):
  For each overexpressed TF with a curated TRRUST regulon (>= MIN_TARGETS targets present in the
  gene panel, self-gene excluded):
    1. FEATURE SELECTION (external, label-free): pick the SAE feature whose top-TOP_K
       associated genes (Pearson gene-association over CONTROL cells) are most enriched for the
       TF's TRRUST targets (min hypergeometric p; Bonferroni-corrected over all features). This
       is "the feature carrying TF X's known regulon". It does NOT use the OE-cell labels, so it
       is invariant under the label shuffle.
    2. CAUSAL ACTIVATION: is that feature's activation higher in TF-OE cells than in control
       cells? Mann-Whitney one-sided greater -> AUC = P(act_OE > act_ctrl); require AUC >= FLOOR.
    3. A TF RECOVERS iff enrichment q < ALPHA (BH over TFs) AND activation q < ALPHA (BH over
       TFs) AND AUC >= FLOOR.
  HEADLINE: n_recovered (real) vs its LABEL-SHUFFLE null (permute the OE labels among non-control
  cells; feature selection stays fixed, only the causal activation is re-tested). Empirical p =
  (1 + #{shuffle >= real}) / (1 + n_shuffle). A positive result is n_recovered above the null.

HONESTY. TRRUST covers only well-studied TFs, and some (e.g. KLF1: 2 targets) are underpowered;
report those as low-power, not as negatives. A non-significant headline is a real negative, not
an artifact of a blind method -- feature selection here is NOT blind (it uses external truth),
which is the whole point of the redesign. Never freeze a recovered table from an inconclusive run.

The analysis core (regulon_recovery) is pure numpy/scipy and is unit-tested against a synthetic
oracle (tests/test_regulon_recovery.py); the CLI trains the SAE (torch) like causal_pipeline or
loads a cached activation matrix.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import hypergeom, mannwhitneyu
from statsmodels.stats.multitest import multipletests

import causal_grounding as cg

TOP_K = 50            # genes per feature considered its "program" for enrichment
MIN_TARGETS = 5       # min TRRUST targets present in panel to test a TF
AUC_FLOOR = 0.55      # activation floor: feature must be up in OE vs control (up-direction)
ALPHA = 0.05          # BH-FDR threshold, both families
MIN_CELLS = 20        # min OE cells to test a TF


def _mw_auc_greater(a, b):
    """AUC = P(a > b) and one-sided Mann-Whitney p (H1: a stochastically greater than b)."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if len(a) == 0 or len(b) == 0:
        return 0.5, 1.0
    # degenerate: no variation across the pooled sample -> undefined, treat as null
    if np.ptp(np.concatenate([a, b])) == 0:
        return 0.5, 1.0
    try:
        u, p = mannwhitneyu(a, b, alternative="greater")
    except ValueError:
        return 0.5, 1.0
    return float(u) / (len(a) * len(b)), float(p)


def _select_feature(prog, targets_idx, n_genes, top_k):
    """Feature whose top-`top_k` associated genes are most enriched for `targets_idx`.

    Returns (feature index, raw hypergeometric p, Bonferroni p over all features)."""
    tset = set(int(i) for i in targets_idx)
    n_draw = len(tset)
    best_f, best_p = -1, 1.0
    F = prog.shape[0]
    for f in range(F):
        top = np.argsort(prog[f])[::-1][:top_k]
        hits = sum(1 for i in top if int(i) in tset)
        p = float(hypergeom.sf(hits - 1, n_genes, n_draw, top_k)) if hits > 0 else 1.0
        if p < best_p:
            best_p, best_f = p, f
    return best_f, best_p, min(1.0, best_p * F)


def regulon_recovery(activations, expression, labels, gene_names, regulons,
                     top_k=TOP_K, min_targets=MIN_TARGETS, auc_floor=AUC_FLOOR,
                     alpha=ALPHA, min_cells=MIN_CELLS, n_shuffle=0, seed=0):
    """See module docstring. Returns a dict with per-TF rows and the shuffle-null headline."""
    activations = np.asarray(activations, float)
    labels = np.asarray(labels).astype(str)
    gname = list(map(str, gene_names))
    gidx = {g: i for i, g in enumerate(gname)}
    n_genes = len(gname)
    is_ctrl = labels == cg.CONTROL_LABEL
    ctrl_idx = np.where(is_ctrl)[0]
    prog = cg.gene_association(activations, expression, is_ctrl)   # [F, G] over control cells

    rows = []
    for tf in sorted(set(labels)):
        if tf == cg.CONTROL_LABEL or "+" in tf or tf not in regulons:
            continue
        targets = sorted({gidx[g] for g in regulons[tf] if g in gidx and g != tf})
        oe_idx = np.where(labels == tf)[0]
        row = {"tf": tf, "n_targets_in_panel": len(targets), "n_oe_cells": int(len(oe_idx))}
        if len(targets) < min_targets or len(oe_idx) < min_cells:
            row.update(low_power=True, feature=-1, recovers=False)
            rows.append(row)
            continue
        f, enrich_p, enrich_p_bonf = _select_feature(prog, targets, n_genes, top_k)
        auc, mw_p = _mw_auc_greater(activations[oe_idx, f], activations[ctrl_idx, f])
        top_genes = [gname[int(i)] for i in np.argsort(prog[f])[::-1][:10]]
        hit_targets = [gname[int(i)] for i in np.argsort(prog[f])[::-1][:top_k]
                       if int(i) in set(targets)]
        row.update(low_power=False, feature=int(f), enrich_p=enrich_p,
                   enrich_p_bonf=enrich_p_bonf, auc=auc, mw_p=mw_p, top_genes=top_genes,
                   regulon_hits=hit_targets[:10], passes_floor=bool(auc >= auc_floor))
        rows.append(row)

    tested = [r for r in rows if not r["low_power"]]

    def _mark_recovers(subset, mw_key="mw_p"):
        if not subset:
            return 0
        eq = multipletests([r["enrich_p_bonf"] for r in subset], method="fdr_bh")[1]
        aq = multipletests([r[mw_key] for r in subset], method="fdr_bh")[1]
        n = 0
        for r, e, a in zip(subset, eq, aq):
            r["enrich_q"] = float(e)
            r[mw_key.replace("_p", "_q")] = float(a)
            r["recovers"] = bool(e < alpha and a < alpha and r["passes_floor"])
            n += int(r["recovers"])
        return n

    n_recovered = _mark_recovers(tested)
    for r in rows:
        r.setdefault("recovers", False)

    out = {"n_tf_tested": len(tested), "n_recovered": int(n_recovered),
           "recovery_rate": (n_recovered / len(tested) if tested else float("nan")),
           "top_k": top_k, "auc_floor": auc_floor, "alpha": alpha,
           "recovered_tfs": [r["tf"] for r in tested if r["recovers"]],
           "per_tf": rows}

    if n_shuffle and tested:
        rng = np.random.default_rng(seed)
        nonctrl = np.where(~is_ctrl)[0]
        feats = {r["tf"]: r["feature"] for r in tested}
        enrich_bonf = {r["tf"]: r["enrich_p_bonf"] for r in tested}
        null_counts = []
        for _ in range(n_shuffle):
            perm = labels.copy()
            perm[nonctrl] = rng.permutation(labels[nonctrl])
            shuf = []
            for r in tested:
                tf = r["tf"]
                f = feats[tf]
                idx = np.where(perm == tf)[0]
                if len(idx) < min_cells:
                    shuf.append({"enrich_p_bonf": enrich_bonf[tf], "mw_p": 1.0,
                                 "passes_floor": False})
                    continue
                auc, mw_p = _mw_auc_greater(activations[idx, f], activations[ctrl_idx, f])
                shuf.append({"enrich_p_bonf": enrich_bonf[tf], "mw_p": mw_p,
                             "passes_floor": bool(auc >= auc_floor)})
            null_counts.append(_mark_recovers(shuf))
        null_counts = np.array(null_counts)
        out["shuffle_null"] = {
            "n_shuffle": n_shuffle, "null_mean": float(null_counts.mean()),
            "null_max": int(null_counts.max()),
            "empirical_p": float((1 + int((null_counts >= n_recovered).sum())) / (1 + n_shuffle))}
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
    ap.add_argument("--trrust", default="config/trrust/trrust_rawdata.human.tsv")
    ap.add_argument("--activations", help="cached (n_cells, latent) .npy to skip SAE training")
    ap.add_argument("--latent", type=int, default=2048)
    ap.add_argument("--k", type=int, default=32)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--n-hvg", type=int, default=4000)
    ap.add_argument("--n-shuffle", type=int, default=50)
    ap.add_argument("--top-k", type=int, default=TOP_K)
    ap.add_argument("--min-targets", type=int, default=MIN_TARGETS)
    ap.add_argument("--out", default="causal_out/regulon_recovery.json")
    args = ap.parse_args()

    from trrust_crossref import load_trrust
    from causal_pipeline import load_perturbseq, train_topk_sae
    regulons = load_trrust(args.trrust)
    X, genes, labels, tested = load_perturbseq(args.adata, args.pert_col, args.control_value,
                                               n_hvg=args.n_hvg)
    print(f"loaded {X.shape[0]} cells x {X.shape[1]} genes; {len(tested)} single-gene perts; "
          f"{len(regulons)} TRRUST regulons")
    if args.activations:
        acts = np.load(args.activations)
    else:
        encode, info = train_topk_sae(X, args.latent, args.k, seed=args.seed)
        acts = encode(X)
        print(f"trained SAE {info}")
    res = regulon_recovery(acts, X, labels, genes, regulons, top_k=args.top_k,
                           min_targets=args.min_targets, n_shuffle=args.n_shuffle, seed=args.seed)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    json.dump(res, open(args.out, "w"), indent=2)
    print(f"[{args.out}] recovered {res['n_recovered']}/{res['n_tf_tested']} TFs "
          f"= {res['recovery_rate']:.3f}"
          + (f"  | shuffle null mean {res['shuffle_null']['null_mean']:.2f} "
             f"p={res['shuffle_null']['empirical_p']:.3f}" if "shuffle_null" in res else ""))
    for r in sorted((r for r in res["per_tf"] if r.get("recovers")),
                    key=lambda r: r.get("enrich_q", 1.0)):
        print(f"    {r['tf']:<8} feat {r['feature']} auc {r['auc']:.3f} "
              f"enrichQ {r.get('enrich_q', 1):.2g} actQ {r.get('mw_q', 1):.2g} "
              f"| regulon hits: {', '.join(r['regulon_hits'][:6])}")


if __name__ == "__main__":
    main()
