#!/usr/bin/env python3
"""Causal head-to-head pipeline (Phase 1, RunPod GPU) — expression vs embedding SAEs.

One command per arm. Trains a TopK sparse autoencoder (the standard in the
foundation-model SAE literature) on CONTROL cells only, encodes all cells, and runs the
pre-registered causal-grounding metric (src/causal_grounding.py) on Perturb-seq data.

  Arm B (contribution): --rep expression   -> SAE on log-norm gene expression
  Arm A (reproduce field): --rep embedding --embedding scgpt_emb.npy
                                            -> SAE on foundation-model embeddings

Both arms compute grounding in gene space (gene-association + DE from the SAME
expression matrix), so the comparison is apples-to-apples. Run each; compare the
`causal_grounding_rate`. Prediction (ELEVATION_PLAN): expression >> embedding, beating
the 6.2-10% embedding ceiling.

Input: an AnnData `.h5ad` with a control (non-targeting) group and per-cell knockdown
labels in an obs column. For Replogle K562: `pertpy.data.replogle_2022_k562_essential()`.

This module imports torch lazily; the metric + loaders are import-safe without a GPU.
Verify wiring with:  python3 src/causal_pipeline.py --synthetic
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import causal_grounding as cg


# --------------------------------------------------------------------------- #
# TopK SAE (lazy torch)
# --------------------------------------------------------------------------- #
def train_topk_sae(rep, latent_dim, k, epochs=200, batch=1024, lr=4e-4, seed=0, device=None):
    """Train a TopK SAE on `rep` (n_cells, d). Returns (encode_fn, info). encode_fn maps
    any (n, d) array to (n, latent_dim) TopK activations."""
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    torch.manual_seed(seed); np.random.seed(seed)
    d = rep.shape[1]
    enc = nn.Linear(d, latent_dim).to(device)
    dec = nn.Linear(latent_dim, d, bias=False).to(device)
    with torch.no_grad():                                   # unit-norm decoder atoms
        dec.weight.copy_(F.normalize(dec.weight, dim=0))
    opt = torch.optim.AdamW(list(enc.parameters()) + list(dec.parameters()), lr=lr)
    mu = rep.mean(0, keepdims=True); Xt = torch.tensor(rep - mu, dtype=torch.float32, device=device)

    def topk(z):
        v, i = torch.topk(z, k, dim=1)
        out = torch.zeros_like(z); out.scatter_(1, i, F.relu(v)); return out

    n = Xt.shape[0]
    for _ in range(epochs):
        perm = torch.randperm(n, device=device)
        for b in range(0, n, batch):
            xb = Xt[perm[b:b + batch]]
            z = topk(enc(xb)); recon = dec(z)
            loss = F.mse_loss(recon, xb)
            opt.zero_grad(); loss.backward(); opt.step()
            with torch.no_grad():
                dec.weight.copy_(F.normalize(dec.weight, dim=0))

    @torch.no_grad()
    def encode(mat):
        xb = torch.tensor(mat - mu, dtype=torch.float32, device=device)
        return topk(enc(xb)).cpu().numpy()

    return encode, {"latent_dim": latent_dim, "k": k, "d": d, "device": str(device)}


# --------------------------------------------------------------------------- #
# Perturb-seq loading
# --------------------------------------------------------------------------- #
def load_perturbseq(h5ad, pert_col, control_value, min_cells=30, n_hvg=0, tf_set=None,
                    exclude_set=None):
    """Return (expression log-norm (n,g), gene_names, pert_labels, tested_perts).

    n_hvg>0 keeps the top-n_hvg highly-variable genes UNION every tested-perturbation
    target gene (so no testable perturbation is lost), and subsets BEFORE densifying, to
    bound memory. A full genome-scale Perturb-seq densified is tens of GB; with n_hvg a
    few thousand it is a few GB. This is also standard practice for SAEs on Perturb-seq.

    tf_set (optional): restrict the SCORED perturbations to these genes (e.g. transcription
    factors). Causal/regulatory grounding is only meaningful for regulators; essential
    housekeeping-gene knockdowns cause broad non-specific effects and are not testable for
    target-specific responses. This also matches the foundation-model benchmark (TF panel).

    exclude_set (optional): drop these genes from the scored perturbations (e.g. Hart CEGv2
    core-essential genes -- knocking down a core-essential gene collapses transcription
    globally, which grounds a broad 'cell-health' feature non-specifically; the first real
    run was ~10x enriched for such genes, see docs/COMPONENT2_RESULTS.md)."""
    import anndata as ad
    import scanpy as sc
    from collections import Counter
    A = ad.read_h5ad(h5ad)
    sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    genes_all = list(map(str, A.var_names))
    raw = A.obs[pert_col].astype(str).values
    labels = np.where(raw == control_value, cg.CONTROL_LABEL, raw)
    counts = Counter(labels)
    present = set(genes_all)
    tested = sorted({p for p in labels if p != cg.CONTROL_LABEL
                     and p in present and counts[p] >= min_cells
                     and (tf_set is None or p in tf_set)
                     and (exclude_set is None or p not in exclude_set)})
    if n_hvg and n_hvg < A.shape[1]:
        sc.pp.highly_variable_genes(A, n_top_genes=n_hvg)
        keep = A.var["highly_variable"].values.copy()
        g2i = {g: i for i, g in enumerate(genes_all)}
        for p in tested:                                   # never drop a testable target
            keep[g2i[p]] = True
        A = A[:, keep].copy()                              # subset while still sparse
        genes = list(map(str, A.var_names))
    else:
        genes = genes_all
    X = np.asarray(A.X.todense() if hasattr(A.X, "todense") else A.X, dtype=np.float32)
    return X, genes, labels, tested


# --------------------------------------------------------------------------- #
# Synthetic wiring check (CPU, no data)
# --------------------------------------------------------------------------- #
def _synthetic(n_prog=8, prog_size=15, n_bg=250, rng=None):
    """Realistic-scale synthetic Perturb-seq: n_prog programs each with a regulator whose
    KD suppresses it, plus n_bg no-effect background perturbations to power the
    specificity null. Validated: the full pipeline grounds ~7/8 true regulators with
    ~0 false positives (see LAB_NOTEBOOK)."""
    rng = rng or np.random.default_rng(0)
    G = n_prog * prog_size
    genes = [f"g{i}" for i in range(G)]

    def cells(activity):
        X = activity[:, np.repeat(np.arange(n_prog), prog_size)] + rng.normal(0, .1, (activity.shape[0], G))
        return np.maximum(X, 0)

    blocks, labels = [], []
    blocks.append(cells(rng.uniform(.8, 1.2, (2000, n_prog)))); labels += ["control"] * 2000
    for f in range(n_prog):
        act = rng.uniform(.8, 1.2, (80, n_prog)); act[:, f] = rng.uniform(0, .1, 80)
        blocks.append(cells(act)); labels += [f"g{f * prog_size}"] * 80
    for j in range(n_bg):
        blocks.append(cells(rng.uniform(.8, 1.2, (40, n_prog)))); labels += [f"bg{j}"] * 40
    X = np.vstack(blocks).astype(np.float32)
    tested = [f"g{f * prog_size}" for f in range(n_prog)] + [f"bg{j}" for j in range(n_bg)]
    return X, genes, np.array(labels), tested


def shuffle_control(acts, expression, labels, genes, tested, n_shuffle, seed,
                    colspec_alpha=cg.COLSPEC_ALPHA):
    """Permute perturbation labels among perturbed cells (controls fixed) and re-run
    grounding, n_shuffle times. If the real grounding is causal signal, the shuffled
    grounding rate collapses toward 0. Returns the null rates + an empirical p."""
    labels = np.asarray(labels)
    is_pert = labels != cg.CONTROL_LABEL
    rng = np.random.default_rng(seed)
    null = []
    for _ in range(n_shuffle):
        lab = labels.copy()
        lab[is_pert] = rng.permutation(lab[is_pert])
        null.append(cg.causal_grounding(acts, expression, lab, genes, tested,
                                        colspec_alpha=colspec_alpha)["causal_grounding_rate"])
    return [float(x) for x in null]


def _tag_breakdown(res, tag_set, tag_name):
    """Split the grounding rate by whether each tested perturbation is in tag_set (e.g.
    core-essential vs not). Establishes whether grounding is driven by the tagged class --
    the essential-gene confound test -- in a single run, without excluding anything."""
    rows = res["per_perturbation"]
    inn = [r for r in rows if r["pert"] in tag_set]
    out = [r for r in rows if r["pert"] not in tag_set]
    def rate(rs):
        g = sum(bool(r.get("grounded")) for r in rs)
        return {"n": len(rs), "grounded": g, "rate": (g / len(rs) if rs else float("nan"))}
    bd = {"tag": tag_name, "in_tag": rate(inn), "not_in_tag": rate(out)}
    res.setdefault("breakdowns", []).append(bd)
    print(f"  breakdown[{tag_name}]: in-tag {bd['in_tag']['grounded']}/{bd['in_tag']['n']} "
          f"= {bd['in_tag']['rate']:.3f} | not-in-tag {bd['not_in_tag']['grounded']}/"
          f"{bd['not_in_tag']['n']} = {bd['not_in_tag']['rate']:.3f}")
    return bd


def run(rep_matrix, expression, genes, labels, tested, latent, k, seed, out, n_shuffle=0,
        auc_floor=cg.AUC_FLOOR, min_active_cells=cg.MIN_ACTIVE_CELLS,
        colspec_alpha=cg.COLSPEC_ALPHA, tag_set=None, tag_name="tagged"):
    encode, info = train_topk_sae(rep_matrix, latent, k, seed=seed)
    acts = encode(rep_matrix)
    res = cg.causal_grounding(acts, expression, labels, genes, tested,
                              auc_floor=auc_floor, min_active_cells=min_active_cells,
                              colspec_alpha=colspec_alpha)
    res["sae"] = info
    if tag_set is not None:
        _tag_breakdown(res, tag_set, tag_name)
    # diagnostic: where in the pipeline do perturbations pass/fail?
    pp = [r for r in res["per_perturbation"] if r.get("matched_feature", -1) >= 0]
    if pp:
        aucs = np.array([r.get("auc", 0.5) for r in pp])
        qs = np.array([r.get("mw_q", 1.0) for r in pp])
        diag = {"n_scored": len(pp),
                "n_pass_match": int(sum(r.get("passes_match", False) for r in pp)),
                "n_pass_floor": int(sum(r.get("passes_floor", False) for r in pp)),
                "median_auc": float(np.median(aucs)), "min_auc": float(aucs.min()),
                "n_mw_q_below_0.05": int((qs < 0.05).sum())}
        res["diagnostic"] = diag
        print(f"  diag: scored {diag['n_scored']} | pass match: {diag['n_pass_match']} "
              f"| pass floor(auc<={res['auc_floor']}): {diag['n_pass_floor']} "
              f"| median auc {diag['median_auc']:.3f} min {diag['min_auc']:.3f} "
              f"| mw_q<0.05: {diag['n_mw_q_below_0.05']}")
    if n_shuffle:
        null = shuffle_control(acts, expression, labels, genes, tested, n_shuffle, seed,
                               colspec_alpha=colspec_alpha)
        real = res["causal_grounding_rate"]
        res["shuffle_control"] = {
            "n_shuffle": n_shuffle, "null_rates_mean": float(np.mean(null)),
            "null_rates_max": float(np.max(null)),
            "empirical_p": float((1 + sum(x >= real for x in null)) / (1 + n_shuffle))}
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    json.dump(res, open(out, "w"), indent=2)
    msg = (f"[{out}] grounded {res['n_grounded']}/{res['n_tested']} "
           f"= {res['causal_grounding_rate']:.3f}")
    if n_shuffle:
        sc = res["shuffle_control"]
        msg += f"  | shuffled null mean {sc['null_rates_mean']:.3f} p={sc['empirical_p']:.3f}"
    print(msg + f"  (SAE {info})")
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--synthetic", action="store_true", help="CPU wiring check, no data")
    ap.add_argument("--adata"); ap.add_argument("--pert-col", default="gene")
    ap.add_argument("--control-value", default="non-targeting")
    ap.add_argument("--rep", choices=["expression", "embedding"], default="expression")
    ap.add_argument("--embedding", help="path to .npy (n_cells, d) foundation-model embedding")
    ap.add_argument("--latent", type=int, default=1024)
    ap.add_argument("--k", type=int, default=32)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="causal_out/grounding.json")
    ap.add_argument("--n-shuffle", type=int, default=0, help="permutation control replicates")
    ap.add_argument("--n-hvg", type=int, default=0, help="keep top-N HVGs + all perturbed genes (bounds memory)")
    ap.add_argument("--tf-list", default=None, help="restrict scored perturbations to genes in this file (e.g. transcription factors)")
    ap.add_argument("--exclude-list", default=None, help="drop genes in this file from scored perturbations (e.g. core-essential genes)")
    ap.add_argument("--tag-list", default=None, help="report grounding rate split by membership in this file (e.g. essential vs not), without excluding")
    ap.add_argument("--auc-floor", type=float, default=cg.AUC_FLOOR, help="matched feature KD-vs-ctrl AUC must be <= this")
    ap.add_argument("--colspec-alpha", type=float, default=cg.COLSPEC_ALPHA, help="column-specificity gate alpha (1.0 disables it)")
    ap.add_argument("--min-active-cells", type=int, default=cg.MIN_ACTIVE_CELLS, help="min control cells a feature must fire in to be matchable")
    args = ap.parse_args()

    if args.synthetic:
        X, genes, labels, tested = _synthetic()
        n_prog = sum(1 for t in tested if not t.startswith("bg"))
        res = run(X, X, genes, labels, tested, latent=128, k=8, seed=args.seed,
                  out=args.out, n_shuffle=20)
        pp = res["per_perturbation"]
        true_g = sum(pp[i]["grounded"] for i in range(n_prog))
        false_g = res["n_grounded"] - true_g
        sc = res["shuffle_control"]
        print(f"  self-test: {true_g}/{n_prog} true regulators grounded, "
              f"{false_g} false positives; shuffled null mean {sc['null_rates_mean']:.3f}")
        assert true_g >= n_prog - 2, f"pipeline should recover most regulators, got {true_g}/{n_prog}"
        assert false_g <= max(3, int(0.05 * (len(tested) - n_prog))), f"too many false positives: {false_g}"
        assert sc["null_rates_max"] < res["causal_grounding_rate"], "shuffle control must collapse below the real rate"
        print("  PIPELINE SELF-TEST PASSED")
        return

    def _load_set(path):
        return {ln.strip() for ln in open(path) if ln.strip()} if path else None
    tf_set = _load_set(args.tf_list)
    exclude_set = _load_set(args.exclude_list)
    tag_set = _load_set(args.tag_list)
    X, genes, labels, tested = load_perturbseq(args.adata, args.pert_col, args.control_value,
                                               n_hvg=args.n_hvg, tf_set=tf_set, exclude_set=exclude_set)
    print(f"loaded {X.shape[0]} cells x {X.shape[1]} genes; {len(tested)} testable perturbations"
          + (f" (excluded {len(exclude_set)}-gene list)" if exclude_set else ""))
    if args.rep == "expression":
        rep = X
    else:
        rep = np.load(args.embedding).astype(np.float32)
        assert rep.shape[0] == X.shape[0], "embedding rows must match cells"
    run(rep, X, genes, labels, tested, args.latent, args.k, args.seed, args.out, n_shuffle=args.n_shuffle,
        auc_floor=args.auc_floor, min_active_cells=args.min_active_cells, colspec_alpha=args.colspec_alpha,
        tag_set=tag_set, tag_name=(Path(args.tag_list).stem if args.tag_list else "tagged"))


if __name__ == "__main__":
    main()
