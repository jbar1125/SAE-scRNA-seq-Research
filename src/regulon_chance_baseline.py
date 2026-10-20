#!/usr/bin/env python3
"""Is the best-matching feature's regulon overlap better than CHANCE? (the v2 verdict)

The regulon-recovery run selects, per TF, the feature (of n_features) whose top-TOP_K genes most
overlap that TF's curated targets. A 0-recovered headline is only interpretable once you know
what overlap the SAME search would produce on random gene sets: the maximum over many features is
much larger than the per-feature expectation, and that is what the Bonferroni correction encodes.

This converts the run's per-TF overlaps into a plain, quantified statement by simulating the
null directly: draw n_features random top-TOP_K gene sets, record the MAX overlap with the TF's
regulon, repeat. If the observed overlap sits at or below that simulated best-of-n_features
expectation, the SAE features are NOT concentrating curated regulons -- a real negative about
feature/ground-truth alignment, not an underpowered test.

CAVEAT (do not overstate): regulon_hits is stored truncated to 10 per TF, so this is valid while
observed overlaps are < 10 (they are: 2-4 in the Norman v2 run). It also says nothing about
whether the features are biologically meaningful in general -- the same SAEs recover a clean
erythroid/hemoglobin program. The claim is specifically about single-TF curated-regulon alignment.

Usage: python3 src/regulon_chance_baseline.py causal_out/rr_v2_seed*.json [--n-genes 5791]
"""
import argparse
import json

import numpy as np


def simulate_max_overlap(n_targets, n_genes, top_k, n_features, n_sim, rng):
    """E[max overlap] when drawing n_features random top_k gene sets against n_targets targets."""
    mx = np.empty(n_sim)
    for i in range(n_sim):
        mx[i] = rng.hypergeometric(n_targets, n_genes - n_targets, top_k, size=n_features).max()
    return float(mx.mean()), float(np.percentile(mx, 95))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("json", nargs="+")
    ap.add_argument("--n-genes", type=int, default=5791, help="genes in the panel for that run")
    ap.add_argument("--n-features", type=int, default=2048, help="SAE latents searched per TF")
    ap.add_argument("--n-sim", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    at_or_below = total = 0
    for path in args.json:
        res = json.load(open(path))
        tested = [r for r in res["per_tf"] if not r.get("low_power")]
        top_k = res.get("top_k", 50)
        print("=" * 74)
        print(f"{path}: recovered {res['n_recovered']}/{res['n_tf_tested']}  "
              f"(panel {args.n_genes} genes, {args.n_features} features, top_k {top_k})")
        print(f"  {'TF':<10}{'targets':>8}{'observed':>10}{'chance max':>12}{'95th pct':>10}"
              f"   verdict")
        for r in sorted(tested, key=lambda r: -r["n_targets_in_panel"]):
            nt = r["n_targets_in_panel"]
            obs = len(r.get("regulon_hits", []))
            if obs >= 10:
                note = "(hits truncated at 10 - inconclusive)"
                print(f"  {r['tf']:<10}{nt:>8}{obs:>10}{'':>12}{'':>10}   {note}")
                continue
            mean_max, p95 = simulate_max_overlap(nt, args.n_genes, top_k, args.n_features,
                                                 args.n_sim, rng)
            verdict = "at/below chance" if obs <= mean_max else "ABOVE chance"
            at_or_below += obs <= mean_max
            total += 1
            print(f"  {r['tf']:<10}{nt:>8}{obs:>10}{mean_max:>12.2f}{p95:>10.1f}   {verdict}")

    print("=" * 74)
    print(f"SUMMARY: {at_or_below}/{total} tested TFs have best-feature regulon overlap at or "
          f"BELOW the\n         best-of-{args.n_features}-random-sets expectation.")
    if total and at_or_below == total:
        print("  => SAE features do not concentrate curated single-TF regulons above chance.")
        print("     This is a REAL negative about feature/ground-truth alignment, not a")
        print("     power failure. It does NOT imply the features are biologically empty.")


if __name__ == "__main__":
    main()
