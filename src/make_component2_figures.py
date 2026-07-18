#!/usr/bin/env python3
"""Deterministic Component-2 DIAGNOSTIC figures (not headline results).

Two panels visualising the two STABLE findings from the first real Arm-B run
(docs/COMPONENT2_RESULTS.md) -- both are diagnostic and survive regardless of the v4
re-run, so they are safe to plot now:
  A. the grounding rate is hyperparameter-DOMINATED (0.136-0.395 from knobs alone), so the
     absolute rate is not quotable; only the matched head-to-head is.
  B. the seed-stable grounded set is ~10x ENRICHED for core-essential genes and DEPLETED of
     nonessential -- the confound that motivated the v4 column-specificity fix.

Numbers are the recorded results (single-seed sweep + 5-seed set), hardcoded so the figure
is reproducible without the ephemeral raw JSONs. Run: python3 src/make_component2_figures.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path("docs/figures")

# --- recorded numbers (docs/COMPONENT2_RESULTS.md) ---
SEED_RATES = [0.136, 0.191, 0.259, 0.191, 0.198]        # seeds 0-4, latent2048/k32
SWEEP = [("512/32", 0.395), ("1024/32", 0.284), ("2048/32", 0.136),
         ("2048/16", 0.191), ("2048/64", 0.358)]         # single-seed sweep
# essential-gene analysis (Hart CEGv2 / NEGv1) on the 26 seed-stable hits
ENRICH = {"core-essential": (0.308, 0.030), "ref-nonessential": (0.0, 0.053)}  # (grounded, pool)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.2))

    # Panel A: hyperparameter sensitivity
    labels = [s[0] for s in SWEEP]; vals = [s[1] for s in SWEEP]
    axA.bar(range(len(vals)), vals, color="#4C78A8")
    mean = sum(SEED_RATES) / len(SEED_RATES)
    axA.axhline(mean, ls="--", color="#333",
                label=f"5-seed mean @2048/32 = {mean:.3f}")
    axA.axhspan(min(SEED_RATES), max(SEED_RATES), color="#4C78A8", alpha=0.12,
                label=f"seed range {min(SEED_RATES):.3f}-{max(SEED_RATES):.3f}")
    axA.set_xticks(range(len(labels))); axA.set_xticklabels(labels)
    axA.set_xlabel("SAE latent / k"); axA.set_ylabel("causal-grounding rate")
    axA.set_title("A. Rate is hyperparameter-dominated (0.14-0.40)")
    axA.legend(fontsize=8, loc="upper right")

    # Panel B: essential-gene enrichment of the grounded set
    cats = list(ENRICH); x = range(len(cats)); w = 0.38
    grounded = [ENRICH[c][0] for c in cats]; pool = [ENRICH[c][1] for c in cats]
    axB.bar([i - w / 2 for i in x], grounded, w, label="grounded set (n=26)", color="#E45756")
    axB.bar([i + w / 2 for i in x], pool, w, label="scored TF pool (n=1839)", color="#BAB0AC")
    axB.set_xticks(list(x)); axB.set_xticklabels(cats)
    axB.set_ylabel("fraction of genes")
    axB.set_title("B. Grounded set is ~10x essential-enriched (the confound)")
    axB.legend(fontsize=8)
    for i, c in enumerate(cats):
        axB.text(i - w / 2, grounded[i] + 0.008, f"{grounded[i]:.0%}", ha="center", fontsize=8)

    fig.suptitle("Component 2 (Arm B, expression) DIAGNOSTICS -- preliminary, pre-v4-rerun",
                 fontsize=11, y=1.02)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(OUT / f"fig_component2_diagnostics.{ext}", bbox_inches="tight", dpi=140)
    print(f"wrote {OUT}/fig_component2_diagnostics.png (+svg)")


if __name__ == "__main__":
    main()
