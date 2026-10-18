#!/usr/bin/env python3
"""Plot historical Component-2 diagnostics, with their evidence limitations visible.

Source: docs/COMPONENT2_RESULTS.md, recorded 2026-07-18, including its later
correction to the essential-gene comparison; see PROJECT_MASTER_SUMMARY.md §8.5.
The original per-perturbation JSONs are not archived in this repository. These
hardcoded, reported aggregates are not independently reproduced results.

A: descriptive v3 seed-0 sensitivity; no sweep significance test was recorded.
B: corrected v4 split within the 162 tested TFs, not the untested 1,839-TF list.
The comparator is outside the core-essential tag, not the NEGv1 reference set.
Wilson intervals are conditional binomial approximations, not uncertainty over
biological replicates, SAE seeds, preprocessing choices or metric versions.

Run: python3 src/make_component2_figures.py
"""
from math import sqrt
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"] = "component2-historical-diagnostics-v2"
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

OUT = Path(__file__).resolve().parents[1] / "docs" / "figures"
N_TESTED = 162

# Reported v3 results: same 162 tested TFs, all-cell SAE training.
SEED_COUNTS = [22, 31, 42, 31, 32]  # seeds 0-4, latent=2048, k=32
SWEEP = [("512 / 32", 64), ("1024 / 32", 46), ("2048 / 32", 22),
         ("2048 / 16", 31), ("2048 / 64", 58)]  # seed 0 only

# Reported later v4 tag breakdown, distinct from the v3 seed-stable set.
# _tag_breakdown in causal_pipeline.py splits membership vs nonmembership.
TAG_SPLIT = [("Core-essential tag", 6, 51), ("Other tested TFs", 14, 111)]

_TITLE = "Historical Replogle K562 diagnostics (all-cell SAE training)"
_DESCRIPTION = (
    "Historical reported aggregates from docs/COMPONENT2_RESULTS.md (2026-07-18, "
    "including the later correction). Original per-perturbation JSONs are not "
    "archived. Panel A: descriptive v3 seed-0 sensitivity. Panel B: corrected v4 "
    "tested-pool split, 6/51 core-essential-tag vs 14/111 other tested TFs; "
    "95% Wilson intervals are binomial approximations. No enrichment, equivalence "
    "or sweep-significance claim. Panels use different metric versions."
)


def wilson_interval(successes, total, z=1.959963984540054):
    """Return a two-sided 95% Wilson score interval (default normal quantile)."""
    proportion = successes / total
    denominator = 1 + z * z / total
    center = (proportion + z * z / (2 * total)) / denominator
    half_width = z * sqrt(proportion * (1 - proportion) / total
                          + z * z / (4 * total * total)) / denominator
    return center - half_width, center + half_width


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(12, 6.0))
    blue, orange, ink, muted = "#4C78A8", "#E66100", "#1f2937", "#4b5563"

    # Panel A deliberately presents observed variation without an inference.
    values = [count / N_TESTED for _, count in SWEEP]
    ax_a.bar(range(len(SWEEP)), values, color=blue, width=0.65, zorder=3)
    seed_rates = [count / N_TESTED for count in SEED_COUNTS]
    ax_a.axhline(sum(seed_rates) / len(seed_rates), ls="--", color=ink,
                 label=f"5-seed mean at 2048 / 32: {sum(seed_rates) / len(seed_rates):.1%}")
    ax_a.axhspan(min(seed_rates), max(seed_rates), color=blue, alpha=0.12,
                 label=f"Observed seed range: {min(seed_rates):.1%}-{max(seed_rates):.1%}")
    for x, (_, count), value in zip(range(len(SWEEP)), SWEEP, values):
        ax_a.text(x, value + 0.013, f"{count}/{N_TESTED}", ha="center", fontsize=9)
    ax_a.set_xticks(range(len(SWEEP)), [label for label, _ in SWEEP], fontsize=9)
    ax_a.set_xlabel("SAE latent dimensions / active features (k)")
    ax_a.set_ylabel("Fraction of tested TFs grounded")
    ax_a.set_ylim(0, 0.54)
    ax_a.set_title("A. v3: recorded seed-0 sensitivity", loc="left", fontsize=11)
    ax_a.legend(fontsize=8, loc="upper right", frameon=False)

    # Panel B uses the tested denominators and does not assert enrichment or equality.
    for x, (label, grounded, total) in enumerate(TAG_SPLIT):
        value = grounded / total
        low, high = wilson_interval(grounded, total)
        ax_b.errorbar(x, value, yerr=[[value - low], [high - value]],
                      fmt="o", color=orange if x == 0 else blue,
                      capsize=7, markersize=9, linewidth=1.8, zorder=3)
        ax_b.text(x, high + 0.012, f"{grounded}/{total} = {value:.1%}",
                  ha="center", fontsize=10, color=ink)
    ax_b.set_xticks(range(len(TAG_SPLIT)), [row[0] for row in TAG_SPLIT], fontsize=10)
    ax_b.set_xlim(-0.65, 1.65)
    ax_b.set_ylim(0, 0.32)
    ax_b.set_ylabel("Fraction grounded within each tested group")
    ax_b.set_title("B. v4: corrected tested-pool split (n=162)", loc="left", fontsize=11)
    ax_b.text(0.5, 0.96, "Points: reported rates; whiskers: 95% Wilson intervals",
              ha="center", va="top", transform=ax_b.transAxes, fontsize=8, color=muted)

    for ax in (ax_a, ax_b):
        ax.yaxis.set_major_formatter(PercentFormatter(xmax=1, decimals=0))
        ax.grid(axis="y", color="#e5e7eb", linewidth=0.8, zorder=0)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.tick_params(axis="both", labelcolor=ink)

    fig.suptitle(_TITLE, fontsize=14, x=0.07, ha="left", y=0.98, color=ink)
    fig.text(0.07, 0.885,
             "Reported historical aggregates; original per-perturbation JSONs are not archived.\n"
             "Descriptive diagnostics, not independently reproduced results.",
             ha="left", fontsize=10, color=muted, linespacing=1.4)
    fig.text(0.07, 0.06,
             "Source: docs/COMPONENT2_RESULTS.md, 2026-07-18, including its later correction. "
             "Panels use different metric versions.\n"
             "B: \"Other tested TFs\" means outside the core-essential tag, not the NEGv1 reference set.\n"
             "Wilson intervals assume binomial sampling; they omit seed/run and biological-replicate uncertainty. "
             "No enrichment or equivalence claim.",
             ha="left", fontsize=8.5, color=muted, linespacing=1.5)
    fig.subplots_adjust(left=0.07, right=0.985, top=0.76, bottom=0.26, wspace=0.32)
    for extension in ("png", "svg"):
        metadata = {"Title": _TITLE, "Description": _DESCRIPTION}
        if extension == "svg":
            metadata.update({"Date": None, "Creator": "make_component2_figures.py"})
        else:
            metadata["Software"] = "make_component2_figures.py"
        path = OUT / f"fig_component2_diagnostics.{extension}"
        fig.savefig(path, bbox_inches="tight", dpi=160, metadata=metadata)
        if extension == "svg":
            # Matplotlib emits trailing blanks in path data; keep the artifact diff clean.
            path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines())
                            + "\n")
    plt.close(fig)
    print("wrote docs/figures/fig_component2_diagnostics.{png,svg}")


if __name__ == "__main__":
    main()
