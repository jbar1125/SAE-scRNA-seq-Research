#!/usr/bin/env python3
"""Generate the Gate-0 result figures from the committed JSON artifacts.

Reproducible, deterministic; reads only the small committed decision artifacts (no
retraining). Writes PNG + SVG to docs/figures/. Design: diverging encoding for the
polarity question (erythroid-more-distributed vs granulocyte-more-distributed), a
CVD-safe orange/blue/gray palette, sign encoded by position AND color (not color
alone), direct value labels, recessive axes.

Run: python3 src/make_figures.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"] = "gate0"      # deterministic SVG element ids
import matplotlib.pyplot as plt

_META = {"Date": None}                              # drop the embedded timestamp -> byte-stable

OUT = Path("docs/figures"); OUT.mkdir(parents=True, exist_ok=True)

# CVD-safe diverging pair + neutral (Wong-style orange / blue, gray midpoint)
WARM, COOL, NEU, INK, MUT = "#E66100", "#1A85FF", "#9CA3AF", "#1f2937", "#6b7280"


def _dec(path, *keys):
    d = json.load(open(path))
    for k in keys:
        d = d[k]
    return d


def load_method_dependence():
    b_m = "data_g0/baselines/baselines_results.json"
    b_h = "data_g0_human/baselines/baselines_results.json"
    rows = {}
    for m, k in [("PCA", "pca"), ("NMF", "nmf"), ("SAE", "sae")]:
        rows[m] = {"mouse": _dec(b_m, "modularity", k, "decision"),
                   "human": _dec(b_h, "modularity", k, "decision")}
    rows["GRN"] = {"mouse": _dec("data_g0/grn/grn_baseline.json", "modularity", "decision"),
                   "human": _dec("data_g0_human/grn/grn_baseline.json", "modularity", "decision")}
    rows["ICA"] = {"mouse": _dec("data_g0/ica/ica_baseline.json", "modularity", "decision"),
                   "human": _dec("data_g0_human/ica/ica_baseline.json", "modularity", "decision")}
    return rows


def fig_method_dependence():
    rows = load_method_dependence()
    order = ["PCA", "ICA", "NMF", "SAE", "GRN"]     # display order
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), sharex=True)
    for ax, sp in zip(axes, ["mouse", "human"]):
        ys = range(len(order))
        for y, m in zip(ys, order):
            d = rows[m][sp]
            diff = float(d["median_ery_n_real"]) - float(d["median_gran_n_real"])
            supp = d["asymmetric_modularity_supported"]
            color = WARM if diff > 0 else (COOL if diff < 0 else NEU)
            ax.barh(y, diff, height=0.62, color=color, zorder=3)
            # direct value label, offset to the open side
            off = 0.06 if diff >= 0 else -0.06
            ha = "left" if diff >= 0 else "right"
            lbl = f"{d['median_ery_n_real']:.1f}/{d['median_gran_n_real']:.1f}"
            ax.text(diff + off, y, lbl, va="center", ha=ha, fontsize=9, color=INK, zorder=4)
            # the lone supporter is carried by the diverging color (only warm bar)
            # plus the subtitle callout; no inline text needed.
        ax.axvline(0, color=INK, lw=1.2, zorder=2)
        ax.set_yticks(list(ys)); ax.set_yticklabels(order, fontsize=10)
        ax.set_xlim(-2.6, 2.0); ax.invert_yaxis()
        ax.set_title(sp.capitalize(), fontsize=11, color=INK, loc="left", pad=8)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.spines["bottom"].set_color(MUT)
        ax.tick_params(length=0)
        ax.set_xticks([-2, -1, 0, 1])
        ax.grid(axis="x", color="#e5e7eb", lw=0.8, zorder=0)
    axes[0].set_ylabel("decomposition", fontsize=10, color=MUT)
    fig.text(0.5, 0.015,
             "erythroid − granulocyte real-program count      "
             "(← granulocyte more distributed    |    erythroid more distributed →)",
             ha="center", fontsize=9.5, color=MUT)
    fig.suptitle("Gene-program modularity is decomposition-dependent",
                 fontsize=13.5, color=INK, x=0.09, ha="left", y=0.98)
    fig.text(0.09, 0.865,
             "median real-program counts (erythroid/granulocyte); the original claim needs "
             "erythroid > granulocyte.\nIt survives in only 1 of 10 method×species cells "
             "(mouse NMF, the single orange bar).",
             ha="left", fontsize=8.8, color=MUT, linespacing=1.35)
    fig.subplots_adjust(left=0.09, right=0.985, top=0.72, bottom=0.14, wspace=0.28)
    for ext in ("png", "svg"):
        fig.savefig(OUT / f"fig_method_dependence.{ext}", dpi=200, metadata=_META)
    plt.close(fig)
    print("wrote docs/figures/fig_method_dependence.{png,svg}")


def fig_robustness():
    """Compact diverging strip: granulocyte-minus-erythroid across every stress test
    (mouse). All >= 0 -> the verdict never flips toward the original claim."""
    pts = []
    r2 = json.load(open("data_g0/robustness/r2_l0_band_summary.json"))["points"]
    for p in r2:
        pts.append((f"L0={p['mean_l0']:.0f}", p["median_gran_n_real"] - p["median_ery_n_real"]))
    r4 = json.load(open("data_g0/robustness/r4_hvg_summary.json"))["points"]
    for p in r4:
        pts.append((f"{p['n_hvg']} HVG", p["median_gran_n_real"] - p["median_ery_n_real"]))
    r3 = json.load(open("data_g0/robustness/r3_20seed_summary.json"))
    pts.append(("20 seeds", r3["median_gran_n_real"] - r3["median_ery_n_real"]))

    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    ys = range(len(pts))
    for y, (lab, v) in zip(ys, pts):
        ax.barh(y, v, height=0.6, color=(COOL if v > 0 else NEU), zorder=3)
        ax.text(v + 0.05, y, f"+{v:.1f}" if v > 0 else f"{v:.1f}",
                va="center", ha="left", fontsize=9, color=INK, zorder=4)
    ax.axvline(0, color=INK, lw=1.2)
    ax.set_yticks(list(ys)); ax.set_yticklabels([p[0] for p in pts], fontsize=10)
    ax.invert_yaxis(); ax.set_xlim(-0.5, 3.0)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(MUT); ax.tick_params(length=0)
    ax.grid(axis="x", color="#e5e7eb", lw=0.8, zorder=0)
    ax.set_title("Mouse verdict is stable across every stress test", fontsize=12.5,
                 color=INK, loc="left", pad=10)
    fig.text(0.125, 0.885, "granulocyte − erythroid real programs; > 0 everywhere "
             "= never the original claim", ha="left", fontsize=9, color=MUT)
    ax.set_xlabel("granulocyte − erythroid real-program count", fontsize=9.5, color=MUT)
    fig.subplots_adjust(left=0.18, right=0.97, top=0.80, bottom=0.12)
    for ext in ("png", "svg"):
        fig.savefig(OUT / f"fig_robustness.{ext}", dpi=200, metadata=_META)
    plt.close(fig)
    print("wrote docs/figures/fig_robustness.{png,svg}")


if __name__ == "__main__":
    fig_method_dependence()
    fig_robustness()
