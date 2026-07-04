#!/usr/bin/env python3
"""V0b: gene-content SAE module definitions (canonical v2).

This script is the single source of truth for V0b. It SUPERSEDES the Colab
notebooks in Drive (`v0b_module_definitions.ipynb`, `_1`, `_2`). An audit of
those notebooks found the canonically-named, most-recently-modified notebook is
actually V1: it uses `CROSS_SEED_THRESHOLD = 3` instead of per-seed
winner-take-all, pads the marker sets with hemoglobin / late-granulocyte genes
that are NOT in the 2,000-HVG set, omits `prob_19Lymph` from the uncommitted
filter (mislabeling ~273 high-Lymph cells), and scores asymmetry with the biased
`H_ery / H_gran` entropy ratio. This file implements the V2 design from
PROJECT_HANDOFF.md section 4.7 and nothing else.

What V0b decides
----------------
Whether "asymmetric modularity" survives a gene-content test: is granulocyte
commitment captured as a unified program (few submodules, one dominant) while
erythroid commitment is distributed across multiple submodules? Every downstream
Phase 2 claim depends on this.

Method (v2)
-----------
1. Per (seed, feature, submodule): hypergeometric enrichment of the feature's
   top-30 decoder genes against canonical, HVG-restricted, submodule-structured
   marker sets. BH-FDR per seed.
2. Per-seed WINNER-TAKE-ALL assignment (NOT cross-seed feature-index matching;
   SAEs do not preserve feature identity across seeds). A feature is assigned to
   the submodule with the smallest q-value, requiring q < 0.05 AND overlap >=
   WINNER_TAKE_ALL_MIN_OVERLAP.
3. Cells are filtered into committed groups by Palantir branch probability. The
   uncommitted test takes the max over ALL SEVEN terminals, including
   prob_19Lymph (the v1 bug fixed here).
4. Asymmetry is scored per seed with two metrics that are robust to the
   different submodule counts (erythroid 3, granulocyte 2):
     - normalized entropy = H(counts) / log(n_submodules) in [0, 1]
       (1.0 = fully distributed, 0.0 = unified)
     - largest-submodule fraction = max(counts) / sum(counts)
       (1.0 = unified)
5. Decision rule (pre-registered): asymmetric modularity is SUPPORTED iff BOTH
   of the following hold in >= 3 of 5 seeds AND both bootstrap CIs exclude 0 in
   the same direction:
     - ery_norm_entropy > gran_norm_entropy   (erythroid more distributed)
     - ery_largest_frac < gran_largest_frac   (erythroid less concentrated)
   If not, the correlation-sign finding does not survive and the claim must be
   softened. That is a real possible outcome and is reported honestly.

Interpretation choices (documented for review)
----------------------------------------------
- "Top-30 decoder genes" = the 30 genes with the largest |decoder weight| for a
  feature's column in W (shape input_dim x latent_dim), matching the notebook.
- Hypergeometric population N = number of HVGs (2000); K = present markers in the
  submodule; n = TOP_K = 30; k = overlap. Upper tail P(X >= k).
- BH-FDR is applied PER SEED across the (feature x submodule) family.
- A submodule must have >= MIN_MARKERS_PRESENT present markers to be testable.
- Bootstrap CIs are computed over the 5 seed-level metric values. With n = 5 this
  is a weak interval; it is reported alongside the per-seed sign-agreement count,
  which is the primary decision criterion. This limitation is stated, not hidden.

Inputs (in --data-dir, default ./data)
  expression_matrix.npy        2730 x 2000 float32, scaled (OPTIONAL: only used
                               for MD5 verification + activation dynamics)
  gene_names.csv               2000 mouse HVGs, MGI symbols, first column
  cell_metadata_palantir.csv   V0a output with prob_<term>_W<id> columns
  sae_seed0.pt ... sae_seed4.pt  5 mouse SAE checkpoints

Outputs (in --output-dir, default ./v0b_outputs)
  module_assignments_v0b.csv      per (seed, feature) assignment  [SHA-256 frozen]
  enrichment_full_v0b.csv         full per-test audit trail
  modularity_metrics_per_seed.csv per-seed normalized entropy + largest fraction
  submodule_dynamics.csv          per-seed binned activations (if X available)
  cell_metadata_v0b.csv           metadata + lineage_group
  v0b_provenance.json             constants, coverage, counts, SHA-256, decision

Run
  python v0b_module_definitions.py --data-dir ./data --output-dir ./v0b_outputs
In Colab, point --data-dir at the Drive folder holding the inputs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

# torch is imported lazily inside the two functions that touch checkpoints
# (load_checkpoints, submodule_dynamics) so the core decision logic can run and
# be unit-tested without a torch install.

# --------------------------------------------------------------------------- #
# Constants (PROJECT_HANDOFF.md 4.7)
# --------------------------------------------------------------------------- #
SEED = 42
N_SEEDS = 5
LATENT_DIM = 128
TOP_K = 30                          # top decoder genes per feature
ENRICHMENT_Q_THRESHOLD = 0.05       # BH-FDR threshold for module assignment
WINNER_TAKE_ALL_MIN_OVERLAP = 2     # require >= 2 marker genes in top-30
MIN_MARKERS_PRESENT = 2             # a submodule needs >= 2 present markers
BRANCH_PROB_COMMITTED = 0.7
BRANCH_PROB_UNCOMMITTED = 0.5
N_BOOTSTRAP = 1000
EXPECTED_X_MD5 = "60183a17983c8b977d036e0f3a58da61"
N_PT_BINS = 40

# v2 marker sets: HVG-present only, organized by lineage submodule.
# Erythroid has 3 submodules, granulocyte has 2 (this asymmetry in submodule
# count is exactly why v2 normalizes entropy by log(n_submodules)).
MARKER_SETS = {
    # ===== ERYTHROID SUBMODULES =====
    "Ery_TF":       ["Gata1", "Klf1", "Tal1", "Lmo2", "Zfpm1", "Stat5a", "Bcl11a", "Myb"],
    # Alas2 is the erythroid-specific heme synthase (the strongest heme gene);
    # it was wrongly omitted before. Verified present in Paul15 2026-06-30.
    "Ery_Heme":     ["Alas2", "Fech", "Hmbs", "Ppox", "Cpox", "Urod", "Alad"],
    # Ermap = erythroid membrane-associated protein, verified present in Paul15.
    "Ery_Membrane": ["Gypc", "Ank1", "Rhag", "Aqp1", "Epor", "Tspo", "Ermap"],
    # Globins under Paul15 (mouse MARS-seq) symbols: Hba-a2, Hbb-b1. The newer
    # Hba-a1/Hbb-bs/bt/y symbols used before are ABSENT from this panel, which is
    # why earlier "erythroid undetected" results were a marker-curation artifact.
    "Ery_Effector": ["Hba-a2", "Hbb-b1"],
    # ===== GRANULOCYTE SUBMODULES =====
    "Gran_TF":      ["Cebpa", "Cebpe", "Runx1"],
    "Gran_Primary": ["Mpo", "Elane", "Prtn3", "Ctsg"],
    # Secondary/specific granule. Also HVG-dropped originally; marker-aware only.
    "Gran_Secondary": ["Ltf", "Lcn2", "Camp", "Ngp", "S100a8", "S100a9", "Mmp8", "Mmp9"],
    # ===== CONTROLS (must not preferentially load on Ery or Gran) =====
    "Progenitor":   ["Kit", "Cd34", "Mllt3", "Eif4ebp1"],
    "Cycling":      ["Top2a", "Pcna", "Mcm2", "Mcm3", "Mcm4", "Mcm5", "Mcm6",
                     "Mcm7", "Cdk1", "Cdk4", "Cdk6", "Cenpe", "Cenpf", "Aurkb",
                     "Birc5", "Ccnb2"],
}
# Lineage submodule lists. Effector/Secondary self-adapt: compute_coverage drops
# them when their genes are absent (original data), includes them under
# marker-aware preprocessing. So the recorded v2/v3 results (globins/late-gran
# absent) are unchanged; the upgraded run gains the effector submodules.
ERY_SUBMODULES = ["Ery_TF", "Ery_Heme", "Ery_Membrane", "Ery_Effector"]
GRAN_SUBMODULES = ["Gran_TF", "Gran_Primary", "Gran_Secondary"]

# The seven Palantir terminal branch-probability columns. prob_19Lymph is
# included in the uncommitted test (the v1 bug fixed in v2).
TERMINAL_COLS = ["prob_1Ery", "prob_14Mo", "prob_16Neu", "prob_13Baso",
                 "prob_11DC", "prob_8Mk", "prob_19Lymph"]
GRAN_BRANCH_COLS = ["prob_14Mo", "prob_16Neu", "prob_13Baso"]


# --------------------------------------------------------------------------- #
# Loading
# --------------------------------------------------------------------------- #
def load_gene_names(data_dir: Path) -> list[str]:
    gene_df = pd.read_csv(data_dir / "gene_names.csv")
    return gene_df.iloc[:, 0].astype(str).tolist()


def load_palantir(data_dir: Path, n_cells_expected: int | None) -> pd.DataFrame:
    meta = pd.read_csv(data_dir / "cell_metadata_palantir.csv")
    meta = meta.loc[:, ~meta.columns.str.startswith("Unnamed")]
    # Strip the V0a _W<id> suffix: prob_1Ery_W37539 -> prob_1Ery
    rename_map = {c: c.split("_W")[0] for c in meta.columns
                 if c.startswith("prob_") and "_W" in c}
    if rename_map:
        meta = meta.rename(columns=rename_map)
    required = ["cell_id", "paul15_clusters", "palantir_pseudotime"] + TERMINAL_COLS
    missing = [c for c in required if c not in meta.columns]
    if missing:
        raise ValueError(f"Palantir metadata missing required columns: {missing}")
    if n_cells_expected is not None and len(meta) != n_cells_expected:
        raise ValueError(f"Row mismatch: meta has {len(meta)} rows, "
                         f"expression matrix has {n_cells_expected}")
    return meta


def load_expression(data_dir: Path) -> np.ndarray | None:
    """Load and MD5-verify the expression matrix. The original matrix must match
    the locked hash; a marker-aware matrix (preprocess_paul15.py) must match its
    OWN recorded hash in preprocess_report.json. Either way it is verified, never
    silently skipped."""
    path = data_dir / "expression_matrix.npy"
    if not path.exists():
        print(f"WARNING: {path} not found. Skipping MD5 verification and "
              f"activation dynamics. Core module decision does not require X.")
        return None
    X = np.load(path)
    actual_md5 = hashlib.md5(X.tobytes()).hexdigest()
    expected, source = EXPECTED_X_MD5, "locked original"
    report = data_dir / "preprocess_report.json"
    if report.exists():
        rec = json.load(open(report)).get("expression_matrix_md5")
        if rec:
            expected, source = rec, "preprocess_report.json (marker-aware)"
    if actual_md5 != expected:
        raise ValueError(f"Expression matrix MD5 mismatch vs {source}! expected "
                         f"{expected}, got {actual_md5}")
    print(f"Expression matrix MD5 verified ({source}): {actual_md5}  shape={X.shape}")
    return X


def load_checkpoints(data_dir: Path, input_dim: int):
    """Load the 5 SAE checkpoints. latent_dim is INFERRED from each checkpoint's
    decoder shape, so overcomplete SAEs (512/1024 latents) load without changes."""
    import torch
    import torch.nn as nn
    import torch.nn.functional as F

    class SparseAutoencoder(nn.Module):
        def __init__(self, input_dim, latent_dim):
            super().__init__()
            self.encoder = nn.Linear(input_dim, latent_dim, bias=True)
            self.decoder = nn.Linear(latent_dim, input_dim, bias=True)

        def forward(self, x):
            z = F.relu(self.encoder(x))
            return self.decoder(z), z

    # Auto-detect all sae_seed*.pt (supports >=10 seeds without editing constants).
    paths = sorted(Path(data_dir).glob("sae_seed*.pt"),
                   key=lambda p: int("".join(c for c in p.stem if c.isdigit()) or "0"))
    if not paths:
        raise FileNotFoundError(f"No sae_seed*.pt in {data_dir}")
    decoder_weights, models = {}, {}
    for seed, path in enumerate(paths):
        ckpt = torch.load(path, map_location="cpu", weights_only=False)
        sd = ckpt.get("model_state", ckpt.get("state_dict", ckpt)) \
            if isinstance(ckpt, dict) else ckpt
        latent_dim = sd["decoder.weight"].shape[1]      # inferred, not hardcoded
        model = SparseAutoencoder(input_dim=input_dim, latent_dim=latent_dim)
        model.load_state_dict(sd)
        model.eval()
        decoder_weights[seed] = model.decoder.weight.detach().numpy()  # (input_dim, latent_dim)
        models[seed] = model
    global N_SEEDS
    N_SEEDS = len(decoder_weights)                      # all downstream loops follow this
    dims = {seed: decoder_weights[seed].shape[1] for seed in decoder_weights}
    print(f"Loaded {N_SEEDS} SAE checkpoints. latent dims: {set(dims.values())}")
    return decoder_weights, models


# --------------------------------------------------------------------------- #
# Marker coverage
# --------------------------------------------------------------------------- #
_MARKERS_SOURCE = "built-in mouse (Paul15)"


def _maybe_override_markers():
    """Optionally replace the built-in MOUSE marker sets with an external panel
    (e.g. human orthologs for the replication arm) when the env var
    V0B_MARKERS_JSON points to a JSON file with keys: marker_sets (dict),
    ery_submodules (list), gran_submodules (list). This leaves the pre-registered
    mouse pipeline byte-for-byte identical when the var is unset, and lets the
    human run reuse the SAME metric/submodule structure with orthologous genes.
    Idempotent: a given file is applied once."""
    global MARKER_SETS, ERY_SUBMODULES, GRAN_SUBMODULES, _MARKERS_SOURCE
    path = os.environ.get("V0B_MARKERS_JSON")
    if not path:
        return
    if _MARKERS_SOURCE == path:
        return
    spec = json.load(open(path))
    MARKER_SETS = {str(k): [str(g) for g in v] for k, v in spec["marker_sets"].items()}
    ERY_SUBMODULES = [str(m) for m in spec["ery_submodules"]]
    GRAN_SUBMODULES = [str(m) for m in spec["gran_submodules"]]
    _MARKERS_SOURCE = path
    print(f"Marker sets OVERRIDDEN from {path} "
          f"({spec.get('species', 'unknown')}); {len(MARKER_SETS)} submodules.")


def compute_coverage(gene_names: list[str]):
    _maybe_override_markers()
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    coverage = {}
    for module, genes in MARKER_SETS.items():
        present = [g for g in genes if g in gene_to_idx]
        coverage[module] = {
            "requested": len(genes),
            "present": len(present),
            "present_genes": present,
            "missing_genes": [g for g in genes if g not in gene_to_idx],
            "present_indices": [gene_to_idx[g] for g in present],
        }
    usable = [m for m, c in coverage.items() if c["present"] >= MIN_MARKERS_PRESENT]
    print("\nMarker coverage (present / requested):")
    for m, c in coverage.items():
        flag = "OK" if c["present"] >= MIN_MARKERS_PRESENT else "DROPPED (low coverage)"
        print(f"  {m:14s} {c['present']:2d}/{c['requested']:2d}  {flag}")
        if c["missing_genes"]:
            print(f"     missing from HVG set: {c['missing_genes']}")
    dropped_lineage = [m for m in ERY_SUBMODULES + GRAN_SUBMODULES if m not in usable]
    if dropped_lineage:
        print(f"\nWARNING: lineage submodules dropped for low HVG coverage: "
              f"{dropped_lineage}. This changes the submodule counts used in the "
              f"asymmetry metric and must be reported.")
    return coverage, usable


# --------------------------------------------------------------------------- #
# Enrichment + assignment
# --------------------------------------------------------------------------- #
def run_enrichment(decoder_weights, coverage, usable, n_genes):
    records = []
    for seed in range(N_SEEDS):
        W = decoder_weights[seed]
        for feat in range(W.shape[1]):            # actual latent dim (supports overcomplete)
            top_idx = set(np.argsort(np.abs(W[:, feat]))[::-1][:TOP_K].tolist())
            for module in usable:
                marker_idx = set(coverage[module]["present_indices"])
                K = len(marker_idx)
                k = len(top_idx & marker_idx)
                p = stats.hypergeom.sf(k - 1, n_genes, K, TOP_K) if k > 0 else 1.0
                records.append({"seed": seed, "feature_idx": feat, "module": module,
                                "overlap": k, "module_size": K, "p_value": p})
    df = pd.DataFrame(records)
    df["q_value"] = np.nan
    for seed in range(N_SEEDS):
        mask = df["seed"] == seed
        _, q, _, _ = multipletests(df.loc[mask, "p_value"].values, method="fdr_bh")
        df.loc[mask, "q_value"] = q
    df["significant"] = (df["q_value"] < ENRICHMENT_Q_THRESHOLD) & \
                        (df["overlap"] >= WINNER_TAKE_ALL_MIN_OVERLAP)
    return df


def winner_take_all(enrich_df):
    """Per (seed, feature): assign to the significant submodule with smallest q,
    tiebreak by larger overlap. UNASSIGNED if no submodule qualifies."""
    rows = []
    sig = enrich_df[enrich_df["significant"]]
    feats = sorted(enrich_df["feature_idx"].unique())     # actual latent dim
    for seed in range(N_SEEDS):
        for feat in feats:
            cand = sig[(sig["seed"] == seed) & (sig["feature_idx"] == feat)]
            if len(cand):
                best = cand.sort_values(["q_value", "overlap"],
                                        ascending=[True, False]).iloc[0]
                rows.append({"seed": seed, "feature_idx": feat,
                             "module": best["module"], "overlap": int(best["overlap"]),
                             "q_value": float(best["q_value"])})
            else:
                rows.append({"seed": seed, "feature_idx": feat, "module": "UNASSIGNED",
                             "overlap": 0, "q_value": np.nan})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Cell filtering
# --------------------------------------------------------------------------- #
def assign_cell_groups(meta: pd.DataFrame) -> pd.DataFrame:
    meta = meta.copy()
    meta["prob_gran_combined"] = meta[GRAN_BRANCH_COLS].sum(axis=1)
    meta["max_branch_prob"] = meta[TERMINAL_COLS].max(axis=1)  # includes prob_19Lymph

    def group(row):
        if row["prob_1Ery"] > BRANCH_PROB_COMMITTED:
            return "committed_erythroid"
        if row["prob_gran_combined"] > BRANCH_PROB_COMMITTED:
            return "committed_granulocyte"
        if row["max_branch_prob"] < BRANCH_PROB_UNCOMMITTED:
            return "uncommitted_progenitor"
        return "intermediate"

    meta["lineage_group"] = meta.apply(group, axis=1)
    return meta


# --------------------------------------------------------------------------- #
# Modularity metrics (v2)
# --------------------------------------------------------------------------- #
def _entropy(counts):
    total = float(sum(counts))
    if total == 0:
        return 0.0
    p = [c / total for c in counts if c > 0]
    return float(-sum(pi * np.log(pi) for pi in p))


def _norm_entropy(counts):
    """Shannon entropy normalized by log(n_submodules) -> [0, 1]."""
    n = len(counts)
    if n <= 1:
        return 0.0
    return _entropy(counts) / np.log(n)


def _largest_fraction(counts):
    total = float(sum(counts))
    if total == 0:
        return np.nan
    return max(counts) / total


def per_seed_modularity(assignments, ery_subs, gran_subs):
    rows = []
    for seed in range(N_SEEDS):
        s = assignments[assignments["seed"] == seed]
        ery_counts = [int((s["module"] == m).sum()) for m in ery_subs]
        gran_counts = [int((s["module"] == m).sum()) for m in gran_subs]
        rows.append({
            "seed": seed,
            "ery_total": sum(ery_counts), "gran_total": sum(gran_counts),
            "ery_counts": ery_counts, "gran_counts": gran_counts,
            "ery_norm_entropy": _norm_entropy(ery_counts),
            "gran_norm_entropy": _norm_entropy(gran_counts),
            "ery_largest_frac": _largest_fraction(ery_counts),
            "gran_largest_frac": _largest_fraction(gran_counts),
        })
    return pd.DataFrame(rows)


def decide(mod_df: pd.DataFrame, rng: np.random.Generator):
    """Apply the pre-registered decision rule. Seeds where either lineage has 0
    assigned features are 'undefined' and excluded from the sign-agreement
    denominator (reported explicitly)."""
    valid = mod_df[(mod_df["ery_total"] > 0) & (mod_df["gran_total"] > 0)].copy()
    n_valid = len(valid)

    entropy_sign = valid["ery_norm_entropy"] > valid["gran_norm_entropy"]
    frac_sign = valid["ery_largest_frac"] < valid["gran_largest_frac"]
    both = entropy_sign & frac_sign
    n_both = int(both.sum())

    def boot_ci(values):
        values = np.asarray(values, dtype=float)
        if len(values) == 0:
            return (np.nan, np.nan)
        means = [rng.choice(values, size=len(values), replace=True).mean()
                 for _ in range(N_BOOTSTRAP)]
        return tuple(np.percentile(means, [2.5, 97.5]))

    d_entropy = (valid["ery_norm_entropy"] - valid["gran_norm_entropy"]).values
    d_frac = (valid["gran_largest_frac"] - valid["ery_largest_frac"]).values
    ci_entropy = boot_ci(d_entropy)
    ci_frac = boot_ci(d_frac)

    entropy_excludes_0_pos = ci_entropy[0] > 0
    frac_excludes_0_pos = ci_frac[0] > 0
    supported = (n_both >= 3) and entropy_excludes_0_pos and frac_excludes_0_pos

    return {
        "n_valid_seeds": n_valid,
        "n_seeds_both_directions_hold": n_both,
        "seed_agreement_threshold": 3,
        "mean_diff_norm_entropy_ery_minus_gran": float(np.mean(d_entropy)) if n_valid else None,
        "bootstrap_ci_norm_entropy_diff": [float(ci_entropy[0]), float(ci_entropy[1])],
        "mean_diff_largest_frac_gran_minus_ery": float(np.mean(d_frac)) if n_valid else None,
        "bootstrap_ci_largest_frac_diff": [float(ci_frac[0]), float(ci_frac[1])],
        "entropy_ci_excludes_0_positive": bool(entropy_excludes_0_pos),
        "frac_ci_excludes_0_positive": bool(frac_excludes_0_pos),
        "asymmetric_modularity_supported": bool(supported),
        "bootstrap_caveat": ("CIs are over 5 seed-level values; with n=5 the "
                             "interval is weak. The >=3/5 sign-agreement count is "
                             "the primary criterion."),
    }


# --------------------------------------------------------------------------- #
# Optional activation dynamics
# --------------------------------------------------------------------------- #
def submodule_dynamics(X, models, meta, assignments, usable):
    import torch
    X_t = torch.tensor(X, dtype=torch.float32)
    acts, recon_mse = {}, {}
    for seed in range(N_SEEDS):
        with torch.no_grad():
            recon, z = models[seed](X_t)
        acts[seed] = z.numpy()
        recon_mse[seed] = float(((recon - X_t) ** 2).mean().item())
    mean_acts = np.stack([acts[s] for s in range(N_SEEDS)]).mean(axis=0)

    meta = meta.copy()
    meta["pt_bin"] = pd.cut(meta["palantir_pseudotime"], bins=N_PT_BINS, labels=False)
    # Features assigned to each submodule in >= 3 of 5 seeds (robust set)
    robust = {}
    for module in usable:
        per_feat = (assignments[assignments["module"] == module]
                    .groupby("feature_idx").size())
        robust[module] = per_feat[per_feat >= 3].index.tolist()

    rows = []
    for grp in ["committed_erythroid", "committed_granulocyte", "uncommitted_progenitor"]:
        cells = meta[meta["lineage_group"] == grp]
        if cells.empty:
            continue
        for module, feats in robust.items():
            if not feats:
                continue
            module_act = mean_acts[cells.index.values][:, feats].mean(axis=1)
            tmp = cells.assign(module_act=module_act)
            for pt_bin, val in tmp.groupby("pt_bin")["module_act"].mean().items():
                rows.append({"lineage_group": grp, "module": module,
                             "pt_bin": pt_bin, "mean_activation": float(val)})
    return pd.DataFrame(rows), recon_mse


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default=os.environ.get("V0B_DATA_DIR", "./data"))
    ap.add_argument("--output-dir", default=os.environ.get("V0B_OUTPUT_DIR", "./v0b_outputs"))
    args = ap.parse_args()

    data_dir = Path(args.data_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    np.random.seed(SEED)
    # torch needs no seeding here: checkpoint load + eval are deterministic, and
    # the only stochastic step (bootstrap) uses the numpy generators above.

    gene_names = load_gene_names(data_dir)
    n_genes = len(gene_names)
    X = load_expression(data_dir)
    n_cells = X.shape[0] if X is not None else None
    if X is not None and X.shape[1] != n_genes:
        raise ValueError(f"gene_names ({n_genes}) != expression cols ({X.shape[1]})")
    meta = load_palantir(data_dir, n_cells)
    decoder_weights, models = load_checkpoints(data_dir, input_dim=n_genes)

    coverage, usable = compute_coverage(gene_names)
    ery_subs = [m for m in ERY_SUBMODULES if m in usable]
    gran_subs = [m for m in GRAN_SUBMODULES if m in usable]

    enrich_df = run_enrichment(decoder_weights, coverage, usable, n_genes)
    assignments = winner_take_all(enrich_df)
    meta = assign_cell_groups(meta)
    group_counts = meta["lineage_group"].value_counts().to_dict()
    print("\nCell group counts:", group_counts)

    mod_df = per_seed_modularity(assignments, ery_subs, gran_subs)
    decision = decide(mod_df, rng)

    # Outputs
    assign_path = out_dir / "module_assignments_v0b.csv"
    assignments.to_csv(assign_path, index=False)
    enrich_df.to_csv(out_dir / "enrichment_full_v0b.csv", index=False)
    mod_df.to_csv(out_dir / "modularity_metrics_per_seed.csv", index=False)
    meta.to_csv(out_dir / "cell_metadata_v0b.csv", index=False)

    recon_mse = None
    if X is not None:
        dyn_df, recon_mse = submodule_dynamics(X, models, meta, assignments, usable)
        dyn_df.to_csv(out_dir / "submodule_dynamics.csv", index=False)

    with open(assign_path, "rb") as f:
        table_sha = hashlib.sha256(f.read()).hexdigest()

    provenance = {
        "version": "v2",
        "constants": {
            "TOP_K": TOP_K, "ENRICHMENT_Q_THRESHOLD": ENRICHMENT_Q_THRESHOLD,
            "WINNER_TAKE_ALL_MIN_OVERLAP": WINNER_TAKE_ALL_MIN_OVERLAP,
            "MIN_MARKERS_PRESENT": MIN_MARKERS_PRESENT,
            "BRANCH_PROB_COMMITTED": BRANCH_PROB_COMMITTED,
            "BRANCH_PROB_UNCOMMITTED": BRANCH_PROB_UNCOMMITTED,
            "N_BOOTSTRAP": N_BOOTSTRAP, "SEED": SEED,
        },
        "expression_matrix_md5_verified": X is not None,
        "expression_matrix_md5_expected": EXPECTED_X_MD5,
        "n_genes": n_genes, "n_cells": n_cells,
        "marker_coverage": {m: {k: v for k, v in c.items() if k != "present_indices"}
                            for m, c in coverage.items()},
        "usable_submodules": usable,
        "ery_submodules_used": ery_subs, "gran_submodules_used": gran_subs,
        "cell_group_counts": group_counts,
        "recon_mse_per_seed": recon_mse,
        "module_assignments_v0b_sha256": table_sha,
        "decision": decision,
    }
    with open(out_dir / "v0b_provenance.json", "w") as f:
        json.dump(provenance, f, indent=2)

    # Summary
    print("\n" + "=" * 72)
    print("V0b PER-SEED MODULARITY METRICS")
    print("=" * 72)
    print(mod_df[["seed", "ery_total", "gran_total", "ery_norm_entropy",
                  "gran_norm_entropy", "ery_largest_frac", "gran_largest_frac"]]
          .to_string(index=False))
    print("\n" + "=" * 72)
    print("DECISION")
    print("=" * 72)
    for k, v in decision.items():
        print(f"  {k}: {v}")
    verdict = ("SUPPORTED" if decision["asymmetric_modularity_supported"]
               else "NOT SUPPORTED -> soften the asymmetric modularity claim")
    print(f"\n  ASYMMETRIC MODULARITY: {verdict}")
    print("\n" + "=" * 72)
    print(f"PRE-REGISTRATION: module_assignments_v0b.csv SHA-256\n  {table_sha}")
    print("  Commit this hash to the repo BEFORE any Phase 2 data inspection.")
    print("=" * 72)


if __name__ == "__main__":
    main()
