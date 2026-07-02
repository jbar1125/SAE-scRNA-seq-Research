# COMPONENT 0 RESULTS (Gate 0: the foundation, EXECUTED)

This file records the ACTUAL executed Gate-0 run on real data, in-container, on
2026-07-02. It supersedes the "needs a GPU run" status in COMPONENT0_STATUS.md for
the mouse arm. Every number here was produced by version-controlled code
(`train_sae.py`, `v0b_v3_1_decision.py`, `baselines_nmf_pca.py`,
`marker_sensitivity.py`) on the marker-aware, corrected-marker Paul15 matrix and,
for the replication arm, on the Setty 2019 human CD34+ marrow atlas.

Read the honest one-liner first:

> **The original "asymmetric modularity" headline (erythroid distributed across
> many programs, granulocyte a single unified program) is NOT SUPPORTED.** On the
> sparsity-corrected, 10-seed SAEs the data show the REVERSE and it is seed-robust:
> granulocyte splits into 3 real above-control programs, erythroid concentrates
> into 1 dominant hemoglobin program. AND — the more important finding for rigor —
> **the modularity verdict is method-dependent**: PCA, NMF, and the SAE give three
> different answers on the identical matrix. No unsupervised decomposition recovers
> a method-invariant "gene program" structure. That is why the causal test
> (Component 2) is the real arbiter, not any single decomposition.

---

## 1. What was run (provenance)

| Arm | Data | Cells x genes | Preprocess | SAE |
|-----|------|---------------|-----------|-----|
| Mouse (primary) | Paul15 myeloid atlas (git-LFS mirror, sha256 bd75eb90...) | 2730 x 2012 | log1p + marker-aware HVG (corrected globins Hba-a2/Hbb-b1 + Alas2 + Ermap), z-scaled. MD5 `184ea5dd...` | 512 latents, l1=0.8, 300 epochs, 10 seeds |
| Human (replication) | Setty 2019 CD34+ marrow (Palantir, sha256 be4f8603...) | 4142 x 2032 | normalize_total(1e4) + log1p + marker-aware HVG (human orthologs), z-scaled. MD5 `62796189...` | 512 latents, l1 retuned for band, 300 epochs, 10 seeds |

Both matrices are MD5-verified at load. Marker panels: mouse = built-in
`v0b.MARKER_SETS`; human = `human_markers.json` (1:1 orthologs of the same
submodules), injected via the `V0B_MARKERS_JSON` override so the pre-registered
metric and submodule structure are byte-for-byte reused, only gene symbols differ.

## 2. Sparsity was fixed (the QC gate that the earlier runs failed)

The pre-registered QC band is mean L0 (active latents/cell) in **20-50**. A key
lesson from this run: **L0 depends on training length**, not just L1. At 150
epochs l1=3.0 gave L0 34.6; at the real 300 epochs the same L1 collapsed to L0 ~17
(too sparse). Retuning L1 at the final epoch count:

Mouse 300-epoch L1 sweep (1 seed each): l1=0.5 -> L0 67.4 (dense); l1=0.75 -> 38.9;
**l1=0.8 -> 37.0 (chosen, band center)**; l1=1.0 -> 28.0; l1=1.5 -> 19.4; l1=2.0 -> 17.4.

Mouse final (l1=0.8, 10 seeds, 300 epochs): **L0 mean 35.75, range [34.6, 37.0],
all 10 seeds in band.** recon_mse 0.4642 (range [0.4634, 0.4651]) - extremely
stable across seeds. This is a properly-regularized run; the earlier under-/over-
regularized runs that produced the v1/v2/v3 numbers are not comparable and are
retired.

## 3. Mouse primary result: pre-registered v3.1 decision

Frozen metric (`v0b_v3_1_decision.py`, spec SHA `fc342829...`): per seed, a
submodule is a REAL program iff its strength exceeds the Progenitor baseline by
> 0.10 AND its abundance-matched permutation q < 0.05. The original claim is
SUPPORTED only if erythroid shows >= 2 real programs and is MORE distributed than
granulocyte in >= 60% of seeds.

Submodule strength (mean over 10 seeds; 1.0 = null) + abundance-matched perm q:

| submodule | strength | perm_q | real program? |
|-----------|---------:|-------:|:-------------:|
| Ery_Effector (globins) | 5.07 | 0.0018 | yes (erythroid's one dominant program) |
| Gran_Primary | 4.30 | 0.0018 | yes |
| Gran_Secondary | 3.36 | 0.0018 | yes |
| Gran_TF | 2.93 | 0.0154 | yes |
| Ery_Heme | 2.40 | 0.0018 | borderline (excess ~0.07-0.18 vs floor 0.10) |
| **Progenitor (control)** | **2.33** | **0.085** | **no (correctly non-significant)** |
| Cycling (2nd control) | 2.01 | 0.0018 | real program (expected: cell cycle is real) |
| Ery_Membrane | 1.99 | 0.0105 | below control baseline |
| Ery_TF | 1.90 | 0.0360 | below control baseline |

Per-seed control-referenced counts: **gran_n_real = 3 in all 10/10 seeds**;
ery_n_real = 2 in 5 seeds, 1 in 5 seeds (median 1.5). gran excess-entropy
(0.82-0.95) exceeds ery excess-entropy (0-0.24) in **all 10/10 seeds**.

**Verdict: asymmetric modularity (original direction) NOT SUPPORTED** (0/10 seeds
meet the pattern; threshold was 6). Plain finding: granulocyte distributed across
3 real programs; erythroid unified into 1 dominant hemoglobin program - the
REVERSE of the v1 headline, consistent with PROJECT_AUDIT.md section J.

Null calibration check: the Progenitor control lands at perm_q = 0.085
(non-significant), exactly where a well-calibrated null should sit. This is the fix
for the earlier mis-calibrated null in which the control tested significant.

## 4. Baselines: the modularity verdict is METHOD-DEPENDENT (the real finding)

`baselines_nmf_pca.py` runs the SAME v3.1 decision on PCA and NMF decoders (rank
512, matched to the SAE) and reports effective dimensionality.

Effective dimensionality (participation ratio, PR = (Σλ)²/Σλ²):
- **SAE: 195.9** | PCA (data): 14.0 | NMF: 5.2.
  The SAE learns a genuinely high-dimensional distributed code; PCA/NMF compress
  variance into a handful of directions.

Modularity verdict per representation (identical data, identical metric):
- **PCA**: ery_n_real=1, gran_n_real=1 -> resolves almost nothing (too coarse).
- **NMF**: ery_n_real=2.5, gran_n_real=2 -> "SUPPORTED" (finds erythroid MORE
  distributed - the OPPOSITE of the SAE; an artifact of parts-based non-negative
  decomposition, and its Gran_Primary q=0.18 is non-significant).
- **SAE**: ery_n_real=1.5, gran_n_real=3 -> NOT supported (granulocyte distributed).

**Three decompositions, three different answers.** This is the artifact-vs-signal
result the project's novelty rests on: "gene program modularity" is not a property
of the data alone; it is co-determined by the chosen decomposition. Any single-
method modularity claim (including the v1 headline and including this SAE result)
is method-contingent, not biological ground truth.

## 5. Marker robustness: verdict is bulletproof to single-marker drops

`marker_sensitivity.py` (leave-one-marker-out, 34 present markers): **0/34 drops
flip the verdict.** gran_n_real stays 3.0 throughout; ery_n_real moves within
1.0-2.0 (dropping Alas2 pulls Ery_Heme below the floor; dropping either globin
retires the 2-gene Ery_Effector submodule) but never enough to support the
original claim. This directly answers the "did you cherry-pick markers?" concern
that the earlier globin-symbol bug made concrete.

## 6. Human replication arm

IN PROGRESS (this commit). The human CD34+ marrow atlas has been obtained
(sha256 be4f8603...), preprocessed with the CELLxGENE pipeline (normalize_total +
log1p + marker-aware HVG on human orthologs; all 9 submodules usable, MD5
`62796189...`), and the human L1 is being retuned to the 20-50 L0 band (human data
are denser: l1=0.4 -> L0 99.5, l1=0.6 -> L0 54.4, higher L1 needed than mouse).
The pre-registered v3.1 decision will then run with `V0B_MARKERS_JSON=human_markers.json`.
This section is filled with the human verdict in the finalizing commit. The
question it answers: does the SAE's granulocyte-distributed / erythroid-unified
direction hold cross-species, and does the PCA/NMF/SAE method-disagreement recur?

## 7. Freeze (provenance hashes)

The Gate-0 mouse decision output is frozen for provenance:
- `data_g0/v0b_outputs/v0b_v3_1_decision.json` sha256 `0e88f982...`
- `data_g0/v0b_outputs/v0b_v3_1_per_seed.csv`   sha256 `683eea9d...`
- combined canonical bundle sha256
  `172861417a923a83705812a96e5ad1e555055e968323289879447e53beb84396`

This freezes a NEGATIVE-plus-method-dependence result from a properly regularized,
10-seed, sparsity-corrected run - not a positive claim from a broken run. That
distinction is the whole point of the gate (CLAUDE.md section 2).

## 8. Honest limitations (do not omit)

- n=10 mouse seeds, one dataset per species. The reversed SAE direction is seed-
  robust but the cross-method disagreement means it is not a species/decomposition-
  invariant biological fact.
- The 20-50 L0 band and the 0.10 excess floor are project conventions; the pre-
  registration froze them before this run, but they are choices.
- Ery_Effector is a 2-gene submodule in mouse (globins); dropping either globin
  retires it. Human uses 3 globins (HBB/HBA1/HBA2), a modest robustness gain.
- "Strength" is decoder-mass concentration, a structural (not causal) measure.
  Whether these programs are causally real is exactly what Component 2 (Replogle
  CRISPRi head-to-head vs Kendiukhov's 6.2% embedding-space null) is designed to
  test. Gate 0 delivers the calibrated, method-aware foundation for that test; it
  does not itself settle biology.
