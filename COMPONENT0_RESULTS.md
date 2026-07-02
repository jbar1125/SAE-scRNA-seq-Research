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

Human CD34+ marrow atlas (Setty 2019, sha256 be4f8603...), CELLxGENE pipeline
(normalize_total(1e4) + log1p + marker-aware HVG on human orthologs; all 9
submodules usable, MD5 `62796189...`). Human L1 retuned to band (human is denser:
l1=0.4->L0 99.5, 0.6->54.4, **0.8->32.5 chosen**, 1.0->21.4). Final 10-seed run:
**L0 mean 32.5, all 10 in band**, recon_mse ~0.666 (extremely tight). Analysis run
with `V0B_MARKERS_JSON=human_markers.json` so the pre-registered metric and
submodule structure are reused with human orthologs.

**v3.1 verdict (human): original claim NOT SUPPORTED (0/10 seeds)** - same as mouse.
median gran_n_real = 2 vs ery_n_real = 1. Submodule strength (1.0=null) + perm_q:

| submodule | strength | perm_q | note |
|-----------|---------:|-------:|------|
| Gran_Primary | 5.44 | 0.0045 | strongest program in human |
| Gran_Secondary | 4.02 | 0.043 | real program |
| Ery_Effector (globins) | 3.75 | 0.216 | present but NOT significant vs abundance-matched null |
| Cycling | 3.67 | 0.0045 | real (expected) |
| Ery_TF | 2.74 | 0.018 | the one erythroid submodule that clears |
| Gran_TF | 2.74 | 0.298 | n.s. |
| Ery_Membrane | 2.51 | 0.304 | n.s. |
| **Progenitor (control)** | **2.22** | **0.748** | **non-significant (null calibrated)** |
| Ery_Heme | 2.21 | 0.736 | at control |

**Cross-species reading (honest):**
- The DIRECTION replicates: granulocyte carries more real above-control programs
  than erythroid in BOTH species; the original "erythroid distributed / granulocyte
  unified" claim is rejected in mouse AND human.
- The MECHANISM differs by species and it is worth stating plainly. In mouse the
  dominant, highly-significant program is erythroid hemoglobin (`Ery_Effector`); in
  human the dominant program is the granulocyte primary granule, and the erythroid
  globin program is present but NOT significant against the abundance-matched null.
  Most likely cause: Setty is CD34+ progenitor-SELECTED marrow, so terminal
  hemoglobin-high erythroblasts are underrepresented relative to Paul15. This is a
  sampling/composition difference, not a contradiction of the direction.
- Control calibration holds cross-species (human Progenitor perm_q 0.748).

**Human baselines** (same v3.1 metric on PCA/NMF decoders, rank 512):
- Effective dimensionality (PR): **SAE 278.6, PCA 270.7, NMF 5.9.** Note the human
  contrast with mouse: here PCA is nearly as high-dimensional as the SAE, while NMF
  is again heavily compressed.
- **PCA**: ery=1, gran=3, NOT supported -> granulocyte-distributed (AGREES with SAE).
- **NMF**: ery=0, gran=0, NOT supported -> detects no above-control program at all.
- **SAE**: ery=1, gran=2, NOT supported -> granulocyte-distributed.

**Human marker sensitivity** (leave-one-out, 36 present markers): **0/36 drops flip
the verdict.** Robust, same as mouse.

## 6b. Cross-species x cross-method synthesis (the money table)

Each cell = does that method+species SUPPORT the ORIGINAL claim (erythroid more
distributed), and the plain direction it points:

| method | mouse | human |
|--------|-------|-------|
| SAE  | NO - granulocyte more distributed (gran 3 / ery 1.5) | NO - granulocyte more distributed (gran 2 / ery 1) |
| PCA  | NO - resolves almost nothing (ery 1 / gran 1) | NO - granulocyte more distributed (gran 3 / ery 1) |
| NMF  | **YES** - erythroid more distributed (ery 2.5 / gran 2) | NO - detects nothing (ery 0 / gran 0) |

Two things to read off this:
1. **The original claim survives in exactly 1 of 6 method x species combinations**
   (mouse NMF), and that one is the method most prone to parts-based artifact
   (inflated strengths, a non-significant Gran_Primary). Every other combination
   rejects it. So the v1 headline is not just unsupported by the SAE; it is
   unsupported by the weight of methods and species.
2. **The SAE direction (granulocyte more distributed) is the only claim that holds
   in both species**, and in human PCA independently agrees. But NMF disagrees in
   both species (opposite in mouse, null in human). So even the reversed finding is
   method-contingent, not decomposition-invariant. This is the artifact-vs-signal
   thesis, now demonstrated across two species and three decompositions.

## 7. Freeze (provenance hashes)

The Gate-0 decision outputs are frozen for provenance.

Mouse (primary):
- combined canonical bundle sha256
  `172861417a923a83705812a96e5ad1e555055e968323289879447e53beb84396`

Human (replication):
- `v0b_v3_1_decision.json` sha256 `28301fec...`; `v0b_v3_1_per_seed.csv` sha256 `5fdc2e8c...`
- combined canonical bundle sha256
  `cc8159c882cf326bca5bbbc7f74d110d4a81338ec0ba07e208eb935bc8f729d3`

(Bundle = `cat <decision.json> <per_seed.csv> | sha256sum`, regenerable from the
committed artifacts.) This freezes a NEGATIVE-plus-method-dependence result from a
properly regularized, 10-seed, sparsity-corrected run in BOTH species - not a
positive claim from a broken run. That distinction is the whole point of the gate
(CLAUDE.md section 2). Gate 0 is now EXECUTED end to end: sparsity fixed, >=10
seeds, pre-registered v3.1, PCA/NMF baselines, marker sensitivity, and human
replication all done in-container.

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
