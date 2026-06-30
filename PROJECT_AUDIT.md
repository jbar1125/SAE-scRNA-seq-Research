# PROJECT AUDIT: consistency + rigor

Authoritative correction and rigor layer over `PROJECT_HANDOFF.md`. The handoff is
a Phase-1-closeout snapshot and is left intact as a historical record; this file
carries every inconsistency found, its resolution status, and the rigor upgrades.
Read both.

Status keys: RESOLVED (fixed in this repo), OPEN (needs data/run), VERIFY (needs
source check before any submission), FLAG (needs a decision/recompute).

---

## A. Inconsistencies found

| # | Location | Issue | Severity | Status |
|---|----------|-------|----------|--------|
| A1 | Handoff 5 (V11), 7 | PCA participation ratio given as both **43.4** and **42.3** | med | FLAG: recompute once with one documented formula; lock the value. No code/data in repo to recompute here. |
| A2 | Handoff 7 | NMF participation ratio: "**1.8x higher** than SAE" vs "**~34.94**". With SAE=24.2, 34.94/24.2 = **1.44x**, not 1.8x (1.8x would be ~43.6). Internal arithmetic conflict. | med | FLAG: compute V11 NMF once, report one number + consistent multiplier. |
| A3 | Handoff 13 | Kendiukhov "**2025**, arXiv **2603.02952**". An arXiv id `2603.xxxxx` is March **2026**, not 2025. Year vs id-date conflict. Same for the two 2026 Kendiukhov ids (2603.01752, 2603.11940) sitting implausibly close to the "2025" one. | med | VERIFY: confirm real ids + years before citing. The 6.2% null depends on this paper; get it exactly right. |
| A4 | Handoff 13 | Beneyto-Calabuig 2023, Lasry 2023, Petti 2019 not confirmed (volume/pages). | med | VERIFY before submission (already flagged at handoff). |
| A5 | Drive | The canonically-named, latest `v0b_module_definitions.ipynb` is **V1**, not v2 (see B). The handoff assumes the user runs v2. | HIGH | RESOLVED: `v0b_module_definitions.py` is the canonical v2 and supersedes all three notebooks. |
| A6 | Drive `v0b_outputs/` | No `module_assignments_v0b.csv` exists anywhere. V0b has **never been completed/frozen**. The linchpin finding is unvalidated. | HIGH | OPEN: run the new script, freeze the SHA-256. |
| A7 | Handoff 4.3 | Phase-1 SAE training config only partially pinned: batch 256 / wd 1e-5 / clip 1.0 / cosine+warmup are attributed to the *Phase-2 embedding-SAE* config, leaving the exact Phase-1 schedule unspecified. | med | OPEN (V3): extract exact Phase-1 hyperparameters from the training notebook. |
| A8 | Handoff 4.2 | CELLxGENE Census release date not pinned. | med | OPEN: record census version + date (handoff action #2). |
| A9 | Handoff 9 (Part 1), 4.6 | The "committed_granulocyte" group is `prob_14Mo + prob_16Neu + prob_13Baso`. Monocytes (14Mo/15Mo) are not granulocytes; basophils are granulocytes. So the group is really "committed myeloid (mono+neu+baso)". The "granulocyte unified module" claim inherits this loose label. | low-med | FLAG: either rename to "myeloid" or justify the grouping explicitly in text. The v2 script keeps the handoff definition but the label is imprecise. |
| A10 | Handoff 7 | Annotation table GATA1/TAL1 absorbed into GFI1B (overlap + short-marker normalization bias). | med | OPEN: regenerate `generate_annotation_table.py` with the curated non-overlapping v2 marker sets; refreeze SHA-256. |

### Cross-checks that PASSED (no inconsistency)
- Recon MSE table (4.4): mean 0.6157 and std 0.0055 are arithmetically correct for the five per-seed values.
- 19 paul15 clusters (4.1) enumerate to exactly 19.
- V0b v2 predicted cell counts sum to 2730 (751+897+700+382), and the v1 variant (751+897+973+109) also sums to 2730. Internally consistent.
- Test-family sizes: 896 = 128x7 and 3072 = 24x128 both check out.
- 640 annotation entries = 128x5. Checks out.

---

## B. The V0b notebook trap (evidence)

The most-recently-modified `v0b_module_definitions.ipynb` (Drive id
`1Jg1fS80...`, 2026-05-29), which also carries the canonical name, fails all of
the handoff's own v2 verification checks:

| Handoff v2 marker | That notebook | Verdict |
|---|---|---|
| `WINNER_TAKE_ALL_MIN_OVERLAP = 2` | `CROSS_SEED_THRESHOLD = 3` | v1 |
| `Ery_TF` = Gata1, Klf1, Tal1, Lmo2, Zfpm1, Stat5a, Bcl11a, Myb (HVG-restricted) | adds Nfe2/Gata2/Mafg/Sox6/Foxo3 + hemoglobin (`Hba-a1`...) and late-granulocyte (`Ltf`, `S100a8`...) submodules not in the HVG set | v1 |
| `terminal_cols` includes `prob_19Lymph` | `['prob_1Ery','prob_14Mo','prob_16Neu','prob_13Baso','prob_11DC','prob_8Mk']` (no `prob_19Lymph`) | v1 |
| metric = normalized entropy + largest-submodule fraction, per seed | metric = `H_ery / H_gran` entropy ratio | v1 |
| output `module_assignments_v0b.csv` | output `module_assignments.csv` | v1 |

Consequence: anyone opening the obvious notebook runs v1, reproducing every bug
v2 was built to fix. The repo's `v0b_module_definitions.py` removes this hazard:
one file, in version control, reviewable, unit-tested.

---

## C. Methodological rigor upgrades

1. **Small-marker-set FDR floor (new finding).** `Gran_TF` has 3 HVG markers and
   `Gran_Primary` 4. Their *best-case* hypergeometric p-values at full overlap are
   ~9e-4 and ~8e-5. Under BH-FDR across ~896 tests/seed, a granulocyte submodule
   only survives q<0.05 when many other tests are co-significant (favorable rank).
   If few features genuinely load on granulocyte programs, granulocyte will look
   under-assigned for a *statistical* reason (tiny marker sets + multiple testing),
   confounding the biological asymmetry signal. This was surfaced by the synthetic
   smoke test. Mitigations to add: report each submodule's minimum achievable q;
   consider per-lineage FDR families; or size-match marker sets. The asymmetry
   verdict must be interpreted against this floor, not in isolation.
2. **n=5 bootstrap is weak.** The asymmetry CI is over 5 seed-level values. The
   script reports it but makes the **>=3/5 per-seed sign-agreement** count the
   primary criterion, and states the caveat in output. Add a label-permutation
   null (shuffle submodule labels within seed) for a calibrated p-value.
3. **Winner-take-all hides co-membership.** Forcing each feature to one submodule
   is conservative for a "distribution" claim but masks features enriched for two
   erythroid submodules. Add a soft/multi-membership robustness pass.
4. **TOP_K sensitivity.** Assignments use top-30 decoder genes. V4's K sweep
   (10/20/30/50) should be run to show the asymmetry verdict is not K-specific.
5. **Granulocyte label precision.** See A9.
6. **Perturbation still circular / buggy.** V14 (OOD self-distance=0), V15
   (decoder/classifier scale), V0d (non-circular marker scoring) remain unaddressed.
   Code is not in the repo. These do not affect V0b but block the perturbation arm.

---

## D. Per V-item rigor status (what "highest rigor" needs)

- **V0a** DONE. Optional: complete the 10-root sensitivity sweep (V7).
- **V0b** Canonical v2 script delivered + logic-tested. OPEN: run on real data, freeze SHA-256, decide the asymmetry verdict honestly.
- **V0c** bimodality (Hartigan dip + BH-FDR) NOT BUILT. Replace the heuristic 49-feature count.
- **V0d / V14 / V15** perturbation NOT BUILT/UNFIXED. See C6.
- **V3** training reproducibility incomplete (A7): pin exact Phase-1 config.
- **V4** marker recovery K-sweep PENDING (C4).
- **V5** cross-seed feature stability NOT DONE (distinct from V0b).
- **V9** logistic trajectory fits NOT DONE: replace the provisional 0.37 crossover.
- **V11/V12** participation ratio + NMF modularity PARTIAL: resolves A1, A2.
- **V13** human concordance PARTIAL: rerun v2-style asymmetry on the 3 human seeds.

---

## E. Reproducibility gaps and what this commit closes

Closed here:
- Canonical, version-controlled V0b code (`v0b_module_definitions.py`).
- Pinned environment (`requirements.txt`, validated by the test run).
- Input manifest with the expression-matrix MD5 and Drive locations (`data/README.md`).
- Logic regression tests (`tests/test_v0b_logic.py`).
- `.gitignore` that still allows committing the small frozen pre-registration artifacts.

Still open:
- Pin CELLxGENE Census version/date (A8).
- Pin exact Phase-1 SAE training config (A7).
- After running V0b: commit `module_assignments_v0b.csv`, `modularity_metrics_per_seed.csv`, `v0b_provenance.json` (the SHA-256 freeze).
- Verify citations A3, A4.

---

## F. Honest bottom line

The headline ("asymmetric modularity") has not survived its own gene-content test,
because that test was never run to a frozen result, and the notebook most likely
to be run is the superseded v1. Until `v0b_module_definitions.py` runs on the real
data and the verdict is recorded, every Phase-2 claim that depends on it is
provisional. The script is built to report a null honestly if the asymmetry does
not hold.

---

## G. V0b first run + diagnostic (2026-06-30)

Ran `v0b_module_definitions.py` on the real data. MD5 verified; cell-group counts
(committed_granulocyte 897, committed_erythroid 751, uncommitted 697,
intermediate 385) match the handoff's v2 predictions; all v2 markers fully present
in the HVG set. Verdict: **asymmetric modularity NOT SUPPORTED**, but the result is
confounded and the metric was uninformative (erythroid 0-1 features assigned per
seed). module_assignments_v0b.csv SHA-256 `f18b088a95e8e4f3248561d72a711f2bd6d1302127afd60b753a4e9d211d96b5`. **NOT frozen as pre-registration** (inconclusive run).

Per-submodule enrichment diagnostic (640 tests = 5 seeds x 128 features each):

| submodule | max overlap | n(overlap>=2) | n(sig, q<.05) | min q |
|-----------|-------------|---------------|---------------|-------|
| Ery_TF | 2 | 7 | 0 | 0.304 |
| Ery_Heme | 3 | 9 | **2** | 0.006 |
| Ery_Membrane | 2 | 6 | 0 | 0.188 |
| Gran_TF | 2 | 2 | 0 | 0.053 |
| Gran_Primary | 4 | 43 | **29** | 0.000 |
| Progenitor | 2 | 1 | 0 | 0.096 |
| Cycling | 4 | 40 | 3 | 0.006 |

Findings:
1. The asymmetry signal is essentially one submodule, **Gran_Primary** (Mpo/Elane/Prtn3/Ctsg). "Granulocyte unified" reduces to "the primary-granule effector co-occurs."
2. **Both** TF programs are undetected (Gran_TF and Ery_TF: 0 significant). The granulocyte TF program is as invisible as the erythroid one, so the lineage asymmetry as stated is absent.
3. Real biological axis is **effector-concentration vs TF-dilution** (Gran_Primary + Ery_Heme concentrate; all TF/membrane sets dilute), orthogonal to granulocyte-vs-erythroid.
4. The test is partly underpowered: Ery_TF/Ery_Membrane/Gran_TF have real overlap-2 co-occurrence killed by BH across ~896 mostly-null tests. Gran_TF min_q=0.053 just misses; dropping controls from the FDR family or per-lineage FDR may change it.
5. The defensible remaining finding is smaller: primary-granule and heme-synthesis effector modules are recovered as concentrated SAE features.

Required redesign (V0b v3): replace thresholded winner-take-all with a continuous
per-feature submodule loading score (e.g., summed |decoder weight| on present
markers, or a rank-enrichment score), measure distribution on the continuous
scores with no significance gate, run the TOP_K sweep (10/20/30/50) and per-lineage
FDR as sensitivity checks, and reframe the hypothesis around effector-vs-TF
concentration. Treat the v1 "supported" result as not replicated under rigorous
markers.
