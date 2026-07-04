# PRE-REGISTRATION (Component 0, Gate 0)

This freezes the single analysis specification BEFORE the Phase 1 finalization
final runs and BEFORE any Phase 2 data inspection. It exists so no metric, marker
set, threshold, or decision rule can be chosen after seeing results. The frozen
verdict history (v1 NO / v3 YES / metric-dependent) is exactly why one metric must
be locked; this is that lock.

## Frozen spec + SHA-256

- Machine-readable spec: [`preregistration_spec.json`](./preregistration_spec.json)
- Frozen SHA-256: **`fc342829191da5a0bb5244975c3f61c8934943cb512818446e5b5f0984f6011a`**
  (also in `preregistration_spec.sha256`)

The spec is generated from the live code constants, so it cannot silently drift
from the implementation. To verify the freeze is intact, regenerate and compare:

```bash
python3 - <<'PY'
import json, hashlib, v0b_module_definitions as v0b, v0b_v3_loading as v3, v0b_v3_1_decision as v31
# (regenerate exactly as in the freeze commit) ...
PY
sha256sum preregistration_spec.json   # must match the SHA above
```

If any constant, marker set, or threshold changes, the SHA changes. After this
freeze that is a NEW pre-registration with an explicit, dated amendment note in
this file — never a silent edit.

## What is frozen (summary; JSON is authoritative)

**Primary modularity metric** — `v0b_v3_1_decision`:
- Continuous decoder-loading strength (no TOP_K, no significance gate at the
  loading step): per (seed, feature, submodule) enrichment = mass on markers /
  total feature mass, size-normalized; submodule strength = mean of top-3 features.
- Control baseline = **Progenitor only** (Cycling excluded from the baseline; it is
  a real program). Reported as a secondary control.
- **Abundance-matched permutation null** (per-gene decoder-mass-matched random gene
  sets, 1000 draws, BH-FDR). A submodule is a REAL program iff excess over the
  Progenitor baseline > 0.10 AND permutation q < 0.05.
- Decision: asymmetric modularity (erythroid more distributed than granulocyte) is
  SUPPORTED iff in >= 3 of 5 seeds: ery_n_real >= 2 AND ery_n_real >= gran_n_real
  AND ery excess-entropy > gran. The plain program-count finding is reported either
  way.

**Marker sets** — frozen in the JSON, including the corrected Paul15 symbols
(`Hba-a2`, `Hbb-b1`, `Alas2`, `Ermap`) that fixed the earlier curation artifact.

**Phase 2 causal test** — head-to-head against the Kendiukhov (2026) 6.2%
foundation-model null: the identical Replogle CRISPRi validation is run on
expression-space AND embedding-space SAE features; one-sided binomial vs 0.062;
Gate 2 thresholds frozen in the JSON.

**Statistics** — BH-FDR per family; cross-seed rule >= 3/5 mouse or >= 2/3 human;
expression-matrix MD5 verified at load.

## What is NOT frozen (and why it is fine)

- Sparsity hyperparameter (L1) and dictionary size are tuned to hit the L0 20-50 QC
  band and reported as a sweep; the metric is applied identically at every setting.
- The number of seeds may increase (>= 10 planned); the >= 3/5 rule scales to a
  >= 60% majority.
- Datasets added in Phase 2 are analyzed with this same frozen spec.

## Freeze status

- Metric, marker sets, thresholds, decision rule: FROZEN at the SHA above.
- Module-assignment output table SHA-256: to be frozen after the finalized run on
  the sparsity-corrected, >=10-seed checkpoints (do NOT freeze an under-regularized
  or human-unreplicated run).
