# COMPONENT 0 STATUS (Gate 0: lock the foundation)

Component 0 is the prerequisite gate before any Phase 2 work: pre-register one
metric, fix sparsity, finish baselines, prove marker robustness, and freeze. All
CODE for this gate is written, tested in-container, committed, and pushed. The
remaining work is EXECUTION on a GPU (Colab), which cannot run in this sandbox
(no data, Paul15 download proxy-blocked). This file states exactly what is done and
gives the single runner that finishes the gate.

## Code-complete and verified in-container

| Item | Deliverable | Verification |
|------|-------------|--------------|
| Pre-register ONE metric | `v0b_v3_1_decision.py` (Progenitor baseline, abundance-matched null, effect-size floor, >=3/5 rule) | `tests/test_v0b_v3_1_logic.py` PASS |
| Freeze the spec | `preregistration_spec.json` + `PREREGISTRATION.md`, SHA `fc342829...` | regenerable from live constants |
| PCA/NMF baselines (V11/V12) | `baselines_nmf_pca.py` (one PR formula; v3.1 modularity on PCA/NMF/SAE) | `tests/test_baselines_logic.py` PASS + synthetic run |
| Marker robustness | `marker_sensitivity.py` (leave-one-marker-out) | synthetic run, 0 verdict flips |
| Overcomplete + marker-aware training | `train_sae.py`, `preprocess_paul15.py` (log1p, corrected globins/Alas2/Ermap) | synthetic train->v3 verified |

Every logic test passes: v0b, v3, v3.1, baselines. All scripts compile.

## Needs a GPU run (the actual Gate 0 execution)

These require the real data + GPU and are the only things left for Gate 0:
1. **Sparsity fix**: L1 sweep so mean L0 lands in 20-50, then retrain **>=10 seeds**.
2. **Run the frozen v3.1** on the sparsity-corrected checkpoints (mouse).
3. **Run baselines** (PCA/NMF PR + modularity) and **marker sensitivity** on the same.
4. **Human replication**: repeat on the 25k human SAEs.
5. Only THEN freeze the module-assignment output table SHA-256.

## The runner (Colab, GPU) — finishes Gate 0

Step A — sparsity sweep (fast, 1 seed each; pick the L1 whose L0 is in 20-50):
```python
import os
if not os.path.isdir('/content/drive/MyDrive'): 
    from google.colab import drive; drive.mount('/content/drive')
!pip -q install scanpy statsmodels scikit-learn
SHA="15d00e4"
for f in ["v0b_module_definitions","v0b_v3_loading","v0b_v3_1_decision","train_sae","preprocess_paul15","baselines_nmf_pca","marker_sensitivity"]:
    !rm -f {f}.py*; !wget -qO {f}.py https://raw.githubusercontent.com/jbar1125/SAE-scRNA-seq-Research/{SHA}/{f}.py
!python preprocess_paul15.py --n-hvg 2000 --out-dir /content/drive/MyDrive/data_g0
for l1 in [0.4, 0.6, 0.8, 1.0]:
    !python train_sae.py --matrix /content/drive/MyDrive/data_g0/expression_matrix.npy --latent-dim 512 --l1 {l1} --seeds 1 --out-dir /tmp/sweep_{l1}
```
Read each print: pick the `--l1` whose L0 is closest to the middle of 20-50.

Step B — final train (>=10 seeds at the chosen L1) + the full frozen analysis:
```python
L1 = 0.6   # <- set to the value chosen in Step A
!python train_sae.py --matrix /content/drive/MyDrive/data_g0/expression_matrix.npy --latent-dim 512 --l1 {L1} --seeds 10 --out-dir /content/drive/MyDrive/data_g0
!python v0b_v3_1_decision.py --data-dir /content/drive/MyDrive/data_g0 --output-dir /content/drive/MyDrive/data_g0/v0b_outputs
!python baselines_nmf_pca.py --data-dir /content/drive/MyDrive/data_g0 --rank 512 --out /content/drive/MyDrive/data_g0/baselines
!python marker_sensitivity.py --data-dir /content/drive/MyDrive/data_g0 --out /content/drive/MyDrive/data_g0/sensitivity
```
Note: the loaders auto-detect all `sae_seed*.pt` files and set the seed count
accordingly (and the decision rule scales as a >=60% majority), so `--seeds 10`
works with no code edits.

## After the run
- Confirm L0 in 20-50 (sparse SAE), controls non-significant (calibrated null).
- Record the v3.1 verdict + the PCA/NMF-vs-SAE modularity comparison.
- Run the same on human; require the direction to hold mouse AND human.
- Freeze `module_assignments`/decision output SHA-256 into the repo. Gate 0 done.

## Honest note
The pre-registration is frozen at the metric/marker level (SHA `fc342829`). The
scientific verdict is NOT yet locked because it must come from a sparsity-corrected,
>=10-seed, human-replicated run — freezing an under-regularized run would repeat the
exact mistake this gate exists to prevent.
