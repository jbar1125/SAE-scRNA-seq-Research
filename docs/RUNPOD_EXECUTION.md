# RUNPOD EXECUTION — the causal head-to-head (Phase 1)

Turnkey steps to run the expression-vs-embedding causal-grounding head-to-head on a
RunPod GPU. The metric and pipeline are already built and CPU-verified
(`src/causal_grounding.py`, `src/causal_pipeline.py`); this is the data + GPU part.

What is VERIFIED here vs. NOT (honesty):
- Verified on CPU: the grounding metric (oracle unit test) and the full pipeline
  (SAE-learning + grounding) on realistic-scale synthetic data (6-7/8 regulators, ~0
  false positives).
- Not verifiable here: the scGPT embedding step (needs the model + GPU) and the real
  Replogle scale. Those run on RunPod.

## 0. RunPod instance

- GPU: 1x A100 40/80GB or L40S (48GB). 16-32 vCPU, 100+ GB disk.
- Image: a PyTorch 2.x CUDA image (RunPod "PyTorch" template is fine).
- Clone: `git clone https://github.com/jbar1125/SAE-scRNA-seq-Research && cd SAE-scRNA-seq-Research`

## 1. Environment

```bash
pip install -r requirements.txt
pip install pertpy scanpy anndata scikit-learn statsmodels
# confirm the built pieces work on this box (GPU auto-detected):
python3 tests/test_causal_grounding.py          # ALL CAUSAL-GROUNDING TESTS PASSED
python3 src/causal_pipeline.py --synthetic       # PIPELINE SELF-TEST PASSED (now on GPU)
```

## 2. Get the Replogle CRISPRi data

```python
import pertpy
adata = pertpy.data.replogle_2022_k562_essential()   # ~300k cells, ~2k perturbations
adata.write("replogle_k562_essential.h5ad")
# obs perturbation column is typically 'gene' (or 'gene_id'); control is 'non-targeting'
print(adata.obs.columns.tolist()); print(adata.obs['gene'].value_counts().head())
```
Note the exact perturbation-column name and the control label; pass them below.

## 3. Arm B — expression-space SAE (the contribution)

```bash
python3 src/causal_pipeline.py \
  --adata replogle_k562_essential.h5ad --pert-col gene --control-value non-targeting \
  --rep expression --latent 2048 --k 32 --seed 0 --n-shuffle 50 \
  --out causal_out/expression_grounding.json
```
Multi-seed: repeat with `--seed 1..4` and average the grounding rate (rigor).

## 4. Arm A — embedding-space SAE (reproduce the field)

Generate scGPT (and/or Geneformer) cell embeddings for the SAME cells, save as an
`(n_cells, d)` `.npy` aligned to `adata` row order, then:

```bash
python3 src/causal_pipeline.py \
  --adata replogle_k562_essential.h5ad --pert-col gene --control-value non-targeting \
  --rep embedding --embedding scgpt_emb.npy --latent 2048 --k 32 --seed 0 \
  --out causal_out/embedding_grounding.json
```

scGPT embedding (the one step needing hands-on setup):
```python
# pip install scgpt; download the whole-human checkpoint per the scGPT repo.
import scgpt, numpy as np, anndata as ad
adata = ad.read_h5ad("replogle_k562_essential.h5ad")
emb = scgpt.tasks.embed_data(adata, model_dir="scGPT_human", ...)   # -> (n_cells, 512)
np.save("scgpt_emb.npy", np.asarray(emb))
```
Target: reproduce the ~6-10% embedding ceiling (validates the pipeline against the
published number). If it lands there, the comparison is trustworthy.

## 5. Compare

```bash
python3 - <<'PY'
import json
e=json.load(open("causal_out/expression_grounding.json"))
a=json.load(open("causal_out/embedding_grounding.json"))
print("expression-space:", e["causal_grounding_rate"], f"({e['n_grounded']}/{e['n_tested']})")
print("embedding-space :", a["causal_grounding_rate"], f"({a['n_grounded']}/{a['n_tested']})")
print("fold improvement:", e["causal_grounding_rate"]/max(a["causal_grounding_rate"],1e-9))
PY
```

**Success = expression-space grounding is clearly and robustly higher** (target: a
multiple of the ~6-10% embedding ceiling), across seeds. That is the positive,
benchmark-beating result. If it is NOT higher, report it honestly (the mechanistic
prior favors a positive result, but the design is internally fair either way).

## 6. After the number

- Pre-register the exact metric config before the final runs (freeze a spec SHA, as in
  Gate 0). A draft spec lives in `config/causal_grounding_spec.json`.
- If positive: build the causally-grounded gene-program atlas and/or the
  causally-supervised SAE (ELEVATION_PLAN sections 4), then reformat the writeup and
  figures around the causal result.

## Compute budget (rough)

- scGPT inference on ~300k cells: ~1-3 h on an A100.
- TopK SAE training (2048 latents) per seed: minutes-to-1h on GPU.
- Grounding metric: minutes (numpy/scipy, CPU-bound).
