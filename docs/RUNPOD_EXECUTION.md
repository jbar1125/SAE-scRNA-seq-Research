# RUNPOD EXECUTION: historical causal head-to-head recipe

> **CURRENT STATUS, 2026-10-18:** this is an archived July recipe, not a verified
> end-to-end setup. Expression-space Replogle K562/RPE1 and Norman analyses later
> ran; their raw JSON archive from rented compute is missing. The scGPT/Geneformer
> embedding arm and the `embed_data(..., ...)` pseudocode below never ran or passed
> a real-data smoke test. Do not paste that stub as executable code.
>
> Current causal spec is v7 plus the training-mode erratum, SHA prefix `aae1a9d9`
> (read the complete hash in `config/causal_grounding_spec.sha256`). Prior causal
> SAE fits used **all cells**, with labels unseen; control-only training failed the
> pipeline synthetic self-test. A future matched comparison must declare training
> cells, preprocessing, gene/TF coverage, model checkpoint and code/spec versions,
> then validate both arms and archive complete outputs. The reported 6.2% from a
> different metric/panel is literature context, not a reproduction target or null.
> See `PROJECT_MASTER_SUMMARY.md` sections 8-10 and 18 for the current workflow.

The steps below record the intended expression-vs-embedding experiment on a RunPod
GPU. Only the metric and synthetic pipeline were verified when this was drafted;
these instructions do not establish a completed head-to-head.

What is VERIFIED here vs. NOT (honesty):
- Verified on CPU: the grounding metric (oracle unit test) and the full pipeline
  (SAE-learning + grounding) on realistic-scale synthetic data (6-7/8 regulators, ~0
  false positives).
- At drafting, unverified: scGPT embedding generation and real Replogle scale.
  Real expression runs were subsequently reported; the embedding step remains unverified.

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
  --n-hvg 2000 --tf-list config/tf_lists/hs_hgnc_tfs.txt \
  --out causal_out/expression_grounding.json
```
Multi-seed: repeat with `--seed 1..4` and average the grounding rate (rigor).

Two flags are REQUIRED at genome scale, both established on 2026-07-18 (see LAB_NOTEBOOK):
- `--n-hvg 2000` bounds memory. The full dense 310k x 8563 matrix (~10 GB) OOM-thrashes a
  24 GB box. This keeps the top 2000 HVGs UNION every tested perturbation gene.
- `--tf-list config/tf_lists/hs_hgnc_tfs.txt` restricts scoring to the 1,839 human TFs.
  Scoring all ~1789 perturbations includes housekeeping/essential knockdowns whose
  effect can reflect housekeeping responses. This broad TF list also includes basal
  machinery, so restriction alone does not establish regulatory specificity.
  The implemented metric is cross-fit Mann-Whitney with the amendment chain in
  current spec v7 + erratum, SHA prefix `aae1a9d9`; v6 `81b44e04` is historical.

The historical command defaults to all-cell SAE training. Record `--train-on all`
explicitly in any newly declared reproduction. For PCA/NMF ladder runs, these flags
transform the SAE input; they do not score direct PCA/NMF factors as causal baselines.

## 4. Arm A — embedding-space SAE (reproduce the field)

Generate scGPT (and/or Geneformer) cell embeddings for the SAME cells, save as an
`(n_cells, d)` `.npy` aligned to `adata` row order, then:

```bash
python3 src/causal_pipeline.py \
  --adata replogle_k562_essential.h5ad --pert-col gene --control-value non-targeting \
  --rep embedding --embedding scgpt_emb.npy --latent 2048 --k 32 --seed 0 \
  --out causal_out/embedding_grounding.json
```

UNEXECUTED scGPT pseudocode (API and output alignment require implementation/testing):
```python
# pip install scgpt; download the whole-human checkpoint per the scGPT repo.
import scgpt, numpy as np, anndata as ad
adata = ad.read_h5ad("replogle_k562_essential.h5ad")
emb = scgpt.tasks.embed_data(adata, model_dir="scGPT_human", ...)   # -> (n_cells, 512)
np.save("scgpt_emb.npy", np.asarray(emb))
```
Do not treat agreement with a literature percentage as validation. Verify cell order,
preprocessing, feature construction and controls independently, then compare the two
within-study arms under one declared metric.

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

The historical hypothesis was higher expression-space grounding. No result currently
supports that comparison. Report the matched difference and uncertainty across declared
seeds/settings, including ties or the opposite ordering. Higher aggregate rates alone
do not establish better lineage-regulatory specificity.

## 6. After the number

- The metric config and dated amendments are recorded in `config/causal_grounding_spec.json`
  (current SHA prefix `aae1a9d9`, v7 plus the 2026-10-18 training erratum; historical
  v6 prefix `81b44e04`). Read the amendment chain to distinguish designs fixed before
  particular runs from corrections made after examining results.
- If positive: build the causally-grounded gene-program atlas (a low-risk formatting of
  the grounded programs + their perturbation certificates; ELEVATION_PLAN section 4),
  then reformat the writeup and figures around the causal result. (The
  causally-supervised SAE was prototyped and did NOT validate on synthetic — shelved.)

## Compute budget (rough)

- scGPT inference on ~300k cells: ~1-3 h on an A100.
- TopK SAE training (2048 latents) per seed: minutes-to-1h on GPU.
- Grounding metric: minutes (numpy/scipy, CPU-bound).
