# REPRODUCIBILITY

How to re-verify this project end to end, and the checks that pass as of the latest
commit. Reproducibility is a first-class deliverable here (and an ISEF/STS judging
criterion), so the verification is scripted, not asserted.

## Environment

- CPU only (no GPU needed for Gate 0 + hardening). ~4 cores, 16 GB RAM is ample.
- Python 3.11; `pip install -r requirements.txt` (numpy, pandas, scipy, statsmodels,
  scikit-learn, torch, scanpy/anndata/h5py for preprocessing, matplotlib for figures).
- Large inputs (matrices, checkpoints, raw h5/h5ad) are gitignored and regenerable;
  the committed evidence is the small JSON/CSV artifacts + figures.

## Verification checks (all PASS at the latest commit)

| # | Check | Command | Expected |
|---|-------|---------|----------|
| 1 | Syntax | `python3 -m py_compile src/*.py tests/*.py` | no error |
| 2 | Logic tests | `python3 tests/test_v0b_logic.py` (and v3, v3_1, baselines) | `ALL ... TESTS PASSED` |
| 3 | Pre-registration freeze | `sha256sum config/preregistration_spec.json` | `fc342829191da5a0bb5244975c3f61c8934943cb512818446e5b5f0984f6011a` (matches `config/preregistration_spec.sha256`) |
| 4 | Frozen mouse decision bundle | `cat data_g0/v0b_outputs/v0b_v3_1_decision.json data_g0/v0b_outputs/v0b_v3_1_per_seed.csv \| sha256sum` | `172861417a923a83705812a96e5ad1e555055e968323289879447e53beb84396` |
| 5 | Frozen human decision bundle | same on `data_g0_human/...` | `cc8159c882cf326bca5bbbc7f74d110d4a81338ec0ba07e208eb935bc8f729d3` |
| 6 | Figures deterministic | `python3 src/make_figures.py` twice | byte-identical PNG/SVG both runs |

## Re-run the analysis from checkpoints (fast, no training)

```bash
# frozen verdict (mouse); add V0B_MARKERS_JSON=config/human_markers.json for human
python3 src/v0b_v3_1_decision.py --data-dir ./data_g0 --output-dir /tmp/check
# baselines, sensitivity, robustness, cross-refs, extra families
python3 src/baselines_nmf_pca.py  --data-dir ./data_g0 --rank 512 --seeds 10 --out /tmp/bl
python3 src/marker_sensitivity.py --data-dir ./data_g0 --out /tmp/ms
python3 src/robustness_sweeps.py  --data-dir ./data_g0 --out /tmp/r1
python3 src/grn_baseline.py --data-dir ./data_g0 --tfs config/tf_lists/mm_mgi_tfs.txt --out /tmp/grn
python3 src/ica_baseline.py --data-dir ./data_g0 --out /tmp/ica     # rank 50 (converges)
python3 src/make_figures.py
```

## Re-run from raw data (retrains; ~10 min/species on CPU)

```bash
python3 src/preprocess_paul15.py --n-hvg 2000 --out-dir ./data_g0
python3 src/train_sae.py --matrix ./data_g0/expression_matrix.npy --latent-dim 512 \
    --l1 0.8 --seeds 10 --epochs 300 --out-dir ./data_g0
# then the analysis block above
```

## Fixed-on-audit note (self-correction)

The pre-registration hash `fc342829` is the SHA-256 of the canonical serialization
`json.dumps(spec, indent=2)`. The committed `preregistration_spec.json` had originally
been written with a trailing newline, so `sha256sum` of the file returned a different
value (`6398c5ed…`) and the documented verification (check #3) did not reproduce the
recorded hash. Caught during a reproducibility pass and fixed by normalizing the file
to the exact canonical bytes; the spec CONTENT is byte-identical to the freeze commit
(`35173e4`) — verified `json.load` equality — so nothing that was frozen changed, only
the trailing byte, and check #3 now passes. This is logged in `LAB_NOTEBOOK.md`.
