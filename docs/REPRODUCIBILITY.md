# Reproducibility

The committed Gate-0 decision bundles and structural summaries can be inspected
and integrity-checked without downloading large data. Reproducing their model
fits additionally requires source data and checkpoints. Original Component 2
result JSONs were lost with ephemeral GPU instances; narrative tables remain,
but exact causal reproduction has not been demonstrated.

## Validated cloud environment

The 2026-10-18 cloud setup uses Python 3.12 in
`/workspace/.venvs/sae-scrna`. Activate it from the repository root:

```bash
source /workspace/.venvs/sae-scrna/bin/activate
export XDG_CACHE_HOME=/workspace/.cache
export MPLCONFIGDIR=/workspace/.cache/matplotlib
export NUMBA_CACHE_DIR=/workspace/.cache/numba
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
```

The reusable install/validation instructions and resolved package list live under
`/workspace/onboarding/sae-scrna/` in this cloud workspace. They are environment
artifacts, not a portable repository lockfile. `requirements.txt` has historical
core pins plus open-ended Torch and Scanpy requirements; do not describe it as a
complete, frozen environment. Record resolved versions, device, CPU thread limits,
and source revision with each new run.

The continuation also archives the exact resolved instance package list and
principal versions under
[`data_program_modulation_audit/environment/`](../data_program_modulation_audit/environment/).
This records the Python 3.12 Linux CPU environment; it does not assert portability
to arbitrary CUDA or operating-system configurations.

Gate-0 analysis and synthetic checks work on CPU. Full Perturb-seq fits can require
substantially more time and memory. Hardware availability and measured run duration
must be checked for the actual run rather than inferred from an old estimate.

## Code checks

The onboarding pass ran package consistency, syntax, all nine existing standalone
test suites, specification hashes, and frozen decision-bundle hashes. These checks
establish the environment and tested code paths, not replication of historical
real-data training or a new discovery.

After the real-data CPU continuation, all **15** standalone suites passed with
zero skips, including counterexample, conditional-identity, checkpoint, duplicate-
preparation, marker-baseline and sensitivity suites. Actual small SAE fitting and the
updated AnnData CLI path executed. The current CLI smoke fixture yielded no
significant named-program matches; it verifies safe unscored behavior, not a
biological result. A later full-data marker baseline and its explicit sensitivities
completed and are archived in [the continuation record](NORMAN_CONTINUATION_RESULTS.md);
full-size SAE training is still running. See [the audit](PROGRAM_MODULATION_AUDIT.md) for the separately
archived 80 design-counterexample runs and 200 deployed-default null simulations.

```bash
python -m pip check
python -m py_compile src/*.py tests/*.py
for test_file in tests/test_*.py; do
  python "$test_file" || exit $?
done
```

`test_pipeline_training.py` trains actual small SAEs. It records that its current
all-cell mode recovers 8/8 planted regulators, while control-only mode recovers
0/8. The latter is a documented power limitation of this synthetic setup.

## Verify preserved evidence without modifying it

Run this read-only check from the repository root:

```bash
python - <<'PY'
import hashlib
from pathlib import Path

for spec in sorted(Path("config").glob("*_spec.json")):
    expected = spec.with_suffix(".sha256").read_text().split()[0]
    actual = hashlib.sha256(spec.read_bytes()).hexdigest()
    assert actual == expected, (str(spec), actual, expected)
    print("SPEC PASS", spec, actual)

bundles = {
    "data_g0": "172861417a923a83705812a96e5ad1e555055e968323289879447e53beb84396",
    "data_g0_human": "cc8159c882cf326bca5bbbc7f74d110d4a81338ec0ba07e208eb935bc8f729d3",
}
for dataset, expected in bundles.items():
    run = Path(dataset) / "v0b_outputs"
    payload = ((run / "v0b_v3_1_decision.json").read_bytes()
               + (run / "v0b_v3_1_per_seed.csv").read_bytes())
    actual = hashlib.sha256(payload).hexdigest()
    assert actual == expected, (dataset, actual, expected)
    print("BUNDLE PASS", dataset, actual)
PY
```

The Gate-0 metric spec is `fc342829...`. Current causal, regulon, and modulation
specifications have separate sidecar hashes. The causal spec includes dated
outcome-informed amendments and an implementation erratum. Matching a hash
establishes integrity of those bytes, not confirmatory status or result validity.

The mouse decision JSON's prose rounds erythroid count to 2, while its numeric
field and per-seed CSV give median 1.5. Cite the numeric evidence. Preserve the
frozen bytes and disclose this wording discrepancy rather than changing the bundle.

## Preserve frozen outputs during regeneration

Do not run training or preprocessing into `data_g0/` or `data_g0_human/` during
development. Use a new run directory and compare its outputs with the reference
artifacts. The committed directories omit large expression matrices, cell metadata,
and checkpoints, so a fresh checkout alone cannot rerun checkpoint analyses.

The supported entrypoints are `preprocess_paul15.py`,
`preprocess_human_marrow.py`, `train_sae.py`, and `v0b_v3_1_decision.py`. Their CLI
help exposes input/output paths. Mouse preprocessing follows its depth-adjusted
source and log policy; human preprocessing requires the Setty raw-count file and
normalizes counts before `log1p`. Their frozen matrix MD5s are:

| Dataset | Expression matrix MD5 |
|---|---|
| Mouse | `184ea5ddfca23db2bbb343b311688c45` |
| Human | `627961897d52a50c7a81f5bbf7e04097` |

Human structural analyses require
`V0B_MARKERS_JSON=config/human_markers.json`. Use each new run's provenance report
for newly generated matrices rather than forcing them to match by editing a hash.
Differences across numerical libraries and devices must be investigated and
reported; seed equality alone does not guarantee byte-identical training outputs.

`src/make_figures.py` reads committed structural summaries and normally writes to
`docs/figures/`. When verifying deterministic rendering, redirect the module's
`OUT` path to a scratch location or use a disposable copy of the checkout so
rendering cannot overwrite tracked figures. A figure reproduced from committed
summaries verifies rendering, not model fitting. The historical causal diagnostic
was built from hardcoded narrative numbers and originally contained a wrong-
background essentiality panel; it is not an independent evidence source.

## Missing causal evidence and publication requirements

The historical Replogle K562/RPE1, effect-size, representation-ladder, Norman,
and regulon-recovery outputs were not committed and are reported as lost in
[the master summary, section 12](PROJECT_MASTER_SUMMARY.md#12-frozen-artifacts-registry).
The current implementation includes later changes; running current defaults is
not automatically a recreation of every historical metric version.

For each causal rerun, retain:

- Code SHA, full command, data source and checksum, preprocessing decisions, and
  the resolved environment.
- Actual SAE training population and held-out design. Historical SAEs used all
  cells. Cross-fitting did not hold test cells out of SAE training.
- Seed, architecture, training diagnostics, per-perturbation results, null outputs,
  uncertainty estimates, and all tested programs including negative/low-power ones.
- Machine-readable results and a manifest in a durable archive; large datasets
  and checkpoints can be deposited separately. Link the archive in the manuscript.

The PCA/NMF ladder trains an SAE on PCA/NMF inputs; plain PCA/NMF feature baselines
remain a separate experiment. A literature grounding rate from another metric is
not a substitute for the unrun matched embedding arm. Program-modulation inference
must address selection, multiplicity, perturbation identity, and the dominant-axis
alternative before any positive is treated as a discovery.

## Historical hash correction

The Gate-0 spec was intended to use the canonical `json.dumps(spec, indent=2)`
serialization. A trailing newline initially caused the file hash to disagree with
the documented `fc342829...` hash. The file was normalized to the canonical bytes;
parsed content was checked against freeze commit `35173e4`. This integrity repair
is documented in `LAB_NOTEBOOK.md` and did not change the analysis specification.
