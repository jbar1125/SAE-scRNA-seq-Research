# SAE / scRNA-seq Hematopoiesis Research

A computational research project by Jacob Barzideh evaluating whether sparse
autoencoder (SAE) features provide reliable descriptions of hematopoietic gene
programs, and whether those descriptions survive structural and perturbation tests.

**Current conclusion:** the measured modularity of erythroid and granulocyte programs
depends on the decomposition. The original asymmetric-modularity hypothesis is not
supported by the frozen SAE analyses in either species. Later perturbation analyses
have not established a new lineage-regulatory discovery or an advantage over
foundation-model features.

Start with [the master summary](docs/PROJECT_MASTER_SUMMARY.md),
[the results](docs/COMPONENT0_RESULTS.md), and
[the reproducibility guide](docs/REPRODUCIBILITY.md). Agents must read
[CLAUDE.md](CLAUDE.md) before working. The master summary records the project through
2026-10-18; its section 17 documents additional historical claims requiring correction.

## Evidence and status

| Component | Status | Evidence and limits |
|---|---|---|
| Gate 0, mouse Paul15 | Executed; frozen decision bundle | 2,730 cells, 2,012 genes, 10 SAE seeds; original claim supported in 0/10 seeds. [Results](docs/COMPONENT0_RESULTS.md). |
| Gate 0, human Setty 2019 CD34+ marrow | Executed; frozen decision bundle | 4,142 cells, 2,032 genes, 10 SAE seeds; original claim supported in 0/10 seeds. |
| Structural comparison | Executed, committed summaries | PCA, NMF, SAE, co-expression GRN, and ICA across two datasets; only mouse NMF supports the original decision rule. ICA uses rank 50 rather than 512. [Hardening](docs/COMPONENT0_HARDENING.md). |
| Structural robustness | Executed, committed summaries | Threshold, sparsity, seed-count, and gene-panel sweeps; marker leave-one-out flips 0/34 mouse and 0/36 human decisions. |
| Replogle CRISPRi, K562 and RPE1 | Historical exploratory results; original JSONs missing | Recorded perturbation-responsive programs are mainly housekeeping/chromatin programs. Original outputs from ephemeral GPU instances must be regenerated and archived. [Causal results](docs/COMPONENT2_RESULTS.md). |
| Causal representation ladder | Historical exploratory results | Every arm trains an SAE: on expression, PCA-reduced input, or NMF-reduced input. This is not a direct SAE-versus-PCA/NMF comparison. |
| Norman CRISPRa grounding and regulon recovery | Historical exploratory results; original JSONs missing | Grounding rate did not exceed the recorded perturbation-identity shuffle null; single-TF regulon recovery was unsuccessful. Neither result establishes absence of all biological signal. |
| Original nine-panel modulation | Three full-size CPU seeds and one control-only sensitivity completed | Three panels match in all primary seeds; KLF1 is the sole repeated legacy candidate for these panels and fails the old background gate. [Completed fits](docs/NORMAN_COMPLETED_MODELS.md). |
| Historical K562 program | Complete exploratory follow-up | 63/61/61 all-cell SAE associations, 57 direct-score calls; 38 survive reference, depth and every single-gene omission. PCA also recovers the program. [Results](docs/NORMAN_K562_CONTEXT.md). |
| Expression versus foundation-model embedding | Not run | The published 6.2% figure uses a different protocol and is not a matched baseline for this project. |

The Gate-0 metric was version-locked before the final runs, following exploratory
metric and marker development. Causal metrics underwent outcome-informed amendments;
report those analyses as exploratory and disclose the amendment chain. A SHA-256
records file integrity, not independent registration or biological validation.

The [continuation audit](docs/PROGRAM_MODULATION_AUDIT.md) found that the historical
quantitative modulation gates pass a uniform-response counterexample in 20/20
simulations. A separate conditional identity-heterogeneity test and a prospective
[real-data validation plan](config/program_modulation_validation_plan.json) now
distinguish that question from generic perturbation response. Full simulated
outputs and calibration records are [archived](data_program_modulation_audit/README.md).
The later fixed boundary audit shows that unequal generic severity can make all
nine identity tests reject despite rank-one group effects, and that correlated
cells invalidate independent-cell permutations. These tests do not establish
selective pathway regulation or biological replication.

The [fresh Norman continuation](docs/NORMAN_CONTINUATION_RESULTS.md) now includes
a completed full-data, held-out-control marker baseline: 37 heme-panel associations,
of which 29 survive target omission, RNA-depth matching and every negative-control
reference. These are conditional marker-score associations. All three full-size CPU
SAE seeds completed their fixed 200 epochs and are archived with a separate
one-seed control-only sensitivity. No candidate combines the old background gate
with two-seed agreement. See [the interpretation](docs/NORMAN_COMPLETED_MODELS.md).
An ordinary signed-PCA comparison also completed, with all annotation failures
and 945 comparisons archived. Its TF and Cycling features recover different
responses from the direct marker scores; unequal capacities and procedures
preclude a superiority claim. Direct NMF completed its fixed 200-iteration cap
without convergence. A [phase diagnostic](docs/NORMAN_PHASE_DIAGNOSTIC.md) found
that 19 jointly retained Cycling PCA increases coexist with decreases in both
external S and G2M scores; this is score discordance, not measured phase arrest.
A [common six-model protocol](docs/NORMAN_SHARED_INFERENCE.md) and
[signed-weight attribution](docs/NORMAN_PCA_ATTRIBUTION.md) now explain the score
means with all hypotheses and gene contributions retained. The
[claim-to-evidence ledger](docs/CLAIM_EVIDENCE_LEDGER.md) maps manuscript wording
to its complete evidence and limits.

The exact historical eight-gene K562 program now recovers CBL, UBASH3B and PTPN12
across all three complete SAE fits and direct scores. All three pass every control
reference, depth matching, individual-transcript criteria and single-gene omissions.
The list mixes biosynthesis and effector genes and was selected from these same
data in earlier work. This is positive program-association evidence, not independent
confirmation, a new regulator mechanism or an SAE advantage. A bounded NMF warm
restart also met its numerical tolerance and changed annotation coverage despite
only a 0.88% reconstruction-error decrease. All original capped results remain.
The [targeted literature review](docs/LITERATURE_REVIEW.md) records prior biology
and unverified novelty rather than presenting known CBL responses as a discovery.

**Training correction:** every historical causal SAE was fit on all cells, including
perturbed cells, without perturbation labels. Feature-gene associations used control
cells. Cross-fitting separated perturbation matching from testing; it did not hold
test cells out of SAE training. `--train-on control` is available as a sensitivity
mode, but its current synthetic test recovers none of the eight planted regulators.
See [the dated erratum](config/causal_grounding_spec.json).

## Development

The repository contains standalone Python scripts rather than a packaged application.
Gate-0 analysis and synthetic tests run on CPU. Full Perturb-seq fits are larger;
choose resources from the actual dataset and record runtime rather than promising a
GPU duration from an estimate.

Use the configured cloud Python environment when available. For a new environment,
`requirements.txt` contains historical core version pins and open-ended Torch/Scanpy
requirements; it is not a complete environment lock. Record the resolved environment
and follow [REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md).

From the repository root, the standalone checks are:

```bash
python3 -m py_compile src/*.py tests/*.py
for test_file in tests/test_*.py; do
  python3 "$test_file" || exit $?
done
```

Synthetic success verifies specified code paths, not a real-data biological result.
See the reproducibility guide for evidence checks that do not retrain or overwrite
the frozen Gate-0 artifacts. Put new matrices, checkpoints, logs, and analyses in a
separate run directory; do not overwrite `data_g0/` or `data_g0_human/`.

## Repository map

| Path | Contents |
|---|---|
| `src/` | Preprocessing, SAE training, structural metrics, causal metrics, diagnostics, figure generators. |
| `tests/` | Standalone logic, synthetic, preprocessing, and actual-training checks. |
| `config/` | Versioned analysis specifications, marker panels, TF lists, reference gene sets. |
| `data_g0/`, `data_g0_human/` | Committed Gate-0 summaries and frozen decision bundles; large matrices/checkpoints are excluded. |
| `docs/figures/` | Structural figures and a corrected, explicitly historical causal diagnostic; its aggregate sources are not archived causal reruns. |
| `LAB_NOTEBOOK.md` | Dated research history, including retractions and corrections. |

## Reading and publication

- [Master summary](docs/PROJECT_MASTER_SUMMARY.md): project history, evidence registry,
  limitations, file index, and remaining directions.
- [Gate-0 results](docs/COMPONENT0_RESULTS.md) and
  [hardening](docs/COMPONENT0_HARDENING.md): committed structural evidence.
- [Causal results](docs/COMPONENT2_RESULTS.md): current interpretation plus a labeled
  historical record; missing raw evidence is explicit.
- [Abstract draft](docs/ABSTRACT.md): evidence-limited manuscript framing.
- [Manuscript draft](docs/MANUSCRIPT_DRAFT.md): preserved structural evidence, methods audit, and explicit submission gaps.
- [Reproducibility](docs/REPRODUCIBILITY.md): artifact verification and regeneration gaps.
- [Research audit](docs/PROJECT_AUDIT.md): earlier methodological corrections.

Older strategy documents and research plans describe proposed work and superseded
interpretations. Use the dated evidence, executable implementation, and corrections
before quoting them. Including historical causal numerical claims requires archived
reruns. Submission also needs audited citations, confirmed authorship and disclosures,
and validation of any claimed biological novelty.
