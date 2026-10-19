# Norman continuation: fresh archived results

The full-data external-marker baseline is complete. The full-size SAE experiment
is running and has no completed result yet. These are separate analyses. All new
marker results below have machine-readable outputs in
[`data_norman_continuation/`](../data_norman_continuation/), unlike the lost
historical Component 2 outputs. No new biological mechanism or SAE advantage is claimed.

## Dataset and design

The raw Norman file was downloaded from the public scverse S3 source, checked at
828,160,302 bytes, and assigned observed SHA-256
`32cc993abad41336eb0f3635f97e802e13fe5721dd9277c1c55e7c13853782e0`.
This is an observed source fingerprint, not an independently supplied checksum.
The existing preparation produced 111,445 cells: 11,855 controls, 57,831 single-
perturbation cells, and 41,759 dual-perturbation cells. All 105 single labels have
enough cells for the marker baseline; the older grounding loader additionally
requires a measured target transcript and therefore calls only 102 testable.

All nine external human marker panels were retained. Every duplicate raw-count
symbol column was summed before normalization: 33,694 columns become 33,660
unique symbols, with exact per-cell count conservation. The 34 duplicate symbols
affect neither the 59 external markers nor the 105 single-target symbols. The
original files remain unchanged. Separate unique-count preparation supplies the
SAE input and retains source Ensembl identifiers.

For the marker baseline, all-gene library normalization to 10,000 and log1p
preceded marker subsetting. Outcome-blind, within-gemgroup control splitting
gave 5,926 discovery and 5,929 inference controls. Discovery controls alone
determined marker detection, means and scales. Three panels had fewer than
three eligible markers and remained unscored.

One-sided stratified Mann–Whitney comparisons used the eight `gemgroup` technical
partitions. All partitions contain controls and every single perturbation.
BH correction covered all 9 × 105 = 945 comparisons, including unscorable p=1
records. A response required q < 0.05 and AUC ≥ 0.55. Removing a perturbation's
own transcript was a separately corrected sensitivity using the same family.
See the [declared plan](../config/marker_program_baseline_plan.json).

## Marker-response results

| Program | Eligible markers | Primary responses | Survive target omission | Also survive depth matching and every control reference |
|---|---:|---:|---:|---:|
| Ery_TF | 7 | 7 | 6 | 3 |
| Ery_Heme | 7 | 37 | 37 | 29 |
| Ery_Membrane | 7 | 25 | 25 | 20 |
| Ery_Effector | 2 | Unscored | Unscored | Unscored |
| Gran_TF | 3 | 6 | 4 | 4 |
| Gran_Primary | 1 | Unscored | Unscored | Unscored |
| Gran_Secondary | 2 | Unscored | Unscored | Unscored |
| Progenitor | 3 | 3 | 3 | 2 |
| Cycling | 16 | 0 | 0 | 0 |

The final column is exploratory. Its plan was written after the primary results,
before sensitivity outcomes. Scores, marker eligibility and the original split
remained fixed. Depth matching used discovery-control UMI quartiles within each
gemgroup; because RNA depth can be downstream of perturbation, this can overadjust
real biology. Control-reference tests used each of four observed negative-control
construct identities, with BH over all 3,780 program × perturbation × reference
comparisons. Intersections are descriptive robustness, not independent replication.
The primary computations were reproduced exactly before running sensitivities.

The full results, all references and descriptive per-partition effects are
archived; no favorable reference or technical partition was selected. See the
[sensitivity plan](../config/marker_baseline_sensitivity_plan.json) and
[945-row comparison table](../data_norman_continuation/marker_baseline/comparisons.csv).

## A direct-target trap and a reference-dependent candidate

KLF1 overexpression has TF-panel AUC 0.615, but the target-omitted TF score has
AUC 0.425 and does not pass the activation criterion. Its heme-panel AUC is 0.626
and membrane-panel AUC 0.580; those panels do not contain KLF1, and both associations
survive all specified sensitivities. Thus a transcript's presence in its own
annotation panel and a downstream marker response must be distinguished.

CEBPA and CEBPE also pass the primary Gran_TF panel, but removing their target
leaves only two eligible markers, below the three-marker rule. Their omission
analyses are **unscorable**, not evidence of absent downstream effects.

CBL and PTPN12 retain heme and membrane associations in all sensitivities. UBASH3B
passes the pooled heme comparison (AUC 0.571) and depth sensitivity, but its lowest
control-reference AUC is 0.548, below the fixed 0.55 effect floor. It therefore
does not belong in the all-reference robust set. This limits the previous
three-candidate narrative rather than confirming it wholesale.

![Archived marker baseline and KLF1 diagnostic](figures/fig_norman_marker_baseline.png)

Panel B is a post hoc diagnostic of an existing candidate. Both plots render the
archived outputs; they are not independent validation experiments.

## Conditional identity test and interpretation

All six adequately covered panels, including Cycling, reject the within-gemgroup
identity-shuffle null at the minimum raw p of 0.001 and Holm p 0.009 over nine
planned panels. Cycling nevertheless has no positive response under the
one-sided activation/effect-floor rule. Identity variation and activation are
different questions; neither finding establishes erythroid pathway specificity.
Three low-power panels retain p=1.

These are conditional cell-level associations in one cell line. Technical
partitions are not biological replicates. The observed `guide_identity` strings
often differ only by vector orientation, and CBL has a single observed construct
identity; do not call these independent sgRNA replications. Generic perturbation
severity, compositional normalization, cell state, viability, construct effects
and unobserved confounding remain possible. Literature novelty, independent
experimental validation and an SAE advantage have not been established.

## Full-size SAE status

The user selected CPU execution. Seed 0 is training on all 111,445 cells with
4,081 selected/force-kept genes, 2,048 latent features, k=32, batch=1,024 and a
200-epoch target. Atomic checkpoints retain input hashes, weights, optimizer,
means and RNG states; small tests verified exact uninterrupted/resumed equality.
Seeds 1 and 2 remain pending. A running or partial seed is not a completed result.
Use [the run instructions](NORMAN_CPU_RUNS.md) and timestamped progress records
for current status rather than an old duration estimate.

### Interim inference readiness check

The [declared diagnostic](../config/norman_interim_diagnostics_plan.json) used a
separate copy of the epoch-20 checkpoint from the 200-epoch parent fit. Inference
and the complete nine-panel analysis finished on all 111,445 cells without further
optimizer steps. This is a partial-model implementation check, not a completed
primary seed, a selected best epoch, or an independent replicate.

Four panels matched an eligible SAE feature: Ery_TF, Ery_Heme, Ery_Effector and
Cycling. Each rejected the conditional identity null (per-seed Holm p=.009);
the other five hypotheses retained p=1. The heme feature had seven exploratory
modulators and a descriptive feature-background percentile of 86. None of the
matched features reached the old 90th-percentile background threshold, which is
itself descriptive and does not establish specificity. Full outputs, negatives
and checkpoint provenance are retained in
[`sae_interim_epoch20/`](../data_norman_continuation/sae_interim_epoch20/).

Annotation coverage differs between analyses: the direct marker baseline requires
detection in at least ten discovery controls, while SAE matching checks marker
presence. HBB was detected in only three discovery controls. An Ery_Effector SAE
match alongside an unscorable direct score therefore does not establish an SAE
advantage. The fixed 200-epoch primary fits remain unchanged.
