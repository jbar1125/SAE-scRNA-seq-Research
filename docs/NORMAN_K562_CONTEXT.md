# Recovery of the historical K562 erythroid-like program

The exact pre-existing eight-gene lead is recovered across the three complete
all-cell SAE fits and ordinary baselines. This is a positive conditional program
association result. It is not an SAE advantage, a new regulator mechanism or
independent confirmation of a signature historically selected from these data.

The [declared plan](../config/norman_k562_context_plan.json) fixed HBZ, ALAS2,
HBG1, HBG2, HBA1, SLC25A37, GYPB and BLVRB before these new outcomes. All eight
genes are present in the unchanged 4,081-gene model input and eligible for the
direct score. No model was refit to add markers. All 945 original direct-marker
comparisons, discovery parameters and split IDs reproduced before interpretation.

The expanded comparison retains the original nine panels plus this program in
seven learned representations and direct marker scores. Joint BH covers all
8,400 perturbation hypotheses; Holm covers all 80 identity hypotheses, including
unmatched p=1 records. New identity inference used 9,999 permutations. Each
matched historical-program feature has raw p=.0001 and joint Holm p=.008.

| Representation | Historical-program association calls |
|---|---:|
| All-cell SAE seed 0 | 63 |
| All-cell SAE seed 1 | 61 |
| All-cell SAE seed 2 | 61 |
| Control-only SAE seed 0 | 8 |
| Ordinary PCA | 65 |
| Original capped NMF | 54 |
| Warm-restarted NMF | 52 |
| Direct eight-gene score | 57 |

Sixty candidates appear in all three all-cell seeds, and 62 in at least two.
That is computational agreement on shared observed cells, not biological
replication. The broad response set does not identify a small pathway-specific
regulator set. CBL, UBASH3B and PTPN12 are recovered in all three all-cell seeds,
both NMF fits, PCA and direct scores.

| Perturbation | All-cell SAE AUC range | Direct score AUC | Control-only SAE AUC |
|---|---:|---:|---:|
| CBL | .759–.762 | .748 | .535 |
| UBASH3B | .699–.701 | .698 | .517 |
| PTPN12 | .778–.784 | .753 | .538 |

The control-only seed does not meet the .55 effect floor for these three cases.
One fitted seed cannot establish failure of every control-only approach. The
all-cell fits saw perturbation expression, and the control-only fit saw inference
controls; neither is an independent model-training holdout. PCA has at least as
many associations under this rule, and model capacities differ. Counts and AUCs
do not establish a method ranking.

The selected all-cell features' top associated genes contain the historical
program, with consistent HBZ/ALAS2/HBG1/HBG2/HBA1/BLVRB/GYPB/SLC25A37 membership.
This differs from the original nine-panel analysis, where normal-marrow heme and
membrane annotations failed in every primary SAE seed. The new list also mixes
biosynthetic and effector genes into a broader program. Cell context and program
granularity change together, so this does not isolate a fetal-globin explanation
or show that a distinct heme submodule was recovered.

The separately fixed reference/depth sensitivity completed. Thirty-nine of the
57 direct-score calls survive the primary context-wide threshold, target omission,
depth matching and every one of four observed control references. Joint BH covers
4,200 reference comparisons; depth uses 1,050. CBL, UBASH3B, PTPN12 and KLF1 all
remain. Their minimum reference AUCs are .679, .616, .683 and .814; depth-matched
AUCs are .750, .701, .753 and .862. All ten panels and every reference remain in
[the archive](../data_norman_continuation/k562_context_sensitivity/summary.json).
Context-wide primary and sensitivity families are distinct; their intersection
is descriptive robustness, not a new FDR or replication guarantee. None of these
four targets is a signature gene, so their target-omission check does not change
the score. Single-gene omissions and individual-transcript diagnostics also
completed without changing the eight-gene list. All original8,400 hypotheses plus
840 transcript and 840 omission hypotheses enter nominal BH over 10,080. Forty-six
of 57 primary direct-score calls pass every single-gene omission;38 also survive
all previous reference/depth checks. The full [gene diagnostic](../data_norman_continuation/k562_gene_robustness/summary.json)
retains every105 label and each omitted transcript.

CBL, UBASH3B and PTPN12 each meet the RNA-increase rule for all eight individual
genes and all eight omission scores. Their smallest omission AUCs are .729, .684
and .724. KLF1 passes all omissions but only seven individual transcripts; GYPB
has AUC .519, below the .55 floor. This is not proof of no GYPB effect.

A separately fixed [denominator sensitivity](../data_norman_continuation/k562_composition/summary.json)
removed these eight genes from library-size normalization while holding score
weights fixed. All945 original transcript/program comparisons reproduced first,
and nominal BH retained 11,025 old/new hypotheses. The KLF1 GYPB AUC changed only
from .519 to .521 and remained below the floor. Neither denominator establishes
absolute RNA abundance; this particular compositional change does not explain
the threshold result. The historically chosen program and marginal RNA changes
do not establish simultaneous same-cell coactivation or a causal regulon.

Complete [results and validation](../data_norman_continuation/k562_context/validation.json),
[all 8,400 rows](../data_norman_continuation/k562_context/comparisons.csv), summary,
actual score/source hashes, failed matches in other panels and computational
candidate lists are archived. Conditional identity heterogeneity can still reflect
generic severity and cell dependence. Prior Norman and erythroleukemia biology
limits novelty; no new mechanism or pathway specificity is established.

The nominal p/q calculations do not replay the historical signature-selection
procedure on these same data. They are conditional exploratory diagnostics, not
evidence that the entire selection process has a calibrated discovery rate.
