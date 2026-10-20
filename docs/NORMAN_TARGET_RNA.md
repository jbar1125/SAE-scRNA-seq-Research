# Target-transcript induction diagnostic

The [declared exploratory diagnostic](../config/norman_target_rna_plan.json) is
complete. All 105 direct target-name comparisons and four alternative transcripts
recorded by the original authors were retained in one 109-test BH family.
Normalization used all genes before extracting transcripts; the original control
split and gemgroup-stratified tests were preserved. A call required q<.05 and
AUC≥.55. No previous phenotype analysis drops labels based on these results.

| Comparison | Planned | Transcript present | RNA-increase calls |
|---|---:|---:|---:|
| Direct target name | 105 | 102 | 90 |
| Author-recorded alternative | 4 | 4 | 2 |

The three absent direct symbols were C19orf26, C3orf72 and KIAA1804. The original
authors' alternatives CBARP and RP5-862P8.2 increased under C19orf26 and KIAA1804
labels (AUC .714 and .630). FOXL2NB under C3orf72 did not meet the effect floor
(AUC .504). RHOXF2B under RHOXF2 likewise fell below it (.539); direct RHOXF2 was
.518. These are separately tested transcript labels, not assumed equivalent genes.
In particular, RHOXF2/RHOXF2B may be distinct paralogous targets; exact guide
sequence specificity has not been established here.

Among 102 observed direct transcripts, 12 failed the combined significance/effect
rule; another three comparisons remained unscored. Several uncalled cases have
small p-values but AUC below .55. A statistically detectable transcript change
does not automatically meet a meaningful effect criterion, and no call does not
prove absent delivery or function. Ninety direct calls plus two alternative-only
calls cover 92 observed labels, without establishing 92 successful functional
gene interventions.

BAK1 and BCL2L11, which the original-author code excluded with an unexplained
“we know these don't work...” comment, show RNA increases here (AUC .793 and .715). We therefore
cannot interpret that comment as failure of RNA induction. It may refer to a
different endpoint; the code comment alone does not identify which. KLF1's own
transcript has AUC .850, consistent with the previously known self-target trap.
These are measurements in the same observed cells, not independent validation.

A post hoc pairing retains all 105 labels and shows that 19 of the 20 fixed
PCA-up/both-phase-down effect-direction cases also have a direct target-RNA
increase. C3orf72 is the exception and has no qualifying alternative-transcript
increase. This does not filter any earlier analysis or add a significance test;
the complete [paired table](../data_norman_continuation/target_rna/cross_diagnostic.csv)
is retained. The score-direction result therefore does not rest entirely on
labels lacking observed target-RNA induction.

All 109 [comparisons](../data_norman_continuation/target_rna/comparisons.csv),
complete [results](../data_norman_continuation/target_rna/results.json), exact
split, source/input hashes, actual run log and validation record are archived.
Machine p/q values that underflow to zero are retained as numerical output; they
do not imply an exact probability of zero. Target RNA induction is not protein
activity, on-target specificity, independent biological replication, mechanism
or novelty. Existing single-guide and cell-dependence limitations remain.
