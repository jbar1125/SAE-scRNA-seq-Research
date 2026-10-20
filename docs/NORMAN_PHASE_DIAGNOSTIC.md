# Exploratory cell-cycle signature diagnostic

The original ordinary PCA analysis called 52 perturbations upward on a control-
annotated Cycling feature, while the broad direct Cycling score called none.
Twenty-three of those 52 had broad-marker AUC≤0.45. This motivated a declared
post-primary diagnostic using external S and G2M signatures, before viewing
phase-score outcomes. No gene subset was selected from those new outcomes.

The maintained Scanpy tutorial specifies first 43 S and remaining 54 G2M genes in
the Regev list. Its checksum-pinned canonical ZIP was inaccessible; the actual
97-symbol list came from a pinned Broad Institute workshop commit, with a Git
blob matching indexed ETH/NBIS copies. The canonical ZIP identity remains
unverified. Source URLs/hashes are in [the marker provenance](../config/cell_cycle/source.json).
The tutorial cites Tirosh et al.; Crossref verifies 2016 for that paper despite the
tutorial's 2015 wording. This is not a full original-paper gene-list audit.

The original control split reproduced. All-gene normalization/log1p preceded
marker subsetting; discovery controls supplied detection eligibility, means and
scales. One S symbol, MLF1IP, was absent; no aliases were guessed. The other 96
markers met the discovery-control rules. Every 105 single perturbation received
both upward and downward stratified Mann–Whitney comparisons. The 945 original
PCA hypotheses and 420 phase/direction hypotheses entered joint Holm correction
(1,365 total), with a separate same-sized target-omission family. A directional
call also required AUC≥0.55 in that tested direction.

| Signature | Eligible/requested | Up calls | Down calls | Down calls retained after target omission |
|---|---:|---:|---:|---:|
| S | 42/43 | 0 | 52 | 52 |
| G2M | 54/54 | 0 | 37 | 37 |

Fifty of the original 52 PCA-up calls pass the new joint correction. Among those
50, 28 have an omission-consistent S decrease, 21 a G2M decrease, and 19 both. For
the original 52-case set, the corresponding counts are 29, 22 and 20. These are
paired score directions on the same observed cells, not independent experiments.

The simple explanation that one phase's score rises while the other falls is
unsupported under this declared rule: neither phase has an upward call. This
does not prove absence of every weak or phase-specific response. A gene-expression
signature is not a direct measurement of progression, arrest or proliferation.
Different projections can respond to different mixtures of genes and cell states;
these results do not by themselves establish PCA failure or method superiority.
Conditional cell-level p-values and technical partitions do not establish
independent biological replication or a causal mechanism.

All 105 paired estimates, complete directional tests, low-power records, fitted
parameters and all PCA hypotheses with the new correction are in
[`phase_diagnostic/`](../data_norman_continuation/phase_diagnostic/). The helper's
auxiliary upward-only BH fields are preserved in its fit record; diagnostic flags
use the explicit joint Holm families, as identified in the output. The next
question was which signed weights of the saved PCA component can raise its score.
That [separately declared attribution](NORMAN_PCA_ATTRIBUTION.md) is now complete,
with explicit checks that the linear projection approximates the original scores.
