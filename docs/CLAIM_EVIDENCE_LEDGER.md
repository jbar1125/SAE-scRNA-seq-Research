# Manuscript claim-to-evidence ledger

Updated 2026-10-20. This maps the reviewable manuscript to complete evidence. A
hash verifies preserved bytes; it does not establish registration, biological
replication or validity of an experimental assumption.

| Claim in the draft | Preserved evidence | Scope required in the wording |
|---|---|---|
| Original asymmetry unsupported in both SAE datasets | `data_g0/v0b_outputs/` and `data_g0_human/v0b_outputs/`, frozen decisions and per-seed CSVs | 0/10 model seeds per dataset; one biological dataset per species; tested operational rule. |
| Only mouse NMF supports the original rule | Both datasets' baseline, GRN and ICA summaries | One of ten method-by-dataset combinations; unequal rank/preprocessing; descriptive comparison. |
| Marker omissions do not restore the claim | Both datasets' robustness/sensitivity summaries | 0/34 mouse and 0/36 human changes; fixed observed cells, not independent biological samples. |
| Generic response passes legacy gates | `data_program_modulation_audit/legacy_counterexamples_20seeds.json` and current reanalysis | 20/20 quantitative-gate passes; unscored manual-coherence criterion excluded. |
| Current null calibration is finite | `data_program_modulation_audit/calibration/` | 6/200 familywise rejections; exact interval 1.11–6.42%; stipulated matched features, no model training. |
| Rank heterogeneity does not establish selective regulation | `data_program_modulation_audit/boundary_audit/` | Graded rank-one generic effects reject all nine in 100/100 worlds. This is correct detection of a narrower quantity. |
| Cell dependence changes null behavior | Same boundary audit | Invalid cell shuffle 100/100 versus separate cluster means 8/100, with wide interval; no claim of general cluster calibration. |
| Fresh direct heme associations survive specified sensitivities | `data_norman_continuation/marker_baseline/` | 37 primary, 29 descriptive intersection; control/detection/normalization rules, depth overadjustment and technical partitions disclosed. |
| KLF1 TF score contains self-target signal | Same archive, original and omission families | AUC .615→.425 after omission; downstream panels persist; not proof of absent regulation. |
| Original full-size SAE endpoint yields no candidate passing both legacy requirements | `data_norman_continuation/sae_primary/` | Three complete 200-epoch seeds; KLF1 repeated but background percentiles63/73/68; seed-2-only TF list. Gates themselves do not prove specificity. |
| Control-only training changes the annotation result | `data_norman_continuation/sae_control_sensitivity/` | One complete seed, full-data preprocessing, trained on inference controls; no independent holdout or general method ranking. |
| Six representations retain different responses under common inference | `data_norman_continuation/shared_inference/` | Full 5,670 association and 54 identity hypotheses; nine panels, all105 singles; unequal capacity and capped nonconverged NMF. |
| PCA-up can coexist with both external phase scores down | `data_norman_continuation/phase_diagnostic/` and `shared_inference/` | 19 jointly significant doubly decreased cases among50 PCA calls in the phase family;20 effect-direction cases among52 in the common family. Different criteria, not conflicting counts. |
| Other weighted genes raise the PCA mean despite falling phase contributions | `data_norman_continuation/pca_attribution/` | All428,505 contributions retained, strict score-fidelity check passes; deterministic mean algebra, not rank or causal mediation. |
| Some dual marker contrasts retain sign on two scales | `data_norman_continuation/dual_marker/` | Joint2,358-test families;21 heme contrasts; pointwise intervals, independent-cell assumption, scale dependence; not validated molecular interactions. |
| CBL/CNN1 erythroid response has prior evidence | `data_literature_review/` and original authors' pinned example code | Original abstract and targeted code reviewed; original full paper/supplements remain unread. No novelty claim. |

The current paper can report the measurements, tested failures and explicitly
exploratory diagnostics. It cannot establish an actual new regulator, causal
lineage mechanism, method superiority or independent biological replication.

Before submission, the author must supply the actual author/affiliation and
contribution statements, funding/conflict disclosures and required assistance
disclosure; verify primary-text claims and select a reuse license; and arrange a
durable archive for large artifacts. No scientific review, approval, funding,
ethics determination, archive identifier or acceptance is implied. These missing
facts are left explicit in [the manuscript](MANUSCRIPT_DRAFT.md).
