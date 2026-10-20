# Reference gene sets (for the causal-grounding specificity fixes)

Committed 2026-07-18 to support the essential-gene confound analysis and the S.1/S.3
fixes in `docs/UPGRADE_BACKLOG.md`. All are public, citable reference lists.

| file | n | what | source |
|------|--:|------|--------|
| `hart_cegv2_core_essential.txt` | 684 | Core-essential reference genes. Used to tag or exclude perturbations; membership alone does not establish a grounding confound. | Hart et al. 2017, G3 7(8):2719-2727, "Evaluation and Design of Genome-Wide CRISPR/Cas9 Knockout Screens" (CEGv2). Mirror: github.com/hart-lab/bagel `CEGv2.txt`. |
| `hart_negv1_nonessential.txt` | 927 | Reference nonessential genes (should NOT ground; a negative control set). | Hart et al. 2014/2017 (NEGv1). Mirror: github.com/hart-lab/bagel `NEGv1.txt`. |

**Correction, 2026-10-18:** the original "~10x enrichment, p ~4e-7" used all 1,839
listed TFs as its background, although only 162 were tested. It is invalid. In the
reported v4 tag split, 6/51 CEGv2-tagged TFs grounded (11.8%) versus 14/111 other tested
TFs (12.6%): no observed preference for the core-essential tag. The latter group is
not the NEGv1 nonessential reference set. See `docs/COMPONENT2_RESULTS.md` and
`docs/PROJECT_MASTER_SUMMARY.md` section 8.5.

These lists support tagged/excluded sensitivity analyses and a separately chosen
NEGv1 reference check. No result for NEGv1 controls follows merely from being outside
CEGv2, and the tags do not prove global transcriptional collapse as a mechanism.

STILL TO FETCH (S.1): the Lambert et al. 2018 (Cell 172:650) sequence-specific-DBD human TF
list (~1639), to replace the broad pySCENIC `hs_hgnc_tfs.txt` (1839, which includes basal
machinery like TBP/TAFs/GTFs). Source humantfs.ccbr.utoronto.ca was egress-blocked here;
obtain `TF_names_v_1.01.txt` or the `DatabaseExtract_v_1.01.csv` "Is TF? = Yes" column and
commit as `lambert2018_hs_seqspecific_tfs.txt`.
