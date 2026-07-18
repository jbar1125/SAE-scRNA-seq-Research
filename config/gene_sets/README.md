# Reference gene sets (for the causal-grounding specificity fixes)

Committed 2026-07-18 to support the essential-gene confound analysis and the S.1/S.3
fixes in `docs/UPGRADE_BACKLOG.md`. All are public, citable reference lists.

| file | n | what | source |
|------|--:|------|--------|
| `hart_cegv2_core_essential.txt` | 684 | Core-essential genes (fitness genes essential across cell lines). Used to detect/exclude the essential-gene grounding confound. | Hart et al. 2017, G3 7(8):2719-2727, "Evaluation and Design of Genome-Wide CRISPR/Cas9 Knockout Screens" (CEGv2). Mirror: github.com/hart-lab/bagel `CEGv2.txt`. |
| `hart_negv1_nonessential.txt` | 927 | Reference nonessential genes (should NOT ground; a negative control set). | Hart et al. 2014/2017 (NEGv1). Mirror: github.com/hart-lab/bagel `NEGv1.txt`. |

Why these are here (see `docs/COMPONENT2_RESULTS.md`): the first real Arm-B run grounded a
set that is ~10x enriched for CEGv2 core-essential genes (8/26 seed-stable hits, vs a 3.0%
base rate among the 1839 scored TFs, hypergeometric p ~ 4e-7). That means much of the
grounding is essential-knockdown transcriptional collapse, not regulatory specificity.
These lists let the pipeline (a) exclude core-essential perturbations and (b) verify that
nonessential-reference genes do not ground.

STILL TO FETCH (S.1): the Lambert et al. 2018 (Cell 172:650) sequence-specific-DBD human TF
list (~1639), to replace the broad pySCENIC `hs_hgnc_tfs.txt` (1839, which includes basal
machinery like TBP/TAFs/GTFs). Source humantfs.ccbr.utoronto.ca was egress-blocked here;
obtain `TF_names_v_1.01.txt` or the `DatabaseExtract_v_1.01.csv` "Is TF? = Yes" column and
commit as `lambert2018_hs_seqspecific_tfs.txt`.
