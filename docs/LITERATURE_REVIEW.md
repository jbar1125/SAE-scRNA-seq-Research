# Literature review and novelty limits

Reviewed 2026-10-18 ET during CPU continuation. This is a targeted review, not a
systematic review or evidence that an unreturned association has never been
published. Search metadata and source fingerprints are in
[`data_literature_review/`](../data_literature_review/).

| Proposed interpretation | Verified source and review scope | Consequence for this project |
|---|---|---|
| CBL-linked erythroid responses are new biology | Norman et al., *Science* 2019, [10.1126/science.aax4438](https://doi.org/10.1126/science.aax4438), abstract read through EuropePMC. It explicitly reports CBL/CNN1 synergy driving erythroid differentiation. | A CBL erythroid association in these same data is contextual replication, not sufficient evidence of novelty. The abstract does not establish our single-perturbation effect size or verify our null. |
| UBASH3B has no prior erythroleukemia connection | Wang et al., *BMC Cancer* 2024;24:326, [10.1186/s12885-024-12075-2](https://doi.org/10.1186/s12885-024-12075-2). Full-text methods/results reviewed. FLI1 perturbation regulates UBASH3B in HEL and K562 contexts; UBASH3B knockdown is investigated for proliferation and signaling. | Prior related biology exists. These different interventions and outcomes neither validate our CRISPRa heme score nor establish erythroid differentiation. Our UBASH3B score also fails one control-reference effect floor. |
| SAE circuit structure establishes cellular causal regulation | Kendiukhov, *Bioinformatics* 2026;42(7):btag379, [10.1093/bioinformatics/btag379](https://doi.org/10.1093/bioinformatics/btag379). Full-text biological-prior, annotation-selection, CRISPRi and limitation sections reviewed. | That work distinguishes model-internal ablations from cellular causal structure and shows the importance of annotation selection and a sign-bias null. We do not claim to originate SAE interpretation or compare its reported percentages with our different estimands. |
| PTPN12 is a new erythroid regulator | Targeted PTPN12/PTP-PEST searches returned related papers, including the original Norman article through indexed search. We have not completed its full text and supplements or established an independent validation experiment. | Novelty remains unestablished. A search hit is not evidence of a specific directional result; absence from a title or abstract is not evidence of novelty. |

The Norman paper is PMID 31395745 / PMC6746554. EuropePMC lists a free HTML copy,
but article content was not retrieved. Initial requests returned 403; later PMC
requests served a browser challenge and EuropePMC served an application shell.
The REST full-text endpoint returned 500 and the no-JavaScript article route 520.
Only its abstract has been reviewed here. Two other open-access full-text XMLs
were retrieved successfully. Their targeted review does not independently
reproduce their experiments. No access restriction was bypassed.

The broader SAE search also returned preprints and papers about foundation-model
interpretability. Metadata or an abstract is insufficient to validate their
numerical claims. These hits are retained for follow-up, not treated as established
comparators. The search returned only its first page (up to 30 results per query),
includes ambiguous gene aliases, and did not exhaustively screen references or
supplements. Its hit counts are retrieval diagnostics, not discovery evidence.

Crossref verified bibliographic fields for Paul, Setty, Benjamini–Hochberg,
Scanpy, Replogle, Norman, TRRUST v2, Wang, Kendiukhov and Tirosh. The phase-diagnostic
Tirosh citation is a 2016 paper, despite the tutorial's 2015 wording; its primary
text was not read here. TRRUST v2 was online in 2017
and appears in the 2018 print issue; the manuscript uses 2018. Crossref metadata
checks do not substitute for primary-text verification of scientific claims.

On 2026-10-20, the original authors' public
[analysis repository](https://github.com/thomasmaxwellnorman/Perturbseq_GI) supplied
additional primary-source code context. Exact Git blobs were verified for its
README, CBL/CNN1 examples and guide-performance notebooks. The latter explicitly
distinguishes the single guide used for Perturb-seq from two-guide averages in
the separate fitness screen, and compares guide cassette positions A/B. The
examples include PTPN12/SNAI1 and CBL/CNN1 cDNA-expression plotting code for CD235a
and CD71. These observations further limit novelty and independent-guide claims;
they do not verify our specific single-target effect sizes. Targeted code/markdown
was read, not executed. This does not complete original-article or supplement
review. [Source hashes and cell references](../data_literature_review/norman_original_code_review.json)
record the scope.

Before proposing a new regulator mechanism, review the original Norman full text
and supplements; distinguish single perturbations from combinations and CRISPRa
from knockdown; examine prior work for the exact gene, direction, cell context and
outcome; and obtain independent experimental evidence. Computational seed agreement
or technical gemgroups do not supply that evidence. The currently defensible
contribution remains the measured representation dependence and the documented
failure modes of the validation design.
