# Targeted literature and citation metadata

`references.json` retains selected Crossref metadata and exact cached-response
hashes. `searches.json` records the EuropePMC query URLs, retrieval times, hit
counts and first-page bibliographic records. `text_review_sources.json` identifies
the two accessed open full texts and their cached XML hashes. Raw API responses
and article XMLs remain outside Git in the environment's literature cache;
their hashes identify what was retrieved but are not a durable full-text archive.

See [the review](../docs/LITERATURE_REVIEW.md) for supported statements and limits.
This is not a systematic search or proof of novelty. Metadata verification,
abstract review, targeted full-text review and experimental reproduction are
different levels of evidence. No funding or ethics statement from another
paper is attributed to this project.

`norman_original_code_review.json` pins verified blobs from the original authors'
repository and identifies reviewed notebook cells. This is code-context review,
not execution or full original-article/supplement review. Ten bibliographic
references now include the Tirosh citation used by the external phase diagnostic.

`norman_fulltext_followup.json` records additional failed public full-text attempts,
the supplementary API's explicit XML error and the linked preprint's metadata.
A public bioRxiv PDF request was denied at CONNECT. These are retrieval diagnostics,
not completed primary-text review.
