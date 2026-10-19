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
