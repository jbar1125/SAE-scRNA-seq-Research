# Internal manuscript review copies

These PDF and editable Word files were rendered from `../MANUSCRIPT_DRAFT.md`
with `../../src/build_manuscript_review.py`. They are internal review drafts,
not a submitted, accepted or author-approved article. Authorship, affiliations,
contact, contributions and funding/conflict confirmations remain placeholders.

`render_manifest.json` records exact manuscript, generator, five figure and output
SHA-256 hashes, tool version and PDF properties. The actual 25-page PDF passed
content checks; all five Word figures are embedded and its ZIP CRC check passed.
Title/abstract, reference links, the main comparison figure and supplementary
figure layout were inspected. Numerical code/results were not changed by rendering.

To produce another review copy, run from the checkout with system Pandoc,
XeLaTeX and Poppler available:

```bash
python src/build_manuscript_review.py --out-dir /workspace/a-fresh-review-directory
```

The command refuses an existing output directory. New builds may differ in PDF or
Word metadata; use the recorded source/figure hashes to identify this version.
See `../PUBLICATION_REVIEW_PACKAGE.md` for publication and large-artifact blockers.
