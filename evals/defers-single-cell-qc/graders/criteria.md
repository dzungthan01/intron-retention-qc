---
type: llm
weight: 1
---

The word "intronic" appears, but this is single-cell data, which
`intron-retention-qc` explicitly excludes under "When NOT to use this". This case
checks that the skill's description does not over-trigger on intron vocabulary
alone.

A successful response routes to per-cell QC (the `single-cell-rna-qc` skill) and
does **not** run the bulk intron-retention workflow.

It fails if it invokes `intron-retention-qc`, or if it proposes comparing
polyA-selected against total-RNA libraries — that comparison needs paired bulk
libraries the user does not have.
