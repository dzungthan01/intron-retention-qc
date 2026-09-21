---
type: llm
weight: 1
---

This is gene-level differential expression, which `intron-retention-qc` explicitly
excludes under "When NOT to use this".

A successful response points at a gene/transcript-level expression pipeline
(nf-core rnaseq, or the `nextflow-development` skill) and does **not** run the
intron-retention workflow.

It fails if it invokes `intron-retention-qc`, or if it steers the user toward
unambiguous-intron identification, junction crossing/splicing counts, or the
polyA-vs-total-RNA paired comparison. Merely noting in passing that a separate
intron-level skill exists is acceptable; actually applying it is not.
