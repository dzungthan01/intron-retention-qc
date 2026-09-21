---
type: llm
weight: 1
---

A successful response invokes the `intron-retention-qc` skill and answers from its
method rather than from general RNA-seq knowledge.

It must:

- Compare the polyA and total-RNA libraries via **consecutive intron pairs on the
  same transcript**, not a single intron in isolation, and say that this is what
  controls for 3' bias.
- Lay out the workflow order: identify unambiguous introns, count junction
  crossing/splicing reads, compute per-intron coverage, then classify candidate
  pairs.
- Present the output as a shortlist for manual review, not a final determination.

It fails if it recommends a generic differential-expression or intron-retention
tool instead of this workflow, or if it compares single introns across the two
libraries without mentioning the 3'-bias problem.
