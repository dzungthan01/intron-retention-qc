---
type: llm
weight: 1
---

The method depends on a paired polyA / total-RNA comparison and cannot run on one
library alone. The skill states this under "When NOT to use this".

A successful response says plainly that this workflow cannot answer the question
without a total-RNA or Ribo-Zero counterpart for the same samples, and explains
why: with only polyA data there is nothing to distinguish genuinely absent signal
from signal that polyA selection removed.

It fails if it walks the user through the workflow anyway, or if it suggests the
single-library result would be merely weaker rather than unobtainable.
