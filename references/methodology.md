# Methodology reference

Source: undergraduate research report, ITMAT Bioinformatics Laboratory, Perelman School of
Medicine, University of Pennsylvania (2024). Validated on Mus musculus mm10 / GENCODE vM23,
in-house liver/heart/lung samples plus an 8-tissue public GEO survey (embryonic stem cells,
pancreatic islet cells, lung, kidney, hippocampus, pituitary, quadriceps, colon).

## Library design this method assumes

Two libraries per sample:
- **PolyA-selected**: standard practice, enriches mature mRNA, presumed to exclude free-floating
  (non-polyadenylated) intron signal almost entirely.
- **"No-select" / total RNA, or Ribo-Zero rRNA-depleted**: captures all sequences, not just
  polyadenylated ones. Since ~80-90% of total RNA is ribosomal, sequence the no-select library at
  roughly 10x the read depth of the polyA library so non-ribosomal signal is comparably powered
  across the two — otherwise a real depth difference will look like a biological difference.

## The three-criterion test, in full

For each pair of consecutive introns (A = more 3', B = the next intron 5' of A) on the same
transcript:

1. **Low polyA signal at A.** Near-zero average coverage for A in the polyA library.
2. **Higher polyA signal at B than at A** (by at least `ratio_threshold`, default 0.3 — i.e., B's
   coverage exceeds A's coverage by that fraction). This is the bias-control step: if A's low
   signal were pure 3' bias, B (further from the 3' end) should be *lower* still, not higher.
   Seeing B higher than A rules that out.
3. **Non-low signal at A in the no-select library.** A is actively transcribed/present in the
   no-select data — the low polyA signal at A isn't simply "nothing is there."

An intron pair passing all three is a **candidate** free-floating-intron locus: signal that
exists in the no-select assay but can't be explained by retained pre-mRNA (which would show up in
polyA too) or by 3' bias (ruled out by criterion 2).

## Manual review after automated flagging

In the source study, an automated candidate list of 59 introns was reviewed by hand in a genome
browser; roughly half were ruled out as false positives. Common rejection reasons to check for:
- Low read depth at the locus generally (candidate is noise, not signal)
- Annotation issues (the "intron" boundary doesn't match visual alignment patterns)
- The apparent signal increase at B is itself explainable by a nearby, unrelated feature

Roughly a third of the reviewed candidates were kept as strong evidence (the source study's final
count: 20 of 59). Set expectations with users accordingly — this skill produces a *shortlist for
review*, not a final determination.

## The cross-tissue/condition summary statistic

For each exon-intron junction: `residual = average_coverage - crossing_read_count`. Sum
`residual` over all junctions in a sample **where the residual is positive**, then normalize by
that sample's sequencing depth. Rationale: crossing reads are direct evidence of intron
*retention* (still attached to flanking exons); coverage in excess of that is more consistent with
signal not explained by retention alone. Comparing this statistic across tissues/conditions
(rather than raw coverage) controls for differences in overall retention baseline.

In the source study, kidney, hippocampus, and quadriceps showed the highest values among eight
surveyed mouse tissues; a matched polyA sample was included as a negative control and scored near
zero, as expected.

## Recommended validation: CAGE-seq

CAGE-seq (Cap Analysis of Gene Expression) captures the 5' cap of transcripts, pinpointing
transcription start sites (TSS). For a candidate free-floating intron, a TSS that falls *within*
the intron (rather than at the gene's annotated TSS) supports independent transcription of that
intron rather than simple retention. This is the recommended orthogonal validation step for any
candidate this skill surfaces before treating it as confirmed.

## Known limitations to disclose to users

- The `ratio_threshold` (default 0.3) was tuned on one study's sequencing depth and library prep;
  it is a starting point, not a validated universal constant.
- "Unambiguous intron" filtering removes cassette exons and alternate-splice-site introns by
  design — this method does not attempt to resolve retention signal in ambiguous/overlapping
  regions.
- This is a candidate-generation and hypothesis-generating method, not a confirmatory one. Always
  present results as "candidates for review," and recommend CAGE-seq or manual genome-browser
  inspection before treating a hit as confirmed.
