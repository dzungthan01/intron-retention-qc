---
name: intron-retention-qc
description: Identifies retained and free-floating intron signal in bulk RNA-seq data by comparing polyA-selected vs. total-RNA ("no-select"/rRNA-depleted) libraries of the same samples, explicitly correcting for 3' sequencing bias. Use whenever a user asks about intron retention, pre-mRNA contamination, splicing QC, free-floating/excised introns, comparing polyA-selected to total-RNA or Ribo-Zero libraries, unambiguous intron-exon junction identification, or wants to survey intron-retention signal across tissues or conditions from bulk RNA-seq. This fills a gap next to the existing nextflow-development skill (gene/transcript-level expression via nf-core rnaseq) and single-cell-rna-qc (per-cell QC) — neither operates at the junction level this skill targets.
compatibility: Requires a Python environment with numpy/pandas/scipy; input is paired-end SAM/BAM alignments, strand-specific bedGraph coverage files, and a GENCODE-style BED/GTF annotation. STAR (or any spliced aligner producing standard CIGAR strings) is assumed upstream.
---

# Intron Retention / Free-Floating Intron QC

Detects whether intron signal observed in RNA-seq is explained by (a) immature pre-mRNA still
attached to flanking exons, or (b) excised introns that have not yet been degraded
("free-floating" introns) — by comparing a polyA-selected library against a total-RNA
("no-select") or rRNA-depleted (Ribo-Zero) library prepared from the *same* biological sample.

This workflow's methodology was developed and validated in the original BIOL3999 report on a
Mus musculus liver/heart/lung study (GENCODE vM23) plus an 8-tissue survey pulled from public GEO
datasets. **The scripts in this repo are a from-scratch reimplementation of that methodology and
have only been smoke-tested on synthetic fixtures, not yet re-run against real sequencing data —
see [Status](README.md#status) for exactly which scripts are direct code adaptations vs.
re-derived from the report's prose.** Nothing below is mouse- or tissue-specific in principle — it
operates on any paired polyA/total-RNA sample pair with a matching GENCODE-style annotation — but
treat outputs as unverified until validated on real data.

## When NOT to use this

- Gene- or transcript-level differential expression → use `nextflow-development`'s `rnaseq`
  pipeline instead.
- Single-cell data → use `single-cell-rna-qc`.
- You only have a polyA-selected library with no total-RNA/no-select counterpart for the same
  sample — the core method here depends on the paired comparison and cannot run on one library
  alone.

## Why 3' bias matters here (read this before running anything)

PolyA selection captures the poly-A tail at the 3' end of a transcript. Any break in the RNA —
native degradation or library-prep shearing — biases sequencing toward the 3' end, so raw
polyA-vs-no-select comparisons at a single intron are confounded: low polyA signal at a 5' intron
can be pure 3' bias, not evidence of anything biological. Every step below that compares the two
library types does so via **consecutive intron pairs on the same transcript**, never via a single
intron in isolation, specifically to cancel this bias out. Don't skip straight to
`classify_free_floating_introns.py` on a single-intron basis — it will produce garbage.

## Workflow

1. **Identify unambiguous introns** — an intron is "unambiguous" if it shares both splice
   junctions with every overlapping intron in every annotated isoform/overlapping gene. This
   excludes cassette exons and alternate splice sites so retention/crossing counts aren't
   confounded by isoform ambiguity.
   ```
   python scripts/identify_unambiguous_introns.py \
     --genes genome.bed --gencode gencode.vM23.annotation.gtf \
     --out unambiguous_introns.bed
   ```

2. **Count junction crossing/splicing reads** — for each unambiguous intron, count reads whose
   CIGAR string crosses the junction (evidence of retention) vs. reads that splice cleanly across
   it. Run once per library (polyA and no-select, each sample).
   ```
   python scripts/count_junction_signal.py \
     --introns unambiguous_introns.bed --sam sample_polyA.sam \
     --out-crossings polyA_crossings.tsv --out-splicings polyA_splicings.tsv
   ```

3. **Compute average per-intron coverage** — from strand-specific forward/reverse bedGraph
   coverage tracks, again once per library.
   ```
   python scripts/compute_coverage.py \
     --introns unambiguous_introns.bed \
     --cov-fwd sample_polyA.fwd.bedgraph --cov-rev sample_polyA.rev.bedgraph \
     --out polyA_avecov.tsv
   ```

4. **Classify candidate free-floating introns** — apply the three-criterion test to every
   consecutive intron pair on a transcript, comparing the polyA and no-select coverage tables from
   steps 2-3:
   - Low signal in the polyA library at the 3'-most intron of the pair
   - *Higher* signal in the polyA library at the subsequent (more 5') intron — this is what rules
     out plain 3' bias, which would predict the opposite direction
   - Non-low signal at that same 3'-most intron in the no-select library
   ```
   python scripts/classify_free_floating_introns.py \
     --polyA-cov polyA_avecov.tsv --noselect-cov noselect_avecov.tsv \
     --ratio-threshold 0.3 --out candidate_introns.tsv
   ```
   `--ratio-threshold` is the minimum fold-increase required in the polyA library from the 3'-most
   intron to the next one (0.3 was the value used in the source study; treat it as a tunable
   knob, not a universal constant — re-derive it from your own negative-control distribution if
   the sequencing depth or library prep differs substantially).

5. **(Optional) Cross-tissue/condition survey** — summarize retention signal per
   tissue/condition so multiple samples can be compared on one axis:
   ```
   python scripts/tissue_survey_summary.py \
     --crossings *_crossings.tsv --avecov *_avecov.tsv --out tissue_summary.tsv
   ```
   The summary statistic is `sum(avg_coverage - crossing_reads)` over junctions where that
   difference is positive, normalized by library read depth. Higher = more retention/free-floating
   signal relative to properly-spliced reads.

## Output

`candidate_introns.tsv` — one row per candidate intron pair, with the polyA/no-select coverage
values that triggered the flag. Treat this as a **candidate list for manual review**, not a final
call: in the source study, roughly half of automatically flagged candidates were ruled out by hand
after inspection in a genome browser. Recommend surfacing IGV/UCSC Genome Browser coordinates
alongside each candidate so a human can do that review quickly — see
`references/methodology.md` for the visual criteria used to accept/reject candidates by eye.

## Validating candidates

The source methodology validates top candidates with CAGE-seq (Cap Analysis of Gene Expression):
if the 5' cap / transcription start site for a "free-floating" candidate intron sits *within* the
intron rather than at the annotated gene TSS, that supports the intron being transcribed
independently rather than merely retained. This skill does not automate CAGE-seq analysis; treat
it as the recommended follow-up wet/computational validation step, not something to skip.

See `references/methodology.md` for the full worked criteria, known failure modes, and the
rationale behind each threshold.
