# Intron Retention / Free-Floating Intron QC

A [Claude Code Skill](https://docs.claude.com/en/docs/claude-code/skills) that identifies retained
and "free-floating" intron signal in bulk RNA-seq data, by comparing polyA-selected vs. total-RNA
("no-select"/rRNA-depleted) libraries of the same biological samples and explicitly correcting for
3' sequencing bias.

See [`SKILL.md`](SKILL.md) for the full method, workflow, and script usage.

## Background

This package packages up the methodology from my BIOL3999 undergraduate research report on
splicing and intron retention (Mus musculus liver/heart/lung, GENCODE vM23, validated against an
8-tissue GEO survey) into a reusable Claude Code Skill: a `SKILL.md` describing when and how to
use it, five Python scripts implementing the workflow steps, and a `references/methodology.md`
with the full worked criteria and known failure modes.

## Status

This is a first draft, not a finished, validated tool. Being upfront about exactly what's been
verified and what hasn't:

- **Direct adaptations of the original report's code** — `scripts/count_junction_signal.py` and
  `scripts/compute_coverage.py` are cleaned-up rewrites of two scripts included in the original
  report's appendix (same core algorithm; `compute_coverage.py` was additionally rewritten as an
  O(n log n) two-pointer sweep instead of the original's per-intron rescan).
- **Reimplemented from the report's prose, not its code** — `scripts/identify_unambiguous_introns.py`,
  `scripts/classify_free_floating_introns.py`, and `scripts/tissue_survey_summary.py` implement
  steps the report described in prose but didn't include runnable code for. Each script's
  docstring says so explicitly.
- **Testing so far**: all five scripts have been smoke-tested end-to-end on small synthetic
  fixtures — they run without errors and produce correctly-shaped output. This is **not** the same
  as biological validation against real SAM/bedGraph files.

### What's needed before this should be considered validated

1. Run it against real data — ideally the same in-house liver/heart/lung samples (or the GEO
   tissue survey) used in the original report, and sanity-check the candidate list against the
   previously-identified validated introns from that report.
2. Re-derive the "low" / "non-low" coverage thresholds in `classify_free_floating_introns.py` from
   an actual negative-control distribution — the current defaults (`1.0`) are placeholders, not
   tuned values.
3. Decide on dependencies — the scripts currently use only the Python stdlib (`re`, `bisect`,
   `csv`, `argparse`). Swapping the hand-rolled SAM/CIGAR parsing for `pysam` would handle edge
   cases (split alignments, clipped mates) more reliably, at the cost of an added dependency.

## License

[MIT](LICENSE)
