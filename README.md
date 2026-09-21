# Intron Retention / Free-Floating Intron QC

[![tests](https://github.com/dzungthan01/intron-retention-qc/actions/workflows/tests.yml/badge.svg)](https://github.com/dzungthan01/intron-retention-qc/actions/workflows/tests.yml)

A [Claude Code Skill](https://docs.claude.com/en/docs/claude-code/skills) that identifies
retained and "free-floating" intron signal in bulk RNA-seq, by comparing polyA-selected
against total-RNA libraries of the same biological sample and explicitly correcting for
3' sequencing bias.

See [`SKILL.md`](SKILL.md) for the workflow and script usage.

## Background

RNA-seq routinely shows signal over introns, and there are two competing explanations for
it. Either the reads come from immature pre-mRNA that still has its introns attached to
the flanking exons, or they come from introns that have already been excised but not yet
degraded — "free-floating" introns. Coverage plots alone cannot tell these apart, because
they show depth without showing whether a molecule is still connected to its neighbours.

The BIOL3999 research this skill packages up separates the two using a paired library
design. PolyA selection captures the poly-A tail, so it enriches mature mRNA and should
almost entirely exclude excised introns, which are not polyadenylated. A "no-select" or
Ribo-Zero library captures everything regardless of tail. Signal present in no-select but
absent in polyA is therefore a free-floating intron candidate, with the polyA library
acting as the control. Because total RNA is 80–90% ribosomal, the no-select libraries were
sequenced roughly ten times deeper so the non-ribosomal fraction was comparably powered.

The complication is that polyA selection introduces a 3' bias: any break in the RNA biases
sequencing toward the poly-A tail, so coverage falls off toward the 5' end. A single intron
compared across the two library types is therefore uninterpretable — low polyA signal at a
5' intron may be pure bias. The method controls for this by only ever comparing
**consecutive intron pairs within the same transcript**. If the 3'-most intron of a pair is
silent in polyA while the next one further 5' is *higher*, that is the opposite of what 3'
bias predicts, and the bias explanation is ruled out.

The work was carried out on in-house mouse liver, heart and lung libraries at the ITMAT
Bioinformatics Laboratory, plus an eight-tissue survey drawn from public GEO data. It
produced a shortlist of 59 candidate introns, of which 20 survived manual review in a
genome browser, and proposed CAGE-seq as the orthogonal validation step.

## This skill

Five scripts implement the five computational steps of that workflow. The sixth step —
CAGE-seq validation of surviving candidates — is deliberately out of scope; this skill
generates the shortlist that step would then test.

Two of the scripts are adaptations of code included in the report's appendix. The other
three implement steps the report described in prose but did not ship runnable code for,
which matters when interpreting their output — see [Future validation steps](#future-validation-steps).

| # | Script | What it does | Provenance |
|---|---|---|---|
| 1 | `identify_unambiguous_introns.py` | Builds a canonical intron set. An intron is "unambiguous" if it shares both splice junctions with every overlapping intron in every annotated isoform, which excludes cassette exons and alternate splice sites so later counts aren't confounded by isoform ambiguity. | reimplemented from prose |
| 2 | `count_junction_signal.py` | Per intron, counts reads whose CIGAR *crosses* an exon-intron junction (evidence of retention) against reads that *splice* cleanly across it (evidence of normal processing). Handles mate-pair de-duplication and skips multi-mappers. | adapted from appendix code |
| 3 | `compute_coverage.py` | Computes average per-base coverage for each intron from strand-specific bedGraph tracks, as an O(n log n) sweep rather than the original's per-intron rescan. | adapted from appendix code |
| 4 | `classify_free_floating_introns.py` | The three-criterion test, applied to consecutive intron pairs (A = 3'-most, B = next toward 5'): A silent in polyA; real signal at B exceeding A; A present in no-select and above B. Emits a candidate shortlist. | reimplemented from prose |
| 5 | `tissue_survey_summary.py` | Reduces each sample to one number so tissues or conditions can be compared on a single axis, from the residual between coverage and crossing-read count. | reimplemented from prose |

`references/methodology.md` carries the full criteria, the thresholds' rationale, and the
known failure modes.

## QA

Validation lives in `tests/` (script behaviour) and `evals/` (skill triggering).

```bash
python3 -m pytest tests/ -q            # 56 passed, 3 xfailed — free, offline
claude plugin eval . --trust-plugin    # optional; launches billed agent runs
```

### Unit tests

56 tests, no sequencing data required, about a second to run.

- **Docs/CLI contract** (`test_docs_contract.py`) — every command in `SKILL.md` is parsed
  and checked against each script's argparse declarations via `ast`, in both directions:
  documented flags must exist, and required flags must be documented. This caught a
  documented command that exited 2.
- **Numerical behaviour** (one file per script) — hand-computed fixtures rather than
  smoke tests: coverage arithmetic and strand routing, junction counts at the overlap
  cutoff boundary (including the off-by-one either side), mate-pair de-duplication,
  multi-mapper rejection, the unambiguous-intron filter, 3'-position ranking on both plus
  and minus strands, and each of the three classification criteria including the case
  where 3' bias should cause rejection.
- **Statistic equivalence** — the summary statistic supports both forms the report
  describes, and a test pins the guarantee that they agree exactly on zero-crossing
  junctions, where a naive ratio would divide by zero.

Known defects are recorded as `strict=True` xfails, so the suite stays green while the
issue is open and turns red the moment one is fixed — the marker has to be deleted
deliberately. Three remain, none affecting the science:

| Defect | Where |
|---|---|
| Mates are taken as consecutive lines with no QNAME check, so a coordinate-sorted SAM pairs unrelated reads and still produces plausible counts | `count_junction_signal.py:137` |
| `--depths liver=…` never matches the derived sample name `liver_crossings`, so depth normalization silently no-ops | `tissue_survey_summary.py:47` |
| The reported intron count ignores the `--chromosomes` allowlist applied at write time | `identify_unambiguous_introns.py:172` |

### Trigger eval design


Four cases under `evals/` test that boundary

| Case | Prompt | Passes when |
|---|---|---|
| `triggers-on-intron-retention` | Mouse liver, polyA + total-RNA from the same samples; introns strong in total-RNA and near-absent in polyA — retained pre-mRNA or free-floating? | The skill fires, and the answer works from consecutive intron pairs, names the 3'-bias control, gives the workflow in order, and calls the output a review shortlist |
| `defers-gene-level-expression` | Differential expression between liver and kidney, three replicates each — what pipeline? | Routes to nf-core `rnaseq` / `nextflow-development`. **Fails if this skill fires.** |
| `defers-single-cell-qc` | 10x single-cell data with high intronic read fractions — how to QC? | Routes to `single-cell-rna-qc`. **Fails if this skill fires.** |
| `requires-a-paired-library` | Only polyA, no total-RNA counterpart — can I still find free-floating introns? | Says the method cannot run, and why. Fails if it walks the workflow anyway. |

`defers-single-cell-qc` is the sharpest of the four: it is dense with intron vocabulary but
is the wrong modality, so it catches a description that over-triggers on keywords alone.

The positive case carries two graders — an `llm` grader scoring the answer's content, and a
`tool_used: Skill` grader that deterministically checks the skill actually loaded rather
than the model answering from memory.

**Runs and arms.** Each case runs three times by default, because model output is
nondeterministic and one pass cannot separate "triggers reliably" from "triggered once";
the score is the pass rate. Each case also runs under two ablation arms — once with the
plugin loaded, once without — and the report gives the delta. That is what attributes the
behaviour to the skill rather than to what the base model already knew, and it matters most
for the two negative cases: if the baseline also routes single-cell work correctly, those
cases say nothing about this skill.

4 cases × 3 runs × 2 arms = 24 agent runs.

## Future validation steps

**1. Identify the dataset behind the published figures.** The
source data exists as several output lineages — `avecov` files keyed on transcript
(`ENSMUST…`, ~121k introns) and older `AVE_COV` files keyed on gene (`ENSMUSG…`, ~139k
introns) — across polyA, no-select and two Ribo-Zero preps, for liver, heart and lung.

**2. Score against the 59-candidate shortlist.** Once the dataset is confirmed, measure
recall against the 59 flagged regions and the 20 that survived manual review.

**3. Calibrate thresholds on the confirmed negative control.** The "low polyA" and "non-low
no-select" cuts are currently 0.0 and 0.1, derived from the liver coverage distributions
rather than from an arbitrary constant, but they need re-deriving once the right libraries
are known.

**4. Resolve the summary statistic.** The report describes it two ways — Methods says
*coverage minus crossing reads*, while the Results text and the Figure 8 caption both say
*ratio of residuals to crossings*. Both are implemented behind `--statistic`. Run both
and compare tissue rankings.

**5. Close the three mechanical defects** listed under [Unit tests](#unit-tests).

**6. Run the trigger eval suite** and record the ablation delta.

**7. CAGE-seq follow-up.** For candidates that survive all of the above, a transcription
start site falling *within* the intron rather than at the gene's annotated TSS would
support independent transcription rather than simple retention. Not automated here.

## [WIP] Validation documentation

Validating this skill against the original research is ongoing. This section logs what each
pass established, so the next one starts from evidence rather than from scratch.

### Attempt 1 — liver, batch B1 (September 2026)

**What the data looks like.** The source outputs survive as two distinct lineages, which
are not interchangeable:

| lineage | keyed on | introns | available |
|---|---|---|---|
| `*_avecov.txt` | transcript (`ENSMUST…`) | ~121k | polyA B1, NoSelect B1, NuRibo B2 |
| `*_AVE_COV.txt`, `*_INTRON_AVE_COV.txt` | gene (`ENSMUSG…`) | ~139k | polyA B1, NoSelect B1, NuRibo B2, OldRibo B2 |

Alongside them, `*_combined_output.csv` files carry junction-level crossing and splicing
counts, two rows per intron for its left and right junction.

The coverage distributions behave exactly as the paired design predicts. The polyA library
is shallower and sparser — 55.6% of introns at zero coverage, median non-zero 0.16, maximum
1,096 — against the no-select library's 48.3% zero, median non-zero 0.47, maximum 21,075.
That spread in maximum depth is the ten-fold-deeper sequencing of the no-select arm showing
up in the data, which confirms the two files are a genuine matched pair.

**What this pass established.**

- **The no-select arm is NoSelect-B1.** It outperforms both Ribo-Zero preparations
  (OldRibo-B2, NuRibo-B2) in every configuration tested, which identifies which library the
  original comparison drew on.
- **Grouping is per gene, not per transcript.** The gene-keyed lineage scores better at
  every threshold — 27 candidate regions recovered against 20 — so the report's
  "consecutive introns on the same transcript" was in practice per gene. Mixing the two
  lineages across a polyA/no-select pair silently compares different intron sets.
- **Thresholds now come from the data.** The "low polyA" and "non-low no-select" cuts were
  placeholder constants at 1.0. The real distributions put them at 0.0 and 0.1; at the
  placeholder, criterion 3 would have discarded 87% of introns before the comparison ran.
- **Crossing reads track expression, not the free-floating property.** Candidate regions
  carry *more* crossing reads than background (median 3.0 against 0.0), so they cannot act
  as a discriminating filter. This closes off an otherwise plausible approach.

**Where it stands.** The best configuration recovers 27 of the report's 59 candidate
regions while emitting roughly 3,100 candidates. Three explanations for that spread were
tested and eliminated — library pairing, threshold choice, and a crossing-read filter —
which leaves provenance as the open variable. The surviving outputs span polyA, no-select
and two Ribo-Zero preparations across liver, heart and lung, and which combination produced
the published figures is not yet confirmed. Until it is, the comparison runs against an
unidentified target, so 27/59 is a baseline to improve on rather than a measurement of the
method itself.

**Next steps:**

1. Confirm which library combination produced the report's Figure 8 and its 59-candidate
   shortlist.
2. Re-run the comparison on the confirmed pair and re-score recall.
3. Score separately against the higher-confidence subset the appendix marks with asterisks.
4. Re-derive both thresholds from the confirmed negative-control distribution.
5. Settle the summary statistic by comparing tissue rankings under both forms.

## License

[MIT](LICENSE)
