#!/usr/bin/env python3
"""Apply the three-criterion free-floating-intron test to consecutive intron pairs.

For each pair of consecutive introns on a transcript (A = closer to the 3' end,
B = the next intron toward the 5' end), flag the pair as a candidate if all three
hold:
  1. Low polyA-library coverage at A (near zero)
  2. Real signal at B in polyA, exceeding A by at least --ratio-threshold. The
     floor on B is what carries this criterion when A is zero -- a bare ratio
     against zero passes on any speck of signal, which is not evidence of
     anything (see references/methodology.md)
  3. Non-low no-select coverage at A, and higher at A than at B

v1 draft, implementing the criteria described in the source report's "Identifying
free-floating introns" methodology section (no runnable code for this exact
comparison step was included in that report's appendix; the coverage/crossing
inputs it consumes come from the two scripts that were, adapted upstream in this
package).

Usage:
    python classify_free_floating_introns.py \
        --introns unambiguous_introns.bed \
        --polyA-cov polyA_avecov.tsv --noselect-cov noselect_avecov.tsv \
        --ratio-threshold 0.3 --out candidate_introns.tsv
"""
import argparse
import csv
from collections import defaultdict

# Derived from the liver no-select coverage distribution: 48.3% of introns sit at
# exactly zero, and the non-zero band under 0.1 is a read or two smeared across a
# long intron. Re-derive against your own polyA library before trusting these.
DEFAULT_LOW_POLYA_THRESHOLD = 0.0
DEFAULT_NONLOW_NOSELECT_THRESHOLD = 0.1
DEFAULT_MIN_POLYA_AT_B = 0.1
DEFAULT_RATIO_THRESHOLD = 0.3

# Keeps the ratio defined when A is exactly zero. Criterion 2's floor on B has
# already ruled out noise by the time this is reached.
RATIO_DENOMINATOR_FLOOR = 1e-9


def load_transcript_structure(introns_path):
    """intron_id -> (transcript_id, position_from_3prime)"""
    structure = {}
    by_transcript = defaultdict(list)
    with open(introns_path) as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 7:
                continue
            intron_id, transcript_id, pos = fields[4], fields[5], int(fields[6])
            structure[intron_id] = (transcript_id, pos)
            by_transcript[transcript_id].append((pos, intron_id))
    for transcript_id in by_transcript:
        by_transcript[transcript_id].sort()
    return structure, by_transcript


def load_coverage(path):
    cov = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            cov[row["intron_id"]] = float(row["avg_coverage"])
    return cov


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--introns", required=True)
    ap.add_argument("--polyA-cov", required=True)
    ap.add_argument("--noselect-cov", required=True)
    ap.add_argument("--ratio-threshold", type=float, default=DEFAULT_RATIO_THRESHOLD)
    ap.add_argument("--low-polyA-threshold", type=float, default=DEFAULT_LOW_POLYA_THRESHOLD)
    ap.add_argument("--min-polyA-at-b", type=float, default=DEFAULT_MIN_POLYA_AT_B)
    ap.add_argument(
        "--nonlow-noselect-threshold", type=float, default=DEFAULT_NONLOW_NOSELECT_THRESHOLD
    )
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    _, by_transcript = load_transcript_structure(args.introns)
    polyA = load_coverage(args.polyA_cov)
    noselect = load_coverage(args.noselect_cov)

    candidates = []
    for transcript_id, ordered_introns in by_transcript.items():
        # ordered_introns is [(position_from_3prime, intron_id), ...] sorted ascending,
        # i.e. index 0 is the 3'-most intron (A), index 1 is the next one in (B), etc.
        for idx in range(len(ordered_introns) - 1):
            _, intron_a = ordered_introns[idx]
            _, intron_b = ordered_introns[idx + 1]
            if intron_a not in polyA or intron_b not in polyA:
                continue
            if intron_a not in noselect or intron_b not in noselect:
                continue

            polyA_a, polyA_b = polyA[intron_a], polyA[intron_b]
            noselect_a, noselect_b = noselect[intron_a], noselect[intron_b]

            low_polyA_at_a = polyA_a <= args.low_polyA_threshold
            signal_at_b = polyA_b >= args.min_polyA_at_b
            rises_toward_b = (
                (polyA_b - polyA_a) / max(polyA_a, RATIO_DENOMINATOR_FLOOR)
                >= args.ratio_threshold
            )
            nonlow_noselect_at_a = noselect_a >= args.nonlow_noselect_threshold
            noselect_falls_toward_b = noselect_a > noselect_b

            criterion_1 = low_polyA_at_a
            criterion_2 = signal_at_b and rises_toward_b
            criterion_3 = nonlow_noselect_at_a and noselect_falls_toward_b

            if criterion_1 and criterion_2 and criterion_3:
                candidates.append(
                    {
                        "transcript_id": transcript_id,
                        "intron_a": intron_a,
                        "intron_b": intron_b,
                        "polyA_a": polyA_a,
                        "polyA_b": polyA_b,
                        "noselect_a": noselect_a,
                        "noselect_b": noselect_b,
                    }
                )

    with open(args.out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "transcript_id", "intron_a_3prime", "intron_b_5prime",
                "polyA_a", "polyA_b", "noselect_a", "noselect_b",
            ]
        )
        for c in candidates:
            w.writerow(
                [
                    c["transcript_id"], c["intron_a"], c["intron_b"],
                    c["polyA_a"], c["polyA_b"], c["noselect_a"], c["noselect_b"],
                ]
            )

    print(f"{len(candidates)} candidate intron pairs written to {args.out}")
    print("Remember: treat this as a shortlist for manual review, not a final call.")


if __name__ == "__main__":
    main()
