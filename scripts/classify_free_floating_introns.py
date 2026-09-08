#!/usr/bin/env python3
"""Apply the three-criterion free-floating-intron test to consecutive intron pairs.

For each pair of consecutive introns on a transcript (A = closer to the 3' end,
B = the next intron toward the 5' end), flag the pair as a candidate if all three
hold:
  1. Low polyA-library coverage at A (near zero)
  2. polyA coverage at B exceeds polyA coverage at A by at least --ratio-threshold
     (this is what rules out plain 3' bias -- see references/methodology.md)
  3. Non-low no-select-library coverage at A (A is present, just not polyadenylated)

v1 draft, implementing the criteria described in the source report's "Identifying
free-floating introns" methodology section (no runnable code for this exact
comparison step was included in that report's appendix; the coverage/crossing
inputs it consumes come from the two scripts that were, adapted upstream in this
package). Needs review of the "near zero" / "low" thresholds below against a real
negative-control distribution before trusting output.

Usage:
    python classify_free_floating_introns.py \
        --introns unambiguous_introns.bed \
        --polyA-cov polyA_avecov.tsv --noselect-cov noselect_avecov.tsv \
        --ratio-threshold 0.3 --out candidate_introns.tsv
"""
import argparse
import csv
from collections import defaultdict

# "Near zero" / "non-low" are necessarily dataset-dependent; these are starting
# points based on the source study and should be re-derived from your own data
# (e.g. the Nth percentile of coverage at introns with no crossing-read support).
DEFAULT_LOW_POLYA_THRESHOLD = 1.0
DEFAULT_NONLOW_NOSELECT_THRESHOLD = 1.0


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
    ap.add_argument("--ratio-threshold", type=float, default=0.3)
    ap.add_argument("--low-polyA-threshold", type=float, default=DEFAULT_LOW_POLYA_THRESHOLD)
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
            if intron_a not in polyA or intron_b not in polyA or intron_a not in noselect:
                continue

            polyA_a, polyA_b = polyA[intron_a], polyA[intron_b]
            noselect_a = noselect[intron_a]

            criterion_1_low_polyA_at_a = polyA_a <= args.low_polyA_threshold
            denom = polyA_a if polyA_a > 0 else 1e-9
            criterion_2_b_exceeds_a = (polyA_b - polyA_a) / denom >= args.ratio_threshold
            criterion_3_nonlow_noselect_at_a = noselect_a >= args.nonlow_noselect_threshold

            if criterion_1_low_polyA_at_a and criterion_2_b_exceeds_a and criterion_3_nonlow_noselect_at_a:
                candidates.append(
                    {
                        "transcript_id": transcript_id,
                        "intron_a": intron_a,
                        "intron_b": intron_b,
                        "polyA_a": polyA_a,
                        "polyA_b": polyA_b,
                        "noselect_a": noselect_a,
                    }
                )

    with open(args.out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            ["transcript_id", "intron_a_3prime", "intron_b_5prime", "polyA_a", "polyA_b", "noselect_a"]
        )
        for c in candidates:
            w.writerow(
                [c["transcript_id"], c["intron_a"], c["intron_b"], c["polyA_a"], c["polyA_b"], c["noselect_a"]]
            )

    print(f"{len(candidates)} candidate intron pairs written to {args.out}")
    print("Remember: treat this as a shortlist for manual review, not a final call.")


if __name__ == "__main__":
    main()
