#!/usr/bin/env python3
"""Summarize intron-retention signal per sample/tissue for cross-sample comparison.

Per exon-intron junction, residual = avg_coverage - crossing_reads. Two summary
forms are supported, because the source report describes both:

  difference  sum(residual) over junctions where residual > 0
              -- the report's Methods section
  ratio       sum(residual / (crossing_reads + pseudocount)) over the same
              -- the report's Results text and Figure 8 caption

Higher = more retention/free-floating signal relative to properly-spliced reads.
Both are normalized by the sample's read depth.

The pseudocount is what keeps `ratio` defined: roughly half of all introns have
zero crossing reads, and those are the most informative cases (coverage with no
evidence of exon attachment), so they must not be dropped. At the default
pseudocount of 1 a zero-crossing junction contributes exactly its coverage, which
is what `difference` gives it too -- the two forms agree precisely where the
ratio would otherwise divide by zero.

Intended for comparing multiple tissues or conditions on one axis, not for
declaring any single sample "positive" or "negative" in isolation -- include a
polyA-only sample as a negative control (expected to score near zero) if you have
one, as in the source study.

v1 draft, adapted from the summary-statistic description in the source report
(this specific aggregation script was not included in that report's appendix;
the crossing-count and coverage inputs it consumes come from the two scripts
that were, adapted upstream in this package).

Usage:
    python tissue_survey_summary.py \
        --crossings liver_crossings.tsv kidney_crossings.tsv \
        --avecov liver_avecov.tsv kidney_avecov.tsv \
        --statistic ratio --depths liver=32000000 kidney=31500000 \
        --out tissue_summary.tsv

--crossings and --avecov must be given in the same sample order. --depths is
optional (sample=read_count pairs); if omitted, raw (unnormalized) sums are
reported and depth normalization is left to the caller.
"""
import argparse
import csv
import os
from collections import defaultdict

DEFAULT_PSEUDOCOUNT = 1.0


def sample_name_from_path(path):
    return os.path.splitext(os.path.basename(path))[0]


def load_counts(path):
    counts = {}
    with open(path) as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 2:
                continue
            counts[fields[0]] = float(fields[1])
    return counts


def load_avecov(path):
    cov = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            cov[row["intron_id"]] = float(row["avg_coverage"])
    return cov


def junction_statistic(coverage, crossings, statistic, pseudocount):
    residual = coverage - crossings
    if statistic == "difference":
        return residual
    return residual / (crossings + pseudocount)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--crossings", nargs="+", required=True)
    ap.add_argument("--avecov", nargs="+", required=True)
    ap.add_argument("--statistic", choices=("difference", "ratio"), default="ratio")
    ap.add_argument("--pseudocount", type=float, default=DEFAULT_PSEUDOCOUNT)
    ap.add_argument("--depths", nargs="*", default=[], help="sample=read_count pairs, optional")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if len(args.crossings) != len(args.avecov):
        raise SystemExit("--crossings and --avecov must have the same number of files, in matching sample order")
    if args.pseudocount <= 0:
        raise SystemExit("--pseudocount must be positive, or the ratio is undefined at zero crossings")

    depth_by_sample = {}
    for pair in args.depths:
        name, value = pair.split("=")
        depth_by_sample[name] = float(value)

    rows = []
    for crossings_path, avecov_path in zip(args.crossings, args.avecov):
        sample = sample_name_from_path(crossings_path)
        crossings = load_counts(crossings_path)
        avecov = load_avecov(avecov_path)

        positive_sum = 0.0
        zero_crossing_junctions = 0
        for intron_id, cov in avecov.items():
            crossing_reads = crossings.get(intron_id, 0.0)
            if crossing_reads == 0:
                zero_crossing_junctions += 1
            value = junction_statistic(cov, crossing_reads, args.statistic, args.pseudocount)
            if value > 0:
                positive_sum += value

        depth = depth_by_sample.get(sample)
        normalized = positive_sum / depth if depth else positive_sum
        rows.append((sample, args.statistic, positive_sum, depth if depth else "", normalized))
        print(f"{sample}: {zero_crossing_junctions} of {len(avecov)} junctions have zero crossing reads")

    with open(args.out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample", "statistic", "raw_sum_positive_residuals", "read_depth", "normalized_statistic"])
        w.writerows(rows)

    print(f"Wrote summary for {len(rows)} samples to {args.out}")
    if not depth_by_sample:
        print("No --depths given: normalized_statistic == raw sum. Pass --depths to normalize.")


if __name__ == "__main__":
    main()
