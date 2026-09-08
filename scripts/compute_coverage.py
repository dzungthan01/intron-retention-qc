#!/usr/bin/env python3
"""Compute average coverage per intron from forward + reverse strand bedGraph files.

v1 draft, adapted from the original average-coverage script in the source report.
Rewritten as an O(n log n) two-pointer sweep over sorted per-chromosome intervals
instead of the original's per-intron rescan of the coverage file (which was
correct but quadratic). Assumes both the intron list and the bedGraph files are
sorted by (chrom, start) -- true for standard bedGraph output; re-sort upstream
if yours isn't.

Usage:
    python compute_coverage.py --introns unambiguous_introns.bed \
        --cov-fwd sample.fwd.bedgraph --cov-rev sample.rev.bedgraph \
        --out avecov.tsv
"""
import argparse
import csv
from collections import defaultdict


def load_introns(path):
    introns = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            chrom, start, end, strand, intron_id = fields[0], int(fields[1]), int(fields[2]), fields[3], fields[4]
            introns.append({"chrom": chrom, "start": start, "end": end, "strand": strand, "id": intron_id})
    introns.sort(key=lambda i: (i["chrom"], i["start"]))
    return introns


def load_bedgraph(path):
    by_chrom = defaultdict(list)
    with open(path) as f:
        first = f.readline()
        if not first.startswith("track"):
            f.seek(0)
        for line in f:
            fields = line.strip().split()
            if len(fields) < 4:
                continue
            chrom, start, end, value = fields[0], int(fields[1]), int(fields[2]), float(fields[3])
            by_chrom[chrom].append((start, end, value))
    for chrom in by_chrom:
        by_chrom[chrom].sort()
    return by_chrom


def total_coverage_in_range(cov_intervals, start, end, cursor_state):
    """Sums coverage*overlap_length for cov_intervals overlapping [start, end],
    advancing a shared cursor since both introns and cov_intervals are sorted."""
    idx = cursor_state["idx"]
    total = 0.0
    # back off in case a previous interval still overlaps this (later-starting) intron
    while idx > 0 and cov_intervals[idx - 1][1] >= start:
        idx -= 1
    while idx < len(cov_intervals) and cov_intervals[idx][0] <= end:
        cov_start, cov_end, value = cov_intervals[idx]
        if cov_end >= start:
            overlap_start, overlap_end = max(start, cov_start), min(end, cov_end)
            if overlap_start <= overlap_end:
                total += (overlap_end - overlap_start + 1) * value
        if cov_end < end:
            idx += 1
        else:
            break
    cursor_state["idx"] = idx
    return total


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--introns", required=True)
    ap.add_argument("--cov-fwd", required=True)
    ap.add_argument("--cov-rev", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    introns = load_introns(args.introns)
    fwd = load_bedgraph(args.cov_fwd)
    rev = load_bedgraph(args.cov_rev)

    cursors = defaultdict(lambda: {"idx": 0})
    rows = []
    for intron in introns:
        cov_source = fwd if intron["strand"] == "+" else rev
        intervals = cov_source.get(intron["chrom"], [])
        cursor = cursors[(intron["chrom"], intron["strand"])]
        total = total_coverage_in_range(intervals, intron["start"], intron["end"], cursor)
        length = intron["end"] - intron["start"] + 1
        avg_cov = total / length if length > 0 else 0.0
        rows.append(
            [intron["id"], intron["chrom"], intron["start"], intron["end"], intron["strand"], total, length, avg_cov]
        )

    with open(args.out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["intron_id", "chrom", "start", "end", "strand", "total_coverage", "length", "avg_coverage"])
        w.writerows(rows)

    print(f"Wrote average coverage for {len(rows)} introns to {args.out}")


if __name__ == "__main__":
    main()
