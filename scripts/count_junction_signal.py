#!/usr/bin/env python3
"""Count junction-crossing and junction-splicing reads per unambiguous intron.

For each intron, a "crossing" read is one whose CIGAR contains an M (match) block
that spans across the junction position without a gap -- evidence the intron is
still attached to (i.e. retained in) the transcript. A "splicing" read is one
whose CIGAR contains an N (skip/intron) block whose boundary lands exactly on the
junction -- evidence of a normal, correctly spliced transcript.

v1 draft, adapted and cleaned up from the original crossing/splicing counter
script in the source report (same core algorithm: CIGAR-block walking + an
overlap cutoff to avoid crediting reads that just barely touch a junction).
Rewritten to use bisect instead of a hand-rolled binary search, and to avoid the
original's numeric sentinel encoding for chrX/chrY. Needs validation against a
real SAM file with known-truth junctions before trusting output at scale.

Usage:
    python count_junction_signal.py --introns unambiguous_introns.bed --sam sample.sam \
        --out-crossings crossings.tsv --out-splicings splicings.tsv
"""
import argparse
import bisect
import csv
import re
import sys
from collections import defaultdict

OVERLAP_CUTOFF = 8  # bp of clearance required on each side of a junction to count a crossing


def load_introns(path):
    """Returns junctions indexed by (chrom, strand, boundary) -> sorted position list,
    plus a lookup from (chrom, strand, position) -> intron_id for each boundary type."""
    starts = defaultdict(list)  # (chrom, strand) -> sorted [start positions]
    ends = defaultdict(list)  # (chrom, strand) -> sorted [end positions]
    start_owner = {}  # (chrom, strand, start) -> intron_id
    end_owner = {}  # (chrom, strand, end) -> intron_id

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            chrom, start, end, strand, intron_id = fields[0], int(fields[1]), int(fields[2]), fields[3], fields[4]
            key = (chrom, strand)
            starts[key].append(start)
            ends[key].append(end)
            start_owner[(chrom, strand, start)] = intron_id
            end_owner[(chrom, strand, end)] = intron_id

    for key in starts:
        starts[key].sort()
    for key in ends:
        ends[key].sort()

    return starts, ends, start_owner, end_owner


def parse_sam_line(line):
    """Returns (chrom, pos, strand, blocks) or None if the line should be skipped.
    blocks is a list of (op, length) for M/N operations only (I/S/D collapsed as in
    the source methodology: deletions treated as matches, insertions/soft-clips dropped)."""
    if line.startswith("@"):
        return None
    fields = line.split("\t")
    if len(fields) < 11:
        return None
    flag = int(fields[1])
    reverse = bool(flag & 16)
    first_in_pair = bool(flag & 64)
    nh_match = re.search(r"NH:i:(\d+)", line)
    if nh_match and int(nh_match.group(1)) > 1:
        return None  # skip multi-mappers, as in the source methodology

    strand = ("-" if reverse else "+") if first_in_pair else ("+" if reverse else "-")
    chrom = fields[2]
    pos = int(fields[3])
    cigar = fields[5]

    cigar = re.sub(r"\d+I", "", cigar)  # drop insertions
    cigar = cigar.replace("D", "M")  # treat deletions as matches
    cigar = re.sub(r"\d+S", "", cigar)  # drop soft clips
    blocks = [(int(n), op) for n, op in re.findall(r"(\d+)([MN])", cigar)]
    return chrom, pos, strand, blocks


def walk_read(chrom, pos, strand, blocks, starts, ends, start_owner, end_owner):
    """Returns (crossing_events, splicing_events) as lists of intron_ids."""
    crossings, splicings = [], []
    key = (chrom, strand)
    chrom_starts, chrom_ends = starts.get(key, []), ends.get(key, [])

    cursor = pos
    for length, op in blocks:
        block_start, block_end = cursor, cursor + length - 1
        if op == "M":
            for owner_map, positions in ((start_owner, chrom_starts), (end_owner, chrom_ends)):
                lo = bisect.bisect_left(positions, block_start)
                hi = bisect.bisect_right(positions, block_end)
                for junction_pos in positions[lo:hi]:
                    if (
                        junction_pos - block_start > OVERLAP_CUTOFF
                        and block_end - junction_pos > OVERLAP_CUTOFF
                    ):
                        crossings.append(owner_map[(chrom, strand, junction_pos)])
        elif op == "N":
            # a skip (N) whose boundaries land exactly on annotated junctions = correct splicing
            if (chrom, strand, block_start) in start_owner:
                splicings.append(start_owner[(chrom, strand, block_start)])
            if (chrom, strand, block_end) in end_owner:
                splicings.append(end_owner[(chrom, strand, block_end)])
        cursor += length

    return crossings, splicings


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--introns", required=True)
    ap.add_argument("--sam", required=True)
    ap.add_argument("--out-crossings", required=True)
    ap.add_argument("--out-splicings", required=True)
    args = ap.parse_args()

    starts, ends, start_owner, end_owner = load_introns(args.introns)
    crossing_counts = defaultdict(int)
    splicing_counts = defaultdict(int)

    with open(args.sam) as sam:
        lines = [l for l in sam if not l.startswith("@")]

    # reads are consumed in mate pairs, as in the source methodology; corrects for
    # double-counting a junction event seen on both reads of an overlapping pair
    i = 0
    n_reads = 0
    while i < len(lines) - 1:
        r1, r2 = parse_sam_line(lines[i]), parse_sam_line(lines[i + 1])
        i += 2
        if r1 is None or r2 is None:
            continue
        n_reads += 2
        c1, s1 = walk_read(*r1, starts, ends, start_owner, end_owner)
        c2, s2 = walk_read(*r2, starts, ends, start_owner, end_owner)
        # de-duplicate crossing/splicing events shared by both mates of the pair
        for intron_id in set(c1) | set(c2):
            crossing_counts[intron_id] += 1 if (intron_id in c1) != (intron_id in c2) else 1
        for intron_id in set(s1) | set(s2):
            splicing_counts[intron_id] += 1 if (intron_id in s1) != (intron_id in s2) else 1
        if n_reads % 1_000_000 == 0:
            print(f"Processed {n_reads} reads", file=sys.stderr)

    with open(args.out_crossings, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        for intron_id, count in sorted(crossing_counts.items()):
            w.writerow([intron_id, count])

    with open(args.out_splicings, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        for intron_id, count in sorted(splicing_counts.items()):
            w.writerow([intron_id, count])

    print(
        f"{len(crossing_counts)} introns with crossing signal, "
        f"{len(splicing_counts)} with splicing signal",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
