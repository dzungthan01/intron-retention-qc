#!/usr/bin/env python3
"""Identify unambiguous introns from a gene-model BED file + GENCODE annotation.

An intron is "unambiguous" if it shares both splice junctions with every other
annotated intron it overlaps (across isoforms / overlapping genes), and does not
overlap any annotated exon. This excludes cassette exons and alternate 5'/3'
splice-site introns, which would otherwise confound junction-level retention
counts with isoform ambiguity.

v1 draft — re-implemented from the described methodology (the original script for
this specific step was not included in the source report's appendix; the crossing/
splicing counter and coverage calculator scripts were, and are faithfully adapted
in count_junction_signal.py / compute_coverage.py). Review against a known-good
gene model before trusting output on a new species/annotation build.

Input BED columns (tab-separated, no header): chrom, start, end, strand, transcript_id
Output columns: chrom, start, end, strand, intron_id, transcript_id, position_from_3prime
  (position_from_3prime is 1 for the intron closest to the transcript's 3' end)
"""
import argparse
import csv
import sys
from collections import defaultdict


def parse_bed_introns(path):
    introns = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            chrom, start, end, strand, transcript_id = line.split("\t")[:5]
            introns.append(
                {
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "transcript_id": transcript_id,
                }
            )
    return introns


def parse_gencode_exons(path):
    """Minimal GTF exon extractor. Expects standard GTF with an 'exon' feature type."""
    exons = defaultdict(list)  # chrom -> list of (start, end)
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5 or fields[2] != "exon":
                continue
            chrom, start, end = fields[0], int(fields[3]), int(fields[4])
            exons[chrom].append((start, end))
    for chrom in exons:
        exons[chrom].sort()
    return exons


def overlaps_any(chrom, start, end, exons_by_chrom):
    exons = exons_by_chrom.get(chrom, [])
    # linear scan is fine for a one-time filter pass; swap for bisect if annotation is huge
    for ex_start, ex_end in exons:
        if ex_start > end:
            break
        if ex_end >= start:
            return True
    return False


def cluster_and_filter_ambiguous(introns):
    """Within each chromosome, cluster overlapping intron intervals. A cluster is
    unambiguous only if every intron in it has identical (start, end, strand) --
    i.e. there is exactly one distinct coordinate set in the cluster."""
    by_chrom = defaultdict(list)
    for intron in introns:
        by_chrom[intron["chrom"]].append(intron)

    kept = []
    for chrom, chrom_introns in by_chrom.items():
        chrom_introns.sort(key=lambda x: (x["start"], x["end"]))
        cluster = []
        cluster_end = None
        for intron in chrom_introns:
            if cluster and intron["start"] <= cluster_end:
                cluster.append(intron)
                cluster_end = max(cluster_end, intron["end"])
            else:
                kept.extend(_unambiguous_members(cluster))
                cluster = [intron]
                cluster_end = intron["end"]
        kept.extend(_unambiguous_members(cluster))
    return kept


def _unambiguous_members(cluster):
    if not cluster:
        return []
    distinct_coords = {(i["start"], i["end"], i["strand"]) for i in cluster}
    if len(distinct_coords) != 1:
        return []  # ambiguous cluster: alternate splice sites present, drop all
    # de-duplicate identical entries (same intron annotated on multiple transcripts
    # is fine -- keep one representative per transcript_id for downstream pairing)
    seen = set()
    result = []
    for i in cluster:
        key = (i["chrom"], i["start"], i["end"], i["strand"], i["transcript_id"])
        if key not in seen:
            seen.add(key)
            result.append(i)
    return result


def assign_position_from_3prime(introns):
    by_transcript = defaultdict(list)
    for intron in introns:
        by_transcript[intron["transcript_id"]].append(intron)

    for transcript_id, tx_introns in by_transcript.items():
        strand = tx_introns[0]["strand"]
        # 3' end is the highest coordinate on + strand, lowest coordinate on - strand
        tx_introns.sort(key=lambda x: x["start"], reverse=(strand == "+"))
        for rank, intron in enumerate(tx_introns, start=1):
            intron["position_from_3prime"] = rank
    return introns


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genes", required=True, help="BED: chrom start end strand transcript_id")
    ap.add_argument("--gencode", required=True, help="GENCODE GTF (exon features used)")
    ap.add_argument("--out", required=True, help="Output BED-like TSV path")
    ap.add_argument(
        "--chromosomes",
        default=None,
        help="Optional comma-separated allowlist, e.g. chr1,chr2,...,chr19 for mouse autosomes",
    )
    args = ap.parse_args()

    introns = parse_bed_introns(args.genes)
    exons = parse_gencode_exons(args.gencode)

    unambiguous = cluster_and_filter_ambiguous(introns)
    unambiguous = [
        i for i in unambiguous if not overlaps_any(i["chrom"], i["start"], i["end"], exons)
    ]
    unambiguous = assign_position_from_3prime(unambiguous)

    allowed = set(args.chromosomes.split(",")) if args.chromosomes else None

    unambiguous.sort(key=lambda x: (x["chrom"], x["start"]))
    with open(args.out, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        for idx, intron in enumerate(unambiguous):
            if allowed and intron["chrom"] not in allowed:
                continue
            intron_id = f"{intron['chrom']}:{intron['start']}-{intron['end']}({intron['strand']})"
            writer.writerow(
                [
                    intron["chrom"],
                    intron["start"],
                    intron["end"],
                    intron["strand"],
                    intron_id,
                    intron["transcript_id"],
                    intron["position_from_3prime"],
                ]
            )
    print(f"Wrote {len(unambiguous)} unambiguous introns to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
