"""Stage 2: crossing/splicing counts on a hand-built SAM.

The intron under test is chr1:1000-1100 on the + strand, so 1000 is the 5' junction
and 1100 the 3' junction. OVERLAP_CUTOFF is 8, and the comparison is strict.
"""
import pytest

from conftest import run_script, write_tsv

INTRON = ["chr1", 1000, 1100, "+", "I1"]

# A proper FR pair: read 1 forward, read 2 reverse. Both resolve to the + strand
# under the script's flag logic, which is what the original appendix script assumed
# when it applied read 1's strand to both mates.
FLAG_READ1 = 99
FLAG_READ2 = 147


def sam_line(qname, flag, position, cigar, alignments=1):
    return "\t".join(
        [qname, str(flag), "chr1", str(position), "255", cigar, "=", "0", "0", "*", "*",
         f"NH:i:{alignments}"]
    )


def count(tmp_path, sam_lines, introns=(INTRON,)):
    intron_file = write_tsv(tmp_path / "introns.bed", list(introns))
    sam = tmp_path / "reads.sam"
    sam.write_text("@HD\tVN:1.6\tSO:queryname\n" + "\n".join(sam_lines) + "\n")
    crossings = tmp_path / "crossings.tsv"
    splicings = tmp_path / "splicings.tsv"
    run_script(
        "count_junction_signal.py",
        "--introns", intron_file,
        "--sam", sam,
        "--out-crossings", crossings,
        "--out-splicings", splicings,
    )
    return _read_counts(crossings), _read_counts(splicings)


def _read_counts(path):
    counts = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        intron_id, value = line.split("\t")
        counts[intron_id] = int(value)
    return counts


def test_read_crossing_a_junction_is_counted(tmp_path):
    # 980-1039 clears the 1000 junction by 20bp on the left and 39bp on the right.
    crossings, _ = count(
        tmp_path,
        [
            sam_line("pair1", FLAG_READ1, 980, "60M"),
            sam_line("pair1", FLAG_READ2, 5000, "50M"),
        ],
    )
    assert crossings == {"I1": 1}


def test_clearance_of_exactly_the_cutoff_is_not_counted(tmp_path):
    # 992-1051 leaves exactly 8bp left of the junction; the comparison is strict.
    crossings, _ = count(
        tmp_path,
        [
            sam_line("pair1", FLAG_READ1, 992, "60M"),
            sam_line("pair1", FLAG_READ2, 5000, "50M"),
        ],
    )
    assert crossings == {}


def test_spliced_read_is_counted_as_splicing_not_crossing(tmp_path):
    # 60M lands on 940-999, then the 101N skip spans exactly 1000-1100.
    crossings, splicings = count(
        tmp_path,
        [
            sam_line("pair1", FLAG_READ1, 940, "60M101N40M"),
            sam_line("pair1", FLAG_READ2, 5000, "50M"),
        ],
    )
    assert crossings == {}
    assert splicings == {"I1": 1}


def test_multi_mapped_pair_is_skipped(tmp_path):
    crossings, _ = count(
        tmp_path,
        [
            sam_line("pair1", FLAG_READ1, 980, "60M", alignments=2),
            sam_line("pair1", FLAG_READ2, 5000, "50M"),
        ],
    )
    assert crossings == {}


def test_both_mates_crossing_counts_once(tmp_path):
    # The de-duplication the original handled with fix_dups.
    crossings, _ = count(
        tmp_path,
        [
            sam_line("pair1", FLAG_READ1, 980, "60M"),
            sam_line("pair1", FLAG_READ2, 985, "60M"),
        ],
    )
    assert crossings == {"I1": 1}


def test_crossing_on_the_second_mate_is_counted(tmp_path):
    # Mate 2 must resolve to the same strand as mate 1, or its evidence is lost.
    crossings, _ = count(
        tmp_path,
        [
            sam_line("pair1", FLAG_READ1, 5000, "50M"),
            sam_line("pair1", FLAG_READ2, 980, "60M"),
        ],
    )
    assert crossings == {"I1": 1}


def test_separate_pairs_accumulate(tmp_path):
    crossings, _ = count(
        tmp_path,
        [
            sam_line("pair1", FLAG_READ1, 980, "60M"),
            sam_line("pair1", FLAG_READ2, 5000, "50M"),
            sam_line("pair2", FLAG_READ1, 981, "60M"),
            sam_line("pair2", FLAG_READ2, 5000, "50M"),
        ],
    )
    assert crossings == {"I1": 2}


@pytest.mark.xfail(
    strict=True,
    reason="count_junction_signal.py:137-138 takes consecutive lines as mates with no "
    "QNAME check, so a coordinate-sorted SAM silently pairs unrelated reads. STAR "
    "emits coordinate-sorted BAM by default; the script should refuse it rather than "
    "produce plausible counts.",
)
def test_unpaired_neighbours_are_rejected(tmp_path):
    # Two reads from different templates sitting next to each other.
    _, _ = count(
        tmp_path,
        [
            sam_line("readA", FLAG_READ1, 980, "60M"),
            sam_line("readB", FLAG_READ1, 5000, "50M"),
        ],
    )
    pytest.fail("expected the script to reject a SAM that is not name-sorted")
