"""Stage 2: the unambiguous-intron filter and 3'-position ranking.

Exons sit at 1-100 and 600-700 on every chromosome, so the 200-300 and 400-500
introns below are clear of them unless a test says otherwise.
"""
import pytest

from conftest import run_script, write_tsv

CHROMOSOMES = ["chr1", "chr2", "chr3", "chr4", "chr5"]


def gtf_exons(tmp_path):
    rows = []
    for chromosome in CHROMOSOMES:
        for start, end in [(1, 100), (600, 700)]:
            rows.append([chromosome, "TEST", "exon", start, end, ".", "+", ".", 'gene_id "g";'])
    return write_tsv(tmp_path / "annotation.gtf", rows)


def identify(tmp_path, bed_rows, extra_arguments=()):
    genes = write_tsv(tmp_path / "genes.bed", bed_rows)
    output = tmp_path / "unambiguous.bed"
    run_script(
        "identify_unambiguous_introns.py",
        "--genes", genes,
        "--gencode", gtf_exons(tmp_path),
        "--out", output,
        *extra_arguments,
    )
    return [line.split("\t") for line in output.read_text().splitlines() if line.strip()]


def test_isolated_introns_are_kept(tmp_path):
    rows = identify(tmp_path, [["chr1", 200, 300, "+", "tx1"], ["chr1", 400, 500, "+", "tx1"]])
    assert {row[4] for row in rows} == {"chr1:200-300(+)", "chr1:400-500(+)"}


def test_overlapping_introns_with_different_coordinates_are_dropped(tmp_path):
    # An alternate splice site: the cluster holds two distinct coordinate sets.
    rows = identify(tmp_path, [["chr3", 200, 300, "+", "txA"], ["chr3", 250, 350, "+", "txB"]])
    assert rows == []


def test_identical_intron_on_two_transcripts_is_kept_for_each(tmp_path):
    rows = identify(tmp_path, [["chr4", 200, 300, "+", "txC"], ["chr4", 200, 300, "+", "txD"]])
    assert sorted(row[5] for row in rows) == ["txC", "txD"]


def test_intron_overlapping_an_exon_is_dropped(tmp_path):
    rows = identify(tmp_path, [["chr5", 50, 150, "+", "txE"]])
    assert rows == []


def test_position_from_three_prime_on_the_plus_strand(tmp_path):
    rows = identify(tmp_path, [["chr1", 200, 300, "+", "tx1"], ["chr1", 400, 500, "+", "tx1"]])
    ranks = {row[4]: int(row[6]) for row in rows}
    # On +, the 3' end is the highest coordinate, so 400-500 ranks first.
    assert ranks == {"chr1:400-500(+)": 1, "chr1:200-300(+)": 2}


def test_position_from_three_prime_on_the_minus_strand(tmp_path):
    rows = identify(tmp_path, [["chr2", 200, 300, "-", "tx2"], ["chr2", 400, 500, "-", "tx2"]])
    ranks = {row[4]: int(row[6]) for row in rows}
    # On -, the 3' end is the lowest coordinate, so 200-300 ranks first.
    assert ranks == {"chr2:200-300(-)": 1, "chr2:400-500(-)": 2}


def test_chromosome_allowlist_filters_output(tmp_path):
    rows = identify(
        tmp_path,
        [["chr1", 200, 300, "+", "tx1"], ["chr2", 200, 300, "-", "tx2"]],
        extra_arguments=("--chromosomes", "chr1"),
    )
    assert {row[0] for row in rows} == {"chr1"}


@pytest.mark.xfail(
    strict=True,
    reason="identify_unambiguous_introns.py:172 reports len(unambiguous), which counts "
    "introns the --chromosomes allowlist excluded at write time. The number printed "
    "does not match the rows written, which misleads the 'did step 1 return a "
    "plausible count' sanity check.",
)
def test_reported_count_matches_rows_written(tmp_path):
    genes = write_tsv(
        tmp_path / "genes.bed",
        [["chr1", 200, 300, "+", "tx1"], ["chr2", 200, 300, "-", "tx2"]],
    )
    output = tmp_path / "unambiguous.bed"
    completed = run_script(
        "identify_unambiguous_introns.py",
        "--genes", genes,
        "--gencode", gtf_exons(tmp_path),
        "--out", output,
        "--chromosomes", "chr1",
    )
    written = len([line for line in output.read_text().splitlines() if line.strip()])
    assert f"Wrote {written} unambiguous introns" in completed.stderr
