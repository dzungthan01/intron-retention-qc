"""Stage 2: per-intron average coverage against hand-computed bedGraph fixtures."""
from conftest import read_tsv, run_script, write_tsv


def compute(tmp_path, intron_rows, forward_rows, reverse_rows):
    introns = write_tsv(tmp_path / "introns.bed", intron_rows)
    forward = write_tsv(tmp_path / "fwd.bedgraph", forward_rows)
    reverse = write_tsv(tmp_path / "rev.bedgraph", reverse_rows)
    output = tmp_path / "avecov.tsv"
    run_script(
        "compute_coverage.py",
        "--introns", introns,
        "--cov-fwd", forward,
        "--cov-rev", reverse,
        "--out", output,
    )
    return {row["intron_id"]: float(row["avg_coverage"]) for row in read_tsv(output)}


def test_flat_coverage_over_whole_intron(tmp_path):
    # 10 bases inclusive at depth 2 -> total 20, average 2.0.
    coverage = compute(
        tmp_path,
        [["chr1", 100, 109, "+", "I1"]],
        [["chr1", 100, 109, 2.0]],
        [],
    )
    assert coverage == {"I1": 2.0}


def test_partial_overlap_is_clipped_to_the_intron(tmp_path):
    # Only 105-109 of the interval falls inside: 5 bases at depth 4 -> 20/10.
    coverage = compute(
        tmp_path,
        [["chr1", 100, 109, "+", "I1"]],
        [["chr1", 105, 120, 4.0]],
        [],
    )
    assert coverage == {"I1": 2.0}


def test_multiple_intervals_are_summed(tmp_path):
    coverage = compute(
        tmp_path,
        [["chr1", 100, 109, "+", "I1"]],
        [["chr1", 100, 104, 1.0], ["chr1", 105, 109, 3.0]],
        [],
    )
    assert coverage == {"I1": 2.0}


def test_strand_selects_the_matching_bedgraph(tmp_path):
    coverage = compute(
        tmp_path,
        [["chr1", 100, 109, "+", "PLUS"], ["chr1", 200, 209, "-", "MINUS"]],
        [["chr1", 100, 109, 2.0], ["chr1", 200, 209, 99.0]],
        [["chr1", 200, 209, 3.0]],
    )
    assert coverage == {"PLUS": 2.0, "MINUS": 3.0}


def test_cursor_advances_across_consecutive_introns(tmp_path):
    coverage = compute(
        tmp_path,
        [
            ["chr1", 100, 109, "+", "I1"],
            ["chr1", 200, 209, "+", "I2"],
            ["chr1", 300, 309, "+", "I3"],
        ],
        [["chr1", 100, 109, 2.0], ["chr1", 200, 209, 4.0], ["chr1", 300, 309, 6.0]],
        [],
    )
    assert coverage == {"I1": 2.0, "I2": 4.0, "I3": 6.0}


def test_one_interval_spanning_every_intron(tmp_path):
    # Exercises the cursor back-off: a single interval overlaps all three introns.
    coverage = compute(
        tmp_path,
        [
            ["chr1", 100, 109, "+", "I1"],
            ["chr1", 200, 209, "+", "I2"],
            ["chr1", 300, 309, "+", "I3"],
        ],
        [["chr1", 50, 400, 1.0]],
        [],
    )
    assert coverage == {"I1": 1.0, "I2": 1.0, "I3": 1.0}


def test_track_header_is_skipped(tmp_path):
    introns = write_tsv(tmp_path / "introns.bed", [["chr1", 100, 109, "+", "I1"]])
    forward = tmp_path / "fwd.bedgraph"
    forward.write_text('track type=bedGraph name="fwd"\nchr1\t100\t109\t2.0\n')
    reverse = write_tsv(tmp_path / "rev.bedgraph", [])
    output = tmp_path / "avecov.tsv"
    run_script(
        "compute_coverage.py",
        "--introns", introns,
        "--cov-fwd", forward,
        "--cov-rev", reverse,
        "--out", output,
    )
    assert {row["intron_id"]: float(row["avg_coverage"]) for row in read_tsv(output)} == {"I1": 2.0}


def test_intervals_are_treated_as_closed_like_the_original_script(tmp_path):
    # Pins the convention rather than endorsing it: the report's appendix script also
    # uses end - start + 1, so the reimplementation is faithful. bedGraph is normally
    # half-open, which would make both the length and the overlap one base too long.
    rows = read_tsv(
        _run_single(tmp_path, ["chr1", 100, 109, "+", "I1"], [["chr1", 100, 109, 2.0]])
    )
    assert int(rows[0]["length"]) == 10
    assert float(rows[0]["total_coverage"]) == 20.0


def _run_single(tmp_path, intron_row, forward_rows):
    introns = write_tsv(tmp_path / "introns.bed", [intron_row])
    forward = write_tsv(tmp_path / "fwd.bedgraph", forward_rows)
    reverse = write_tsv(tmp_path / "rev.bedgraph", [])
    output = tmp_path / "avecov.tsv"
    run_script(
        "compute_coverage.py",
        "--introns", introns,
        "--cov-fwd", forward,
        "--cov-rev", reverse,
        "--out", output,
    )
    return output
