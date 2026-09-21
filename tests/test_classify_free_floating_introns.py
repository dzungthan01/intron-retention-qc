"""Stage 2: the three-criterion test, on coverage tables with known answers.

A = the 3'-most intron of a consecutive pair (position_from_3prime 1),
B = the next intron toward the 5' end (position 2).

Cases are (polyA_a, polyA_b, noselect_a, noselect_b).
"""
from conftest import read_tsv, run_script, write_coverage_table, write_tsv


def build_pairs(tmp_path, cases):
    intron_rows = []
    polyA_coverage = {}
    noselect_coverage = {}

    for index, (transcript_id, values) in enumerate(sorted(cases.items())):
        polyA_a, polyA_b, noselect_a, noselect_b = values
        start = 1000 + index * 1000
        intron_a = f"{transcript_id}.A"
        intron_b = f"{transcript_id}.B"
        intron_rows.append(["chr1", start, start + 100, "+", intron_a, transcript_id, 1])
        intron_rows.append(["chr1", start + 200, start + 300, "+", intron_b, transcript_id, 2])
        polyA_coverage[intron_a] = polyA_a
        polyA_coverage[intron_b] = polyA_b
        noselect_coverage[intron_a] = noselect_a
        noselect_coverage[intron_b] = noselect_b

    introns = write_tsv(tmp_path / "introns.bed", intron_rows)
    polyA = write_coverage_table(tmp_path / "polyA_avecov.tsv", polyA_coverage)
    noselect = write_coverage_table(tmp_path / "noselect_avecov.tsv", noselect_coverage)
    return introns, polyA, noselect


def classify(tmp_path, cases, **overrides):
    introns, polyA, noselect = build_pairs(tmp_path, cases)
    output = tmp_path / "candidates.tsv"
    arguments = [
        "--introns", introns,
        "--polyA-cov", polyA,
        "--noselect-cov", noselect,
        "--out", output,
    ]
    for flag, value in overrides.items():
        arguments += [f"--{flag.replace('_', '-')}", value]
    run_script("classify_free_floating_introns.py", *arguments)
    return {row["transcript_id"] for row in read_tsv(output)}


def test_true_positive_is_flagged(tmp_path):
    # A silent in polyA, real signal at B, A present in no-select and above B.
    flagged = classify(tmp_path, {"TX_TRUE": (0.0, 5.0, 10.0, 2.0)})
    assert flagged == {"TX_TRUE"}


def test_three_prime_bias_is_rejected(tmp_path):
    # Signal falling toward the 5' end is what 3' bias predicts. Criterion 1 is
    # relaxed here so criterion 2 is the only thing under test.
    flagged = classify(
        tmp_path,
        {"TX_BIAS": (0.5, 0.1, 10.0, 2.0)},
        low_polyA_threshold=1.0,
    )
    assert flagged == set()


def test_negligible_signal_at_b_is_rejected(tmp_path):
    # B is indistinguishable from zero, so the rise from A is not evidence of
    # anything. The floor on B is what catches this; a bare ratio against zero
    # would pass it.
    flagged = classify(tmp_path, {"TX_NOISE": (0.0, 0.0001, 10.0, 2.0)})
    assert flagged == set()


def test_signal_floor_at_b_is_enforced(tmp_path):
    # Just under and just over the default 0.1 floor.
    assert classify(tmp_path, {"TX_UNDER": (0.0, 0.09, 10.0, 2.0)}) == set()
    assert classify(tmp_path, {"TX_OVER": (0.0, 0.11, 10.0, 2.0)}) == {"TX_OVER"}


def test_absent_in_no_select_is_rejected(tmp_path):
    # Nothing at A in either library: criterion 3's floor must reject it.
    flagged = classify(tmp_path, {"TX_ABSENT": (0.0, 5.0, 0.05, 0.0)})
    assert flagged == set()


def test_no_select_must_be_higher_at_a_than_at_b(tmp_path):
    # The report's criterion 3 is relational: no-select at A must surpass B.
    flagged = classify(tmp_path, {"TX_NS": (0.0, 5.0, 10.0, 40.0)})
    assert flagged == set()


def test_high_polya_at_a_is_rejected(tmp_path):
    # A is not silent in polyA, so it is ordinary retention: criterion 1 rejects it.
    flagged = classify(tmp_path, {"TX_RETAINED": (5.0, 9.0, 10.0, 2.0)})
    assert flagged == set()


def test_cases_are_judged_independently(tmp_path):
    flagged = classify(
        tmp_path,
        {
            "TX_TRUE": (0.0, 5.0, 10.0, 2.0),
            "TX_NOISE": (0.0, 0.0001, 10.0, 2.0),
            "TX_ABSENT": (0.0, 5.0, 0.05, 0.0),
            "TX_NS": (0.0, 5.0, 10.0, 40.0),
        },
    )
    assert flagged == {"TX_TRUE"}


def test_candidate_row_reports_both_no_select_values(tmp_path):
    introns, polyA, noselect = build_pairs(tmp_path, {"TX_TRUE": (0.0, 5.0, 10.0, 2.0)})
    output = tmp_path / "candidates.tsv"
    run_script(
        "classify_free_floating_introns.py",
        "--introns", introns,
        "--polyA-cov", polyA,
        "--noselect-cov", noselect,
        "--out", output,
    )
    row = read_tsv(output)[0]
    assert float(row["noselect_a"]) == 10.0
    assert float(row["noselect_b"]) == 2.0
