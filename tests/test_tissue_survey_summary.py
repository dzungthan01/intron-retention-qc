"""Stage 2: the cross-tissue summary statistic, in both forms, and normalization."""
import subprocess
import sys

import pytest

from conftest import SCRIPTS_DIR, read_tsv, run_script, write_coverage_table, write_tsv


def summarize(tmp_path, samples, depths=(), statistic=None, pseudocount=None):
    """samples: {stem: {intron_id: (avg_coverage, crossings)}}."""
    crossing_files = []
    coverage_files = []
    for stem, introns in sorted(samples.items()):
        crossing_files.append(
            write_tsv(
                tmp_path / f"{stem}_crossings.tsv",
                [[intron_id, crossings] for intron_id, (_, crossings) in sorted(introns.items())],
            )
        )
        coverage_files.append(
            write_coverage_table(
                tmp_path / f"{stem}_avecov.tsv",
                {intron_id: coverage for intron_id, (coverage, _) in introns.items()},
            )
        )

    output = tmp_path / "summary.tsv"
    arguments = ["--crossings", *crossing_files, "--avecov", *coverage_files]
    if statistic:
        arguments += ["--statistic", statistic]
    if pseudocount is not None:
        arguments += ["--pseudocount", pseudocount]
    if depths:
        arguments += ["--depths", *depths]
    arguments += ["--out", output]
    run_script("tissue_survey_summary.py", *arguments)
    return {row["sample"]: row for row in read_tsv(output)}


def raw_sum(rows, sample="liver_crossings"):
    return float(rows[sample]["raw_sum_positive_residuals"])


def test_difference_sums_only_positive_residuals(tmp_path):
    # I1: 10 - 4 = 6 counts; I2: 2 - 5 = -3 is dropped.
    rows = summarize(tmp_path, {"liver": {"I1": (10.0, 4), "I2": (2.0, 5)}}, statistic="difference")
    assert raw_sum(rows) == 6.0


def test_ratio_divides_by_crossings_plus_pseudocount(tmp_path):
    # I1: (10 - 4) / (4 + 1) = 1.2; I2: (2 - 5) / (5 + 1) is negative, dropped.
    rows = summarize(tmp_path, {"liver": {"I1": (10.0, 4), "I2": (2.0, 5)}}, statistic="ratio")
    assert raw_sum(rows) == pytest.approx(1.2)


def test_ratio_is_the_default(tmp_path):
    rows = summarize(tmp_path, {"liver": {"I1": (10.0, 4)}})
    assert rows["liver_crossings"]["statistic"] == "ratio"
    assert raw_sum(rows) == pytest.approx(1.2)


def test_zero_crossing_junction_is_kept_and_matches_the_difference(tmp_path):
    # The graceful-division guarantee: at pseudocount 1 a junction with no crossing
    # reads contributes exactly its coverage under both forms. These are the most
    # informative junctions, so neither form may drop them.
    sample = {"liver": {"I1": (10.0, 0)}}
    ratio_dir = tmp_path / "ratio"
    difference_dir = tmp_path / "difference"
    ratio_dir.mkdir()
    difference_dir.mkdir()
    as_ratio = summarize(ratio_dir, sample, statistic="ratio")
    as_difference = summarize(difference_dir, sample, statistic="difference")
    assert raw_sum(as_ratio) == 10.0
    assert raw_sum(as_difference) == 10.0


def test_missing_crossing_counts_are_treated_as_zero(tmp_path):
    # I2 has coverage but no crossing row, so its whole coverage counts as residual.
    rows = summarize(
        tmp_path, {"liver": {"I1": (10.0, 4), "I2": (3.0, 0)}}, statistic="difference"
    )
    assert raw_sum(rows) == 9.0


def test_samples_are_summarized_independently(tmp_path):
    rows = summarize(
        tmp_path,
        {"liver": {"I1": (10.0, 4)}, "kidney": {"I1": (20.0, 5)}},
        statistic="difference",
    )
    assert raw_sum(rows, "liver_crossings") == 6.0
    assert raw_sum(rows, "kidney_crossings") == 15.0


def test_depth_normalization_divides_the_sum(tmp_path):
    # Keyed on the file stem, which is what the script actually derives.
    rows = summarize(
        tmp_path,
        {"liver": {"I1": (10.0, 4)}},
        depths=["liver_crossings=2"],
        statistic="difference",
    )
    assert float(rows["liver_crossings"]["normalized_statistic"]) == 3.0


def run_expecting_failure(tmp_path, extra_arguments):
    crossings = write_tsv(tmp_path / "a_crossings.tsv", [["I1", 1]])
    coverage = write_coverage_table(tmp_path / "a_avecov.tsv", {"I1": 5.0})
    return subprocess.run(
        [
            sys.executable, str(SCRIPTS_DIR / "tissue_survey_summary.py"),
            "--crossings", str(crossings),
            "--avecov", str(coverage),
            "--out", str(tmp_path / "summary.tsv"),
            *extra_arguments,
        ],
        capture_output=True,
        text=True,
    )


def test_mismatched_file_counts_are_rejected(tmp_path):
    crossings = write_tsv(tmp_path / "a_crossings.tsv", [["I1", 1]])
    coverage_a = write_coverage_table(tmp_path / "a_avecov.tsv", {"I1": 5.0})
    coverage_b = write_coverage_table(tmp_path / "b_avecov.tsv", {"I1": 5.0})
    completed = subprocess.run(
        [
            sys.executable, str(SCRIPTS_DIR / "tissue_survey_summary.py"),
            "--crossings", str(crossings),
            "--avecov", str(coverage_a), str(coverage_b),
            "--out", str(tmp_path / "summary.tsv"),
        ],
        capture_output=True,
        text=True,
    )
    assert completed.returncode != 0


def test_non_positive_pseudocount_is_rejected(tmp_path):
    completed = run_expecting_failure(tmp_path, ["--pseudocount", "0"])
    assert completed.returncode != 0


@pytest.mark.xfail(
    strict=True,
    reason="The documented --depths key is the tissue name ('liver=32000000'), but "
    "sample_name_from_path (tissue_survey_summary.py:47) derives 'liver_crossings' "
    "from the filename. The lookup misses, normalization silently no-ops, and the "
    "warning only fires when --depths is omitted entirely -- so a no-select library "
    "sequenced 10x deeper compares against polyA unnormalized.",
)
def test_documented_depth_key_normalizes(tmp_path):
    rows = summarize(
        tmp_path, {"liver": {"I1": (10.0, 4)}}, depths=["liver=2"], statistic="difference"
    )
    assert float(rows["liver_crossings"]["normalized_statistic"]) == 3.0
