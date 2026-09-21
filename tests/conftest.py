"""Shared helpers for the fixture-based workflow tests."""
import csv
import subprocess
import sys
from pathlib import Path

SKILL_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = SKILL_ROOT / "scripts"


def run_script(script_name, *arguments):
    """Run a workflow script; fail with its stderr if it exits non-zero."""
    completed = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / script_name), *map(str, arguments)],
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, f"{script_name} exited {completed.returncode}:\n{completed.stderr}"
    return completed


def write_tsv(path, rows, header=None):
    lines = []
    if header:
        lines.append("\t".join(header))
    lines.extend("\t".join(str(field) for field in row) for row in rows)
    path.write_text("\n".join(lines) + "\n")
    return path


def write_coverage_table(path, coverage_by_intron):
    """A minimal avecov table: classify/summary read intron_id and avg_coverage."""
    return write_tsv(
        path,
        sorted(coverage_by_intron.items()),
        header=["intron_id", "avg_coverage"],
    )


def read_tsv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))
