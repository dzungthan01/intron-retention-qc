"""Stage 1: SKILL.md's commands must match what the scripts actually accept.

Reads argparse declarations with `ast` rather than running the scripts, so a
broken command line is caught without needing fixtures or real data.
"""
import ast
import re
from pathlib import Path

import pytest

SKILL_ROOT = Path(__file__).resolve().parent.parent
SKILL_MD = SKILL_ROOT / "SKILL.md"
SCRIPTS_DIR = SKILL_ROOT / "scripts"


def fenced_blocks(text):
    return re.findall(r"```(?:\w+)?\n(.*?)```", text, re.DOTALL)


def documented_invocations():
    """[(script_name, [flags]), ...] for every `python scripts/X.py` in SKILL.md."""
    invocations = []
    for block in fenced_blocks(SKILL_MD.read_text()):
        joined = re.sub(r"\\\n\s*", " ", block)
        for line in joined.splitlines():
            match = re.search(r"python\s+scripts/(\S+\.py)", line)
            if not match:
                continue
            flags = re.findall(r"(--[A-Za-z0-9][-A-Za-z0-9]*)", line)
            invocations.append((match.group(1), flags))
    return invocations


def declared_arguments(script_name):
    """flag -> required, from the script's ap.add_argument(...) calls."""
    tree = ast.parse((SCRIPTS_DIR / script_name).read_text())
    declared = {}
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        if not isinstance(node.func, ast.Attribute) or node.func.attr != "add_argument":
            continue
        flags = [
            arg.value
            for arg in node.args
            if isinstance(arg, ast.Constant) and str(arg.value).startswith("--")
        ]
        if not flags:
            continue
        required = any(
            keyword.arg == "required"
            and isinstance(keyword.value, ast.Constant)
            and keyword.value.value is True
            for keyword in node.keywords
        )
        for flag in flags:
            declared[flag] = required
    return declared


INVOCATIONS = documented_invocations()


def invocation_id(invocation):
    return invocation[0]


def test_skill_md_documents_every_script():
    documented = {name for name, _ in INVOCATIONS}
    on_disk = {path.name for path in SCRIPTS_DIR.glob("*.py")}
    assert documented == on_disk


@pytest.mark.parametrize("invocation", INVOCATIONS, ids=invocation_id)
def test_documented_script_exists(invocation):
    script_name, _ = invocation
    assert (SCRIPTS_DIR / script_name).is_file()


@pytest.mark.parametrize("invocation", INVOCATIONS, ids=invocation_id)
def test_documented_flags_are_accepted(invocation):
    script_name, flags = invocation
    declared = declared_arguments(script_name)
    unknown = [flag for flag in flags if flag not in declared]
    assert unknown == []


@pytest.mark.parametrize("invocation", INVOCATIONS, ids=invocation_id)
def test_required_flags_are_documented(invocation):
    script_name, flags = invocation
    declared = declared_arguments(script_name)
    missing = [flag for flag, required in declared.items() if required and flag not in flags]
    assert missing == []
