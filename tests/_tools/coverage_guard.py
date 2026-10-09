"""Coverage guard: no source line the baseline executed may go missed.

Compares two coverage.py JSON reports (``--cov-report=json``) of
``src/tanglepack`` and fails when a statement that the baseline executed is
missed in the candidate report.

Usage::

    env/bin/python tests/_tools/coverage_guard.py BASE.json CANDIDATE.json [--base-ref REF]

The baseline source is read from git (``git show REF:path``, ``REF`` defaults
to ``HEAD``), the candidate source from the working tree. Prints the
per-module coverage change and every newly missed line (at its candidate line
number), and exits with status 1 on any new miss (0 otherwise).

Note:
    A line is matched across the two versions by its KEY, never by its number
    (TEST-6, WP11): ``(enclosing qualname, stripped source text, occurrence of
    that text inside the qualname)``, with ``<module>`` for module-level lines.
    An edit above a covered statement shifts its number but not its key, so it
    is neither a false miss nor a mask for a real one. A deleted statement has
    no candidate key and is never reported; dead API is deleted, not
    whitelisted.
"""

from __future__ import annotations

import argparse
import ast
import json
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Callable, Iterable

REPO_ROOT = Path(__file__).resolve().parents[2]

LineKey = tuple[str, str, int]


def _qualname_spans(source: str) -> list[tuple[str, int, int]]:
    """Every function/class qualname of ``source`` with its (first, last) line, outermost first."""
    spans: list[tuple[str, int, int]] = []

    def visit(node: ast.AST, prefix: str) -> None:
        for child in ast.iter_child_nodes(node):
            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                name = f"{prefix}{child.name}"
                first = min([child.lineno] + [d.lineno for d in child.decorator_list])
                spans.append((name, first, child.end_lineno or child.lineno))
                visit(child, f"{name}.")

    visit(ast.parse(source), "")
    return spans


def line_keys(source: str, lines: Iterable[int]) -> dict[int, LineKey]:
    """Key each of ``lines`` by (enclosing qualname, stripped text, occurrence index).

    Args:
        source: The module's source text.
        lines: 1-based line numbers to key.

    Returns:
        ``{line: key}``; the occurrence index counts the same stripped text
        earlier in the same qualname, so repeated statements stay distinct.
    """
    text = source.splitlines()
    owner = ["<module>"] * (len(text) + 1)
    for name, first, last in _qualname_spans(source):  # outermost first: inner spans overwrite
        owner[first:last + 1] = [name] * (last + 1 - first)
    seen: Counter[tuple[str, str]] = Counter()
    keys: dict[int, LineKey] = {}
    wanted = set(lines)
    for number, line in enumerate(text, start=1):
        pair = (owner[number], line.strip())
        if number in wanted:
            keys[number] = (*pair, seen[pair])
        seen[pair] += 1
    return keys


def _load(path: str) -> dict[str, dict]:
    """Load the ``files`` table of a coverage.py JSON report, keyed by repo-relative path."""
    with open(path) as handle:
        data = json.load(handle)
    files: dict[str, dict] = {}
    for name, payload in data["files"].items():
        key = name
        if Path(name).is_absolute():
            try:
                key = str(Path(name).relative_to(REPO_ROOT))
            except ValueError:
                pass
        files[key] = payload
    return files


def git_source(ref: str) -> Callable[[str], str]:
    """A reader of repo-relative paths at the git revision ``ref``."""
    def read(path: str) -> str:
        return subprocess.run(["git", "show", f"{ref}:{path}"], cwd=REPO_ROOT,
                              capture_output=True, text=True, check=True).stdout
    return read


def _working_tree(path: str) -> str:
    return (REPO_ROOT / path).read_text()


def compare(base: dict[str, dict], cand: dict[str, dict], base_source: Callable[[str], str],
            cand_source: Callable[[str], str] = _working_tree) -> int:
    """Compare two coverage ``files`` tables and print the result.

    Args:
        base: Baseline report's ``files`` table, keyed by repo-relative path.
        cand: Candidate report's ``files`` table, keyed the same way.
        base_source: Reads a path's baseline source.
        cand_source: Reads a path's candidate source (default: the working tree).

    Returns:
        The process exit status: 1 when a baseline-executed statement is missed.
    """
    failures = 0
    for path in sorted(set(base) | set(cand)):
        if not path.startswith("src/tanglepack/"):
            continue
        b, c = base.get(path), cand.get(path)
        if b is None or c is None:
            print(f"{path}: present only in {'candidate' if b is None else 'baseline'}")
            continue
        executed = set(line_keys(base_source(path), b["executed_lines"]).values())
        missing = line_keys(cand_source(path), c["missing_lines"])
        new_missed = sorted(line for line, key in missing.items() if key in executed)
        b_pct, c_pct = b["summary"]["percent_covered"], c["summary"]["percent_covered"]
        if abs(c_pct - b_pct) > 1e-9 or new_missed:
            print(f"{path}: {b_pct:6.2f}% -> {c_pct:6.2f}% ({c_pct - b_pct:+.2f})")
        if new_missed:
            print(f"    NEWLY MISSED: {new_missed}")
            failures += len(new_missed)
    if failures:
        print(f"coverage guard FAILED: {failures} newly missed line(s)")
        return 1
    print("coverage guard OK: no newly missed line")
    return 0


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point.

    Args:
        argv: Arguments (defaults to ``sys.argv[1:]``).

    Returns:
        The exit status.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("base")
    parser.add_argument("candidate")
    parser.add_argument("--base-ref", default="HEAD", help="git revision of the baseline source")
    args = parser.parse_args(argv)
    return compare(_load(args.base), _load(args.candidate), git_source(args.base_ref))


if __name__ == "__main__":
    sys.exit(main())
