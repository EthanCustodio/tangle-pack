"""Coverage guard for the 2026-10-05 test-suite refactor.

Compares two coverage.py JSON reports (``--cov-report=json``) of
``src/tanglepack`` and fails when a source line that the baseline executed is
missed in the candidate report, unless that line lies inside a whitelisted
dead-API qualname (whose tests the refactor drops while the code stays).

Usage::

    env/bin/python tests/_tools/coverage_guard.py BASE.json CANDIDATE.json

Prints the per-module coverage change and every newly missed line, and exits
with status 1 on any non-whitelisted new miss (0 otherwise).

This file lives under ``tests/_tools`` and is not collected by pytest (its name
does not match ``test_*.py``).

Dev Notes:
    The whitelist is resolved by qualname through ``ast`` on the CURRENT source
    tree, so it survives line shifts. A whitelisted qualname that no longer
    exists is reported (it is a stale entry), but does not fail the guard.
"""

from __future__ import annotations

import argparse
import ast
import json
import sys
from pathlib import Path
from typing import Iterable

# Dead API (plan decision 4): the code stays, its tests are dropped.
WHITELIST: dict[str, tuple[str, ...]] = {
    "src/tanglepack/numerics/IntersectionRegistry.py": (
        "IntersectionRegistry.on_interval",
        "IntersectionRegistry._get_lambda_u",
    ),
    "src/tanglepack/numerics/ManifoldMachine.py": ("ManifoldMachine._curvature_area",),
    "src/tanglepack/loom/TangleSession.py": ("TangleSession.invalidate_trellises",),
    "src/tanglepack/topology/Pseudoneighbor.py": ("forward_unstable_branch_cycle",),
    "src/tanglepack/topology/StrongPip.py": ("forward_stable_branch_cycle",),
}

REPO_ROOT = Path(__file__).resolve().parents[2]


def _qualname_spans(source: str) -> dict[str, tuple[int, int]]:
    """Map every function/class qualname in ``source`` to its (first, last) line span."""
    spans: dict[str, tuple[int, int]] = {}

    def visit(node: ast.AST, prefix: str) -> None:
        for child in ast.iter_child_nodes(node):
            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                name = f"{prefix}{child.name}"
                first = min(
                    [child.lineno] + [d.lineno for d in getattr(child, "decorator_list", [])]
                )
                spans[name] = (first, child.end_lineno or child.lineno)
                visit(child, f"{name}.")

    visit(ast.parse(source), "")
    return spans


def whitelisted_lines(path: str, qualnames: Iterable[str]) -> tuple[set[int], list[str]]:
    """Return the line numbers covered by ``qualnames`` in ``path`` and the stale names.

    Args:
        path: Repo-relative source path.
        qualnames: Dotted qualnames to whitelist.

    Returns:
        The set of whitelisted line numbers and the list of qualnames not found.
    """
    file = REPO_ROOT / path
    if not file.exists():
        return set(), list(qualnames)
    spans = _qualname_spans(file.read_text())
    lines: set[int] = set()
    stale: list[str] = []
    for name in qualnames:
        if name not in spans:
            stale.append(name)
            continue
        lo, hi = spans[name]
        lines.update(range(lo, hi + 1))
    return lines, stale


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


def compare(base_path: str, cand_path: str) -> int:
    """Compare two coverage reports and print the result.

    Args:
        base_path: Baseline coverage JSON.
        cand_path: Candidate coverage JSON.

    Returns:
        The process exit status: 1 when a non-whitelisted line is newly missed.
    """
    base = _load(base_path)
    cand = _load(cand_path)
    failures = 0
    for path in sorted(set(base) | set(cand)):
        if not path.startswith("src/tanglepack/"):
            continue
        b = base.get(path)
        c = cand.get(path)
        if b is None or c is None:
            print(f"{path}: present only in {'candidate' if b is None else 'baseline'}")
            continue
        b_pct = b["summary"]["percent_covered"]
        c_pct = c["summary"]["percent_covered"]
        base_executed = set(b["executed_lines"])
        new_missed = sorted(base_executed & set(c["missing_lines"]))
        allowed, stale = whitelisted_lines(path, WHITELIST.get(path, ()))
        bad = [line for line in new_missed if line not in allowed]
        ok = [line for line in new_missed if line in allowed]
        if abs(c_pct - b_pct) > 1e-9 or new_missed or stale:
            print(f"{path}: {b_pct:6.2f}% -> {c_pct:6.2f}% ({c_pct - b_pct:+.2f})")
        if ok:
            print(f"    newly missed, whitelisted dead API: {ok}")
        if bad:
            print(f"    NEWLY MISSED: {bad}")
            failures += len(bad)
        for name in stale:
            print(f"    stale whitelist entry (not found): {name}")
    if failures:
        print(f"coverage guard FAILED: {failures} newly missed line(s) outside the whitelist")
        return 1
    print("coverage guard OK: no newly missed line outside the whitelist")
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
    args = parser.parse_args(argv)
    return compare(args.base, args.candidate)


if __name__ == "__main__":
    sys.exit(main())
