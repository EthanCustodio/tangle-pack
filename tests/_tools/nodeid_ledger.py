"""Account for every test node id the 2026-10-05 refactor removed.

Usage::

    env/bin/python tests/_tools/nodeid_ledger.py RUN_DIR LEDGER.md

``RUN_DIR`` holds the per-phase collections ``nodeids_<phase>.txt`` (the
baseline is ``nodeids_base.txt``), the per-phase deletion lists
``<phase>_deleted.txt`` and, for a move-only phase, ``<phase>_moves.tsv``
(``old_path<TAB>new_path`` per line, as ``git diff --name-status -M`` prints a
rename). The phases are walked in :data:`PHASES` order; a missing snapshot is
skipped (its transition folds into the next one).

For each transition ``prev -> next`` every node id in ``prev`` but not in
``next`` must be either listed in ``next``'s deletion list or carried by a
recorded move (its path renamed, the rest of the id unchanged and present in
``next``). Every listed deletion must also name its test function somewhere in
the ledger markdown, so no deletion goes unexplained. Exit status 1 on any
unaccounted id.

Dev Notes:
    This is a one-off audit tool for the refactor's final verification (plan
    "Verification (end to end)"); it is not collected by pytest (no ``test_``
    prefix) and has no runtime dependency on the library.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

PHASES: tuple[str, ...] = (
    "base", "p2", "p3a", "p3b", "p4", "p5", "p6", "p7a", "p7b", "p8", "p9", "p10a", "p10",
)


def _read_ids(path: Path) -> set[str]:
    """The non-empty node ids of one listing, one per line."""
    return {line.strip() for line in path.read_text().splitlines() if "::" in line}


def _read_moves(path: Path) -> dict[str, str]:
    """``old path -> new path`` of a move-only phase, or empty when none."""
    if not path.exists():
        return {}
    moves = {}
    for line in path.read_text().splitlines():
        if line.strip():
            old, new = line.split("\t")
            moves[old] = new
    return moves


def _function_name(node_id: str) -> str:
    """The test function of a node id, without its parameter brackets."""
    return re.sub(r"\[.*\]$", "", node_id.split("::")[-1])


def audit(run_dir: Path, ledger: Path) -> int:
    """Walk the phases and report every unaccounted removal.

    Args:
        run_dir: The refactor's run folder.
        ledger: The ledger markdown that must explain every deletion.

    Returns:
        The number of problems found (0 = every removed id is accounted for).
    """
    ledger_text = ledger.read_text()
    snapshots = [(p, run_dir / f"nodeids_{p}.txt") for p in PHASES]
    snapshots = [(p, path) for p, path in snapshots if path.exists()]
    problems = 0
    total_removed = 0
    for (prev_name, prev_path), (name, path) in zip(snapshots, snapshots[1:]):
        prev, current = _read_ids(prev_path), _read_ids(path)
        removed = prev - current
        deleted_path = run_dir / f"{name}_deleted.txt"
        deleted = _read_ids(deleted_path) if deleted_path.exists() else set()
        moves = _read_moves(run_dir / f"{name}_moves.tsv")
        moved = set()
        for node_id in removed:
            file_path, _, rest = node_id.partition("::")
            new_path = moves.get(file_path)
            if new_path is not None and f"{new_path}::{rest}" in current:
                moved.add(node_id)
        unaccounted = removed - deleted - moved
        unexplained = sorted(
            node_id for node_id in deleted if _function_name(node_id) not in ledger_text
        )
        stale = sorted(deleted - removed)
        total_removed += len(removed)
        print(
            f"{prev_name:>5} -> {name:<5} removed {len(removed):4d} "
            f"(deleted {len(removed & deleted):4d}, moved {len(moved):4d}), "
            f"added {len(current - prev):4d}, collected {len(current):5d}"
        )
        for node_id in sorted(unaccounted):
            print(f"    UNACCOUNTED {node_id}")
        for node_id in unexplained:
            print(f"    NOT IN LEDGER {node_id}")
        for node_id in stale:
            print(f"    LISTED BUT NOT REMOVED {node_id}")
        problems += len(unaccounted) + len(unexplained) + len(stale)
    base = _read_ids(snapshots[0][1])
    final = _read_ids(snapshots[-1][1])
    print(
        f"baseline {len(base)} -> final {len(final)} ({snapshots[-1][0]}); "
        f"{len(base - final)} baseline ids gone, {total_removed} removals walked, "
        f"{problems} problem(s)"
    )
    return problems


def main(argv: list[str]) -> int:
    """Command-line entry point.

    Args:
        argv: ``[RUN_DIR, LEDGER.md]``.

    Returns:
        The process exit status.
    """
    if len(argv) != 2:
        print(__doc__)
        return 2
    return 1 if audit(Path(argv[0]), Path(argv[1])) else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
