"""Registry insertion stays near-linear (opt-in: ``pytest -m perf``).

The sorted key arrays plus the collision prefilter replace the old O(N) rescan
per insert. This is a wall-clock measurement, so it is excluded from the
default run (``addopts = -m 'not perf'``) and selected with ``-m perf``.
"""

from __future__ import annotations

import gc
import time

import pytest

from tanglepack.numerics.IntersectionRegistry import IntersectionRegistry


@pytest.mark.perf
def test_registry_insert_is_near_linear() -> None:
    """20k inserts complete quickly and the per-insert cost does not blow up."""
    registry = IntersectionRegistry()

    def add_block(start: int, count: int) -> float:
        t0 = time.perf_counter()
        for i in range(start, start + count):
            registry.add_synthetic(
                coords=(float(i), 0.0),
                unstable_cdist=float(i),
                stable_cdist=float(20000 - i),
            )
        return time.perf_counter() - t0

    # Fix for a PRE-EXISTING flake (it failed roughly one full-suite run in three
    # at Phase 5's HEAD, while always passing standalone; more tests in the suite
    # made it worse, which is what surfaced it). The ratio below compares two
    # wall-clock blocks, so it is only meaningful against a settled heap: run
    # inside the full suite, with several hundred other tests' garbage still
    # uncollected, the FIRST block is measured on a fragmented allocator and comes
    # out artificially fast, and the ratio blows past any bound. Collecting first
    # makes the two blocks comparable. The bound itself is unchanged.
    gc.collect()

    total_start = time.perf_counter()
    first_block = add_block(0, 1000)
    add_block(1000, 18000)
    last_block = add_block(19000, 1000)
    total = time.perf_counter() - total_start

    assert len(registry) == 20000
    assert total < 2.0, f"20k inserts took {total:.2f}s"
    # The bisect is O(log N) but the list insert behind it is O(N), so the
    # per-insert cost does grow with N by construction: this bound only has to
    # catch a return to the O(N) RESCAN (which was ~50x across this range), and
    # it has to survive a loaded machine, so it is deliberately loose.
    assert last_block < 15 * first_block, (
        f"per-insert cost grew from {first_block:.4f}s/1k to {last_block:.4f}s/1k"
    )
