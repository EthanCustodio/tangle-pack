"""
Regression (PK-3, D5): no branch of a period-k orbit is left silently without
its window holes.

A reference trajectory appears on each of the k stable branches of its orbit
at ``+1 .. +(k-1)``. A period-k chain is grown one map step at a time, so its
unstable branches are staggered (each reaches beta further than the one before
it in the cycle) and a reference on the outermost lap of the longest branch has
no registered ``+j`` image: the exact iterate table (WP2) has no link, and the
image crossing does not exist in the trellis (or lies past its pip-orbit cut,
D12: a pair never leaves T(q0)). Nothing is invented for it, but
each stable branch left without that appearance is named in one WARNING of the
pseudoneighbor logger. Which branches those are is read off the table by rule;
no count of a build is pinned.
"""

from __future__ import annotations

import logging

import pytest

from cases import BUILDERS
from helpers.logs import assert_logged

_LOGGER = "tanglepack.topology.Pseudoneighbor"


def _gaps(trellis) -> set[tuple[tuple[int, int], tuple]]:
    """``(reference, stable branch)`` for every appearance the trellis does not hold.

    A reference whose first member image off the trellis is at ``+j`` (``j <
    k_value``: unregistered, or past the pip-orbit cut, D12) has no appearance
    on ``advance_key(branch, i)`` for ``i = j .. k_value - 1``.
    """
    gaps = set()
    for ref in trellis.reference_pseudoneighbors:
        fp = ref.branch_key[0]
        j = next(
            (
                j
                for j in range(1, fp.k_value)
                if any(trellis.iterate(end, j) not in trellis for end in ref.as_tuple())
            ),
            fp.k_value,
        )
        gaps |= {(ref.as_tuple(), fp.advance_key(ref.branch_key, i)) for i in range(j, fp.k_value)}
    return gaps


def _rerun_pseudoneighbors(caplog, trellis) -> None:
    """Recompute the pairs of a built trellis with the pseudoneighbor logger captured."""
    with caplog.at_level(logging.WARNING, logger=_LOGGER):
        trellis.compute_pseudoneighbors()


@pytest.mark.parametrize("name", ["inversion_period2", "p3_6_blasts"])
def test_missing_window_appearances_are_named_per_branch(caplog, name):
    case = BUILDERS[name]()
    (fp,) = case.fixed_points
    trellis = case.session.trellis(fp)
    gaps = _gaps(trellis)
    assert gaps, "the builder no longer leaves a +j appearance unregistered"

    # Nothing invented: no direct hole of a gapped reference sits on its gap branch.
    direct = {
        (hole.origin, hole.pair.branch_key)
        for hole in trellis.holes
        if hole.side_source == "signs"
    }
    assert not direct & gaps

    caplog.clear()
    _rerun_pseudoneighbors(caplog, trellis)
    assert_logged(caplog, logging.WARNING, _LOGGER, count=len({branch for _, branch in gaps}))


@pytest.mark.parametrize("name", ["p3", "k28_one_blast", "inversion"])
def test_complete_trajectories_log_nothing(caplog, name):
    case = BUILDERS[name]()
    trellis = case.session.trellis(case.fixed_points[0])
    assert not _gaps(trellis), "the builder now leaves a +j appearance unregistered"

    caplog.clear()
    _rerun_pseudoneighbors(caplog, trellis)
    assert_logged(caplog, logging.WARNING, _LOGGER, count=0)
