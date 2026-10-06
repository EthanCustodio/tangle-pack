"""The 2026-09-30 cross-branch rules of the iterated (empty-stretch) cut.

A hole's image lobe sits on the stable branch of its BASE (the image chain's
``(x_0, x_m)``) and folds through it along chords of consecutive interior
crossings. Three rules keep other stable branches out of that reading:

* only the base branch's crossings are paired into chords (an excursion of the
  image across another orbit branch or fixed point folds nothing there);
* a chain ending on another branch than it starts marks nothing (WARNING);
* an empty stretch whose far end lies on another stable branch is no flank of
  an image crossing: it is skipped at WARNING and cuts nothing.

The first two read ``IteratedHomotopyPartition._empty_stretches`` (an allowed
kernel) on synthetic chains; the third runs the public ``from_minimal`` on the
nested tangle with one foreign far end injected per image crossing.
"""

from __future__ import annotations

import importlib
import logging
from types import SimpleNamespace

import pytest

from cases import Case, build_nested
from helpers.logs import assert_logged
from tanglepack.topology.PartitionFamily import IteratedHomotopyPartition
from tanglepack.topology.TopologyResults import PartitionInterval

#: The module (the package attribute of the same name is the class).
family_module = importlib.import_module("tanglepack.topology.PartitionFamily")
LOGGER = "tanglepack.topology.PartitionFamily"


class _FakeTrellis:
    """Just enough trellis for ``_empty_stretches``: a crossing's stable branch."""

    def __init__(self, branches: dict[int, tuple]) -> None:
        self._branches = branches

    def intersection(self, iid: int) -> SimpleNamespace:
        """The crossing, carrying only its stable branch key."""
        return SimpleNamespace(manifold_b_key=self._branches[iid])


def _hole_homotopy(lo_id: int, hi_id: int, branch: tuple) -> SimpleNamespace:
    """A homotopy family with one open-open interval, the hole ``(lo, hi)``."""
    hole = PartitionInterval(
        lo_id=lo_id, hi_id=hi_id, lo_cdist=1.0, hi_cdist=2.0,
        closed_lo=False, closed_hi=False,
    )
    return SimpleNamespace(results={(branch, "right"): SimpleNamespace(intervals=[hole])})


def test_chords_pair_only_the_base_branchs_crossings(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    own, foreign = ("fp", "stable", 1, 0), ("fp", "stable", 2, 0)
    branches = {1: own, 2: own, 10: own, 11: own, 12: own, 13: own, 98: foreign, 99: foreign}
    minimal = SimpleNamespace(trellis=_FakeTrellis(branches), hole_bridge_ids=[(1, 2)])
    # The image of the hole bridge crosses the foreign branch twice between
    # its two folds on its own branch.
    chain = [(10, 11), (11, 99), (99, 98), (98, 12), (12, 13)]
    monkeypatch.setattr(family_module, "_image_chain", lambda trellis, bid: chain)
    with caplog.at_level(logging.INFO, logger=LOGGER):
        empty = IteratedHomotopyPartition._empty_stretches(
            minimal, _hole_homotopy(1, 2, own)
        )
    assert (13, "base") in empty[10]
    assert (12, "chord") in empty[11]
    assert 99 not in empty and 98 not in empty
    assert_logged(caplog, logging.INFO, LOGGER)


def test_a_lobe_ending_on_another_branch_marks_nothing(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    own, foreign = ("fp", "stable", 1, 0), ("fp", "stable", 2, 0)
    branches = {1: own, 2: own, 10: own, 11: own, 12: foreign}
    minimal = SimpleNamespace(trellis=_FakeTrellis(branches), hole_bridge_ids=[(1, 2)])
    monkeypatch.setattr(family_module, "_image_chain", lambda trellis, bid: [(10, 11), (11, 12)])
    with caplog.at_level(logging.WARNING, logger=LOGGER):
        empty = IteratedHomotopyPartition._empty_stretches(
            minimal, _hole_homotopy(1, 2, own)
        )
    assert set(empty) == {1, 2}  # only the hole itself
    assert_logged(caplog, logging.WARNING, LOGGER)


@pytest.fixture
def nested() -> Case:
    """The nested period-1 + period-3 tangle, outer zone blasted twice."""
    return build_nested()


def test_a_far_end_on_another_branch_is_no_flank(
    nested: Case, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    """An empty stretch reaching another stable branch is skipped at WARNING and
    cuts nothing: the partition is the one without it."""
    session = nested.session
    minimal, homotopy = session.minimal_trellis(), session.homotopy_partition()
    unpatched = session.iterated_partition()
    trellis = minimal.trellis
    assert minimal.image_bridge_ids, "the nested case carries image bridges"
    stable_ids = {
        key: list(branch.intersection_ids)
        for key, branch in trellis.branches.items()
        if key[1] == "stable" and branch.intersection_ids
    }
    assert len(stable_ids) > 1

    def foreign_to(iid: int) -> int:
        own = trellis.intersection(iid).manifold_b_key
        return next(ids[-1] for key, ids in stable_ids.items() if key != own)

    with caplog.at_level(logging.WARNING, logger=LOGGER):
        IteratedHomotopyPartition.from_minimal(minimal, homotopy)
    baseline = len(
        [r for r in caplog.records if r.name == LOGGER and r.levelno == logging.WARNING]
    )
    caplog.clear()
    original = family_module.IteratedHomotopyPartition._empty_stretches

    def with_foreign_far_ends(minimal_, homotopy_):
        empty = {iid: list(v) for iid, v in original(minimal_, homotopy_).items()}
        for bid in minimal_.image_bridge_ids:
            for iid in bid:
                other = foreign_to(iid)
                empty.setdefault(iid, []).append((other, "base"))
        return empty

    monkeypatch.setattr(
        family_module.IteratedHomotopyPartition,
        "_empty_stretches",
        staticmethod(with_foreign_far_ends),
    )
    with caplog.at_level(logging.WARNING, logger=LOGGER):
        iterated = IteratedHomotopyPartition.from_minimal(minimal, homotopy)
    warnings = assert_logged(caplog, logging.WARNING, LOGGER)
    assert len(warnings) > baseline, "each foreign far end is skipped at WARNING"
    assert iterated.signature() == unpatched.signature()
    applied = [cut for cut in iterated.cuts if cut.opened is not None]
    assert applied
    assert all(
        trellis.intersection(cut.partner).manifold_b_key
        == trellis.intersection(cut.intersection_id).manifold_b_key
        for cut in applied
    )
