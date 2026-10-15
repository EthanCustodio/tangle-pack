"""The cross-branch rules of the iterated (empty-stretch) cut (C19, WP10b).

A hole's image lobe sits on the stable branch of its BASE (the image chain's
``(x_0, x_m)``) and folds through stable branches along chords:

* on every stable branch, the base branch included, the interior crossings pair
  in STABLE order: walking a branch from its anchor, which no lobe image holds,
  the image is entered and left in turn (Jordan alternation; C19, B7b);
* the base lies on ONE branch (its ends are images of one branch's crossings), read
  at a non-anchor end, since an anchor lies on every branch of its point (D4); a
  base across two branches asserts;
* every empty stretch lies on one branch, so only an anchor can abut one along
  another branch (asserted).

The first three read ``IteratedHomotopyPartition._empty_stretches`` (an allowed
kernel) on synthetic chains; the last runs the public ``from_minimal`` on the
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
from tanglepack.topology.Trellis import Trellis

#: The module (the package attribute of the same name is the class).
family_module = importlib.import_module("tanglepack.topology.PartitionFamily")
LOGGER = "tanglepack.topology.PartitionFamily"

UNSTABLE = ("fp", "unstable", 0, 0)
OWN, FOREIGN = ("fp", "stable", 0, 1), ("fp", "stable", 1, 0)


def _crossing(branch: tuple, stable_cdist: float = 1.0) -> SimpleNamespace:
    """A crossing of ``UNSTABLE`` with ``branch``."""
    return SimpleNamespace(
        manifold_a_key=UNSTABLE, manifold_b_key=branch, stable_cdist=stable_cdist, is_anchor=False
    )


#: The anchor of ``UNSTABLE``'s point: on every branch there, ``s(0, 0)`` by C7.
ANCHOR = SimpleNamespace(
    manifold_a_key=UNSTABLE, manifold_b_key=("fp", "stable", 0, 0), stable_cdist=0.0, is_anchor=True
)


class _FakeTrellis:
    """Just enough trellis for ``_empty_stretches``: the crossings themselves."""

    #: The production anchor rules.
    unstable_key_of = Trellis.unstable_key_of
    stable_key_at = Trellis.stable_key_at

    def __init__(self, crossings: dict[int, SimpleNamespace]) -> None:
        self.registry = crossings

    def intersection(self, iid: int) -> SimpleNamespace:
        return self.registry[iid]


def _hole_homotopy(lo_id: int, hi_id: int, branch: tuple) -> SimpleNamespace:
    """A homotopy family with one open-open interval, the hole ``(lo, hi)``."""
    hole = PartitionInterval(
        lo_id=lo_id, hi_id=hi_id, lo_cdist=1.0, hi_cdist=2.0,
        closed_lo=False, closed_hi=False,
    )
    return SimpleNamespace(results={(branch, "right"): SimpleNamespace(intervals=[hole])})


def _empty_stretches(monkeypatch, crossings: dict, chain: list) -> dict:
    """``_empty_stretches`` of the hole ``(1, 2)`` on ``OWN`` whose image chain is ``chain``."""
    crossings = {1: _crossing(OWN), 2: _crossing(OWN), **crossings}
    minimal = SimpleNamespace(trellis=_FakeTrellis(crossings), hole_bridge_ids=[(1, 2)])
    monkeypatch.setattr(family_module, "_image_chain", lambda trellis, bid: chain)
    return IteratedHomotopyPartition._empty_stretches(minimal, _hole_homotopy(1, 2, OWN))


def _chain(crossings: list[int]) -> list[tuple[int, int]]:
    return list(zip(crossings, crossings[1:]))


def test_foreign_crossings_pair_in_stable_order(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    """The p3_6_blasts excursion: 41..44 on a foreign branch at stable cdists 1.158,
    1.025, 0.973, 1.285 pair as (43, 42) and (41, 44), never as the chain visits them."""
    crossings = {
        **{iid: _crossing(OWN) for iid in (10, 11, 12, 13)},
        41: _crossing(FOREIGN, 1.158),
        42: _crossing(FOREIGN, 1.025),
        43: _crossing(FOREIGN, 0.973),
        44: _crossing(FOREIGN, 1.285),
    }
    with caplog.at_level(logging.WARNING, logger=LOGGER):
        empty = _empty_stretches(monkeypatch, crossings, _chain([10, 11, 41, 42, 43, 44, 12, 13]))
    assert (13, "base") in empty[10]
    assert empty[43] == [(42, "chord")] and empty[41] == [(44, "chord")]
    assert_logged(caplog, logging.WARNING, LOGGER, count=0)


def test_base_crossings_pair_in_stable_order(monkeypatch: pytest.MonkeyPatch) -> None:
    """The k32_series hole (6, 7)'s lobe: its chain 8..13 spirals back across its
    own base branch, interior stable order 11, 10, 9, 12, so the chords are (11, 10)
    and (9, 12), the gaps inside the lobe image, not the chain's (9, 10), (11, 12)."""
    cdists = {8: 0.176, 9: 6.786, 10: 6.665, 11: 6.571, 12: 7.418, 13: 0.036}
    crossings = {iid: _crossing(OWN, cdist) for iid, cdist in cdists.items()}
    empty = _empty_stretches(monkeypatch, crossings, _chain(list(range(8, 14))))
    assert (13, "base") in empty[8]
    assert (10, "chord") in empty[11] and (12, "chord") in empty[9]
    assert (10, "chord") not in empty[9] and (12, "chord") not in empty[11]


def test_the_base_is_read_at_a_non_anchor_end(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    """A chain leaving the anchor whose C7 branch is the OTHER half (the
    inversion_period2 pattern) still has its base on the branch of its far end."""
    crossings = {0: ANCHOR, **{iid: _crossing(OWN) for iid in (11, 12, 13)}}
    with caplog.at_level(logging.WARNING, logger=LOGGER):
        empty = _empty_stretches(monkeypatch, crossings, _chain([0, 11, 12, 13]))
    assert (13, "base") in empty[0]
    assert (12, "chord") in empty[11]
    assert_logged(caplog, logging.WARNING, LOGGER, count=0)


def test_a_base_across_two_branches_asserts(monkeypatch: pytest.MonkeyPatch) -> None:
    crossings = {10: _crossing(OWN), 11: _crossing(OWN), 12: _crossing(FOREIGN)}
    with pytest.raises(AssertionError):
        _empty_stretches(monkeypatch, crossings, _chain([10, 11, 12]))


@pytest.fixture
def nested() -> Case:
    """The nested period-1 + period-3 tangle, outer zone blasted twice."""
    return build_nested()


def test_only_an_anchor_abuts_a_stretch_on_another_branch(
    nested: Case, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Every empty stretch lies on one branch: a far end on another stable branch
    at a non-anchor image crossing asserts (it used to be skipped at WARNING)."""
    session = nested.session
    minimal, homotopy = session.minimal_trellis(), session.homotopy_partition()
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

    original = family_module.IteratedHomotopyPartition._empty_stretches

    def with_foreign_far_ends(minimal_, homotopy_):
        empty = {iid: list(v) for iid, v in original(minimal_, homotopy_).items()}
        for bid in minimal_.image_bridge_ids:
            for iid in bid:
                empty.setdefault(iid, []).append((foreign_to(iid), "base"))
        return empty

    monkeypatch.setattr(
        family_module.IteratedHomotopyPartition,
        "_empty_stretches",
        staticmethod(with_foreign_far_ends),
    )
    with pytest.raises(AssertionError):
        IteratedHomotopyPartition.from_minimal(minimal, homotopy)
