"""The trellis of the chosen pip (D10, D12 reversed 2026-10-31): choosing a pip cuts N at its orbit.

``T(q0) = (W^u computed, W^s[P, q0])``: every stable branch ends at its
pip-orbit point, and no unstable branch is cut, so a computed crossing deeper
than the pip's registered orbit stays in N. Choosing a pip cuts; un-choosing restores the
computed trellis. The synthetic tests hand-build the registry (explicit
canonical distances, keys and iterate links); one real build checks that the
strong-pip classification never sees the cut.
"""

from __future__ import annotations

import logging
from types import SimpleNamespace

import numpy as np
import pytest

from cases import build_k10
from helpers.fakes import bare_fixed_point, synthetic_trellis
from helpers.logs import assert_logged
from tanglepack.numerics.IntersectionRegistry import IntersectionRegistry
from tanglepack.topology.Trellis import Trellis

_LOGGER = "tanglepack.topology.Trellis"


def _window():
    """A period-1 (lambda = 4) window: the pip q0 at (s, u) = (4, 1), a middle
    crossing, and f(q0) = (1, 4), linked; plus one crossing past the stable cut
    and one past f(q0) on the unstable branch (f^m(q0), m = 1)."""
    fp = bare_fixed_point(1, beta=0.25)
    registry = IntersectionRegistry()
    stable, unstable = (fp, "stable", 0, 0), (fp, "unstable", 0, 0)

    def crossing(s: float, u: float) -> int:
        return registry.add_synthetic(
            (s, u), unstable_cdist=u, stable_cdist=s, manifold_a_key=unstable, manifold_b_key=stable
        )

    pip, middle, image = crossing(4.0, 1.0), crossing(3.0, 2.0), crossing(1.0, 4.0)
    registry.register_iterate(pip, 1, image)
    past_stable, past_unstable = crossing(8.0, 0.5), crossing(2.0, 6.0)
    return synthetic_trellis(registry, fp), pip, (middle, image), (past_stable, past_unstable)


def test_choosing_a_pip_cuts_the_stable_branches_only_and_unchoosing_restores():
    trellis, pip, inside, (past_stable, past_f_m) = _window()
    computed = {key: list(branch.intersection_ids) for key, branch in trellis.branches.items()}

    trellis.strong_pip = pip
    assert set(trellis.own_intersection_ids) == {pip, *inside, past_f_m}
    for branch in trellis.branches.values():
        assert past_stable not in branch.intersection_ids
        assert past_f_m in branch.intersection_ids
    assert pip in trellis and past_f_m in trellis and past_stable not in trellis

    trellis.strong_pip = None
    assert {key: branch.intersection_ids for key, branch in trellis.branches.items()} == computed
    assert past_stable in trellis


def test_iterates_stay_raw_lookups_past_the_cut():
    """X is every iterate of N, so the table still answers for points outside T."""
    trellis, pip, (_middle, image), _outside = _window()
    trellis.strong_pip = pip
    assert trellis.iterate(pip, 1) == image
    assert trellis.iterate(image, -1) == pip


def test_a_placed_cut_is_one_warning(caplog):
    """A period-3 pip with no registered +1 / +2 iterate: its cuts on the other two
    branches are placed by image_cdist, in ONE WARNING of the Trellis logger."""
    fp = bare_fixed_point(3, beta=0.5)
    registry = IntersectionRegistry()
    pip = registry.add_synthetic(
        (4.0, 1.0), unstable_cdist=1.0, stable_cdist=4.0,
        manifold_a_key=(fp, "unstable", 0, 0), manifold_b_key=(fp, "stable", 0, 0),
    )
    trellis = synthetic_trellis(registry, fp)
    with caplog.at_level(logging.WARNING, logger=_LOGGER):
        trellis.strong_pip = pip
        trellis.strong_pip = pip  # the same cut twice: no second warning
    assert_logged(caplog, logging.WARNING, _LOGGER, count=1)


def test_classification_never_sees_the_cut():
    """Every candidate is classified over the whole computed crossing set: a trellis
    cut at any candidate classifies the same candidates as a fresh snapshot."""
    case = build_k10(through="pips")
    session, fp = case.session, case.fixed_point
    fresh = Trellis.from_workbench(session.workbench, fp).classify_strong_pips(choose_default=False)
    for candidate in fresh:
        trellis = session.trellis(fp)
        trellis.set_strong_pip(candidate)
        assert trellis.classify_strong_pips(choose_default=False) == fresh
        assert trellis.strong_pip is None


def test_a_cached_bridge_polyline_is_read_only():
    """The memo shares one array between callers, so nobody can write to it."""
    trellis, *_ = _window()
    bridge = SimpleNamespace(version=0)
    polyline = trellis.cached_bridge_polyline(bridge, lambda: np.zeros((3, 2)))

    assert trellis.cached_bridge_polyline(bridge, lambda: None) is polyline
    assert not polyline.flags.writeable
    with pytest.raises(ValueError):
        polyline[0, 0] = 1.0
