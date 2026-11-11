"""Strong-pip classification is restricted to q0's own pair of manifolds.

A strong pip is a property of one pair, the unstable manifold of q0's unstable
fixed point and the stable branch cycle of its stable one, so only crossings of
that pair may disqualify a candidate: a homoclinic q0's own tangle, or for a
heteroclinic q0 in ``W^u(A) x W^s(B)`` the other ``W^u(A) x W^s(B)`` crossings
(2026-11-10). A heteroclinic q0 needs A and B of one period. These tests pin
that behaviour with a hand-built registry so the disqualifier is present and
unambiguous (the heavy Henon fixture happens to grow no such crossing).
"""

from __future__ import annotations

import inspect

import pytest

from helpers.fakes import bare_fixed_point
from tanglepack.numerics.FixedPoint import FixedPoint
from tanglepack.numerics.IntersectionRegistry import IntersectionRegistry
from tanglepack.topology.StrongPip import classify_strong_pips, is_strong_pip
from tanglepack.topology.Trellis import Trellis

#: The default slack of the strong-pip cdist-collision fallback.
COLLISION_RTOL = inspect.signature(is_strong_pip).parameters["collision_rtol"].default


def _fixed_point(period: int, lambda_u: float) -> FixedPoint:
    """A bare no-inversion fixed point whose full-cycle unstable eigenvalue is ``lambda_u``."""
    return bare_fixed_point(period, beta=lambda_u ** (-1.0 / period))


def _trellis(registry: IntersectionRegistry, *fixed_points: FixedPoint) -> Trellis:
    return Trellis(
        fixed_points=list(fixed_points),
        registry=registry,
        branches={},
        bridges=[],
    )


@pytest.mark.parametrize(
    "owner, disqualifies",
    [("other_fixed_point", False), ("same_fixed_point", True)],
    ids=["heteroclinic_ignored", "homoclinic_blocks"],
)
def test_only_the_own_tangle_disqualifies_a_strong_pip(owner: str, disqualifies: bool):
    """The one point inside q0's open box decides by whose tangle it belongs to.

    q0 is a homoclinic crossing of the period-3 inner tangle; the only point
    inside its open box ``(0, 1) x (0, 1)`` sits on the SAME stable branch. Its
    unstable side is:

    * ``other_fixed_point`` -- the period-1 saddle's: a heteroclinic crossing,
      not part of the period-3 tangle, so it is skipped and q0 stays strong;
    * ``same_fixed_point`` -- fp3's own: a homoclinic crossing, which blocks q0.
    """
    fp3 = _fixed_point(3, 4.0)
    fp1 = _fixed_point(1, 6.0)
    reg = IntersectionRegistry()

    stable_branch = (fp3, "stable", 0, 0)
    q0 = reg.add_synthetic(
        (0.0, 0.0), unstable_cdist=1.0, stable_cdist=1.0,
        manifold_a_key=(fp3, "unstable", 0, 0), manifold_b_key=stable_branch,
    )
    unstable_side = {
        "other_fixed_point": (fp1, "unstable", 0, 0),
        "same_fixed_point": (fp3, "unstable", 0, 0),
    }[owner]
    inside = reg.add_synthetic(
        (0.1, 0.1), unstable_cdist=0.5, stable_cdist=0.5,
        manifold_a_key=unstable_side, manifold_b_key=stable_branch,
    )

    result = is_strong_pip(_trellis(reg, fp3, fp1), q0)

    assert result.is_strong_pip is not disqualifies
    assert result.blocking_intersection_id == (inside if disqualifies else None)


def test_table_linked_own_iterate_does_not_disqualify_strong_pip():
    """q0's own forward image on the next stable branch is recognised by
    iterate-table lookup, not by the cdist collision: registered
    2 * collision_rtol short on the stable side and slightly short on the unstable
    side, it maps back strictly inside q0's box by scaling and would block q0
    — unless the table says it IS q0's image. Unlinked, the collision fallback
    still blocks (unchanged behaviour)."""
    fp3 = _fixed_point(3, 4.0)
    beta = fp3.per_step_beta("unstable")

    def build(link: bool):
        reg = IntersectionRegistry()
        q0 = reg.add_synthetic(
            (0.0, 0.0), unstable_cdist=1.0, stable_cdist=1.0,
            manifold_a_key=(fp3, "unstable", 0, 0),
            manifold_b_key=(fp3, "stable", 0, 0),
        )
        image = reg.add_synthetic(
            (0.1, 0.1),
            unstable_cdist=0.995 * beta,
            stable_cdist=(1.0 - 2.0 * COLLISION_RTOL) / beta,
            manifold_a_key=(fp3, "unstable", 1, 0),
            manifold_b_key=(fp3, "stable", 1, 0),
        )
        if link:
            reg.register_iterate(q0, 1, image)
        return reg, q0, image

    reg, q0, _image = build(link=True)
    result = is_strong_pip(_trellis(reg, fp3), q0)
    assert result.is_strong_pip
    assert result.blocking_intersection_id is None

    reg, q0, image = build(link=False)
    result = is_strong_pip(_trellis(reg, fp3), q0)
    assert not result.is_strong_pip
    assert result.blocking_intersection_id == image


@pytest.mark.parametrize(
    "q0_pair, inside_pair, disqualifies",
    [
        (("A", "B"), ("A", "B"), True),
        (("A", "B"), ("B", "B"), False),
        (("B", "B"), ("A", "B"), False),
        (("B", "B"), ("B", "B"), True),
    ],
    ids=["heteroclinic_own_pair_blocks", "heteroclinic_ignores_homoclinic",
         "homoclinic_ignores_heteroclinic", "homoclinic_own_pair_blocks"],
)
def test_only_q0s_own_pair_disqualifies(q0_pair, inside_pair, disqualifies):
    """Two period-1 saddles A and B; q0 and the one point inside its open box
    ``(0, 1) x (0, 1)`` both lie on B's stable branch, each named by its
    ``(unstable, stable)`` fixed points. The point blocks q0 iff it joins q0's
    own pair: a heteroclinic q0 in ``W^u(A) x W^s(B)`` is blocked by another
    ``W^u(A) x W^s(B)`` crossing and not by a homoclinic one of B, and a
    homoclinic q0 is unaffected by heteroclinic crossings (today's rule)."""
    fps = {"A": _fixed_point(1, 4.0), "B": _fixed_point(1, 6.0)}
    reg = IntersectionRegistry()

    def crossing(coords, cdist, pair):
        unstable, stable = (fps[name] for name in pair)
        return reg.add_synthetic(
            coords, unstable_cdist=cdist, stable_cdist=cdist,
            manifold_a_key=(unstable, "unstable", 0, 0), manifold_b_key=(stable, "stable", 0, 0),
        )

    q0 = crossing((0.0, 0.0), 1.0, q0_pair)
    inside = crossing((0.1, 0.1), 0.5, inside_pair)

    result = is_strong_pip(_trellis(reg, fps["A"], fps["B"]), q0)

    assert result.is_strong_pip is not disqualifies
    assert result.blocking_intersection_id == (inside if disqualifies else None)


def test_heteroclinic_pip_across_periods_raises():
    """A heteroclinic q0 between a period-1 and a period-3 fixed point is no
    strong pip: the direct test raises, and the batch classification skips it."""
    fp1, fp3 = _fixed_point(1, 6.0), _fixed_point(3, 4.0)
    reg = IntersectionRegistry()
    q0 = reg.add_synthetic(
        (0.0, 0.0), unstable_cdist=1.0, stable_cdist=1.0,
        manifold_a_key=(fp1, "unstable", 0, 0), manifold_b_key=(fp3, "stable", 0, 0),
    )
    trellis = _trellis(reg, fp1, fp3)

    with pytest.raises(ValueError, match="periods"):
        is_strong_pip(trellis, q0)
    assert classify_strong_pips(trellis, [q0]) == {}
