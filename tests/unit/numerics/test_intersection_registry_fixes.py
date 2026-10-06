"""Regression tests for the Phase 1 Intersection / IntersectionRegistry fixes.

Pins:
    1.9  ``Intersection.fixed_points`` must not raise IndexError when only
         ``manifold_b_key`` is set (the normal shape for crossings born on an
         iterated bridge).
    1.10 ``Intersection.synthetic`` must build with keyword arguments so the
         label does not land in the ``id`` slot, and must forward manifold keys.
"""

from __future__ import annotations

from typing import Optional

import pytest

from helpers.fakes import bare_fixed_point
from tanglepack.numerics.Intersection import Intersection
from tanglepack.numerics.IntersectionRegistry import IntersectionRegistry


# ── 1.9 Intersection.fixed_points ──────────────────────────────────────────

#: (a-key fixed point, b-key fixed point, expected ``fixed_points``), with
#: ``"p"`` / ``"q"`` naming two distinct fixed points and ``None`` an unset key.
FIXED_POINTS_TABLE: list[tuple[Optional[str], Optional[str], tuple[str, ...]]] = [
    (None, "p", ("p",)),      # only the b key: the IndexError case (1.9)
    ("p", None, ("p",)),      # only the a key
    ("p", "p", ("p",)),       # one fixed point on both keys: deduplicated
    ("p", "q", ("p", "q")),   # two fixed points: a-key first
    (None, None, ()),         # no keys at all
]


@pytest.mark.parametrize(
    ("a_fp", "b_fp", "expected"),
    FIXED_POINTS_TABLE,
    ids=["only_b", "only_a", "same_deduped", "distinct_ab_order", "no_keys"],
)
def test_fixed_points_table(
    a_fp: Optional[str], b_fp: Optional[str], expected: tuple[str, ...]
) -> None:
    """``Intersection.fixed_points`` reads whichever keys are set, a-key first."""
    points = {"p": bare_fixed_point(beta=1 / 3), "q": bare_fixed_point(beta=1 / 5)}
    keys = {}
    if a_fp is not None:
        keys["manifold_a_key"] = (points[a_fp], "unstable", 0, 0)
    if b_fp is not None:
        keys["manifold_b_key"] = (points[b_fp], "stable", 0, 0)

    ix = Intersection(coords=(0.0, 0.0), unstable_cdist=1.0, stable_cdist=1.0, **keys)

    assert ix.fixed_points == tuple(points[name] for name in expected)


# ── 1.10 Intersection.synthetic ────────────────────────────────────────────


def test_synthetic_keeps_the_label_out_of_the_id_slot_and_forwards_keys() -> None:
    """The label lands in ``label`` (not ``id``) and both manifold keys are kept."""
    fp = bare_fixed_point(beta=1 / 3)
    a_key = (fp, "unstable", 0, 0)
    b_key = (fp, "stable", 0, 1)

    ix = Intersection.synthetic(
        (1.0, 2.0),
        0.5,
        0.25,
        label="anchor",
        manifold_a_key=a_key,
        manifold_b_key=b_key,
    )

    assert ix.id is None
    assert ix.label == "anchor"
    assert ix.is_synthetic
    assert ix.coords == (1.0, 2.0)
    assert ix.unstable_cdist == 0.5
    assert ix.stable_cdist == 0.25
    assert ix.manifold_a_key == a_key
    assert ix.manifold_b_key == b_key


def test_registry_add_synthetic_keeps_label_and_registry_id() -> None:
    """A registered synthetic crossing carries the registry id and its label."""
    fp = bare_fixed_point(beta=1 / 3)
    registry = IntersectionRegistry()
    iid = registry.add_synthetic(
        (0.0, 0.0),
        0.0,
        0.0,
        label="anchor",
        manifold_a_key=(fp, "unstable", 0, 0),
        manifold_b_key=(fp, "stable", 0, 0),
    )

    stored = registry[iid]
    assert stored.id == iid
    assert stored.label == "anchor"
    assert stored.is_synthetic
