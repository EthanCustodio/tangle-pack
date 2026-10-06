"""The physical-law checks of the law tier, one function per law.

Every check has the signature ``check_<law>(case) -> int``: it asserts its law
on a built :class:`cases.Case` and returns HOW MANY ITEMS it checked, so the
law tier can insist that a law was not vacuous (a count of 0 is a failure
unless the ``(law, case)`` pair is listed in ``cases.NOT_APPLICABLE``). The
checks only read the case; they never mutate the session.

The checks are grouped into :data:`LAYERS` (manifolds, crossings, anchors,
bridges, ...). :func:`run_layer` runs every check of a layer, collects the
failures and raises ONE ``AssertionError`` naming each failing check.

Dev Notes:

* Phase 1 of the 2026-10-05 test-suite refactor lands the numerics layers
  (manifolds, crossings, anchors, bridges). The topology layers (partition,
  arrangement, classes, iterated, dual graph, symbolic) are added with the
  law-tier modules in Phase 4, which also deletes the scattered twins.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Callable

import pytest

from helpers.invariants import (
    assert_cdist_monotonic,
    assert_iterate_relation,
    assert_no_geometric_spikes,
    assert_one_to_one,
)

if TYPE_CHECKING:  # pragma: no cover - typing only
    from cases import Case

#: Relative slack of the per-link area check: the library's own iterate
#: identification slack (``collision_rtol`` default of ``StrongPip`` and
#: ``Pseudoneighbor``). Measured worst link at baseline: 1.4e-3 (p3, nested),
#: 2.3e-4 (k28 two blasts), 1.8e-4 (k10).
AREA_LINK_RTOL: float = 1e-2

#: A law check: asserts on a case and returns the number of items checked.
LawCheck = Callable[["Case"], int]


def _anchor_ids(case: "Case") -> list[int]:
    """Registry ids of the synthetic anchors (both cdists exactly 0, no segments)."""
    return [
        intersection_id
        for intersection_id, ix in case.registry
        if ix.seg_ids is None
        and ix.unstable_cdist is not None
        and ix.stable_cdist is not None
        and float(ix.unstable_cdist) == 0.0
        and float(ix.stable_cdist) == 0.0
    ]


# --------------------------------------------------------------------------- #
# Manifolds
# --------------------------------------------------------------------------- #
def check_manifold_cdist_monotone(case: "Case") -> int:
    """cdist is non-decreasing along every manifold (ties allowed at folds)."""
    manifolds = list(case.workbench.manifolds.values())
    for manifold in manifolds:
        assert_cdist_monotonic(manifold, strict=False)
    return len(manifolds)


def check_manifold_no_spikes(case: "Case") -> int:
    """No manifold carries a geometric spike (a scrambled node)."""
    manifolds = list(case.workbench.manifolds.values())
    for manifold in manifolds:
        assert_no_geometric_spikes(manifold)
    return len(manifolds)


def check_manifold_iterate_law(case: "Case") -> int:
    """``c_iterate = stretch_param * c`` along every manifold's iterate chain."""
    manifolds = list(case.workbench.manifolds.values())
    for manifold in manifolds:
        assert_iterate_relation(manifold)
    return len(manifolds)


def check_manifold_one_to_one(case: "Case") -> int:
    """The geometric and iterate lists of every manifold are acyclic and consistent."""
    manifolds = list(case.workbench.manifolds.values())
    for manifold in manifolds:
        assert_one_to_one(manifold)
    return len(manifolds)


# --------------------------------------------------------------------------- #
# Crossings
# --------------------------------------------------------------------------- #
def check_crossings_unstable_by_stable(case: "Case") -> int:
    """Every registered crossing is one unstable key (a) times one stable key (b)."""
    count = 0
    for intersection_id, ix in case.registry:
        assert ix.manifold_a_key is not None and ix.manifold_b_key is not None, intersection_id
        assert ix.manifold_a_key[1] == "unstable", (intersection_id, ix.manifold_a_key)
        assert ix.manifold_b_key[1] == "stable", (intersection_id, ix.manifold_b_key)
        count += 1
    return count


def check_crossings_cdists_defined(case: "Case") -> int:
    """Every registered crossing has two non-negative canonical distances."""
    count = 0
    for intersection_id, ix in case.registry:
        assert ix.unstable_cdist is not None and float(ix.unstable_cdist) >= 0.0, intersection_id
        assert ix.stable_cdist is not None and float(ix.stable_cdist) >= 0.0, intersection_id
        count += 1
    return count


def check_crossings_area_along_chains(
    case: "Case", *, rtol: float = AREA_LINK_RTOL
) -> int:
    """``unstable_cdist * stable_cdist`` is preserved by every registered ``+1`` link.

    Compared LINK BY LINK, not against the head of the chain, so the error
    of a long chain (period 3 runs a dozen links) does not accumulate. The
    anchors (product 0) are skipped.
    """
    registry = case.registry
    anchors = set(_anchor_ids(case))
    count = 0
    for intersection_id, ix in registry:
        image_id = registry.iterate_table[intersection_id, 1]
        if image_id is None or intersection_id in anchors:
            continue
        image = registry[image_id]
        before = float(ix.unstable_cdist) * float(ix.stable_cdist)
        after = float(image.unstable_cdist) * float(image.stable_cdist)
        assert after == pytest.approx(before, rel=rtol), (
            f"area not preserved on the link {intersection_id} -> {image_id}: "
            f"{before!r} -> {after!r}"
        )
        count += 1
    return count


# --------------------------------------------------------------------------- #
# Anchors
# --------------------------------------------------------------------------- #
def check_one_anchor_per_unstable_branch(case: "Case") -> int:
    """Every unstable branch with a registered anchor has exactly one (author rule)."""
    per_branch: dict = {}
    for intersection_id in _anchor_ids(case):
        key = case.registry[intersection_id].manifold_a_key
        per_branch.setdefault(key, []).append(intersection_id)
    for key, ids in per_branch.items():
        assert len(ids) == 1, f"{len(ids)} anchors on unstable branch {key}: {ids}"
    return len(per_branch)


# --------------------------------------------------------------------------- #
# Bridges
# --------------------------------------------------------------------------- #
def check_bridge_partial_iff_no_id(case: "Case") -> int:
    """A bridge is partial exactly when it has no ``BridgeId``."""
    bridges = case.workbench.bridges
    for bridge in bridges:
        assert bridge.partial == (bridge.id is None), bridge
    return len(bridges)


def check_bridge_id_in_unstable_order(case: "Case") -> int:
    """``BridgeId = (first, second)`` is ordered by unstable cdist, on the bridge's branch."""
    registry = case.registry
    count = 0
    for bridge in case.workbench.bridges:
        if bridge.id is None:
            continue
        first, second = (registry[i] for i in bridge.id)
        assert float(first.unstable_cdist) <= float(second.unstable_cdist), bridge.id
        assert first.manifold_a_key == bridge.manifold_key, bridge.id
        assert second.manifold_a_key == bridge.manifold_key, bridge.id
        count += 1
    return count


def check_bridges_at_exact(case: "Case") -> int:
    """``bridges_at(x)`` lists exactly the identified bridges ending at ``x``."""
    workbench = case.workbench
    expected: dict[int, set] = {}
    for bridge in workbench.bridges:
        if bridge.id is None:
            continue
        for end in bridge.id:
            expected.setdefault(end, set()).add(bridge.id)
    count = 0
    for intersection_id, _ix in case.registry:
        assert set(workbench.bridges_at(intersection_id)) == expected.get(
            intersection_id, set()
        ), intersection_id
        count += 1
    return count


def check_bridge_single_copy(case: "Case") -> int:
    """Each ``BridgeId`` is registered once and ``bridge(id)`` returns that copy."""
    workbench = case.workbench
    identified = [bridge for bridge in workbench.bridges if bridge.id is not None]
    ids = [bridge.id for bridge in identified]
    assert len(ids) == len(set(ids)), "a BridgeId is registered twice"
    for bridge in identified:
        assert workbench.bridge(bridge.id) is bridge, bridge.id
    return len(identified)


#: The law checks by layer.
LAYERS: dict[str, list[LawCheck]] = {
    "manifolds": [
        check_manifold_cdist_monotone,
        check_manifold_no_spikes,
        check_manifold_iterate_law,
        check_manifold_one_to_one,
    ],
    "crossings": [
        check_crossings_unstable_by_stable,
        check_crossings_cdists_defined,
        check_crossings_area_along_chains,
    ],
    "anchors": [check_one_anchor_per_unstable_branch],
    "bridges": [
        check_bridge_partial_iff_no_id,
        check_bridge_id_in_unstable_order,
        check_bridges_at_exact,
        check_bridge_single_copy,
    ],
}


def run_layer(case: "Case", layer: str) -> dict[str, int]:
    """
    Run every check of one layer and report the item count of each.

    Args:
        case: The built case.
        layer: A key of :data:`LAYERS`.

    Returns:
        ``{check name: items checked}``.

    Raises:
        AssertionError: One or more checks failed; the message names each
            failing check with its own message.
    """
    counts: dict[str, int] = {}
    failures: list[str] = []
    for check in LAYERS[layer]:
        try:
            counts[check.__name__] = check(case)
        except AssertionError as error:
            failures.append(f"{check.__name__}: {error}")
    if failures:
        raise AssertionError(
            f"{len(failures)} law check(s) of layer {layer!r} failed on case "
            f"{case.name!r}:\n" + "\n".join(failures)
        )
    return counts
