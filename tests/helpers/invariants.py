"""Reusable geometric/numeric invariant checks for the numerics layer.

These helpers encode the hard physical/numerical laws every manifold must obey
and are shared across the test suite. They deliberately have no ``test_`` prefix
so pytest does not collect them as tests.

The fundamental invariants (see CLAUDE.md):

0. Smooth curve   - a manifold is a simple non-self-intersecting curve, sampled
                    densely; no node juts off it (``assert_no_geometric_spikes``).
                    This is the *observable* correctness property -- a scrambled
                    manifold shows it as a spike even when cdist looks fine.
                    A short scramble (two near points swapped) shows instead as a
                    hairpin reversal of the polyline (``assert_no_reversals``).
1. Monotonicity   - along a manifold's geometric ordering, cdist is
                    non-decreasing (``assert_cdist_monotonic``). It is *not*
                    guaranteed strictly increasing: at a high-stretch fold the
                    refiner bridges sub-ULP gaps with equal-cdist points, so ties
                    are legitimate, so ``strict`` defaults to False. Pass
                    ``strict=True`` only for low-stretch growth that is known
                    to stay injective.
2. Iterate law    - ``c_iterate = stretch_param * c`` along the iterate chain, in
                    the growth direction (``assert_iterate_relation``).
3. One-to-one     - the geometric and iterate linked lists are acyclic and
                    mutually consistent (``assert_one_to_one``).
4. Area preserved - along one iterate chain of crossings in the registry,
                    ``unstable_cdist * stable_cdist`` is invariant
                    (``assert_area_preserved_along_chain``). Only *within* a chain:
                    equal products never imply shared chain membership.

Note:
    A ``Point`` stores a single scalar ``cdist``; a ``BranchPoint`` stores a
    ``(unstable, stable)`` tuple and resolves it via ``get_cdist(stability)``. The
    root fixed point is a ``BranchPoint`` whose ``cdists`` is ``None`` (it is the
    fixed point, distance zero / undefined) and is skipped by these checks. The
    refiner also caches *phoney* pre-iterates (a ``Point`` with ``cdist is None``)
    on the non-growth side; those are skipped too.
"""

from __future__ import annotations

from typing import Optional, TYPE_CHECKING

import numpy as np

from tanglepack import BaseManifold, BranchPoint

if TYPE_CHECKING:
    from tanglepack import IntersectionRegistry


def walk_nodes(manifold: BaseManifold) -> list:
    """Return the manifold's nodes in geometric order (root -> tail), once each."""
    return manifold.get_point_array(return_nodes=True)


def node_cdist(node, stability: str) -> Optional[float]:
    """Resolve a node's cdist for ``stability``, or ``None`` if undefined.

    Returns ``None`` for the root fixed point (a ``BranchPoint`` with no cdists)
    and for phoney cached pre-iterates (a ``Point`` with ``cdist is None``).
    """
    if isinstance(node, BranchPoint):
        if node.cdists is None:
            return None
        value = node.get_cdist(stability)
    else:
        value = node.cdist
    return None if value is None else float(value)


def manifold_cdists(manifold: BaseManifold, stability: str) -> list[float]:
    """Ordered list of defined cdists along the manifold (skips ``None`` nodes)."""
    out = []
    for node in walk_nodes(manifold):
        c = node_cdist(node, stability)
        if c is not None:
            out.append(c)
    return out


def assert_cdist_monotonic(
    manifold: BaseManifold, *, strict: bool = False
) -> None:
    """Assert cdist is non-decreasing (or strictly increasing) along the curve.

    Args:
        manifold: The manifold (or bridge) to check.
        strict: If True, require strictly increasing (no equal neighbours).
            Defaults to False: cdist ties are legitimate at a high-stretch fold
            (see the module docstring). Pass True only for low-stretch growth
            known to stay injective.

    Raises:
        AssertionError: Reporting the first offending index and the two cdists.
    """
    stability = manifold.stability
    cdists = manifold_cdists(manifold, stability)
    for i in range(len(cdists) - 1):
        a, b = cdists[i], cdists[i + 1]
        ok = (b > a) if strict else (b >= a)
        assert ok, (
            f"cdist not {'strictly ' if strict else ''}increasing at index {i}: "
            f"{a!r} -> {b!r} (stability={stability!r}, n={len(cdists)})"
        )


def assert_no_cdist_collision(
    manifold: BaseManifold, *, atol: float = 1e-12
) -> None:
    """Assert no two distinct nodes share a cdist (within ``atol``).

    Sanctioned for low-stretch growth only (``strict`` cdist checks are too).

    Note:
        cdist ties are *legitimate* at a high-stretch fold, where the refiner
        bridges a sub-ULP gap with equal-cdist points (they are spliced
        geometrically, never re-sorted, so they cannot scramble the curve). This
        check therefore applies only to low-stretch growth that is expected to
        stay injective; for the general correctness property use
        ``assert_no_geometric_spikes`` instead.
    """
    stability = manifold.stability
    cdists = sorted(manifold_cdists(manifold, stability))
    for i in range(len(cdists) - 1):
        gap = cdists[i + 1] - cdists[i]
        assert gap > atol, (
            f"cdist collision: two nodes within {atol} at cdist {cdists[i]!r} "
            f"(gap={gap!r}, stability={stability!r})"
        )


def find_geometric_spikes(
    manifold: BaseManifold,
    *,
    length_ratio: float = 30.0,
    edge_skip: int = 3,
) -> list[tuple[int, float]]:
    """Find the over-long polyline segments that mark a scrambled manifold.

    A manifold is a smooth simple curve sampled densely enough that refinement
    keeps every step short. The high-stretch bug instead leaves two geometrically
    distant points adjacent in the linked list (their cdists collapsed to within a
    float ULP and could no longer order them), so the polyline darts across a long
    segment to a misplaced point and back -- the visible "spike". This flags every
    segment whose length exceeds ``length_ratio`` times the median segment length.

    The first/last ``edge_skip`` segments are ignored: the manifold's seed near the
    fixed point is legitimately a single long-but-straight step (the initial
    coarse segment before any refinement), not a scramble.

    Args:
        manifold: The manifold (or bridge) to inspect.
        length_ratio: A segment longer than this many median segment lengths is a
            spike.
        edge_skip: Number of segments at each end to ignore (the coarse seed).

    Returns:
        A list of ``(segment_index, segment_length)`` for each flagged segment.
    """
    nodes = walk_nodes(manifold)
    pts = [np.asarray(n.get_point(), dtype=float) for n in nodes]
    if len(pts) < 2 * edge_skip + 2:
        return []

    seglen = [float(np.linalg.norm(pts[i + 1] - pts[i])) for i in range(len(pts) - 1)]
    positive = [s for s in seglen if s > 0]
    if not positive:
        return []
    median_seg = float(np.median(positive))
    threshold = length_ratio * median_seg

    spikes = []
    for i in range(edge_skip, len(seglen) - edge_skip):
        if seglen[i] > threshold:
            spikes.append((i, seglen[i]))
    return spikes


def assert_no_geometric_spikes(
    manifold: BaseManifold,
    *,
    length_ratio: float = 30.0,
    edge_skip: int = 3,
) -> None:
    """Assert the manifold has no over-long (scrambled) segments (see
    :func:`find_geometric_spikes`).

    Raises:
        AssertionError: Reporting how many spikes were found and the worst one.
    """
    spikes = find_geometric_spikes(
        manifold, length_ratio=length_ratio, edge_skip=edge_skip
    )
    if spikes:
        worst = max(spikes, key=lambda s: s[1])
        pts = [np.asarray(n.get_point(), dtype=float) for n in walk_nodes(manifold)]
        seglen = [float(np.linalg.norm(b - a)) for a, b in zip(pts, pts[1:])]
        median_seg = float(np.median([s for s in seglen if s > 0]))
        raise AssertionError(
            f"{len(spikes)} geometric spike(s) on the {manifold.stability!r} "
            f"manifold: segment {worst[0]} is {worst[1]:.3e} long, "
            f"{worst[1] / median_seg:.1f}x the median segment {median_seg:.3e} "
            f"(threshold {length_ratio:g}x; a misplaced point the polyline "
            f"darts out to and back)."
        )


def find_reversals(
    manifold: BaseManifold, *, max_cos: float = -0.5
) -> list[tuple[int, float]]:
    """Find the hairpin turns that mark two adjacent points stored out of order.

    A refined manifold turns gently from one segment to the next (refinement
    subdivides wherever the curvature area exceeds ``area_cutoff``). When two
    nearby points are spliced in the wrong order -- their cdists collapsed to
    within a float ULP and could no longer order them -- the polyline steps past
    the second point, doubles back to it, and doubles back again: a zig-zag of
    consecutive segments pointing almost opposite ways. This flags every joint
    whose turning cosine is below ``max_cos`` (the default -0.5 is a turn sharper
    than 120 degrees).

    Args:
        manifold: The manifold (or bridge) to inspect.
        max_cos: A joint whose cosine between consecutive segments is below this
            is a reversal.

    Returns:
        A list of ``(joint_index, cosine)`` for each flagged joint (the joint at
        node ``joint_index + 1``).
    """
    pts = np.asarray(
        [n.get_point() for n in walk_nodes(manifold)], dtype=float
    ).reshape(-1, 2)
    if len(pts) < 3:
        return []
    seg = np.diff(pts, axis=0)
    length = np.hypot(seg[:, 0], seg[:, 1])
    keep = (length[1:] > 0) & (length[:-1] > 0)
    cos = np.full(len(seg) - 1, 1.0)
    cos[keep] = np.sum(seg[1:][keep] * seg[:-1][keep], axis=1) / (
        length[1:][keep] * length[:-1][keep]
    )
    return [(int(i), float(cos[i])) for i in np.flatnonzero(cos < max_cos)]


def assert_no_reversals(manifold: BaseManifold, *, max_cos: float = -0.5) -> None:
    """Assert the manifold polyline never doubles back (see :func:`find_reversals`).

    Raises:
        AssertionError: Reporting how many reversals were found and the sharpest.
    """
    reversals = find_reversals(manifold, max_cos=max_cos)
    if reversals:
        joint, cos = min(reversals, key=lambda r: r[1])
        raise AssertionError(
            f"{len(reversals)} hairpin reversal(s) on the {manifold.stability!r} "
            f"manifold; the sharpest at joint {joint} has turning cosine {cos:.3f} "
            f"(threshold {max_cos:g}: two nearby points stored out of order)."
        )


def assert_iterate_relation(
    manifold: BaseManifold, *, rtol: float = 1e-6, atol: float = 1e-12
) -> None:
    """Assert ``c_growth_iterate = stretch_param * c`` along the iterate chain.

    For an unstable manifold the growth iterate is ``next_iterate``; for a stable
    manifold it is ``prev_iterate``. Nodes whose growth iterate is missing, is a
    phoney cached point (``cdist is None``), or that lack a ``stretch_param`` are
    skipped.

    Raises:
        AssertionError: On the first node whose iterate cdist does not match.
    """
    stability = manifold.stability
    grow_attr = "next_iterate" if stability == "unstable" else "prev_iterate"

    for node in walk_nodes(manifold):
        c = node_cdist(node, stability)
        if c is None:
            continue
        stretch = getattr(node, "stretch_param", None)
        if stretch is None:
            continue
        iterate = getattr(node, grow_attr, None)
        if iterate is None:
            continue
        c_iter = node_cdist(iterate, stability)
        if c_iter is None:
            continue
        expected = float(stretch) * c
        assert np.isclose(c_iter, expected, rtol=rtol, atol=atol), (
            f"iterate law violated ({stability}): {grow_attr}.cdist={c_iter!r} "
            f"but stretch_param*cdist={expected!r} "
            f"(stretch={stretch!r}, cdist={c!r}, rtol={rtol})"
        )


def assert_one_to_one(manifold: BaseManifold) -> None:
    """Assert the geometric and iterate linked lists are acyclic and consistent.

    - The geometric walk visits each node at most once (acyclic).
    - Each node's growth iterate links back to it (``next_iterate.prev_iterate``
      is the node for unstable; ``prev_iterate.next_iterate`` for stable), for
      real (cdist-bearing) iterates only.
    """
    stability = manifold.stability
    grow_attr = "next_iterate" if stability == "unstable" else "prev_iterate"
    back_attr = "prev_iterate" if stability == "unstable" else "next_iterate"

    nodes = walk_nodes(manifold)
    ids = [id(n) for n in nodes]
    assert len(ids) == len(set(ids)), (
        f"geometric list is not acyclic: {len(ids) - len(set(ids))} repeated "
        f"node(s) (stability={stability!r})"
    )

    for node in nodes:
        if node_cdist(node, stability) is None:
            continue
        iterate = getattr(node, grow_attr, None)
        if iterate is None or node_cdist(iterate, stability) is None:
            continue
        assert getattr(iterate, back_attr, None) is node, (
            f"iterate back-link broken ({stability}): node's {grow_attr} does not "
            f"point back via {back_attr}"
        )


def assert_area_preserved_along_chain(
    registry: IntersectionRegistry,
    start_id: int,
    *,
    rtol: float = 1e-3,
    max_steps: int = 64,
) -> None:
    """Assert ``unstable_cdist * stable_cdist`` is invariant along an iterate chain.

    Walks the ``n=1`` links of ``registry.iterate_table`` from ``start_id`` and
    checks that the canonical-area product is constant. One forward map step
    stretches the unstable canonical distance and contracts the stable one by the
    same per-step factor, so the product is preserved — this is the registry-level
    expression of area preservation.

    Args:
        registry: The ``IntersectionRegistry`` holding the crossings and the
            iterate table.
        start_id: Registry id of the crossing to start the walk at.
        rtol: Relative tolerance on the product. The default 1e-3 accommodates the
            tolerance of the canonical-distance match that records the links.
        max_steps: Stop after this many forward steps (guards a cyclic table).

    Raises:
        AssertionError: If ``start_id`` has no usable product, or on the first
            step whose product disagrees.

    Note:
        Per CLAUDE.md this only validates invariance *within* a known chain. Two
        different iterate chains generally carry different products, so equal
        products do NOT imply two crossings share a chain — never use the product
        for membership; compare the stable and unstable cdists individually.
    """

    def product(intersection_id: int) -> Optional[float]:
        ix = registry[intersection_id]
        u, s = ix.unstable_cdist, ix.stable_cdist
        if u is None or s is None:
            return None
        return float(u) * float(s)

    p0 = product(start_id)
    assert p0 is not None, (
        f"starting crossing {start_id} has no canonical-area product "
        f"(a cdist is None)"
    )

    node_id = registry.iterate_table[start_id, 1]
    steps = 0
    while node_id is not None and steps < max_steps:
        p = product(node_id)
        if p is not None:
            assert np.isclose(p, p0, rtol=rtol), (
                f"canonical area not preserved along chain from {start_id}: "
                f"crossing {node_id} has product {p!r} vs {p0!r} "
                f"(step {steps + 1}, rtol={rtol})"
            )
        node_id = registry.iterate_table[node_id, 1]
        steps += 1
