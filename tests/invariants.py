"""Reusable geometric/numeric invariant checks for the numerics layer.

These helpers encode the hard physical/numerical laws every manifold must obey
and are shared across the test suite. They deliberately have no ``test_`` prefix
so pytest does not collect them as tests.

The fundamental invariants (see CLAUDE.md):

0. Smooth curve   - a manifold is a simple non-self-intersecting curve, sampled
                    densely; no node juts off it (``assert_no_geometric_spikes``).
                    This is the *observable* correctness property -- a scrambled
                    manifold shows it as a spike even when cdist looks fine.
1. Monotonicity   - along a manifold's geometric ordering, cdist is
                    non-decreasing (``assert_cdist_monotonic``). It is *not*
                    guaranteed strictly increasing: at a high-stretch fold the
                    refiner bridges sub-ULP gaps with equal-cdist points, so ties
                    are legitimate. Pass ``strict=True`` only for low-stretch
                    growth that is known to stay injective.
2. Iterate law    - ``c_iterate = stretch_param * c`` along the iterate chain, in
                    the growth direction (``assert_iterate_relation``).
3. One-to-one     - the geometric and iterate linked lists are acyclic and
                    mutually consistent (``assert_one_to_one``).
4. Area preserved - for an intersection chain ``unstable_cdist * stable_cdist`` is
                    invariant (``assert_area_preserved_along_chain``).

Note:
    A ``Point`` stores a single scalar ``cdist``; a ``BranchPoint`` stores a
    ``(unstable, stable)`` tuple and resolves it via ``get_cdist(stability)``. The
    root fixed point is a ``BranchPoint`` whose ``cdists`` is ``None`` (it is the
    fixed point, distance zero / undefined) and is skipped by these checks. The
    refiner also caches *phoney* pre-iterates (a ``Point`` with ``cdist is None``)
    on the non-growth side; those are skipped too.
"""

from __future__ import annotations

from typing import Literal, Optional

import numpy as np

from tanglepack import BaseManifold, BranchPoint, Point


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
    manifold: BaseManifold, *, strict: bool = True
) -> None:
    """Assert cdist is (strictly) increasing along the geometric ordering.

    Args:
        manifold: The manifold (or bridge) to check.
        strict: If True, require strictly increasing (no equal neighbours). This
            is what catches cdist collisions that leaked out of a merge.

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
        raise AssertionError(
            f"{len(spikes)} geometric spike(s) on the {manifold.stability!r} "
            f"manifold: segment {worst[0]} is {worst[1]:.3e} long, "
            f"{worst[1] / (length_ratio):.3e} = {length_ratio:g}x the median "
            f"(a misplaced point the polyline darts out to and back)."
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
    branch_point: BranchPoint, *, rtol: float = 1e-4, max_steps: int = 64
) -> None:
    """Assert ``unstable_cdist * stable_cdist`` is invariant along an iterate chain.

    Walks the ``next_iterate`` chain of intersection ``BranchPoint``s starting at
    ``branch_point`` and checks the canonical-area product is constant. Per
    CLAUDE.md this only validates invariance *within* a known chain; equal
    products do NOT imply two points share a chain.
    """
    def product(bp: BranchPoint) -> Optional[float]:
        if not isinstance(bp, BranchPoint) or bp.cdists is None:
            return None
        u = bp.get_cdist("unstable")
        s = bp.get_cdist("stable")
        if u is None or s is None:
            return None
        return float(u) * float(s)

    p0 = product(branch_point)
    assert p0 is not None, "starting branch point has no canonical-area product"

    node = branch_point.next_iterate
    steps = 0
    while node is not None and steps < max_steps:
        p = product(node)
        if p is not None:
            assert np.isclose(p, p0, rtol=rtol), (
                f"canonical area not preserved along chain: {p!r} vs {p0!r} "
                f"(step {steps + 1}, rtol={rtol})"
            )
        node = node.next_iterate
        steps += 1
