"""
Regression (NUM-9, WP10): a bridge whose existing images lie on two curves.

In the nested build, two inner blasts iterate a period-3 bridge that shares its
two straddle nodes with an annulus bridge; an unfiltered outer blast then maps
that annulus bridge, whose images are the grown manifold's tip followed by the
start of its neighbour's image. Merging them as ONE chain walked off the first
curve's tail into ``None`` (``AttributeError`` in ``merge_manifolds``), which the
blast used to skip as a "rare under-resolved iterate". Every bridge the blast
feeds must now be iterated, and its image must run along one curve from the
image of its root to the image of its tail with non-decreasing cdist.
"""

from __future__ import annotations

from cases import build_nested
from tanglepack import BaseManifold


def _zones(case):
    """``(inner, outer)`` zones of a nested build."""
    zones = sorted(case.zones, key=lambda zone: zone.area)
    return zones[0], zones[-1]


def _image_is_one_curve(bridge) -> bool:
    """Whether the walk from ``f(root)`` reaches ``f(tail)``, cdist non-decreasing."""
    first, last = bridge.root.next_iterate, bridge.tail.next_iterate
    branch = bridge.fixed_point.advance_key(bridge.manifold_key, 1)[3]
    walker = BaseManifold(
        first, "unstable", bridge.stretch_param, bridge.fixed_point, branch_index=branch, manifold_key=None
    )
    cdists, previous, node = [], None, first
    while node is not None and node is not last:
        cdists.append(float(node.cdist))
        previous, node = node, walker.walk_fwd(previous, node)
    return node is last and cdists == sorted(cdists) and cdists[-1] <= last.cdist


def test_an_unfiltered_outer_blast_iterates_every_annulus_bridge():
    case = build_nested(outer_blasts=0, through="zones")
    session = case.session
    inner, outer = _zones(case)
    for _ in range(2):
        session.blast_zone(inner, num_iterations=1, fixed_point=None, min_separation=1e-4, diagnose=False)

    result = session.blast_zone(outer, num_iterations=1, fixed_point=None, min_separation=1e-4, diagnose=False)

    parents = result.steps[0].interior_parents
    assert any(bridge.fixed_point is not outer.key for bridge in parents), "an annulus bridge must be fed"
    for bridge in parents:
        assert bridge.iterated, bridge.id
        assert _image_is_one_curve(bridge), bridge.id
