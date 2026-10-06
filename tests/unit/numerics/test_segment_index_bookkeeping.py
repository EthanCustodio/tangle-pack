"""Tangle segment-index bookkeeping and point-insertion side effects.

Two regression pins:
1. ``Tangle.add_manifold`` must never record ``None`` as a segment id.
   ``_insert_segment`` returns ``None`` for an already-indexed edge, and a
   ``None`` in ``_manifold_segs`` poisons every later ``_seg_lookup[sid]``
   walk (``_segments_touching``, the trim paths) with a ``KeyError``.
2. Inserting a point into a linked list must never overwrite the host
   point's ``stretch_param``. The inserted node inherits a missing value
   instead (a fresh crossing separator carries ``stretch_param=None``).
"""

from __future__ import annotations

from tanglepack import Point


def test_add_manifold_twice_records_no_none_ids(small_tangle):
    workbench, fp = small_tangle
    tangle = workbench.Tangle

    # Re-adding purges and re-indexes; a second add of the same manifold is
    # exactly the case where _insert_segment starts returning None.
    for manifold in workbench.manifolds.values():
        tangle.add_manifold(manifold)
        tangle.add_manifold(manifold)

    for manifold, seg_ids in tangle._manifold_segs.items():
        assert None not in seg_ids, f"None segment id recorded for {manifold}"
        # every recorded id must resolve, which is what downstream walks rely on
        for sid in seg_ids:
            assert sid in tangle._seg_lookup


def test_insert_point_keeps_host_stretch_param():
    host = Point(0.0, 0.0, cdist=0.0, stretch_param=2.5)
    separator = Point(1.0, 0.0, cdist=1.0)  # no stretch_param, like a crossing separator

    host.insert_point_forward(separator)

    assert host.stretch_param == 2.5, "host stretch_param was overwritten"
    assert separator.stretch_param == 2.5, "inserted node did not inherit"


def test_insert_point_does_not_clobber_existing_node_param():
    host = Point(0.0, 0.0, cdist=0.0, stretch_param=2.5)
    node = Point(1.0, 0.0, cdist=1.0, stretch_param=3.5)

    host.insert_point_backward(node)

    assert host.stretch_param == 2.5
    assert node.stretch_param == 3.5


def test_insert_iterate_touches_no_stretch_param():
    """Iterate-list inserts leave stretch_param alone on both sides.

    Callers construct iterate points with their stretch already set (or fill
    it explicitly afterwards), and the host may even be a BranchPoint, which
    has no scalar stretch_param at all.
    """
    host = Point(0.0, 0.0, cdist=1.0, stretch_param=2.5)
    image = Point(0.5, 0.5, cdist=2.0)

    host.insert_next_iterate(image)

    assert host.stretch_param == 2.5
    assert image.stretch_param is None

    preimage = Point(-0.5, -0.5, cdist=0.5)
    host.insert_prev_iterate(preimage)

    assert host.stretch_param == 2.5
    assert preimage.stretch_param is None
