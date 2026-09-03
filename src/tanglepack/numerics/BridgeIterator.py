"""Forward iteration of bridges, and the derived bridge genealogy.

A bridge is a stretch of unstable manifold between two crossings, so its forward
image is a stretch of the branch one map step forward. This module owns both
halves of that statement: :meth:`BridgeIterator.iterate_bridge`, which computes
the image and cuts it, and :meth:`BridgeIterator.image_bridges` /
:meth:`BridgeIterator.preimage_bridges`, which DERIVE the genealogy from the
crossing iterate table rather than storing it on the bridge.

:class:`~.TangleWorkbench.TangleWorkbench` builds one of these in its constructor
and forwards the four public methods to it, so no caller changes.

Dev Notes:

Cutting an iterated bridge (plan row 2.2) only ever cuts at the crossings just
born on that image -- the ones whose ``unstable_manifold`` IS the image. A crossing
found on an edge the image SHARES with an already-indexed curve is registered but
not cut at here, because its bracketing points lie on the other curve; the cut that
uses it is the one made on that curve. Widening the cut set to every crossing whose
cdist falls in the image's span would reintroduce the cross-curve bridges plan row
1.1 removed, so it must not be done without a per-curve span test.

Partial bridges (an image piece bounded by fewer than two crossings) have no
endpoint pair, hence no ``BridgeId``, hence no identity to dedupe on. They are
kept out of ``_bridges`` in a separate ``_partial_bridges`` list: they are real
stretches of unstable manifold the blast frontier must carry forward, but they are
not bridges in the topological sense and never appear in a genealogy answer. Two
consequences to keep in mind:

* repeated blasting of the same region can accumulate several Bridge objects over
  the same partial arc (nothing identifies them as the same arc), which is why the
  nested period-3 blast script registers 141 bridges where the pre-2.2 code
  reported 140;
* ``image_bridges`` can therefore lose a mid-arc leading/trailing piece of an
  image, so the derived frontier is the bridges of the image, not the whole image.

Both are by definition, not a gap to be papered over with a cdist-only signature --
that guesswork is exactly what plan row 2.2 removed.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Optional

from .Bridge import Bridge, BridgeId
from .Intersection import ManifoldKey

if TYPE_CHECKING:  # pragma: no cover - import cycle guard
    from .TangleWorkbench import TangleWorkbench

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class BridgeIterator:
    """
    The bridge-iteration half of a :class:`~.TangleWorkbench.TangleWorkbench`.

    Holds no state of its own: the bridge registry, the intersection registry and
    the manifold machine all live on the workbench, and this object is the set of
    operations that read and extend them. Exactly one is constructed per workbench,
    which forwards its four public methods here.

    Args:
        workbench: The workbench whose bridges this iterates.
    """

    def __init__(self, workbench: "TangleWorkbench") -> None:
        """
        Args:
            workbench (TangleWorkbench): The workbench whose bridges, registry
                and manifolds this iterator reads and extends.
        """
        self._workbench = workbench

    @property
    def workbench(self) -> "TangleWorkbench":
        """The workbench this iterator operates on."""
        return self._workbench

    def image_bridges(
        self, bridge_id: BridgeId, n: int = 1
    ) -> Optional[list[BridgeId]]:
        """
        The bridges tiling the ``n``-th forward image of one bridge.

        Genealogy is DERIVED, not stored: a bridge spans the arc between crossings
        ``a`` and ``b``, so its image spans the arc between ``f^n(a)`` and
        ``f^n(b)`` on the branch ``n`` map steps forward
        (:meth:`FixedPoint.advance_key`). The bridges of that image are exactly the
        registered bridges of the image branch whose two endpoints both lie inside
        that span. Both iterates are read from the registry's iterate table -- no
        cdist guesswork -- so the answer is only as complete as the table.

        Args:
            bridge_id: The bridge whose image is wanted.
            n: Number of map steps; negative walks backward.

        Returns:
            The image bridges' ids in increasing unstable canonical distance, or
            ``None`` when either endpoint has no registered ``n``-iterate. The list
            can be empty when the image arc is grown but not yet cut.

        Raises:
            KeyError: If ``bridge_id`` is not registered.

        Note:
            Partial pieces have no id and are never returned, so an image whose
            leading (or trailing) stretch runs past its last crossing loses that
            stretch here. That is by definition: a partial arc is not a bridge.
        """
        workbench = self._workbench
        return self._image_bridge_ids(workbench.bridge(bridge_id), n)

    def preimage_bridges(
        self, bridge_id: BridgeId, n: int = 1
    ) -> Optional[list[BridgeId]]:
        """
        The bridges tiling the ``n``-th BACKWARD image of one bridge.

        Exactly :meth:`image_bridges` with the sign of ``n`` flipped; see there for
        the semantics and the ``None`` case.

        Args:
            bridge_id: The bridge whose preimage is wanted.
            n: Number of backward map steps (positive means ``f^-n``).

        Returns:
            The preimage bridges' ids in increasing unstable canonical distance, or
            ``None`` when either endpoint has no registered ``-n``-iterate.

        Raises:
            KeyError: If ``bridge_id`` is not registered.
        """
        workbench = self._workbench
        return self._image_bridge_ids(workbench.bridge(bridge_id), -n)

    def _image_span(
        self, bridge: Bridge, n: int
    ) -> Optional[tuple[tuple[int, int], tuple[float, float], ManifoldKey]]:
        """The ``n``-iterates of a bridge's endpoints, their cdists, and the image key.

        The two ids come back ORDERED BY UNSTABLE CDIST (low, high), matching the
        cdists returned alongside them. Which of ``f^n(a)``, ``f^n(b)`` is the lower
        is not assumed: an inversion branch or a bookkeeping slip could swap them,
        and a caller comparing an ordered chain against them must not silently invert.
        """
        workbench = self._workbench
        if bridge.id is None:
            return None
        table = workbench._intersection_registry.iterate_table
        images = [table[endpoint, n] for endpoint in bridge.id]
        if any(image is None for image in images):
            return None

        registry = workbench._intersection_registry
        images.sort(key=lambda image: float(registry[image].unstable_cdist))
        lo, hi = (float(registry[image].unstable_cdist) for image in images)
        image_key = bridge.fixed_point.advance_key(bridge.manifold_key, n)
        return (images[0], images[1]), (lo, hi), image_key

    def _image_bridge_ids(self, bridge: Bridge, n: int) -> Optional[list[BridgeId]]:
        """The ids of the registered bridges inside a bridge's n-th image span."""
        workbench = self._workbench
        if n == 0:
            return [bridge.id] if bridge.id is not None else None

        span = self._image_span(bridge, n)
        if span is None:
            return None
        _images, (lo, hi), image_key = span

        registry = workbench._intersection_registry
        tol = registry.cdist_tol
        inside: list[tuple[float, BridgeId]] = []
        for bid, other in workbench._bridges.items():
            if other.manifold_key != image_key:
                continue
            first = float(registry[bid[0]].unstable_cdist)
            second = float(registry[bid[1]].unstable_cdist)
            if first >= lo - tol and second <= hi + tol:
                inside.append((first, bid))

        return [bid for _cdist, bid in sorted(inside)]

    def _existing_image_bridges(self, bridge: Bridge) -> Optional[list[Bridge]]:
        """
        The already-computed bridges that tile this bridge's forward image, if the
        image is entirely covered by them.

        A bridge spanning crossings ``a -> b`` images onto the arc
        ``f(a) -> f(b)`` of the branch one map step forward
        (:meth:`FixedPoint.advance_key`). If that arc has already been grown and
        cut, those bridges ARE the forward image -- there is no need (and it is
        actively harmful) to re-map and re-cut it: re-mapping lays a second,
        slightly different polyline over curve that already exists (the "zig-zag"),
        spawns overlapping duplicate bridges, and trips the same-stability
        (unstable x unstable) detector along the near-coincident pair.

        The image is only reusable when the existing bridges TILE it exactly: the
        first starts at ``f(a)``, the last ends at ``f(b)``, and each shares an
        endpoint with the next. Anything less means part of the image is curve that
        has not been computed, which is precisely the case that must fall through
        to real iteration.

        Args:
            bridge: The bridge about to be iterated.

        Returns:
            The existing image bridges in unstable cdist order, or ``None`` to fall
            through to real iteration.
        """
        workbench = self._workbench
        span = self._image_span(bridge, 1)
        if span is None:
            logger.debug(
                "bridge %s: no reusable image -- an endpoint has no registered n=1 "
                "iterate, so the image will be recomputed by mapping it forward",
                bridge.id,
            )
            return None
        (lo_image, hi_image), _cdists, _key = span

        image_ids = self._image_bridge_ids(bridge, 1)
        if not image_ids:
            logger.debug(
                "bridge %s: image arc %s -> %s carries no registered bridge; "
                "recomputing",
                bridge.id,
                lo_image,
                hi_image,
            )
            return None

        # The chain must run from f(a) to f(b) with no gap: the endpoint ids come
        # from _image_span already ordered by cdist, as image_ids is.
        if image_ids[0][0] != lo_image or image_ids[-1][1] != hi_image:
            logger.debug(
                "bridge %s: image arc %s -> %s is only partly cut (%s); recomputing",
                bridge.id,
                lo_image,
                hi_image,
                image_ids,
            )
            return None
        if any(lhs[1] != rhs[0] for lhs, rhs in zip(image_ids, image_ids[1:])):
            logger.debug(
                "bridge %s: image bridges %s do not tile the arc; recomputing",
                bridge.id,
                image_ids,
            )
            return None

        return [workbench._bridges[bid] for bid in image_ids]

    def iterate_bridge(self, bridge: Bridge) -> list[Bridge]:
        """
        Map a bridge forward one iterate, add the result to the tangle, detect new
        intersections with the stable manifold, cut the result into new bridges,
        and return those bridges.

        Marks the original bridge as iterated. The genealogy is not stored on the
        bridge: ``workbench.image_bridges(bridge.id)`` derives it from the crossing
        iterate table afterwards.

        Args:
            bridge: A bridge created by create_bridges() or a previous iterate_bridge().

        Returns:
            List of new Bridge objects from cutting the iterated result. If the
            image crosses the stable manifold fewer than twice, returns a
            single-element list holding the unsplit image, marked
            :attr:`Bridge.partial` (it is bounded by fewer than two crossings, so
            it is not a bridge in the topological sense).

        Raises:
            ValueError: If bridge has already been iterated.
            ValueError: If create_bridges() has not been called yet.
        """
        workbench = self._workbench
        if bridge.iterated:
            raise ValueError(
                "This bridge has already been iterated. Its image is "
                "workbench.image_bridges(bridge.id)."
            )
        if not workbench._bridges and not workbench._partial_bridges:
            raise ValueError(
                "No bridges registered. Call create_bridges() before iterate_bridge()."
            )

        # 0. If this bridge's forward image is a stretch of manifold that has already
        #    been grown, that section is already cut into bridges -- reuse those
        #    existing objects instead of mapping the points forward again. Re-mapping
        #    would lay a second, slightly-different polyline over curve that already
        #    exists (the "zig-zag"), spawn overlapping duplicate bridges, and trip the
        #    same-stability (unstable x unstable) detector along the near-coincident
        #    pair. Only when the image runs past the grown extent is there genuinely
        #    new curve to compute (handled below).
        existing_image = self._existing_image_bridges(bridge)
        if existing_image is not None:
            bridge.iterated = True
            workbench._bump_generation()
            return existing_image

        # 1. map forward
        iterated = workbench._man_machine.iterate_bridge(bridge)

        # The iterated bridge is M(bridge): a segment of the unstable manifold one
        # orbit step forward of the parent's branch. Recording that key -- before the
        # image is indexed -- is what lets every crossing detected on it carry its
        # unstable branch identity (and lets its children advance the key again).
        image_key = bridge.fixed_point.advance_key(bridge.manifold_key, 1)
        iterated.manifold_key = image_key

        # 2. register with the tangle. The bridge is unstable manifold, so every
        #    curve it can legitimately cross is stable and already in the rtree
        #    from compute_intersections -- querying finds all its crossings, and
        #    skipping the rtree inserts keeps repeated blasting fast.
        workbench.Tangle.add_manifold(iterated, index_segments=False)

        # 3. resolve only new crossings involving the iterated bridge
        new_intersections = workbench.Tangle.populate_intersections_for_manifold(
            iterated
        )

        for ix in new_intersections:
            # A crossing that collides with one already registered IS that crossing;
            # take the registry's id for it so the cut below records the canonical id.
            ix.id = workbench._intersection_registry.add(ix)

        # 4. cut at crossings. Only the crossings just born on this image are cut at:
        #    they alone carry bracketing points that lie on the image's own polyline.
        new_bridges = (
            workbench.Tangle.create_bridges(
                new_intersections,
                for_manifold=iterated,
                cdist_tol=workbench._intersection_registry.cdist_tol,
            )
            if new_intersections
            else []
        )

        # A bridge is bounded by two crossings, so create_bridges yields nothing when
        # the image crosses the stable manifold fewer than twice -- either zero times
        # (no crossing) or exactly once (a single new homoclinic point cannot bound a
        # bridge on its own). Either way the forward image is still a real piece of
        # unstable manifold whose dynamics must be carried forward, so keep it as one
        # PARTIAL bridge (both endpoint ids stay None) rather than dropping it. (Any
        # new crossing was already registered above; it is picked up as a bridge
        # boundary on later iterates, once a second crossing appears alongside it.)
        # ManifoldMachine.iterate_bridge may return a BaseManifold, so wrap it as a
        # Bridge -- downstream consumers (uniiterated_bridges, genealogy) need the
        # Bridge attributes.
        if not new_bridges:
            if isinstance(iterated, Bridge):
                new_bridges = [iterated]
            else:
                new_bridges = [
                    Bridge(
                        root=iterated.root,
                        stability=iterated.stability,
                        stretch_param=iterated.stretch_param,
                        fixed_point=iterated.fixed_point,
                        tail=iterated.tail,
                        branch_index=iterated.branch_index,
                        manifold_key=image_key,
                    )
                ]

        # 5. record the forward iterate of the parent's two endpoint crossings
        workbench._iterate_inference.register_endpoint_iterates(
            bridge, new_intersections
        )

        # 6. Single-copy invariant: a bridge is uniquely defined by the two
        #    intersections it connects, so a freshly cut child whose BridgeId is
        #    already registered IS that bridge -- _register_bridge hands back the
        #    stored object instead of a duplicate. The caller still sees every piece
        #    of the image; it just points at the one stored copy. Iterating the
        #    fixed-point bridge, for instance, re-traces curve already held, so its
        #    children all resolve to existing bridges and nothing is added.
        children = [workbench._register_bridge(child) for child in new_bridges]
        bridge.iterated = True
        # Mapping forward can refine (and cut at) curve that already exists --
        # the image of a bridge whose points all have iterates is laid over
        # grown manifold -- so every memoised walk is dropped, not just the
        # children's.
        workbench._bump_generation(invalidate_walks=True)

        return children

    def iterate_all_bridges(self) -> list[Bridge]:
        """
        Iterate all bridges that have not yet been mapped forward.

        Returns:
            All new bridges produced across all iterations.
        """
        workbench = self._workbench
        pending = list(
            workbench.uniiterated_bridges
        )  # snapshot before loop mutates _bridges

        all_new: list[Bridge] = []
        for bridge in pending:
            all_new.extend(workbench.iterate_bridge(bridge))

        return all_new
