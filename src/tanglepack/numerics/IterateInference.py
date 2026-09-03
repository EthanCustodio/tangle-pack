"""Filling the crossing iterate table.

The forward image of a crossing is not searched for by phase-space coordinates --
those drift under the nonlinear map. It is PREDICTED: the image lies on the
branches one orbit step forward (:meth:`~.FixedPoint.FixedPoint.advance_key`) at
canonical distances scaled by each side's own per-map-step factor
(:meth:`~.FixedPoint.FixedPoint.per_step_beta`), and the registry is then asked
for the crossing that sits there.

:class:`~.TangleWorkbench.TangleWorkbench` builds one of these in its constructor
and forwards :meth:`~.TangleWorkbench.TangleWorkbench.infer_iterates` and
:meth:`~.TangleWorkbench.TangleWorkbench.infer_iterate_table` to it, so no caller
changes; :class:`~.BridgeIterator.BridgeIterator` calls
:meth:`IterateInference.register_endpoint_iterates` directly.

Dev Notes:

The two canonical distances are always compared INDIVIDUALLY. Their PRODUCT is
preserved along an iterate chain (area preservation, CLAUDE.md), which makes it
worthless as a discriminator: two different chains can share it, so a
product-based match would happily cross chains.

Each side advances on -- and is scaled by -- its OWN fixed point. On a
heteroclinic crossing the unstable branch belongs to a different periodic orbit
than the stable one, with a different ``k_value``, so one point's per-step factor
cannot serve for both cdists.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Callable, Iterable, Optional

import numpy as np
from numpy.typing import NDArray

from .Intersection import Intersection
from .IntersectionRegistry import IntersectionRegistry

if TYPE_CHECKING:  # pragma: no cover - import cycle guard
    from .Bridge import Bridge
    from .TangleWorkbench import TangleWorkbench

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class IterateInference:
    """
    The iterate-table half of a :class:`~.TangleWorkbench.TangleWorkbench`.

    Holds no state of its own: the intersection registry, the bridge set and the
    dynamical system all live on the workbench, and this object is the set of
    operations that fill the registry's iterate table. Exactly one is constructed
    per workbench, which forwards its two public inference methods here.

    Args:
        workbench: The workbench whose iterate table this fills.
    """

    def __init__(self, workbench: "TangleWorkbench") -> None:
        """
        Args:
            workbench (TangleWorkbench): The workbench whose registry and
                iterate table this fills in.
        """
        self._workbench = workbench

    @property
    def workbench(self) -> "TangleWorkbench":
        """The workbench this reads and writes."""
        return self._workbench

    def infer_iterate_table(self, cdist_rtol: float = 0.05) -> int:
        """
        Scan all iterated bridges and record the n=1 forward iterate relationship for
        each boundary intersection.

        Bridge topology identifies *which* intersections to process (only the two
        endpoints of each iterated bridge). The image f(i_src) is then identified by
        canonical distance and branch identity rather than phase-space coordinates,
        which drift under the nonlinear map and are only approximate. Under one
        application of M the image lies on the branches one orbit step forward of the
        source (see FixedPoint.advance_key), with the unstable cdist stretched and
        the stable cdist contracted by the per-map-step factor
        FixedPoint.per_step_beta.

        The stable branch (manifold_b_key) is the reliable discriminator: the stable
        manifold is always indexed, so every intersection carries it. The unstable
        branch (manifold_a_key) is missing on intersections born from an iterated
        bridge (those bridges are not indexed manifolds), so it is used only as an
        extra constraint when both source and candidate have it. Stable branch plus
        canonical distance resolves the arc-length ambiguity that motivated the old
        coordinate-based match.

        Args:
            cdist_rtol: Maximum relative canonical-distance error between the predicted
                image and a candidate (on the correct stable branch) for the match to
                be accepted. Defaults to 0.05, comfortably above per-step scaling noise.

        Returns:
            Number of new iterate relationships recorded.
        """
        workbench = self._workbench
        registry = workbench._intersection_registry
        recorded = 0

        for bridge in workbench._bridges.values():
            if not bridge.iterated:
                continue
            for src_id in (bridge.first_intersection, bridge.second_intersection):
                if src_id is None:
                    continue
                if self._register_forward_iterate(src_id, registry, cdist_rtol):
                    recorded += 1

        return recorded

    def infer_iterates(self, cdist_rtol: float = 0.05) -> int:
        """
        Fill the iterate table for *every* intersection by canonical-distance mapping.

        For each registered intersection this records its n=1 forward image M(i) — the
        intersection one orbit step forward, found by predicting the image's branches
        (:meth:`FixedPoint.advance_key`) and canonical distances (each scaled by its
        own branch's :meth:`FixedPoint.per_step_beta`) and matching against the
        registry. This generalizes :meth:`infer_iterate_table`,
        which records the same relationship but only for bridge-boundary intersections.

        :meth:`~.TangleWorkbench.TangleWorkbench.compute_intersections` calls this
        automatically, so the table is dense
        as soon as intersections are computed; iterating bridges then keeps it current
        through the bridge machinery. Idempotent: it skips intersections that already
        have an n=1 entry, so repeated calls are cheap.

        Args:
            cdist_rtol: Maximum relative canonical-distance error for a match (default
                0.05, comfortably above per-step scaling noise).

        Returns:
            Number of new iterate relationships recorded.
        """
        workbench = self._workbench
        registry = workbench._intersection_registry
        recorded = 0
        for src_id in registry.all_ids():
            if self._register_forward_iterate(src_id, registry, cdist_rtol):
                recorded += 1
        return recorded

    def _register_forward_iterate(
        self, src_id: int, registry: IntersectionRegistry, cdist_rtol: float
    ) -> bool:
        """
        Record the n=1 forward iterate of one intersection, if it can be identified.

        The image f(src) lies on the branches one map step forward (see
        :meth:`FixedPoint.advance_key`), with each canonical distance scaled by its
        own side's per-map-step factor (:meth:`FixedPoint.per_step_beta`): the
        unstable cdist stretched on the unstable branch's fixed point, the stable
        cdist contracted on the stable branch's. The two coincide on a homoclinic
        crossing; on a HETEROCLINIC one the two branches belong to different
        periodic orbits with different k_values, and using one point's factor for
        both cdists is simply wrong.

        The image is matched by stable branch (always present) plus canonical
        distance, with the unstable branch as an extra constraint when both source
        and candidate carry it (intersections born from iterated bridges have no
        unstable branch). Coordinates are deliberately not used — they drift under
        the nonlinear map, whereas canonical distances and branch keys are exact.

        Returns:
            True if a new iterate edge was recorded, else False.
        """
        workbench = self._workbench
        if (src_id, 1) in registry.iterate_table:
            return False

        source = registry[src_id]
        if self._register_anchor_iterate(src_id, source, registry):
            return True

        # The image's stable branch and stable canonical distance are both
        # predicted, so the candidates are one bisected slice of that branch's
        # own cdist order rather than the whole registry. A crossing outside the
        # window fails the relative stable-cdist test by construction, so the
        # match this narrows to is the same one the full scan would return.
        prediction = self._image_prediction(source)
        if prediction is None:
            return False
        _a_key, b_key, _u_pred, s_pred = prediction
        # The same tolerance _match_image applies below, read once: the window
        # is exactly the set of candidates whose relative stable-cdist error can
        # still come in under cdist_rtol.
        cdist_tol = workbench._intersection_registry.cdist_tol
        window = cdist_rtol * (abs(s_pred) + cdist_tol)
        candidates = registry.candidates_near_stable_cdist(b_key, s_pred, window)

        best_id, _best_err = self._match_image(
            source, candidates, cdist_rtol, exclude={src_id}
        )
        if best_id is None:
            return False

        registry.register_iterate(src_id, 1, best_id)
        return True

    def _register_anchor_iterate(
        self, src_id: int, source: Intersection, registry: IntersectionRegistry
    ) -> bool:
        """
        Record the forward image of an ANCHOR -- a periodic point's own crossing.

        A periodic point sits at canonical distance ``(0, 0)`` on the two branches
        it anchors, and the map carries it to the next point of its orbit, whose
        anchor sits at ``(0, 0)`` on the two ADVANCED branches. So an anchor's
        image is not something to search for by canonical distance: it is fully
        determined by the branch bookkeeping, and
        :meth:`IntersectionRegistry.find` resolves it exactly.

        Doing it here rather than in the generic matcher closes the one case that
        matcher provably cannot handle. On a period-1 orbit the anchor is its own
        image (the periodic point is a fixed point of the map), and the
        "an image is not its source" guard in :meth:`_match_image` -- correct for
        every crossing at cdist ``c > 0``, whose predicted image sits at
        ``beta * c`` far from ``c`` -- excludes the only candidate that can match.
        On a period > 1 orbit the generic search does find the next branch's
        anchor, but by a cdist comparison of ``0`` against ``0``, which is a
        coincidence of the prediction rather than a statement about the orbit.

        Args:
            src_id: Registry id of the crossing being iterated.
            source: The crossing itself.
            registry: The registry to read and write.

        Returns:
            True if the source is an anchor and its image was recorded, else
            False (including for an anchor whose image branch has no anchor yet,
            which happens while a period > 1 orbit is only partly initialized).
        """
        tol = registry.cdist_tol
        if abs(source.unstable_cdist) > tol or abs(source.stable_cdist) > tol:
            return False

        prediction = self._image_prediction(source)
        if prediction is None:
            return False
        a_key, b_key, _u_pred, _s_pred = prediction

        probe = Intersection.synthetic(
            coords=source.coords,
            unstable_cdist=0.0,
            stable_cdist=0.0,
            manifold_a_key=a_key if a_key is not None else source.manifold_a_key,
            manifold_b_key=b_key,
        )
        target = registry.find(probe)
        if target is None:
            logger.debug(
                "Anchor %d has no image anchor on branches (%s, %s) yet",
                src_id,
                a_key,
                b_key,
            )
            return False

        registry.register_iterate(src_id, 1, target)
        return True

    @staticmethod
    def _image_prediction(
        source: Intersection,
    ) -> Optional[tuple[Optional[tuple], tuple, float, float]]:
        """
        What the n=1 image of ``source`` must look like: its two branch keys and
        its two canonical distances.

        Each side advances on -- and is scaled by -- its OWN fixed point: on a
        heteroclinic crossing the unstable branch belongs to a different periodic
        orbit than the stable one, with a different ``k_value``, so using one
        point's per-step factor for both cdists is simply wrong.

        Returns:
            ``(advanced unstable key or None, advanced stable key, predicted
            unstable cdist, predicted stable cdist)``, or None when the crossing
            has no stable branch or its fixed point carries no eigendata.
        """
        if source.manifold_b_key is None:
            return None
        stable_fp = source.manifold_b_key[0]
        if (
            not getattr(stable_fp, "unstable_eigenvalues", None)
            or stable_fp.k_value is None
        ):
            return None

        b_key = stable_fp.advance_key(source.manifold_b_key, 1)
        a_key = (
            source.manifold_a_key[0].advance_key(source.manifold_a_key, 1)
            if source.manifold_a_key is not None
            else None
        )
        unstable_fp = (
            stable_fp if source.manifold_a_key is None else source.manifold_a_key[0]
        )
        u_pred = source.unstable_cdist * unstable_fp.per_step_beta("unstable")
        s_pred = source.stable_cdist * stable_fp.per_step_beta("stable")
        return a_key, b_key, u_pred, s_pred

    def _match_image(
        self,
        source: Intersection,
        candidates: Iterable[tuple[int, Intersection]],
        cdist_rtol: float,
        *,
        exclude: "set[int] | frozenset[int]" = frozenset(),
        accept: Optional[Callable[[Intersection], bool]] = None,
    ) -> tuple[Optional[int], float]:
        """
        Pick the crossing among ``candidates`` that is the n=1 image of ``source``.

        The single place the branch/beta prediction of an iterate is turned into a
        match, shared by the registry-wide heuristic
        (:meth:`_register_forward_iterate`) and the explicit bridge-endpoint
        registration (:meth:`register_endpoint_iterates`). A candidate must sit on
        the advanced STABLE branch -- canonical distance restarts at 0 on every
        branch, so on a period > 1 or heteroclinic tangle the cdist window alone
        would happily pick a crossing on the wrong branch -- and, when both sides
        carry one, on the advanced unstable branch too. The two canonical distances
        are then compared INDIVIDUALLY; their product is never used, since two
        different iterate chains share it (CLAUDE.md, area preservation).

        Args:
            source: The crossing whose forward image is wanted.
            candidates: ``(id, Intersection)`` pairs to choose from.
            cdist_rtol: Maximum relative canonical-distance error for a match.
            exclude: Ids that may not be chosen (the source itself, and any image
                already claimed by another source in the same call).
            accept: Optional extra predicate a candidate must satisfy -- used to
                gate on the real map's image coordinates.

        Returns:
            ``(matched id or None, best relative error seen)``.
        """
        workbench = self._workbench
        prediction = self._image_prediction(source)
        if prediction is None:
            return None, float("inf")

        a_key, b_key, u_pred, s_pred = prediction
        cdist_tol = workbench._intersection_registry.cdist_tol

        best_id, best_err = None, float("inf")
        for tgt_id, tgt in candidates:
            if tgt_id is None or tgt_id in exclude:
                continue
            if tgt.manifold_b_key != b_key:
                continue
            if (
                a_key is not None
                and tgt.manifold_a_key is not None
                and tgt.manifold_a_key != a_key
            ):
                continue
            if accept is not None and not accept(tgt):
                continue
            u_rel = abs(tgt.unstable_cdist - u_pred) / (abs(u_pred) + cdist_tol)
            s_rel = abs(tgt.stable_cdist - s_pred) / (abs(s_pred) + cdist_tol)
            err = max(u_rel, s_rel)
            if err < best_err:
                best_err, best_id = err, tgt_id

        if best_id is not None and best_err <= cdist_rtol:
            return best_id, best_err
        return None, best_err

    def register_endpoint_iterates(
        self,
        bridge: Bridge,
        image_crossings: list[Intersection],
        cdist_rtol: float = 0.05,
    ) -> int:
        """
        Record the n=1 forward iterate of a bridge's two endpoint crossings.

        The image of a bridge's endpoint crossing is a crossing of the image bridge,
        so it is one of the crossings just born on that image -- no registry-wide
        search and no heuristic is needed. The match is made by
        :meth:`_match_image`: the advanced stable (and, when present, unstable)
        branch key, then both canonical distances individually, each scaled by its
        own branch's :meth:`FixedPoint.per_step_beta`.

        Two further guards apply here, where the registry-wide heuristic cannot use
        them:

        * The REAL map is a sanity gate. ``f(source.coords)`` is where the image
          must be, so a candidate further than ``1e-3 * max(1, |f|)`` from it is
          rejected outright. That both catches a cdist coincidence and stops a false
          registration when the true image lies beyond the computed manifold and so
          is not among the crossings at all.
        * The map is injective, so the two endpoints of one bridge cannot share an
          image: a candidate claimed by the first endpoint is off limits to the
          second.

        Args:
            bridge: The parent bridge that was just iterated.
            image_crossings: The crossings born on its forward image, ids assigned.
            cdist_rtol: Maximum relative canonical-distance error for a match.

        Returns:
            Number of new iterate relationships recorded.
        """
        workbench = self._workbench
        if not image_crossings:
            return 0

        registry = workbench._intersection_registry
        candidates = [(ix.id, ix) for ix in image_crossings if ix.id is not None]
        claimed: set[int] = set()
        recorded = 0

        for src_id in (bridge.first_intersection, bridge.second_intersection):
            if src_id is None or src_id not in registry:
                continue
            if (src_id, 1) in registry.iterate_table:
                continue

            source = registry[src_id]
            if source.coords is None:
                continue

            image = np.asarray(
                workbench.dynamical_system.map(np.asarray(source.coords, dtype=float)),
                dtype=float,
            ).ravel()
            gate = 1e-3 * max(1.0, float(np.linalg.norm(image)))

            # image/gate bound as defaults: the predicate is consumed inside this
            # iteration, and binding makes that independent of the loop variable.
            def near_the_true_image(
                candidate: Intersection,
                image: NDArray[np.float64] = image,
                gate: float = gate,
            ) -> bool:
                if candidate.coords is None:
                    return False
                offset = np.asarray(candidate.coords, dtype=float) - image
                return float(np.linalg.norm(offset)) <= gate

            best_id, best_err = self._match_image(
                source,
                candidates,
                cdist_rtol,
                exclude=claimed | {src_id},
                accept=near_the_true_image,
            )

            if best_id is None:
                logger.debug(
                    "no image of crossing %s among the %d crossings born on the "
                    "image of its bridge (best relative cdist error %.3g, %d "
                    "candidate(s) already claimed)",
                    src_id,
                    len(candidates),
                    best_err,
                    len(claimed),
                )
                continue

            registry.register_iterate(src_id, 1, best_id)
            claimed.add(best_id)
            recorded += 1

        return recorded
