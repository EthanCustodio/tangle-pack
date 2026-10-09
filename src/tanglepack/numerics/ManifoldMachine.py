"""
Growing manifolds one iteration of the map at a time.

:class:`ManifoldMachine` is the numerical engine of the library. It maps a
manifold's points forward (backward, for a stable manifold) one step, splices
the images into the geometric linked list in canonical-distance order, merges
the new arc with what was already there, and refines the result wherever the
curvature area between adjacent points exceeds ``area_cutoff``.

Invariants it maintains:

* The geometric list is ordered by canonical distance. That order is
  non-decreasing, not strictly increasing: at a high-stretch fold adjacent
  cdists can collapse to within a float ULP, and a tie is harmless because
  points are spliced geometrically and never re-sorted by cdist.
* Every point it inserts keeps its iterate links, so a grown curve can always
  be walked back to the seed segment.
* Refinement is breadth-first and batched: one map evaluation per layer, not
  one per inserted point.

Dev Notes:

NOT A BUG -- refinement "explosion" at high development is just exponential
manifold growth into the map's escape region, and we never need to grow that far.

Symptom: at the k=2.1 inner period-3 fp ``construct_fixed_point([[0,1],[-1,0],
[-1,1]])`` with a tight area_cutoff (1e-7), a manifold arm grows gradually for a
while and then the point count and runtime blow up in a single iteration (in one
seen case to ~2e5 points; pushing one further hit millions and tens of GB).

Diagnosis (instrumented, 2026-06): the blow-up is NOT a tight-fold refinement
artifact (the earlier guess in this note). It is the arm ESCAPING to infinity. Two
signals settle it, traced per growth iteration on both the stable [1,1] and the
unstable [0,1] arm of that fp:
  * The minimum segment length never shrinks -- it sits at the ~1.4e-6 seed
    spacing the whole time. A tight fold would drive it toward zero; there is no
    tight fold.
  * The point count tracks coordinate divergence exactly. For the unstable [0,1]
    arm: iter 8 -> 478 pts, max|coord| ~ 3; iter 9 -> ~9.9k pts, max|coord| ~ 2e3;
    iter 10 -> ~1.5e6 pts, max|coord| ~ 2e26. The arm flies off to infinity, and
    refinement faithfully (and pointlessly) tries to resolve segments that span an
    ever-larger, diverging curve -- each inserted preiterate-midpoint maps even
    farther out, a divergence feedback, not fold resolution.

So the manifold is simply growing exponentially and escaping, which is the
expected behaviour of these maps for arms that leave the trapped region (cf.
``test_high_stretch_period3_growth``, which stops at iter 4 noting "the arm escapes
by the fifth iterate"). There is nothing to fix in ``refine_manifold`` /
``_refine_layer``: a fully developed tangle needs only a handful of iterations, and
arms are meant to be grown with a stop condition (``grow_until_turnaround`` /
``grow_until_arclength`` / ``grow_until_intersection``), NOT a large fixed
iteration count that walks an escaping arm out to infinity. The ``< 1e-8`` spatial
guard in ``_refine_layer`` stays as-is; it is fine for the regime we actually use.
See project memory ``solver-rootcause-and-seed-step``.

The point budget (2026-10-07, after the numerics stress campaign): a fixed
iteration count has no notion of escape, and past the escape step the point
count can jump 5-8x in ONE refinement step -- the E3 last-factor prediction
underestimated it by an order of magnitude -- so a step could run out of memory
before anything outside it got a say. Every top-level call (one growth step of
a whole ``k_value`` chain, one bridge image, one initializer image) therefore
counts the points it creates and raises
:class:`~.budget.PointBudgetExceeded` BEFORE the map evaluation that would
cross :attr:`ManifoldMachine.max_points_per_step`. The workbench and the
initializer share ONE machine, so its ``area_cutoff`` and budget reach every
refinement, the initial fundamental segments of a ``k_value > 1`` orbit
included. Prefer the ``grow_until_*`` drivers (``grow_until_bounded`` stops on
``max|x|``, the only condition shown to stop before take-off) over
``grow_n_times`` with a large count.

The default cutoff is 1e-7 (2026-10-07, was 1e-4): the stress campaign's E4
found it the cheapest cutoff converged on every case tested (at 1e-4 the
period-3 crossing count still moves with the cutoff at fixed depth; a coarse
polygon both misses true crossings and makes spurious ones). It costs ~8-12x
the points of 1e-4 per step, so the escape cliff sits ~10x closer in memory:
the point budget is what makes the finer default safe. A recipe tuned at 1e-4
(the k=10 tests and scripts grow 9-10 steps) passes ``area_cutoff=1e-4``
explicitly; the topology there is the same at both.

Open questions:

* ``iterate_manifold`` maps exactly one step and takes no ``num_iterations``
  input. An old TODO asked for one; repeated iteration is the caller's loop
  (``grow_x_times``), so the extra parameter is not wanted.
"""

from __future__ import annotations

from contextlib import contextmanager
from typing import Iterator, Literal, Optional
import logging

import numpy as np

from .Point import Point
from .FixedPoint import FixedPoint
from .DynamicalSystem import DynamicalSystem
from .BaseManifold import BaseManifold
from .ManifoldView import ManifoldView
from .BranchPoint import BranchPoint
from .Bridge import Bridge
from .Intersection import ManifoldKey, Stability
from .budget import DEFAULT_MAX_POINTS_PER_STEP, PointBudgetExceeded
from .geometry import curvature_area_batch

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

#: Library default for :attr:`ManifoldMachine.area_cutoff` (1e-4 before
#: 2026-10-07; see the module Dev Notes).
DEFAULT_AREA_CUTOFF = 1e-7


class ManifoldMachine:
    """
    Grows manifolds: one iteration of the map, spliced in and refined.

    A machine holds no manifold state: apart from the system it maps with, the
    refinement threshold and the point budget, it only keeps the tally of the
    current top-level call. Every method takes the manifold it acts on. The
    three things it does are:

    * :meth:`iterate_manifold` -- map every point of a curve one step in its own
      dynamical direction (forward for an unstable manifold, backward for a
      stable one), splice the images into the geometric linked list, and merge
      the result with the part of the curve that already existed.
    * :meth:`refine_manifold` -- insert points wherever the curvature area
      between two neighbours exceeds :attr:`area_cutoff`, breadth-first and one
      batched map call per layer.
    * :meth:`grow_x_times` -- drive the two above around a periodic orbit.

    Invariants every method preserves:

    * The geometric list stays ordered by canonical distance (non-decreasing --
      see the module Dev Notes on ties at a high-stretch fold).
    * Inserted points keep their iterate links, so a grown curve can always be
      walked back to the seed segment.
    * A refined point is spliced geometrically between its neighbours; the list
      is never re-sorted.

    Attributes:
        system (DynamicalSystem): The map, its inverse and their batched forms.
        area_cutoff (float): Curvature-area threshold above which a pair of
            adjacent points gets a new point between them.
        max_points_per_step (Optional[int]): Points one top-level call may
            create before it raises :class:`~.budget.PointBudgetExceeded`;
            ``None`` disables the budget.
        points_created (int): Points created by the current (or last)
            top-level call.
        max_abs_coord (float): Largest ``|coordinate|`` among the points the
            current (or last) top-level call mapped or created.
    """

    def __init__(
        self,
        system: DynamicalSystem,
        area_cutoff: float = DEFAULT_AREA_CUTOFF,
        max_points_per_step: Optional[int] = DEFAULT_MAX_POINTS_PER_STEP,
    ) -> None:
        """
        Bind a machine to the system it will map with.

        Args:
            system (DynamicalSystem): The dynamical system to grow with.
            area_cutoff (float): Refinement threshold (see the class docs);
                :data:`DEFAULT_AREA_CUTOFF` (1e-7) by default.
            max_points_per_step (Optional[int]): Per-call point budget; ``None``
                disables it.
        """

        self.system = system
        self.area_cutoff = area_cutoff
        self.max_points_per_step = max_points_per_step
        self.points_created = 0
        self.max_abs_coord = 0.0
        self._call_depth = 0

    @contextmanager
    def _top_level_call(self) -> Iterator[None]:
        """Reset the tally on entry to the OUTERMOST machine call only, so a
        growth chain is one budget and consecutive calls never share one."""
        if self._call_depth == 0:
            self.points_created = 0
            self.max_abs_coord = 0.0
        self._call_depth += 1
        try:
            yield
        finally:
            self._call_depth -= 1

    def _note_coords(self, coords: np.ndarray) -> None:
        """Fold coordinates into :attr:`max_abs_coord` (NaN ignored, inf kept)."""
        self.max_abs_coord = float(
            np.fmax.reduce(np.abs(coords).ravel(), initial=self.max_abs_coord)
        )

    def _charge(
        self, num_points: int, stage: Literal["map", "refine"], manifold: BaseManifold
    ) -> None:
        """Count ``num_points`` about-to-be-created points, raising first if they
        would take the call over its budget."""
        budget = self.max_points_per_step
        if budget is not None and self.points_created + num_points > budget:
            raise PointBudgetExceeded(
                budget,
                self.points_created,
                stage,
                self.max_abs_coord,
                manifold.manifold_key,
            )
        self.points_created += num_points

    def _branch_view(
        self,
        fixed_point: FixedPoint,
        key: ManifoldKey,
        stretch_param: float,
    ) -> BaseManifold:
        """A transient BaseManifold over the piece a manifold key names."""
        view = BaseManifold(
            fixed_point.branch_points[key[2]],
            key[1],
            stretch_param=stretch_param,
            fixed_point=fixed_point,
            branch_index=key[3],
            manifold_key=key,
        )

        view.root = view.first_node()

        return view

    def new_grow_manifold(
        self,
        fixed_point: FixedPoint,
        stability: Stability,
        branch_index: Optional[int] = None,
    ) -> None:
        """
        Grow every manifold of the given stability by one map iteration.

        The manifold pieces attached to a periodic orbit form ONE chain under the
        map: the piece at ``(orbit, branch)`` maps onto the piece one
        :meth:`FixedPoint.advance_key` step along. Growing them all by one
        iteration is therefore a single walk of that chain, iterating each piece
        in turn -- ``k_value`` steps, which closes the chain (``k_value`` is the
        period without inversion and twice the period with it).

        Growth happens by side effect: :meth:`iterate_manifold` splices the images
        into the next orbit point's linked list, so the manifold it returns is
        never needed here -- the next step re-reads the spliced-in points through a
        fresh view.

        The walk starts one step BEHIND ``(0, b)`` so that the piece whose image
        lands on orbit point 0 is iterated first; it keeps every piece growing by
        exactly one fundamental segment per call.

        Args:
            fixed_point (FixedPoint): The fixed point whose manifolds are grown.
            stability (Stability): Stability to grow. Unstable
                pieces advance along the chain (+1 map step), stable pieces are
                grown by the inverse map and so walk it backwards (-1).
            branch_index (int, optional): Which eigenvector branch to grow. None
                grows every allocated branch. On an inversion point every branch
                index names the same chain, which covers both branches.

        Note:
            A branch index names an EIGENDIRECTION, and each one is a chain of
            ``k_value`` pieces walked by ``advance_key``. An inversion chain
            visits both branches, so it is walked once from branch 0 whatever is
            asked.

        Raises:
            PointBudgetExceeded: If the step would create more than
                :attr:`max_points_per_step` points over the whole chain. The
                chain is then partly mapped and under-refined, not resumable.
        """
        step = 1 if stability == "unstable" else -1

        if fixed_point.check_inversion():
            start_branches = [0]
        else:
            start_branches = [0, 1] if branch_index is None else [branch_index]

        with self._top_level_call():
            for b in start_branches:
                key = fixed_point.advance_key((fixed_point, stability, 0, b), -step)

                view = self._branch_view(fixed_point, key, stretch_param=1)
                if view.root is None:
                    continue  # branch not initialized

                # the whole chain shares the stretch parameter of its first piece
                stretch_param = view.root.stretch_param
                view.stretch_param = stretch_param

                for _ in range(fixed_point.k_value):

                    self.iterate_manifold(view)

                    key = fixed_point.advance_key(key, step)
                    view = self._branch_view(fixed_point, key, stretch_param)
                    assert view.root is not None, (
                        f"manifold piece {key[2], key[3]} of a chain that started "
                        f"at branch {b} is not initialized"
                    )

    def iterate_manifold(self, manifold: BaseManifold) -> BaseManifold:
        """
        Iterates all uniterated points in a manifold
        Returns a new manifold pointing only to the iterate

        Parameters:
            manifold: manifold to iterate

        Raises:
            PointBudgetExceeded: If mapping or refining the image would take the
                current top-level call over :attr:`max_points_per_step`. The
                parent's points may already be linked to a half-refined image
                (the links are made before refinement); that state is not
                resumable.
        """
        with self._top_level_call():
            return self._iterate_manifold(manifold)

    @staticmethod
    def _manifold_from_points(
        points: list[Point], stability: Stability, stretch_param: float, fixed_point: FixedPoint
    ) -> BaseManifold:
        """A key-less manifold through fresh points in cdist order, linked along the
        geometric list (forward on an unstable branch, backward on a stable one)."""
        manifold = BaseManifold(points[0], stability, stretch_param, fixed_point, manifold_key=None)
        for point, following in zip(points, points[1:]):
            if stability == "unstable":
                point.insert_point_forward(following)
            else:
                point.insert_point_backward(following)
        manifold.tail = points[-1]
        return manifold

    def _iterate_manifold(self, manifold: BaseManifold) -> BaseManifold:
        """The body of :meth:`iterate_manifold`, inside the caller's budget scope."""
        viewer = ManifoldView(manifold, self.system)

        non_iterated_coords = manifold.get_non_iterated_point_array()
        non_iterated_cdists = manifold.get_non_iterated_cdist_array()
        non_iterated_points = manifold.get_non_iterated_point_array(return_nodes=True)

        old_points = manifold.get_iterated_point_array(return_nodes=True)

        if len(non_iterated_coords):

            self._note_coords(non_iterated_coords)
            self._charge(len(non_iterated_coords), "map", manifold)
            iterated_points = viewer.map_fwd_batch(non_iterated_coords)
            self._note_coords(iterated_points)
            distances = np.ravel(manifold.stretch_param * non_iterated_cdists)

            xvals = iterated_points[:, 0]
            yvals = iterated_points[:, 1]

            new_points = [
                Point(x, y, float(cdist), stretch_param=manifold.stretch_param)
                for x, y, cdist in zip(xvals, yvals, distances)
            ]

            for i, point in enumerate(non_iterated_points):
                if manifold.stability == "unstable":
                    point.insert_next_iterate(new_points[i])
                else:
                    point.insert_prev_iterate(new_points[i])

            # The image of a piece lives on the branch one map step forward, which
            # only the caller can name (see iterate_bridge); the machine leaves the
            # key unset rather than copying the un-advanced one.
            new_iterated_points = self._manifold_from_points(
                new_points, manifold.stability, manifold.stretch_param, manifold.fixed_point
            )

        # The images that already exist may lie on several curves (NUM-9): merge
        # each contiguous chain of them in turn.
        chains = self._image_chains(manifold, old_points)

        if len(non_iterated_coords):
            mapped_manifold = new_iterated_points
            for chain in chains:
                mapped_manifold = self.merge_manifolds(chain, mapped_manifold)
            self.refine_manifold(mapped_manifold)

            # cdist is the manifold's ordering key, so the geometric list must be
            # cdist-sorted. It is only *non-decreasing*, not strictly increasing:
            # where a high-stretch fold has collapsed adjacent cdists to within a
            # float ULP the refiner bridges the gap with equal-cdist points (a tie
            # is harmless -- such points are spliced geometrically, never re-sorted).
            assert (np.diff(mapped_manifold.get_cdist_array(), axis=0) >= 0).all()

            return mapped_manifold

        # Every point already has its image (a bridge cut from the middle of a grown
        # manifold): the image is bracketed by those images directly; do NOT step one
        # more iterate forward (that image does not exist yet). Images on two curves
        # are joined first, as above.
        if not chains:
            raise ValueError("No points in manifold to iterate")
        image = chains[0]
        for chain in chains[1:]:
            image = self.merge_manifolds(image, chain)
        if len(chains) > 1:
            self.refine_manifold(image)

        piece = dict(
            root=image.root,
            stability=manifold.stability,
            stretch_param=manifold.stretch_param,
            fixed_point=manifold.fixed_point,
            tail=image.tail,
            branch_index=manifold.branch_index,
            manifold_key=None,
        )
        if isinstance(manifold, Bridge):
            return Bridge(name=manifold.name, **piece)
        return BaseManifold(**piece)

    def _image_chains(self, manifold: BaseManifold, images: list[Point]) -> list[BaseManifold]:
        """
        The existing images of a curve's points, as geometrically contiguous chains.

        A bridge's images are usually one run of one curve, but a bridge whose
        straddle nodes it shares with a neighbour that was iterated first has its
        images on TWO curves (NUM-9: the grown manifold's tip, then the start of the
        neighbour's image). Walking the curve from the first image, the images it
        reaches (to the end of that curve, or past the last image's cdist) form one
        chain from the first to the last reached; the next chain starts at the
        first image after it. An image the walk passes over without reaching (a
        node a merge dropped as a duplicate cdist keeps its preimage's link) is
        skipped, as one merged chain always skipped it.

        Args:
            manifold: The curve being iterated (its stability and branch walk the images).
            images: The existing images, in the curve's order.

        Returns:
            One manifold per chain, in cdist order; empty when there are no images.
        """
        def chain(first: Point, last: Point) -> BaseManifold:
            return BaseManifold(
                first,
                manifold.stability,
                manifold.stretch_param,
                manifold.fixed_point,
                tail=last,
                branch_index=manifold.branch_index,
                manifold_key=None,
            )

        if not len(images):
            return []
        walker = chain(images[0], images[-1])
        bound = images[-1].cdist
        chains, i = [], 0
        while i < len(images):
            reached, previous, node = set(), None, images[i]
            while node is not None and node.cdist <= bound:
                reached.add(id(node))
                previous, node = node, walker.walk_fwd(previous, node)
            last = max(j for j in range(i, len(images)) if id(images[j]) in reached)
            chains.append(chain(images[i], images[last]))
            i = last + 1
        return chains

    def iterate_bridge(self, manifold: Bridge) -> BaseManifold:
        """
        Iterates a bridge forward and returns another bridge.

        Args:
            manifold (Bridge): Bridge to iterate forward.

        Returns:
            The forward image, with ``manifold_key`` left unset. The image lies on
            the branch one map step forward of the parent's, so writing the
            parent's (un-advanced) key here would mislabel it; the caller
            (:meth:`TangleWorkbench.iterate_bridge`) advances the key and sets it
            before the image is indexed.

        Raises:
            PointBudgetExceeded: As :meth:`iterate_manifold`; one bridge image
                is one budget.
        """

        return self.iterate_manifold(manifold)

    def grow_x_times(
        self,
        fixed_point: FixedPoint,
        stability: Stability,
        num_times: int = 1,
        branch_index: Optional[int] = None,
    ) -> None:
        """
        Grows all the manifolds of the given stability from the fixed point
        by iterating them num_times and merging them back together

        Parameters:
            fixed_point: fixed point to grow the manifolds from
            stability: stability of the manifolds to grow
            num_times: number of times iterating the manifolds
        """

        for _ in range(num_times):
            self.new_grow_manifold(fixed_point, stability, branch_index)

    def merge_manifolds(
        self, manifold_1: BaseManifold, manifold_2: BaseManifold
    ) -> BaseManifold:
        """
        O(n) inplace method for merging two linked lists (manifold_1 and manifold_2)
        Merges two manifolds based on their cdists.
        Does not return a 3rd manifold, returns the whichever manifold has the smaller
        starting cdist with the new tail set

        Parameters:
            manifold_1: manifold to merge
            manifold_2: manifold to merge
        """

        head_1 = manifold_1.root
        head_2 = manifold_2.root

        # NOTE If there are negative cdists this merge method will fail
        # Choose which manifold to start with based on the lower cdist
        if head_1.cdist <= head_2.cdist:
            current_point = head_1
            head_1 = manifold_1.walk_fwd(None, head_1)
            output_manifold = manifold_1
        else:
            current_point = head_2
            head_2 = manifold_2.walk_fwd(None, head_2)
            output_manifold = manifold_2

        # set the termination point right after the manifold tail
        over_one_1 = manifold_1.walk_fwd(None, manifold_1.tail)
        over_one_2 = manifold_2.walk_fwd(None, manifold_2.tail)

        while head_1 is not over_one_1 and head_2 is not over_one_2:

            if head_1.cdist < head_2.cdist:
                next_head = manifold_1.walk_fwd(None, head_1)
                if manifold_1.walk_fwd(None, current_point) is not head_1:
                    self._insert_point_geometrically(current_point, head_1, manifold_1)
                current_point = head_1
                head_1 = next_head

            elif head_2.cdist < head_1.cdist:
                next_head = manifold_2.walk_fwd(None, head_2)
                if manifold_2.walk_fwd(None, current_point) is not head_2:
                    self._insert_point_geometrically(current_point, head_2, manifold_2)
                current_point = head_2
                head_2 = next_head

            else:  # they share a cdist
                if head_1 is head_2:
                    # The same physical node appears in both lists.
                    if manifold_1.walk_fwd(None, current_point) is not head_1:
                        self._insert_point_geometrically(
                            current_point, head_1, manifold_1
                        )
                    head_1 = manifold_1.walk_fwd(None, head_1)
                    head_2 = manifold_2.walk_fwd(None, head_2)
                    current_point = manifold_1.walk_fwd(None, current_point)
                else:
                    # Two different points collided onto the same cdist. cdist is
                    # the only ordering key, so once it saturates (one ULP at large
                    # magnitude exceeds the true separation near a fold) it can no
                    # longer order the pair -- keeping both and ordering by the tied
                    # key scrambles them into a spike. The newly generated point has
                    # no iterate links yet; keep the established (iterate-linked)
                    # point and skip the duplicate.
                    logger.debug(
                        "Duplicate cdist collision at cdist=%s; coords %s vs %s",
                        head_1.cdist,
                        head_1.get_point(),
                        head_2.get_point(),
                    )
                    head_1_has_iter = (
                        head_1.next_iterate is not None
                        or head_1.prev_iterate is not None
                    )
                    head_2_has_iter = (
                        head_2.next_iterate is not None
                        or head_2.prev_iterate is not None
                    )
                    if head_1_has_iter or not head_2_has_iter:
                        # head_1 is the keeper (or neither has iterates — pick head_1)
                        next_head_1 = manifold_1.walk_fwd(None, head_1)
                        if manifold_1.walk_fwd(None, current_point) is not head_1:
                            self._insert_point_geometrically(
                                current_point, head_1, manifold_1
                            )
                        current_point = head_1
                        head_1 = next_head_1
                        head_2 = manifold_2.walk_fwd(None, head_2)
                    else:
                        # head_2 is the keeper
                        next_head_2 = manifold_2.walk_fwd(None, head_2)
                        if manifold_2.walk_fwd(None, current_point) is not head_2:
                            self._insert_point_geometrically(
                                current_point, head_2, manifold_2
                            )
                        current_point = head_2
                        head_2 = next_head_2
                        head_1 = manifold_1.walk_fwd(None, head_1)

        # One list ran out first: splice the exhausted list's last kept node onto
        # the remaining head and take that manifold's tail as the merged tail, so
        # the whole remainder (not just the next node) stays on the curve.
        if head_1 is not over_one_1:
            self._insert_point_geometrically(
                current_point, head_1, manifold_1, only_forward=True
            )
            output_manifold.tail = manifold_1.tail
        else:
            self._insert_point_geometrically(
                current_point, head_2, manifold_2, only_forward=True
            )
            output_manifold.tail = manifold_2.tail

        return output_manifold

    def refine_manifold(
        self,
        manifold: BaseManifold,
        branch_index: Optional[int] = None,
        final_node: Optional[Point] = None,
    ) -> set[Point]:
        """
        Adds additional points in areas of the manifold with high curvature.

        Refines breadth-first: every consecutive pair of points is a candidate.
        For each pair a left and right parabola plus a chord are fit and the
        bounded area is measured; pairs whose area exceeds ``self.area_cutoff``
        get a new point inserted at the (forward image of the) midpoint of their
        preiterates. A pair that is refined spawns its two child pairs, which are
        checked in the next layer. All of a layer's midpoints are mapped forward
        in a single batched call, so the dynamical map is evaluated once per
        layer rather than once per inserted point.

        Parameters:
            manifold (BaseManifold): current manifold
            branch_index: Optional if starting at a fixed point
            final_node: Unused; kept for backwards-compatible signature.

        Note:
            The per-pair decision depends only on local geometry and preiterates,
            so the breadth-first traversal reaches the same ``area < cutoff``
            stopping condition as the old depth-first scheme.
        """
        logger.debug("Refining manifold: %r", manifold.get_point_array())

        num_initial_points = (
            len(manifold.get_point_array())
            if logger.isEnabledFor(logging.INFO)
            else None
        )

        viewer = ManifoldView(manifold, self.system)

        # Layer 0: every consecutive geometric pair from root to tail.
        layer: list[tuple] = []
        previous_point = manifold.root
        current_point = manifold.walk_fwd(None, previous_point)
        while current_point is not None:
            layer.append((previous_point, current_point))
            if current_point is manifold.tail:
                break
            next_point = manifold.walk_fwd(previous_point, current_point)
            previous_point, current_point = current_point, next_point

        modified_points: set = set()
        while layer:
            layer = self._refine_layer(
                layer, manifold, viewer, branch_index, modified_points
            )

        if num_initial_points is not None:
            num_final_points = len(manifold.get_point_array())
            logger.debug(
                "%d points added during refinement (%d total)",
                num_final_points - num_initial_points,
                num_final_points,
            )

        return modified_points

    def _refine_layer(
        self,
        layer: list[tuple],
        manifold: BaseManifold,
        viewer: ManifoldView,
        branch_index,
        modified_points: set,
    ) -> list[tuple]:
        """
        Refine one breadth-first layer of pairs and return the next layer.

        Computes curvature areas for the whole layer at once, inserts a new
        point into each pair that needs refining (mapping all midpoints forward
        in a single batched call), and returns the child pairs of the refined
        pairs as the next layer.

        Args:
            layer: list of ``(p0, p1)`` node pairs to consider.
            manifold: the manifold being refined.
            viewer: view binding the manifold to the dynamical system.
            branch_index: branch index for inserting after a BranchPoint.
            modified_points: set accumulating every touched/created point.

        Returns:
            The next layer of ``(p0, p1)`` pairs (empty when refinement is done).
        """
        stability = manifold.stability

        # Gather node references and per-pair geometry, dropping pairs that are
        # degenerate or missing the data needed to refine (mirrors the old
        # per-pair guards: dead endpoints, near-coincident points, and missing
        # preiterates on bridge boundary points).
        p0_nodes, p1_nodes = [], []
        p0_xy, p1_xy, left_xy, right_xy = [], [], [], []
        pre_mid, c0_list, c1_list = [], [], []

        nan_row = np.array([np.nan, np.nan])

        for p0, p1 in layer:
            if p0 is None or p1 is None:
                continue

            a = p0.get_point()
            b = p1.get_point()
            if np.hypot(b[0] - a[0], b[1] - a[1]) < 1e-8:
                continue

            pre0 = ManifoldMachine._get_preiterate(p0, stability, 1)
            pre1 = ManifoldMachine._get_preiterate(p1, stability, 1)
            if pre0 is None or pre1 is None:
                logger.debug(
                    "Skipping refinement for pair at cdist %.3g-%.3g: "
                    "preiterate missing on one or both points.",
                    p0.cdist,
                    p1.cdist,
                )
                continue

            left = manifold.walk_back(p1, p0)
            right = manifold.walk_fwd(p0, p1)

            p0_nodes.append(p0)
            p1_nodes.append(p1)
            p0_xy.append(a)
            p1_xy.append(b)
            left_xy.append(left.get_point() if left is not None else nan_row)
            right_xy.append(right.get_point() if right is not None else nan_row)
            pre_mid.append(0.5 * (pre0.get_point() + pre1.get_point()))
            c0_list.append(float(p0.cdist))
            c1_list.append(float(p1.cdist))

        if not p0_nodes:
            return []

        areas = self._curvature_area_batch(
            np.array(p0_xy),
            np.array(p1_xy),
            np.array(left_xy),
            np.array(right_xy),
        )

        # Refine every pair whose curvature area exceeds the cutoff. The new point
        # takes the cdist midpoint of its neighbours; at a high-stretch fold the
        # curve has expanded so much that adjacent points are far apart in space
        # while their cdists have collapsed to within a float ULP, so that midpoint
        # may round exactly onto an endpoint. That tie is harmless here -- the point
        # is spliced geometrically between its neighbours, not re-sorted by cdist --
        # and bridging the gap is exactly what keeps the curve smooth instead of
        # leaving a coarse spike. The ``< 1e-8`` spacing guard above guarantees the
        # subdivision still terminates once the two points effectively coincide.
        c0_arr = np.array(c0_list)
        c1_arr = np.array(c1_list)
        mid_cd = 0.5 * (c0_arr + c1_arr)

        flagged = np.flatnonzero(areas >= self.area_cutoff)
        if flagged.size == 0:
            return []

        self._charge(flagged.size, "refine", manifold)
        pre_mid = np.array(pre_mid)
        new_coords = viewer.map_fwd_batch(pre_mid[flagged])
        self._note_coords(new_coords)
        new_cdists = mid_cd[flagged]

        next_layer: list[tuple] = []
        for out_idx, pair_idx in enumerate(flagged):
            p0 = p0_nodes[pair_idx]
            p1 = p1_nodes[pair_idx]
            x, y = new_coords[out_idx]
            new_point = Point(
                x, y, float(new_cdists[out_idx]), stretch_param=manifold.stretch_param
            )
            # Cache the preiterate so the child pairs can be refined next layer.
            self._cache_preiterate(new_point, pre_mid[pair_idx], stability)

            self._insert_point_geometrically(p0, new_point, manifold, branch_index)
            modified_points.update((p0, p1, new_point))

            next_layer.append((p0, new_point))
            next_layer.append((new_point, p1))

        return next_layer

    def _get_refined_point(
        self,
        p0: Point | BranchPoint,
        p1: Point | BranchPoint,
        viewer: ManifoldView,
        stability: Stability,
    ):
        """
        Reference single-pair refinement: the midpoint of the two preiterates
        mapped one step forward. The batched _refine_layer is the live path;
        this scalar version is kept as the pinned reference implementation.
        """

        p0_preiterate = self._get_preiterate(p0, stability)
        p1_preiterate = self._get_preiterate(p1, stability)

        new_point_coords_back = 0.5 * (
            p1_preiterate.get_point() + p0_preiterate.get_point()
        )

        new_distance = 0.5 * (float(p0.cdist) + float(p1.cdist))

        new_point_coords = self._get_iterate(new_point_coords_back, viewer)

        x = new_point_coords[0]
        y = new_point_coords[1]
        new_point = Point(x, y, float(new_distance), stretch_param=p0.stretch_param)

        self._cache_preiterate(new_point, new_point_coords_back, stability)

        return new_point

    @staticmethod
    def _get_iterate(point: Point | BranchPoint, viewer: ManifoldView):
        """Map a point one step forward in the manifold's stability direction."""
        return viewer.map_fwd(point)

    @staticmethod
    def _get_preiterate(
        point: Point | BranchPoint,
        stability: Stability,
        num_iterates: int = 1,
    ):
        """
        Helper function to get a point's preiterate based on stability

        Parameters:
            point (Point or BranchPoint): point to get preiterate
            stability (str): stability of the manifold
        """

        if stability == "unstable":
            return point.get_prev_iterate(num_iterates)
        else:
            return point.get_next_iterate(num_iterates)

    @staticmethod
    def _cache_preiterate(
        point: Point | BranchPoint,
        preiterate_coords,
        stability: Stability,
    ):
        """
        Helper function to cache preiterate based on stability.
        Only stores the coordinates to be accessed for refinement.

        Parameters:
            point (Point or BranchPoint): Current Points
            preiterate_coords: Preiterate coordinates to cache
            stability (str): stability of the manifold
        """

        # phoney cached preiterate point without cdist
        cached_preiterate = Point(preiterate_coords[0], preiterate_coords[1])

        if stability == "unstable":
            point.prev_iterate = cached_preiterate
        else:
            point.next_iterate = cached_preiterate

    @staticmethod
    def _insert_point_geometrically(
        p0: Point | BranchPoint,
        new_point: Point | BranchPoint,
        manifold: BaseManifold,
        branch_index: int = None,
        only_forward: bool = False,
    ):
        """
        Takes in a point (p0) on the manifold and inserts new_point after it.
        In this case we are inserting away from the fixed point. On the
        unstable manifold that is 'forward'. On the stable manifold that
        is 'backward'.

        Parameters:
            p0 (Point): Anchor point we wil insert after
            new_point (Point): Point to insert into the manifold
            manifold (BaseManifold): Current manifold
            branch_index (int): Optional parameter if inserting after fixed point

        Note:
            Use this method to smoothly handle insertion regardless of stability
            and BranchPoints

        REFACTOR:
            branch_index may not be necessary if that information is already in manifold
        """

        if isinstance(p0, BranchPoint):
            if manifold.stability == "unstable":
                if p0.forward_branches[branch_index] is new_point:
                    return
                p0.insert_point_forward(
                    new_point, branch_index=branch_index, only_forward=only_forward
                )
            else:
                if p0.backward_branches[branch_index] is new_point:
                    return
                p0.insert_point_backward(
                    new_point, branch_index=branch_index, only_forward=only_forward
                )

        else:
            if manifold.stability == "unstable":
                if p0.forward is new_point:
                    return
                p0.insert_point_forward(new_point, only_forward=only_forward)
            else:
                if p0.backward is new_point:
                    return
                p0.insert_point_backward(new_point, only_forward=only_forward)

    @staticmethod
    def _curvature_area_batch(
        p0: np.ndarray,
        p1: np.ndarray,
        left: np.ndarray,
        right: np.ndarray,
    ) -> np.ndarray:
        """The curvature area of each segment (:func:`.geometry.curvature_area_batch`)."""
        return curvature_area_batch(p0, p1, left, right)

