"""
Growing manifolds one iteration of the map at a time.

:class:`ManifoldMachine` is the numerical engine of the library. It maps a
manifold's points forward (backward, for a stable manifold) one step, puts the
images on a chain of the point store (:class:`~.PointStore.Chain`), merges that
chain with the run of the curve that was already there by canonical distance,
and refines the result wherever the curvature area between adjacent points
exceeds ``area_cutoff``. Every step works on chain positions: nothing walks a
curve point by point.

Invariants it maintains:

* Every chain is ordered by canonical distance. That order is
  non-decreasing, not strictly increasing: at a high-stretch fold adjacent
  cdists can collapse to within a float ULP, and a tie is harmless because
  points are inserted at a position and never re-sorted by cdist.
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

Diagnosis (instrumented): the blow-up is NOT a tight-fold refinement artifact.
It is the arm ESCAPING to infinity. Two signals settle it, traced per growth
iteration on both the stable [1,1] and the unstable [0,1] arm of that fp:
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

The point budget: a fixed
iteration count has no notion of escape, and past the escape step the point
count can jump 5-8x in ONE refinement step -- the last-factor prediction
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

The default cutoff is 1e-7: the numerics stress campaign
found it the cheapest cutoff converged on every case tested (at 1e-4 the
period-3 crossing count still moves with the cutoff at fixed depth; a coarse
polygon both misses true crossings and makes spurious ones). It costs ~8-12x
the points of 1e-4 per step, so the escape cliff sits ~10x closer in memory:
the point budget is what makes the finer default safe. A recipe tuned at 1e-4
(the k=10 tests and scripts grow 9-10 steps) passes ``area_cutoff=1e-4``
explicitly; the topology there is the same at both.

``iterate_manifold`` maps exactly one step and takes no ``num_iterations``:
repeated iteration is the caller's loop (``grow_n_times``), settled (an old
TODO asked for the parameter; it is not wanted).
"""

from __future__ import annotations

from contextlib import contextmanager
from typing import Iterator, Literal, Optional
import logging

import numpy as np
from numpy.typing import NDArray

from .Point import Point
from .PointStore import LEFT_OFF, Chain, PointStore
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

#: Library default for :attr:`ManifoldMachine.area_cutoff` (see the module
#: Dev Notes).
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
      stable one), put the images on a chain, and merge it with the part of
      the curve that already existed.
    * :meth:`refine_manifold` -- insert points wherever the curvature area
      between two neighbours exceeds :attr:`area_cutoff`, breadth-first and one
      batched map call per layer.
    * :meth:`new_grow_manifold` -- drive the two above around a periodic orbit.

    Invariants every method preserves:

    * Every chain stays ordered by canonical distance (non-decreasing -- see
      the module docstring on ties at a high-stretch fold).
    * Inserted points keep their iterate links, so a grown curve can always be
      walked back to the seed segment.
    * A refined point is inserted at its position between its neighbours; a
      chain is never re-sorted.

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
        chain = view._slot_chain()
        view.root = None if chain is None else int(chain.ids[0])
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

        Growth happens by side effect: :meth:`iterate_manifold` merges the images
        into the next orbit point's chain, so the manifold it returns is
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
                if view._root is None:
                    continue  # branch not initialized

                # the whole chain shares the stretch parameter of its first piece
                stretch_param = view._span()[1].stretch_param
                view.stretch_param = stretch_param

                for _ in range(fixed_point.k_value):

                    self.iterate_manifold(view)

                    key = fixed_point.advance_key(key, step)
                    view = self._branch_view(fixed_point, key, stretch_param)
                    assert view._root is not None, (
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
    def _manifold_from_rows(
        store: PointStore,
        rows: NDArray[np.int64],
        stability: Stability,
        stretch_param: float,
        fixed_point: FixedPoint,
    ) -> BaseManifold:
        """A key-less manifold over a new chain of fresh rows, in cdist order."""
        store.new_chain(rows, stretch_param)
        return BaseManifold(
            int(rows[0]), stability, stretch_param, fixed_point, tail=int(rows[-1]),
            manifold_key=None, store=store,
        )

    @staticmethod
    def _manifold_from_points(
        points: list[Point], stability: Stability, stretch_param: float, fixed_point: FixedPoint
    ) -> BaseManifold:
        """:meth:`_manifold_from_rows` over hand-built points (an adapter for the
        tests, which N7b deletes)."""
        return ManifoldMachine._manifold_from_rows(
            points[0].store, [point.row for point in points], stability, stretch_param, fixed_point
        )

    def _iterate_manifold(self, manifold: BaseManifold) -> BaseManifold:
        """The body of :meth:`iterate_manifold`, inside the caller's budget scope."""
        viewer = ManifoldView(manifold, self.system)

        non_iterated_coords = manifold.get_non_iterated_point_array()
        non_iterated_cdists = manifold.get_non_iterated_cdist_array()
        _head, store, sources = manifold._select("non_iterated")
        image_head, image_store, image_rows = manifold._select("iterated")

        if len(non_iterated_coords):

            self._note_coords(non_iterated_coords)
            self._charge(len(non_iterated_coords), "map", manifold)
            iterated_points = viewer.map_fwd_batch(non_iterated_coords)
            self._note_coords(iterated_points)
            distances = np.ravel(manifold.stretch_param * non_iterated_cdists)

            # The images go into the store of the curve (a chain's rows share one
            # store), in one append; each pair is linked both ways (the image on
            # an unstable curve, the preimage on a stable one).
            rows = store.append(iterated_points, distances)
            unstable = manifold.stability == "unstable"
            if ((store.image if unstable else store.preimage)[sources] >= 0).any():
                raise ValueError(f"{'Iterate' if unstable else 'Pre-iterate'} 1 already exists")
            if unstable:
                store.link(sources, rows)
            else:
                store.link(rows, sources)

            # The image of a piece lives on the branch one map step forward, which
            # only the caller can name (see iterate_bridge); the machine leaves the
            # key unset rather than copying the un-advanced one.
            new_iterated_points = self._manifold_from_rows(
                store, rows, manifold.stability, manifold.stretch_param, manifold.fixed_point
            )

        # The images that already exist may lie on several curves: merge
        # each contiguous chain of them in turn.
        chains = self._image_chains(manifold, image_head, image_store, image_rows)

        if len(non_iterated_coords):
            mapped_manifold = new_iterated_points
            for chain in chains:
                mapped_manifold = self.merge_manifolds(chain, mapped_manifold)
            self.refine_manifold(mapped_manifold)

            # cdist is the manifold's ordering key, so the chain must be
            # cdist-sorted. It is only *non-decreasing*, not strictly increasing:
            # where a high-stretch fold has collapsed adjacent cdists to within a
            # float ULP the refiner bridges the gap with equal-cdist points (a tie
            # is harmless -- such points are inserted at a position, never re-sorted).
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
            root=image._root,
            store=image.store,
            stability=manifold.stability,
            stretch_param=manifold.stretch_param,
            fixed_point=manifold.fixed_point,
            tail=image._tail,
            branch_index=manifold.branch_index,
            manifold_key=None,
        )
        if isinstance(manifold, Bridge):
            return Bridge(name=manifold.name, **piece)
        return BaseManifold(**piece)

    def _image_chains(
        self,
        manifold: BaseManifold,
        head: Optional[BranchPoint],
        store: Optional[PointStore],
        rows: NDArray[np.int64],
    ) -> list[BaseManifold]:
        """
        The existing images of a curve's points, as runs of their chains.

        A bridge's images are usually one run of one chain, but a bridge whose
        straddle nodes it shares with a neighbour that was iterated first has its
        images on TWO chains (the grown manifold's tip, then the start of the
        neighbour's image). The run starting at an image takes every later image
        on the same chain at or after its position whose cdist is within the last
        image's (the extent a walk from it along its curve reached); the next run
        starts at the first image after the last one taken. An image on no chain
        inside a run (a node a merge dropped as a duplicate cdist keeps its
        preimage's link) is skipped, as one merged chain always skipped it; one
        that starts a run (an initializer's fictitious point) is a chain of one.
        A periodic point (the image of an anchor bridge's root, always first)
        starts a run at its slot on the curve's branch.

        Args:
            manifold: The curve being iterated (its stability, stretch and branch
                go to the runs).
            head: The periodic point heading the images, or None.
            store: The store of the images' rows.
            rows: The images' rows after the head, in the curve's order.

        Returns:
            One manifold per run, in cdist order; empty when there are no images.

        Raises:
            AssertionError: A run would start at a point a merge left off its chain.
        """
        def run(first: int | BranchPoint, last: int) -> BaseManifold:
            return BaseManifold(
                first,
                manifold.stability,
                manifold.stretch_param,
                manifold.fixed_point,
                tail=last,
                branch_index=manifold.branch_index,
                manifold_key=None,
                store=store,
            )

        if head is None and not len(rows):
            return []
        images = ([] if head is None else [head]) + rows.tolist()
        chains, positions, cdists = store._chain[rows], store._position[rows], store.cdist[rows]
        if head is not None:
            # The image of an anchor bridge's root is the next periodic point: its
            # run starts at the slot of the bridge's branch, position -1 of the
            # chain hanging there.
            slots = head.forward_branches if manifold.stability == "unstable" else head.backward_branches
            slot_chain = slots[manifold.branch_index]
            rows = np.concatenate(([-1], rows))
            chains = np.concatenate(([slot_chain.index], chains))
            positions = np.concatenate(([-1], positions))
            cdists = np.concatenate(([head.cdist], cdists))
        bound = cdists[-1]
        runs, i = [], 0
        while i < len(images):
            if chains[i] < 0:
                assert positions[i] != LEFT_OFF, "a run of images starts at a point a merge left off"
                store.new_chain(rows[i : i + 1], manifold.stretch_param)
                last = i
            else:
                reached = (
                    (chains[i:] == chains[i]) & (positions[i:] >= positions[i]) & (cdists[i:] <= bound)
                )
                last = i + int(np.flatnonzero(reached)[-1])
            runs.append(run(images[i], images[last]))
            i = last + 1
        return runs

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

    def merge_manifolds(
        self, manifold_1: BaseManifold, manifold_2: BaseManifold
    ) -> BaseManifold:
        """
        Merges two runs of points by canonical distance into one chain, in place.

        ``manifold_1`` is a run of existing points (from :meth:`_image_chains`),
        ``manifold_2`` the new images (or an earlier merge's output). The merged
        run starts at the lower root (``manifold_1``'s on a tie) and ends at the
        tail of the run that outlasts the other. Where the two runs lie on
        different chains and share no cdist (the usual merge) the order is one ``searchsorted``; otherwise the
        zipper below decides it point by point: tied ROOTS are both kept, the
        existing first; a point on both runs is kept once; of two other points
        on one cdist the existing one is kept unless only the new one has an
        iterate link, and then the existing one leaves its chain.

        Parameters:
            manifold_1: The existing run.
            manifold_2: The run merged into it.

        Returns:
            Whichever manifold starts lower, its tail set to the merged run's end.

        Raises:
            AssertionError: The merged run cannot be one chain (no case reaches
                these): the run placed second starts after the beginning of its
                chain, the run that ran out first has points past its tail, or
                both runs end on a tie.
        """
        # A periodic point heads only an existing run (an anchor bridge's image);
        # at cdist 0 it comes first, and the zipper starts after it.
        root_1, chain_1, lo_1, hi_1 = manifold_1._span()
        root_2, chain_2, lo_2, hi_2 = manifold_2._span()
        assert root_2 is None, "a periodic point heads only the existing run"
        store = chain_1.store
        ids_1, ids_2 = chain_1.ids[lo_1:hi_1], chain_2.ids[lo_2:hi_2]
        c1, c2 = store.cdist[ids_1], store.cdist[ids_2]
        ranges = f"cdists [{c1[0]!r}, {c1[-1]!r}] and [{c2[0]!r}, {c2[-1]!r}]"
        first_is_1 = root_1 is not None or c1[0] <= c2[0]

        below = np.minimum(np.searchsorted(c1, c2), len(c1) - 1)
        if chain_1 is not chain_2 and not (c1[below] == c2).any():
            merged = np.insert(ids_1, np.searchsorted(c1, c2, side="right"), ids_2)
            rest_is_1 = c1[-1] > c2[-1]
        else:
            if root_1 is not None:
                order, i, j = [], 0, 0
            else:
                order = [ids_1[0] if first_is_1 else ids_2[0]]
                i, j = (1, 0) if first_is_1 else (0, 1)
            while i < len(ids_1) and j < len(ids_2):
                if c1[i] < c2[j]:
                    order.append(ids_1[i])
                    i += 1
                elif c2[j] < c1[i]:
                    order.append(ids_2[j])
                    j += 1
                elif ids_1[i] == ids_2[j]:
                    order.append(ids_1[i])
                    i, j = i + 1, j + 1
                else:
                    # Two different points collided onto the same cdist. cdist is
                    # the only ordering key, so once it saturates (one ULP at large
                    # magnitude exceeds the true separation near a fold) it can no
                    # longer order the pair -- keeping both and ordering by the tied
                    # key scrambles them into a spike. Keep the established
                    # (iterate-linked) point and skip the duplicate.
                    pair = [ids_1[i], ids_2[j]]
                    logger.debug(
                        "Duplicate cdist collision at cdist=%s; coords %s vs %s",
                        float(c1[i]),
                        store._xy_read[pair[0]],
                        store._xy_read[pair[1]],
                    )
                    # An image or a preimage row, a phantom's included.
                    linked = (store._image[pair] >= 0) | (store._preimage[pair] >= 0)
                    keep_1 = linked[0] or not linked[1]
                    order.append(ids_1[i] if keep_1 else ids_2[j])
                    i, j = i + 1, j + 1
            rest_is_1 = i < len(ids_1)
            assert chain_1 is chain_2 or rest_is_1 or j < len(ids_2), f"both runs end on a tie: {ranges}"
            merged = np.concatenate((order, ids_1[i:] if rest_is_1 else ids_2[j:])).astype(np.int64)

        assert root_1 is None or merged[0] == ids_1[0], (
            f"a new point before a periodic point's first point: {ranges}"
        )
        if chain_1 is chain_2:
            lo = min(lo_1, lo_2)
            assert np.array_equal(chain_1.ids[lo : lo + len(merged)], merged), (
                f"two runs of one chain that do not join: {ranges}"
            )
        else:
            first, start = (chain_1, lo_1) if first_is_1 else (chain_2, lo_2)
            second, second_lo = (chain_2, lo_2) if first_is_1 else (chain_1, lo_1)
            assert second_lo == 0 and second.anchor is None, (
                f"the run placed second starts after the beginning of its chain: {ranges}"
            )
            done, done_hi = (chain_2, hi_2) if rest_is_1 else (chain_1, hi_1)
            assert done_hi == len(done), f"the run that ran out goes on past its tail: {ranges}"
            rest, rest_hi = (chain_1, hi_1) if rest_is_1 else (chain_2, hi_2)
            first.absorb(second, np.concatenate((first.ids[:start], merged, rest.ids[rest_hi:])), start)

        output = manifold_1 if first_is_1 else manifold_2
        output.tail = manifold_1._tail if rest_is_1 else manifold_2._tail
        return output

    def refine_manifold(
        self,
        manifold: BaseManifold,
        branch_index: Optional[int] = None,
        final_node: Optional[Point] = None,
    ) -> None:
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
            branch_index: Unused since a slot holds its chain (N7a): a point
                inserted right after a BranchPoint root needs no slot update.
            final_node: Unused; kept for backwards-compatible signature.

        Note:
            The per-pair decision depends only on local geometry and preiterates,
            so the traversal order does not matter: breadth-first reaches the
            same ``area < cutoff`` stopping condition as any other order.
        """
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("Refining manifold: %r", manifold.get_point_array())

        num_initial_points = (
            len(manifold.get_point_array())
            if logger.isEnabledFor(logging.INFO)
            else None
        )

        viewer = ManifoldView(manifold, self.system)
        if branch_index is None:
            branch_index = manifold.branch_index

        # Layer 0: every consecutive pair from root to tail, by the position of
        # its first point (-1: the BranchPoint root, the chain's anchor).
        head, chain, lo, hi = manifold._span()
        if chain is None:
            return
        layer = np.arange(lo - (head is not None), hi - 1)
        while len(layer):
            layer = self._refine_layer(chain, layer, manifold, viewer, branch_index)

        if num_initial_points is not None:
            num_final_points = len(manifold.get_point_array())
            logger.debug(
                "%d points added during refinement (%d total)",
                num_final_points - num_initial_points,
                num_final_points,
            )

    def _refine_layer(
        self,
        chain: Chain,
        layer: NDArray[np.int64],
        manifold: BaseManifold,
        viewer: ManifoldView,
        branch_index: Optional[int],
    ) -> NDArray[np.int64]:
        """
        Refine one breadth-first layer of pairs and return the next layer.

        Computes curvature areas for the whole layer at once, inserts a new
        point into each pair that needs refining (mapping all midpoints forward
        in a single batched call, inserting all of them with one chain insert),
        and returns the child pairs of the refined pairs as the next layer.

        A pair is its first point's position: the pair is that point and the
        next one on the chain. Position -1 is the chain's anchor, which is also
        the left neighbour of position 0; the right neighbour of the last
        position is missing (NaN), and so is the left of the anchor.

        Args:
            chain: the chain the pairs lie on.
            layer: ascending positions of the pairs' first points.
            manifold: the manifold being refined.
            viewer: view binding the manifold to the dynamical system.
            branch_index: Unused since a slot holds its chain (N7a): a point
                inserted right after the anchor needs no slot update.

        Returns:
            The next layer, the positions of the child pairs (empty when
            refinement is done).
        """
        stability = manifold.stability
        store, ids, anchor = chain.store, chain.ids, chain.anchor
        xy, cdist = store.xy, store.cdist
        preiterate = store.preimage if stability == "unstable" else store.image
        nan_row = np.array([np.nan, np.nan])

        def coords(positions: NDArray[np.int64]) -> NDArray[np.float64]:
            """The points at chain positions: the anchor at -1, NaN off the chain."""
            inside = (positions >= 0) & (positions < len(ids))
            out = np.where(inside[:, None], xy[ids[np.clip(positions, 0, len(ids) - 1)]], nan_row)
            if anchor is not None:
                out[positions == -1] = anchor.get_point()
            return out

        p0_xy, p1_xy = coords(layer), coords(layer + 1)
        rows_0 = ids[np.maximum(layer, 0)]
        rows_1 = ids[layer + 1]
        pre_rows_0, pre_rows_1 = preiterate[rows_0], preiterate[rows_1]
        pre0 = xy[np.maximum(pre_rows_0, 0)]
        pre1 = xy[np.maximum(pre_rows_1, 0)]
        c0, c1 = cdist[rows_0], cdist[rows_1]
        at_anchor = layer == -1
        if at_anchor.any():
            anchor_pre = self._get_preiterate(anchor, stability)
            pre_rows_0 = np.where(at_anchor, -1 if anchor_pre is None else 0, pre_rows_0)
            if anchor_pre is not None:
                pre0[at_anchor] = anchor_pre
            c0 = np.where(at_anchor, float(anchor.cdist), c0)

        # Pairs that are degenerate (near-coincident points) or missing the data
        # needed to refine (a missing preiterate on a bridge boundary point) are
        # dropped.
        apart = np.hypot(p1_xy[:, 0] - p0_xy[:, 0], p1_xy[:, 1] - p0_xy[:, 1]) >= 1e-8
        has_pre = (pre_rows_0 >= 0) & (pre_rows_1 >= 0)
        if logger.isEnabledFor(logging.DEBUG):
            for k in np.flatnonzero(apart & ~has_pre):
                logger.debug(
                    "Skipping refinement for pair at cdist %.3g-%.3g: "
                    "preiterate missing on one or both points.",
                    c0[k],
                    c1[k],
                )
        keep = np.flatnonzero(apart & has_pre)
        if not keep.size:
            return np.empty(0, dtype=np.int64)
        layer = layer[keep]

        areas = self._curvature_area_batch(
            p0_xy[keep], p1_xy[keep], coords(layer - 1), coords(layer + 2)
        )

        # Refine every pair whose curvature area exceeds the cutoff. The new point
        # takes the cdist midpoint of its neighbours; at a high-stretch fold the
        # curve has expanded so much that adjacent points are far apart in space
        # while their cdists have collapsed to within a float ULP, so that midpoint
        # may round exactly onto an endpoint. That tie is harmless here -- the point
        # is inserted at its position between its neighbours, not re-sorted by cdist
        # -- and bridging the gap is exactly what keeps the curve smooth instead of
        # leaving a coarse spike. The ``< 1e-8`` spacing guard above guarantees the
        # subdivision still terminates once the two points effectively coincide.
        mid_cd = 0.5 * (c0[keep] + c1[keep])
        pre_mid = 0.5 * (pre0[keep] + pre1[keep])

        flagged = np.flatnonzero(areas >= self.area_cutoff)
        if flagged.size == 0:
            return np.empty(0, dtype=np.int64)

        self._charge(flagged.size, "refine", manifold)
        new_coords = viewer.map_fwd_batch(pre_mid[flagged])
        self._note_coords(new_coords)

        # The new points, then their preiterates (preimages on an unstable curve,
        # images on a stable one): phantom rows (no cdist) the child pairs refine
        # from next layer, named by a one-way link. A new first point (position 0)
        # needs no slot update: the slot holds the chain.
        rows = store.append(new_coords, mid_cd[flagged])
        phantoms = store.append(pre_mid[flagged], None)
        store.link_phantoms(rows, phantoms, stability)

        refined = layer[flagged]
        chain.insert(refined + 1, rows)

        # After the insert the k-th refined pair's first point sits k places
        # further on, and its new point right after it.
        shifted = refined + np.arange(len(refined))
        return np.column_stack((shifted, shifted + 1)).ravel()

    @staticmethod
    def _get_preiterate(point: BranchPoint, stability: Stability) -> Optional[np.ndarray]:
        """
        The coordinates of a periodic point's preiterate (its preimage on an
        unstable curve, its image on a stable one), or None if it has none.

        Parameters:
            point (BranchPoint): The chain's anchor (a row's preiterate is read
                off the store's link columns).
            stability (str): stability of the manifold
        """
        pre = point.prev_iterate if stability == "unstable" else point.next_iterate
        return None if pre is None else pre.get_point()

    @staticmethod
    def _curvature_area_batch(
        p0: np.ndarray,
        p1: np.ndarray,
        left: np.ndarray,
        right: np.ndarray,
    ) -> np.ndarray:
        """The curvature area of each segment (:func:`.geometry.curvature_area_batch`)."""
        return curvature_area_batch(p0, p1, left, right)

