from __future__ import annotations

from typing import Literal, Optional, Tuple
from collections import deque

import numpy as np
from .BranchPoint import BranchPoint

ManifoldKey = Tuple["FixedPoint", Literal["unstable", "stable"], int, int]

# A planar saddle has exactly two eigendirections per manifold, so a BranchPoint
# always carries two branch slots. This is the number of DIRECTIONS available,
# not the number of pieces the dynamical chain visits (FixedPoint.num_branches).
_EIGENDIRECTIONS = 2

# Relative slack on |lambda_u * lambda_s| == 1 before a map is treated as
# something other than area preserving (see FixedPoint.per_step_beta).
_DET_TOL = 1e-9

"""
Dev Notes:

Ensure that the shape of the coordinates array is consisntent everywhere
in the package. In this class it is (2, 1). I know in BasePoint it is (2,).

Check how the unstable eigenvectors is implemented when it comes to points 
with inversion. We may or may not need to store two vectors, the method
we use elsewhere I think just takes the initial vector and maps it around
k times so we don't need to be able to walk in both eigendirections likely.

Specify why we store the partial_jacobians. What are they used for?
"""


class FixedPoint:
    """
    Implements the data structure to store all information about a fixed point.

    Note:
        Does not include the methods to compute nor construct fixed points.
        Those are located in FixedPointSolver.py.

    Attributes:

        period (int): Period of the fixed point.
        branch_points (List[BranchPoint]): List of the individual BranchPoints that
            make up the fixed point.
        coordinates (List[np.ndarray]): Array storing the coordinates of the
            fixed point. The shape is (2, 1).

        accuracy (float): Computed as the norm of the difference between the
            fixed point and it's iterate. {latex}
        k_value (int): Number of iterations it takes to return to the segment
            of manifold you started from.

        unstable_eigenvectors (List[np.ndarray]): List of arrays storing the
            eigenvectors corresponding to the direction of the unstable manifold.
        unstable_eigenvalues (List[float]): List of arrays storing the eigenvalues
            corresponding to the unstable manifold. If an eigenvalue is < 0 the
            fixed point has inversion.
        stable_eigenvectors (List[np.ndarray]): List of arrays storing the eigenvectors
            corresponding to the direction fo the stable manifold.
        stable_eigenvalues (List[float]): List of arrays storing the eigenvalues
            corresponding to the stable manifold. If an eigenvalue is < 0 the
            fixed point has inversion.

        jacobians (List[np.ndarray]): List of the full cycle Jacobians at each
            iterate of the fixed point. These Jacobians are used to compute the
            eigenstuffs.
        partial_jacobians (List[np.ndarray]): List of the single step Jacobians
            at each iterate of the fixed point. These Jacobians are used to
    """

    def __init__(self, period: int) -> None:
        """
        Allocates the data structures which store fixed point information.

        Note:
            Does not include the methods to compute nor construct fixed points.
            Those are located in FixedPointSolver.py.

        Args:
            period (int): Period of the fixed point.

        Note:
            The branch count is NOT an input: it is derived from the eigenvalues
            by :attr:`num_branches` once :meth:`set_k_value` has run. Each
            BranchPoint is allocated with a slot per EIGENDIRECTION (always two),
            which is what the previous callers passed anyway; how many of those
            slots the k-chain visits is the derived :attr:`num_branches`.
        """

        self.period = period

        self.branch_points = [
            BranchPoint(_EIGENDIRECTIONS, (0.0, 0.0)) for _ in range(period)
        ]
        self.coordinates = [np.empty((2, 1)) for _ in range(period)]

        self.unstable_eigenvectors = [np.empty((2, 1)) for _ in range(period)]
        self.unstable_eigenvalues = [0.0] * period

        self.stable_eigenvectors = [np.empty((2, 1)) for _ in range(period)]
        self.stable_eigenvalues = [0.0] * period

        self.accuracy = 0.0

        self.jacobians = [np.empty((2, 2)) for _ in range(period)]
        self.partial_jacobians = [np.empty((2, 2)) for _ in range(period)]

        # Filled in by set_k_value() once the eigenvalues are known.
        self.k_value: Optional[int] = None

        # Memoised branch bookkeeping. Both are functions of (k_value, period)
        # alone, and k_value only ever changes in set_k_value(), which clears
        # them. Keyed by stability.
        self._branch_cycle_cache: dict[str, list[ManifoldKey]] = {}
        self._branch_position_cache: dict[str, tuple[dict[ManifoldKey, int], int]] = {}

    @property
    def num_branches(self) -> int:
        """
        Number of manifold branches per stability the dynamical chain visits.

        Two if the fixed point has inversion (the negative eigenvalue swaps the
        two eigendirections every ``period`` map steps, so one chain of
        ``k_value`` pieces covers both), one if it does not.

        This is DERIVED, never supplied: it is ``k_value // period``. A caller
        may still initialize a non-inversion saddle on both eigendirections --
        those are two independent chains of one branch each, not this count.

        Raises:
            ValueError: If called before set_k_value().
        """
        return 2 if self.check_inversion() else 1

    def check_inversion(self) -> bool:
        """
        Returns True if the Fixed Point has inversion
        False if it does not.

        Raises:
            ValueError: If called before set_k_value().
        """

        if self.k_value is None:
            raise ValueError(
                "k_value has not been set yet; call set_k_value() after the "
                "eigenvalues are computed."
            )

        return self.period != self.k_value

    def set_k_value(self) -> None:
        """
        Sets the value 'k' which describes how many iterations it takes a nearby
        point to get back to that neighborhoods. If the fixed point has inversion
        this is double the period.

        Raises:
            ValueError: If the unstable and stable eigenvalues have DISAGREEING
                signs. That is an orientation-reversing map (det J < 0): exactly
                one of the two manifolds inverts, so the unstable and stable
                sides have different inversion status and a single ``k_value``
                cannot describe both. Eigenvalues still at their ``0.0``
                placeholder carry no sign and are skipped.

        Note:
            Only the sign matters here; the magnitudes are used by
            :meth:`per_step_beta`.
        """

        for orbit_index, (lambda_u, lambda_s) in enumerate(
            zip(self.unstable_eigenvalues, self.stable_eigenvalues)
        ):
            u = float(np.asarray(lambda_u).ravel()[0])
            s = float(np.asarray(lambda_s).ravel()[0])
            if u == 0.0 or s == 0.0:
                continue  # eigenvalue not computed yet -- no sign to compare
            if (u < 0) != (s < 0):
                raise ValueError(
                    f"orbit point {orbit_index} has eigenvalues of disagreeing "
                    f"sign (unstable {u!r}, stable {s!r}): the map is "
                    "orientation reversing there, so only one of the two "
                    "manifolds inverts and a single k_value cannot model it."
                )

        multiplier = 2 if any(x < 0 for x in self.unstable_eigenvalues) else 1

        self.k_value = self.period * multiplier

        # k_value is the only input to the branch chain, so this is the one
        # place the memos below can go stale.
        self._branch_cycle_cache.clear()
        self._branch_position_cache.clear()

    def per_step_beta(self, stability: Literal["unstable", "stable"]) -> float:
        """
        The canonical-distance factor of ONE application of the map.

        A point at canonical distance ``c`` on a branch images to ``c *
        per_step_beta(stability)`` on the branch one map step along the chain
        (:meth:`advance_key`). The stored eigenvalues are those of the
        FULL-CYCLE Jacobian DM^period, so the factor is the ``period``-th root
        of the eigenvalue -- NOT the ``k_value``-th root: the unstable factor is
        greater than 1 (expansion), the stable factor less than 1 (contraction).

        The distinction only bites on an inversion point. There M^period already
        carries a point all the way onto the OPPOSITE branch, growing its
        canonical distance by the full ``|lambda_u|`` in ``period`` steps; a
        branch RETURN is ``k_value = 2 * period`` steps and multiplies the
        canonical distance by ``lambda_u ** 2``. In general
        ``per_step_beta(s) ** k_value == |lambda| ** num_branches``.

        Args:
            stability: Which manifold's canonical distance is being scaled.

        Returns:
            float: The per-map-step factor.

        Raises:
            ValueError: If ``stability`` is not "unstable"/"stable", or if
                ``set_k_value()`` has not been called.

        Note:
            The stable factor is taken as ``1 / per_step_beta("unstable")``
            whenever ``lambda_u * lambda_s == 1`` within 1e-9 -- i.e. for the
            area-preserving (det J = +1) maps this package is written for. That
            is exactly what every call site did before this method existed
            (``stable_cdist / beta_unstable``), so results are unchanged, and it
            is the more accurate of the two on such a map because the two
            eigenvalues are then one number, not two independent estimates.
            Otherwise the stable eigenvalue's own magnitude is used. A stable
            eigenvalue still at its ``0.0`` placeholder (eigendata incomplete)
            also falls back to the reciprocal.

        Note:
            ``k_value`` is still required (hence the guard below): it is what
            says whether the point inverts, and a caller counting branch returns
            must go through it.
        """
        if stability not in ("unstable", "stable"):
            raise ValueError(
                f"stability must be 'unstable' or 'stable', got {stability!r}"
            )
        if self.k_value is None:
            raise ValueError(
                "k_value has not been set yet; call set_k_value() after the "
                "eigenvalues are computed."
            )

        lambda_u = float(np.abs(np.asarray(self.unstable_eigenvalues[0]).ravel()[0]))
        beta_u = lambda_u ** (1.0 / self.period)
        if stability == "unstable":
            return beta_u

        lambda_s = float(np.abs(np.asarray(self.stable_eigenvalues[0]).ravel()[0]))
        if lambda_s == 0.0 or abs(lambda_u * lambda_s - 1.0) <= _DET_TOL:
            return 1.0 / beta_u
        return lambda_s ** (1.0 / self.period)

    def per_step_factor(self, k_step_factor: float) -> float:
        """
        Convert a MEASURED ``k_value``-step growth ratio into a per-map-step one.

        The initializer measures how far a seed point moves over one full return
        to its own branch (``ManifoldInitializer.get_first_point_back`` applies
        the map ``k_value`` times) rather than reading the eigenvalue, because
        that ratio is the stretch actually applied to the points it creates.
        This is the one place that conversion happens, so no call site takes a
        ``k_value``-th root by hand.

        Args:
            k_step_factor: The growth ratio over ``k_value`` map steps.

        Returns:
            float: The equivalent single-map-step factor.

        Raises:
            ValueError: If ``set_k_value()`` has not been called.

        Note:
            The exponent here is ``1 / k_value`` while :meth:`per_step_beta`
            uses ``1 / period``, and the two nevertheless agree: they are roots
            of DIFFERENT quantities. This one takes the ratio measured over
            ``k_value`` actual map applications, which on an inversion point is
            ``lambda ** 2`` (a branch return), so its ``k_value``-th root is
            ``|lambda| ** (1 / period)`` -- exactly what per_step_beta computes
            from the full-cycle eigenvalue. They differ only by the
            linearization error of the seed step.
        """
        if self.k_value is None:
            raise ValueError(
                "k_value has not been set yet; call set_k_value() after the "
                "eigenvalues are computed."
            )
        return k_step_factor ** (1 / self.k_value)

    def branch_cycle(
        self, stability: Literal["unstable", "stable"]
    ) -> list[ManifoldKey]:
        """
        The ``k_value`` manifold pieces of one stability, in map-step order.

        Element ``j`` maps to element ``j + 1`` under one application of M, and
        the last maps back to the first: this is exactly the orbit of
        ``(self, stability, 0, 0)`` under :meth:`advance_key` with ``n = +1``.
        Without inversion it is the ``period`` orbit points on branch 0; with
        inversion it is all ``2 * period`` pieces, the branch flipping on each
        wrap.

        Args:
            stability: Which manifold's branches to list.

        Returns:
            list[ManifoldKey]: The ``k_value`` keys in forward-map order. A
            fresh list each call; the chain itself is memoised and only
            recomputed after :meth:`set_k_value`.

        Raises:
            ValueError: If ``stability`` is invalid or ``k_value`` is unset.

        Note:
            The order is FORWARD for both stabilities: M maps W^s(z_i) onto
            W^s(z_{i+1}) exactly as it maps W^u(z_i) onto W^u(z_{i+1}), so the
            stable branches cycle in the same orbit order the unstable ones do
            (the stable manifold is *grown* by M^-1, which is a different
            question). The topological algorithms index this list and take
            differences of positions modulo ``k_value`` as counts of forward map
            steps, so reversing it for the stable side would silently negate
            every one of those counts.
        """
        if stability not in ("unstable", "stable"):
            raise ValueError(
                f"stability must be 'unstable' or 'stable', got {stability!r}"
            )
        if self.k_value is None:
            raise ValueError(
                "k_value has not been set yet; call set_k_value() after the "
                "eigenvalues are computed."
            )

        cached = self._branch_cycle_cache.get(stability)
        if cached is None:
            key: ManifoldKey = (self, stability, 0, 0)
            cached = [key]
            for _ in range(self.k_value - 1):
                key = self.advance_key(key, 1)
                cached.append(key)
            self._branch_cycle_cache[stability] = cached
        return list(cached)

    def branch_position_map(
        self, stability: Literal["unstable", "stable"]
    ) -> tuple[dict[ManifoldKey, int], int]:
        """
        Where each branch of one stability sits in the forward map-step chain.

        The lookup form of :meth:`branch_cycle`: ``position[key]`` is the number
        of forward map steps from ``(self, stability, 0, 0)`` to ``key``, and
        differences of positions modulo the returned length are counts of map
        steps -- which is how the topological algorithms compare two branches.

        Args:
            stability: Which manifold's branches to map.

        Returns:
            ``(position-of-each-branch, k_value)``.

        Raises:
            ValueError: If ``stability`` is invalid or ``k_value`` is unset.

        Note:
            The mapping is the memo itself, shared by every caller and only
            rebuilt after :meth:`set_k_value`. Treat it as read-only.
        """
        cached = self._branch_position_cache.get(stability)
        if cached is None:
            cycle = self.branch_cycle(stability)
            cached = ({key: i for i, key in enumerate(cycle)}, len(cycle))
            self._branch_position_cache[stability] = cached
        return cached

    def advance_key(self, key: ManifoldKey, n: int = 1) -> ManifoldKey:
        """
        Advance a manifold key by ``n`` applications of the map.

        This is the single source of truth for "which manifold piece does one map
        step send this one to". One application of M sends the manifold piece
        attached to orbit point ``i`` onto the piece at orbit point ``(i + 1) mod
        period``, for BOTH stabilities: M maps W^s(x_i) onto W^s(x_{i+1}) exactly
        as it maps W^u(x_i) onto W^u(x_{i+1}). ``stability`` is therefore carried
        through untouched -- it is the caller who decides whether it wants the
        forward step (``n = +1``, how the unstable fundamental segment is grown)
        or the backward one (``n = -1``, how the stable one is).

        The branch index labels the eigenvector direction at orbit point 0; the
        pieces at the other orbit points inherit it by propagation. Without
        inversion the branch never changes. With inversion (negative eigenvalue,
        ``k_value == 2 * period``) M^period returns to orbit point 0 on the
        OPPOSITE branch, so the branch flips exactly when the orbit index wraps
        (period-1 -> 0 going forward, 0 -> period-1 going backward) and two full
        orbits -- ``k_value`` steps -- return to the start.

        The map is a group action on the ``k_value`` pieces: ``advance_key`` is
        additive in ``n`` and has period ``k_value``.

        Args:
            key: A manifold key ``(fixed_point, stability, orbit_index,
                branch_index)``. Its fixed point must be this one.
            n: Number of map steps. Negative means applications of M^-1.
                Defaults to 1.

        Returns:
            ManifoldKey: The key ``n`` map steps along the chain.

        Raises:
            ValueError: If the key belongs to another fixed point, if
                ``set_k_value()`` has not been called, or if the branch index is
                neither of the two eigendirections. Note that branch 1 of a
                point WITHOUT inversion is a valid key -- the two
                eigendirections of a simple saddle are two independent chains,
                each advancing its orbit index with the branch carried along --
                even though such a point's :attr:`num_branches` is 1.

        Note:
            The ordering this induces on an inversion point is
            ``(0, 0), (1, 0), ..., (p-1, 0), (0, 1), ..., (p-1, 1)``: the piece
            after ``(p-1, 0)`` is ``(0, 1)``. The pre-2026-09 initializer instead
            walked ``for orbit_index: for branch_index:``, visiting ``(0, 1)``
            before ``(1, 0)``, which is not the physical chain.
        """

        fp, stability, orbit_index, branch_index = key

        if fp is not self:
            raise ValueError(
                "advance_key was given a key belonging to a different fixed point."
            )

        if self.k_value is None:
            raise ValueError(
                "k_value has not been set yet; call set_k_value() after the "
                "eigenvalues are computed."
            )

        if branch_index not in range(_EIGENDIRECTIONS):
            raise ValueError(
                f"branch_index {branch_index} is outside this fixed point's branch "
                f"range (a planar saddle has {_EIGENDIRECTIONS} eigendirections)."
            )

        if self.check_inversion():
            assert self.k_value == 2 * self.period, (
                "an inversion point's chain must run over both branches of the "
                f"orbit; got k_value={self.k_value} for period {self.period}"
            )
            # The k_value pieces are laid out as orbit_index + period * branch_index,
            # so one map step is +1 on that single cyclic index and the branch flip
            # falls out of the wrap.
            position = (orbit_index + self.period * branch_index + n) % self.k_value
            return (fp, stability, position % self.period, position // self.period)

        return (fp, stability, (orbit_index + n) % self.period, branch_index)

    def get_iterable_array(
        self, stability: Literal["unstable", "stable"], shift: int = 0
    ) -> list[int]:
        """
        Construct a list of orbit indices based off the stability.
        Starting at the zeroth index this is the order a manifold will move
        between the indices of the fixed points

        Ex:
            unstable [0, 1, 2, 3]
            stable [0, 3, 2, 1]

        Parameters:
            stability: The stability of the manifold.
            shift: Amount to cyclically shift the resulting array
                    "takes the last element and moves it to the front".
                    This is equivalent to starting at another point in the cycle.
        """

        orbit_indices = [i for i in range(self.period)]
        if stability == "stable":
            temp = orbit_indices[1:]
            orbit_indices = [0] + temp[::-1]

        orbit_indices = deque(orbit_indices)
        orbit_indices.rotate(shift)

        return list(orbit_indices)

    def get_branch_array(self) -> list[int]:
        """
        Get an array to iterate over the branch indices.

        Returns:
            list[int]: A list of branch indices [0] or [0, 1], depending on
                whether the fixed point has inversion (see :attr:`num_branches`).
        """

        return list(range(self.num_branches))
