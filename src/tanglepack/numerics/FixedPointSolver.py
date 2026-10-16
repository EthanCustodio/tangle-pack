"""
Locating periodic orbits and building :class:`~.FixedPoint.FixedPoint` objects.

:class:`FixedPointSolver` runs the shooting solve for a period-p orbit, checks
that what it found is a saddle, and fills in the eigendata, k-value and branch
points of the resulting fixed point.

Dev Notes:

The branch count is derived from the eigenvalues by FixedPoint.num_branches (2
with inversion, 1 without), never passed to construct_fixed_point. Growing a
NON-inversion saddle on both of its eigendirections is supported, but it is a
request made per manifold at initialization time
(ManifoldInitializer.construct_kevin_way / TangleWorkbench.initialize_manifold),
not a property of the fixed point -- the two directions are two independent
one-branch chains there, not the two halves of one chain.

FixedPoint.set_k_value rejects a point whose unstable and stable eigenvalues
disagree in sign. det DM^period = (det DM)^period, so on an orientation-reversing
map (det J < 0) that product is negative for every ODD-period orbit: the guard
therefore closes the door on odd-period orbits of det < 0 maps entirely, not just
on the occasional awkward one. Even-period orbits of such a map are fine (the
product is positive). Supporting the odd case needs a per-stability inversion
flag rather than one k_value, which is a data-model change, not a solver one.
"""

from __future__ import annotations

import logging
from typing import Callable, Optional

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import fsolve
from scipy.differentiate import jacobian as jacob

from .DynamicalSystem import DynamicalSystem
from .FixedPoint import FixedPoint

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Hook signature: (orbit_index, unstable_eigenvector, stable_eigenvector) ->
# (unstable_eigenvector, stable_eigenvector), each a (2, 1) array.
OrientHook = Callable[
    [int, NDArray[np.float64], NDArray[np.float64]],
    tuple[NDArray[np.float64], NDArray[np.float64]],
]

# Relative tolerance on the imaginary part of an eigenvalue before the pair is
# considered complex rather than a real pair with round-off.
_IMAG_TOL = 1e-9

# Relative tolerance under which two orbit points are the same point. A solve that
# slid to a lower-period orbit repeats it to ~1e-14; genuine orbit points sit
# >= 1e-2 apart across the Henon family, so anything in 1e-12..1e-4 separates.
_PERIOD_RTOL = 1e-9


class FixedPointSolver:
    """
    Toolbox for computing fixed points of any period and initializing
    FixedPoint objects.

    Attributes:
        dynamical_map (func): Function for the map of the system.
        dynamical_map_inverse (func): Function for the inverse map of the system.
        jacobian_function (func, optional): Function to compute the Jacobian at
            any point in the system. Without an explicit Jacobian function finite
            difference will be used to compute the Jacobians.
        orient (OrientHook, optional): Per-orbit-index eigenvector orientation
            hook. See :meth:`__init__`.
    """

    def __init__(
        self, system: DynamicalSystem, orient: Optional[OrientHook] = None
    ) -> None:
        """
        Initializes the functions associated with the dynamical system.

        Args:
            system (DynamicalSystem): Object storing all the system maps.
            orient (OrientHook, optional): Called once per orbit index with
                ``(orbit_index, unstable_eigenvector, stable_eigenvector)`` right
                after the eigen-decomposition, returning the (possibly re-signed)
                pair. Defaults to None, meaning the eigenvectors are used exactly
                as ``np.linalg.eig`` returns them.

        Note:
            ``np.linalg.eig`` fixes each eigenvector only up to sign, so which of
            the two branches of a manifold gets grown is otherwise arbitrary. The
            intended, supported way to choose directions is
            :meth:`TangleWorkbench.orient_eigenvectors`, which re-signs the
            eigenvectors post hoc, per fixed point, against directions the caller
            names. ``orient`` exists only for the narrow case of a map where a
            fixed sign rule is known ahead of time.
        """

        self.dynamical_map = system.map
        self.dynamical_map_inverse = system.map_inv
        self.jacobian_function = system.jacobian
        self.orient = orient

    def construct_fixed_point(
        self, initial_guess: NDArray[np.float64]
    ) -> FixedPoint:
        """
        Computes the fixed point from an initial guess using a multipoint
        shooting Newton's method and initializes a FixedPoint object containing
        all the information.

        Args:
            initial_guess (np.ndarray): A (period, 2) array.
                Each row is an initial guess for one iterate.

        Returns:
            FixedPoint: The fully constructed fixed point. Its branch count is
                derived from the eigenvalues by ``FixedPoint.num_branches``.

        Raises:
            ValueError: If the orbit does not converge, if it converged to an
                orbit whose minimal period is a proper divisor of the requested
                one (checked before the saddle test, so a slid orbit is reported
                as a period problem and not as "elliptic"), if the eigenvalues
                are not a real saddle pair, or if the unstable and stable
                eigenvalues disagree in sign (``FixedPoint.set_k_value``).

        Note:
            ``period`` is the MINIMAL period of the orbit. A saddle with
            inversion has ``k_value = 2 * period``, and ``k_value`` is never a
            period to request: asking period 2 for a period-1 saddle with
            inversion repeats the point and is rejected here.
        """

        period, _ = np.shape(np.atleast_2d(initial_guess))

        fixed_point = self.compute_fixed_point(initial_guess)

        minimal = self.minimal_period(fixed_point)
        if minimal < period:
            raise ValueError(
                f"the solve converged to an orbit of minimal period {minimal}, not "
                f"the requested period {period}. Pass the {minimal}-point orbit "
                f"{fixed_point[:minimal].tolist()} instead: period is the minimal "
                f"period, and k_value (2 * period for a saddle with inversion) is "
                f"never a period to request."
            )

        difference = np.abs(self.multipoint_shoot(fixed_point) - fixed_point)
        accuracy = (np.average(np.linalg.norm(difference, axis=1))) ** (1 / 3)

        jacobians = self.compute_jacobian(fixed_point)
        partial_jacobians = self.compute_partial_jacobians(fixed_point)

        eigenvalues, eigenvectors = self.compute_eigenvectors(fixed_point, jacobians)

        point = FixedPoint(period)

        point.coordinates = fixed_point
        point.accuracy = accuracy

        for i in range(period):

            p_x, p_y = point.coordinates[i]
            point.branch_points[i].set_x(p_x)
            point.branch_points[i].set_y(p_y)

            point.branch_points[i].cdist = 0.0

            point.unstable_eigenvectors[i] = eigenvectors[i][0]
            point.stable_eigenvectors[i] = eigenvectors[i][1]

            point.unstable_eigenvalues[i] = eigenvalues[i][0]
            point.stable_eigenvalues[i] = eigenvalues[i][1]

            point.jacobians[i] = jacobians[i]
            point.partial_jacobians[i] = partial_jacobians[i]

            point.branch_points[i].next_iterate = point.branch_points[(i + 1) % period]
            point.branch_points[i].prev_iterate = point.branch_points[(i - 1) % period]

        point.set_k_value()

        return point

    @staticmethod
    def minimal_period(
        orbit: NDArray[np.float64], rtol: float = _PERIOD_RTOL
    ) -> int:
        """
        The smallest period ``d`` with which a converged orbit repeats itself.

        Args:
            orbit: A ``(p, 2)`` array, the points of a period-p cycle in map order.
            rtol: Two points are the same when they differ by at most
                ``rtol * max(1, max|orbit|)`` in every coordinate.

        Returns:
            int: The smallest divisor ``d`` of ``p`` with ``orbit`` invariant
                under a shift by ``d``; ``p`` when no proper divisor works.

        Note:
            Only divisors of ``p`` need checking: on a cycle, ``x_i = x_j`` forces
            ``x_{i+n} = x_{j+n}`` for every ``n``, so the orbit repeats with period
            ``j - i`` as well as ``p``, hence with their gcd, which divides ``p``.
        """

        orbit = np.atleast_2d(orbit)
        period = len(orbit)
        tol = rtol * max(1.0, float(np.max(np.abs(orbit))))
        return next(
            d for d in range(1, period + 1)
            if period % d == 0
            and np.max(np.abs(np.roll(orbit, -d, axis=0) - orbit)) <= tol
        )

    def compute_fixed_point(
        self, initial_guess: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """
        Computes the coordinates of a fixed point based on the initial guess. Uses
        a multipoint shooting Newton's method.

        Args:
            initial_guess (np.ndarray): A (period, 2) array representing an
                initial guess for the orbit.

        Returns:
            np.ndarray: A (period, 2) array representing the converged fixed point
                coordinates.

        Raises:
            ValueError: If MINPACK reports anything other than convergence
                (``ier != 1``); the message carries MINPACK's own diagnosis.
                Returning the non-converged vector would give a "fixed point"
                that is not on any orbit.
        """

        initial_guess_flattened = self.flatten_trajectory(initial_guess)

        # Multipoint shooting couples every coordinate (residual_i = map(x_{i-1}) -
        # x_i), so this is a multivariate root-find. ``scipy.optimize.newton`` on an
        # array solves each coordinate as an INDEPENDENT scalar problem, which only
        # "converges" when the guess is already the answer (e.g. the exact orbit at
        # k=2) and silently returns a non-orbit otherwise -- the source of the
        # k-sensitive, eigenvector-flip-requiring behaviour. ``fsolve`` (MINPACK
        # hybrd) solves the coupled system properly and, staying near the guess,
        # also keeps the orbit labelling consistent with ``initial_guess``.
        shape = np.shape(initial_guess_flattened)
        solution, _info, ier, message = fsolve(
            lambda x: np.ravel(
                self.multipoint_shoot_flattened_difference(np.reshape(x, shape))
            ),
            np.ravel(initial_guess_flattened),
            xtol=1e-13,
            maxfev=10000,
            full_output=True,
        )
        if ier != 1:
            raise ValueError(
                "fsolve did not converge to a fixed point from the initial guess "
                f"{np.ravel(initial_guess).tolist()} (MINPACK ier={ier}): "
                f"{' '.join(str(message).split())}"
            )

        fixed_point_flattened = solution.reshape(shape)

        fixed_point_full = self.unflatten_trajectory(fixed_point_flattened)

        return np.array(fixed_point_full)

    def compute_eigenvectors(
        self,
        fixed_point: NDArray[np.float64],
        jacobians: Optional[list[NDArray[np.float64]]] = None,
    ) -> tuple[list[NDArray[np.float64]], list[list[NDArray[np.float64]]]]:
        """
        Computes the eigenvectors for each iterate of the fixed point.

        Args:
            fixed_point (np.ndarray): (period, 2) array of orbit coordinates.
            jacobians (list of np.ndarray, optional): Full-cycle Jacobians at
                each iterate. Automatically computes the Jacobians if not provided.

        Returns:
            tuple: (eigenvalues, eigenvectors), where:
                eigenvalues: Per iterate, the floats ``[unstable, stable]``.
                eigenvectors: Per iterate, the two ``(2,)`` arrays
                    ``[unstable, stable]``.

        Raises:
            ValueError: If any iterate is not a saddle (see
                :meth:`_validate_saddle`).

        Note:
            ``np.linalg.eig`` fixes the eigenvectors only up to sign. If the
            solver was built with an ``orient`` hook it is called here, once per
            orbit index, to re-sign the pair; otherwise the eigenvectors are
            returned exactly as the decomposition produced them and directions
            are chosen post hoc by
            :meth:`TangleWorkbench.orient_eigenvectors`.
        """

        if jacobians is None:
            jacobians = self.compute_jacobian(fixed_point)

        period, _ = np.shape(np.atleast_2d(fixed_point))

        eigenvector_list: list[list[NDArray[np.float64]]] = []
        eigenvalue_list: list[list[float]] = []

        for i in range(period):

            jacobian = jacobians[i]

            eigenvalues, eigenvectors = np.linalg.eig(jacobian)

            self._validate_saddle(i, eigenvalues)
            eigenvalues = np.real(eigenvalues)
            eigenvectors = np.real(eigenvectors)

            unstable_index = int(np.argmax(np.abs(eigenvalues)))

            unstable_vector = eigenvectors[:, unstable_index].reshape(2)
            stable_vector = eigenvectors[:, 1 - unstable_index].reshape(2)

            if self.orient is not None:
                unstable_vector, stable_vector = self.orient(
                    i, unstable_vector, stable_vector
                )

            eigenvector_list.append([unstable_vector, stable_vector])

            logger.debug("eigenvalues at orbit index %d: %s", i, eigenvalues)
            eigenvalue_list.append(
                [float(eigenvalues[unstable_index]), float(eigenvalues[1 - unstable_index])]
            )

        return eigenvalue_list, eigenvector_list

    @staticmethod
    def _validate_saddle(
        orbit_index: int, eigenvalues: NDArray[np.complex128]
    ) -> None:
        """
        Check that one full-cycle Jacobian belongs to a saddle.

        A saddle of an area-preserving map has two real eigenvalues, one of
        modulus above 1 and one below. Anything else (a complex pair at an
        elliptic point, or a degenerate pair at a parabolic one) means the guess
        converged onto a periodic orbit that has no manifolds to grow.

        Args:
            orbit_index (int): Which iterate of the orbit these eigenvalues
                belong to; used in the error message.
            eigenvalues (np.ndarray): The two eigenvalues of the full-cycle
                Jacobian at that iterate.

        Raises:
            ValueError: If the eigenvalues are not those of a saddle. This is a
                ValueError rather than an assertion because it is an input
                error — the caller's initial guess converged to a non-saddle —
                not a broken internal invariant.
        """
        values = np.ravel(np.asarray(eigenvalues))

        scale = np.maximum(np.abs(values), 1.0)
        if np.any(np.abs(np.imag(values)) > _IMAG_TOL * scale):
            raise ValueError(
                f"fixed point iterate {orbit_index} is not a saddle: eigenvalues "
                f"{values.tolist()} are complex, so the orbit is elliptic and has "
                "no stable/unstable manifolds"
            )

        moduli = np.abs(np.real(values))
        if not (np.max(moduli) > 1.0 and np.min(moduli) < 1.0):
            raise ValueError(
                f"fixed point iterate {orbit_index} is not a saddle: eigenvalues "
                f"{np.real(values).tolist()} do not have exactly one modulus above "
                "1 and one below"
            )

    def compute_partial_jacobians(
        self, fixed_point: NDArray[np.float64]
    ) -> list[NDArray[np.float64]]:
        """
        Computes the partial step jacobians for each step of the
        fixed point.

        Args:
            fixed_point (np.ndarray): A (period, 2) array representing the fixed point.

        Returns:
            list[np.ndarray]: List of (2, 2) Jacobians for single steps at each iterate.
        """

        period, _ = np.shape(np.atleast_2d(fixed_point))

        difference = np.abs(self.multipoint_shoot(fixed_point) - fixed_point)
        initial_step = (np.average(np.linalg.norm(difference, axis=1))) ** (1 / 3)

        partial_jacobians = [np.empty((2, 2)) for i in range(period)]

        for i in range(period):

            # compute the single step jacobian at each iterate in the fixed point
            x_i = fixed_point[i]

            if self.jacobian_function is None:
                # if a jacobian function isn't provided, compute numerically. The
                # finite-difference step reuses the orbit residual^(1/3); when the
                # orbit is located exactly (residual 0) that is 0, which makes the
                # difference stencil degenerate (nan). Fall back to scipy's own
                # default step in that case.
                if initial_step > 0:
                    jacobian = jacob(self.dynamical_map, x_i, initial_step=initial_step)
                else:
                    jacobian = jacob(self.dynamical_map, x_i)
                jacobian = jacobian.df

            else:
                jacobian = self.jacobian_function(x_i)

            partial_jacobians[i] = jacobian

        return partial_jacobians

    def compute_jacobian(
        self, fixed_point: NDArray[np.float64]
    ) -> list[NDArray[np.float64]]:
        """
        Computes the Jacobian for the fixed point.
        Computes p matrices for a period p fixed point.
        This routine computes the Jacobian for each step and then
        returns the product of those matrices for the full cycle jacobians.
        Will use a jacobian function if provided
        otherwise it will switch to a finite difference.

        Args:
            fixed_point (np.ndarray): A (period, 2) array representing the fixed point.

        Returns:
            list[np.ndarray]: List of (2, 2) full-cycle Jacobians, one per iterate.
        """

        period, _ = np.shape(np.atleast_2d(fixed_point))

        partial_jacobians = self.compute_partial_jacobians(fixed_point)

        true_jacobians = []
        # compute each cyclic permutation to get the full cycle jacobians
        for shift in range(period):

            factors = [partial_jacobians[(shift + i) % period] for i in range(period)]
            factors = factors[::-1]

            if period == 1:
                product = factors[0].copy()
            else:
                product = np.linalg.multi_dot(factors)
            true_jacobians.append(product)

        return true_jacobians

    def multipoint_shoot_flattened_difference(
        self, trajectory: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """
        Takes the difference between an iterate and the current trajectory.

        Args:
            trajectory (np.ndarray): (2 * period, 1) array of points to map forward.

        Returns:
            np.ndarray: (2 * period, 1) array.
        """

        shoot = self.multipoint_shoot_flattened(trajectory)

        difference = shoot - trajectory

        return difference

    def multipoint_shoot_flattened(
        self, trajectory: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """
        Takes a flattened trajectory and iterates it forward.

        Args:
            trajectory (np.ndarray): (2 * period, 1) array of points to shoot forward.

        Returns:
            np.ndarray: (2 * period, 1) array of the mapped trajectory.
        """

        trajectory_full = self.unflatten_trajectory(trajectory)

        trajectory_full_mapped = self.multipoint_shoot(trajectory_full)

        trajectory_flattened_mapped = self.flatten_trajectory(trajectory_full_mapped)

        return trajectory_flattened_mapped

    @staticmethod
    def flatten_trajectory(trajectory: NDArray[np.float64]) -> NDArray[np.float64]:
        """
        Takes a trajectory and makes it a 2n x 1 vector where n is
        the number of iterates.

        Args:
            trajectory (np.ndarray): (period, 2) array of points representing
                a trajectory.

        Returns
            np.ndarray: (2 * period, 1) array of points representing a
                flattened trajectory.
        """

        number_iterates, _ = np.shape(np.atleast_2d(trajectory))

        trajectory_reshaped = np.reshape(trajectory, (2 * number_iterates, 1))

        return trajectory_reshaped

    @staticmethod
    def unflatten_trajectory(trajectory: NDArray[np.float64]) -> NDArray[np.float64]:
        """
        Takes a flattened trajectory and reformats it back into a p x 2 matrix.

        Args:
            trajectory (np.ndarray): (2 * period, 1) array of points representing a
                flattened trajectory.

        Returns:
            np.ndarray: (period, 2) array of points representing an
                unflattened trajectory.
        """

        number_iterates_doubled, _ = np.shape(np.atleast_2d(trajectory))

        number_iterates = number_iterates_doubled // 2

        trajectory_reshaped = np.reshape(trajectory, (number_iterates, 2))

        return trajectory_reshaped

    def multipoint_shoot(self, trajectory: NDArray[np.float64]) -> NDArray[np.float64]:
        """
        Takes a trajectory and dynamical_maps it forward once.

        Args:
            trajectory (np.ndarray): (period, 2) array of points representing
                a trajectory.

        Returns:
            np.ndarray: (period, 2) array of the resulting trajectory.
        """

        period, _ = np.shape(np.atleast_2d(trajectory))

        maped_trajectory = np.array(
            [self.dynamical_map(trajectory[i, :]) for i in range(period)]
        )

        maped_trajectory = np.roll(maped_trajectory, 1, axis=0)

        return maped_trajectory
