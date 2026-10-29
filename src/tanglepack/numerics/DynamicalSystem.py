"""
The map, its inverse and (optionally) its Jacobian.

:class:`DynamicalSystem` is the one place the library touches the user's map.
It stores the forward and inverse maps, derives batched versions of both (a
whole refinement layer is mapped in one call), and enforces the point/matrix
shapes the numerical layer relies on.

With ``workers > 1`` a large batch is split into ``workers`` chunks mapped on
that many cores, in a thread pool (a vectorised map whose numpy work releases
the GIL) or a forked process pool (a per-point Python map, e.g. ``solve_ivp``).
The chunks are mapped in order and concatenated, so the result is the serial
one exactly for an ELEMENTWISE map, one whose image of a point does not depend
on the other points of the batch. A batch-adaptive integrator (one step size
for the whole batch) is not elementwise: that is a property of the user's map,
not of tanglepack.
"""

from __future__ import annotations

import multiprocessing
import weakref
from concurrent.futures import BrokenExecutor, Executor, ProcessPoolExecutor, ThreadPoolExecutor
from functools import partial
from typing import Annotated, Callable, Literal, Optional

import numpy as np
import numpy.typing as npt

from numpy import float64

# A 2D point: exactly shape (2,) and dtype float64
Point2D = Annotated[npt.NDArray[float64], (2,)]

# Function type: takes a 2D point, returns a 2D point
MapFunc = Callable[[Point2D], Point2D]

Matrix2D = Annotated[npt.NDArray[float64], (2, 2)]

# Function type for jacobian function: takes a 2D point, returns a 2D matrix
JacFunc = Callable[[Point2D], Matrix2D]

MapExecutor = Literal["thread", "process"]

#: A batch of a batch-capable map in a thread pool with fewer points runs in the
#: caller: below it the hand-off costs more than the split saves (2026-10-29, 4-8
#: cores: the Henon map breaks even near 1e5 points, the leapfrog Duffing map
#: near 2e4).
MIN_PARALLEL_POINTS: int = 32768
#: The same for an expensive call: a map evaluated point by point, or any map in
#: a process pool (chosen for an expensive map). The forked pool's hand-off, about
#: 0.5 ms, is repaid from a dozen points of an ODE map (``solve_ivp``, ~2 ms a
#: point); a cheap map in a process pool loses at most that per call.
MIN_PARALLEL_POINTS_EXPENSIVE: int = 16

# The maps a forked worker process applies, installed by its initializer.
_CHILD_MAPS: dict[str, MapFunc] = {}


def _install_maps(maps: dict[str, MapFunc]) -> None:
    """Process-pool initializer: the forked child inherits ``maps``, never pickles them."""
    _CHILD_MAPS.update(maps)


def _map_points(
    fn: MapFunc, batchable: bool, coords: npt.NDArray[float64]
) -> npt.NDArray[float64]:
    """Map ``(N, 2)`` points in one batched call, or point by point."""
    if batchable:
        # fn takes/returns (2, N) (coordinate on axis 0); work in (N, 2).
        return np.asarray(fn(coords.T), dtype=float).T
    return np.vstack([np.asarray(fn(p), dtype=float) for p in coords])


def _child_map_points(
    name: str, batchable: bool, coords: npt.NDArray[float64]
) -> npt.NDArray[float64]:
    """:func:`_map_points` in a worker process, with the map its initializer installed."""
    return _map_points(_CHILD_MAPS[name], batchable, coords)


class DynamicalSystem:
    """
    Object which contains the map functions for a dynamical system. This object also
    stores and enforces the types that the mapping functions should handle.

    Attributes:
        map (MapFunc): The dynamical map of the system. Takes a point on the plane and
            returns a point on the plane.
        map_inv (MapFunc): The inverse dynamical map of the system. Takes a point on
            the plane and returns a point on the plane. map(map_inv()) is the identity.
        jacobian (JacFunc): Optional function which takes a point on the plane and
            returns the jacobian at that point as a 2x2 matrix. Useful to include if
            you have a fast method to compute the jacobian.
        name (str): Optional name of the dynamical system.
        workers (int): Cores a large batch is split across (1: serial).
        map_executor (MapExecutor): ``"thread"`` or ``"process"`` pool.
    """

    def __init__(
        self,
        dynamical_map: MapFunc,
        dynamical_map_inverse: MapFunc,
        jacobian_function: JacFunc | None = None,
        name: str = "unnamed",
        *,
        workers: int = 1,
        map_executor: MapExecutor = "thread",
    ) -> None:
        """
        Initalizes the system with the mapping functions.

        Args:
            dynamical_map (MapFunc): The dynamical map of the system.
            dynamical_map_inverse (MapFunc): The inverse dynamical map of the system.
            jacobian_function (JacFunc | None, optional): Optional function which takes
                a point on the plane and returns the jacobian at that point as a
                2x2 matrix. Defaults to None.
            name (str, optional): Optional name of the dynamical system.
                Defaults to "unnamed".
            workers (int, optional): See :attr:`workers`. Defaults to 1.
            map_executor (MapExecutor, optional): See :attr:`map_executor`.
                Defaults to ``"thread"``.
        """

        self.map = dynamical_map
        self.map_inv = dynamical_map_inverse
        self.jacobian = jacobian_function
        self.name = name

        # Cached "does this callable accept an (N, 2) batch?" flags, resolved
        # lazily on the first batched call (None = not yet probed).
        self._map_batchable: bool | None = None
        self._map_inv_batchable: bool | None = None

        # The pool is built on the first split batch, for the maps it was built with.
        self._pool: Optional[Executor] = None
        self._pool_maps: tuple = ()
        self._pool_finalizer: Optional[weakref.finalize] = None
        self._workers = 1
        self._map_executor: MapExecutor = "thread"
        self.workers = workers
        self.map_executor = map_executor

    @property
    def workers(self) -> int:
        """
        Cores a large batch is split across (at least :data:`MIN_PARALLEL_POINTS`
        points for a batch-capable map in threads, else
        :data:`MIN_PARALLEL_POINTS_EXPENSIVE`); 1 maps every batch in the caller.

        Setting it shuts the current pool down; the next split builds a new one.
        It is a runtime resource: ``save`` does not store it.
        """
        return self._workers

    @workers.setter
    def workers(self, value: int) -> None:
        if int(value) != value or value < 1:
            raise ValueError(f"workers must be a positive integer, got {value!r}")
        self._shutdown_pool()
        self._workers = int(value)

    @property
    def map_executor(self) -> MapExecutor:
        """
        The pool a split batch runs in: ``"thread"`` (default) for a vectorised
        map whose numpy work releases the GIL, ``"process"`` for a per-point
        Python map (a forked pool: the map is inherited, only the coordinate
        chunks are pickled). Setting it shuts the current pool down.

        Raises:
            ValueError: On an unknown executor, or ``"process"`` where the
                fork start method is unavailable.
        """
        return self._map_executor

    @map_executor.setter
    def map_executor(self, value: MapExecutor) -> None:
        if value not in ("thread", "process"):
            raise ValueError(f"map_executor must be 'thread' or 'process', got {value!r}")
        if value == "process" and "fork" not in multiprocessing.get_all_start_methods():
            raise ValueError("map_executor='process' needs the fork start method")
        self._shutdown_pool()
        self._map_executor = value

    def _shutdown_pool(self) -> None:
        """Shut the pool down (it is rebuilt lazily)."""
        if self._pool_finalizer is not None:
            self._pool_finalizer()
        self._pool, self._pool_maps, self._pool_finalizer = None, (), None

    def _executor(self) -> Executor:
        """The pool, built for the current maps (a map replaced since rebuilds it)."""
        maps = (self.map, self.map_inv)
        if self._pool is not None and self._pool_maps == maps:
            return self._pool
        self._shutdown_pool()
        if self.map_executor == "thread":
            pool: Executor = ThreadPoolExecutor(max_workers=self.workers)
        else:
            pool = ProcessPoolExecutor(
                max_workers=self.workers,
                mp_context=multiprocessing.get_context("fork"),
                initializer=_install_maps,
                initargs=({"map": self.map, "map_inv": self.map_inv},),
            )
        self._pool, self._pool_maps = pool, maps
        self._pool_finalizer = weakref.finalize(self, pool.shutdown)
        return pool

    def check_area_preserving(
        self, point: npt.ArrayLike, *, rtol: float = 1e-6
    ) -> bool | None:
        """
        Whether the map preserves area at ``point``: ``|det J| == 1`` within ``rtol``.

        Area preservation is the fundamental property every topological result
        rests on (CLAUDE.md, "Area Preservation"). The absolute value is
        compared, so an orientation-reversing map (``det J = -1``) passes.

        Args:
            point: A ``(2,)`` point at which to evaluate :attr:`jacobian`.
            rtol: Allowed relative deviation of ``|det J|`` from 1.

        Returns:
            True or False from the Jacobian, or None when the system has no
            Jacobian (nothing to check; finite differences would only test
            their own step size).
        """
        if self.jacobian is None:
            return None
        determinant = np.linalg.det(np.asarray(self.jacobian(np.asarray(point, float))))
        return bool(abs(abs(determinant) - 1.0) <= rtol)

    def map_batch(self, coords: npt.NDArray[float64]) -> npt.NDArray[float64]:
        """
        Apply :attr:`map` to a batch of points ``(N, 2) -> (N, 2)``.

        Maps follow the columns-of-points convention (coordinate on axis 0), the
        same one ``scipy.differentiate`` uses, so the batch is handed to the map
        transposed as ``(2, N)`` and the result transposed back. A map written
        as ``x, y = point; ...`` is therefore batch-capable for free. If the map
        does not accept the batch it falls back to a per-point loop so legacy
        scalar-only maps still work (just without the speed-up). When a GPU
        backend has wrapped :attr:`map`, this single batched call is where the
        device evaluation happens. Otherwise, with :attr:`workers` above 1, a
        large batch is split into that many chunks mapped in the pool and
        concatenated in order: the serial result for an elementwise map (see
        the module docstring).

        Args:
            coords: ``(N, 2)`` array of points.

        Returns:
            ``(N, 2)`` array of images.
        """
        return self._apply_batch("map", coords, "_map_batchable")

    def map_inv_batch(self, coords: npt.NDArray[float64]) -> npt.NDArray[float64]:
        """Apply :attr:`map_inv` to a batch of points ``(N, 2) -> (N, 2)``.

        See :meth:`map_batch`.
        """
        return self._apply_batch("map_inv", coords, "_map_inv_batchable")

    def _apply_batch(
        self, name: str, coords: npt.NDArray[float64], cache_attr: str
    ) -> npt.NDArray[float64]:
        coords = np.asarray(coords, dtype=float)
        if coords.shape[0] == 0:
            return coords.reshape(0, 2)

        fn = getattr(self, name)
        batchable = getattr(self, cache_attr)
        if batchable is None:
            batchable = self._detect_batchable(fn, coords)
            setattr(self, cache_attr, batchable)

        threads = self.map_executor == "thread"
        smallest = MIN_PARALLEL_POINTS if batchable and threads else MIN_PARALLEL_POINTS_EXPENSIVE
        # A map routed by enable_gpu gets the whole batch on the device.
        if self.workers == 1 or len(coords) < smallest or getattr(self, "_gpu_map", False):
            return _map_points(fn, batchable, coords)
        task = (
            partial(_map_points, fn, batchable)
            if threads
            else partial(_child_map_points, name, batchable)
        )
        chunks = np.array_split(coords, min(self.workers, len(coords)))
        try:
            return np.concatenate(list(self._executor().map(task, chunks)))
        except BrokenExecutor:
            # A worker died (killed, segfault): drop the pool so the next call builds a new one.
            self._shutdown_pool()
            raise

    @staticmethod
    def _detect_batchable(fn: MapFunc, coords: npt.NDArray[float64]) -> bool:
        """
        Decide whether ``fn`` accepts a ``(2, N)`` batch of column points.

        Calls ``fn`` on the transposed batch and verifies the first few columns
        against the per-point result. The column check is what stops a scalar map
        from being mistaken for batch-capable on an ``N == 2`` input, where it
        could coincidentally return a ``(2, 2)`` array.
        """
        cols = coords.T  # (2, N)
        try:
            out = np.asarray(fn(cols), dtype=float)
        except Exception:
            return False
        if out.shape != cols.shape:
            return False
        k = min(3, coords.shape[0])
        try:
            expected = np.vstack(
                [np.asarray(fn(coords[i]), dtype=float) for i in range(k)]
            )  # (k, 2)
        except Exception:
            return False
        return bool(np.allclose(out[:, :k].T, expected, rtol=1e-9, atol=1e-12))
