from __future__ import annotations

from typing import Callable

import numpy as np
import numpy.typing as npt
from typing_extensions import Annotated

from numpy import float64

# A 2D point: exactly shape (2,) and dtype float64
Point2D = Annotated[npt.NDArray[float64], (2,)]

# Function type: takes a 2D point, returns a 2D point
MapFunc = Callable[[Point2D], Point2D]

Matrix2D = Annotated[npt.NDArray[float64], (2, 2)]

# Function type for jacobian function: takes a 2D point, returns a 2D matrix
JacFunc = Callable[[Point2D], Matrix2D]


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
    """

    def __init__(
        self,
        dynamical_map: MapFunc,
        dynamical_map_inverse: MapFunc,
        jacobian_function: JacFunc | None = None,
        name: str = "unnamed",
    ):
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
        """

        self.map = dynamical_map
        self.map_inv = dynamical_map_inverse
        self.jacobian = jacobian_function
        self.name = name

        # Cached "does this callable accept an (N, 2) batch?" flags, resolved
        # lazily on the first batched call (None = not yet probed).
        self._map_batchable: bool | None = None
        self._map_inv_batchable: bool | None = None

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
        device evaluation happens.

        Args:
            coords: ``(N, 2)`` array of points.

        Returns:
            ``(N, 2)`` array of images.
        """
        return self._apply_batch(self.map, coords, "_map_batchable")

    def map_inv_batch(self, coords: npt.NDArray[float64]) -> npt.NDArray[float64]:
        """Apply :attr:`map_inv` to a batch of points ``(N, 2) -> (N, 2)``.

        See :meth:`map_batch`.
        """
        return self._apply_batch(self.map_inv, coords, "_map_inv_batchable")

    def _apply_batch(
        self, fn: MapFunc, coords: npt.NDArray[float64], cache_attr: str
    ) -> npt.NDArray[float64]:
        coords = np.asarray(coords, dtype=float)
        if coords.shape[0] == 0:
            return coords.reshape(0, 2)

        batchable = getattr(self, cache_attr)
        if batchable is None:
            batchable = self._detect_batchable(fn, coords)
            setattr(self, cache_attr, batchable)

        if batchable:
            # fn takes/returns (2, N) (coordinate on axis 0); work in (N, 2).
            return np.asarray(fn(coords.T), dtype=float).T
        return np.vstack([np.asarray(fn(p), dtype=float) for p in coords])

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
