"""The Hénon family, the one map every tanglepack script and test uses (notebooks still carry their own copies).

Dev Notes:

Before this module the same three functions were re-typed at the top of roughly
thirty files -- every script in ``scripts/``, most test modules, both
visualisation drivers -- each with its own ``k, b`` literals and its own choice
of ``np.array`` versus ``np.stack``. The two spellings are numerically identical
(``np.array([u, v])`` and ``np.stack([u, v], axis=0)`` build the same ``(2,)``
point and the same ``(2, N)`` batch), but only one of them advertises that the
map is batch-capable, so the copies quietly disagreed about the contract. They
are factories here instead: :func:`henon_map` closes over ``(k, b)`` and returns
a function with the axis-0 batch convention the whole numerics layer assumes
(coordinate on axis 0, so a single call maps a whole refinement layer, on the
GPU when one is enabled).

The map is

    (x, y) -> (y - k + x**2, -b * x)

which is area-preserving for ``|b| = 1`` -- the determinant of
:func:`henon_jacobian` is ``b``.
"""

from __future__ import annotations

import copy
from typing import Callable, Union

import numpy as np
from numpy.typing import NDArray

#: Parameters of the binary-horseshoe Hénon map: a single clean saddle at
#: ``[4, -4]`` whose tangle is a full binary horseshoe. The default fixture map.
HENON_K10: tuple[float, float] = (10, 1)

#: Parameters of the nested period-3 Hénon map: the period-3 saddle orbit inside
#: the period-1 tangle, the fixture for every nested / higher-period test.
HENON_P3: tuple[float, float] = (2, 1)

#: The fixed-point seeds each parameter pair is known to converge from. Keys are
#: the ``(k, b)`` pairs; values map a name to the guess
#: :meth:`~tanglepack.TangleWorkbench.construct_fixed_point` takes.
_SADDLE_GUESSES: dict[tuple[float, float], dict[str, list]] = {
    HENON_K10: {
        # (4.3166, -4.3166): both eigenvalues positive, no inversion.
        "saddle": [4, -4],
        # (-2.3166, 2.3166): both eigenvalues negative, so k_value = 2 * period
        # and the two manifold branches are a single dynamical chain.
        "inversion": [-2.3166, 2.3166],
    },
    HENON_P3: {
        "period_3": [[0, 1], [-1, 0], [-1, 1]],
        # (2.732, -2.732) -- the outer period-1 saddle of the nested tangle.
        "saddle": [4, -4],
    },
}

#: A map, its inverse, or a Jacobian: takes a ``(2,)`` point or a ``(2, N)``
#: batch of points (coordinate on axis 0).
MapFunction = Callable[[NDArray[np.float64]], NDArray[np.float64]]


def henon_map(k: float, b: float = 1) -> MapFunction:
    """
    Build the Hénon map ``(x, y) -> (y - k + x**2, -b * x)``.

    Args:
        k (float): The nonlinearity parameter.
        b (float): The Jacobian determinant. Defaults to 1 (area-preserving).

    Returns:
        MapFunction: The map, batch-capable on axis 0 -- it takes a ``(2,)``
        point or a ``(2, N)`` batch and returns the same shape.
    """

    def _map(point: NDArray[np.float64]) -> NDArray[np.float64]:
        x, y = point[0], point[1]
        return np.stack([y - k + x**2, -b * x], axis=0)

    return _map


def henon_map_inverse(k: float, b: float = 1) -> MapFunction:
    """
    Build the inverse Hénon map ``(x, y) -> (-y / b, x + k - y**2 / b**2)``.

    Args:
        k (float): The nonlinearity parameter.
        b (float): The Jacobian determinant. Defaults to 1 (area-preserving).

    Returns:
        MapFunction: The inverse map, batch-capable on axis 0.
    """

    def _map_inverse(point: NDArray[np.float64]) -> NDArray[np.float64]:
        x, y = point[0], point[1]
        return np.stack([-y / b, x + k - (y**2) / (b**2)], axis=0)

    return _map_inverse


def henon_jacobian(k: float, b: float = 1) -> MapFunction:
    """
    Build the Jacobian ``[[2x, 1], [-b, 0]]`` of the Hénon map.

    Args:
        k (float): The nonlinearity parameter. Unused -- the Jacobian does not
            depend on it; accepted so the three factories share a signature.
        b (float): The Jacobian determinant. Defaults to 1 (area-preserving).

    Returns:
        MapFunction: The Jacobian at a single ``(2,)`` point. Unlike the map and
        its inverse this one is NOT batched: the solver evaluates it point by
        point.
    """

    def _jacobian(point: NDArray[np.float64]) -> NDArray[np.float64]:
        x = point[0]
        return np.array([[2 * x, 1], [-b, 0]])

    return _jacobian


def saddle_guesses(
    k: float, b: float = 1
) -> dict[str, Union[list[float], list[list[float]]]]:
    """
    The fixed-point seeds known to converge for this parameter pair.

    Args:
        k (float): The nonlinearity parameter.
        b (float): The Jacobian determinant. Defaults to 1.

    Returns:
        dict[str, list]: A fresh copy of the named guesses. ``HENON_K10`` has
        ``"saddle"`` (the plain saddle at ``(4.3166, -4.3166)``) and
        ``"inversion"`` (the saddle at ``(-2.3166, 2.3166)``, whose eigenvalues
        are both negative). ``HENON_P3`` has ``"period_3"`` (the three-point
        orbit, as the multipoint guess ``construct_fixed_point`` expects) and
        ``"saddle"`` (the outer period-1 saddle of the nested tangle).

    Raises:
        KeyError: No guesses are recorded for ``(k, b)``.
    """
    try:
        known = _SADDLE_GUESSES[(k, b)]
    except KeyError:
        raise KeyError(
            f"no recorded saddle guesses for (k, b) = ({k}, {b}); "
            f"known pairs are {sorted(_SADDLE_GUESSES)}"
        ) from None

    return copy.deepcopy(known)
