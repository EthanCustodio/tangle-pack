"""The forced Duffing oscillator: an expensive map, for the ``workers`` option.

The double-well Duffing oscillator driven periodically,

    q'' = q - q**3 + eps * cos(omega * t),

has a saddle at the origin when unforced, whose two homoclinic loops wind
round the wells at ``q = +-1``. Under the forcing the saddle becomes a
hyperbolic periodic orbit near the origin and the loops break into a
homoclinic tangle (Melnikov), as Hénon's. Its stroboscopic map (one forcing
period ``T = 2 pi / omega``) is area-preserving, and evaluating it is an ODE
integration: the case the multicore map split is for.

Two forms of the map:

* :func:`duffing_map` / :func:`duffing_map_inverse` integrate with leapfrog
  (kick-drift-kick), a composition of shears, so the map is EXACTLY symplectic
  (``det J = 1`` to rounding) and time-reversible: the inverse is the same
  scheme run with ``-h`` from ``T`` back to 0. Both are batch-capable on
  axis 0, as :mod:`.henon`'s. Handed a CuPy array (``enable_gpu``) they run
  the whole leapfrog loop as ONE device kernel, bit-identical to the numpy
  loop: the worked example of a device map in ``tanglepack.numerics.gpu``.
* :func:`duffing_map_ivp` / :func:`duffing_map_ivp_inverse` integrate one
  point at a time with ``scipy.integrate.solve_ivp`` (DOP853), a batch in a
  Python loop: the GIL-bound form ``map_executor="process"`` is for. They are
  area-preserving only to their tolerance.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from ..numerics.namespace import array_namespace
from .henon import MapFunction

#: Forcing amplitude: four steps of each manifold cross at ``area_cutoff = 1e-5``.
DUFFING_EPS: float = 0.2
#: Forcing frequency: the saddle stretches by about ``e^T ~ 66`` per period
#: ``T = 2 pi / omega``, and MINPACK converges on its orbit (at ``omega = 1``,
#: ``e^T ~ 535``, it stalls at rounding).
DUFFING_OMEGA: float = 1.5
#: Leapfrog steps per forcing period: a step of ``T / 200`` costs about 2 ms per
#: batch call, an expensive map, and still resolves the flow.
DUFFING_STEPS: int = 200

#: The leapfrog loop as one CUDA kernel, in numpy's operation order (``0.5 * h``
#: first, ``(q - q*q*q) + F``); compiled with ``--fmad=false`` so no
#: multiply-add is fused, which makes it bit-identical to the numpy loop.
_LEAPFROG_KERNEL = """
double q = q0, p = p0;
for (int i = 0; i < steps; i++) {
    p = p + (0.5 * h) * ((q - q * q * q) + forcing[i]);
    q = q + h * p;
    p = p + (0.5 * h) * ((q - q * q * q) + forcing[i + 1]);
}
q1 = q;
p1 = p;
"""


def _leapfrog(
    eps: float, omega: float, steps: int, direction: int
) -> MapFunction:
    """The leapfrog flow over one forcing period, forward (+1) or backward (-1)."""
    h = direction * 2 * np.pi / omega / steps
    # Forcing at the grid times i * T / steps, read from the same integer index
    # forward and backward, so the inverse retraces the forward map's floats.
    times = range(steps + 1) if direction > 0 else range(steps, -1, -1)
    forcing = [eps * np.cos(2 * np.pi * i / steps) for i in times]
    device = {}  # the compiled kernel and the forcing table, built on the first CuPy call

    def device_flow(cp, point):
        """The flow of a CuPy batch, by :data:`_LEAPFROG_KERNEL`: one launch per call."""
        if not device:
            device["kernel"] = cp.ElementwiseKernel(
                "raw float64 forcing, float64 h, int32 steps, float64 q0, float64 p0",
                "float64 q1, float64 p1",
                _LEAPFROG_KERNEL,
                "duffing_leapfrog",
                options=("--fmad=false",),
            )
            device["forcing"] = cp.asarray(forcing, dtype=cp.float64)
        q, p = device["kernel"](device["forcing"], h, steps, point[0], point[1])
        return cp.stack([q, p], axis=0)

    def flow(point: NDArray[np.float64]) -> NDArray[np.float64]:
        xp = array_namespace(point)
        if xp is not np:
            return device_flow(xp, point)
        q, p = point[0], point[1]
        for before, after in zip(forcing, forcing[1:]):
            p = p + 0.5 * h * (q - q * q * q + before)
            q = q + h * p
            p = p + 0.5 * h * (q - q * q * q + after)
        return np.stack([q, p], axis=0)

    return flow


def duffing_map(
    eps: float = DUFFING_EPS, omega: float = DUFFING_OMEGA, steps: int = DUFFING_STEPS
) -> MapFunction:
    """
    Build the stroboscopic map of the forced Duffing oscillator, by leapfrog.

    Args:
        eps (float): The forcing amplitude.
        omega (float): The forcing frequency; one map step is ``2 pi / omega``.
        steps (int): Leapfrog steps per forcing period.

    Returns:
        MapFunction: ``(q, p)`` at ``t = 0`` to ``(q, p)`` at ``t = T``,
        batch-capable on axis 0 and exactly symplectic.
    """
    return _leapfrog(eps, omega, steps, +1)


def duffing_map_inverse(
    eps: float = DUFFING_EPS, omega: float = DUFFING_OMEGA, steps: int = DUFFING_STEPS
) -> MapFunction:
    """
    Build the inverse of :func:`duffing_map`: the same leapfrog with ``-h``.

    Args:
        eps (float): The forcing amplitude.
        omega (float): The forcing frequency.
        steps (int): Leapfrog steps per forcing period.

    Returns:
        MapFunction: The inverse map, batch-capable on axis 0; it undoes
        :func:`duffing_map` to rounding.
    """
    return _leapfrog(eps, omega, steps, -1)


def _ivp(eps: float, omega: float, rtol: float, direction: int) -> MapFunction:
    """The flow over one forcing period by ``solve_ivp``, one point per integration."""
    period = 2 * np.pi / omega
    span = (0.0, period) if direction > 0 else (period, 0.0)

    def field(t: float, state: NDArray[np.float64]) -> list[float]:
        q, p = state
        return [p, q - q**3 + eps * np.cos(omega * t)]

    def flow(point: NDArray[np.float64]) -> NDArray[np.float64]:
        point = np.asarray(point, dtype=float)
        images = [
            solve_ivp(field, span, column, method="DOP853", rtol=rtol, atol=rtol).y[:, -1]
            for column in point.reshape(2, -1).T
        ]
        return np.stack(images, axis=1).reshape(point.shape)

    return flow


def duffing_map_ivp(
    eps: float = DUFFING_EPS, omega: float = DUFFING_OMEGA, rtol: float = 1e-10
) -> MapFunction:
    """
    Build the stroboscopic map by ``solve_ivp`` (DOP853), per point.

    Args:
        eps (float): The forcing amplitude.
        omega (float): The forcing frequency.
        rtol (float): The integrator's relative (and absolute) tolerance.

    Returns:
        MapFunction: The map on a ``(2,)`` point or a ``(2, N)`` batch, each
        point integrated on its own in Python (GIL-bound). Area-preserving
        only to ``rtol``.
    """
    return _ivp(eps, omega, rtol, +1)


def duffing_map_ivp_inverse(
    eps: float = DUFFING_EPS, omega: float = DUFFING_OMEGA, rtol: float = 1e-10
) -> MapFunction:
    """
    Build the inverse of :func:`duffing_map_ivp`: the flow from ``T`` back to 0.

    Args:
        eps (float): The forcing amplitude.
        omega (float): The forcing frequency.
        rtol (float): The integrator's tolerance.

    Returns:
        MapFunction: The inverse map, per point as :func:`duffing_map_ivp`.
    """
    return _ivp(eps, omega, rtol, -1)


def saddle_guess(eps: float = DUFFING_EPS, omega: float = DUFFING_OMEGA) -> list[float]:
    """
    A seed for the saddle's periodic orbit at ``t = 0``.

    The linearised response at the origin, ``q'' - q = eps cos(omega t)``, has
    the periodic solution ``q = -eps cos(omega t) / (1 + omega**2)``.

    Args:
        eps (float): The forcing amplitude.
        omega (float): The forcing frequency.

    Returns:
        list[float]: ``[q, p]`` for :meth:`~tanglepack.TangleWorkbench.construct_fixed_point`.
    """
    return [-eps / (1 + omega**2), 0.0]
