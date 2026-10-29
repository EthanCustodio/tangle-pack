"""Optional GPU acceleration of the batched map evaluation and the same-stability sweep.

:func:`enable_gpu` does two things. It wraps a
:class:`~tanglepack.numerics.DynamicalSystem`'s ``map``/``map_inv`` callables
(unless ``route_map=False``) so that a batch of points is evaluated on the GPU
via `CuPy <https://cupy.dev>`_ and the result copied back to host memory, and it
marks the system so that ``numerics/diagnostics.py``'s same-stability sweep
builds and pre-tests its candidate pairs on the device (:func:`gpu_enabled`).
The point store, refinement's bookkeeping, the R-tree and crossing resolution
stay on the CPU.

Only the *vectorized* map call is accelerated: every per-layer refinement sweep
and every manifold iteration maps a whole ``(2, N)`` batch in a single call (see
:meth:`DynamicalSystem.map_batch`), and that single call is what runs on the
device. Per-point evaluations (a single ``(2,)`` point, as used by the
initializer and the fixed-point solver) stay on the CPU to avoid paying the
host<->device transfer for two numbers. While the map is wrapped,
``DynamicalSystem``'s ``workers`` split is skipped: the device gets the whole
batch.

What runs where
---------------
* **Device-capable.** The user's map, by NEP-18 under :func:`enable_gpu` (the
  only place the library moves data between host and device), and
  :func:`~tanglepack.numerics.geometry.curvature_area_batch` /
  :func:`~tanglepack.numerics.geometry.chord_deviation_batch`, which compute in
  their inputs' namespace (``numerics/namespace.py``'s ``array_namespace``):
  CuPy in, CuPy out, a mix raises ``TypeError``. No library caller hands them
  device arrays.
* **Host-only**, each for its own reason:

  - the merge (``ManifoldMachine.merge_manifolds``): it writes chains and the
    store;
  - refinement's layer loop (``ManifoldMachine._refine_layer``): it reads the
    store and appends to it between kernel calls, so on a device every layer
    would upload its four ``(M, 2)`` arrays and download its flags;
  - the segment build (``Tangle.add_manifolds`` and ``Tangle._query_curve``):
    libspatialindex runs only on the CPU;
  - crossing resolution (``Tangle._find_true_intersection``): a 2x2 LAPACK
    solve per candidate pair, not a batched kernel.

  The point store and the intersection registry stay host arrays.

* **The same-stability sweep** (``numerics/diagnostics.py``'s ``_hits``, the
  ``diagnose`` count of every blast) runs on the device while
  :func:`gpu_enabled`: the segment rows are gathered on the host and uploaded
  once, and the edge dedup, the sweep axis and sort, the strips and the
  candidate pairs are built there, pre-tested by one fused kernel whose margins
  are loosened so it keeps every pair the exact test keeps. Only the surviving
  pairs come back, and the host re-runs its exact test on them, so the listing
  and the count are byte-identical to the CPU's.

* **A library kernel on the device is parity to tolerance, never
  bit-identical.** CuPy's ``pow`` and ``hypot`` differ from numpy's in the
  last bit, so the curvature kernel on the device differs by up to 2.8e-7 (relative) on real refinement
  inputs (k=2.8 after two blasts) and moves no chord across ``area_cutoff``.
  Hénon under NEP-18 likewise leaves a different saved state with equal
  same-stability counts. A map written as below is bit-identical.

Measured result
---------------
Two things gain: the sweep, on any map, and a map written as one device kernel.

* **The sweep, 6.5x** (2026-10-29, the 1e-7 k=2.5 scale series, 32 blasts with
  ``diagnose=True``, cores 16-31, RTX 5060 Ti, medians of 3): with
  ``route_map=False`` the ``diagnose`` time falls from 38.7 s to 5.9 s and the
  blast series from 50.1 s to 17.4 s (2.9x), the saved state and every
  ``SameStabilityCount`` byte-equal to the CPU's. With the Hénon map routed it
  takes 18.0 s and the state is not byte-equal (see above): by 32 blasts the
  same-stability counts differ too.
* **Plain NEP-18 loses on an expensive map.** Every numpy call of the map body
  becomes its own kernel launch: leapfrog Duffing at 200 steps launches about
  1600 per call, about 10 ms whatever the batch size, and its growth takes
  2.35 s on the GPU against 0.27 s on the CPU (2026-10-28).
* **One fused kernel wins by about 100x** (:mod:`tanglepack.examples.duffing`,
  2026-10-28): 0.07 ms at 100 points, 1.3 ms at 1e5 and 12 ms at 1e6, against
  0.7, 54 and 1180 ms on the CPU, and bit-identical to the numpy map (2e5
  points, 0 differing; the whole saved Duffing state is byte-equal CPU against
  GPU, 2026-10-29). Duffing growth at 2000 steps per map call takes 0.56 s
  routed against 1.93 s on the CPU and 19.6 s under plain NEP-18. A fused map breaks even at about 10 points, so
  ``min_batch_points`` defaults to 16 (it was 64, when no batch broke even).
* **A cheap map gains nothing from routing.** The Hénon map is 0.06% of the 1e-7 scale
  series; the map, refinement, merge, segment build and crossing resolution together are
  3.5% of it, so even infinitely fast they would cap it at 1.037x (re-measured
  2026-10-27: 0.98-1.00x end to end). The numerics stress campaign
  (``stress/numerics/report/report.pdf``, 2026-10-07) measured 0.96-1.00x, the
  device path 88% transfer at ``N = 10^7`` and plain-CuPy break-even between
  ~6e5 and ~4e6 points per batch.

Turn it on whenever blasts run with ``diagnose=True``. Route the map for an
integration (an ODE flow) written as below; pass ``route_map=False`` for a
closed-form map such as Hénon.

Map requirements
----------------
For a map to run on the GPU at all, its body must use the NumPy *function* API
so the calls dispatch to CuPy under the array-function protocol (NEP-18):

* index the coordinate axis (``x, y = point`` or ``point[0]`` / ``point[1]``),
* build the result with ``np.stack([...], axis=0)`` -- **not** ``np.array([...])``
  (constructing an array from a Python list of CuPy scalars does not dispatch to
  the GPU and will raise).

The bundled Hénon maps follow this form: correct, and slow on a map of many
operations (see above).

Writing a device map
--------------------
To be fast, a map handed a CuPy array runs its WHOLE body as one kernel; to
leave the state byte-equal, that kernel must give the numpy map's floats:

* **one fused kernel**: a ``cupy.ElementwiseKernel`` holding the whole loop, one
  thread per point, with any table (a forcing term) uploaded once as a ``raw``
  argument;
* **``options=("--fmad=false",)``**: nvcc otherwise fuses ``a * b + c`` into one
  multiply-add, rounded once where numpy rounds twice;
* **numpy's operation order**, parenthesised as Python evaluates it:
  ``p + 0.5 * h * (q - q * q * q + f)`` is
  ``p + (0.5 * h) * ((q - q * q * q) + f)``.

Plain ``+ - * /`` round exactly alike on both sides, so such a kernel is
bit-identical; a transcendental (``exp``, ``sin``, ``pow``) is not (CuPy's
differ from numpy's in the last bit), so a map needing one computes it on the
host into a table, as the Duffing forcing does. Branch on the input's namespace
(``numerics/namespace.py``'s ``array_namespace``) and keep the numpy path as it
was. :func:`tanglepack.examples.duffing.duffing_map` is the worked example.

Usage
-----
Declare it once at the top of a script, before growing anything::

    import tanglepack

    session = tanglepack.TangleSession(my_map, my_map_inv, my_jac)
    tanglepack.enable_gpu(session)                  # a device map (Duffing): map + sweep
    tanglepack.enable_gpu(session, route_map=False) # Hénon: the sweep only, state byte-equal
    ...

``enable_gpu`` accepts a :class:`DynamicalSystem`, a ``TangleWorkbench`` (uses its
``dynamical_system``), or a ``TangleSession`` (uses its ``workbench``).
"""

from __future__ import annotations

import logging
from typing import Any, Callable

import numpy as np

from .DynamicalSystem import DynamicalSystem

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def _resolve_system(target: Any) -> DynamicalSystem:
    """Return the DynamicalSystem held by ``target`` (system/workbench/session)."""
    if isinstance(target, DynamicalSystem):
        return target
    system = getattr(target, "dynamical_system", None)
    if isinstance(system, DynamicalSystem):
        return system
    workbench = getattr(target, "workbench", None)
    if workbench is not None:
        system = getattr(workbench, "dynamical_system", None)
        if isinstance(system, DynamicalSystem):
            return system
    raise TypeError(
        "enable_gpu expected a DynamicalSystem, TangleWorkbench, or TangleSession; "
        f"got {type(target).__name__!r} with no reachable dynamical_system."
    )


def gpu_enabled(target: Any) -> bool:
    """
    Whether :func:`enable_gpu` put ``target``'s same-stability sweep on the GPU.

    Args:
        target: A :class:`DynamicalSystem`, ``TangleWorkbench`` or ``TangleSession``.

    Returns:
        True between :func:`enable_gpu` and :func:`disable_gpu`, whether or not
        the map was routed too (``route_map``).

    Raises:
        TypeError: ``target`` holds no reachable :class:`DynamicalSystem`.
    """
    return getattr(_resolve_system(target), "_gpu_enabled", False)


def _import_cupy():
    try:
        import cupy as cp  # type: ignore
    except Exception as exc:  # pragma: no cover - depends on environment
        raise ImportError(
            "GPU acceleration requires CuPy. Install it for your CUDA version, "
            "e.g. `pip install tanglepack[gpu]` or `pip install cupy-cuda12x`."
        ) from exc
    return cp


def _wrap_for_gpu(fn: Callable, cp, min_batch_points: int) -> Callable:
    """Wrap ``fn`` so large batches run on the GPU and small inputs stay on CPU."""

    def gpu_fn(point):
        arr = np.asarray(point)
        # A single point (ndim <= 1) or a tiny batch is cheaper on the CPU than
        # paying for a host<->device round trip. Coordinate is on axis 0, so the
        # number of points is the product of the trailing dimensions.
        n_points = int(arr.size // 2) if arr.ndim >= 2 else 1
        if arr.ndim <= 1 or n_points < min_batch_points:
            return fn(point)
        out = fn(cp.asarray(arr))  # numpy-API ops dispatch to CuPy
        return cp.asnumpy(out)

    return gpu_fn


def enable_gpu(
    target: Any, min_batch_points: int = 16, *, route_map: bool = True
) -> DynamicalSystem:
    """Put the same-stability sweep, and unless told not to the map, on the GPU.

    Wraps ``map`` and ``map_inv`` in place when ``route_map``. Idempotent and
    reversible via :func:`disable_gpu`; a second call with another
    ``route_map`` replaces the first. The original CPU callables are stashed on
    the system so they can be restored.

    Args:
        target: a :class:`DynamicalSystem`, ``TangleWorkbench``, or
            ``TangleSession``.
        min_batch_points: batches with fewer points than this run on the CPU
            (the transfer overhead is not worth it for small inputs; a fused
            device map breaks even at about 10). Single points always run on
            the CPU.
        route_map: hand batches of the map to the device. False keeps the map
            (and ``workers``' split of it) on the host, so the state stays
            byte-equal to the CPU's for a map whose CuPy evaluation is not
            bit-identical (Hénon) while the sweep still runs on the device.

    Returns:
        The underlying :class:`DynamicalSystem` (for chaining).

    Raises:
        ImportError: CuPy cannot be imported.
        TypeError: ``target`` holds no reachable :class:`DynamicalSystem`.

    Note:
        Only a map written as one device kernel gains (see the module
        docstring's "Writing a device map"); for a closed-form map such as
        Hénon pass ``route_map=False``.
    """
    system = _resolve_system(target)
    cp = _import_cupy()

    _restore_cpu_map(system)
    if route_map:
        system._cpu_map = system.map
        system._cpu_map_inv = system.map_inv
        system.map = _wrap_for_gpu(system._cpu_map, cp, min_batch_points)
        system.map_inv = _wrap_for_gpu(system._cpu_map_inv, cp, min_batch_points)
        system._gpu_map = True
    system._gpu_enabled = True

    logger.info(
        "GPU acceleration enabled for system %r (route_map=%s, min_batch_points=%d); "
        "only a map written as one device kernel gains, see the "
        "tanglepack.numerics.gpu module docstring.",
        system.name,
        route_map,
        min_batch_points,
    )
    return system


def disable_gpu(target: Any) -> DynamicalSystem:
    """Restore the CPU map callables and the host sweep undone by :func:`enable_gpu`."""
    system = _resolve_system(target)
    _restore_cpu_map(system)
    system._gpu_enabled = False
    return system


def _restore_cpu_map(system: DynamicalSystem) -> None:
    """Put back the callables a routed :func:`enable_gpu` stashed, re-probing batchability."""
    if getattr(system, "_gpu_map", False):
        system.map = system._cpu_map
        system.map_inv = system._cpu_map_inv
        system._gpu_map = False
    system._map_batchable = None
    system._map_inv_batchable = None
