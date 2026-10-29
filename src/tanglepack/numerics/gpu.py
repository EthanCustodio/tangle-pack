"""Optional GPU acceleration for the batched map evaluation.

This module sits *beside* the numerics code: enabling it changes nothing in
``ManifoldMachine``, ``Tangle``, or any other file. It only wraps a
:class:`~tanglepack.numerics.DynamicalSystem`'s ``map``/``map_inv`` callables so
that, when they are handed a batch of points, the batch is evaluated on the GPU
via `CuPy <https://cupy.dev>`_ and the result is copied back to host memory. The
point store, refinement's bookkeeping, the R-tree and crossing resolution stay on
the CPU.

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

* **A device result is parity to tolerance, never bit-identical.** CuPy's
  ``pow`` and ``hypot`` differ from numpy's in the last bit, so the curvature
  kernel on the device differs by up to 2.8e-7 (relative) on real refinement
  inputs (k=2.8 after two blasts) and moves no chord across ``area_cutoff``.

Measured result
---------------
GPU acceleration is opt-in and gives **no gain** (re-measured 2026-10-27, two
runs):

* end to end the GPU runs at 0.98-1.00x the CPU speed (growth, intersections
  and bridges, median of three interleaved repeats, equal point counts: k=2.8
  at 11 steps 0.990x and 1.000x, period 3 0.981x and 0.977x, k=10 0.983x and
  0.983x);
* the map is 0.56-0.60% of k=2.8 unstable-growth time (cProfile) and 0.06% of
  the 1e-7 scale series (``stress/profile/profile_numerics.py --cases scale``);
* the map, refinement, merge, segment build and crossing resolution together are
  3.5% of that series, so even infinitely fast they would cap it at 1.037x; the
  time is in the blast's same-stability count (``diagnose``, 78%);
* the break-even batch is out of reach: k=2.8 at 12 steps already raises
  ``PointBudgetExceeded`` (2,000,000 points) at the escape.

The numerics stress campaign (``stress/numerics/report/report.pdf``, 2026-10-07)
measured 0.96-1.00x end to end, the map at most 0.84% of growth, the device path
88% transfer at ``N = 10^7`` and break-even between ~6e5 and ~4e6 points per
batch, so ``min_batch_points`` has nothing to tune. Leave it off.

Map requirements
----------------
For a map to run on the GPU its body must use the NumPy *function* API so the
calls dispatch to CuPy under the array-function protocol (NEP-18). In practice:

* index the coordinate axis (``x, y = point`` or ``point[0]`` / ``point[1]``),
* build the result with ``np.stack([...], axis=0)`` -- **not** ``np.array([...])``
  (constructing an array from a Python list of CuPy scalars does not dispatch to
  the GPU and will raise).

The bundled Hénon maps already follow this form.

Usage
-----
Declare it once at the top of a script, before growing anything::

    import tanglepack

    session = tanglepack.TangleSession(my_map, my_map_inv, my_jac)
    tanglepack.enable_gpu(session)      # everything downstream batches on GPU
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


def _import_cupy():
    try:
        import cupy as cp  # type: ignore
    except Exception as exc:  # pragma: no cover - depends on environment
        raise ImportError(
            "GPU acceleration requires CuPy. Install it for your CUDA version, "
            "e.g. `pip install tanglepack[gpu]` or `pip install cupy-cuda12x`."
        ) from exc
    return cp


def gpu_enabled(system: DynamicalSystem) -> bool:
    """Whether :func:`enable_gpu` routes ``system`` through the GPU."""
    return getattr(system, "_gpu_enabled", False)


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


def enable_gpu(target: Any, min_batch_points: int = 64) -> DynamicalSystem:
    """Route the system's batched map evaluations through the GPU.

    Wraps ``map`` and ``map_inv`` in place. Idempotent and reversible via
    :func:`disable_gpu`. The original CPU callables are stashed on the system so
    they can be restored.

    Args:
        target: a :class:`DynamicalSystem`, ``TangleWorkbench``, or
            ``TangleSession``.
        min_batch_points: batches with fewer points than this run on the CPU
            (the transfer overhead is not worth it for small inputs). Single
            points always run on the CPU.

    Returns:
        The underlying :class:`DynamicalSystem` (for chaining).

    Raises:
        ImportError: CuPy cannot be imported.
        TypeError: ``target`` holds no reachable :class:`DynamicalSystem`.

    Note:
        Measured as no end-to-end gain (see the module docstring); it stays
        opt-in and is not deprecated.
    """
    system = _resolve_system(target)
    cp = _import_cupy()

    # Stash the CPU callables once so repeated enable/disable stays correct.
    if not gpu_enabled(system):
        system._cpu_map = system.map
        system._cpu_map_inv = system.map_inv

    system.map = _wrap_for_gpu(system._cpu_map, cp, min_batch_points)
    system.map_inv = _wrap_for_gpu(system._cpu_map_inv, cp, min_batch_points)
    system._gpu_enabled = True

    # Force the batchable probe to re-run against the freshly wrapped callables.
    system._map_batchable = None
    system._map_inv_batchable = None

    logger.info(
        "GPU acceleration enabled for system %r (min_batch_points=%d); measured as "
        "no end-to-end gain (map <= 0.84%% of growth time), see the "
        "tanglepack.numerics.gpu module docstring.",
        system.name,
        min_batch_points,
    )
    return system


def disable_gpu(target: Any) -> DynamicalSystem:
    """Restore the original CPU map callables wrapped by :func:`enable_gpu`."""
    system = _resolve_system(target)
    if not gpu_enabled(system):
        return system

    system.map = system._cpu_map
    system.map_inv = system._cpu_map_inv
    system._gpu_enabled = False
    system._map_batchable = None
    system._map_inv_batchable = None
    return system
