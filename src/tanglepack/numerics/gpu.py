"""Optional GPU acceleration for the batched map evaluation.

This module sits *beside* the numerics code: enabling it changes nothing in
``ManifoldMachine``, ``Tangle``, or any other file. It only wraps a
:class:`~tanglepack.numerics.DynamicalSystem`'s ``map``/``map_inv`` callables so
that, when they are handed a batch of points, the batch is evaluated on the GPU
via `CuPy <https://cupy.dev>`_ and the result is copied back to host memory. The
linked-list / refinement / intersection logic stays on the CPU.

Only the *vectorized* map call is accelerated: every per-layer refinement sweep
and every manifold iteration maps a whole ``(2, N)`` batch in a single call (see
:meth:`DynamicalSystem.map_batch`), and that single call is what runs on the
device. Per-point evaluations (a single ``(2,)`` point, as used by the
initializer and the fixed-point solver) stay on the CPU to avoid paying the
host<->device transfer for two numbers.

Measured result
---------------
GPU acceleration is opt-in and, as measured by the numerics stress campaign
(2026-10-07, experiments E1/E2 in ``stress/numerics/report/report.pdf``), gives
**no gain**:

* end to end the GPU runs at 0.96-1.00x the CPU speed (median of three repeats:
  k=2.8 0.999x, period 3 0.975x, k=10 0.958x);
* the map is at most 0.84% of unstable-growth CPU time (cProfile), so Amdahl's
  law caps any map accelerator at 1.008x; refinement and the linked list
  (Python ``Point`` objects, pointer walks) take the rest;
* the device path is transfer-dominated (88% of GPU time at ``N = 10^7``), and
  through ``map_batch`` it only breaks even somewhere between ~6e5 and ~4e6
  points per batch, so ``min_batch_points`` has nothing to tune.

It could pay only once refinement and the linked list themselves live on the
device as arrays, so that a growth step no longer round-trips through host
memory batch by batch. Until then leave it off.

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
    if not getattr(system, "_gpu_enabled", False):
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
    if not getattr(system, "_gpu_enabled", False):
        return system

    system.map = system._cpu_map
    system.map_inv = system._cpu_map_inv
    system._gpu_enabled = False
    system._map_batchable = None
    system._map_inv_batchable = None
    return system
