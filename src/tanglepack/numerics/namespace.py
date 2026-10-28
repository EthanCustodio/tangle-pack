"""The array namespace of a kernel's inputs: numpy itself, or CuPy for device arrays.

A kernel that calls :func:`array_namespace` computes where its inputs live and
returns an array of that namespace; it never moves data between host and device.
Every library caller hands such kernels host arrays, so on the library's own
paths the namespace is always ``numpy`` (see the ``tanglepack.numerics.gpu``
module docstring for which stages could run on a device and which are host-only).
"""

from __future__ import annotations

import sys
from types import ModuleType

import numpy as np


def array_namespace(*arrays: object) -> ModuleType:
    """
    The module a kernel computes with: ``numpy`` itself unless its inputs are CuPy arrays.

    Args:
        arrays: The kernel's array inputs.

    Returns:
        The ``cupy`` module when every input is a CuPy array, else ``numpy``.

    Raises:
        TypeError: Host and device arrays are mixed (a kernel never moves data).
    """
    cp = sys.modules.get("cupy")  # never imported here: no CuPy array can exist without it
    if cp is None:
        return np
    on_device = [isinstance(a, cp.ndarray) for a in arrays]
    if not any(on_device):
        return np
    if not all(on_device):
        raise TypeError("array_namespace: host and CuPy arrays mixed; move them explicitly")
    return cp
