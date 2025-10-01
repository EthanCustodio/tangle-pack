# tanglepack_webdash/maps.py
from __future__ import annotations
import numpy as np
from typing import Callable, Dict, Tuple

MapPair = Tuple[Callable[[np.ndarray], np.ndarray], Callable[[np.ndarray], np.ndarray]]


def henon_binary() -> MapPair:
    """Your working example from the script."""
    k, b = 10.0, 1.0

    def f(point: np.ndarray) -> np.ndarray:
        x, y = float(point[0]), float(point[1])
        return np.array([y - k + x * x, -b * x], dtype=float)

    def finv(point: np.ndarray) -> np.ndarray:
        x, y = float(point[0]), float(point[1])
        return np.array([-y / b, x + k - (y * y) / (b * b)], dtype=float)

    return f, finv


# Add more presets here as you like:
# def standard_map(...): ...
# def baker_map(...): ...

PRESETS: Dict[str, tuple[str, Callable[[], MapPair]]] = {
    "henon_binary": ("Hénon (binary horseshoe)", henon_binary),
    # "standard_map": ("Standard map", standard_map),
}
