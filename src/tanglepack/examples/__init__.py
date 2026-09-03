"""Reference dynamical systems shared by the tests, scripts and notebooks.

Only the Hénon family lives here today. Import the factories rather than
re-typing a map: a script that defines its own copy is a map the test suite
cannot pin.
"""

from __future__ import annotations

from .henon import (
    HENON_K10,
    HENON_P3,
    MapFunction,
    henon_jacobian,
    henon_map,
    henon_map_inverse,
    saddle_guesses,
)

__all__ = [
    "HENON_K10",
    "HENON_P3",
    "MapFunction",
    "henon_jacobian",
    "henon_map",
    "henon_map_inverse",
    "saddle_guesses",
]
