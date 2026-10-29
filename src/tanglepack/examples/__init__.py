"""Reference dynamical systems shared by the tests, scripts and notebooks.

The Hénon family is exported here; the forced Duffing oscillator, the
expensive map of the ``workers`` option, is imported from
:mod:`.duffing`. Import the factories rather than
re-typing a map: a script that defines its own copy is a map the test suite
cannot pin.
"""

from __future__ import annotations

from .henon import (
    HENON_K10,
    HENON_K4,
    HENON_P3,
    MapFunction,
    henon_jacobian,
    henon_map,
    henon_map_inverse,
    saddle_guesses,
)

__all__ = [
    "HENON_K10",
    "HENON_K4",
    "HENON_P3",
    "MapFunction",
    "henon_jacobian",
    "henon_map",
    "henon_map_inverse",
    "saddle_guesses",
]
