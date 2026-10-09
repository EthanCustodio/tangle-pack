"""NUM-20: on a heteroclinic link each side of the crossing moves on its OWN saddle.

A crossing of ``W^u(A)`` and ``W^s(B)`` maps to a crossing of the same two
manifolds: its unstable cdist grows by A's per-step factor and its stable cdist
shrinks by B's, so the product scales by ``beta_A / beta_B`` (1 only when A is
B). The p1(u) x p3(s) recipes are the fixtures with such links (13u / 15s and
13u / 16s, blasted; the unblasted 13u / 16s crossings below already hold 12); the
iterate inference was right all along, and the two link laws used to read A's
factor for both sides (links ``7 -> 215``, ``5 -> 9`` and, here, ``4 -> 8``
failed).
"""

from __future__ import annotations

import cases
from helpers import laws


def test_heteroclinic_links_scale_by_their_own_saddles() -> None:
    case = cases.build_nested(
        outer_blasts=0, p3_unstable_steps=13, p3_stable_steps=16, through="intersected"
    )
    registry = case.registry
    heteroclinic = [
        source
        for source, n, _ in registry.iterate_table.items()
        if n == 1
        and not registry[source].is_synthetic
        and registry[source].manifold_a_key[0] is not registry[source].manifold_b_key[0]
    ]
    assert heteroclinic, "the recipe must register a heteroclinic iterate link"
    assert laws.check_area_along_iterate_links(case) > 0
    assert laws.check_iterate_link_scales_by_beta(case) > 0
