"""The deep period-3 runs: 15 or 16 unstable steps, or 6 blasts.

Past the reference switch (BCG-4) the hole set is deeper, and the word of an
image-bridge class names the class of ITS image, which the minimal trellis
does not keep: that pair used to be a virtual ``new1`` (PK-1). Every pair of a
walk is the class of one piece of ``f(B)`` between two crossings, a real bridge
of ``W^u``, so the symbolic dynamics closes the alphabet under the word map
(WP7, D3). What holds there:

* invariant I1 -- all holes of one origin share their side of their bridge.
  A bridge's unstable branch is carried by ``Bridge.manifold_key``; inferring it
  from the endpoint intersections instead was unreliable on a period > 1 orbit
  (every branch's anchor sits at unstable cdist 0), so the backward propagation
  picked a container on the wrong branch and an orbit's holes flipped
  ``bridge_side`` partway down the chain (fixed 2026-10-02).
* every class walks (D3), and the alphabet is closed: every symbol of every
  word names a class of the table, and the dynamics is reliable.
"""

from __future__ import annotations

import pytest

from cases import BUILDERS, Case
from tanglepack.topology import check_holes_share_bridge_side

#: The deep runs (``tests/cases.py``): 15 and 16 unstable steps, 6 blasts.
DEEP_RUNS: tuple[str, ...] = ("p3_15", "p3_16", "p3_6_blasts")


def _build(name: str) -> Case:
    """A fresh build of one deep run."""
    return BUILDERS[name]()


@pytest.mark.parametrize("name", DEEP_RUNS)
def test_deep_p3_holes_share_bridge_side(name: str) -> None:
    """I1 holds deep into the tangle, and the symbolic dynamics runs."""
    case = _build(name)
    (fp,) = case.fixed_points
    trellis = case.session.trellis(fp)
    assert trellis.holes, f"{name} must punch holes"
    check_holes_share_bridge_side(
        trellis.holes, orientation_preserving=trellis.orientation_preserving
    )
    assert case.session.symbolic_dynamics().classes


@pytest.mark.parametrize("name", DEEP_RUNS)
def test_deep_p3_alphabet_is_closed(name: str) -> None:
    """Every symbol of every word is a class of the table, every class resolves, and
    the dynamics is reliable (no virtual class: D3)."""
    dyn = _build(name).session.symbolic_dynamics()
    for letter, word in dyn.rules.items():
        for symbol in word:
            assert symbol.bridge_class in dyn.table, (letter, symbol.text)
    assert dyn.unresolved == []
    assert dyn.is_reliable


@pytest.mark.parametrize("name", DEEP_RUNS)
def test_deep_p3_every_class_walks(name: str) -> None:
    """Every class gets an itinerary (D3): the walk exists once the openings follow
    the crossing signs (D6). The geometric openings left an element boundary out
    and made one class unreachable."""
    dyn = _build(name).session.symbolic_dynamics()
    assert not [cd for cd in dyn.classes.values() if cd.unresolved_reason]
