"""The deep period-3 runs: 15 or 16 unstable steps, or 6 blasts.

Past the reference switch (BCG-4) the hole set is deeper, and the word of an
image-bridge class names the class of ITS image, which the minimal trellis
does not keep: that pair used to be a virtual ``new1`` (PK-1). The symbolic
dynamics closes the alphabet under the word map (WP7, C9), and a pair no
computed bridge spans is a class all the same, walked from its elements'
landings and flagged "no computed bridge yet" (D11, WP8b). What holds there:

* invariant I1 -- all holes of one origin share their side of their bridge.
  A bridge's unstable branch is carried by ``Bridge.manifold_key``; inferring it
  from the endpoint intersections instead was unreliable on a period > 1 orbit
  (every branch's anchor sits at unstable cdist 0), so the backward propagation
  picked a container on the wrong branch and an orbit's holes flipped
  ``bridge_side`` partway down the chain (fixed 2026-10-02).
* every class walks (D3), and the alphabet is closed: every symbol of every
  word names a class of the table, every class the closure added is spanned by
  a registered bridge or flagged as spanned by none, and the dynamics is
  reliable.
* every non-mutating law of the law tier (WP10), on these runs and on
  ``p3_16_12`` and ``inversion_period2_10``.
"""

from __future__ import annotations

import logging

import pytest

from cases import BUILDERS, Case
from helpers.law_tier import build_products, run_layer
from helpers.laws import LAYERS
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
def test_deep_p3_alphabet_is_closed(name: str, caplog) -> None:
    """Every symbol of every word is a class of the table; every class the closure
    added is spanned by a registered bridge or flagged ``spanned`` False (D11), and
    the deep runs name at least one such class. Nothing is unresolved, no
    "no bridge spans" WARNING is logged, and the dynamics is reliable."""
    session = _build(name).session
    caplog.clear()
    with caplog.at_level(logging.WARNING):
        dyn = session.symbolic_dynamics(rebuild=True)
    evidence = session.minimal_trellis().classes
    for letter, word in dyn.rules.items():
        for symbol in word:
            assert symbol.bridge_class in dyn.table, (letter, symbol.text)
    for cd in dyn.classes.values():
        assert cd.named_by is None or (cd.bridge_class in evidence) == cd.entry.spanned, cd.letter
    assert any(not cd.entry.spanned for cd in dyn.classes.values()), "no class lacks a bridge"
    assert dyn.unresolved == [] and dyn.is_reliable
    assert not [r for r in caplog.records if "no bridge spans" in r.getMessage()]


@pytest.mark.parametrize("name", DEEP_RUNS)
def test_deep_p3_every_class_walks(name: str) -> None:
    """Every class gets a walk (D3): the walk exists once the openings follow the
    crossing signs (D6). The geometric openings left an element boundary out and
    made one class unreachable."""
    dyn = _build(name).session.symbolic_dynamics()
    for cd in dyn.classes.values():
        assert cd.itinerary is not None, cd.letter


@pytest.mark.parametrize("name", DEEP_RUNS + ("p3_16_12", "inversion_period2_10"))
def test_every_law_holds_on_the_deep_runs(name: str) -> None:
    """Every non-mutating law holds deep into the tangle (WP10: bridges bracketed
    across a shared segment, the partition and face-side laws restated)."""
    case = _build(name)
    build_products(case)
    for layer in LAYERS:
        run_layer(layer, case)
