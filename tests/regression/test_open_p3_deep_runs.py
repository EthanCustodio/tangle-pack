"""The deep period-3 runs: 15 or 16 unstable steps, or 6 blasts.

Past the reference switch (BCG-4) the hole set is deeper, and the word of an
image-bridge class names the class of ITS image, which the minimal trellis
does not keep: that pair used to be a virtual ``new1`` (PK-1). The symbolic
dynamics closes the alphabet under the word map over the classes registered
bridges span (WP7, C9); a pair no bridge spans is a construction error (D3),
never a class, and whether it must get its bridge is DEC-21. What holds there:

* invariant I1 -- all holes of one origin share their side of their bridge.
  A bridge's unstable branch is carried by ``Bridge.manifold_key``; inferring it
  from the endpoint intersections instead was unreliable on a period > 1 orbit
  (every branch's anchor sits at unstable cdist 0), so the backward propagation
  picked a container on the wrong branch and an orbit's holes flipped
  ``bridge_side`` partway down the chain (fixed 2026-10-02).
* every class walks (D3), and the alphabet is closed: every symbol of every
  word names a class of the table, every class the closure added has a
  registered bridge, and a walk pair no bridge spans unresolves its class with
  a WARNING, so the dynamics is not reliable.
"""

from __future__ import annotations

import logging

import pytest

from cases import BUILDERS, Case
from helpers.logs import assert_logged
from tanglepack.topology import check_holes_share_bridge_side
from tanglepack.topology.BridgeClass import oriented_class

_LOGGER = "tanglepack.topology.SymbolicDynamics"

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


def _unspanned_pairs(session, dyn, cd) -> list:
    """The pairs of a class's walk whose class no registered bridge spans."""
    if cd.search is None or cd.search.walk is None:
        return []
    known = set(dyn.table.classes) | set(session.minimal_trellis().classes.classes)
    walk = cd.search.walk.itinerary
    pairs = [
        oriented_class(dyn.naming.parent_of(x), dyn.naming.parent_of(y), dyn.fixed_points)
        for x, y in zip(walk[::2], walk[1::2])
    ]
    return [pair for pair, direction in pairs if direction and pair not in known]


@pytest.mark.parametrize("name", DEEP_RUNS)
def test_deep_p3_alphabet_is_closed(name: str, caplog) -> None:
    """Every symbol of every word is a class of the table, and every class the closure
    added is spanned by a registered bridge. A walk pair no bridge spans is a
    construction error (D3): its class is unresolved with one WARNING, and the
    dynamics is not reliable."""
    session = _build(name).session
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger=_LOGGER):
        dyn = session.symbolic_dynamics(rebuild=True)
    evidence = session.minimal_trellis().classes
    for letter, word in dyn.rules.items():
        for symbol in word:
            assert symbol.bridge_class in dyn.table, (letter, symbol.text)
    for cd in dyn.classes.values():
        assert cd.named_by is None or cd.bridge_class in evidence, cd.letter
    unspanned = [cd for cd in dyn.classes.values() if _unspanned_pairs(session, dyn, cd)]
    assert all(cd.itinerary is None for cd in unspanned)
    if unspanned:
        assert len(assert_logged(caplog, logging.WARNING, _LOGGER)) >= len(unspanned)
        assert not dyn.is_reliable


@pytest.mark.parametrize("name", DEEP_RUNS)
def test_deep_p3_every_class_walks(name: str) -> None:
    """Every class gets a walk (D3): the walk exists once the openings follow the
    crossing signs (D6). The geometric openings left an element boundary out and
    made one class unreachable. A walk naming a pair no bridge spans is kept on the
    search though the class is unresolved."""
    dyn = _build(name).session.symbolic_dynamics()
    for cd in dyn.classes.values():
        assert cd.itinerary is not None or (cd.search is not None and cd.search.walk), cd.letter
