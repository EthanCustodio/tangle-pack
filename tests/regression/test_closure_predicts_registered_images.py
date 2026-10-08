"""Regression (PK-1, WP7): the classes the word map adds are real bridges' classes.

A pair of a walk is the class of one piece of ``f(B)`` between two crossings, a
real bridge of ``W^u``. On the 16-step period-3 run some of those pieces are
registered (the image of an image bridge): the classes the closure adds
(``named_by`` set) carry them as members, and their walked words must agree
with the registered images. Rule-based: no count, letter or word is pinned.
"""

from __future__ import annotations

from cases import BUILDERS
from tanglepack.topology.BridgeClass import _image_chain


def test_closure_classes_verify_against_their_registered_images() -> None:
    session = BUILDERS["p3_16"]().session
    trellis = session.trellis()
    dyn = session.symbolic_dynamics()
    closure = [cd for cd in dyn.classes.values() if cd.named_by is not None]
    evidenced = [
        cd
        for cd in closure
        if any(
            not member.is_loop and _image_chain(trellis, member.bridge_id)
            for member in cd.entry.members
        )
    ]
    assert evidenced, "no closure class has a member with a registered image chain"
    for cd in evidenced:
        assert cd.evidence, cd.letter
        assert cd.verified is True, (cd.letter, cd.word, [e.word for e in cd.evidence])
