"""
Regression (WP9, NEST-1 / NEST-2, D7): a heteroclinic pair is a pseudoneighbor pair
of the tangle whose stable branch it sits on.

The ``heteroclinic`` case (period 1 + period 3, period 3 grown 16 unstable steps)
has a thin lobe of the inner orbit's unstable manifold A poking across the outer
saddle's stable manifold B inside B's fundamental window. Its two crossings are
consecutive on B's stable branch and on A's unstable branch, so they are a
reference pair of B: tested against A's crossings with A's branch cycle (the arc
is A's), holed in B's partition, and its backward chain runs along A's unstable
cycle into A's bridges. Such a pair exists only when both tangles are selected.

WP9c (D13, D14): a heteroclinic bridge is a bridge like any other. A class joining
two tangles has no member lobe, but the images of arcs on other unstable branches
are walls its walk may not cross, and they decide its tie; one unstable step further
(``heteroclinic_17u``, perf) the registered images of the connecting bridges agree.
"""

from __future__ import annotations

import logging

import pytest

from cases import BUILDERS


def test_heteroclinic_pairs_hole_the_stable_owners_partition(caplog):
    with caplog.at_level(logging.WARNING, logger="tanglepack.topology.Pseudoneighbor"):
        case = BUILDERS["heteroclinic"]()
        trellis = case.session.topology().trellis
    assert not [r for r in caplog.records if "missing from the branch cycle" in r.message]

    def unstable_owner(pair):
        return trellis.intersection(pair.as_tuple()[0]).manifold_a_key[0]

    heteroclinic = [
        pair for pair in trellis.reference_pseudoneighbors
        if unstable_owner(pair) is not pair.branch_key[0]
    ]
    assert heteroclinic, "no heteroclinic reference pair"
    for pair in heteroclinic:
        owner, stable_owner = unstable_owner(pair), pair.branch_key[0]
        assert pair.hole is not None
        assert {opening[3] for opening in pair.hole.openings} == {pair.branch_key}
        backward = [h for h in trellis.holes if h.origin == pair.as_tuple() and h.iterate < 0]
        assert backward, "the heteroclinic reference propagates no hole"
        for hole in backward:
            assert trellis.bridge_between(*hole.bounding_ids).fixed_point is owner

    alone = case.session.topology([stable_owner]).trellis
    assert all(
        trellis.intersection(p.as_tuple()[0]).manifold_a_key[0] is stable_owner
        for p in alone.pseudoneighbors
    )


def _connecting(dyn) -> list:
    """The class records whose two elements lie on two fixed points."""
    return [
        cd for cls, cd in dyn.classes.items()
        if cls.source.fixed_point is not cls.target.fixed_point
    ]


def test_a_tie_no_lobe_decides_is_broken_by_the_images_of_disjoint_arcs():
    """WP9c (D13, D14): a heteroclinic bridge is a bridge like any other. A class whose
    arcs join two stable branches has no member lobe, but the image of its arcs still
    cannot cross the images of arcs on other unstable branches, so the decided walks of
    those classes break its tie: every connecting class walks one way, and some class
    without a lobe owes that to the walls."""
    dyn = BUILDERS["heteroclinic"]().session.symbolic_dynamics()
    connecting = _connecting(dyn)
    assert connecting
    for cd in connecting:
        assert cd.search is not None and cd.search.status == "resolved", (
            cd.letter, cd.status)
    assert any(
        cd.search is not None and cd.search.rejected and not cd.search.constrained
        for cd in dyn.classes.values()
    )
    assert not [cd.letter for cd in dyn.classes.values() if cd.ambiguous]


@pytest.mark.perf
def test_registered_connecting_images_agree_with_the_wall_rule():
    """One unstable step further the connecting bridges' images are computed: every
    connecting class with registered evidence walks one way, and its registered images
    never contradict the walk (the walls pick the side growth confirms)."""
    dyn = BUILDERS["heteroclinic_17u"]().session.symbolic_dynamics()
    judged = [cd for cd in _connecting(dyn) if cd.evidence]
    assert judged
    for cd in judged:
        assert cd.search is not None and cd.search.status == "resolved", (
            cd.letter, cd.status)
        assert cd.verified is not False, cd.letter
