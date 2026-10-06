"""Every cached session product equals the same product built directly.

:class:`~tanglepack.TangleSession` is a facade: it gathers partitions, holes
and strong pips from its per-fixed-point trellises and hands them to the
topology functions. This module rebuilds each product WITHOUT the session's
caches -- from the public topology constructors over the partitions, holes and
pips read off the per-fixed-point trellises -- and checks the session's result
against it, on the ``k10`` case and on the two-fixed-point ``nested`` case
(author decision 6, 2026-10-05). This is the one place in the suite that
assembles the pipeline by hand.

Comparisons are structural and letter-free: ids within ONE build, element
names, and words as oriented element pairs. Hit/miss behaviour is
``test_session_caches.py``.
"""

from __future__ import annotations

from typing import Any, Callable

import pytest

from cases import Case, build_k10, build_nested
from helpers.names import (
    homotopy_pair,
    itinerary_pairs,
    matrix_by_classes,
    refined_words_in_names,
    word_in_names,
)
from tanglepack.topology.BridgeClass import bridge_classes
from tanglepack.topology.DualGraph import DualGraph
from tanglepack.topology.MinimalTrellis import minimal_trellis
from tanglepack.topology.PartitionFamily import (
    HomotopyPartition,
    IteratedHomotopyPartition,
)
from tanglepack.topology.SymbolicDynamics import symbolic_dynamics
from tanglepack.topology.Trellis import Trellis

_BUILDERS: dict[str, Callable[[], Case]] = {
    "k10": build_k10,
    "nested": build_nested,
}


# --------------------------------------------------------------------------- #
# the direct build
# --------------------------------------------------------------------------- #
class _Direct:
    """The pipeline of one partitioned case, assembled without the session caches.

    Attributes:
        case: The built case.
        trellis: A fresh all-fixed-points trellis of the workbench.
        partitions: Every per-fixed-point stable partition, one per
            ``(branch_key, side)``.
        holes: Every per-fixed-point hole.
        strong_pips: The chosen strong pip of every fixed point that has one.
    """

    def __init__(self, case: Case) -> None:
        """Read the per-fixed-point results off the case's session.

        Args:
            case: A case built through the ``partitioned`` stage.
        """
        session = case.session
        self.case = case
        self.trellis = Trellis.from_workbench(session.workbench, None)
        seen: dict[tuple, Any] = {}
        for fp in case.fixed_points:
            for result in session.trellis(fp).stable_partitions:
                seen.setdefault((result.branch_key, result.side), result)
        self.partitions = list(seen.values())
        self.holes = [hole for fp in case.fixed_points for hole in session.trellis(fp).holes]
        self.strong_pips = [
            session.trellis(fp).strong_pip
            for fp in case.fixed_points
            if session.trellis(fp).strong_pip is not None
        ]

    def homotopy(self) -> HomotopyPartition:
        """The homotopy partition family."""
        return HomotopyPartition.from_results(self.partitions, trellis=self.trellis)

    def classes(self):
        """The bridge-class table (unlettered: letters are a session annotation)."""
        return bridge_classes(self.trellis, self.partitions)

    def minimal(self):
        """The minimal trellis."""
        return minimal_trellis(self.trellis, self.holes, self.classes())

    def iterated(self) -> IteratedHomotopyPartition:
        """The iterated homotopy partition."""
        return IteratedHomotopyPartition.from_minimal(self.minimal(), self.homotopy())

    def dual_graph(self) -> DualGraph:
        """The dual graph over the minimal trellis and the iterated partition."""
        minimal = self.minimal()
        iterated = IteratedHomotopyPartition.from_minimal(minimal, self.homotopy())
        return DualGraph(minimal, iterated, strong_pips=self.strong_pips)

    def dynamics(self):
        """The symbolic dynamics over the direct dual graph and the session's
        lettered table (the table's own equivalence is checked separately)."""
        return symbolic_dynamics(self.dual_graph(), self.case.session.bridge_classes())


# --------------------------------------------------------------------------- #
# structural views
# --------------------------------------------------------------------------- #
def _trellis_view(trellis: Trellis) -> tuple:
    """Branch orderings and bridge ids."""
    return (
        sorted((repr(b.key), tuple(b.intersection_ids)) for b in trellis.branches.values()),
        sorted(b.id for b in trellis.bridges if b.id is not None),
    )


def _arrangement_view(arrangement) -> list:
    """The regions by corner ids."""
    return sorted(tuple(region.corners) for region in arrangement.regions)


def _classes_view(table) -> list:
    """Each class with its members, without the session's annotations."""
    return [(entry.bridge_class, tuple(entry.members)) for entry in table]


def _family_view(family) -> list:
    """A partition family's signature, order-free (the session gathers its
    per-fixed-point partitions in trellis-cache order)."""
    return sorted(family.signature(), key=lambda entry: (repr(entry[0]), entry[1]))


def _minimal_view(minimal) -> tuple:
    """Kept, hole and image bridge ids."""
    return (
        sorted(minimal.kept_bridge_ids),
        sorted(minimal.hole_bridge_ids),
        sorted(minimal.image_bridge_ids),
    )


def _dual_view(graph: DualGraph) -> tuple:
    """Stable nodes with their elements, face corners and fundamental segments."""
    return (
        {key: node.elements for key, node in graph.stable_nodes.items()},
        [face.corners for face in graph.face_nodes],
        graph.fundamental_segments,
    )


def _dynamics_view(dyn) -> tuple:
    """Every class's itinerary, word and refined words, plus both matrices,
    spelled in element names."""
    classes = sorted(
        (
            homotopy_pair(dyn, cd),
            tuple(itinerary_pairs(dyn, cd) or ()),
            tuple(word_in_names(dyn, cd)),
            tuple(sorted((k, tuple(v)) for k, v in refined_words_in_names(dyn, cd).items())),
        )
        for cd in dyn.classes.values()
    )
    return (
        classes,
        matrix_by_classes(dyn, refined=True),
        matrix_by_classes(dyn, refined=False),
    )


# --------------------------------------------------------------------------- #
# product -> (session result view, direct result view)
# --------------------------------------------------------------------------- #
_PRODUCTS: dict[str, Callable[[Case, _Direct], tuple[Any, Any]]] = {
    "trellis": lambda c, d: (_trellis_view(c.session.trellis()), _trellis_view(d.trellis)),
    "arrangement": lambda c, d: (
        _arrangement_view(c.session.arrangement()),
        _arrangement_view(d.trellis.arrangement),
    ),
    "homotopy_partition": lambda c, d: (
        _family_view(c.session.homotopy_partition()),
        _family_view(d.homotopy()),
    ),
    "bridge_classes": lambda c, d: (
        _classes_view(c.session.bridge_classes()),
        _classes_view(d.classes()),
    ),
    "minimal_trellis": lambda c, d: (
        _minimal_view(c.session.minimal_trellis()),
        _minimal_view(d.minimal()),
    ),
    "iterated_partition": lambda c, d: (
        _family_view(c.session.iterated_partition()),
        _family_view(d.iterated()),
    ),
    "dual_graph": lambda c, d: (_dual_view(c.session.dual_graph()), _dual_view(d.dual_graph())),
    "symbolic_dynamics": lambda c, d: (
        _dynamics_view(c.session.symbolic_dynamics()),
        _dynamics_view(d.dynamics()),
    ),
}


@pytest.mark.parametrize("case_name", list(_BUILDERS))
@pytest.mark.parametrize("product", list(_PRODUCTS))
def test_session_product_equals_a_direct_build(product: str, case_name: str) -> None:
    """The session's cached product is the direct build's, structurally."""
    case = _BUILDERS[case_name]()
    direct = _Direct(case)

    from_session, from_direct = _PRODUCTS[product](case, direct)

    assert from_session == from_direct
    assert from_session, f"{product} of {case_name} is empty"

