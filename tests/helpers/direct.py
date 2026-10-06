"""Direct builds of every cached session product, and structural views of them.

:class:`~tanglepack.TangleSession` is a facade: it gathers partitions, holes
and strong pips from its per-fixed-point trellises and hands them to the
topology functions. :class:`DirectBuild` assembles the same products WITHOUT
any session cache -- from the public topology constructors over the
partitions, holes and pips read off the per-fixed-point trellises (the
current state) -- and :func:`session_view` / :func:`direct_view` reduce a
product to a structural, letter-free signature so two builds can be compared
by content, never by object identity.

Shared by ``tests/facade/test_session_equivalence.py`` (a session product
equals a direct build) and ``tests/facade/test_session_caches.py`` (it still
does after every cache-relevant event). This is the one place in the suite
that assembles the pipeline by hand.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Callable, Iterable, Optional

from helpers.names import (
    homotopy_pair,
    itinerary_pairs,
    matrix_by_classes,
    refined_words_in_names,
    word_in_names,
)
from tanglepack.loom.BridgeAlphabet import BridgeAlphabet
from tanglepack.topology.BridgeClass import bridge_classes
from tanglepack.topology.DualGraph import DualGraph
from tanglepack.topology.MinimalTrellis import minimal_trellis
from tanglepack.topology.PartitionFamily import (
    HomotopyPartition,
    IteratedHomotopyPartition,
)
from tanglepack.topology.SymbolicDynamics import symbolic_dynamics
from tanglepack.topology.Trellis import Trellis

if TYPE_CHECKING:  # pragma: no cover - typing only
    from cases import Case
    from tanglepack.numerics.FixedPoint import FixedPoint


class DirectBuild:
    """The pipeline of one partitioned case, assembled without the session caches.

    Attributes:
        case: The built case.
        selected: The fixed points of the selection.
        trellis: A fresh trellis of the workbench over the selection.
        partitions: Every per-fixed-point stable partition (of EVERY fixed
            point, as the session gathers them), one per ``(branch_key, side)``.
        holes: Every hole of the selected fixed points.
        strong_pips: The chosen strong pip of every selected fixed point that
            has one.
    """

    def __init__(
        self,
        case: "Case",
        fixed_points: Optional["FixedPoint | Iterable[FixedPoint]"] = None,
    ) -> None:
        """Read the per-fixed-point results off the case's session.

        Args:
            case: A case built through the ``partitioned`` stage.
            fixed_points: The selection (None for every fixed point), as the
                session accessors take it.
        """
        session = case.session
        self.case = case
        if fixed_points is None:
            self.selected = list(case.fixed_points)
        elif isinstance(fixed_points, Iterable):
            self.selected = list(fixed_points)
        else:
            self.selected = [fixed_points]
        self.trellis = Trellis.from_workbench(session.workbench, fixed_points)
        seen: dict[tuple, Any] = {}
        for fp in case.fixed_points:
            for result in session.trellis(fp).stable_partitions:
                seen.setdefault((result.branch_key, result.side), result)
        self.partitions = list(seen.values())
        self.holes = [hole for fp in self.selected for hole in session.trellis(fp).holes]
        self.strong_pips = [
            session.trellis(fp).strong_pip
            for fp in self.selected
            if session.trellis(fp).strong_pip is not None
        ]

    def homotopy(self) -> HomotopyPartition:
        """The homotopy partition family."""
        return HomotopyPartition.from_results(self.partitions, trellis=self.trellis)

    def classes(self):
        """The bridge-class table, lettered from a FRESH alphabet (letters are a
        session annotation; every view here is letter-free)."""
        table = bridge_classes(self.trellis, self.partitions)
        alphabet = BridgeAlphabet()
        for entry in table:
            if entry.active:
                entry.letter = alphabet.letter_for(entry.bridge_class)
        return table

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
        """The symbolic dynamics over the direct dual graph and direct classes."""
        return symbolic_dynamics(self.dual_graph(), self.classes())


# --------------------------------------------------------------------------- #
# structural views
# --------------------------------------------------------------------------- #
def trellis_view(trellis: Trellis) -> tuple:
    """Branch orderings and bridge ids."""
    return (
        sorted((repr(b.key), tuple(b.intersection_ids)) for b in trellis.branches.values()),
        sorted(b.id for b in trellis.bridges if b.id is not None),
    )


def arrangement_view(arrangement) -> list:
    """The regions by corner ids."""
    return sorted(tuple(region.corners) for region in arrangement.regions)


def classes_view(table) -> list:
    """Each class with its members, without the session's annotations."""
    return [(entry.bridge_class, tuple(entry.members)) for entry in table]


def family_view(family) -> list:
    """A partition family's signature, order-free (the session gathers its
    per-fixed-point partitions in trellis-cache order)."""
    return sorted(family.signature(), key=lambda entry: (repr(entry[0]), entry[1]))


def minimal_view(minimal) -> tuple:
    """Kept, hole and image bridge ids."""
    return (
        sorted(minimal.kept_bridge_ids),
        sorted(minimal.hole_bridge_ids),
        sorted(minimal.image_bridge_ids),
    )


def dual_view(graph: DualGraph) -> tuple:
    """Stable nodes with their elements, face corners and fundamental segments."""
    return (
        {key: node.elements for key, node in graph.stable_nodes.items()},
        [face.corners for face in graph.face_nodes],
        graph.fundamental_segments,
    )


def dynamics_view(dyn) -> tuple:
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


#: product -> (view of a product object, the direct build of that product).
_PRODUCTS: dict[str, tuple[Callable[[Any], Any], Callable[[DirectBuild], Any]]] = {
    "trellis": (trellis_view, lambda d: d.trellis),
    "arrangement": (arrangement_view, lambda d: d.trellis.arrangement),
    "homotopy_partition": (family_view, DirectBuild.homotopy),
    "bridge_classes": (classes_view, DirectBuild.classes),
    "minimal_trellis": (minimal_view, DirectBuild.minimal),
    "iterated_partition": (family_view, DirectBuild.iterated),
    "dual_graph": (dual_view, DirectBuild.dual_graph),
    "symbolic_dynamics": (dynamics_view, DirectBuild.dynamics),
}

#: Every product :func:`view_of` / :func:`direct_view` know.
PRODUCT_NAMES: tuple[str, ...] = tuple(_PRODUCTS)


def view_of(product: str, obj: Any) -> Any:
    """
    The structural view of one product object.

    Args:
        product: A name from :data:`PRODUCT_NAMES` (the session accessor's name).
        obj: The product, from the session or a direct build.

    Returns:
        A comparable, letter-free signature.
    """
    return _PRODUCTS[product][0](obj)


def direct_view(
    case: "Case",
    product: str,
    fixed_points: Optional["FixedPoint | Iterable[FixedPoint]"] = None,
) -> Any:
    """
    The view of ``product`` built directly from the case's current state.

    Args:
        case: A partitioned case.
        product: A name from :data:`PRODUCT_NAMES`.
        fixed_points: The selection (None for every fixed point).

    Returns:
        The view of a fresh, cache-free build.
    """
    return view_of(product, _PRODUCTS[product][1](DirectBuild(case, fixed_points)))
