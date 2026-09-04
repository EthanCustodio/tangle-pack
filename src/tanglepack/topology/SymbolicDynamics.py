"""
The symbolic dynamics of a trellis: one symbol per bridge class, and its image.

A bridge class is a symbol (see :mod:`~tanglepack.topology.BridgeClass`), and one
map step carries a bridge onto a chain of bridges — so a symbol rewrites to a
WORD in the same alphabet. Collect one rule per symbol and the tangle's
substitution is complete: ``a -> ab``, a transition graph, a transition matrix,
and from those the usual growth rate and connectivity questions.

The word itself is read off the dual graph
(:meth:`~tanglepack.topology.DualGraph.DualGraph.walk_bridge`), which handles
the three shapes a bridge's image takes in a finite picture: fully registered,
registered-then-walked, or fully walked.

Dev Notes:

* Symbols are :class:`~tanglepack.topology.BridgeClass.BridgeClass` OBJECTS, and
  labels (``a``, ``b``, ``a1``) are presentation only. A label keys on the fixed
  point's PERIOD, so two distinct saddles of the same period would print the
  same element name; :meth:`SymbolicDynamics.display_label` disambiguates them
  by the fixed point's POSITION in the trellis instead, and that is what every
  report here prints. Never group or compare by a label.
* The order of the symbols is :func:`~tanglepack.topology.BridgeClass.class_sort_key`
  — fixed point, orbit, branch, side, element id — so two runs over the same
  trellis name the symbols identically and :meth:`SymbolicDynamics.describe` is
  string-comparable across builds. New symbols sort after the known ones.
* A **new symbol** is a letter that no computed bridge realises. It is not an
  error: it says the trellis is not closed under the map at this extent (the
  image of some bridge runs through a face pairing that the computed manifolds
  do not yet contain a bridge for). It is kept, labelled and logged, and it has
  no rule of its own — a rule needs a bridge to walk.
* D.4 (one word per class) is enforced strictly. Two bridges of one class whose
  images spell different words mean the partition is too coarse to be a symbolic
  alphabet, and that is a diagnostic worth an exception rather than a silent
  choice.

Open question: the substitution is under ONE map step, so on a period-k orbit
the symbols cycle between branches and the matrix is not primitive by
construction. Composing to ``f^k`` is a matrix power, not new machinery, but the
right alphabet for that is the per-branch one and nothing here chooses it.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Mapping, Sequence, TYPE_CHECKING, Union

import numpy as np
from numpy.typing import NDArray

from .BridgeClass import BridgeClass, class_sort_key
from .DualGraph import DEFAULT_MAX_PATHS, DualGraph, Walk, _position_of
from .TopologyResults import ElementRef

if TYPE_CHECKING:  # pragma: no cover - typing only
    import networkx

    from ..numerics.Bridge import BridgeId
    from ..numerics.FixedPoint import FixedPoint

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

#: A word in the tangle's alphabet: bridge classes in traversal order.
Word = tuple[BridgeClass, ...]


def symbol_label(index: int) -> str:
    """
    The label of the ``index``-th symbol: ``a``..``z``, then ``a1``..``z1``, ...

    Args:
        index: Position of the symbol in the alphabet, from 0.

    Returns:
        The label. Past 26 symbols a tier number is appended, so the alphabet
        never runs out and the labels stay short and sortable.

    Raises:
        ValueError: If ``index`` is negative.
    """
    if index < 0:
        raise ValueError(f"symbol index must be non-negative, got {index}")
    letter = chr(ord("a") + index % 26)
    tier = index // 26
    return letter if tier == 0 else f"{letter}{tier}"


@dataclass
class SymbolicDynamics:
    """
    The substitution a trellis induces on its bridge classes.

    Build one with :func:`symbolic_dynamics`; the fields are all derived from
    the dual graph and the bridge classes and nothing here recomputes them.

    Attributes:
        symbols: Label -> class, in symbol order (known classes first, then the
            new symbols).
        words: Class -> the word its image spells. Only the classes that have a
            bridge to walk appear, so a new symbol is absent.
        rules: Label -> the spelled-out image word (``"a": "ab"``). Labels are
            joined by :attr:`separator`.
        members: Class -> the bridge ids that realise it, as handed in.
        walks: Bridge id -> the :class:`~tanglepack.topology.DualGraph.Walk`
            that produced its word. Keeps the per-bridge diagnosis (which of the
            three cases applied, how much was registered, the route).
        new_symbols: The letters that no computed bridge realises, in symbol
            order.
        transition_graph: A ``networkx.DiGraph`` on the labels, with an edge
            ``u -> v`` weighted by how often ``v`` occurs in ``u``'s word.
        fixed_points: The trellis's fixed points, in the order that defines the
            ``fp{index}`` part of :meth:`display_label`.
        separator: What :attr:`rules` joins labels with — empty while every
            label is a single character, a space once the alphabet passes 26.
    """

    symbols: dict[str, BridgeClass]
    words: dict[BridgeClass, Word]
    rules: dict[str, str]
    members: dict[BridgeClass, list["BridgeId"]]
    walks: dict["BridgeId", Walk]
    new_symbols: list[BridgeClass]
    transition_graph: "networkx.DiGraph"
    fixed_points: list["FixedPoint"]
    separator: str = ""
    _label_of: dict[BridgeClass, str] = field(default_factory=dict)

    # ── lookups ─────────────────────────────────────────────────────────────

    @property
    def alphabet(self) -> list[str]:
        """The symbol labels in symbol order."""
        return list(self.symbols)

    @property
    def classes(self) -> list[BridgeClass]:
        """The symbols themselves, in symbol order."""
        return list(self.symbols.values())

    def label(self, bridge_class: BridgeClass) -> str:
        """
        The label standing for one class.

        Args:
            bridge_class: The class to name.

        Returns:
            Its label (``"a"``, ``"b1"``, ...).

        Raises:
            KeyError: If the class is not one of this alphabet's symbols.
        """
        try:
            return self._label_of[bridge_class]
        except KeyError:
            raise KeyError(
                f"{bridge_class.label} is not a symbol of this alphabet"
            ) from None

    def word_of(self, bridge_class: BridgeClass) -> Word:
        """
        The word one symbol rewrites to.

        Args:
            bridge_class: The class whose image word is wanted.

        Returns:
            The word, as a tuple of classes in traversal order.

        Raises:
            KeyError: If the class has no word — either it is not a symbol at
                all, or it is a NEW symbol, which by definition has no bridge to
                walk.
        """
        try:
            return self.words[bridge_class]
        except KeyError:
            known = bridge_class in self._label_of
            raise KeyError(
                f"{bridge_class.label} has no word because it is "
                + (
                    "a new symbol: no computed bridge realises it, so there is "
                    "nothing to map forward"
                    if known
                    else "not a symbol of this alphabet"
                )
            ) from None

    def display_label(
        self, item: Union[ElementRef, BridgeClass]
    ) -> str:
        """
        Render an element or a class with its fixed point's POSITION, not period.

        :attr:`~tanglepack.topology.TopologyResults.ElementRef.label` names the
        fixed point by its period, so two saddles of the same period print
        alike. This form uses the fixed point's index in :attr:`fixed_points`
        instead — ``fp0@1.0/L#2`` — which is unique within one trellis and
        stable across runs.

        Args:
            item: An :class:`~tanglepack.topology.TopologyResults.ElementRef` or
                a :class:`~tanglepack.topology.BridgeClass.BridgeClass`.

        Returns:
            The rendered name (a class as ``"<first> -> <second>"``).

        Raises:
            TypeError: If ``item`` is neither of those.
        """
        if isinstance(item, ElementRef):
            position = _position_of(item.fixed_point, self.fixed_points)
            return (
                f"fp{position}@{item.orbit_index}.{item.branch_index}"
                f"/{item.side[0].upper()}#{item.element_id}"
            )
        if isinstance(item, BridgeClass):
            return (
                f"{self.display_label(item.first)} -> "
                f"{self.display_label(item.second)}"
            )
        raise TypeError(
            f"expected an ElementRef or a BridgeClass, got {type(item).__name__}"
        )

    # ── matrices and reports ────────────────────────────────────────────────

    def transition_matrix(self) -> NDArray[np.int_]:
        """
        The substitution matrix, in symbol order.

        Returns:
            An ``(n, n)`` integer array whose ``[i, j]`` entry counts how often
            symbol ``j`` occurs in symbol ``i``'s word. Row ``i`` therefore sums
            to the length of that word, and a new symbol's row is all zeros.
        """
        order = list(self.symbols.values())
        index = {bridge_class: i for i, bridge_class in enumerate(order)}
        matrix = np.zeros((len(order), len(order)), dtype=np.int_)
        for row, bridge_class in enumerate(order):
            for letter in self.words.get(bridge_class, ()):
                matrix[row, index[letter]] += 1
        return matrix

    def case_counts(self) -> dict[str, int]:
        """
        How many bridges fell into each of the three walk cases.

        Returns:
            A dict with keys ``"i"``, ``"ii"`` and ``"iii"`` (see
            :data:`~tanglepack.topology.DualGraph.WalkCase`), always all three.
        """
        counts = {"i": 0, "ii": 0, "iii": 0}
        for walk in self.walks.values():
            counts[walk.case] += 1
        return counts

    def describe(self) -> str:
        """
        A deterministic text report: the symbol table, the rules, the matrix.

        Returns:
            The report. It contains no addresses, areas or canonical distances,
            so two builds of the same trellis produce byte-identical text.
        """
        lines: list[str] = []
        counts = self.case_counts()
        lines.append(
            f"Symbolic dynamics: {len(self.symbols)} symbols "
            f"({len(self.new_symbols)} new), {len(self.walks)} bridges "
            f"(cases i/ii/iii: {counts['i']}/{counts['ii']}/{counts['iii']})"
        )

        lines.append("")
        lines.append("Symbols")
        width = max((len(name) for name in self.symbols), default=1)
        for name, bridge_class in self.symbols.items():
            members = len(self.members.get(bridge_class, ()))
            tail = (
                f"{members} bridge{'' if members == 1 else 's'}"
                if bridge_class in self.words
                else "new symbol, no bridge realises it"
            )
            lines.append(
                f"  {name:>{width}} = {self.display_label(bridge_class)}  [{tail}]"
            )

        lines.append("")
        lines.append("Rules")
        for name in self.symbols:
            image = self.rules.get(name)
            lines.append(
                f"  {name:>{width}} -> " + (image if image else "(no rule)")
            )

        lines.append("")
        lines.append("Transition matrix (row = symbol, column = letter of its word)")
        matrix = self.transition_matrix()
        cell = max(width, max((len(str(int(v))) for v in matrix.ravel()), default=1))
        header = " ".join(f"{name:>{cell}}" for name in self.symbols)
        lines.append(f"  {'':>{width}}   {header}")
        for row, name in enumerate(self.symbols):
            body = " ".join(f"{int(value):>{cell}}" for value in matrix[row])
            lines.append(f"  {name:>{width}} | {body}")
        return "\n".join(lines)

    def __repr__(self) -> str:
        counts = self.case_counts()
        return (
            f"<SymbolicDynamics: {len(self.symbols)} symbols "
            f"({len(self.new_symbols)} new), {len(self.walks)} bridges, "
            f"cases {counts['i']}/{counts['ii']}/{counts['iii']}>"
        )


def word_of_bridge(
    dual_graph: DualGraph,
    bridge_id: "BridgeId",
    classes: Mapping[BridgeClass, Sequence["BridgeId"]],
    *,
    max_paths: int = DEFAULT_MAX_PATHS,
) -> Word:
    """
    The word one bridge's image spells.

    A thin front for
    :meth:`~tanglepack.topology.DualGraph.DualGraph.walk_bridge`, kept for
    diagnostics and tests that want the word without the route.

    Args:
        dual_graph: The dual graph to walk.
        bridge_id: The bridge whose image to read.
        classes: The trellis's bridge classes (see
            :func:`~tanglepack.topology.BridgeClass.bridge_classes`).
        max_paths: Hard cap on the walk's path enumeration.

    Returns:
        The word.

    Raises:
        AmbiguousWalkError: If the image spells more than one word.
        ValueError: See :meth:`~tanglepack.topology.DualGraph.DualGraph.walk_bridge`.
    """
    return dual_graph.walk_bridge(bridge_id, classes, max_paths=max_paths).word


def symbolic_dynamics(
    dual_graph: DualGraph,
    classes: Mapping[BridgeClass, Sequence["BridgeId"]],
    *,
    max_paths: int = DEFAULT_MAX_PATHS,
) -> SymbolicDynamics:
    """
    Read a trellis's substitution off its dual graph.

    Every bridge is walked forward one map step and its image spelled as a word
    (D.1-D.3); the words are then grouped by class (D.4). All bridges of one
    class must spell the SAME word — they are one symbol, so their images must
    be one word — and a disagreement is raised rather than resolved, because it
    means the partition is not fine enough to be an alphabet.

    Args:
        dual_graph: The dual graph of the ALL-fixed-points arrangement, carrying
            the stable partitions the classes were built from.
        classes: Class -> its bridge ids, as
            :func:`~tanglepack.topology.BridgeClass.bridge_classes` returns.
        max_paths: Hard cap on each walk's path enumeration.

    Returns:
        The :class:`SymbolicDynamics`.

    Raises:
        ValueError: If two bridges of one class spell different words (the
            message names the class, both bridges and both words), or for any of
            the walk failures
            :meth:`~tanglepack.topology.DualGraph.DualGraph.walk_bridge` raises.
        AmbiguousWalkError: If one bridge's image spells more than one word.

    Note:
        The bridge-id -> class inverse of ``classes`` is built ONCE here and
        handed to every walk, rather than rebuilt per bridge.

    Note:
        Letters that no computed bridge realises are kept as NEW SYMBOLS and
        logged at WARNING, one line per symbol: the trellis is not closed under
        the map at this extent, which is information, not a failure.
    """
    import networkx as nx

    fixed_points = list(dual_graph.trellis.fixed_points)
    walks: dict["BridgeId", Walk] = {}
    words: dict[BridgeClass, Word] = {}
    source: dict[BridgeClass, "BridgeId"] = {}
    members: dict[BridgeClass, list["BridgeId"]] = {}

    lookup = DualGraph._class_of_bridge(classes)
    for bridge_class, bridge_ids in classes.items():
        members[bridge_class] = list(bridge_ids)
        for bridge_id in bridge_ids:
            walk = dual_graph.walk_bridge(
                bridge_id,
                classes,
                max_paths=max_paths,
                class_of_bridge=lookup,
            )
            walks[bridge_id] = walk
            previous = words.get(bridge_class)
            if previous is None:
                words[bridge_class] = walk.word
                source[bridge_class] = bridge_id
            elif previous != walk.word:
                raise ValueError(
                    f"class {bridge_class.label} is realised by bridges "
                    f"{source[bridge_class]} and {bridge_id}, whose images "
                    f"spell different words: "
                    f"{_spell(previous)} vs {_spell(walk.word)}. One class is "
                    "one symbol, so the stable partition is not fine enough to "
                    "be an alphabet at this extent"
                )

    ordered = _ordered_symbols(classes, words, fixed_points)
    known = set(classes)
    new_symbols = [item for item in ordered if item not in known]
    for symbol in new_symbols:
        logger.warning(
            "new symbol %s: no computed bridge realises this class, so the "
            "trellis is not closed under the map at this extent",
            symbol.label,
        )

    labels = {
        bridge_class: symbol_label(index)
        for index, bridge_class in enumerate(ordered)
    }
    symbols = {labels[bridge_class]: bridge_class for bridge_class in ordered}
    separator = "" if all(len(name) == 1 for name in symbols) else " "
    rules = {
        labels[bridge_class]: separator.join(
            labels[letter] for letter in word
        )
        for bridge_class, word in words.items()
    }

    graph = nx.DiGraph()
    for name, bridge_class in symbols.items():
        graph.add_node(name, symbol=bridge_class, is_new=bridge_class in new_symbols)
    for bridge_class, word in words.items():
        counts: dict[BridgeClass, int] = {}
        for letter in word:
            counts[letter] = counts.get(letter, 0) + 1
        for letter, weight in counts.items():
            graph.add_edge(labels[bridge_class], labels[letter], weight=weight)

    result = SymbolicDynamics(
        symbols=symbols,
        words=words,
        rules=rules,
        members=members,
        walks=walks,
        new_symbols=new_symbols,
        transition_graph=graph,
        fixed_points=fixed_points,
        separator=separator,
        _label_of=labels,
    )
    logger.debug("%r", result)
    return result


def _ordered_symbols(
    classes: Mapping[BridgeClass, Sequence["BridgeId"]],
    words: Mapping[BridgeClass, Word],
    fixed_points: Sequence["FixedPoint"],
) -> list[BridgeClass]:
    """Known classes in class-sort order, then the new symbols in the same order."""
    known = sorted(classes, key=lambda item: class_sort_key(item, fixed_points))
    seen = set(known)
    unknown: list[BridgeClass] = []
    for word in words.values():
        for letter in word:
            if letter in seen:
                continue
            seen.add(letter)
            unknown.append(letter)
    unknown.sort(key=lambda item: class_sort_key(item, fixed_points))
    return known + unknown


def _spell(word: Word) -> str:
    """A word as its classes' labels, for an error message."""
    return " ".join(letter.label for letter in word) or "<empty>"
