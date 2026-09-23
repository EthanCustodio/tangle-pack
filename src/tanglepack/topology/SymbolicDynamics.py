"""
Symbolic dynamics of a trellis, read off bridge classes and dual-graph walks.

Each bridge class ``{X, Y}`` (:mod:`~tanglepack.topology.BridgeClass`) is an
arc between two homotopy elements of the stable manifold. One map step sends
``X`` and ``Y`` into the iterated elements ``X'`` and ``Y'``
(:func:`~tanglepack.topology.DualWalk.land_element`), and the image of the
class is the shortest walk of the dual graph from ``X'`` to ``Y'``
(:func:`~tanglepack.topology.DualWalk.shortest_walks`). The walk records both
sides of every unified stable node it crosses, so its itinerary has even
length and splits into DISJOINT consecutive pairs of iterated elements; each
pair, read through its homotopy parents, is a bridge class traversed
forward or backward, i.e. one SYMBOL of the class's word. The words of all
classes are the rewriting rules; their transition matrix is the symbolic
dynamics.

Dev Notes:

* Directions are always shown, inert classes included (author's decision,
  2026-09-21): k=10 reads ``a -> a u^-1 a^-1``, the walk crossing the inert
  class ``u = {L_1, L_3}`` against its anchor-outward orientation. Inert
  classes are lettered ``u, v, w, x, y, z, uu, uv, ...`` by THIS layer only —
  a bijective numeral over the six letters ``u..z`` (:func:`_inert_letter`),
  skipping any letter an active class already carries (the session alphabet
  can reach ``u``). The session's ``BridgeAlphabet`` stays active-only and is
  not imported (topology never imports loom).
* Refinement occurrences come ONLY from the walked (or trellis) itineraries of
  all classes, never from registered member images: on k=10 a registered
  image endpoint sits exactly on a closed cut boundary and names
  ``(R_1^2, R_3^2)``, a superscript-level mismatch that is reported at INFO
  through :attr:`SymbolicDynamics.unmatched_members` and nothing else. A
  class is split when its symbol occurs with two or more distinct iterated
  pairs; every refined child inherits the parent's word rewritten in refined
  symbols (``a_1 -> a_1 u^-1 a_2^-1`` AND ``a_2 -> a_1 u^-1 a_2^-1``). Members
  are matched to refined children only for classes that actually split.
* The singleton path: a landing that is a singleton element owns no dual-graph
  node, so such a class walks the REGULAR trellis instead — the image chain
  (:func:`~tanglepack.topology.BridgeClass._image_chain`, iterate table only)
  of its first non-loop member that has one, read by
  :func:`~tanglepack.topology.DualWalk.trellis_itinerary`. No registered chain
  → unresolved ("grow or blast").
* Registered member images are EVIDENCE, the walk is the answer: every
  non-loop member with a registered image chain is translated with the same
  function and compared at word level (``verified``) and at exact-itinerary
  level (``itinerary_verified``, INFO only). A length mismatch between the
  registered and the walked word is a WARNING — the shortest walk can be
  topologically wrong where no image is registered.
* A pair whose class is in no table entry is a VIRTUAL class (an arc no bridge
  spans yet); it is named ``new1, new2, ...`` in first-seen order with a
  WARNING and enters the graph as a sink. A pair whose two elements sit on
  different sides is flagged on the symbol (``cross_side``) with a WARNING and
  not asserted: merged outer faces or several stable branches could produce
  one.
* Ambiguous classes (several shortest walks with different itineraries) feed
  the first candidate, in the walk search's deterministic order, into the
  refinement and the graph; :attr:`SymbolicDynamics.is_reliable` is False.
* Inert and virtual symbols are sinks of the transition graph (an inert class
  with a non-empty word is reported at WARNING when the graph is built) and
  their matrix rows are zero, so the graph and the matrix agree.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field, replace
from typing import Iterable, Literal, Optional, Sequence, TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from .BridgeClass import (
    BridgeClass,
    BridgeClassEntry,
    BridgeClassTable,
    _image_chain,
    element_sort_key,
    oriented_class,
)
from .DualWalk import ElementLanding, WalkSearch, land_element, shortest_walks, trellis_itinerary
from .ElementNaming import ElementNaming
from .StablePartition import row_of_end
from .TopologyResults import ElementRef

if TYPE_CHECKING:
    import networkx

    from ..numerics.Bridge import BridgeId
    from ..numerics.FixedPoint import FixedPoint
    from .DualGraph import DualGraph
    from .PartitionFamily import PartitionFamily
    from .Trellis import Trellis

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

SymbolKind = Literal["active", "inert", "virtual"]
ItinerarySource = Literal["walk", "trellis"]

#: The digits of the inert letter series.
_INERT_DIGITS = "uvwxyz"
#: The prefix of virtual class names.
_VIRTUAL_PREFIX = "new"


# ── symbols ─────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class Symbol:
    """
    One token of a word: a bridge class traversed in one direction.

    Attributes:
        bridge_class: The class the token names.
        letter: The class's letter (session letter for an active class, inert
            letter for an inert one, ``new<n>`` for a virtual one).
        direction: ``+1`` for ``source -> target`` (``a``), ``-1`` for
            ``target -> source`` (``a^-1``).
        kind: ``"active"``, ``"inert"`` or ``"virtual"`` (no table entry).
        occurrence: The two ITERATED elements the token's pair consists of,
            class-source end first (swapped for direction ``-1``).
        refined_index: The 1-based index of the refined child this
            occurrence belongs to, ``None`` when the class is unrefined or the
            token has not been refined yet.
        cross_side: True when the pair's two elements sit on different sides
            of the stable manifold — flagged, never asserted.
    """

    bridge_class: BridgeClass
    letter: str
    direction: int
    kind: SymbolKind
    occurrence: tuple[ElementRef, ElementRef]
    refined_index: Optional[int] = None
    cross_side: bool = False

    @property
    def base(self) -> str:
        """The undirected name: ``a`` or ``a_2`` (the refined child's name)."""
        if self.refined_index is None:
            return self.letter
        return f"{self.letter}_{self.refined_index}"

    @property
    def text(self) -> str:
        """The token text: ``a``, ``a^-1``, ``a_2``, ``a_2^-1``, ``u^-1``."""
        return f"{self.base}^-1" if self.direction < 0 else self.base

    @property
    def mathtext(self) -> str:
        """The token as matplotlib mathtext, e.g. ``$a_{2}^{-1}$``."""
        body = self.letter
        if self.refined_index is not None:
            body += f"_{{{self.refined_index}}}"
        if self.direction < 0:
            body += "^{-1}"
        return f"${body}$"

    @property
    def is_inverse(self) -> bool:
        """True for a ``target -> source`` traversal."""
        return self.direction < 0

    def __str__(self) -> str:
        return self.text


@dataclass
class ImageEvidence:
    """
    The registered forward image of one member bridge, as a word.

    Attributes:
        bridge_id: The member.
        direction: Its direction in its class.
        itinerary: The iterated elements along its image chain, in the same
            even-length, disjoint-pair format as a dual-graph walk, running
            from the image of the class source to the image of the target.
        symbols: The itinerary translated (loops dropped).
        word_agrees: True when the non-loop ``(class, direction)`` sequence
            equals the class's primary word; ``None`` when the class has no
            itinerary to compare with.
        itinerary_agrees: True when the iterated refs agree exactly.
    """

    bridge_id: "BridgeId"
    direction: int
    itinerary: tuple[ElementRef, ...]
    symbols: list[Symbol]
    word_agrees: Optional[bool] = None
    itinerary_agrees: Optional[bool] = None

    @property
    def word(self) -> str:
        """The registered word as text."""
        return " ".join(symbol.text for symbol in self.symbols)


@dataclass
class ClassDynamics:
    """
    Everything the layer found out about one bridge class.

    Attributes:
        entry: The class's table entry.
        letter: Its symbol letter.
        kind: ``"active"`` or ``"inert"``.
        landings: Where the class's source and target elements land under one
            map step.
        search: The dual-graph walk search, ``None`` on the singleton/trellis
            path or when a landing failed.
        source: How the itinerary was obtained: ``"walk"`` (dual graph),
            ``"trellis"`` (image chain of a member, singleton path) or ``None``
            when unresolved.
        itinerary: The primary itinerary (even length), or ``None``.
        symbols: The itinerary's non-loop tokens — the class's word.
        loops: The itinerary pairs that fell into one element (dropped from
            the word).
        ambiguous: True when several shortest walks disagree; the first
            candidate is the primary itinerary.
        unresolved_reason: Why no itinerary was found, else ``None``.
        evidence: The registered member images, one record each.
        verified: Word-level agreement of every piece of evidence; ``None``
            without evidence or without an itinerary.
        itinerary_verified: Exact-itinerary agreement, reported only.
    """

    entry: BridgeClassEntry
    letter: str
    kind: Literal["active", "inert"]
    landings: tuple[Optional[ElementLanding], Optional[ElementLanding]]
    search: Optional[WalkSearch] = None
    source: Optional[ItinerarySource] = None
    itinerary: Optional[tuple[ElementRef, ...]] = None
    symbols: list[Symbol] = field(default_factory=list)
    loops: list[tuple[ElementRef, ElementRef]] = field(default_factory=list)
    ambiguous: bool = False
    unresolved_reason: Optional[str] = None
    evidence: list[ImageEvidence] = field(default_factory=list)
    verified: Optional[bool] = None
    itinerary_verified: Optional[bool] = None

    @property
    def bridge_class(self) -> BridgeClass:
        """The class."""
        return self.entry.bridge_class

    @property
    def resolved(self) -> bool:
        """True when an itinerary was found."""
        return self.itinerary is not None

    @property
    def word(self) -> str:
        """The word as text (empty for a trivial loop or an unresolved class)."""
        return " ".join(symbol.text for symbol in self.symbols)

    @property
    def status(self) -> str:
        """One word: ``unresolved``, ``ambiguous``, ``trivial`` or ``resolved``."""
        if self.itinerary is None:
            return "unresolved"
        if self.ambiguous:
            return "ambiguous"
        if not self.symbols:
            return "trivial"
        return "resolved"


@dataclass
class RefinedClass:
    """
    One refined child of a bridge class: the class restricted to one
    iterated-element pair.

    Attributes:
        parent: The class refined.
        letter: The parent's letter.
        index: The 1-based child index, anchor outward.
        occurrence: The iterated pair this child stands for (source end
            first).
        name: ``letter_index``, e.g. ``a_2``.
        members: The member bridges whose own endpoint pair matches
            ``occurrence``.
    """

    parent: BridgeClass
    letter: str
    index: int
    occurrence: tuple[ElementRef, ElementRef]
    name: str
    members: list["BridgeId"] = field(default_factory=list)


@dataclass(frozen=True)
class _SymbolNode:
    """One node of the transition graph."""

    name: str
    kind: SymbolKind
    bridge_class: BridgeClass
    parent: str
    inert: bool
    unresolved: bool


# ── letters ─────────────────────────────────────────────────────────────────


def _bijective(index: int, digits: str) -> str:
    """The ``index``-th (zero-based) numeral of the bijective base ``len(digits)``."""
    if index < 0:
        raise ValueError(f"letter index must be non-negative, got {index}")
    base = len(digits)
    chars = []
    index += 1
    while index > 0:
        index, remainder = divmod(index - 1, base)
        chars.append(digits[remainder])
    return "".join(reversed(chars))


def _inert_letter(index: int) -> str:
    """The ``index``-th inert letter: ``u, v, w, x, y, z, uu, uv, ...``."""
    return _bijective(index, _INERT_DIGITS)


def _active_letter(index: int) -> str:
    """The ``index``-th active fallback letter: ``a, b, ..., z, aa, ...``."""
    return _bijective(index, "abcdefghijklmnopqrstuvwxyz")


def _letter_key(letter: str) -> tuple[int, str]:
    """Order letters ``a < b < ... < z < aa`` (and ``new1 < new2 < ... < new10``)."""
    if letter.startswith(_VIRTUAL_PREFIX) and letter[len(_VIRTUAL_PREFIX):].isdigit():
        return (int(letter[len(_VIRTUAL_PREFIX):]), letter)
    return (len(letter), letter)


def inert_letters(
    table: BridgeClassTable,
    *,
    used: Iterable[str] = (),
) -> dict[BridgeClass, str]:
    """
    Letter the inert classes of a table: ``u, v, w, x, y, z, uu, ...``.

    Args:
        table: The bridge classes; the active entries' ``letter`` s are the
            letters to avoid. The inert entries are lettered in TABLE order
            (a table built by
            :func:`~tanglepack.topology.BridgeClass.bridge_classes` is sorted
            by ``min_unstable_cdist``, the same order the active letters
            follow).
        used: Extra letters to skip (e.g. fallback letters handed to
            unlettered active classes by the caller).

    Returns:
        ``{inert class: letter}``; a letter an active class already carries is
        skipped.
    """
    taken = {entry.letter for entry in table.active if entry.letter is not None}
    taken.update(used)
    entries = list(table.inert)
    letters: dict[BridgeClass, str] = {}
    index = 0
    for entry in entries:
        candidate = _inert_letter(index)
        index += 1
        while candidate in taken:
            candidate = _inert_letter(index)
            index += 1
        taken.add(candidate)
        letters[entry.bridge_class] = candidate
    return letters


def _active_letters(table: BridgeClassTable) -> dict[BridgeClass, str]:
    """The active classes' letters, with local ``a, b, ...`` fallbacks for unlettered ones."""
    letters: dict[BridgeClass, str] = {}
    taken = {entry.letter for entry in table if entry.letter is not None}
    index = 0
    for entry in table.active:
        if entry.letter is not None:
            letters[entry.bridge_class] = entry.letter
            continue
        candidate = _active_letter(index)
        index += 1
        while candidate in taken:
            candidate = _active_letter(index)
            index += 1
        taken.add(candidate)
        letters[entry.bridge_class] = candidate
        logger.debug(
            "active class %s came unlettered; using local fallback letter %r",
            entry.bridge_class.label,
            candidate,
        )
    return letters


# ── itineraries → words ─────────────────────────────────────────────────────


def _virtual_name(letters: dict[BridgeClass, str]) -> str:
    """The next ``new<n>`` name, counting the virtual names already handed out."""
    count = sum(
        1
        for name in letters.values()
        if name.startswith(_VIRTUAL_PREFIX) and name[len(_VIRTUAL_PREFIX):].isdigit()
    )
    return f"{_VIRTUAL_PREFIX}{count + 1}"


def translate_itinerary(
    itinerary: Sequence[ElementRef],
    naming: ElementNaming,
    table: BridgeClassTable,
    letters: dict[BridgeClass, str],
    fixed_points: Sequence["FixedPoint"],
) -> tuple[list[Symbol], list[tuple[ElementRef, ElementRef]]]:
    """
    Read an even-length itinerary of iterated elements as a word.

    The itinerary splits into disjoint consecutive pairs ``(e_0, e_1),
    (e_2, e_3), ...``; each pair's homotopy parents
    (:meth:`~tanglepack.topology.ElementNaming.ElementNaming.parent_of`) name
    a bridge class and a direction through
    :func:`~tanglepack.topology.BridgeClass.oriented_class`. A pair with both
    parents equal is a loop and contributes no token.

    Args:
        itinerary: The iterated element refs, even in number.
        naming: Resolves each iterated ref to its homotopy parent.
        table: The bridge classes (decides ``active`` / ``inert`` / virtual).
        letters: ``{class: letter}`` for the classes of ``table``. EXTENDED IN
            PLACE with ``new1, new2, ...`` for virtual classes met here, so
            one dict shared over every itinerary names them consistently. A
            table class missing from it falls back to the entry's own letter,
            then to the class label.
        fixed_points: The ordering context of
            :func:`~tanglepack.topology.BridgeClass.oriented_class`.

    Returns:
        ``(symbols, loops)``: the non-loop tokens in itinerary order, and the
        loop pairs (iterated refs) dropped from the word.

    Raises:
        ValueError: If the itinerary has odd length.
    """
    itinerary = tuple(itinerary)
    if len(itinerary) % 2:
        raise ValueError(
            f"an itinerary has even length (entry, exit per crossing); got "
            f"{len(itinerary)} elements"
        )
    symbols: list[Symbol] = []
    loops: list[tuple[ElementRef, ElementRef]] = []
    for position in range(0, len(itinerary), 2):
        first, second = itinerary[position], itinerary[position + 1]
        parent_first = naming.parent_of(first)
        parent_second = naming.parent_of(second)
        bridge_class, direction = oriented_class(parent_first, parent_second, fixed_points)
        if direction == 0:
            loops.append((first, second))
            continue
        cross_side = first.side != second.side
        if cross_side:
            logger.warning(
                "itinerary pair (%s, %s) connects elements on different sides of "
                "the stable manifold; flagged, not asserted",
                first.label,
                second.label,
            )
        if bridge_class in table:
            entry = table[bridge_class]
            kind: SymbolKind = "active" if entry.active else "inert"
            letter = letters.get(bridge_class)
            if letter is None:
                letter = entry.letter if entry.letter is not None else bridge_class.label
        else:
            kind = "virtual"
            letter = letters.get(bridge_class)
            if letter is None:
                letter = _virtual_name(letters)
                letters[bridge_class] = letter
                logger.warning(
                    "itinerary pair (%s, %s) names class %s, which no bridge spans; "
                    "naming it %r (virtual class)",
                    first.label,
                    second.label,
                    bridge_class.label,
                    letter,
                )
        occurrence = (first, second) if direction > 0 else (second, first)
        symbols.append(
            Symbol(bridge_class, letter, direction, kind, occurrence, None, cross_side)
        )
    return symbols, loops


def _signature(symbols: Iterable[Symbol]) -> tuple[tuple[BridgeClass, int], ...]:
    """The ``(class, direction)`` sequence of a word — what word-level agreement compares."""
    return tuple((symbol.bridge_class, symbol.direction) for symbol in symbols)


# ── refinement ──────────────────────────────────────────────────────────────


def refine(
    words: dict[BridgeClass, list[Symbol]],
    fixed_points: Sequence["FixedPoint"],
) -> dict[BridgeClass, list[RefinedClass]]:
    """
    Split every class that occurs with more than one iterated-element pair.

    Args:
        words: The primary words of every resolved class, ``{class: tokens}``.
            The occurrences counted are those of the TOKENS (a class's own
            key contributes nothing by itself).
        fixed_points: The ordering context of
            :func:`~tanglepack.topology.BridgeClass.element_sort_key`; children
            are numbered anchor outward by their occurrence's source then
            target.

    Returns:
        ``{class: [child_1, child_2, ...]}`` for the classes with two or more
        distinct occurrences; a class with a single occurrence is absent
        (unrefined). The children carry no members yet.
    """
    occurrences: dict[BridgeClass, dict[tuple[ElementRef, ElementRef], str]] = {}
    for symbols in words.values():
        for symbol in symbols:
            occurrences.setdefault(symbol.bridge_class, {}).setdefault(
                symbol.occurrence, symbol.letter
            )
    refined: dict[BridgeClass, list[RefinedClass]] = {}
    for bridge_class, seen in occurrences.items():
        if len(seen) < 2:
            continue
        ordered = sorted(
            seen.items(),
            key=lambda item: element_sort_key(item[0][0], fixed_points)
            + element_sort_key(item[0][1], fixed_points),
        )
        letter = ordered[0][1]
        refined[bridge_class] = [
            RefinedClass(bridge_class, letter, index, occurrence, f"{letter}_{index}")
            for index, (occurrence, _letter) in enumerate(ordered, start=1)
        ]
        logger.debug(
            "class %s (%s) refines into %d symbols",
            bridge_class.label,
            letter,
            len(ordered),
        )
    return refined


def _refined_symbol(
    symbol: Symbol, refined: dict[BridgeClass, list[RefinedClass]]
) -> Symbol:
    """The token with its refined index attached (unchanged for an unrefined class)."""
    children = refined.get(symbol.bridge_class)
    if not children:
        return symbol
    for child in children:
        if child.occurrence == symbol.occurrence:
            return replace(symbol, refined_index=child.index)
    return symbol


# ── registered evidence ─────────────────────────────────────────────────────


def registered_evidence(
    trellis: "Trellis",
    iterated: "PartitionFamily",
    naming: ElementNaming,
    table: BridgeClassTable,
    letters: dict[BridgeClass, str],
    entry: BridgeClassEntry,
) -> list[ImageEvidence]:
    """
    The registered forward images of a class's members, as words.

    Args:
        trellis: The full trellis (iterate table, crossings).
        iterated: The iterated partition family owning the image crossings.
        naming: Resolves iterated refs to their homotopy parents.
        table: The bridge classes.
        letters: The letter map (see :func:`translate_itinerary`).
        entry: The class whose members to read.

    Returns:
        One :class:`ImageEvidence` per non-loop member with a registered image
        chain, in member order, with ``word_agrees`` / ``itinerary_agrees``
        left ``None`` (the driver fills them). Loop members and members without
        a registered chain contribute nothing.
    """
    evidence: list[ImageEvidence] = []
    for member in entry.members:
        if member.is_loop:
            continue
        chain = _image_chain(trellis, member.bridge_id)
        if not chain:
            continue
        try:
            itinerary = tuple(trellis_itinerary(trellis, iterated, chain, member.direction))
        except ValueError as error:
            logger.warning(
                "registered image of member %s of class %s cannot be read as an "
                "itinerary; it contributes no evidence (%s)",
                member.bridge_id,
                entry.letter if entry.letter is not None else entry.bridge_class.label,
                error,
            )
            continue
        symbols, _loops = translate_itinerary(
            itinerary, naming, table, letters, trellis.fixed_points
        )
        evidence.append(ImageEvidence(member.bridge_id, member.direction, itinerary, symbols))
    return evidence


def _member_pair(
    trellis: "Trellis", iterated: "PartitionFamily", bridge_id: "BridgeId", direction: int
) -> tuple[ElementRef, ElementRef]:
    """A member's own endpoint pair in iterated owners, class-source end first."""
    first = iterated.owner_of_intersection(bridge_id[0], row_of_end(trellis, bridge_id, "first"))
    second = iterated.owner_of_intersection(
        bridge_id[1], row_of_end(trellis, bridge_id, "second")
    )
    return (first, second) if direction > 0 else (second, first)


# ── the result ──────────────────────────────────────────────────────────────


class SymbolicDynamics:
    """
    The words of every bridge class, their refinement and transition structure.

    Built by :func:`symbolic_dynamics`; constructible by hand from finished
    :class:`ClassDynamics` records (the refinement and the rules are derived
    in the constructor, member matching by :meth:`match_members`).

    Attributes:
        naming: The element naming (homotopy and iterated families).
        table: The bridge classes.
        fixed_points: The ordering context.
        classes: ``{class: ClassDynamics}`` for every table entry, in table
            order.
        inert_letters: ``{inert class: letter}``.
        virtual_classes: ``{class: "new<n>"}`` for the classes no bridge spans
            that some itinerary named.
        refined: ``{class: children}`` for every class that split.
        member_refinement: ``{bridge id: child or None}`` for the non-loop
            members of the split classes.
        unmatched_members: ``{bridge id: its own iterated pair}`` for the
            members of split classes matching no child.
        rules: ``{letter: word}`` over unrefined symbols, one per resolved
            class.
        refined_rules: ``{name: word}`` over refined symbols; a split class
            contributes one rule per child (all with the same word), an
            unsplit one keeps its letter.
    """

    def __init__(
        self,
        naming: ElementNaming,
        table: BridgeClassTable,
        fixed_points: Sequence["FixedPoint"],
        classes: dict[BridgeClass, ClassDynamics],
        *,
        inert_letters: Optional[dict[BridgeClass, str]] = None,
        virtual_classes: Optional[dict[BridgeClass, str]] = None,
    ) -> None:
        """
        Args:
            naming: The element naming.
            table: The bridge classes.
            fixed_points: The ordering context (normally ``trellis.fixed_points``).
            classes: The per-class records, itineraries and symbols filled in.
            inert_letters: ``{inert class: letter}``; derived from ``classes``
                when omitted.
            virtual_classes: ``{class: name}`` of the virtual classes; derived
                from the symbols when omitted.
        """
        self.naming = naming
        self.table = table
        self.fixed_points = list(fixed_points)
        self.classes = dict(classes)
        self.inert_letters: dict[BridgeClass, str] = (
            dict(inert_letters)
            if inert_letters is not None
            else {cls: cd.letter for cls, cd in self.classes.items() if cd.kind == "inert"}
        )
        self.virtual_classes: dict[BridgeClass, str] = (
            dict(virtual_classes) if virtual_classes is not None else {}
        )
        if virtual_classes is None:
            for cd in self.classes.values():
                for symbol in cd.symbols:
                    if symbol.kind == "virtual":
                        self.virtual_classes.setdefault(symbol.bridge_class, symbol.letter)
        self.member_refinement: dict["BridgeId", Optional[RefinedClass]] = {}
        self.unmatched_members: dict["BridgeId", tuple[ElementRef, ElementRef]] = {}

        words = {
            cls: cd.symbols for cls, cd in self.classes.items() if cd.itinerary is not None
        }
        self.refined: dict[BridgeClass, list[RefinedClass]] = refine(words, self.fixed_points)
        self.rules: dict[str, list[Symbol]] = {}
        self.refined_rules: dict[str, list[Symbol]] = {}
        for cls, cd in self.classes.items():
            if cd.itinerary is None:
                continue
            self.rules[cd.letter] = list(cd.symbols)
            refined_word = [_refined_symbol(symbol, self.refined) for symbol in cd.symbols]
            children = self.refined.get(cls)
            if children:
                for child in children:
                    self.refined_rules[child.name] = list(refined_word)
            else:
                self.refined_rules[cd.letter] = refined_word

    # ── members ─────────────────────────────────────────────────────────

    def match_members(self, trellis: "Trellis", iterated: "PartitionFamily") -> None:
        """
        Attach each split class's non-loop members to the child they belong to.

        A member's own endpoint pair — the iterated owners of its two
        crossings on their rows, swapped for direction ``-1`` — is compared
        with the children's occurrences. No match leaves the member in
        :attr:`unmatched_members` (INFO). Under the footprint cut rule the
        fixtures match every member; an endpoint on an EXISTING hole
        boundary, which the cut never re-opens, can still be owned by the
        neighbour element.

        Args:
            trellis: The full trellis.
            iterated: The iterated partition family.
        """
        for cls, children in self.refined.items():
            cd = self.classes.get(cls)
            if cd is None:
                continue
            for member in cd.entry.members:
                if member.is_loop:
                    continue
                try:
                    pair = _member_pair(trellis, iterated, member.bridge_id, member.direction)
                except ValueError as error:
                    logger.info(
                        "member %s of class %s has no iterated owner pair (%s)",
                        member.bridge_id,
                        cd.letter,
                        error,
                    )
                    continue
                match = next((child for child in children if child.occurrence == pair), None)
                self.member_refinement[member.bridge_id] = match
                if match is None:
                    self.unmatched_members[member.bridge_id] = pair
                    logger.info(
                        "member %s of class %s sits on (%s, %s), which is none of its "
                        "refined occurrences %s",
                        member.bridge_id,
                        cd.letter,
                        self._name(pair[0]),
                        self._name(pair[1]),
                        [self._pair_text(child.occurrence) for child in children],
                    )
                else:
                    match.members.append(member.bridge_id)

    # ── lookups ─────────────────────────────────────────────────────────

    def dynamics_of(self, letter: str) -> ClassDynamics:
        """
        The record of the class a letter names.

        Args:
            letter: An active or inert letter.

        Returns:
            Its :class:`ClassDynamics`.

        Raises:
            KeyError: If no class carries the letter.
        """
        for cd in self.classes.values():
            if cd.letter == letter:
                return cd
        raise KeyError(f"no bridge class carries the letter {letter!r}")

    def word(self, letter: str, *, refined: bool = True) -> str:
        """
        The right-hand side of one rule as text.

        Args:
            letter: A rule key: a class letter, or a refined child's name
                (``a_2``) when ``refined``.
            refined: Read :attr:`refined_rules` (default) or :attr:`rules`.

        Returns:
            The word, e.g. ``"a_1 u^-1 a_2^-1"``; empty for a trivial loop.

        Raises:
            KeyError: If no rule has that key (an unresolved class has none).
        """
        rules = self.refined_rules if refined else self.rules
        if letter not in rules:
            raise KeyError(f"no {'refined ' if refined else ''}rule for symbol {letter!r}")
        return " ".join(symbol.text for symbol in rules[letter])

    @property
    def is_reliable(self) -> bool:
        """True when no class is ambiguous, unresolved or contradicted by evidence."""
        return all(
            cd.itinerary is not None and not cd.ambiguous and cd.verified is not False
            for cd in self.classes.values()
        )

    @property
    def unresolved(self) -> list[ClassDynamics]:
        """The classes with no itinerary."""
        return [cd for cd in self.classes.values() if cd.itinerary is None]

    # ── transition structure ────────────────────────────────────────────

    def symbol_nodes(self, *, refined: bool = True) -> list[_SymbolNode]:
        """
        The transition graph's nodes in their deterministic order.

        Active symbols first (letter order, a split class's children in
        index order), then inert, then virtual.

        Args:
            refined: Use refined children (default) or the bare letters.

        Returns:
            The nodes, in row order of :meth:`transition_matrix`.
        """
        nodes: list[_SymbolNode] = []

        def add(cls: BridgeClass, letter: str, kind: SymbolKind, unresolved: bool) -> None:
            children = self.refined.get(cls) if refined else None
            if children:
                for child in children:
                    nodes.append(
                        _SymbolNode(child.name, kind, cls, letter, kind != "active", unresolved)
                    )
            else:
                nodes.append(_SymbolNode(letter, kind, cls, letter, kind != "active", unresolved))

        for kind in ("active", "inert"):
            records = sorted(
                (cd for cd in self.classes.values() if cd.kind == kind),
                key=lambda cd: _letter_key(cd.letter),
            )
            for cd in records:
                add(cd.bridge_class, cd.letter, kind, cd.itinerary is None)
        for cls, name in sorted(self.virtual_classes.items(), key=lambda item: _letter_key(item[1])):
            add(cls, name, "virtual", True)
        return nodes

    def _rules_for(self, refined: bool) -> dict[str, list[Symbol]]:
        return self.refined_rules if refined else self.rules

    def _edges(self, refined: bool) -> dict[tuple[str, str], tuple[int, int]]:
        """``{(source, target): (weight, inverse count)}`` over the sink-respecting rules."""
        nodes = {node.name: node for node in self.symbol_nodes(refined=refined)}
        edges: dict[tuple[str, str], tuple[int, int]] = {}
        for name, symbols in self._rules_for(refined).items():
            node = nodes.get(name)
            if node is None:
                continue
            if node.kind != "active":
                if symbols:
                    logger.warning(
                        "%s class %s has the non-empty word %r; it is a sink of the "
                        "transition graph regardless",
                        node.kind,
                        name,
                        " ".join(symbol.text for symbol in symbols),
                    )
                continue
            for symbol in symbols:
                target = symbol.base if refined else symbol.letter
                weight, inverse = edges.get((name, target), (0, 0))
                edges[(name, target)] = (weight + 1, inverse + (symbol.direction < 0))
        return edges

    def transition_graph(self, *, refined: bool = True) -> "networkx.DiGraph":
        """
        The transition graph: one node per symbol, one edge per token.

        Args:
            refined: Over refined symbols (default) or bare letters.

        Returns:
            A directed graph whose nodes carry ``kind``, ``bridge_class``,
            ``parent``, ``inert`` and ``unresolved``; each edge carries
            ``weight`` (how many times the target occurs in the source's word)
            and ``inverse`` (how many of those are ``^-1``). Inert and virtual
            symbols are sinks.
        """
        import networkx

        graph = networkx.DiGraph()
        for node in self.symbol_nodes(refined=refined):
            graph.add_node(
                node.name,
                kind=node.kind,
                bridge_class=node.bridge_class,
                parent=node.parent,
                inert=node.inert,
                unresolved=node.unresolved,
            )
        for (source, target), (weight, inverse) in self._edges(refined).items():
            if target not in graph:
                graph.add_node(
                    target,
                    kind="virtual",
                    bridge_class=None,
                    parent=target,
                    inert=True,
                    unresolved=True,
                )
            graph.add_edge(source, target, weight=weight, inverse=inverse)
        return graph

    def transition_matrix(self, *, refined: bool = True) -> tuple[list[str], NDArray[np.int64]]:
        """
        The transition matrix over the symbols in :meth:`symbol_nodes` order.

        Args:
            refined: Over refined symbols (default) or bare letters.

        Returns:
            ``(symbols, matrix)``: row ``i`` counts the occurrences of symbol
            ``j`` in the word of symbol ``i`` (zero rows for the inert and
            virtual sinks).
        """
        names = [node.name for node in self.symbol_nodes(refined=refined)]
        index = {name: position for position, name in enumerate(names)}
        matrix = np.zeros((len(names), len(names)), dtype=np.int64)
        for (source, target), (weight, _inverse) in self._edges(refined).items():
            if source in index and target in index:
                matrix[index[source], index[target]] = weight
        return names, matrix

    def spectral_radius(self, *, refined: bool = True) -> float:
        """
        The largest eigenvalue modulus of the transition matrix.

        Args:
            refined: Over refined symbols (default) or bare letters.

        Returns:
            The spectral radius (``0.0`` for an empty matrix).
        """
        _names, matrix = self.transition_matrix(refined=refined)
        if matrix.size == 0:
            return 0.0
        return float(np.max(np.abs(np.linalg.eigvals(matrix.astype(np.float64)))))

    # ── reporting ───────────────────────────────────────────────────────

    def _name(self, ref: Optional[ElementRef]) -> str:
        """An element's name text (its label when the naming cannot name it)."""
        if ref is None:
            return "-"
        try:
            return str(self.naming.name(ref).text)
        except Exception:  # noqa: BLE001 - a report never fails on a name
            return ref.label

    def _pair_text(self, pair: tuple[ElementRef, ElementRef]) -> str:
        return f"({self._name(pair[0])}, {self._name(pair[1])})"

    def itinerary_text(self, itinerary: Optional[Sequence[ElementRef]]) -> str:
        """
        An itinerary in element names, ``" | "`` between its disjoint pairs.

        Args:
            itinerary: The iterated refs, or ``None``.

        Returns:
            E.g. ``"R_1^1 R_3^3 | L_3 L_1^2 | R_3^1 R_1^3"``; ``"-"`` for ``None``.
        """
        if itinerary is None:
            return "-"
        pairs = [
            f"{self._name(itinerary[i])} {self._name(itinerary[i + 1])}"
            for i in range(0, len(itinerary) - 1, 2)
        ]
        if len(itinerary) % 2:
            pairs.append(self._name(itinerary[-1]))
        return " | ".join(pairs)

    def describe(self) -> str:
        """
        A multi-line report: naming, every class, the refinement, the matrix.

        Returns:
            The report text.
        """
        lines = [self.naming.describe(), ""]
        lines.append(
            f"{len(self.classes)} bridge class(es): "
            f"{sum(cd.kind == 'active' for cd in self.classes.values())} active, "
            f"{len(self.inert_letters)} inert, {len(self.virtual_classes)} virtual"
            f"{'' if self.is_reliable else ' (NOT reliable)'}"
        )
        for cd in self.classes.values():
            landing_names = ", ".join(
                self._name(landing.target) if landing is not None else "-"
                for landing in cd.landings
            )
            lines.append(
                f"  {cd.letter}: {cd.bridge_class.label} [{cd.kind}] lands in "
                f"({landing_names}); source {cd.source or '-'}; {cd.status}"
                + (f" ({cd.unresolved_reason})" if cd.unresolved_reason else "")
            )
            lines.append(f"    itinerary: {self.itinerary_text(cd.itinerary)}")
            if cd.itinerary is not None:
                refined_word = " ".join(
                    symbol.text
                    for symbol in (
                        _refined_symbol(symbol, self.refined) for symbol in cd.symbols
                    )
                )
                lines.append(f"    word: {cd.letter} -> {cd.word or '(empty)'}")
                lines.append(f"    refined: {refined_word or '(empty)'}")
                if cd.loops:
                    lines.append(
                        "    loops: " + ", ".join(self._pair_text(pair) for pair in cd.loops)
                    )
            lines.append(
                f"    verified: {cd.verified}; itinerary verified: "
                f"{cd.itinerary_verified}; evidence: {len(cd.evidence)}"
            )
        if self.virtual_classes:
            lines.append("virtual classes:")
            for cls, name in self.virtual_classes.items():
                lines.append(f"  {name}: {cls.label}")
        if self.refined:
            lines.append("refined symbols:")
            for cls, children in self.refined.items():
                for child in children:
                    lines.append(
                        f"  {child.name} = {self._pair_text(child.occurrence)} "
                        f"members {child.members}"
                    )
        if self.unmatched_members:
            lines.append("members matching no refined symbol:")
            for bridge_id, pair in self.unmatched_members.items():
                lines.append(f"  {bridge_id}: {self._pair_text(pair)}")
        names, matrix = self.transition_matrix()
        lines.append("transition matrix (rows = source symbol):")
        width = max((len(name) for name in names), default=1)
        lines.append("  " + " " * width + "  " + " ".join(f"{name:>{width}}" for name in names))
        for name, row in zip(names, matrix):
            lines.append(
                f"  {name:>{width}}  " + " ".join(f"{int(value):>{width}}" for value in row)
            )
        if names:
            lines.append(f"spectral radius: {self.spectral_radius():.6g}")
        return "\n".join(lines)

    def __repr__(self) -> str:
        return (
            f"<SymbolicDynamics classes={len(self.classes)} rules={len(self.refined_rules)} "
            f"reliable={self.is_reliable}>"
        )


# ── the driver ──────────────────────────────────────────────────────────────


def _homotopy_of(partition: "PartitionFamily") -> "PartitionFamily":
    """The homotopy family behind an iterated one (the family itself when it is not iterated)."""
    homotopy = getattr(partition, "homotopy", None)
    return homotopy if homotopy is not None else partition


def _resolve_itinerary(
    dual_graph: "DualGraph",
    iterated: "PartitionFamily",
    cd: ClassDynamics,
    max_walks: int,
) -> None:
    """Fill a class's itinerary: landing check, then the trellis or the walk path."""
    trellis = dual_graph.trellis
    source_landing, target_landing = cd.landings
    for landing in (source_landing, target_landing):
        if landing is None or landing.target is None:
            reason = landing.reason if landing is not None else "no landing"
            cd.unresolved_reason = (
                f"element {landing.source.label if landing else '?'} has no image "
                f"element ({reason})"
            )
            logger.warning("class %s unresolved: %s", cd.letter, cd.unresolved_reason)
            return
    assert source_landing is not None and target_landing is not None
    start, goal = source_landing.target, target_landing.target
    assert start is not None and goal is not None

    if source_landing.singleton or target_landing.singleton:
        for member in cd.entry.members:
            if member.is_loop:
                continue
            chain = _image_chain(trellis, member.bridge_id)
            if not chain:
                continue
            try:
                cd.itinerary = tuple(
                    trellis_itinerary(trellis, iterated, chain, member.direction)
                )
            except ValueError as error:
                cd.unresolved_reason = (
                    f"a landing is a singleton element and the image chain of "
                    f"member {member.bridge_id} cannot be read as an itinerary "
                    f"({error})"
                )
                logger.warning("class %s unresolved: %s", cd.letter, cd.unresolved_reason)
                return
            cd.source = "trellis"
            logger.debug(
                "class %s: singleton landing; itinerary read off the image chain of "
                "member %s (%d pairs)",
                cd.letter,
                member.bridge_id,
                len(chain),
            )
            return
        cd.unresolved_reason = (
            "a landing is a singleton element and no member has a registered image "
            "chain; grow or blast"
        )
        logger.warning("class %s unresolved: %s", cd.letter, cd.unresolved_reason)
        return

    search = shortest_walks(dual_graph, start, goal, max_walks=max_walks)
    cd.search = search
    walk = search.walk
    if walk is None or search.status not in ("unique", "trivial", "ambiguous"):
        cd.unresolved_reason = f"dual-graph walk {search.status}"
        logger.warning(
            "class %s unresolved: walk from %s to %s is %s",
            cd.letter,
            start.label,
            goal.label,
            search.status,
        )
        return
    cd.itinerary = tuple(walk.itinerary)
    cd.source = "walk"
    cd.ambiguous = search.status == "ambiguous"


def _compare_evidence(cd: ClassDynamics) -> None:
    """Fill the evidence agreement flags and the class's ``verified`` fields."""
    if not cd.evidence or cd.itinerary is None:
        cd.verified = None
        cd.itinerary_verified = None
        return
    signature = _signature(cd.symbols)
    for evidence in cd.evidence:
        evidence.word_agrees = _signature(evidence.symbols) == signature
        evidence.itinerary_agrees = tuple(evidence.itinerary) == tuple(cd.itinerary)
        if len(evidence.symbols) != len(cd.symbols):
            logger.warning(
                "class %s: registered image of member %s reads %r (%d symbol(s)) but "
                "the %s itinerary reads %r (%d); the shortest walk may be wrong where "
                "no image is registered",
                cd.letter,
                evidence.bridge_id,
                evidence.word,
                len(evidence.symbols),
                cd.source,
                cd.word,
                len(cd.symbols),
            )
    cd.verified = all(evidence.word_agrees for evidence in cd.evidence)
    cd.itinerary_verified = all(evidence.itinerary_agrees for evidence in cd.evidence)
    if not cd.verified:
        disagreeing = [e.bridge_id for e in cd.evidence if not e.word_agrees]
        logger.warning(
            "class %s: word %r disagrees with the registered image(s) of member(s) %s",
            cd.letter,
            cd.word,
            disagreeing,
        )
    elif not cd.itinerary_verified:
        logger.info(
            "class %s: word verified but the registered itinerary differs from the "
            "%s one at element level",
            cd.letter,
            cd.source,
        )


def symbolic_dynamics(
    dual_graph: "DualGraph",
    classes: BridgeClassTable,
    *,
    max_walks: int = 64,
) -> SymbolicDynamics:
    """
    Compute the symbolic dynamics of a trellis from its dual graph and classes.

    Per class: land both elements
    (:func:`~tanglepack.topology.DualWalk.land_element`); a failed landing
    leaves the class unresolved; a singleton landing takes the trellis path
    (image chain of the first non-loop member that has one); otherwise the
    shortest dual-graph walk between the two landings is the itinerary. The
    itinerary is translated to a word, checked against the registered member
    images, and the words of all classes drive the refinement and the rules.

    Args:
        dual_graph: The dual graph of the minimal trellis, built over the
            iterated homotopy partition (its ``partition.homotopy`` is the
            homotopy family).
        classes: The bridge classes of the same trellis, lettered by the
            session (unlettered active classes get local ``a, b, ...``
            fallbacks).
        max_walks: The cap on shortest walks enumerated per class.

    Returns:
        The :class:`SymbolicDynamics`.
    """
    trellis = dual_graph.trellis
    iterated = dual_graph.partition
    homotopy = _homotopy_of(iterated)
    naming = ElementNaming(homotopy, iterated)
    fixed_points = trellis.fixed_points

    active = _active_letters(classes)
    inert = inert_letters(classes, used=active.values())
    letters: dict[BridgeClass, str] = {**active, **inert}

    records: dict[BridgeClass, ClassDynamics] = {}
    for entry in classes:
        cls = entry.bridge_class
        kind: Literal["active", "inert"] = "active" if entry.active else "inert"
        landings = tuple(
            land_element(trellis, homotopy, iterated, ref) for ref in (cls.source, cls.target)
        )
        cd = ClassDynamics(entry, letters[cls], kind, (landings[0], landings[1]))
        _resolve_itinerary(dual_graph, iterated, cd, max_walks)
        if cd.itinerary is not None:
            cd.symbols, cd.loops = translate_itinerary(
                cd.itinerary, naming, classes, letters, fixed_points
            )
        cd.evidence = registered_evidence(trellis, iterated, naming, classes, letters, entry)
        _compare_evidence(cd)
        records[cls] = cd
        logger.debug(
            "class %s: %s, word %r, verified %s",
            cd.letter,
            cd.status,
            cd.word,
            cd.verified,
        )

    virtual = {cls: name for cls, name in letters.items() if cls not in classes}
    dynamics = SymbolicDynamics(
        naming,
        classes,
        fixed_points,
        records,
        inert_letters=inert,
        virtual_classes=virtual,
    )
    dynamics.match_members(trellis, iterated)
    logger.info(
        "symbolic dynamics: %d class(es), %d resolved, %d refined, reliable=%s",
        len(records),
        sum(cd.itinerary is not None for cd in records.values()),
        len(dynamics.refined),
        dynamics.is_reliable,
    )
    return dynamics
