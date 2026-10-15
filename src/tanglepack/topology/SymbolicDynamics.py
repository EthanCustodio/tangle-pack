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

The refined level (alignment plan A.8d). Each token also carries the pair of
elements of the iterated partition ``P1`` it lands on, its OCCURRENCE; the
occurrences of a class across all words are its refined symbols ``c_1, c_2,
...`` (:func:`refine`). **Theorem T**: every member of a class has the same
``P1``-level freely reduced itinerary (the strip two members bound with the
stable manifold holds no hole, so no ``P1`` cut separates their images). The
refined word of ``c_k`` is the class word with every token written as the
child its occurrence names; by T it is the word of every member of ``c_k``,
not an inheritance. The refined matrix ``M_1`` is then the IN-SPLITTING of the
class matrix ``M``: siblings have equal rows and the children's columns sum to
the parent's, a conjugacy with the same entropy. T and the in-splitting are
laws on the law cases and every blast series (``members_share_the_P1_itinerary``,
``refined_matrix_is_an_in_splitting``); the base-branch chords pair in stable
order since B7b (``PartitionFamily``).

Dev Notes:

* Directions are always shown, inert classes included (author's decision,
  2026-09-21): k=10 reads ``a -> a u^-1 a^-1``, the walk crossing the inert
  class ``u = {L_1, L_3}`` against its anchor-outward orientation. Letters are
  ONE bijective numeral over two disjoint digit sets (:func:`letter`, SIMP-6,
  WP6): active classes ``a..t, aa, ab, ...``, inert ones ``u..z, uu, uv, ...``,
  so the two series never collide and nothing is skipped. Inert classes are
  lettered by THIS layer only; the session's ``BridgeAlphabet`` letters the
  active ones with the same :func:`letter` (topology does not depend on loom).
* The ALPHABET is the table handed in -- normally ``MinimalTrellis.alphabet``,
  the classes of the minimal trellis's bridges (D2), lettered here once judged --
  closed under the word map (WP7, C9). Every pair of a walk is the class of one
  piece of ``f(B)``, and a word naming a class outside the alphabet adds it,
  walked in turn and recording the class whose word first named it
  (``named_by``). When some registered bridge spans the pair (the evidence
  table, ``minimal.classes``) the class is a copy of its evidence entry
  (members); when none
  does it is a class all the same (D11, author 2026-10-08: "we do not need
  computed bridges to find words, we strictly use the dual graph unless it is
  a singleton"): a member-less entry,
  walked from its elements' landings and flagged "no computed bridge yet"
  (``BridgeClassEntry.spanned`` False), informational only. The closure runs
  until no word names a new pair, which ends because element pairs are
  finite, and it is a function of the dual graph alone, never of which
  bridges a blast happened to compute. Only a member-less class landing on a
  singleton, which cannot walk, stays unresolved ("grow or blast"). The
  minimal trellis keeps the FIRST images of the hole bridges, and the image of
  an image bridge lies one layer further out (BCG-4), so this follows the
  hole-bridge transport chain. Every other bridge is image evidence only.
* Refinement occurrences come ONLY from the walked (or trellis) itineraries of
  all classes, never from registered member images. A class is split when
  its symbol occurs with two or more distinct iterated pairs; the refined
  word of every child is the word of every member, by T, written in refined
  symbols (``a_1 -> a_1 u^-1 a_2^-1`` AND ``a_2 -> a_1 u^-1 a_2^-1``). Members are
  matched to refined children only for classes that actually split, by their
  OWN endpoint owners: under the empty-stretch cut of ``PartitionFamily`` an
  image crossing is owned by the element of the stable edge the walk crosses
  there, so a child's occurrence names exactly the elements its members land
  in (k=10's ``(8, 9)`` sits on ``a_1 = (R_1^1, R_3^3)``). A member matching
  no child goes to :attr:`SymbolicDynamics.unmatched_members` at INFO and
  nothing else.
* The singleton path: a landing that is a singleton element owns no dual-graph
  node, so such a class walks the REGULAR trellis instead — the image chain
  (:func:`~tanglepack.topology.BridgeClass._image_chain`, iterate table only)
  of its first non-loop member that has one, read by
  :func:`~tanglepack.topology.DualWalk.trellis_itinerary`. No registered chain
  → unresolved ("grow or blast").
* Registered member images are EVIDENCE, the walk is the answer: every
  non-loop member with a registered image chain is translated with the same
  function and compared at word level (``verified``) and at exact-itinerary
  level (``itinerary_verified``, INFO only). A word-level disagreement carries
  an :class:`EvidenceDiagnosis`: ``"stable end"`` when two consecutive
  crossings of the registered chain on ONE stable branch share a crossing sign
  (the image passed around the end of the computed stable branch, so its
  crossings with the trimmed-away part are not registered and the pairing
  shifts: k=2.8 at 6 blasts); ``"order flip"`` when a member end and a bound
  of its own element have registered images in the REVERSE stable order (the
  map preserves the order along a branch, so the registered cdists there are
  noise: a deep lambda-lemma copy hugging the hole-bridge end it sits 1e-6
  beside, k=2.8 at 9-13 blasts and k=3.2 at 9, INT-11); else ``"walk"`` (no
  numerical explanation). A ``"walk"`` disagreement is a WARNING and counts
  against :attr:`SymbolicDynamics.is_reliable`; ``"stable end"`` evidence is
  INCOMPLETE and ``"order flip"`` evidence is NOISE, and neither gives a
  verdict (C11, WP10; C21, WP12: ``word_agrees`` None, INFO), since the walk
  is the word and the chain is only its cross-check (D3). Both words are
  freely reduced (:func:`free_reduce`), so a word is compared with a word.
  Evidence never adds a symbol: a class outside the closed table is named by
  its label (D2).
* A pair whose two elements sit on
  different sides of ONE stable branch is flagged on the symbol
  (``cross_side``) with a WARNING and not asserted. Across two branches there
  is no side relation (D4: the left of an inversion point's branch 0.0 is the
  right of 0.1), so such a pair is never flagged.
* The non-crossing rule (D13, WP9b): the walk runs between the landing SPANS
  and, for a class with a member lobe (a kept non-loop member with both ends
  on one stable branch), it is the shortest walk whose closed curve holds
  exactly the images of the markers the lobe holds
  (:func:`~tanglepack.topology.DualWalk.lobe_image_markers`). It applies to
  EVERY such class, not only to ties (C13): a unique shortest walk can break
  it (the heteroclinic ``bh``). No walk up to the detour cap holding the
  markers leaves the class unresolved (``"no_admissible_walk"``). A walk is a
  claim only when it is established (B6a): a search the budget cut short
  (reason ``"truncated"``, even with candidates), a landing span overlapping
  no stable node (``"span_meets_no_node"``) and two member lobes holding
  different markers (:class:`~tanglepack.topology.DualWalk.MembersDisagree`,
  caught before walking) all leave the class unresolved, and an unresolved
  class is no wall.
* Walls (D13's own words, "f(bridge) cannot cross the images of the unstable
  arcs"; WP9c, C17). Every class's arcs lie on known unstable branches: its
  members' (``Trellis.unstable_key_of``), or, for a class no computed bridge
  spans yet (D11), one map step past the arcs of the class that named it.
  Arcs on different unstable branches are disjoint, so are their images, so
  the decided walk (a ``resolved`` search) of a class walked EARLIER whose
  branches are disjoint from this class's is a wall: a tie loses the walks
  whose chords cross one (``shortest_walks(walls=)``). Every bridge is the
  same (D14): a class joining two stable branches or two tangles has no member
  lobe, but it has arcs, so the walls apply to it like to any other; its
  ``tangle`` follows the table's rule with the arcs' branches for the
  members'. A tie no rule decides stays ambiguous, feeds the first candidate
  in the walk search's deterministic order into the refinement and the
  graph, and :attr:`SymbolicDynamics.is_reliable` is False.
* INERT = TRANSIENT (author, Q3 I-a, 2026-10-09). Once the alphabet is
  closed, a class is inert iff no cycle is reachable from it in the graph of
  the class-level words: its word is empty or names only inert classes
  (:func:`inert_classes`). An unresolved class (no word) is active, and so is
  every class whose word names one. Only then are the classes lettered, in
  table order: active ones by the session's ``letter_for``, inert ones in the
  inert series. The minimal trellis decided before the words which hole
  bridges to map, by the author's sufficient condition (a zeroth-hole bridge
  is inert, 2026-10-14; ``MinimalTrellis``).
* Inert classes are NOT nodes of the transition graph or the matrix (author,
  2026-09-23): they contribute nothing to the dynamics, and their symbols
  are still visible in the words. A token naming an inert class is dropped
  from the edges; an inert word names only inert classes by definition.
  Every token of an active word names a class of the closed table, so every
  edge target is a node.
* A landing that straddles a cut (``contained`` False) still resolves to its
  midpoint owner, so the walk exists (D3); but under a correct partition a
  homotopy element's image lies in ONE iterated element, so a straddle is a
  construction inconsistency and counts against ``is_reliable`` (INT-8).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field, replace
from itertools import groupby
from typing import Callable, Iterable, Literal, Optional, Sequence, TYPE_CHECKING

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
from .DualWalk import (
    ElementLanding,
    MembersDisagree,
    WalkSearch,
    land_element,
    lobe_image_markers,
    shortest_walks,
    trellis_itinerary,
    walk_chords,
)
from .ElementNaming import ElementNaming
from .TopologyResults import ElementRef

if TYPE_CHECKING:
    import networkx

    from ..numerics.Bridge import BridgeId
    from ..numerics.FixedPoint import FixedPoint
    from ..numerics.Intersection import ManifoldKey
    from .DualGraph import Chord, DualGraph, FaceNode
    from .PartitionFamily import PartitionFamily
    from .Trellis import Trellis

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

SymbolKind = Literal["active", "inert"]
ItinerarySource = Literal["walk", "trellis"]

#: Active classes are lettered over ``a..t``, inert ones over ``u..z``: one
#: bijective numeral over two disjoint digit sets (SIMP-6, WP6).
_ACTIVE_DIGITS, _INERT_DIGITS = "abcdefghijklmnopqrst", "uvwxyz"


# ── symbols ─────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class Symbol:
    """
    One token of a word: a bridge class traversed in one direction.

    Attributes:
        bridge_class: The class the token names.
        letter: The class's letter (session letter for an active class, inert
            letter for an inert one; the class label for a class outside the
            table, which only registered evidence can name).
        direction: ``+1`` for ``source -> target`` (``a``), ``-1`` for
            ``target -> source`` (``a^-1``).
        kind: ``"active"`` or ``"inert"`` (``"active"`` outside the table).
        occurrence: The two ITERATED elements the token's pair consists of,
            class-source end first (swapped for direction ``-1``).
        refined_index: The 1-based index of the refined child this
            occurrence belongs to, ``None`` when the class is unrefined or the
            token has not been refined yet.
        cross_side: True when the pair's two elements sit on different sides
            of ONE stable branch (D4) — flagged, never asserted.
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

    def __str__(self) -> str:
        return self.text


@dataclass(frozen=True)
class EvidenceDiagnosis:
    """
    Why a registered image disagrees with the walked word (INT-1).

    Attributes:
        position: The first token where the two words differ.
        walked: The walked token's occurrence there (iterated pair, class
            source end first), ``None`` past the end of the walked word.
        registered: The registered token's occurrence there, ``None`` past the
            end of the registered word.
        cause: ``"stable end"`` when two consecutive crossings of the
            registered chain on one stable branch share a crossing sign: the
            image passed around the end of the computed stable branch, so its
            crossings with the trimmed-away part are not registered and the
            registered word is incomplete. ``"order flip"`` when a member
            end and a bound of its own element have registered images in the
            reverse stable order: the map preserves the order along a branch,
            so the registered cdists contradict the dynamics there and the
            evidence is noise (C21). ``"walk"`` otherwise: no numerical
            explanation, the shortest walk is not the image's class.
        crossings: The equal-sign run behind ``"stable end"``, the
            ``(end, bound)`` behind ``"order flip"``, else ``()``.
    """

    position: int
    walked: Optional[tuple[ElementRef, ElementRef]]
    registered: Optional[tuple[ElementRef, ElementRef]]
    cause: Literal["stable end", "order flip", "walk"]
    crossings: tuple[int, ...] = ()


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
        symbols: The itinerary translated (loops dropped), freely reduced.
        chain: The registered image chain the itinerary was read off.
        word_agrees: True when the non-loop ``(class, direction)`` sequence
            equals the class's primary word, False when it does not; ``None``
            when the class has no itinerary to compare with, or when the
            evidence is incomplete or noise (a ``"stable end"`` or
            ``"order flip"`` diagnosis, C11 / C21).
        itinerary_agrees: True when the iterated refs agree exactly.
        diagnosis: Why the word disagrees, when it does.
    """

    bridge_id: "BridgeId"
    direction: int
    itinerary: tuple[ElementRef, ...]
    symbols: list[Symbol]
    chain: tuple["BridgeId", ...] = ()
    word_agrees: Optional[bool] = None
    itinerary_agrees: Optional[bool] = None
    diagnosis: Optional[EvidenceDiagnosis] = None

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
        symbols: The itinerary's non-loop tokens, freely reduced — the class's word.
        loops: The itinerary pairs that fell into one element (dropped from
            the word).
        ambiguous: True when several shortest walks disagree (after the
            non-crossing rule, so only for a class with no member lobe); the
            first candidate is the primary itinerary.
        unresolved_reason: Why the class has no itinerary (a failed landing,
            no walk, a walk not established -- a truncated search, a span
            meeting no node, member lobes holding different markers (B6a) --
            or a singleton landing with no registered member image), else
            ``None``. A pair no computed bridge spans unresolves nothing
            (D11): it is a class of its own, ``entry.spanned`` False.
        evidence: The registered member images, one record each.
        verified: Word-level agreement of every piece of evidence that gives a
            verdict; ``None`` without such evidence or without an itinerary.
        itinerary_verified: Exact-itinerary agreement, reported only.
        named_by: ``None`` for a class of the alphabet handed in; for a class
            the closure added, the class whose word first named it.
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
    named_by: Optional[BridgeClass] = None

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


def letter(index: int, *, inert: bool = False) -> str:
    """
    The ``index``-th letter of the active or the inert series.

    Args:
        index: Zero-based position in the series.
        inert: The inert series ``u..z, uu, uv, ...`` instead of the active
            ``a..t, aa, ab, ...``.

    Returns:
        The bijective numeral of ``index`` over the series' digits.

    Raises:
        ValueError: If ``index`` is negative.
    """
    if index < 0:
        raise ValueError(f"letter index must be non-negative, got {index}")
    digits = _INERT_DIGITS if inert else _ACTIVE_DIGITS
    chars = []
    index += 1
    while index > 0:
        index, remainder = divmod(index - 1, len(digits))
        chars.append(digits[remainder])
    return "".join(reversed(chars))


def _letter_key(text: str) -> tuple[int, str]:
    """Order letters ``a < b < ... < z < aa``."""
    return (len(text), text)


def inert_letters(table: BridgeClassTable) -> dict[BridgeClass, str]:
    """
    Letter the inert classes of a table: ``u, v, w, x, y, z, uu, ...`` in table order.

    Args:
        table: The bridge classes (a table built by
            :func:`~tanglepack.topology.BridgeClass.bridge_classes` is in the
            order the active letters follow too).

    Returns:
        ``{inert class: letter}``.
    """
    return {entry.bridge_class: letter(i, inert=True) for i, entry in enumerate(table.inert)}


def _active_letters(
    table: BridgeClassTable, letter_for: Optional[Callable[[BridgeClass], str]]
) -> dict[BridgeClass, str]:
    """The active classes' letters: an entry's own, else ``letter_for``'s, else
    ``letter(i)`` by its table position."""
    return {
        entry.bridge_class: entry.letter
        or (letter_for(entry.bridge_class) if letter_for is not None else letter(i))
        for i, entry in enumerate(table.active)
    }


def inert_classes(words: dict[BridgeClass, Optional[list[BridgeClass]]]) -> set[BridgeClass]:
    """
    The transient classes of a word graph (Q3, I-a): no cycle is reachable from them.

    Args:
        words: ``{class: the classes its word names}``, loops dropped; ``None``
            for an unresolved class.

    Returns:
        The fixed point "inert iff the word is empty or names only inert
        classes". An unresolved class is active, and so is every class whose
        word reaches one.
    """
    inert: set[BridgeClass] = set()
    while grown := {
        cls
        for cls, word in words.items()
        if cls not in inert and word is not None and all(token in inert for token in word)
    }:
        inert |= grown
    return inert


# ── itineraries → words ─────────────────────────────────────────────────────


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
        table: The bridge classes (decides ``active`` / ``inert``). The
            driver closes it under the word map before translating, so only
            registered evidence can name a class outside it; such a token is
            named by its class label and counted active.
        letters: ``{class: letter}`` for the classes of ``table``; never
            changed. A table class missing from it falls back to the entry's
            own letter, then to the class label.
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
        cross_side = first.branch_key == second.branch_key and first.side != second.side
        if cross_side:
            logger.warning(
                "itinerary pair (%s, %s) connects elements on different sides of "
                "one stable branch; flagged, not asserted",
                first.label,
                second.label,
            )
        entry = table[bridge_class] if bridge_class in table else None
        kind: SymbolKind = "inert" if entry is not None and entry.inert else "active"
        letter = letters.get(bridge_class)
        if letter is None:
            letter = entry.letter if entry is not None and entry.letter else bridge_class.label
        occurrence = (first, second) if direction > 0 else (second, first)
        symbols.append(
            Symbol(bridge_class, letter, direction, kind, occurrence, None, cross_side)
        )
    return symbols, loops


def _signature(symbols: Iterable[Symbol]) -> tuple[tuple[BridgeClass, int], ...]:
    """The ``(class, direction)`` sequence of a word — what word-level agreement compares."""
    return tuple((symbol.bridge_class, symbol.direction) for symbol in symbols)


def free_reduce(tokens: Sequence[tuple[BridgeClass, int]]) -> list[int]:
    """
    The positions of a word surviving free reduction (A.8c).

    A token next to its own inverse cancels, ``c c^-1`` as ``c^-1 c``; one
    stack pass reaches the fixed point. The class word (``ClassDynamics.symbols``),
    the closure and inertness words and the registered evidence are all reduced
    by this one function, so a word is always compared with a word.
    :func:`translate_itinerary` stays unreduced: the T law reads the finer
    same-occurrence reduction at the ``P1`` level.

    Args:
        tokens: The ``(class, direction)`` sequence of a word, loops dropped.

    Returns:
        The indices of the tokens kept, in order.
    """
    kept: list[int] = []
    for i, (cls, direction) in enumerate(tokens):
        if kept and tokens[kept[-1]] == (cls, -direction):
            logger.debug("free reduction cancels %s^%+d against its inverse", cls.label, direction)
            kept.pop()
        else:
            kept.append(i)
    return kept


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
        word = [symbols[i] for i in free_reduce(_signature(symbols))]
        evidence.append(
            ImageEvidence(member.bridge_id, member.direction, itinerary, word, tuple(chain))
        )
    return evidence


def _member_pair(
    trellis: "Trellis", iterated: "PartitionFamily", bridge_id: "BridgeId", direction: int
) -> tuple[ElementRef, ElementRef]:
    """A member's own endpoint pair in iterated owners, class-source end first."""
    first = iterated.owner_at_end(trellis, bridge_id, "first")
    second = iterated.owner_at_end(trellis, bridge_id, "second")
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
        table: The bridge classes: the alphabet closed under the word map.
        fixed_points: The ordering context.
        classes: ``{class: ClassDynamics}`` for every table entry, in table
            order.
        inert_letters: ``{inert class: letter}``.
        refined: ``{class: children}`` for every class that split.
        member_refinement: ``{bridge id: child or None}`` for the non-loop
            members of the split classes.
        unmatched_members: ``{bridge id: its own iterated pair}`` for the
            members of split classes matching no child.
        rules: ``{letter: word}`` over unrefined symbols, one per resolved
            class.
        refined_rules: ``{name: word}`` over refined symbols; a split class
            contributes one rule per child, each the word of every member,
            by Theorem T (so all children share it); an unsplit one keeps
            its letter.
    """

    def __init__(
        self,
        naming: ElementNaming,
        table: BridgeClassTable,
        fixed_points: Sequence["FixedPoint"],
        classes: dict[BridgeClass, ClassDynamics],
        *,
        inert_letters: Optional[dict[BridgeClass, str]] = None,
    ) -> None:
        """
        Args:
            naming: The element naming.
            table: The bridge classes.
            fixed_points: The ordering context (normally ``trellis.fixed_points``).
            classes: The per-class records, itineraries and symbols filled in.
            inert_letters: ``{inert class: letter}``; derived from ``classes``
                when omitted.
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
        :attr:`unmatched_members` (INFO). Under the empty-stretch cut every
        image crossing is owned by the element of the edge the walk crosses
        there, so the fixtures match every member.

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
        """True when no class is ambiguous, unresolved, contradicted by evidence or
        landed across a cut (a straddling landing, INT-8).

        Unresolved (``itinerary is None``) covers every walk not established
        (B6a): a search the ``max_walks`` budget cut short, a landing span
        meeting no stable node, and member lobes holding different markers.
        """
        return all(
            cd.itinerary is not None
            and not cd.ambiguous
            and cd.verified is not False
            and all(landing.contained for landing in cd.landings if landing is not None)
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

        Active symbols in letter order, a split class's children in index
        order. Inert classes are not nodes (author,
        2026-09-23): they contribute nothing to the dynamics, and their
        symbols remain visible in the words.

        Args:
            refined: Use refined children (default) or the bare letters.

        Returns:
            The nodes, in row order of :meth:`transition_matrix`.
        """
        nodes: list[_SymbolNode] = []

        records = sorted(
            (cd for cd in self.classes.values() if cd.kind == "active"),
            key=lambda cd: _letter_key(cd.letter),
        )
        for cd in records:
            children = self.refined.get(cd.bridge_class) if refined else None
            names = [child.name for child in children] if children else [cd.letter]
            nodes.extend(
                _SymbolNode(name, "active", cd.bridge_class, cd.letter, False, cd.itinerary is None)
                for name in names
            )
        return nodes

    def _rules_for(self, refined: bool) -> dict[str, list[Symbol]]:
        return self.refined_rules if refined else self.rules

    def _edges(self, refined: bool) -> dict[tuple[str, str], tuple[int, int]]:
        """
        ``{(source, target): (weight, inverse count)}`` over the rules.

        Only active symbols are sources; tokens naming an inert class are
        dropped (inert classes are not in the graph).
        """
        nodes = {node.name for node in self.symbol_nodes(refined=refined)}
        edges: dict[tuple[str, str], tuple[int, int]] = {}
        for name, symbols in self._rules_for(refined).items():
            if name not in nodes:
                continue
            for symbol in symbols:
                if symbol.kind == "inert":
                    continue
                target = symbol.base if refined else symbol.letter
                weight, inverse = edges.get((name, target), (0, 0))
                edges[(name, target)] = (weight + 1, inverse + (symbol.direction < 0))
        return edges

    def transition_graph(self, *, refined: bool = True) -> "networkx.DiGraph":
        """
        The transition graph: one node per active symbol, one edge per token
        naming such a symbol.

        Inert classes are not part of the graph (their tokens in a word are
        skipped). The alphabet is closed under the word map, so every edge
        target is a node.

        Args:
            refined: Over refined symbols (default) or bare letters.

        Returns:
            A directed graph whose nodes carry ``kind``, ``bridge_class``,
            ``parent``, ``inert`` and ``unresolved``; each edge carries
            ``weight`` (how many times the target occurs in the source's word)
            and ``inverse`` (how many of those are ``^-1``).
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
            assert target in graph, f"symbol {target} of the word of {source} is no node"
            graph.add_edge(source, target, weight=weight, inverse=inverse)
        return graph

    def transition_matrix(self, *, refined: bool = True) -> tuple[list[str], NDArray[np.int64]]:
        """
        The transition matrix over the symbols in :meth:`symbol_nodes` order.

        Over refined symbols it is ``M_1``, the in-splitting of the class
        matrix ``M`` (A.8d): sibling rows are equal, and the children's
        columns sum to the parent's column of ``M``.

        Args:
            refined: Over refined symbols (default, ``M_1``) or bare letters
                (``M``).

        Returns:
            ``(symbols, matrix)``: row ``i`` counts the occurrences of symbol
            ``j`` in the word of symbol ``i`` (inert classes have no row or
            column).
        """
        names = [node.name for node in self.symbol_nodes(refined=refined)]
        index = {name: position for position, name in enumerate(names)}
        matrix = np.zeros((len(names), len(names)), dtype=np.int64)
        for (source, target), (weight, _inverse) in self._edges(refined).items():
            if source in index and target in index:
                matrix[index[source], index[target]] = weight
        return names, matrix

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
            f"{len(self.inert_letters)} inert"
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
                + ("" if cd.entry.spanned else "; no computed bridge yet")
            )
            lines.append(f"    iterated itinerary: {self.itinerary_text(cd.itinerary)}")
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
        closure = [cd for cd in self.classes.values() if cd.named_by is not None]
        if closure:
            lines.append(
                f"closure: {len(closure)} class(es) named by words "
                f"({sum(not cd.entry.spanned for cd in closure)} with no computed bridge "
                "yet): "
                + ", ".join(
                    f"{cd.letter} (by {self.classes[cd.named_by].letter})" for cd in closure
                )
            )
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


def _arc_branches(
    trellis: "Trellis", entry: BridgeClassEntry, namer_branches: frozenset["ManifoldKey"]
) -> frozenset["ManifoldKey"]:
    """The unstable branches a class's arcs lie on: its members', or, for a class no
    computed bridge spans yet (D11), one map step past the arcs of the class that named it."""
    if entry.members:
        return frozenset(trellis.unstable_key_of(member.bridge_id) for member in entry.members)
    return frozenset(key[0].advance_key(key, 1) for key in namer_branches)


def _resolve_itinerary(
    dual_graph: "DualGraph",
    iterated: "PartitionFamily",
    cd: ClassDynamics,
    max_walks: int,
    walls: Sequence[tuple["FaceNode", "Chord"]] = (),
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

    try:
        wanted = lobe_image_markers(dual_graph, cd.entry)
    except MembersDisagree as error:
        # Never walk unconstrained instead: that word would rest on fewer markers.
        cd.unresolved_reason = (
            f"members hold different markers: the class is not one homotopy class "
            f"(resolution; {error})"
        )
        logger.warning("class %s unresolved: %s", cd.letter, cd.unresolved_reason)
        return
    search = shortest_walks(
        dual_graph,
        start,
        goal,
        max_walks=max_walks,
        landings=(source_landing, target_landing),
        wanted=wanted,
        walls=walls,
    )
    cd.search = search
    if search.status == "unresolved":
        cd.unresolved_reason = (
            "search budget reached before uniqueness was established"
            if search.reason == "truncated"
            else f"dual-graph walk {search.reason}"
        )
        logger.warning(
            "class %s unresolved: walk from %s to %s: %s",
            cd.letter,
            start.label,
            goal.label,
            cd.unresolved_reason,
        )
        return
    cd.itinerary = tuple(search.walk.itinerary)
    cd.source = "walk"
    cd.ambiguous = search.status == "ambiguous"


def _diagnose(
    trellis: "Trellis",
    homotopy: "PartitionFamily",
    cls: BridgeClass,
    evidence: ImageEvidence,
    symbols: Sequence[Symbol],
) -> EvidenceDiagnosis:
    """Where a registered word parts from the walked one, and why (see :class:`EvidenceDiagnosis`)."""
    pairs = list(zip(_signature(symbols), _signature(evidence.symbols)))
    position = next((i for i, (w, r) in enumerate(pairs) if w != r), len(pairs))
    walked = symbols[position].occurrence if position < len(symbols) else None
    registered = (
        evidence.symbols[position].occurrence if position < len(evidence.symbols) else None
    )
    chain = evidence.chain
    crossings = [chain[0][0], *(pair[1] for pair in chain)] if chain else []
    # The first two or more consecutive crossings on one stable branch with one sign.
    def mark(crossing: int) -> tuple:
        ix = trellis.intersection(crossing)
        return ix.manifold_b_key, ix.crossing_sign

    runs = (tuple(run) for _, run in groupby(crossings, key=mark))
    run = next((run for run in runs if len(run) > 1), None)
    if run:
        return EvidenceDiagnosis(position, walked, registered, "stable end", run)
    flip = _order_flip(trellis, homotopy, cls, evidence)
    if flip:
        return EvidenceDiagnosis(position, walked, registered, "order flip", flip)
    return EvidenceDiagnosis(position, walked, registered, "walk")


def _order_flip(
    trellis: "Trellis", homotopy: "PartitionFamily", cls: BridgeClass, evidence: ImageEvidence
) -> Optional[tuple[int, int]]:
    """A member end and a bound of its own element whose registered images run in
    the reverse stable order (C21): the map preserves the order along a branch, so
    the registered cdists there are noise. Lookup only."""
    ends = evidence.bridge_id if evidence.direction >= 0 else evidence.bridge_id[::-1]
    for end, ref in zip(ends, (cls.source, cls.target)):
        element = homotopy.element(ref)
        for bound in (element.lo_id, element.hi_id):
            if bound is None or bound == end:
                continue
            images = trellis.iterate(end, 1), trellis.iterate(bound, 1)
            if None in images:
                continue
            x, y = (trellis.intersection(i) for i in (end, bound))
            fx, fy = (trellis.intersection(i) for i in images)
            if x.manifold_b_key == y.manifold_b_key and fx.manifold_b_key == fy.manifold_b_key and (
                (x.stable_cdist < y.stable_cdist) != (fx.stable_cdist < fy.stable_cdist)
            ):
                return end, bound
    return None


def _compare_evidence(cd: ClassDynamics, trellis: "Trellis", homotopy: "PartitionFamily") -> None:
    """Fill the evidence agreement flags, the diagnoses and the class's ``verified`` fields."""
    if not cd.evidence or cd.itinerary is None:
        cd.verified = None
        cd.itinerary_verified = None
        return
    signature = _signature(cd.symbols)
    for evidence in cd.evidence:
        evidence.itinerary_agrees = tuple(evidence.itinerary) == tuple(cd.itinerary)
        if _signature(evidence.symbols) == signature:
            evidence.word_agrees = True
            continue
        evidence.diagnosis = diagnosis = _diagnose(
            trellis, homotopy, cd.entry.bridge_class, evidence, cd.symbols
        )
        incomplete = diagnosis.cause in ("stable end", "order flip")
        evidence.word_agrees = None if incomplete else False
        if diagnosis.cause == "stable end":
            why = (
                f"incomplete: the image passes the end of the computed stable branch at "
                f"crossings {list(diagnosis.crossings)} (same branch, same sign), so it "
                "gives no verdict"
            )
        elif diagnosis.cause == "order flip":
            end, bound = diagnosis.crossings
            gaps = [
                trellis.intersection(a).stable_cdist - trellis.intersection(b).stable_cdist
                for a, b in ((end, bound), (trellis.iterate(end, 1), trellis.iterate(bound, 1)))
            ]
            why = (
                f"noise: member end {end} and its element's bound {bound} are "
                f"{gaps[0]:+.2e} apart in stable cdist but their registered images "
                f"{gaps[1]:+.2e}, an order the map cannot reverse, so it gives no verdict"
            )
        else:
            why = "no numerical explanation: the shortest walk is not the image's class"
        (logger.info if incomplete else logger.warning)(
            "class %s: the registered image of member %s reads %r but the %s itinerary "
            "reads %r; they part at token %d: %s",
            cd.letter,
            evidence.bridge_id,
            evidence.word,
            cd.source,
            cd.word,
            diagnosis.position,
            why,
        )
    verdicts = [evidence.word_agrees for evidence in cd.evidence if evidence.word_agrees is not None]
    cd.verified = all(verdicts) if verdicts else None
    cd.itinerary_verified = all(evidence.itinerary_agrees for evidence in cd.evidence)
    if cd.verified and not cd.itinerary_verified:
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
    letter_for: Optional[Callable[[BridgeClass], str]] = None,
) -> SymbolicDynamics:
    """
    Compute the symbolic dynamics of a trellis from its dual graph and classes.

    Per class: land both elements
    (:func:`~tanglepack.topology.DualWalk.land_element`); a failed landing
    leaves the class unresolved; a singleton landing takes the trellis path
    (image chain of the first non-loop member that has one); otherwise the
    shortest dual-graph walk between the two landings is the itinerary. A pair
    of the itinerary whose class is not in the alphabet yet adds that class
    (the closure under the word map, see the module Dev Notes), flagged "no
    computed bridge yet" when no registered bridge spans it (D11), and it is
    walked in turn. Once the alphabet is closed the class-level words decide
    which classes are inert (:func:`inert_classes`), the classes are lettered,
    every itinerary is translated to a word and checked against the registered
    member images, and the words of all classes drive the refinement and the
    rules.

    Args:
        dual_graph: The dual graph of the minimal trellis, built over the
            iterated homotopy partition (its ``partition.homotopy`` is the
            homotopy family; its ``minimal.classes`` is the evidence table: a
            word may only name a class it holds).
        classes: The ALPHABET: normally ``MinimalTrellis.alphabet`` (D2). Its
            entries' ``inert`` and, on the active ones, ``letter`` are set here.
        max_walks: The cap on shortest walks enumerated per class.
        letter_for: Letters an ACTIVE class that carries no letter yet
            (normally the session's ``BridgeAlphabet.letter_for``, so the
            letter persists); without it ``letter(i)`` by its position among
            the active classes. Inert classes take the inert series.

    Returns:
        The :class:`SymbolicDynamics`, its ``table`` the closed alphabet (the
        classes handed in first, then the closure in discovery order).
    """
    trellis = dual_graph.trellis
    iterated = dual_graph.partition
    homotopy = _homotopy_of(iterated)
    naming = ElementNaming(homotopy, iterated)
    fixed_points = trellis.fixed_points
    evidence_table = dual_graph.minimal.classes

    entries = list(classes)
    listed = {entry.bridge_class for entry in entries}
    closure: dict[BridgeClass, BridgeClass] = {}  # added class -> the class that named it

    def add(named: BridgeClass, named_by: BridgeClass) -> BridgeClassEntry:
        """A class a word names: a copy of its evidence entry, or a member-less
        entry when no computed bridge spans it yet (D11)."""
        if named in evidence_table:
            entry = replace(evidence_table[named], letter=None, inert=None)
        else:
            # The table's tangle rule, its arcs (one step past its namer's, on the
            # same fixed points) standing in for the members.
            owners = {id(named.source.fixed_point), id(named.target.fixed_point)}
            owners |= {id(key[0]) for key in branches[named_by]}
            tangle = next(
                (i for i, fp in enumerate(fixed_points) if {id(fp)} == owners), None
            )
            entry = BridgeClassEntry(named, tangle=tangle)
        listed.add(named)
        closure[named] = named_by
        return entry

    records: dict[BridgeClass, ClassDynamics] = {}
    words: dict[BridgeClass, Optional[list[BridgeClass]]] = {}
    # D13 / D14: images of disjoint arcs are disjoint, so the decided walk of a class
    # whose arcs lie on other unstable branches is a wall for every later class.
    branches: dict[BridgeClass, frozenset] = {}
    decided: list[tuple[frozenset, list]] = []
    for entry in entries:  # the list grows while the words name new classes
        cls = entry.bridge_class
        namer = branches.get(closure.get(cls), frozenset())
        branches[cls] = mine = _arc_branches(trellis, entry, namer)
        landings = tuple(
            land_element(trellis, homotopy, iterated, ref) for ref in (cls.source, cls.target)
        )
        # Named by its label until the verdicts and letters are in.
        cd = ClassDynamics(
            entry, cls.label, "active", (landings[0], landings[1]), named_by=closure.get(cls)
        )
        walls = [wall for theirs, chords in decided if not theirs & mine for wall in chords]
        _resolve_itinerary(dual_graph, iterated, cd, max_walks, walls)
        # A walk standing for several face paths names none of them, and an unresolved
        # class (a truncated search, B6a) established none: no wall.
        search = cd.search
        if search and search.status == "resolved" and search.walk.multiplicity == 1:
            decided.append((mine, walk_chords(search.walk)))
        itinerary = cd.itinerary or ()
        named = [
            oriented_class(naming.parent_of(x), naming.parent_of(y), fixed_points)
            for x, y in zip(itinerary[::2], itinerary[1::2])
        ]
        tokens = [(pair, direction) for pair, direction in named if direction]
        word = [tokens[i][0] for i in free_reduce(tokens)]
        words[cls] = None if cd.itinerary is None else word
        new = {pair: None for pair in word if pair not in listed}
        entries.extend(add(pair, cls) for pair in new)
        records[cls] = cd

    table = BridgeClassTable(entries)
    inert_set = inert_classes(words)
    for entry in table:
        entry.inert = entry.bridge_class in inert_set
    inert = inert_letters(table)
    letters = {**_active_letters(table, letter_for), **inert}
    assert len(set(letters.values())) == len(letters), "a letter names two classes"
    for entry in table.active:
        entry.letter = letters[entry.bridge_class]
    for cls, cd in records.items():
        cd.letter, cd.kind = letters[cls], "inert" if cls in inert_set else "active"
    if closure:
        bare = [letters[named] for named in closure if named not in evidence_table]
        logger.info(
            "the words name %d class(es) outside the alphabet; closed it under the word "
            "map (%s), %d with no computed bridge yet (%s)",
            len(closure),
            ", ".join(letters[named] for named in closure),
            len(bare),
            ", ".join(bare) or "-",
        )

    for cd in records.values():
        if cd.itinerary is not None:
            symbols, cd.loops = translate_itinerary(
                cd.itinerary, naming, table, letters, fixed_points
            )
            cd.symbols = [symbols[i] for i in free_reduce(_signature(symbols))]
        cd.evidence = registered_evidence(trellis, iterated, naming, table, letters, cd.entry)
        _compare_evidence(cd, trellis, homotopy)
        logger.debug(
            "class %s: %s, word %r, verified %s",
            cd.letter,
            cd.status,
            cd.word,
            cd.verified,
        )

    dynamics = SymbolicDynamics(naming, table, fixed_points, records, inert_letters=inert)
    dynamics.match_members(trellis, iterated)
    logger.info(
        "symbolic dynamics: %d class(es), %d resolved, %d refined, reliable=%s",
        len(records),
        sum(cd.itinerary is not None for cd in records.values()),
        len(dynamics.refined),
        dynamics.is_reliable,
    )
    return dynamics
