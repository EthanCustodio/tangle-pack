"""Letter-free spellers for bridge classes, words, partition rows and matrices.

Letters (``a``, ``b``, ``u``, ...) are assigned by table order and by a
persistent session alphabet; they are presentation, not structure. A test
that pins a fixture fact spells it in ELEMENT NAMES instead (author decision 2,
2026-10-05): a class is its oriented homotopy element pair (source anchor
outward first), a word is the sequence of those pairs with each inverse token
reversed, a refined child is its oriented iterated pair.

Every speller takes ``short``: ``False`` (default) uses the full
``ElementName.text`` (``R_(0.0;1)^2``, with the fixed-point letter when there
are several fixed points), ``True`` the short form (``R_1^2``) that
CLAUDE.md's single-branch fixture facts are written in.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:  # pragma: no cover - typing only
    from tanglepack.topology.BridgeClass import BridgeClass
    from tanglepack.topology.ElementNaming import ElementName, ElementNaming
    from tanglepack.topology.SymbolicDynamics import ClassDynamics, Symbol, SymbolicDynamics
    from tanglepack.topology.TopologyResults import ElementRef, StablePartitionResult

#: A spelled pair of element names, ``(source, target)`` in traversal order.
NamePair = tuple[str, str]


def _text(name: "ElementName", short: bool) -> str:
    """The name's short or full text."""
    return str(name.short_text if short else name.text)


def homotopy_name(naming: "ElementNaming", ref: "ElementRef", *, short: bool = False) -> str:
    """
    The text of one HOMOTOPY element's name.

    Args:
        naming: The element naming.
        ref: A homotopy-family ref.
        short: Drop the branch code (and fixed-point letter).

    Returns:
        E.g. ``"R_1"`` (short) or ``"R_(0.0;1)"``.
    """
    return _text(naming.homotopy_name(ref), short)


def iterated_name(naming: "ElementNaming", ref: "ElementRef", *, short: bool = False) -> str:
    """
    The text of one ITERATED element's name.

    Args:
        naming: The element naming.
        ref: An iterated-family ref.
        short: Drop the branch code (and fixed-point letter).

    Returns:
        E.g. ``"R_1^2"`` (short) or ``"R_(0.0;1)^2"``.
    """
    return _text(naming.name(ref), short)


def class_pair(
    naming: "ElementNaming", bridge_class: "BridgeClass", *, short: bool = False
) -> NamePair:
    """
    A bridge class as its oriented homotopy pair ``(source, target)``.

    Args:
        naming: The element naming.
        bridge_class: The class (``source`` is the anchor-nearer element).
        short: Use the short names.

    Returns:
        The two names, source first.
    """
    return (
        homotopy_name(naming, bridge_class.source, short=short),
        homotopy_name(naming, bridge_class.target, short=short),
    )


def homotopy_pair(
    dyn: "SymbolicDynamics", cd: "ClassDynamics", *, short: bool = False
) -> NamePair:
    """
    One class's oriented homotopy pair, read off its symbolic-dynamics record.

    Args:
        dyn: The symbolic dynamics.
        cd: The class record.
        short: Use the short names.

    Returns:
        ``(source, target)`` in homotopy names.
    """
    return class_pair(dyn.naming, cd.bridge_class, short=short)


def symbol_pair(naming: "ElementNaming", symbol: "Symbol", *, short: bool = False) -> NamePair:
    """
    One word token as its oriented homotopy pair, inverse tokens reversed.

    Args:
        naming: The element naming.
        symbol: The token.
        short: Use the short names.

    Returns:
        ``(source, target)`` for direction ``+1``, ``(target, source)`` for ``-1``.
    """
    pair = class_pair(naming, symbol.bridge_class, short=short)
    return pair if symbol.direction > 0 else (pair[1], pair[0])


def word_in_names(
    dyn: "SymbolicDynamics", cd: "ClassDynamics", *, short: bool = False
) -> list[NamePair]:
    """
    One class's image word as oriented homotopy pairs (loops dropped, as in the word).

    Args:
        dyn: The symbolic dynamics.
        cd: The class record.
        short: Use the short names.

    Returns:
        One pair per token, in word order.
    """
    return [symbol_pair(dyn.naming, symbol, short=short) for symbol in cd.symbols]


def occurrence_in_names(
    naming: "ElementNaming", occurrence: tuple["ElementRef", "ElementRef"], *, short: bool = False
) -> NamePair:
    """
    An iterated element pair (a refined occurrence) in iterated names.

    Args:
        naming: The element naming.
        occurrence: Two iterated refs.
        short: Use the short names.

    Returns:
        The two names in the given order.
    """
    return (
        iterated_name(naming, occurrence[0], short=short),
        iterated_name(naming, occurrence[1], short=short),
    )


def refined_children(
    dyn: "SymbolicDynamics", cd: "ClassDynamics", *, short: bool = False
) -> list[NamePair]:
    """
    A class's refined children as their oriented iterated pairs, anchor outward.

    Args:
        dyn: The symbolic dynamics.
        cd: The class record.
        short: Use the short names.

    Returns:
        One pair per child; empty when the class did not split.
    """
    return [
        occurrence_in_names(dyn.naming, child.occurrence, short=short)
        for child in dyn.refined.get(cd.bridge_class, [])
    ]


def refined_symbol_pair(
    dyn: "SymbolicDynamics", symbol: "Symbol", *, short: bool = False
) -> NamePair:
    """
    One REFINED word token as an oriented pair, inverse tokens reversed.

    A token of a split class is spelled by its refined child's iterated pair
    (the child's ``occurrence``, source end first); a token of an unsplit class
    by its homotopy pair, exactly as :func:`symbol_pair`.

    Args:
        dyn: The symbolic dynamics.
        symbol: A token of :attr:`SymbolicDynamics.refined_rules`.
        short: Use the short names.

    Returns:
        ``(source, target)`` for direction ``+1``, ``(target, source)`` for ``-1``.
    """
    if symbol.refined_index is None:
        return symbol_pair(dyn.naming, symbol, short=short)
    child = dyn.refined[symbol.bridge_class][symbol.refined_index - 1]
    pair = occurrence_in_names(dyn.naming, child.occurrence, short=short)
    return pair if symbol.direction > 0 else (pair[1], pair[0])


def refined_words_in_names(
    dyn: "SymbolicDynamics", cd: "ClassDynamics", *, short: bool = False
) -> dict[NamePair, list[NamePair]]:
    """
    A class's refined rule(s), keyed and spelled letter-free.

    A split class has one rule per refined child, keyed by the child's
    oriented iterated pair; an unsplit class has one rule keyed by its
    homotopy pair. Each word token is spelled by :func:`refined_symbol_pair`.

    Args:
        dyn: The symbolic dynamics.
        cd: The class record (resolved; an unresolved class has no rule).
        short: Use the short names.

    Returns:
        ``{rule key: [token pair, ...]}``; empty when the class has no rule.
    """
    children = dyn.refined.get(cd.bridge_class, [])
    if children:
        keyed = [
            (occurrence_in_names(dyn.naming, child.occurrence, short=short), child.name)
            for child in children
        ]
    else:
        keyed = [(homotopy_pair(dyn, cd, short=short), cd.letter)]
    return {
        key: [refined_symbol_pair(dyn, symbol, short=short) for symbol in dyn.refined_rules[name]]
        for key, name in keyed
        if name in dyn.refined_rules
    }


def itinerary_pairs(
    dyn: "SymbolicDynamics", cd: "ClassDynamics", *, short: bool = False
) -> Optional[list[NamePair]]:
    """
    A class's itinerary as its disjoint consecutive pairs of iterated names.

    Args:
        dyn: The symbolic dynamics.
        cd: The class record.
        short: Use the short names.

    Returns:
        The pairs, or None when the class has no itinerary.

    Raises:
        AssertionError: The itinerary has odd length.
    """
    if cd.itinerary is None:
        return None
    refs = list(cd.itinerary)
    assert len(refs) % 2 == 0, f"odd itinerary of length {len(refs)}"
    return [
        occurrence_in_names(dyn.naming, (refs[i], refs[i + 1]), short=short)
        for i in range(0, len(refs), 2)
    ]


def brackets(
    naming: "ElementNaming",
    result: "StablePartitionResult",
    *,
    iterated: bool = True,
    short: bool = False,
) -> list[tuple[str, str]]:
    """
    One partition row, anchor outward, as ``(name, closedness)``.

    Args:
        naming: The element naming.
        result: One ``(branch, side)`` partition result of the family named.
        iterated: Name the elements as iterated (default) or homotopy elements.
        short: Use the short names.

    Returns:
        ``[(name, "[ ]" | "( )" | "[ )" | "( ]"), ...]``.
    """
    from tanglepack.topology.TopologyResults import ElementRef

    row: list[tuple[str, str]] = []
    for interval in result.intervals:
        ref = ElementRef(result.branch_key, result.side, interval.element_id)
        name = (
            iterated_name(naming, ref, short=short)
            if iterated
            else homotopy_name(naming, ref, short=short)
        )
        shape = ("[" if interval.closed_lo else "(") + " " + ("]" if interval.closed_hi else ")")
        row.append((name, shape))
    return row


def _symbol_key(
    dyn: "SymbolicDynamics", symbol_name: str, refined: bool, short: bool
) -> NamePair:
    """The letter-free key of one transition-matrix symbol name."""
    if refined:
        for children in dyn.refined.values():
            for child in children:
                if child.name == symbol_name:
                    return occurrence_in_names(dyn.naming, child.occurrence, short=short)
    for cls, cd in dyn.classes.items():
        if cd.letter == symbol_name:
            return class_pair(dyn.naming, cls, short=short)
    for cls, letter in dyn.virtual_classes.items():
        if letter == symbol_name:
            return class_pair(dyn.naming, cls, short=short)
    raise KeyError(f"no class or refined child named {symbol_name!r}")


def matrix_by_classes(
    dyn: "SymbolicDynamics", *, refined: bool = True, short: bool = False
) -> dict[NamePair, dict[NamePair, int]]:
    """
    The transition matrix keyed by element pairs instead of letters.

    A row or column of an unrefined class is keyed by its oriented homotopy
    pair; a refined child's by its oriented iterated pair.

    Args:
        dyn: The symbolic dynamics.
        refined: Over refined symbols (default) or bare classes.
        short: Use the short names.

    Returns:
        ``{row key: {column key: count}}`` with every row and column present.
    """
    names, matrix = dyn.transition_matrix(refined=refined)
    keys = [_symbol_key(dyn, name, refined, short) for name in names]
    return {
        row_key: {col_key: int(matrix[i, j]) for j, col_key in enumerate(keys)}
        for i, row_key in enumerate(keys)
    }


def class_by_pair(
    dyn: "SymbolicDynamics", pair: NamePair, *, short: bool = False
) -> "ClassDynamics":
    """
    The one class whose oriented homotopy pair is ``pair``.

    Args:
        dyn: The symbolic dynamics.
        pair: ``(source, target)`` in homotopy names, anchor outward.
        short: ``pair`` is spelled in short names.

    Returns:
        The class record.

    Raises:
        AssertionError: No class, or more than one, carries that pair.
    """
    matches = [cd for cd in dyn.classes.values() if homotopy_pair(dyn, cd, short=short) == pair]
    assert len(matches) == 1, f"expected one class {pair}, found {len(matches)}"
    return matches[0]
