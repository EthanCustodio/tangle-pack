"""Letter-free spellers for bridge classes, words, itineraries and matrices.

Letters (``a``, ``b``, ``u``, ...) are assigned by table order and by a
persistent session alphabet; they are presentation, not structure. When a test
compares two builds of one case (a session product against a direct build),
it spells both in ELEMENT NAMES instead: a class is its oriented homotopy
element pair (source anchor outward first), a word is the sequence of those
pairs with each inverse token reversed, a refined child is its oriented
iterated pair. Names are the full ``ElementName.text`` (``R_(0.0;1)^2``, with
the fixed-point letter when there are several fixed points).

These spellers never pin a case-specific fact (author, 2026-10-05): they only
make two results of the same case comparable.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:  # pragma: no cover - typing only
    from tanglepack.topology.BridgeClass import BridgeClass
    from tanglepack.topology.ElementNaming import ElementNaming
    from tanglepack.topology.SymbolicDynamics import ClassDynamics, Symbol, SymbolicDynamics
    from tanglepack.topology.TopologyResults import ElementRef

#: A spelled pair of element names, ``(source, target)`` in traversal order.
NamePair = tuple[str, str]


def homotopy_name(naming: "ElementNaming", ref: "ElementRef") -> str:
    """
    The text of one HOMOTOPY element's name.

    Args:
        naming: The element naming.
        ref: A homotopy-family ref.

    Returns:
        E.g. ``"R_(0.0;1)"``.
    """
    return str(naming.homotopy_name(ref).text)


def iterated_name(naming: "ElementNaming", ref: "ElementRef") -> str:
    """
    The text of one ITERATED element's name.

    Args:
        naming: The element naming.
        ref: An iterated-family ref.

    Returns:
        E.g. ``"R_(0.0;1)^2"``.
    """
    return str(naming.name(ref).text)


def class_pair(naming: "ElementNaming", bridge_class: "BridgeClass") -> NamePair:
    """
    A bridge class as its oriented homotopy pair ``(source, target)``.

    Args:
        naming: The element naming.
        bridge_class: The class (``source`` is the anchor-nearer element).

    Returns:
        The two names, source first.
    """
    return (
        homotopy_name(naming, bridge_class.source),
        homotopy_name(naming, bridge_class.target),
    )


def homotopy_pair(dyn: "SymbolicDynamics", cd: "ClassDynamics") -> NamePair:
    """
    One class's oriented homotopy pair, read off its symbolic-dynamics record.

    Args:
        dyn: The symbolic dynamics.
        cd: The class record.

    Returns:
        ``(source, target)`` in homotopy names.
    """
    return class_pair(dyn.naming, cd.bridge_class)


def symbol_pair(naming: "ElementNaming", symbol: "Symbol") -> NamePair:
    """
    One word token as its oriented homotopy pair, inverse tokens reversed.

    Args:
        naming: The element naming.
        symbol: The token.

    Returns:
        ``(source, target)`` for direction ``+1``, ``(target, source)`` for ``-1``.
    """
    pair = class_pair(naming, symbol.bridge_class)
    return pair if symbol.direction > 0 else (pair[1], pair[0])


def word_in_names(dyn: "SymbolicDynamics", cd: "ClassDynamics") -> list[NamePair]:
    """
    One class's image word as oriented homotopy pairs (loops dropped, as in the word).

    Args:
        dyn: The symbolic dynamics.
        cd: The class record.

    Returns:
        One pair per token, in word order.
    """
    return [symbol_pair(dyn.naming, symbol) for symbol in cd.symbols]


def occurrence_in_names(
    naming: "ElementNaming", occurrence: tuple["ElementRef", "ElementRef"]
) -> NamePair:
    """
    An iterated element pair (a refined occurrence) in iterated names.

    Args:
        naming: The element naming.
        occurrence: Two iterated refs.

    Returns:
        The two names in the given order.
    """
    return (
        iterated_name(naming, occurrence[0]),
        iterated_name(naming, occurrence[1]),
    )


def refined_symbol_pair(dyn: "SymbolicDynamics", symbol: "Symbol") -> NamePair:
    """
    One REFINED word token as an oriented pair, inverse tokens reversed.

    A token of a split class is spelled by its refined child's iterated pair
    (the child's ``occurrence``, source end first); a token of an unsplit class
    by its homotopy pair, exactly as :func:`symbol_pair`.

    Args:
        dyn: The symbolic dynamics.
        symbol: A token of :attr:`SymbolicDynamics.refined_rules`.

    Returns:
        ``(source, target)`` for direction ``+1``, ``(target, source)`` for ``-1``.
    """
    if symbol.refined_index is None:
        return symbol_pair(dyn.naming, symbol)
    child = dyn.refined[symbol.bridge_class][symbol.refined_index - 1]
    pair = occurrence_in_names(dyn.naming, child.occurrence)
    return pair if symbol.direction > 0 else (pair[1], pair[0])


def refined_words_in_names(
    dyn: "SymbolicDynamics", cd: "ClassDynamics"
) -> dict[NamePair, list[NamePair]]:
    """
    A class's refined rule(s), keyed and spelled letter-free.

    A split class has one rule per refined child, keyed by the child's
    oriented iterated pair; an unsplit class has one rule keyed by its
    homotopy pair. Each word token is spelled by :func:`refined_symbol_pair`.

    Args:
        dyn: The symbolic dynamics.
        cd: The class record (resolved; an unresolved class has no rule).

    Returns:
        ``{rule key: [token pair, ...]}``; empty when the class has no rule.
    """
    children = dyn.refined.get(cd.bridge_class, [])
    if children:
        keyed = [
            (occurrence_in_names(dyn.naming, child.occurrence), child.name)
            for child in children
        ]
    else:
        keyed = [(homotopy_pair(dyn, cd), cd.letter)]
    return {
        key: [refined_symbol_pair(dyn, symbol) for symbol in dyn.refined_rules[name]]
        for key, name in keyed
        if name in dyn.refined_rules
    }


def itinerary_pairs(dyn: "SymbolicDynamics", cd: "ClassDynamics") -> Optional[list[NamePair]]:
    """
    A class's itinerary as its disjoint consecutive pairs of iterated names.

    Args:
        dyn: The symbolic dynamics.
        cd: The class record.

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
        occurrence_in_names(dyn.naming, (refs[i], refs[i + 1]))
        for i in range(0, len(refs), 2)
    ]


def _symbol_key(dyn: "SymbolicDynamics", symbol_name: str, refined: bool) -> NamePair:
    """The letter-free key of one transition-matrix symbol name."""
    if refined:
        for children in dyn.refined.values():
            for child in children:
                if child.name == symbol_name:
                    return occurrence_in_names(dyn.naming, child.occurrence)
    for cls, cd in dyn.classes.items():
        if cd.letter == symbol_name:
            return class_pair(dyn.naming, cls)
    for cls, letter in dyn.virtual_classes.items():
        if letter == symbol_name:
            return class_pair(dyn.naming, cls)
    raise KeyError(f"no class or refined child named {symbol_name!r}")


def matrix_by_classes(
    dyn: "SymbolicDynamics", *, refined: bool = True
) -> dict[NamePair, dict[NamePair, int]]:
    """
    The transition matrix keyed by element pairs instead of letters.

    A row or column of an unrefined class is keyed by its oriented homotopy
    pair; a refined child's by its oriented iterated pair.

    Args:
        dyn: The symbolic dynamics.
        refined: Over refined symbols (default) or bare classes.

    Returns:
        ``{row key: {column key: count}}`` with every row and column present.
    """
    names, matrix = dyn.transition_matrix(refined=refined)
    keys = [_symbol_key(dyn, name, refined) for name in names]
    return {
        row_key: {col_key: int(matrix[i, j]) for j, col_key in enumerate(keys)}
        for i, row_key in enumerate(keys)
    }

