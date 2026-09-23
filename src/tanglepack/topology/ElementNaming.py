"""
Human-readable names for the elements of the stable partitions.

The homotopy partition names its elements side-tagged, anchor outward and
1-based: ``L_1, L_2, ...`` on the left and ``R_1, R_2, ...`` on the right of
each stable branch. The iterated homotopy partition refines those elements,
and each child carries a superscript counting the children of its parent
anchor outward: ``R_1^1, R_1^2``. An element with a single child (an unsplit
element) prints plain, ``R_3``. When more than one stable branch is
partitioned every name is prefixed with its branch tag, ``p3@1.0:R_2^1``, in
the same ``p{period}@{orbit}.{branch}`` form as
:attr:`~tanglepack.topology.TopologyResults.ElementRef.label`.

:class:`ElementNaming` is the ONE place that says which family an
:class:`~tanglepack.topology.TopologyResults.ElementRef` belongs to. A ref
carries no family tag, so the refs of the homotopy family and of the iterated
family collide by value; the naming keeps the two in separate tables and its
methods say which family each argument is read from.

Dev Notes:

* Subscripts are homotopy element ids plus one, so a child's subscript is its
  parent's. Refined bridge-class symbols (``a_1, a_2``) subdivide a class the
  same way and live in ``SymbolicDynamics``, not here.
* Singletons (``[x, x]``) are numbered like any other element; nothing in a
  name says whether an element is a singleton.
* An iterated element whose ``parent_element_id`` is None (only a hand-built
  family; ``IteratedHomotopyPartition.from_minimal`` always stamps it) is read
  as its own parent: the homotopy element with the same id on the same branch
  and side, which must exist.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Iterator, Optional, TYPE_CHECKING

from .TopologyResults import ElementRef, Side

if TYPE_CHECKING:  # pragma: no cover - typing only
    from ..numerics.Intersection import ManifoldKey
    from .PartitionFamily import PartitionFamily

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

#: The side letter of a name, per partition side.
SIDE_LETTER: dict[Side, str] = {"left": "L", "right": "R"}


@dataclass(frozen=True)
class ElementName:
    """
    The name of one partition element: ``R_1^2`` and friends.

    Attributes:
        branch_key: Manifold key of the stable branch the element lies on.
        side: Which side's partition the element belongs to.
        subscript: The 1-based anchor-outward index of the HOMOTOPY element
            (its ``element_id + 1``); a child shares its parent's subscript.
        superscript: The 1-based index of the element among the children of
            its homotopy parent, anchor outward; 1 for an unsplit element.
        siblings: How many children the parent has; 1 means unsplit and the
            name prints without a superscript.
        tagged: True when the branch tag is shown (more than one stable
            branch is partitioned in the family).
    """

    branch_key: "ManifoldKey"
    side: Side
    subscript: int
    superscript: int
    siblings: int
    tagged: bool

    @property
    def side_letter(self) -> str:
        """``"L"`` or ``"R"``."""
        return SIDE_LETTER[self.side]

    @property
    def is_split(self) -> bool:
        """True when the parent has more than one child (the superscript prints)."""
        return self.siblings > 1

    @property
    def tag(self) -> str:
        """
        The branch tag, ``p{period}@{orbit}.{branch}``, shown or not.

        Returns:
            The same branch text
            :attr:`~tanglepack.topology.TopologyResults.ElementRef.label` uses.
        """
        fixed_point, _stability, orbit, branch = self.branch_key
        return f"p{fixed_point.period}@{orbit}.{branch}"

    @property
    def text(self) -> str:
        """
        The plain-text name: ``R_1^2``, ``R_3`` when unsplit, ``p3@1.0:R_2^1``
        when tagged.
        """
        body = f"{self.side_letter}_{self.subscript}"
        if self.is_split:
            body += f"^{self.superscript}"
        return f"{self.tag}:{body}" if self.tagged else body

    @property
    def mathtext(self) -> str:
        """
        The name as matplotlib mathtext: ``$R_{1}^{2}$`` / ``$R_{3}$``; a
        branch tag is plain prefix text outside the math (``p3@1.0:$R_{2}^{1}$``).
        """
        body = f"{self.side_letter}_{{{self.subscript}}}"
        if self.is_split:
            body += f"^{{{self.superscript}}}"
        math = f"${body}$"
        return f"{self.tag}:{math}" if self.tagged else math

    def __str__(self) -> str:
        return self.text

    def __repr__(self) -> str:
        return f"ElementName({self.text!r})"


class ElementNaming:
    """
    Names for the elements of a homotopy partition and, optionally, of the
    iterated partition refining it.

    Two separate tables are kept, one per family, because a ref carries no
    family tag. Every method says which family its argument is read from.
    With no iterated family every homotopy element is its own single child,
    so :meth:`name`, :meth:`parent_of`, :meth:`children_of` and :meth:`ref_of`
    act on the homotopy refs themselves.

    Attributes:
        homotopy: The homotopy family.
        iterated: The iterated family refining it, or None.
        tagged: Whether names carry their branch tag (more than one stable
            branch is partitioned in the homotopy family).
    """

    def __init__(
        self,
        homotopy: "PartitionFamily",
        iterated: Optional["PartitionFamily"] = None,
    ) -> None:
        """
        Args:
            homotopy: The homotopy partition family (one result per
                ``(branch, side)``; element ids anchor outward).
            iterated: The iterated homotopy partition family refining it. Its
                elements point at their parents through
                ``PartitionInterval.parent_element_id``.

        Raises:
            ValueError: If an iterated result covers a ``(branch, side)`` the
                homotopy family does not, or an iterated element names a
                parent the homotopy result does not have.
        """
        self.homotopy = homotopy
        self.iterated = iterated
        self.tagged: bool = len(homotopy.branch_keys) > 1

        # Homotopy refs -> names (superscript 1 / siblings 1 until children
        # are counted), and the ordered list of homotopy refs.
        self._homotopy_refs: list[ElementRef] = []
        self._homotopy_names: dict[ElementRef, ElementName] = {}
        # Iterated refs -> names / parents; homotopy refs -> children (anchor
        # outward). Never mixed with the homotopy tables.
        self._iterated_refs: list[ElementRef] = []
        self._iterated_names: dict[ElementRef, ElementName] = {}
        self._parents: dict[ElementRef, ElementRef] = {}
        self._children: dict[ElementRef, list[ElementRef]] = {}
        # Name -> the ref it names (iterated when there is an iterated family,
        # else homotopy).
        self._ref_of_name: dict[ElementName, ElementRef] = {}
        self._ref_of_text: dict[str, ElementRef] = {}

        for result in homotopy.as_list():
            for interval in result.intervals:
                ref = result.ref(interval.element_id)
                self._homotopy_refs.append(ref)
                self._children[ref] = []

        if iterated is None:
            for ref in self._homotopy_refs:
                self._children[ref] = [ref]
                self._parents[ref] = ref
            self._iterated_refs = list(self._homotopy_refs)
        else:
            self._link_iterated(iterated)

        for parent in self._homotopy_refs:
            children = self._children[parent]
            siblings = len(children)
            self._homotopy_names[parent] = ElementName(
                parent.branch_key, parent.side, parent.element_id + 1, 1, 1, self.tagged
            )
            for position, child in enumerate(children, start=1):
                name = ElementName(
                    parent.branch_key,
                    parent.side,
                    parent.element_id + 1,
                    position,
                    siblings,
                    self.tagged,
                )
                self._iterated_names[child] = name
                self._ref_of_name[name] = child
                self._ref_of_text[name.text] = child

        logger.debug(
            "element naming: %d homotopy element(s), %d named element(s), tagged=%s",
            len(self._homotopy_refs),
            len(self._iterated_refs),
            self.tagged,
        )

    def _link_iterated(self, iterated: "PartitionFamily") -> None:
        """Fill the iterated tables and the children lists, validating parents."""
        homotopy = self.homotopy
        for result in iterated.as_list():
            key = (result.branch_key, result.side)
            if key not in homotopy.results:
                raise ValueError(
                    f"the iterated partition covers branch {result.branch_key[1:]} "
                    f"side {result.side!r}, which the homotopy partition does not; "
                    "every iterated element needs a homotopy parent"
                )
            parent_result = homotopy.results[key]
            for interval in result.intervals:
                ref = result.ref(interval.element_id)
                parent_id = interval.parent_element_id
                if parent_id is None:
                    parent_id = interval.element_id
                    logger.debug(
                        "iterated element %s carries no parent id; it is read as "
                        "its own parent",
                        ref.label,
                    )
                if not 0 <= parent_id < len(parent_result.intervals):
                    raise ValueError(
                        f"iterated element {ref.label} names parent element "
                        f"#{parent_id}, but the {result.side!r} homotopy partition "
                        f"of branch {result.branch_key[1:]} has only "
                        f"{len(parent_result.intervals)} element(s)"
                    )
                parent = parent_result.ref(parent_id)
                self._iterated_refs.append(ref)
                self._parents[ref] = parent
                self._children[parent].append(ref)
        for parent, children in self._children.items():
            children.sort(key=lambda child: child.element_id)

    # ── the two families ────────────────────────────────────────────────────

    @property
    def has_iterated(self) -> bool:
        """True when an iterated family is named (else the homotopy refs are
        the named refs)."""
        return self.iterated is not None

    @property
    def homotopy_refs(self) -> list[ElementRef]:
        """The homotopy refs, result by result, anchor outward."""
        return list(self._homotopy_refs)

    @property
    def refs(self) -> list[ElementRef]:
        """The named refs (iterated, or homotopy without an iterated family),
        result by result, anchor outward."""
        return list(self._iterated_refs)

    @property
    def names(self) -> list[ElementName]:
        """The names of :attr:`refs`, in the same order."""
        return [self._iterated_names[ref] for ref in self._iterated_refs]

    @property
    def homotopy_names(self) -> list[ElementName]:
        """The names of :attr:`homotopy_refs`, in the same order."""
        return [self._homotopy_names[ref] for ref in self._homotopy_refs]

    def items(self) -> Iterator[tuple[ElementRef, ElementName]]:
        """Iterate ``(named ref, name)`` pairs in :attr:`refs` order."""
        for ref in self._iterated_refs:
            yield ref, self._iterated_names[ref]

    # ── lookups ─────────────────────────────────────────────────────────────

    def homotopy_name(self, ref: ElementRef) -> ElementName:
        """
        The name of one HOMOTOPY element (superscript 1, siblings 1).

        Args:
            ref: A ref of the homotopy family.

        Returns:
            Its :class:`ElementName`.

        Raises:
            KeyError: If the ref names no homotopy element.
        """
        try:
            return self._homotopy_names[ref]
        except KeyError:
            raise KeyError(
                f"{ref.label} is not an element of the homotopy partition"
            ) from None

    def name(self, iterated_ref: ElementRef) -> ElementName:
        """
        The name of one ITERATED element (a homotopy element without an
        iterated family).

        Args:
            iterated_ref: A ref of the iterated family.

        Returns:
            Its :class:`ElementName`, with the parent's subscript and the
            child's superscript.

        Raises:
            KeyError: If the ref names no element of that family.
        """
        try:
            return self._iterated_names[iterated_ref]
        except KeyError:
            raise KeyError(
                f"{iterated_ref.label} is not an element of the "
                f"{self._named_kind} partition"
            ) from None

    def parent_of(self, iterated_ref: ElementRef) -> ElementRef:
        """
        The homotopy element one iterated element refines.

        Args:
            iterated_ref: A ref of the iterated family.

        Returns:
            The HOMOTOPY ref of its parent (the ref itself without an
            iterated family).

        Raises:
            KeyError: If the ref names no element of the iterated family.
        """
        try:
            return self._parents[iterated_ref]
        except KeyError:
            raise KeyError(
                f"{iterated_ref.label} is not an element of the "
                f"{self._named_kind} partition"
            ) from None

    def children_of(self, homotopy_ref: ElementRef) -> list[ElementRef]:
        """
        The iterated elements refining one homotopy element.

        Args:
            homotopy_ref: A ref of the homotopy family.

        Returns:
            Its children's ITERATED refs, anchor outward (the ref itself, alone,
            without an iterated family).

        Raises:
            KeyError: If the ref names no homotopy element.
        """
        try:
            return list(self._children[homotopy_ref])
        except KeyError:
            raise KeyError(
                f"{homotopy_ref.label} is not an element of the homotopy partition"
            ) from None

    def ref_of(self, name: ElementName) -> ElementRef:
        """
        The element a name names.

        Args:
            name: An :class:`ElementName` this naming produced (a name built
                by hand with the same fields also resolves).

        Returns:
            The ITERATED ref (the homotopy ref without an iterated family).

        Raises:
            KeyError: If no element carries that name.
        """
        try:
            return self._ref_of_name[name]
        except KeyError:
            raise KeyError(f"no element is named {name.text}") from None

    def lookup(self, text: str) -> ElementRef:
        """
        The element whose plain-text name is ``text``.

        Args:
            text: A name as :attr:`ElementName.text` prints it (``"R_1^2"``,
                ``"L_3"``, ``"p3@1.0:R_2^1"``).

        Returns:
            The ITERATED ref (the homotopy ref without an iterated family).

        Raises:
            KeyError: If no element prints that name.
        """
        try:
            return self._ref_of_text[text]
        except KeyError:
            raise KeyError(f"no element is named {text!r}") from None

    # ── reporting ───────────────────────────────────────────────────────────

    @property
    def _named_kind(self) -> str:
        """The kind word of the family :meth:`name` reads (for messages)."""
        return "iterated" if self.has_iterated else "homotopy"

    def describe(self) -> str:
        """
        A human-readable report: one line per homotopy element with its
        children's names, grouped by ``(branch, side)``.

        Returns:
            The report.
        """
        lines = [
            f"element naming: {len(self._homotopy_refs)} homotopy element(s), "
            f"{len(self._iterated_refs)} {self._named_kind} element(s)"
            + (", branch-tagged" if self.tagged else "")
        ]
        current: Optional[tuple["ManifoldKey", Side]] = None
        for parent in self._homotopy_refs:
            key = (parent.branch_key, parent.side)
            if key != current:
                current = key
                name = self._homotopy_names[parent]
                lines.append(f"{parent.side} partition of {name.tag}:")
            children = self._children[parent]
            child_names = ", ".join(self._iterated_names[c].text for c in children)
            parent_name = self._homotopy_names[parent].text
            if len(children) == 1 and children[0] == parent and not self.has_iterated:
                lines.append(f"  {parent_name}  (#{parent.element_id})")
            else:
                lines.append(
                    f"  {parent_name}  (#{parent.element_id}) -> {child_names or '(no children)'}"
                )
        return "\n".join(lines)

    def __len__(self) -> int:
        """The number of named (iterated, else homotopy) elements."""
        return len(self._iterated_refs)

    def __iter__(self) -> Iterator[ElementRef]:
        return iter(self._iterated_refs)

    def __contains__(self, ref: object) -> bool:
        """Whether a ref names an element of the NAMED family (see :meth:`name`)."""
        return ref in self._iterated_names

    def __repr__(self) -> str:
        return (
            f"<ElementNaming homotopy={len(self._homotopy_refs)} "
            f"{self._named_kind}={len(self._iterated_refs)} tagged={self.tagged}>"
        )
