"""Plumbing of the physical-law tier: shared-build products, the fingerprint, test factories.

* :func:`build_products` makes every topology product of a case once, so the
  laws only ever read cached products.
* :func:`fingerprint` snapshots everything a law could change on a shared
  build (generations, crossings, bridges, holes, pips, partition signatures,
  the session alphabet and the identity of every cached product).
* :func:`law_test` turns one ``helpers.laws`` check into a test function
  parametrized over the law cases (``cases.law_params``: known issues are
  ``xfail(strict=True)``, not-applicable pairs are skipped) that also insists
  the law checked at least one item.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Callable

import pytest

from cases import law_params
from helpers.laws import LawCheck, law_id

if TYPE_CHECKING:  # pragma: no cover - typing only
    from cases import Case

#: The session products the law tier reads, in build order.
PRODUCTS: tuple[str, ...] = (
    "arrangement",
    "bridge_classes",
    "minimal_trellis",
    "iterated_partition",
    "dual_graph",
    "symbolic_dynamics",
)


def build_products(case: "Case") -> None:
    """
    Build (and so cache) every product the laws read.

    Args:
        case: The case whose session to fill.
    """
    for fp in case.fixed_points:
        case.session.trellis(fp)
    case.session.trellis()
    for product in PRODUCTS:
        getattr(case.session, product)()


def fingerprint(case: "Case") -> dict[str, Any]:
    """
    A snapshot of the parts of a build a law could change.

    Args:
        case: The built case.

    Returns:
        ``{aspect: value}``; two snapshots compare equal exactly when the
        build was left alone.
    """
    session = case.session
    registry = case.registry
    trellises = [session.trellis(fp) for fp in case.fixed_points]
    return {
        "workbench generation": case.workbench.generation,
        "registry generation": registry.generation,
        "crossings": tuple(
            (iid, float(ix.unstable_cdist), float(ix.stable_cdist)) for iid, ix in registry
        ),
        "bridge ids": tuple(sorted(str(bridge.id) for bridge in case.workbench.bridges)),
        "holes": tuple(
            tuple(
                (h.bounding_ids, h.iterate, h.bridge_side, tuple(h.openings or ()))
                for h in trellis.holes
            )
            for trellis in trellises
        ),
        "pips": tuple(trellis.strong_pip for trellis in trellises),
        "partition signature": session.homotopy_partition().signature(),
        "alphabet": tuple(session.bridge_alphabet.assigned.values()),
        "products": tuple(id(getattr(session, product)()) for product in PRODUCTS),
    }


def assert_non_vacuous(check: LawCheck, case: "Case") -> int:
    """
    Run one law on a case and insist it checked something.

    Args:
        check: The law check.
        case: The built case.

    Returns:
        The number of items checked.

    Raises:
        AssertionError: The law failed, or checked nothing (a ``(law, case)``
            pair with nothing to check belongs in ``cases.NOT_APPLICABLE``).
    """
    count = check(case)
    assert count > 0, (
        f"law {law_id(check)!r} checked nothing on case {case.name!r}; list the pair in "
        "cases.NOT_APPLICABLE if that is expected"
    )
    return count


def law_test(check: LawCheck, *, fixture: str = "law_case") -> Callable[..., None]:
    """
    A test function for one law, parametrized over the law cases.

    Args:
        check: The law check (``helpers.laws.check_<law>``).
        fixture: ``"law_case"`` (the module's shared read-only build) or
            ``"fresh_law_case"`` (a fresh build, for a law that mutates).

    Returns:
        The test function, named ``test_<law>``, documented by the check.
    """
    name = law_id(check)

    if fixture == "law_case":

        def test(law_case: "Case") -> None:
            assert_non_vacuous(check, law_case)

    else:

        def test(fresh_law_case: "Case") -> None:
            assert_non_vacuous(check, fresh_law_case)

    test.__name__ = f"test_{name}"
    test.__qualname__ = test.__name__
    test.__doc__ = check.__doc__
    return pytest.mark.parametrize(fixture, law_params(name), indirect=True)(test)
