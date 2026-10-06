"""Plumbing of the physical-law tier: shared-build products, the fingerprint, test factories.

* :func:`build_products` makes every topology product of a case once, so the
  laws only ever read cached products.
* :func:`fingerprint` snapshots everything a law could change on a shared
  build (generations, crossings, bridges, holes, pips, partition signatures,
  the session alphabet and the identity of every cached product).
* :func:`run_layer` runs every check of one layer on one case, insists each
  applicable check was not vacuous, and raises ONE ``AssertionError`` naming
  every failing check and its message.
* :func:`layer_test` turns one layer into a test function parametrized over
  :data:`cases.LAW_CASES` (one test per (layer, case)); :func:`known_issue_test`
  gives every ``cases.KNOWN_ISSUES`` pair of the layer a test of its own,
  ``xfail(strict=True)``, and :func:`run_layer` leaves that pair out, so a
  known issue stays precise and a fix still flips to XPASS.
"""

from __future__ import annotations

import traceback
from typing import TYPE_CHECKING, Any, Callable, Optional

import pytest

from cases import KNOWN_ISSUES, LAW_CASES, NOT_APPLICABLE, issue_marks
from helpers.laws import LAYERS, MUTATING_LAYERS, LawCheck, law_id

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


def layer_checks(layer: str) -> list[LawCheck]:
    """
    The checks of one layer, read-only or mutating.

    Args:
        layer: A key of ``helpers.laws.LAYERS`` or ``MUTATING_LAYERS``.

    Returns:
        The layer's checks, in run order.
    """
    if layer in MUTATING_LAYERS:
        return list(MUTATING_LAYERS[layer])
    return list(LAYERS[layer])


def _failure_text(exc: BaseException) -> str:
    """One line per failure: the exception, its message and the failing source line."""
    message = str(exc).strip() or "(no message)"
    frames = traceback.extract_tb(exc.__traceback__)
    where = ""
    if frames:
        frame = frames[-1]
        where = f" [{frame.filename.rsplit('/', 1)[-1]}:{frame.lineno}: {(frame.line or '').strip()}]"
    return f"{type(exc).__name__}: {message}{where}"


def run_layer(layer: str, case: "Case") -> dict[str, int]:
    """
    Run every check of a layer on a case and report every failure at once.

    A check whose ``(law, case)`` pair is in ``cases.NOT_APPLICABLE`` is
    skipped and not counted; one in ``cases.KNOWN_ISSUES`` is left out (it
    has its own ``xfail(strict=True)`` test, :func:`known_issue_test`). Every
    other check must pass AND check at least one item.

    Args:
        layer: The layer name.
        case: The built case.

    Returns:
        ``{law_id: items checked}`` of the checks that ran.

    Raises:
        AssertionError: One error naming every failing or vacuous check with
            its message.
    """
    counts: dict[str, int] = {}
    failures: list[str] = []
    for check in layer_checks(layer):
        law = law_id(check)
        if (law, case.name) in NOT_APPLICABLE or (law, case.name) in KNOWN_ISSUES:
            continue
        try:
            count = check(case)
        except Exception as exc:  # noqa: BLE001 - every failure is reported together
            failures.append(f"{law}: {_failure_text(exc)}")
            continue
        if count <= 0:
            failures.append(
                f"{law}: checked nothing; list ({law!r}, {case.name!r}) in "
                "cases.NOT_APPLICABLE if that is expected"
            )
            continue
        counts[law] = count
    if failures:
        raise AssertionError(
            f"{len(failures)} law(s) of layer {layer!r} failed on case {case.name!r}:\n  "
            + "\n  ".join(failures)
        )
    return counts


def _runnable(layer: str, case_name: str) -> bool:
    """True when at least one check of the layer is neither not applicable nor a known issue."""
    return any(
        (law_id(check), case_name) not in NOT_APPLICABLE
        and (law_id(check), case_name) not in KNOWN_ISSUES
        for check in layer_checks(layer)
    )


def _layer_params(layer: str) -> list:
    """The ``pytest.param`` list of law cases for one layer (skip when nothing runs)."""
    params = []
    for name in LAW_CASES:
        marks = []
        if not _runnable(layer, name):
            marks.append(pytest.mark.skip(reason=f"no applicable law of {layer!r} on {name!r}"))
        params.append(pytest.param(name, marks=marks, id=name))
    return params


def layer_test(layer: str) -> Callable[..., None]:
    """
    The test function of one layer, parametrized over the law cases.

    A read-only layer reads the module's shared ``law_case`` build; a layer
    of ``helpers.laws.MUTATING_LAYERS`` takes ``fresh_law_case``.

    Args:
        layer: The layer name.

    Returns:
        The test function, named ``test_<layer>_laws``.
    """
    fixture = "fresh_law_case" if layer in MUTATING_LAYERS else "law_case"
    if fixture == "law_case":

        def test(law_case: "Case") -> None:
            run_layer(layer, law_case)

    else:

        def test(fresh_law_case: "Case") -> None:
            run_layer(layer, fresh_law_case)

    test.__name__ = f"test_{layer}_laws"
    test.__qualname__ = test.__name__
    test.__doc__ = f"Every law of the {layer!r} layer: " + "; ".join(
        (check.__doc__ or law_id(check)).strip().splitlines()[0]
        for check in layer_checks(layer)
    )
    return pytest.mark.parametrize(fixture, _layer_params(layer), indirect=True)(test)


def known_issue_pairs(layer: str) -> list[tuple[str, str]]:
    """The ``(law_id, case)`` pairs of ``cases.KNOWN_ISSUES`` that belong to a layer."""
    laws_here = {law_id(check) for check in layer_checks(layer)}
    return sorted(pair for pair in KNOWN_ISSUES if pair[0] in laws_here)


def known_issue_test(layer: str) -> Optional[Callable[..., None]]:
    """
    One small ``xfail(strict=True)`` test per known issue of a layer.

    Args:
        layer: The layer name.

    Returns:
        The test function ``test_<layer>_known_issue`` parametrized over the
        layer's known ``(law, case)`` pairs, or None when the layer has none.
    """
    pairs = known_issue_pairs(layer)
    if not pairs:
        return None
    checks = {law_id(check): check for check in layer_checks(layer)}
    fixture = "fresh_law_case" if layer in MUTATING_LAYERS else "law_case"
    params = [
        pytest.param(law, case, marks=issue_marks(law, case), id=f"{law}-{case}")
        for law, case in pairs
    ]

    if fixture == "law_case":

        def test(law: str, law_case: "Case") -> None:
            assert_non_vacuous(checks[law], law_case)

    else:

        def test(law: str, fresh_law_case: "Case") -> None:
            assert_non_vacuous(checks[law], fresh_law_case)

    test.__name__ = f"test_{layer}_known_issue"
    test.__qualname__ = test.__name__
    test.__doc__ = f"A known issue of the {layer!r} layer, kept apart so a fix flips to XPASS."
    return pytest.mark.parametrize(("law", fixture), params, indirect=[fixture])(test)
