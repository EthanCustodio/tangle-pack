"""Every cached session product equals the same product built directly.

:class:`~tanglepack.TangleSession` is a facade: it gathers partitions, holes
and strong pips from its per-fixed-point trellises and hands them to the
topology functions. ``helpers/direct.py`` rebuilds each product WITHOUT the
session's caches -- from the public topology constructors over the
partitions, holes and pips read off the per-fixed-point trellises -- and this
module checks the session's result against it, on the ``k10`` case and on the
two-fixed-point ``nested`` case (author decision 6, 2026-10-05).

Comparisons are structural and letter-free: ids within ONE build, element
names, and words as oriented element pairs. That the session stays equal to a
direct build after every cache-relevant event is ``test_session_caches.py``.
"""

from __future__ import annotations

from typing import Callable

import pytest

from cases import Case, build_k10, build_nested
from helpers.direct import PRODUCT_NAMES, direct_view, view_of

_BUILDERS: dict[str, Callable[[], Case]] = {
    "k10": build_k10,
    "nested": build_nested,
}


@pytest.mark.parametrize("case_name", list(_BUILDERS))
@pytest.mark.parametrize("product", PRODUCT_NAMES)
def test_session_product_equals_a_direct_build(product: str, case_name: str) -> None:
    """The session's cached product is the direct build's, structurally."""
    case = _BUILDERS[case_name]()

    from_session = view_of(product, getattr(case.session, product)())
    from_direct = direct_view(case, product)

    assert from_session == from_direct
    assert from_session, f"{product} of {case_name} is empty"
