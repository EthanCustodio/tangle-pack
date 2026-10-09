"""Fixtures of the physical-law tier (``tests/invariants/``).

The law tier is the ONE place the suite shares a build between tests (author
decision 7, 2026-10-05): each case is built once per RUN (``law_case``,
session scope since WP10, C12; it was once per module, which rebuilt every
case for each of the ~12 law modules) indirect over ``cases.BUILDERS`` names,
with every topology product made inside the fixture, and every test only
reads it. An autouse FINGERPRINT GUARD snapshots the shared build before each test
and compares after it, so a test that mutates the build fails loudly instead
of leaking state into the next law.

The two laws that must mutate a session (a recompute of the crossings) take
``fresh_law_case`` instead: a function-scoped build of their own.
"""

from __future__ import annotations

from typing import Iterator

import pytest

from cases import BUILDERS, Case
from helpers.law_tier import build_products, fingerprint


@pytest.fixture(scope="session")
def law_case(request: pytest.FixtureRequest) -> Case:
    """One read-only build per run of the case named by ``request.param``.

    Every product the laws read (arrangement, bridge classes, minimal
    trellis, iterated partition, dual graph, symbolic dynamics) is built
    here, so a law never triggers a build.
    """
    case = BUILDERS[request.param]()
    build_products(case)
    return case


@pytest.fixture
def fresh_law_case(request: pytest.FixtureRequest) -> Case:
    """A fresh, function-scoped build for a law that mutates its session."""
    return BUILDERS[request.param]()


@pytest.fixture(autouse=True)
def _fingerprint_guard(request: pytest.FixtureRequest) -> Iterator[None]:
    """Fail a test that changed the shared ``law_case`` build."""
    if "law_case" not in request.fixturenames:
        yield
        return
    case = request.getfixturevalue("law_case")
    before = fingerprint(case)
    yield
    after = fingerprint(case)
    assert after == before, (
        f"the test mutated the shared {case.name!r} build: "
        + ", ".join(key for key in before if before[key] != after.get(key))
    )
