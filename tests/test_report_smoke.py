"""Smoke test of the human-readable reports (``summary()``, ``describe()``, repr).

Author decision 5 (2026-10-05): report text is never pinned word for word. A
report is only required to be non-empty and to mention each class or element it
reports on; the structure behind it is asserted elsewhere.
"""

from __future__ import annotations

from tanglepack.topology.DualWalk import Walk, WalkSearch


def test_reports_are_non_empty_and_mention_what_they_report(k10_partitioned) -> None:
    """Every report on the k=10 tangle is non-empty text naming its subjects."""
    session, _fp = k10_partitioned
    dual = session.dual_graph()
    dynamics = session.symbolic_dynamics()

    texts = [
        dual.summary(),
        repr(dual),
        repr(dual.unified_nodes[0]),
        repr(dual.face_nodes[0]),
    ]

    naming = dual.naming
    assert naming is not None
    naming_report = naming.describe()
    texts.append(naming_report)
    for parent in naming.homotopy_refs:
        assert naming.homotopy_name(parent).text in naming_report

    classes_report = session.describe_bridge_classes()
    texts.append(classes_report)
    for entry in session.bridge_classes():
        if entry.letter is not None:
            assert entry.letter in classes_report

    searches = [cd.search for cd in dynamics.classes.values() if cd.search is not None]
    assert searches, "the k=10 active class is resolved by a dual-graph walk"
    for search in searches:
        assert isinstance(search, WalkSearch)
        texts.append(repr(search))
        if search.walk is not None:
            assert isinstance(search.walk, Walk)
            texts.append(repr(search.walk))

    texts.append(dynamics.describe())
    texts.append(repr(dynamics))

    assert all(isinstance(text, str) and text.strip() for text in texts)
