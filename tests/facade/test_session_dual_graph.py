"""
TangleSession.iterated_partition: describe.

The plot delegates are ``tests/plotting/``; the cache contract is ``tests/facade/test_session_caches.py``; session = direct
build is ``tests/facade/test_session_equivalence.py``.
"""

from __future__ import annotations


def test_describe_iterated_partition(k10_partitioned):
    session, fp = k10_partitioned
    assert session.describe_iterated_partition() == session.iterated_partition().describe()
