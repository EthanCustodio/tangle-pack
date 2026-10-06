"""
TangleSession.symbolic_dynamics / describe_symbolic_dynamics.

Mirrors ``test_session_dual_graph.py``. Nothing here pins a registry id, a
bridge count or a fixture word: those are pinned once in ``tests/golden/``.
The plot delegates are ``tests/plotting/``; the cache contract and session =
direct build are ``tests/facade/``.
"""

from __future__ import annotations


# --------------------------------------------------------------------------- #
# k=10: the report
# --------------------------------------------------------------------------- #
def test_describe_symbolic_dynamics(k10_partitioned):
    session, fp = k10_partitioned
    text = session.describe_symbolic_dynamics([fp])
    assert text and isinstance(text, str)
    assert text == session.symbolic_dynamics([fp]).describe()
    for cd in session.symbolic_dynamics([fp]).classes.values():
        assert cd.letter in text
