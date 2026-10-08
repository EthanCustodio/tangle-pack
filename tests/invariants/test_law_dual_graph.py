"""Law tier: the dual graph.

Each stable edge has two open wall nodes or one unified node; one face node
per merged face with a single unbounded outer node; a face attaches to the
side node of the edge it geometrically lies on; the unified edges are
exactly the strong pip's fundamental segment ``(f^k(q0), q0]`` on the pip's
OWN branch (author, 2026-09-30), and that is where ``f(S)`` stops: on every
stable branch of a pip's fixed point the unified span is ``S \\ f(S)`` (the
walk's precondition, PK-5 / DEC-20); the graph is bipartite with walls of
degree 1 and unified nodes of degree 2.
"""

from __future__ import annotations

from helpers.law_tier import layer_tests

test_dual_graph_laws, test_dual_graph_set_aside = layer_tests("dual_graph")
