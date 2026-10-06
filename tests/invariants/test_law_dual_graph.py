"""Law tier: the dual graph.

Each stable edge has two open wall nodes or one unified node; one face node
per merged face with a single unbounded outer node; a face attaches to the
side node of the edge it geometrically lies on; the unified edges are
exactly the strong pip's fundamental segment ``(f^k(q0), q0]`` on the pip's
OWN branch (author, 2026-09-30); the graph is bipartite with walls of degree
1 and unified nodes of degree 2.
"""

from __future__ import annotations

from helpers.law_tier import layer_test

test_dual_graph_laws = layer_test("dual_graph")
