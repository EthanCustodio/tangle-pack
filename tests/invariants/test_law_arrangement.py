"""Law tier: the planar arrangement and its regions.

Euler's formula holds per connected component; every region is a closed
cycle whose representative point survives an independent ray cast; regions
are pairwise interior-disjoint (a containing face swallows several); a
region's combinatorial image is where the real map sends its representative
point, with the same area (the maps are area preserving); and image and
preimage undo each other where both are recorded.
"""

from __future__ import annotations

from helpers.law_tier import known_issue_test, layer_test

test_arrangement_laws = layer_test("arrangement")
test_arrangement_known_issue = known_issue_test("arrangement")
