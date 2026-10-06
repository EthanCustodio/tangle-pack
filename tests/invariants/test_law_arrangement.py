"""Law tier: the planar arrangement and its regions.

Euler's formula holds per connected component; every region is a closed
cycle whose representative point survives an independent ray cast; regions
are pairwise interior-disjoint (a containing face swallows several); a
region's combinatorial image is where the real map sends its representative
point, with the same area (the maps are area preserving); and image and
preimage undo each other where both are recorded.
"""

from __future__ import annotations

from helpers import laws
from helpers.law_tier import law_test

test_arrangement_euler = law_test(laws.check_arrangement_euler)
test_arrangement_regions_sound = law_test(laws.check_arrangement_regions_sound)
test_arrangement_regions_disjoint = law_test(laws.check_arrangement_regions_disjoint)
test_arrangement_image_of = law_test(laws.check_arrangement_image_of)
test_arrangement_preimage_inverts_image = law_test(
    laws.check_arrangement_preimage_inverts_image
)
