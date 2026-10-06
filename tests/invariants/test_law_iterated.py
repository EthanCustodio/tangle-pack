"""Law tier: the iterated homotopy partition, element names and the minimal trellis.

An iterated element lies inside its homotopy parent; every homotopy boundary
survives the cut with its closedness and every new boundary is an image
bridge endpoint on that row; every crossing (and every stable arc midpoint of
the minimal arrangement) has one owner per side; an element under a
same-branch image bridge records the bridge covering it. Names agree with the
structure (side letter, parent subscript, round trips). The minimal trellis
keeps the hole bridges and the images of the active ones and nothing else,
keeps the right nodes, and its sparse arrangement carries exactly the kept
bridges and closes up.
"""

from __future__ import annotations

from helpers import laws
from helpers.law_tier import law_test

test_iterated_child_inside_parent = law_test(laws.check_iterated_child_inside_parent)
test_iterated_keeps_homotopy_boundaries = law_test(laws.check_iterated_keeps_homotopy_boundaries)
test_iterated_unique_owner = law_test(laws.check_iterated_unique_owner)
test_iterated_cut_provenance = law_test(laws.check_iterated_cut_provenance)
test_names_agree_with_structure = law_test(laws.check_names_agree_with_structure)
test_minimal_trellis_bridges = law_test(laws.check_minimal_trellis_bridges)
test_minimal_trellis_nodes_and_arrangement = law_test(
    laws.check_minimal_trellis_nodes_and_arrangement
)
