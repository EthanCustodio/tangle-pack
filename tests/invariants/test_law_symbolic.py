"""Law tier: symbolic dynamics.

Every itinerary has even length and splits into same-side pairs; each pair,
read through its homotopy parents, is a loop or one class or its inverse --
exactly the word's tokens; each class element lands in one iterated element;
a split class's children keep its ``BridgeClass`` and inherit its word;
member matching is consistent; inert classes are outside the transition
graph and matrix; and the matrix counts the word tokens.
"""

from __future__ import annotations

from helpers import laws
from helpers.law_tier import law_test

test_itineraries_even = law_test(laws.check_itineraries_even)
test_itinerary_pairs_same_side = law_test(laws.check_itinerary_pairs_same_side)
test_itinerary_pairs_are_classes = law_test(laws.check_itinerary_pairs_are_classes)
test_landings_contained = law_test(laws.check_landings_contained)
test_refined_children_inherit_word = law_test(laws.check_refined_children_inherit_word)
test_member_matching_consistent = law_test(laws.check_member_matching_consistent)
test_inert_classes_outside_transitions = law_test(laws.check_inert_classes_outside_transitions)
test_matrix_is_token_counts = law_test(laws.check_matrix_is_token_counts)
