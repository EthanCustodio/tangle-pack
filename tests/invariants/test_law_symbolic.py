"""Law tier: symbolic dynamics.

Every itinerary has even length; a pair on one stable branch joins two
elements of one side, and only such a pair may be flagged cross-side (D4: at
an inversion point L of 0.0 is R of 0.1); every active class that does not
land on a singleton has a walk, and no itinerary names a virtual class (D3:
the shortest walk IS the word); each pair, read through its homotopy parents,
is a loop or one class or its inverse -- exactly the word's tokens; each class
element lands in one iterated element; a split class's children keep its
``BridgeClass`` and inherit its word; member matching is consistent; inert
classes are outside the transition graph and matrix; and the matrix counts
the word tokens.
"""

from __future__ import annotations

from helpers.law_tier import layer_tests

test_symbolic_laws, test_symbolic_set_aside = layer_tests("symbolic")
