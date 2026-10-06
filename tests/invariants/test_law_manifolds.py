"""Law tier: the manifold invariants, on every manifold (and bridge) of every case.

cdist is non-decreasing along the curve (ties are legitimate at a
high-stretch fold, so never strict here), no node juts off the curve (no
geometric spike), ``c_iterate = stretch_param * c`` along the iterate chain,
and the geometric and iterate linked lists are acyclic and consistent.
"""

from __future__ import annotations

from helpers import laws
from helpers.law_tier import law_test

test_manifold_cdist_monotone = law_test(laws.check_manifold_cdist_monotone)
test_manifold_no_spikes = law_test(laws.check_manifold_no_spikes)
test_manifold_iterate_law = law_test(laws.check_manifold_iterate_law)
test_manifold_one_to_one = law_test(laws.check_manifold_one_to_one)
