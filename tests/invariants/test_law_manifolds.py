"""Law tier: the manifold invariants, on every manifold (and bridge) of every case.

cdist is non-decreasing along the curve (ties are legitimate at a
high-stretch fold, so never strict here), no node juts off the curve (no
geometric spike), ``c_iterate = stretch_param * c`` along the iterate chain,
and the curve's rows and their iterate links are acyclic and consistent.
"""

from __future__ import annotations

from helpers.law_tier import layer_tests

test_manifolds_laws, test_manifolds_set_aside = layer_tests("manifolds")
