"""Law tier: bridges and their identity.

A bridge is partial exactly when it has no ``BridgeId``; its id is its two
crossings in increasing unstable cdist, both on the bridge's OWN unstable
branch (the period-3 regression of the 2026-09 codebase audit) and inside its
span; ``bridges_at`` indexes exactly the two ends; each id has one registered
copy; two bridges of one branch only touch; and the derived image genealogy
lands on the advanced branch and round-trips through the preimage.
"""

from __future__ import annotations

from helpers.law_tier import layer_test

test_bridges_laws = layer_test("bridges")
