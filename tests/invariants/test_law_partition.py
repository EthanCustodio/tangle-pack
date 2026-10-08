"""Law tier: holes and the stable (homotopy) partition.

The partition covers each branch on each side in contiguous positional
elements, every crossing has exactly one owner per side, pinched singletons
own exactly their point. Every hole is classified. THE HOLE RULE (author,
2026-10-05, I1): a hole maps backward onto the same side of the BRIDGE it is
punched in -- every hole of an origin shares the reference hole's
``bridge_side``, and the reference hole carried back by the real inverse map
lies on that side of the propagated hole's own bridge. It is NOT "the same
side of the stable manifold": a backward hole's stable-manifold row may
differ from its origin's, and no law compares rows along an orbit. A direct
hole's side is the side of its bridge its coordinates sit on. The openings
law: a direct hole opens its inward pair on its own bridge's row; every
opening sits on the row of ITS OWN bridge at that bound, on the branch the
record names; at an anchor (one point on every branch, D4) the signs pick the
branch, and the opened branch's first node lies on the hole's side of its
bridge (D6, checked by geometry); a bound opens nothing only at an anchor or
the branch end. Propagated holes land on the
branch the cycle predicts and stop at periodicity; nothing at iterate
``>= k_value`` is punched (the firm half of holes-backward-only; the
``+1 .. +(k-1)`` exemption is provisional and untested); reference pairs are
consecutive on both branches inside the reference window.
"""

from __future__ import annotations

from helpers.law_tier import layer_tests

test_partition_laws, test_partition_set_aside = layer_tests("partition")
