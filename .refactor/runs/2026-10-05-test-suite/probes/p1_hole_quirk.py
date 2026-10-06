"""P1: the hole quirk -- do rows really flip between a direct hole and its backward images?

For every hole in build_period3() and the k=2.8 two-blast build, print
bounding_ids, near/far bounds, whether the near bound is an anchor, openings,
the row of the hole's bridge at each bound (crossing-sign rule and geometric
_row_at), the hole's own GEOMETRIC side of the stable manifold, and for each
bound with no opening whether its halves face an anchor (no node below) or a
removed tail (no node above).

Decision rule (plan Phase 0, P1):
 (a) every difference from the origin's (which, row) set is an opening dropped
     at an anchor or tail -> openings law;
 (b) rows really flip between the origin and its backward images -> xfail.
"""

from __future__ import annotations

import logging
import sys
from collections import defaultdict

import numpy as np

from tanglepack.examples.henon_cases import build_period3
from tanglepack.topology.StablePartition import (
    _near_far, _row_at, _row_polyline, _stable_frame, _side_of, row_of_end, rows_of_bridge,
)

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from _builds import build_k28  # noqa: E402

logging.basicConfig(level=logging.WARNING)


def geometric_row(trellis, point):
    """The side of the stable manifold (looking anchorward) the point sits on, nearest segment."""
    best = None
    for key, manifold in trellis.manifolds.items():
        if key[1] != "stable":
            continue
        pts = manifold.get_point_array()
        if pts is None or len(pts) < 2:
            continue
        pts = np.asarray(pts, dtype=float)
        a, b = pts[:-1], pts[1:]
        ab = b - a
        t = np.clip(np.einsum("ij,ij->i", point - a, ab) / np.maximum(np.einsum("ij,ij->i", ab, ab), 1e-300), 0, 1)
        proj = a + t[:, None] * ab
        d = np.linalg.norm(proj - point, axis=1)
        i = int(np.argmin(d))
        if best is None or d[i] < best[0]:
            # manifold point arrays run root(anchor)->tail, so anchorward = a - b
            best = (float(d[i]), key, _side_of(a[i] - b[i], point - proj[i]))
    return best


def bridge_by_ends(trellis):
    """Map frozenset(endpoint ids) -> Bridge."""
    return {frozenset((b.first_intersection, b.second_intersection)): b for b in trellis.bridges if not b.partial}


def report(name, trellis):
    """Print every hole and return the per-origin rows."""
    print(f"\n===== {name}: {len(trellis.holes)} holes, orientation_preserving={trellis.orientation_preserving}")
    ends = bridge_by_ends(trellis)
    per_origin = defaultdict(list)
    for hole in sorted(trellis.holes, key=lambda h: (str(h.origin), -(h.iterate or 0))):
        near, far = _near_far(trellis, *hole.bounding_ids)
        bridge = ends.get(frozenset(hole.bounding_ids))
        poly, _ = _row_polyline(trellis, bridge)
        rows_geo = {iid: _row_at(trellis, poly, iid) for iid in (near, far)}
        try:
            rows_comb = rows_of_bridge(trellis, bridge)
        except Exception as exc:  # noqa: BLE001
            rows_comb = f"err {exc!r}"
        nix, fix = trellis.intersection(near), trellis.intersection(far)
        geo = geometric_row(trellis, np.asarray(hole.coords, dtype=float))
        opened_ids = {iid for iid, _w, _r in hole.openings}
        missing = []
        for iid in (near, far):
            if iid in opened_ids:
                continue
            _aw, _ow, below, above = _stable_frame(trellis, trellis.intersection(iid))
            missing.append((iid, "anchor(no node below)" if below is None else ("tail(no node above)" if above is None else "both halves present: faces other branch/side test failed")))
        print(
            f"origin={hole.origin} it={hole.iterate:+d} direct={hole.pair is not None} side={hole.bridge_side} "
            f"bounds near={near}(s={nix.stable_cdist:.3g},label={nix.label},b={nix.manifold_b_key[2:] if nix.manifold_b_key else None}) "
            f"far={far}(s={fix.stable_cdist:.3g},b={fix.manifold_b_key[2:] if fix.manifold_b_key else None}) "
            f"bridge_u={bridge.manifold_key[2:] if bridge and bridge.manifold_key else None}"
        )
        print(f"     openings={hole.openings}")
        print(f"     row@bound geo={rows_geo} comb={rows_comb}  hole geometric stable-side={geo[2]} (dist {geo[0]:.2e}, branch {geo[1][2:]})")
        if missing:
            print(f"     missing={missing}")
        per_origin[hole.origin].append((hole.iterate, hole.bridge_side, sorted((w, r) for _i, w, r in hole.openings), geo[2]))
    return per_origin


def verdict(per_origin):
    """Compare each backward image's rows against its origin's direct hole rows."""
    flips = []
    for origin, rows in per_origin.items():
        direct = [r for r in rows if r[0] >= 0]
        if not direct:
            continue
        ref_rows = {r for _w, r in direct[-1][2]} if direct else set()
        for it, side, ops, geo in rows:
            if it >= 0:
                continue
            rows_here = {r for _w, r in ops}
            if rows_here and rows_here != ref_rows:
                flips.append((origin, it, sorted(ref_rows), sorted(rows_here), geo))
    return flips


if __name__ == "__main__":
    build = build_period3()
    (fp3,) = build.fixed_points
    p3 = report("p3 (build_period3 defaults)", build.session.trellis(fp3))
    session, fp = build_k28(2)
    k28 = report("k28 two blasts", session.trellis(fp))
    for name, po in (("p3", p3), ("k28_two_blasts", k28)):
        print(f"\n{name} row flips (origin, iterate, origin rows, image rows, image geometric side):")
        for f in verdict(po):
            print("   ", f)


def linked_bound_law(trellis):
    """Check: a propagated hole's bound that is the registered iterate of its
    origin's direct-hole bound opens the same (which, row) there.

    Returns (checked, violations, unlinked) counts/lists.
    """
    direct = {}
    for hole in trellis.holes:
        if hole.pair is not None and hole.iterate == 0:
            direct[hole.origin] = hole
    checked, violations, unlinked = 0, [], 0
    for hole in trellis.holes:
        if hole.pair is not None or hole.iterate is None or hole.iterate >= 0:
            continue
        ref = direct.get(hole.origin)
        if ref is None:
            continue
        ref_open = {iid: (w, r) for iid, w, r in ref.openings}
        here_open = {iid: (w, r) for iid, w, r in hole.openings}
        for ref_id in ref.bounding_ids:
            img = trellis.iterate(ref_id, hole.iterate)
            if img is None or img not in hole.bounding_ids:
                unlinked += 1
                continue
            if ref_id in ref_open and img in here_open:
                checked += 1
                if ref_open[ref_id] != here_open[img]:
                    violations.append((hole.origin, hole.iterate, ref_id, img, ref_open[ref_id], here_open[img]))
    return checked, violations, unlinked


def law_sweep():
    """Run the linked-bound law on every law case that builds."""
    from tanglepack.examples.henon_cases import build_nested
    from _builds import build_k10
    cases = {}
    s, fp = build_k10(); cases["k10"] = [s.trellis(fp)]
    s, fp = build_k28(1); cases["k28_one_blast"] = [s.trellis(fp)]
    s, fp = build_k28(2); cases["k28_two_blasts"] = [s.trellis(fp)]
    b = build_period3(); cases["p3"] = [b.session.trellis(f) for f in b.fixed_points]
    b = build_nested(); cases["nested"] = [b.session.trellis(f) for f in b.fixed_points]
    print("\nlinked-bound law (propagated hole at a registered iterate of its origin's bound keeps (which,row)):")
    for name, trellises in cases.items():
        tot = [0, [], 0]
        for T in trellises:
            c, v, u = linked_bound_law(T)
            tot[0] += c; tot[1] += v; tot[2] += u
        print(f"   {name}: checked={tot[0]} violations={tot[1]} unlinked_bounds={tot[2]}")


if __name__ == "__main__" and "--law" in sys.argv:
    law_sweep()
