"""P2: the k=10 inversion saddle (-2.3166, 2.3166) through the FULL pipeline.

Runs every stage in order, catching the first exception per stage (later
stages still run when possible), then evaluates the physical laws on whatever
completed, and counts anchors per unstable branch (the author's rule: exactly
one; today 2 are expected on an inversion point).
"""

from __future__ import annotations

import logging
import sys
import time
import traceback
from collections import Counter

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
sys.path.insert(0, "tests")
from invariants import (  # noqa: E402
    assert_cdist_monotonic, assert_iterate_relation, assert_no_geometric_spikes, assert_one_to_one,
    assert_area_preserved_along_chain,
)
from tanglepack import TangleSession  # noqa: E402
from tanglepack.examples import HENON_K10, henon_jacobian, henon_map, henon_map_inverse, saddle_guesses  # noqa: E402
from tanglepack.topology import check_bridge_rows_consistent, check_holes_share_bridge_side  # noqa: E402

logging.basicConfig(level=logging.ERROR)

RESULTS: list[tuple[str, str, str]] = []


def record(kind: str, name: str, status: str) -> None:
    """Store and print one stage/law outcome."""
    RESULTS.append((kind, name, status))
    print(f"{kind:5s} {name:48s} {status}")


def stage(name, fn):
    """Run one stage; record OK or the exception type and message head."""
    t0 = time.perf_counter()
    try:
        out = fn()
    except Exception as exc:  # noqa: BLE001
        tb = traceback.extract_tb(exc.__traceback__)[-1]
        record("stage", name, f"RAISES {type(exc).__name__}: {str(exc)[:140]} @ {tb.filename.rsplit('/',1)[-1]}:{tb.lineno}")
        return None, False
    record("stage", name, f"ok ({time.perf_counter() - t0:.2f}s)")
    return out, True


def law(name, fn):
    """Run one law check; record PASS (with item count), FAIL or ERROR."""
    try:
        n = fn()
    except AssertionError as exc:
        record("law", name, f"FAIL: {str(exc)[:160]}")
        return
    except Exception as exc:  # noqa: BLE001
        record("law", name, f"ERROR {type(exc).__name__}: {str(exc)[:160]}")
        return
    record("law", name, f"PASS ({n} items)" if n else "VACUOUS (0 items)")


def main() -> None:
    session = TangleSession(henon_map(*HENON_K10), henon_map_inverse(*HENON_K10), henon_jacobian(*HENON_K10))
    fp = session.construct_fixed_point(saddle_guesses(*HENON_K10)["inversion"])
    print(f"fixed point {np.ravel(fp.coordinates[0])}, period={fp.period}, k_value={fp.k_value}, inversion={fp.check_inversion()}, "
          f"lambda_u={np.ravel(fp.unstable_eigenvalues[0])}, lambda_s={np.ravel(fp.stable_eigenvalues[0])}")
    wb = session.workbench
    ok = {}
    _, ok["init"] = stage("initialize_both_manifolds", lambda: session.initialize_both_manifolds(fp))
    _, ok["grow_u"] = stage("grow_n_times unstable 6", lambda: session.grow_n_times(fp, "unstable", num_iterations=6))
    _, ok["grow_s"] = stage("grow_n_times stable 5", lambda: session.grow_n_times(fp, "stable", num_iterations=5))
    _, ok["intersect"] = stage("compute_intersections", lambda: session.compute_intersections([fp]))
    _, ok["trim"] = stage("trim_stable_manifolds", lambda: session.trim_stable_manifolds(fp))
    _, ok["bridges"] = stage("create_bridges", lambda: session.create_bridges(fp))
    _, ok["infer"] = stage("infer_iterate_table", lambda: session.infer_iterate_table())
    _, ok["pips"] = stage("classify_strong_pips", lambda: session.classify_strong_pips())
    _, ok["pn"] = stage("compute_pseudoneighbors", lambda: session.compute_pseudoneighbors())
    _, ok["holes"] = stage("punch_holes", lambda: session.punch_holes())
    _, ok["partition"] = stage("partition_stable_manifold", lambda: session.partition_stable_manifold())
    arr, ok["arrangement"] = stage("arrangement", lambda: session.arrangement())
    classes, ok["classes"] = stage("bridge_classes", lambda: session.bridge_classes())
    minimal, ok["minimal"] = stage("minimal_trellis", lambda: session.minimal_trellis())
    iterated, ok["iterated"] = stage("iterated_partition", lambda: session.iterated_partition())
    dual, ok["dual"] = stage("dual_graph", lambda: session.dual_graph())
    dyn, ok["symbolic"] = stage("symbolic_dynamics", lambda: session.symbolic_dynamics())

    T = session.trellis(fp)
    reg = wb.intersection_registry
    print(f"\ncounts: crossings={len(reg.all_ids())}, bridges={len(wb.bridges)}, holes={len(T.holes)}, "
          f"pips={T.strong_pip}, pseudoneighbors={len(T.pseudoneighbors)}, partitions={len(T.stable_partitions)}")

    # ---- anchors ----
    anchors = [reg[i] if hasattr(reg, "__getitem__") else T.intersection(i) for i in reg.all_ids()]
    anchors = [T.intersection(i) for i in reg.all_ids() if T.intersection(i).label == "anchor"]
    per_u = Counter(a.manifold_a_key[2:] for a in anchors)
    per_pair = Counter((a.manifold_a_key[2:], a.manifold_b_key[2:]) for a in anchors)
    print(f"\nanchors: total={len(anchors)} per unstable branch={dict(per_u)} per (u,s) pair={dict(per_pair)}")
    law("one anchor per unstable branch", lambda: (_ for _ in ()).throw(AssertionError(f"per branch {dict(per_u)}")) if any(v != 1 for v in per_u.values()) else len(per_u))

    # ---- manifolds ----
    mans = [(k, m) for k, m in wb.manifolds.items() if k[0] is fp]
    def manifolds_law(check):
        def run():
            for _k, m in mans:
                check(m)
            return len(mans)
        return run
    law("manifold cdist monotone (non-strict)", manifolds_law(lambda m: assert_cdist_monotonic(m, strict=False)))
    law("manifold cdist strictly monotone (info)", manifolds_law(lambda m: assert_cdist_monotonic(m, strict=True)))
    law("manifold no geometric spikes", manifolds_law(assert_no_geometric_spikes))
    law("manifold iterate law", manifolds_law(assert_iterate_relation))
    law("manifold one-to-one", manifolds_law(assert_one_to_one))

    # ---- crossings ----
    ids = reg.all_ids()
    def keys_law():
        for i in ids:
            x = T.intersection(i)
            assert x.manifold_a_key[1] == "unstable" and x.manifold_b_key[1] == "stable", (i, x.manifold_a_key, x.manifold_b_key)
        return len(ids)
    law("every crossing is unstable x stable", keys_law)
    def area_law():
        n = 0
        for i in ids:
            x = T.intersection(i)
            if x.label == "anchor" or T.iterate(i, 1) is None:
                continue
            assert_area_preserved_along_chain(reg, i)
            n += 1
        return n
    law("area product preserved along iterate chains", area_law)

    # ---- bridges ----
    bridges = [b for b in wb.bridges if b.manifold_key is not None and b.manifold_key[0] is fp]
    def order_law():
        n = 0
        for b in bridges:
            if b.partial:
                assert b.id is None
                continue
            assert b.id is not None
            u1 = T.intersection(b.id[0]).unstable_cdist; u2 = T.intersection(b.id[1]).unstable_cdist
            assert u1 <= u2, (b.id, u1, u2)
            n += 1
        return n
    law("BridgeId in unstable order / partial <=> id None", order_law)
    def own_branch_law():
        n = 0
        for b in bridges:
            if b.partial:
                continue
            for i in b.id:
                assert T.intersection(i).manifold_a_key == b.manifold_key, (b.id, i, T.intersection(i).manifold_a_key[2:], b.manifold_key[2:])
                n += 1
        return n
    law("bridge endpoints on the bridge's own unstable branch", own_branch_law)
    def bridges_at_law():
        n = 0
        for i in ids:
            expect = sorted(b.id for b in bridges if b.id is not None and i in b.id)
            got = sorted(wb.bridges_at(i))
            assert got == expect, (i, got, expect)
            n += 1
        return n
    law("bridges_at exact", bridges_at_law)
    def single_copy_law():
        c = Counter(b.id for b in bridges if b.id is not None)
        dup = {k: v for k, v in c.items() if v > 1}
        assert not dup, dup
        return len(c)
    law("single copy per BridgeId", single_copy_law)

    # ---- holes / partition ----
    law("I1 holes of one origin share bridge_side", lambda: (check_holes_share_bridge_side(T.holes, orientation_preserving=T.orientation_preserving), len(T.holes))[1])
    def i2():
        n = 0
        for b in T.bridges:
            if b.partial:
                continue
            check_bridge_rows_consistent(T, b); n += 1
        return n
    law("I2 bridge rows consistent", i2)
    def owner_law():
        n = 0
        for res in T.stable_partitions.values() if isinstance(T.stable_partitions, dict) else T.stable_partitions:
            for iid, eid in res.element_of_intersection.items():
                owners = [iv for iv in res.intervals if iv.contains_id(iid)] if hasattr(res.intervals[0], "contains_id") else [eid]
                assert len(owners) == 1, (iid, owners)
                n += 1
        return n
    law("partition: one owner per crossing per side", owner_law)
    law("direct hole at iterate >= k_value absent", lambda: (lambda bad: (_ for _ in ()).throw(AssertionError(bad)) if bad else len(T.holes))([h for h in T.holes if h.pair is not None and (h.iterate or 0) >= fp.k_value]))

    # ---- arrangement ----
    if arr is not None:
        def image_law():
            n = 0
            for region in arr.regions[:20]:
                img = arr.image_of(region)
                if img is None:
                    continue
                assert abs(img.area - region.area) <= 1e-3 * max(region.area, 1e-12), (region.area, img.area)
                n += 1
            return n
        print(f"arrangement: {len(arr.regions)} bounded regions")
        law("arrangement image_of preserves area", image_law)
    # ---- classes ----
    if classes is not None:
        def class_law():
            members = Counter()
            for e in classes.entries:
                for m in e.members:
                    members[m.bridge_id] += 1
            dup = {k: v for k, v in members.items() if v > 1}
            assert not dup, dup
            return len(members)
        law("every bridge in exactly one class", class_law)
    if dyn is not None:
        def even_law():
            n = 0
            for cd in dyn.classes.values():
                if cd.itinerary:
                    assert len(cd.itinerary) % 2 == 0
                    n += 1
            return n
        law("itineraries even", even_law)
        law("no virtual symbols (every class spanned)", lambda: (_ for _ in ()).throw(AssertionError(f"virtual {list(dyn.virtual_classes.values())}")) if dyn.virtual_classes else len(dyn.classes))
        law("is_reliable", lambda: (_ for _ in ()).throw(AssertionError("is_reliable False")) if not dyn.is_reliable else 1)
        print("is_reliable:", dyn.is_reliable)
        print(dyn.describe())


if __name__ == "__main__":
    t0 = time.perf_counter()
    main()
    print(f"\nwall {time.perf_counter() - t0:.1f}s")
