# Higher-Period Bridge Iteration Bug — Root Cause & Fix Plan

## The crash

```
File "ManifoldMachine.py", line 872, in _get_refined_point
    p1_preiterate.get_point() + p0_preiterate.get_point()
AttributeError: 'NoneType' object has no attribute 'get_point'
```

`p1_preiterate` (or `p0_preiterate`) is `None`, meaning `point.prev_iterate is None` for some
point in the manifold being refined during `iterate_bridge`.

---

## Root cause: `create_bridges` mixes segments from different orbit branches

### What happens for period-1

For a period-1 fixed point there is exactly one orbit branch. Every intersecting unstable segment
belongs to the same manifold. `create_bridges` sorts them by cdist and creates bridges between
consecutive pairs. Root and tail of every bridge land on the **same** geometric list.

### What breaks for period-N (N > 1)

For a period-3 fixed point there are three orbit branches, each stored as a separate
`BaseManifold` in `wb.manifolds` with keys `(fp, "unstable", 0, 0)`, `(fp, "unstable", 1, 0)`,
`(fp, "unstable", 2, 0)`.

`compute_intersections` indexes **all three** branches into the R-tree.  Intersections are
detected between each orbit branch and the stable manifold independently, so the registry may
contain intersections from orbit 0 at cdist 1.4, orbit 1 at cdist 1.5, orbit 2 at cdist 1.3,
and so on — interleaved by cdist value.

`create_bridges` (`Tangle.py:524–598`) then:

1. Collects every intersecting unstable segment into a single flat list.
2. **Sorts the entire list by cdist** — across all orbit branches simultaneously.
3. Creates bridges between every consecutive pair `(seg_i, seg_{i+1})`, regardless of which
   orbit branch each segment belongs to.

When `seg_i` belongs to orbit 0 and `seg_{i+1}` belongs to orbit 1:

```python
seg1.p0.insert_point_forward(root_point, ...)   # root_point → orbit 0 geometric list
seg2.p0_seg1.insert_point_backward(tail_point, ...)  # tail_point → orbit 1 geometric list

bridge = Bridge(root=root_point, ..., tail=tail_point)  # root and tail on DIFFERENT lists
```

`root_point` and `tail_point` are on geometrically disconnected manifolds.  Walking forward
from `root_point` through orbit 0 never reaches `tail_point` on orbit 1.

### How the None propagates to the crash

When `iterate_bridge(bridge)` is called on a cross-orbit bridge:

1. `iterate_manifold(bridge)` calls `get_non_iterated_point_array()`, which walks from
   `root_point` through orbit 0 until `current is None` (never reaching `tail_point`).

2. `old_iterated_points = BaseManifold(old_points[0], ..., tail=old_points[-1])` creates a
   view of orbit 1 (the `next_iterate`s of orbit 0 interior points), walking from
   `old_points[0]` to `old_points[-1]` through orbit 1's **real** geometric list.

3. Because `tail_point` was inserted **into orbit 1's geometric list** by
   `seg2.p0_seg1.insert_point_backward(tail_point, ...)`, it now sits **between**
   `old_points[0]` and `old_points[-1]` in orbit 1.  The walk through `old_iterated_points`
   therefore **includes `tail_point`**.

4. `tail_point` was created fresh in `Tangle._get_nearby_point` with no iterate links:
   `tail_point.prev_iterate = None`.

5. `refine_manifold(mapped_manifold)` walks adjacent pairs and calls
   `_get_refined_point(p_before_tail, tail_point, ...)`.
   `_get_preiterate(tail_point, "unstable", 1)` → `tail_point.prev_iterate` → **`None`** → crash.

The same logic applies to `root_point` if it too ends up inside the walked range.

---

## Secondary issue: dangling tail iterates after the last grow

After `grow_n_times(fp3, "unstable", 7)`, the final call to `new_grow_manifold` ends with
step `i=2` (orbit 2 → orbit 0).  The new orbit 0 tail points created in this step have
`prev_iterate = orbit 2 point` but `next_iterate = None` (there is no 8th grow).

These points appear in `non_iterated_points` during `iterate_bridge` and are mapped forward
correctly.  Their images' `prev_iterate` is set to these orbit 0 tail points, which have
coordinates, so `get_point()` works.  **This does NOT directly cause the current crash** but
it means the cdist ordering assumption in `merge_manifolds` may be violated for those tail
points if they land in an unusual part of the orbit 1 manifold.

---

## All locations that need changes

| File | Line range | Issue |
|---|---|---|
| `Tangle.py` | `524–598` (`create_bridges`) | Mixes orbit-branch segments; creates cross-orbit bridges |
| `Tangle.py` | `601–628` (`_get_nearby_point`) | Boundary points created with no iterate links |
| `ManifoldMachine.py` | `831–843` (commented guard) | Safety guard was removed; needs to be reinstated as a fallback |
| `ManifoldMachine.py` | `856–888` (`_get_refined_point`) | No fallback when `prev_iterate is None` |

---

## Fix 1 (core): make `create_bridges` orbit-aware  — `Tangle.py`

Group intersecting segments **by their parent manifold** before sorting.  Create bridges only
between consecutive intersections on the **same** manifold.

```python
def create_bridges(self, for_manifold=None):
    from collections import defaultdict

    # --- 1. Collect intersecting unstable segments, grouped by manifold ---
    manifold_segs: dict[BaseManifold, list[tuple[float, _Segment]]] = defaultdict(list)
    seen_ids: set[int] = set()

    for sid_pair in self._intersecting_segments:
        if for_manifold is not None and not (sid_pair & self._manifold_segs.get(for_manifold, set())):
            continue
        for sid in sid_pair:
            if sid in seen_ids:
                continue
            seen_ids.add(sid)
            seg = self._seg_lookup[sid]
            if seg.manifold.stability != "unstable":
                continue
            cdist = 0.5 * (
                seg.p0.get_cdist("unstable") + seg.p0_seg1.get_cdist("unstable")
            )
            manifold_segs[seg.manifold].append((cdist, seg))

    # --- 2. Within each manifold, sort by cdist and create bridges ---
    all_bridges: list[Bridge] = []

    for manifold, segs in manifold_segs.items():
        segs.sort(key=lambda x: x[0])

        for i in range(len(segs) - 1):
            _, seg1 = segs[i]
            _, seg2 = segs[i + 1]

            root_point = self._get_nearby_point(seg1, "root")
            tail_point = self._get_nearby_point(seg2, "tail")

            # cache pre-iterates on boundary points (see Fix 2)
            self._cache_boundary_preiterate(root_point, seg1, manifold.stability)
            self._cache_boundary_preiterate(tail_point, seg2, manifold.stability)

            seg1.p0.insert_point_forward(root_point, manifold.branch_index)
            seg1.p0 = root_point
            seg2.p0_seg1.insert_point_backward(tail_point, manifold.branch_index)
            seg2.p0_seg1 = tail_point

            bridge = Bridge(
                root=root_point,
                stability=manifold.stability,
                stretch_param=manifold.stretch_param,
                fixed_point=manifold.fixed_point,
                tail=tail_point,
                branch_index=manifold.branch_index,
            )
            all_bridges.append(bridge)

    # wire next_bridge / prev_bridge doubly-linked list
    for i in range(len(all_bridges) - 1):
        all_bridges[i].next_bridge = all_bridges[i + 1]
        all_bridges[i + 1].prev_bridge = all_bridges[i]

    return all_bridges
```

Key change: `manifold_segs[seg.manifold].append(...)` ensures segments are bucketed by which
orbit branch they live on.  The inner loop never pairs a segment from orbit 0 with one from
orbit 1.

---

## Fix 2: cache pre-iterates on bridge boundary points — `Tangle.py`

Add a helper `_cache_boundary_preiterate` that interpolates the preiterate for a freshly
created boundary point from the known preiterates of its geometric neighbors:

```python
def _cache_boundary_preiterate(
    self,
    boundary_point: Point,
    seg: _Segment,
    stability: str,
) -> None:
    """
    Approximate and cache `prev_iterate` (unstable) or `next_iterate` (stable)
    on a freshly created bridge boundary point.

    Uses the same weighted interpolation as `_linear_interpolation`: the boundary
    point sits alpha=0.1 away from the intersection toward seg_point, so the
    preiterate is interpolated with the same weight.
    """
    p0, p1 = seg.p0, seg.p0_seg1

    if stability == "unstable":
        pre0 = p0.prev_iterate
        pre1 = p1.prev_iterate
    else:
        pre0 = p0.next_iterate
        pre1 = p1.next_iterate

    if pre0 is None or pre1 is None:
        return  # can't interpolate; leave None and rely on guard in refine

    # 0.9 * intersection_preiterate + 0.1 * seg_point_preiterate
    # (mirrors the 0.1 alpha used in _linear_interpolation for the point itself)
    coords = 0.9 * 0.5 * (pre0.get_point() + pre1.get_point()) + 0.1 * pre0.get_point()

    cached = Point(coords[0], coords[1])
    if stability == "unstable":
        boundary_point.prev_iterate = cached
    else:
        boundary_point.next_iterate = cached
```

This gives boundary points an approximate preiterate that is geometrically consistent
with their position, allowing `_get_refined_point` to proceed without crashing.

---

## Fix 3: reinstate the safety guard — `ManifoldMachine.py`

Uncomment and tighten the guard in `refine_two_points` (lines 831–843):

```python
# Skip refinement if either point's preiterate is missing.
if (
    ManifoldMachine._get_preiterate(p0, manifold.stability, 1) is None
    or ManifoldMachine._get_preiterate(p1, manifold.stability, 1) is None
):
    logger.debug(
        "Skipping refinement for pair at cdist %.3g–%.3g: "
        "preiterate missing on one or both points.",
        p0.cdist,
        p1.cdist,
    )
    continue
```

This is a **defensive fallback**, not the root fix.  With Fix 1 in place, this guard should
rarely trigger for correct inputs.  It prevents a crash in any edge case where a boundary
point's preiterate could not be set (e.g., because its neighbors also lacked preiterates and
Fix 2 returned early).

---

## Fix 4: `TangleWorkbench.compute_intersections` must index each orbit separately — `TangleWorkbench.py`

Currently `index_manifolds` adds all orbit branches into the **same** `Tangle._manifold_segs`
dict.  This is fine; the real fix is in `create_bridges`.  No change needed here beyond Fix 1.

However, `plot_intersections` should be verified to still work after Fix 1 since `_intersecting_coords`
is still populated per-segment regardless of orbit.

---

## Fix 5: `iterate_bridge` manifold-key propagation — `TangleWorkbench.py`

After `iterate_bridge`, the child bridges returned by `Tangle.create_bridges(for_manifold=iterated)`
will now be orbit-aware.  But `iterated.manifold_key` may not be set (it's a fresh
`BaseManifold`/`Bridge` returned by `iterate_manifold`).

Ensure `iterate_manifold` propagates `manifold_key`:

```python
# in iterate_manifold, when building new_iterated_points:
new_iterated_points.manifold_key = manifold.manifold_key
```

And in `iterate_bridge` in `ManifoldMachine.py`, propagate the key to the returned manifold:

```python
def iterate_bridge(self, manifold: Bridge):
    iterated_manifold = self.iterate_manifold(manifold)
    iterated_manifold.manifold_key = manifold.manifold_key  # propagate for Tangle lookup
    return iterated_manifold
```

Without this, `Tangle._key_of(seg)` returns `None` for iterated bridge segments, which
breaks `manifold_a_key`/`manifold_b_key` on new intersections (silently — no crash, but
registry lookups by fixed point and branch will miss those intersections).

---

## Fix 6: `new_grow_manifold` — close the dangling-tail iterate gap — `ManifoldMachine.py`

After `grow_n_times(fp, "unstable", N)`, the last orbit in the final grow cycle has new tail
points with `next_iterate = None`.  This does not currently crash, but it means those tail
points cannot participate in the iterate table (graph and `infer_iterate_table` skip them).

**Option A (simple):** After `grow_n_times`, do one extra half-step: iterate only orbit 0
without growing the others.  This ensures every orbit 0 tail point has `next_iterate` set.
Cost: one extra `iterate_manifold` call that refines but doesn't extend orbit 1.

**Option B (architectural):** Redesign `new_grow_manifold` so that all three orbits are
grown symmetrically — grow orbit 0, 1, 2 all by one step before looping.  This ensures
after N grows every orbit is fully connected in the iterate chain.

Option A is simpler and lower-risk.  Implement as a post-grow helper in `TangleWorkbench`:

```python
def _close_iterate_gaps(self, fixed_point: FixedPoint, stability: str) -> None:
    """One extra iterate pass on orbit 0 to wire up dangling tail points."""
    orbit_indices = fixed_point.get_iterable_array(stability, shift=1)
    current_manifold = BaseManifold(
        fixed_point.branch_points[orbit_indices[0]],
        stability, stretch_param=1, fixed_point=fixed_point, branch_index=0,
    )
    temp_root = current_manifold.root
    if isinstance(current_manifold.root, BranchPoint):
        current_manifold.root = current_manifold.walk_fwd(None, temp_root)
    if current_manifold.root is None:
        return
    current_manifold.stretch_param = current_manifold.root.stretch_param
    self._man_machine.iterate_manifold(current_manifold)  # side-effect: sets next_iterate
```

Call this inside `grow_n_times` after the main grow loop.

---

## Implementation order

1. **Fix 3** — uncomment the safety guard (`ManifoldMachine.py:831–843`).
   Zero risk; immediately stops the crash.  Can ship alone.

2. **Fix 1** — orbit-aware `create_bridges` (`Tangle.py`).
   This is the root fix.  Required for topologically correct bridge topology on period-N orbits.

3. **Fix 2** — cache preiterates on boundary points (`Tangle.py`).
   Needed for smooth refinement at bridge edges even after Fix 1, because the
   bridge-boundary/interior boundary pair still involves one point without a real preiterate.

4. **Fix 5** — `manifold_key` propagation (`ManifoldMachine.py`, `TangleWorkbench.py`).
   Low-risk, ensures registry correctness for iterated bridges.

5. **Fix 6** — close iterate gaps after grow (`TangleWorkbench.py`).
   Lower priority; does not crash, but affects iterate-table completeness.

---

## Testing

After fixes, verify:

- `henon_new_intersection_structures.py` (period-1, k=10) still runs without regression.
- `henon_period_3_new_intersection_structures.py` (period-3, k=2) completes without crash.
- Each bridge returned by `create_bridges` has its `root` and `tail` on the **same** manifold
  (checkable by asserting `bridge.root` is reachable from `bridge.root` by forward walk).
- `registry.by_unstable_cdist` is non-empty and `manifold_a_key` is set on all intersections.
- `build_intersection_graph()` produces a connected graph with iterate edges.

Add a unit test asserting bridge connectivity:
```python
def _bridge_is_connected(bridge: Bridge) -> bool:
    current = bridge.root
    prev = None
    while current is not None:
        if current is bridge.tail:
            return True
        nxt = bridge.walk_fwd(prev, current)
        prev, current = current, nxt
    return False
```

---

## Summary

The crash is not a simple off-by-one — it's a structural mismatch: **`create_bridges` was
written for period-1 (single orbit) and silently produces geometrically invalid bridges for
period-N (N>1) orbits by mixing segments from different orbit branches in a global cdist sort.**
The downstream crash in `_get_refined_point` is a symptom of that invalid bridge structure,
not the root problem.
