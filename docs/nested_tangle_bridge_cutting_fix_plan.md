# Nested-Tangle Bridge Rendering & Cutting — Diagnosis & Patch

Follow-up to `nested_tangle_intersection_fix_plan.md` (which fixed intersection *detection*).
This addresses the report that bridges in the nested run looked wrong — **"one color crosses
multiple intersections"** — while the single-tangle script looked fine.

> **Correction to an earlier draft of this file.** A previous version of this doc concluded the
> visual was "just accumulation overlap, not a bug." That was **wrong**. There was a real bug —
> in the bridge **coloring** — plus a separate metadata bug in `create_bridges`. Both are fixed
> below. All changes are applied and verified in the source.

---

## TL;DR — what was actually wrong

Three findings, in order of how much they explain the symptom:

1. **Coloring bug (the main visual symptom).** `TangleWorkbench.plot_all_bridges` did
   `cm.get_cmap("tab20", n)` with `n = number of bridges`. For `n > 20` that **resamples** the
   colormap, interpolating between tab20's anchor colors, so consecutive bridges receive
   nearly-identical hues. A run of ~4 adjacent bridges then reads as a single colour band — and
   since each bridge spans one intersection gap, that band visibly crosses several intersections.
   This matches your words exactly ("the colors are not correct; one color will cross multiple
   intersections"), and it shows up even in the sparse outer tangle. **Fix: cycle the 20 discrete
   tab20 colors by `index % 20` instead of resampling, so neighbouring bridges are always
   distinct.**

2. **Same-segment metadata corruption in `create_bridges` (real, from the line-696 edge case).**
   Two consecutive crossings landing on the *same* unstable segment is common here (12–17 of
   ~60–73 crossing-bearing segments per manifold). The old per-crossing code mutated the shared
   `_Segment` (`seg.p0 = root_point` / `seg.p0_seg1 = tail_point`), corrupting every later
   crossing that reused it and producing **41 bridges** (38 fp3 / 3 fp1) with degenerate/wrong
   **cdists**. That cdist feeds `nearest_by_unstable_cdist` → bridge↔intersection assignment →
   the graph, so it is worth fixing even though it does not change the drawn curve. **Fix:
   capture the original endpoints once and splice all boundary points in a single ordered pass —
   no mutation.**

3. **Latent crash in `iterate_bridge` (defensive fix).** When an iterated bridge has no new
   crossings, `TangleWorkbench.iterate_bridge` stored the raw `BaseManifold` returned by
   `ManifoldMachine.iterate_bridge` into `self._bridges`; downstream code then accessed
   `.iterated` / `.children` on it and crashed. **Fix: wrap that manifold as a `Bridge`.** (Not
   triggered by the current scripts, but a real bug.)

### What is NOT a bug — and why the bridge geometry was left alone

A bridge's `root`/`tail` are deliberately placed **just outside** their two bounding
intersections (root toward the lower endpoint, tail toward the higher). This **bracketing is by
design** — `Bridge`'s own docstring: *"the root and tail … we want every bridge to have a point
on either side of the intersection point so that when they are mapped forward we can find the
new intersection."* The forward image of a bracketing bridge straddles the image intersections,
which is how `iterate_bridge` detects them.

I tried "fixing" this by placing boundary points strictly inside each gap. It made the static
picture marginally cleaner but **broke the iterate machinery**: the single-tangle script's
iterate relationships dropped 4 → 0 and its intersection count 12 → 6, because non-bracketing
bridges no longer straddle a crossing when mapped forward. That change was **reverted**. So each
bridge does touch its two endpoint intersections by design — with distinct colors (fix #1) that
reads correctly as "this bridge spans the gap between these two intersections," not as one colour
running through many.

---

## Evidence

- **Coloring**: re-rendering the inner/outer tangle with a 20-colour cycle makes consecutive
  bridges visibly distinct and the colour clearly changes at each intersection marker; the
  crossing-free outer excursions remain single bridges (correct). With the old resampled map the
  same bridges blended into ~6 bands.
- **Edge case pervasiveness** (nested run): crossings vs. segments hosting ≥2 crossings —
  fp3-A 105/73 → 12; fp1 92/72 → 17; fp3-B 81/61 → 9; fp3-C 79/59 → 9. The plan's "rare" remark
  was wrong for this regime.
- **Edge-case corruption**: zero-span (degenerate-cdist) bridges 41 → 0 after the splice fix; no
  inverted bridges; all bridges' `root → tail` walk reaches its tail.
- **Cutting ratio is correct**: bridges = crossings − (#manifolds), i.e. one bridge per
  consecutive-intersection gap (outer alone 85/86; inner alone 258/259; nested 262+91).

---

## Patches (all applied in source)

### Patch 1 — `TangleWorkbench.plot_all_bridges`: cycle discrete colors

```python
        if bridges is None:
            bridges = self._bridges
        n = len(bridges)
        if n == 0:
            return
        # Cycle the 20 discrete tab20 colors by index (mod 20) rather than
        # resampling the colormap across all n bridges. Resampling makes adjacent
        # bridges nearly identical in hue, so a run of consecutive bridges reads as
        # a single colour spanning several intersections; cycling guarantees
        # neighbouring bridges are always visually distinct.
        palette = cm.get_cmap("tab20").colors
        for i, bridge in enumerate(bridges):
            bridge.plot(color=palette[i % len(palette)])
```

### Patch 2 — `Tangle.create_bridges`: original-endpoint capture + single ordered splice

The boundary geometry is unchanged (still the by-design bracketing via `_boundary_point`). The
only change vs. the previous version is that the original segment endpoints are captured at
collection time and all boundary points for a segment are spliced in one ascending-`t` pass, so
a segment hosting several crossings is never corrupted by a prior crossing's insertion. Supporting
helpers added/changed: `_cdist_between`, `_fractional_position`, `_boundary_point` (endpoint-based,
returns `(point, t)`), and `_cache_boundary_preiterate` (now takes `p0, p1` instead of a
`_Segment`). See the source for the full bodies — they are already in place.

Key points of the splice (step 3 of `create_bridges`):

```python
        # Splice every segment's boundary points in one ordered pass. An original
        # segment [p0, p1] is an adjacent pair (p0.forward is p1), so the queued
        # points link straight between them in ascending t.
        for p0, p1, branch_index, inserts in pending_inserts.values():
            inserts.sort(key=lambda x: x[0])
            prev = p0
            for _, point in inserts:
                prev.insert_point_forward(point, branch_index)
                prev = point
```

`branch_index` is threaded through because a `BranchPoint` root requires it; on a plain `Point`
it lands in the `only_forward` slot exactly as the original per-crossing call passed it.

### Patch 3 — `TangleWorkbench.iterate_bridge`: wrap a crossing-free iterate as a Bridge

```python
        # 4. cut at crossings
        if new_intersections:
            new_bridges = self.Tangle.create_bridges(for_manifold=iterated)
        else:
            # No crossings: the iterated manifold is itself one unsplit bridge.
            # ManifoldMachine.iterate_bridge returns a BaseManifold, so wrap it as
            # a Bridge (carrying its manifold_key) — downstream consumers
            # (uniiterated_bridges, infer_iterate_table, genealogy) require Bridge
            # attributes such as .iterated / .children.
            if isinstance(iterated, Bridge):
                new_bridges = [iterated]
            else:
                wrapped = Bridge(
                    root=iterated.root,
                    stability=iterated.stability,
                    stretch_param=iterated.stretch_param,
                    fixed_point=iterated.fixed_point,
                    tail=iterated.tail,
                    branch_index=iterated.branch_index,
                )
                wrapped.manifold_key = getattr(iterated, "manifold_key", None)
                new_bridges = [wrapped]
```

---

## Verification

- **Coloring**: inner-core render now shows distinct alternating bridge colours changing at each
  intersection marker; outer crossing-free lobes stay single bridges.
- **Nested script** end-to-end: 354 intersections, graph 354 nodes / 706 edges, no exceptions.
- **Single-tangle regression**: 12 intersections, **4** iterate relationships, graph 12 / 24 —
  identical to before (confirms the geometry revert restored the iterate machinery).
- **Edge case**: degenerate-cdist bridges 41 → 0; no inverted bridges.

If, after the colour fix, you still see a single bridge run *through* an intersection it should
stop at, that would point back at intersection **detection resolution** in the tightest folds
(coarse polygonal segments stepping over a tangency) — Symptom 1 territory from the previous
plan — not bridge cutting or coloring.
```
