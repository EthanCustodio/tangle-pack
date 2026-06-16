# Nested-Tangle Intersection Bug — Root Cause & Fix Plan

Debugging target: `scripts/henon_period_3_new_intersection_structures.py`
(k=2, b=1 Hénon map; inner period-3 fixed point `fp3` nested inside outer period-1 `fp1`.)

This builds on `higher_period_bridge_bug_plan.md`, which fixed the cross-orbit bridge
crash. That fix (orbit-aware `create_bridges`) is in place and working. The problems
described here are **different** and were exposed once the bridge crash was resolved.

---

## Symptoms

1. **A couple of intersections are missing in the *outer* period-1 tangle figure** —
   most are present, a few are not drawn.
2. **The *inner* period-3 tangle shows no intersections at all.**
3. **No heteroclinic (inner↔outer) intersections are computed.**

The user's hypothesis for (1): *"my R-tree implementation does not allow a single
segment to be a member of two different intersection pairs; this should be allowed."*
That hypothesis is essentially correct, but the limitation is **not** in the R-tree or in
`_intersecting_segments` (which is a `set[frozenset[int]]` and already stores a segment in
many pairs). It lives in the **per-segment lookup dicts** and in the **per-segment cdist**.
Details below.

---

## Root cause #2 & #3 (the big one): `compute_intersections` is single-fixed-point and destroys prior state

```python
# TangleWorkbench.compute_intersections (TangleWorkbench.py:291)
def compute_intersections(self, fixed_point):
    self.Tangle.clear_all()                               # ← wipes EVERYTHING
    self._intersection_registry = IntersectionRegistry()  # ← throws away the registry
    self.index_manifolds(fixed_point, "unstable")         # ← only THIS fixed point
    self.index_manifolds(fixed_point, "stable")
    self.Tangle.populate_intersection_dict()
    for intersection in self.Tangle._intersections:
        self._intersection_registry.add(intersection)
```

The script does:

```python
wb.compute_intersections(fp3)   # computes fp3 homoclinic, builds registry
wb.trim_stable_manifolds(fp3)

wb.compute_intersections(fp1)   # clear_all() DESTROYS every fp3 segment + intersection,
wb.trim_stable_manifolds(fp1)   # resets the registry, indexes ONLY fp1

bridges  = wb.create_bridges(fp3)   # Tangle now holds only fp1 → no fp3 bridges
bridges1 = wb.create_bridges(fp1)
```

Consequences, exactly matching the symptoms:

- **Symptom 2** (inner period-3 has no intersections): they *were* computed by the first
  call, then `clear_all()` in the second call deleted every fp3 segment from the R-tree and
  `IntersectionRegistry()` discarded the fp3 intersections. By the time `create_bridges(fp3)`
  runs, the Tangle contains only fp1 geometry. fp3 has nothing.

- **Symptom 3** (no heteroclinic intersections): `compute_intersections(fp)` only ever indexes
  **one** fixed point's manifolds. fp3's manifolds and fp1's manifolds are **never in the
  R-tree at the same time**, so a crossing between fp3's unstable manifold and fp1's stable
  manifold can never be detected — there is no moment when both segments coexist in the index.

The detection machinery itself is fine for cross-fixed-point crossings: `_insert_segment`
(`Tangle.py:420`) and `intersections_for_segment` only skip pairs where
`other.manifold is seg.manifold`, so two **different** manifolds (different fp, or different
orbit branch) would be paired correctly — *if they were ever indexed together.*

### There is also a secondary correctness gap once everything is co-indexed

**The fundamental invariant of the system (see `CLAUDE.md`):** two unstable manifolds can
*never* intersect, and two stable manifolds can *never* intersect — whether they belong to the
same fixed point or to different fixed points. A shared point of two stable manifolds would have
to converge forward to *both* fixed points at once; a shared point of two unstable manifolds
would have to converge backward to both. Only an unstable manifold may cross a stable one
(homoclinic within one fixed point, heteroclinic across two). So **co-indexing fp3 and fp1 does
not introduce any new *real* same-stability crossings** — none can exist geometrically.

The gap is therefore not that "u×u crossings now appear and must be modeled," but that the
detection layer is purely geometric on *polygonal* approximations. `_insert_segment`
(`Tangle.py:420-430`) only skips pairs where `other.manifold is seg.manifold`; it registers any
straddling segment pair into `_intersecting_segments`. Two distinct same-stability manifolds that
run near-parallel and very close (e.g. near a tangency, or where two stable manifolds crowd
toward a fixed point) can produce a **spurious** straddle of their straight-line segments even
though the true smooth curves do not cross. Fed through `populate_intersection_dict`
(`Tangle.py:189-197`), such an artifact pair mis-assigns cdists:

```python
unstable_cdist = seg_1_cdist if seg_1.manifold.stability == "unstable" else seg_2_cdist
stable_cdist   = seg_1_cdist if seg_1.manifold.stability == "stable"   else seg_2_cdist
```

For a u×u artifact this stores one segment's unstable cdist *as if it were a stable cdist* —
garbage. With a single homoclinic tangle the artifact had no second same-stability manifold to
straddle against, so it never bit; with nested co-indexed tangles a near-parallel pair can.
**Every pair must be filtered to exactly one unstable + one stable segment** before becoming an
`Intersection`. This is a defensive guard against numerical noise — *not* the handling of a real
geometric case — so a same-stability pair that fires the filter should be logged at
`logger.debug`/`warning` level (it flags a near-tangency or under-resolved manifold), never
silently modeled as a crossing.

---

## Root cause #1: per-segment storage collapses multiple crossings on one segment

Detection is complete — `_intersecting_segments` holds `{A, B}` and `{A, C}` as distinct
frozensets when one coarse segment `A` crosses two segments `B` and `C`. The loss happens in
`populate_intersection_dict` (`Tangle.py:213-216`), which writes into dicts **keyed by a
single segment id**:

```python
self._intersecting_coords[seg1_id] = point   # dict[int, coord]
self._intersecting_coords[seg2_id] = point
self._intersecting_points[seg1_id] = branch_point
self._intersecting_points[seg2_id] = branch_point
# and earlier: self._intersection_by_seg[seg1_id] = intersection
```

When segment `A` participates in two crossings, the second overwrites the first under key `A`.

Why only "*a couple*" disappear from the figure (not half of them): each crossing is stored
under **both** of its segment ids, so a crossing is only lost from `plot_intersections`
(which draws `_intersecting_coords.values()`, `TangleWorkbench.py:314`) when **both** of its
host segments get overwritten by *other* crossings. That happens only in tight folds /
near-tangencies where several crossings cluster onto a few coarse segments — i.e. exactly the
busiest parts of a dense outer tangle. Hence a *few* points vanish rather than a systematic
fraction.

### A related accuracy bug: per-segment-midpoint cdist

`populate_intersection_dict` computes each crossing's cdist as the **midpoint of the segment**:

```python
seg_1_cdist = 0.5 * (seg_1.p0.get_cdist(...) + seg_1.p0_seg1.get_cdist(...))
```

Two crossings sharing unstable segment `A` therefore receive the **same** `unstable_cdist`.
This:
- corrupts `by_unstable_cdist` ordering (ties at identical cdist),
- makes `nearest_by_unstable_cdist` (used for bridge-endpoint assignment,
  `TangleWorkbench._assign_bridge_intersections`) ambiguous,
- and can trip `IntersectionRegistry._find_collision` (`IntersectionRegistry.py:555`) into
  treating two genuinely distinct crossings as one if their stable cdists also happen to be
  within `cdist_tol`.

The true intersection coordinate is already computed (`_find_true_intersection`); the cdist
should be interpolated **at that true point**, not taken as the segment midpoint.

---

## Affected locations

| File | Location | Issue |
|---|---|---|
| `TangleWorkbench.py` | `compute_intersections` (291) | `clear_all()` + single-fp index + registry reset; cannot co-index fixed points |
| `TangleWorkbench.py` | `index_manifolds` (282) | No "index everything" path |
| `TangleWorkbench.py` | `create_bridges` (327) | Ignores its `fixed_point` arg; always bridges all manifolds → double-adds when called per-fp |
| `Tangle.py` | `populate_intersection_dict` (144) | Per-segment dicts collapse multiple crossings; segment-midpoint cdist; no u×s filter |
| `Tangle.py` | `_intersecting_coords` / `_intersecting_points` / `_intersection_by_seg` | Keyed by single seg id |
| `Tangle.py` | `_get_nearby_point` (666) | Reads `_intersecting_coords[seg.id]` — assumes one crossing per segment |
| `Tangle.py` | `populate_intersections_for_manifold` (298) | Same per-segment-dict + midpoint-cdist issues for iterated bridges |

---

## Fix plan

Implement in this order. Fix A unblocks symptoms 2 & 3. Fix B unblocks symptom 1. Fix C makes
`create_bridges` correct once A co-indexes multiple manifolds and adds the per-fixed-point
filter you asked for. Every block below is **complete drop-in code** with its exact insertion
point; nothing is pseudocode.

A note on bridge ↔ fixed-point linkage: each `Bridge` is constructed with
`fixed_point=manifold.fixed_point` (it lives on `BaseManifold`, `Bridge.py:33-65`). So a single
global `create_bridges()` call already tags every bridge with the correct fixed point — calling
it globally does **not** lose that information. Fix C also adds an optional `fixed_point=`
filter so you can build the two tangles' bridges separately when you want to.

---

### Fix A — make intersection computation cumulative and multi-fixed-point

**A1. Replace `TangleWorkbench.compute_intersections` (`TangleWorkbench.py:291-307`) in full.**

```python
def compute_intersections(self, fixed_points, *, reset: bool = True):
    """
    Compute intersections among the manifolds of one or more fixed points.

    All manifolds of every supplied fixed point are indexed into the SAME Tangle
    before crossings are resolved, so homoclinic crossings (within one fixed
    point) and heteroclinic crossings (between two fixed points) are detected
    together.

    Args:
        fixed_points: A single FixedPoint or an iterable of FixedPoints whose
            manifolds should be co-indexed and intersected.
        reset: If True (default) the Tangle and registry are cleared first. Pass
            False to accumulate further manifolds into an existing computation.

    Returns:
        List of (x, y) coordinates, one per detected crossing.
    """
    if isinstance(fixed_points, FixedPoint):
        fixed_points = [fixed_points]

    if reset:
        self.Tangle.clear_all()
        self._intersection_registry = IntersectionRegistry()

    for fp in fixed_points:
        self.index_manifolds(fp, "unstable")
        self.index_manifolds(fp, "stable")

    self.Tangle.populate_intersection_dict()

    for intersection in self.Tangle._intersections:
        self._intersection_registry.add(intersection)

    return self.Tangle.iter_intersection_coords()
```

`FixedPoint` is already imported (`TangleWorkbench.py:16`), so the `isinstance` guard works and
preserves the old single-fixed-point call signature for the period-1 regression script.

**A2. Update the script** `scripts/henon_period_3_new_intersection_structures.py` (lines 56-63)
to a single co-indexed compute, then trim, then bridge:

```python
# ── Intersection phase (homoclinic + heteroclinic, both tangles) ───────────
wb.compute_intersections([fp3, fp1])
wb.trim_stable_manifolds(fp3)
wb.trim_stable_manifolds(fp1)

bridges = wb.create_bridges(fp3)
bridges1 = wb.create_bridges(fp1)
```

Trim **after** the single combined compute (not between two computes). `trim_stable_manifolds`
already filters to each manifold's own segments, so trimming the two fixed points' stable
manifolds independently still works with both indexed at once.

---

### Fix B — store crossings per pair, not per segment; interpolate cdist at the true point

**B0. Add a module logger to `Tangle.py`** (the same-stability filter in B3/B4 logs through it).
`Tangle.py` currently has no logger. Add to the import header (`Tangle.py:1-6`):

```python
import logging
```

and just below the imports / Dev Notes block (`Tangle.py:19`):

```python
logger = logging.getLogger(__name__)
```

This matches the `logging`-over-`print` convention in `CLAUDE.md`; the library leaves handler
configuration to the application (`ManifoldMachine.py` already attaches a `NullHandler`).

**B1. Change the dict declarations in `Tangle.__init__` (`Tangle.py:109-114`).** Replace:

```python
        self._intersecting_segments: set[frozenset[int]] = set()
        self._intersecting_coords: dict[int, tuple[float, float]] = {}
        self._intersecting_points: dict[int, BranchPoint] = {}

        self._intersections: list[Intersection] = []
        self._intersection_by_seg: dict[int, Intersection] = {}
        self._processed_pairs: set[frozenset[int]] = set()
```

with (now keyed by the segment **pair**; `_intersection_by_seg` becomes a list per segment so a
segment that hosts several crossings keeps all of them):

```python
        self._intersecting_segments: set[frozenset[int]] = set()
        # keyed by the crossing's segment pair, so one segment can host many crossings
        self._intersecting_coords: dict[frozenset[int], tuple[float, float]] = {}
        self._intersecting_points: dict[frozenset[int], BranchPoint] = {}

        self._intersections: list[Intersection] = []
        self._intersection_by_seg: defaultdict[int, list[Intersection]] = defaultdict(list)
        self._processed_pairs: set[frozenset[int]] = set()
```

`clear_all` (`Tangle.py:137-141`) already calls `.clear()` on each of these, which keeps the
`defaultdict` a `defaultdict` — no change needed there. (`defaultdict` and `Intersection` are
already imported at the top of `Tangle.py`.)

**B2. Add two helpers to `Tangle`** — an iterator for plotting and the true-point cdist. Insert
them just above `_get_nearby_point` (`Tangle.py:666`):

```python
    def iter_intersection_coords(self) -> list[tuple[float, float]]:
        """Return the (x, y) of every detected crossing, exactly one per crossing."""
        return list(self._intersecting_coords.values())

    def _cdist_at_point(self, seg: _Segment, point: np.ndarray) -> float:
        """
        Interpolate the cdist of `seg`'s manifold at the true crossing `point`.

        Uses the fractional projection of `point` onto the segment, so two
        crossings that share one segment receive distinct cdists (unlike the old
        segment-midpoint value).
        """
        a = seg.p0.get_point()
        b = seg.p0_seg1.get_point()
        ab = b - a
        denom = float(ab @ ab)
        t = 0.0 if denom == 0.0 else float((np.asarray(point) - a) @ ab) / denom
        t = min(1.0, max(0.0, t))
        ca = seg.p0.get_cdist(seg.manifold.stability)
        cb = seg.p0_seg1.get_cdist(seg.manifold.stability)
        return (1 - t) * ca + t * cb
```

**B3. Replace `populate_intersection_dict` (`Tangle.py:144-218`) in full** — per-pair storage,
true-point cdist, and a u×s filter:

```python
    def populate_intersection_dict(self):
        """
        Resolve every detected crossing pair into an Intersection.

        Each crossing is keyed by its segment PAIR (frozenset of the two segment
        ids), so a single segment may participate in many crossings without any
        being overwritten. Only unstable x stable pairs are kept. A same-stability
        pair (u x u or s x s) is geometrically impossible (see CLAUDE.md's
        fundamental invariant) -- if one appears it is a polygonal/numerical
        artifact of two near-parallel manifolds straddling, so it is logged and
        discarded, never turned into an Intersection.
        """
        for seg_id_pair in list(self._intersecting_segments):

            if seg_id_pair in self._processed_pairs:
                continue  # already processed this pair

            seg1_id, seg2_id = tuple(seg_id_pair)

            if seg1_id not in self._seg_lookup or seg2_id not in self._seg_lookup:
                # drop stale pair and skip
                self._intersecting_segments.discard(seg_id_pair)
                continue

            seg_1, seg_2 = self._seg_lookup[seg1_id], self._seg_lookup[seg2_id]

            # A real crossing is always one unstable + one stable segment. A
            # same-stability pair cannot be a true crossing (fundamental
            # invariant); it is a near-tangency artifact -- log and drop it.
            if seg_1.manifold.stability == seg_2.manifold.stability:
                logger.warning(
                    "Discarding impossible %s x %s segment pair %s as a "
                    "numerical artifact (near-tangency / under-resolved manifold)",
                    seg_1.manifold.stability,
                    seg_2.manifold.stability,
                    tuple(seg_id_pair),
                )
                self._processed_pairs.add(seg_id_pair)
                continue

            point = self._find_true_intersection(seg_1, seg_2)

            # cdist interpolated at the true crossing (not the segment midpoint)
            seg_1_cdist = self._cdist_at_point(seg_1, point)
            seg_2_cdist = self._cdist_at_point(seg_2, point)

            unstable_cdist = (
                seg_1_cdist if seg_1.manifold.stability == "unstable" else seg_2_cdist
            )
            stable_cdist = (
                seg_1_cdist if seg_1.manifold.stability == "stable" else seg_2_cdist
            )

            branch_point = BranchPoint(
                2, (unstable_cdist, stable_cdist), point[0], point[1]
            )

            if seg_1.manifold.stability == "unstable":
                manifold_a_key, cdist_a = Tangle._key_of(seg_1), seg_1_cdist
                manifold_b_key, cdist_b = Tangle._key_of(seg_2), seg_2_cdist
            else:
                manifold_a_key, cdist_a = Tangle._key_of(seg_2), seg_2_cdist
                manifold_b_key, cdist_b = Tangle._key_of(seg_1), seg_1_cdist

            intersection = Intersection.from_segments(
                coords=tuple(point),
                unstable_cdist=cdist_a,
                stable_cdist=cdist_b,
                seg1_id=seg1_id,
                seg2_id=seg2_id,
                manifold_a_key=manifold_a_key,
                manifold_b_key=manifold_b_key,
            )

            self._intersections.append(intersection)
            self._intersection_by_seg[seg1_id].append(intersection)
            self._intersection_by_seg[seg2_id].append(intersection)

            self._intersecting_coords[seg_id_pair] = tuple(point)
            self._intersecting_points[seg_id_pair] = branch_point

            self._processed_pairs.add(seg_id_pair)
```

**B4. Apply the same three changes to `populate_intersections_for_manifold`** (used by
`iterate_bridge`). Replace `Tangle.py:298-394` in full:

```python
    def populate_intersections_for_manifold(
        self, manifold: BaseManifold
    ) -> list[Intersection]:
        """
        Resolve only crossing pairs that involve a segment from `manifold`.

        Use this after add_manifold() when incrementally adding an iterated bridge
        rather than rebuilding the full tangle from scratch.

        Args:
            manifold: The newly added manifold to resolve intersections for.

        Returns:
            List of newly created Intersection objects.
        """
        manifold_seg_ids = self._manifold_segs.get(manifold, set())
        new_intersections: list[Intersection] = []

        for seg_id_pair in list(self._intersecting_segments):
            if seg_id_pair in self._processed_pairs:
                continue
            # only process pairs where at least one segment belongs to manifold
            if not (seg_id_pair & manifold_seg_ids):
                continue

            seg1_id, seg2_id = tuple(seg_id_pair)  # frozenset unpack

            if seg1_id not in self._seg_lookup or seg2_id not in self._seg_lookup:
                self._intersecting_segments.discard(seg_id_pair)
                continue

            seg_1 = self._seg_lookup[seg1_id]
            seg_2 = self._seg_lookup[seg2_id]

            # A same-stability pair cannot be a true crossing (fundamental
            # invariant); it is a near-tangency artifact -- log and drop it.
            if seg_1.manifold.stability == seg_2.manifold.stability:
                logger.warning(
                    "Discarding impossible %s x %s segment pair %s as a "
                    "numerical artifact (near-tangency / under-resolved manifold)",
                    seg_1.manifold.stability,
                    seg_2.manifold.stability,
                    tuple(seg_id_pair),
                )
                self._processed_pairs.add(seg_id_pair)
                continue

            point = self._find_true_intersection(seg_1, seg_2)

            seg_1_cdist = self._cdist_at_point(seg_1, point)
            seg_2_cdist = self._cdist_at_point(seg_2, point)

            unstable_cdist = (
                seg_1_cdist if seg_1.manifold.stability == "unstable" else seg_2_cdist
            )
            stable_cdist = (
                seg_1_cdist if seg_1.manifold.stability == "stable" else seg_2_cdist
            )

            branch_point = BranchPoint(
                2, (unstable_cdist, stable_cdist), point[0], point[1]
            )
            self._intersecting_coords[seg_id_pair] = tuple(point)
            self._intersecting_points[seg_id_pair] = branch_point

            if seg_1.manifold.stability == "unstable":
                manifold_a_key, cdist_a = Tangle._key_of(seg_1), seg_1_cdist
                manifold_b_key, cdist_b = Tangle._key_of(seg_2), seg_2_cdist
            else:
                manifold_a_key, cdist_a = Tangle._key_of(seg_2), seg_2_cdist
                manifold_b_key, cdist_b = Tangle._key_of(seg_1), seg_1_cdist

            intersection = Intersection.from_segments(
                coords=tuple(point),
                unstable_cdist=cdist_a,
                stable_cdist=cdist_b,
                seg1_id=seg1_id,
                seg2_id=seg2_id,
                manifold_a_key=manifold_a_key,
                manifold_b_key=manifold_b_key,
            )

            self._intersections.append(intersection)
            self._intersection_by_seg[seg1_id].append(intersection)
            self._intersection_by_seg[seg2_id].append(intersection)
            new_intersections.append(intersection)

            self._processed_pairs.add(seg_id_pair)

        return new_intersections
```

**B5. Replace `_get_nearby_point` (`Tangle.py:666-693`)** so the crossing coordinate is passed
in (it can no longer be looked up by a single `seg.id`):

```python
    def _get_nearby_point(
        self,
        seg: _Segment,
        intersection_coords: tuple[float, float],
        side: Literal["root", "tail"],
    ) -> Point:
        """
        Create a Point just inside `seg`, offset 10% from the crossing toward the
        segment endpoint on `side`.

        The crossing coordinate is passed in explicitly (rather than looked up by
        seg.id) so a single segment may host several crossings, each producing its
        own boundary point.
        """
        if side == "root":
            seg_point = seg.p0
        elif side == "tail":
            seg_point = seg.p0_seg1
        else:
            raise ValueError(f"Invalid side: {side}")

        new_point = self._linear_interpolation(
            intersection_coords, seg_point.get_point(), 0.1
        )

        new_cdist = self._cdist_at_point(seg, np.asarray(intersection_coords))

        return Point(new_point[0], new_point[1], new_cdist)
```

**B6. Route plotting through the new iterator.** In `TangleWorkbench.plot_intersections`
(`TangleWorkbench.py:314`), change:

```python
        pts = np.array(list(self.Tangle._intersecting_coords.values()))
```

to:

```python
        pts = np.array(self.Tangle.iter_intersection_coords())
```

(The web dashboard reader `tanglepack_webdash/utils/figures.py:204` consumes `.values()` too;
it still returns one coord per crossing, so it now plots every crossing — no edit required, but
worth a glance when verifying.)

---

### Fix C — `create_bridges`: per-crossing cutting + per-fixed-point filter

**C1. Replace `Tangle.create_bridges` (`Tangle.py:601-664`) in full.** It now iterates
crossings (not deduped segments), groups by the parent unstable manifold, and accepts an
optional `fixed_point` filter:

```python
    def create_bridges(self, for_manifold=None, fixed_point=None):
        """
        Cut every indexed unstable manifold into bridges at its crossings with the
        stable manifold(s).

        Crossings are grouped by their parent unstable manifold and sorted by the
        true crossing cdist, so bridges never span two manifolds and a single
        segment that hosts several crossings is handled correctly (each crossing is
        a distinct cut, not deduped by segment id).

        Args:
            for_manifold: If given, only build bridges for crossings that involve a
                segment of this specific manifold.
            fixed_point: If given, only build bridges on unstable manifolds that
                emanate from this fixed point. Each Bridge still records its own
                fixed_point, so a global call (fixed_point=None) followed by
                filtering on bridge.fixed_point is equivalent to per-fixed-point
                calls.

        Returns:
            List of Bridge objects, doubly linked via next_bridge / prev_bridge.
        """
        from collections import defaultdict

        for_manifold_segs = (
            self._manifold_segs.get(for_manifold, set())
            if for_manifold is not None
            else None
        )

        # --- 1. Collect crossings grouped by their parent unstable manifold ---
        # each entry: (true_unstable_cdist, crossing_coords, unstable_segment)
        manifold_crossings: dict[
            BaseManifold, list[tuple[float, tuple[float, float], _Segment]]
        ] = defaultdict(list)

        for sid_pair in self._intersecting_segments:
            if for_manifold_segs is not None and not (sid_pair & for_manifold_segs):
                continue
            if sid_pair not in self._intersecting_coords:
                continue  # not a resolved unstable x stable crossing

            sid1, sid2 = tuple(sid_pair)
            seg_1, seg_2 = self._seg_lookup[sid1], self._seg_lookup[sid2]

            if seg_1.manifold.stability == "unstable":
                u_seg = seg_1
            elif seg_2.manifold.stability == "unstable":
                u_seg = seg_2
            else:
                continue  # no unstable segment (shouldn't happen post-filter)

            if (
                fixed_point is not None
                and u_seg.manifold.fixed_point is not fixed_point
            ):
                continue

            coords = self._intersecting_coords[sid_pair]
            cdist = self._cdist_at_point(u_seg, np.asarray(coords))
            manifold_crossings[u_seg.manifold].append((cdist, coords, u_seg))

        # --- 2. Within each manifold, sort by cdist and cut consecutive crossings ---
        all_bridges: list[Bridge] = []

        for manifold, crossings in manifold_crossings.items():
            crossings.sort(key=lambda c: c[0])

            for i in range(len(crossings) - 1):
                _, coords1, seg1 = crossings[i]
                _, coords2, seg2 = crossings[i + 1]

                root_point = self._get_nearby_point(seg1, coords1, "root")
                tail_point = self._get_nearby_point(seg2, coords2, "tail")

                # cache pre-iterates on boundary points
                self._cache_boundary_preiterate(root_point, seg1, manifold.stability)
                self._cache_boundary_preiterate(tail_point, seg2, manifold.stability)

                # insert boundary points just inside each crossing
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

**C2. Pass the filter through `TangleWorkbench.create_bridges` (`TangleWorkbench.py:327-332`).**
Replace in full:

```python
    def create_bridges(self, fixed_point: Optional[FixedPoint] = None):
        """
        Cut indexed unstable manifolds into bridges.

        Args:
            fixed_point: If given, only build bridges for that fixed point's
                unstable manifolds. If None, build bridges for every indexed
                unstable manifold at once. Each bridge keeps its own fixed_point
                linkage either way.

        Returns:
            The newly created bridges.
        """
        bridges = self.Tangle.create_bridges(fixed_point=fixed_point)
        self._bridges.extend(bridges)
        self._assign_bridge_intersections(bridges)
        return bridges
```

Now `wb.create_bridges(fp3)` and `wb.create_bridges(fp1)` build the two tangles' bridges
separately with no double-counting, and `wb.create_bridges()` (no arg) builds both at once —
your choice per call site.

> **Known remaining edge case (document, don't block on it):** if two *consecutive* crossings
> on one manifold fall on the **same** unstable segment object, the `seg.p0 = root_point` /
> `seg.p0_seg1 = tail_point` reassignments mutate that shared segment, so a third crossing
> reusing it interpolates from the moved endpoint. This only occurs when a single coarse
> segment carries ≥2 crossings that also become bridge boundaries — rare for well-resolved
> manifolds. If it shows up, the clean fix is to stop mutating `seg.p0` / `seg.p0_seg1` and
> instead anchor each insertion on the original endpoint captured at collection time. Leave the
> reassignments in for now to preserve parity with the working period-1 path.

---

## Verification

After the fixes, re-run and check:

1. **Outer tangle, no missing points.** `henon_period_3_new_intersection_structures.py` figure
   shows every fp1 crossing. Cross-check: number of points drawn ==
   `len(wb.Tangle._intersecting_coords)` (now keyed by crossing pair, so one entry per crossing)
   == number of unstable×stable pairs in `_intersecting_segments`.
2. **Inner tangle present.** `registry.from_fixed_point(fp3)` is non-empty; fp3 intersections
   are drawn by `plot_intersections(fp3)`.
3. **Heteroclinic present.** There exist registered intersections whose `manifold_a_key[0] is fp3`
   and `manifold_b_key[0] is fp1` (or vice versa) — i.e. `len(ix.fixed_points) == 2`.
   Add a quick assert in the script:
   ```python
   hetero = [ix for _, ix in registry if len(ix.fixed_points) == 2]
   print(f"{len(hetero)} heteroclinic intersections")
   ```
4. **Distinct cdists.** No two distinct crossings share an identical `unstable_cdist` unless
   they are truly the same point (guards against the midpoint-tie regression).
4b. **No same-stability crossings survive.** Every entry in `_intersecting_coords` is a u×s
   pair (`{seg.manifold.stability for seg in pair} == {"unstable", "stable"}`). If the
   `logger.warning` from the filter fires at all, treat it as a *diagnostic* — it means two
   same-stability manifolds came close enough to straddle as polygons, i.e. that region wants a
   finer manifold resolution — not as a normal expected event.
5. **Regression.** `henon_new_intersection_structures.py` (period-1, k=10) still produces the
   same intersection count and bridges as before.
6. **Bridges.** Each bridge's `root`/`tail` are on the same manifold (existing invariant from
   the prior plan) and `first_intersection`/`second_intersection` resolve to distinct ids.

---

## Summary

Two independent defects, both surfacing now that the period-N bridge crash is fixed:

- **Symptoms 2 & 3** are caused by `compute_intersections` being **single-fixed-point and
  destructive**: the second call wipes the first, and no two fixed points are ever co-indexed,
  so the inner tangle is erased and heteroclinic crossings can never form. Fix A makes the
  computation cumulative and multi-fixed-point. The accompanying u×s filter is **not** there to
  model new same-stability crossings (those are geometrically impossible — see the fundamental
  invariant) but to discard the rare polygonal artifact a near-tangency can produce once two
  same-stability manifolds are co-indexed; it logs whenever it fires.

- **Symptom 1** is caused by **per-segment storage** (`_intersecting_coords` etc. keyed by a
  single seg id) collapsing multiple crossings that share a segment, compounded by a
  **segment-midpoint cdist** that ties distinct crossings together. Fix B keys crossings by
  the segment **pair**, interpolates cdist at the true intersection, and exposes a
  per-crossing iterator for plotting. Fix C makes bridge-cutting robust to the
  multiple-crossings-per-segment case.
</content>
</invoke>
