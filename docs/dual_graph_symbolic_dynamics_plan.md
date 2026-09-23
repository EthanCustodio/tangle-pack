# Dual Graph, Bridge Classes, and Symbolic Dynamics — Implementation Plan

> **Retracted in part, 2026-09-14.** The symbolic dynamics did not work
> properly and was removed together with the walk machinery it alone used:
> `topology/SymbolicDynamics.py`, `DualGraph.walk_bridge` (with `Walk`,
> `WalkCase`, `AmbiguousWalkError` and the private walk helpers),
> `Arrangement.face_at_stub`, `plotting.plot_transition_graph`, the
> `TangleSession.symbolic_dynamics` / `plot_transition_graph` wrappers, their
> tests, and `scripts/henon_symbolic_dynamics.py`. Phase A (element identity,
> combinatorial row, bridge classes), Phase B minus the walks (the dual graph:
> arc nodes, merged face nodes, fill, `element_images`, `plot_dual_graph`) and
> Phase C (`image_cdist`, the `image_of_element` fallback) survive.
>
> Two changes landed with the retraction. **Fill:** the strong pip now DRIVES
> the fill instead of being log-checked against a segment derived from the
> outermost crossings — on each pip's own stable branch the arcs meeting
> `(cdist(f^k(q0)), cdist(q0)]` (`k = k_value`) are filled, every other branch
> stays hollow, no pip warns and fills nothing, and `TangleSession.dual_graph`
> keys its cache on the gathered pips too. (A `plot_dual_graph_curved` that
> drew edges along the unstable manifold was tried the same day and dropped;
> the straight-edge `plot_dual_graph` is the one visualisation.) Scripts:
> `scripts/henon_dual_graph.py`
> (both fixtures) and `scripts/henon_k10_dual_graph.py` (chosen pip). Phase D
> and the D-related Deferred items below are historical.

> **A new, different symbolic dynamics landed on 2026-09-21.** It maps each
> bridge class's two homotopy elements one step into the iterated homotopy
> partition and reads the class's word off the shortest dual-graph walk between
> the two landings (both sides recorded at every crossing, so the itinerary
> splits into disjoint pairs, each a class or its inverse; singleton landings
> read the same itinerary off the regular trellis instead), then refines
> classes from the itineraries alone. It revives none of the retracted code:
> `topology/ElementNaming.py`, `topology/DualWalk.py`,
> `topology/SymbolicDynamics.py`, with `TangleSession.symbolic_dynamics` and
> `scripts/henon_symbolic_itineraries.py`. See CLAUDE.md, "Symbolic dynamics".

Date: 2026-09-03. Baseline: branch `regions-and-identity` at `66f5688`, suite 485
passed / 1 skipped. Builds directly on `regions_bridge_identity_and_cleanup_plan.md`
(the arrangement, `Region`, `Arc`, partition elements, `image_of_element`).

## Goal

Turn the arrangement of a trellis into a dual graph and read the symbolic dynamics
of the tangle off it, without growing anything:

1. **Bridge classes**: the pairs of stable partition elements that the trellis's
   bridges connect. One class per distinct element pair; a class is a symbol.
2. **Dual graph**: one node per face of the plane (regions, open faces, the merged
   outer faces, the unbounded face — exactly one per face), one node per stable arc
   (an "open circle" carrying the element on each side), and edges joining a face to
   every stable arc on its boundary.
3. **Fill**: the arc nodes on which an image bridge may cross the stable manifold
   at a crossing the trellis has not seen are marked passable ("filled circles").
4. **Element images**: every element mapped forward to the element(s) covering its
   image, from the iterate table where it is populated and from canonical-distance
   scaling otherwise.
5. **Words**: for every bridge class, the walk through the dual graph from the
   image of its first element to the image of its second, recorded as the sequence
   of bridge classes it traverses. `a -> ab`.
6. **Transition graph**: symbol table, substitution rules, and a networkx digraph
   plus a matrix, with a text report and a plot of the dual graph.

Out of scope: any GUI, blasting, growing manifolds to verify (used only in tests).

## Decisions taken with the user (2026-09-03)

| Question | Decision |
|---|---|
| Which side's elements define a bridge's class? | The side the bridge itself lies on at each end (its *row*). One class per bridge. |
| Several walks between the image elements? | Enumerate; require exactly one distinct **word**; raise with all candidates otherwise. |
| How is an element's (and a crossing's) image found? | Iterate table first; when an endpoint's iterate is unregistered, scale its cdist by `per_step_beta("stable")` onto the advanced branch. |
| Which nodes are filled? | Derived from last crossings: on stable branch `j`, arcs with cdist in `(cdist(f(T_{j-1})), cdist(T_j)]`, `T` the outermost crossing of each branch and `j-1` the branch one step back along `advance_key`. With every branch trimmed at the strong pip orbit this is `(f^k(q0), q0]` on the pip's branch and empty elsewhere. Logged when the two disagree. |

Assumptions stated to the user and not contested:

- The face point is a hub. The stored graph is bipartite (face nodes, arc nodes); a
  "through-face edge" (arc, face, arc) is derived, never stored, so edge count is
  linear in the boundary length.
- Anchor bridges (periodic point out to the first crossing) are bridge classes; their
  first element is the one whose `lo_id` is None.
- Partial bridges are ignored.
- Symbolic dynamics is under ONE map step. On a period-k orbit the classes cycle
  between branches; nothing is composed to `f^k`.

## Definitions

- **ElementRef**: the global identity of a partition element,
  `(branch_key, side, element_id)`. A frozen dataclass so it hashes and prints.
- **Row** of a bridge end: the side of the stable branch the bridge approaches that
  crossing from, in the stable dynamical orientation (`_row_at` today).
- **BridgeClass**: `(first: ElementRef, second: ElementRef)`, ordered like a
  `BridgeId` (unstable dynamical direction). Frozen, hashable.
- **Face node**: one merged face of the plane (see B.2). **Arc node**: one stable
  `Arc` of the arrangement (identity = `Arc.edge_key`).
- **Face side of an arc**: which side of the stable branch, in dynamical orientation,
  a face lies on. Combinatorial: the arrangement traverses every face with the face on
  the RIGHT of each half-edge (bounded faces clockwise). A stable arc traversed
  `hi -> lo` (`Arc.reverse`) runs in the dynamical direction, so that face is on the
  `right`; traversed `lo -> hi` the face is on the `left`.
- **Row of a bridge end, combinatorially**: at crossing `p` with sign `s = crossing_sign(p)`,
  the ray `u+` is on the `left` of the dynamical stable direction iff `s > 0`
  (`cross(s-, u+) = cross(u+, s+) = s`), and `u-` is on the opposite side. A bridge
  leaves its first endpoint along `u+` and enters its second along `u-`, so
  `row(first) = left iff sign(first) > 0`, `row(second) = left iff sign(second) < 0`.
  For a same-branch bridge this reproduces invariant I2 (same row both ends) exactly
  when the two signs differ. This rule needs no geometry and works at an anchor.
- **Word**: a tuple of `BridgeClass` in traversal order.

## Phase A — Element identity and bridge classes

**A.1 `ElementRef`** in `TopologyResults.py`. `PartitionInterval.ref` property returns
it (raises if unstamped). `StablePartitionResult.ref(element_id)`.

**A.2 Combinatorial row.** `StablePartition.row_of_end(trellis, bridge_id, endpoint) -> Side`
from `crossing_sign` as defined above. Test: agrees with the geometric `_row_at` on
every non-partial bridge of the k=10 and nested p3 fixtures where `_row_at` is not
None; on same-branch bridges it never produces a mismatch (I2).

**A.3 `BridgeClass`** in new `topology/BridgeClass.py`:
```python
@dataclass(frozen=True)
class BridgeClass:
    first: ElementRef
    second: ElementRef

def bridge_classes(trellis, partitions: Iterable[StablePartitionResult],
                   *, bridge_ids: Optional[Iterable[BridgeId]] = None
                   ) -> dict[BridgeClass, list[BridgeId]]
```
For each non-partial bridge id of the trellis (default: every bridge with an id, so
the all-fixed-points trellis covers nested and heteroclinic cases), take the row at
each end, resolve the element with `element_of_intersection` of the partition result
for `(branch of that endpoint, row)`. A missing partition for a branch is a
`ValueError` naming the branch. Classes keep their members sorted by `BridgeId`.

**A.4 Session gathering.** `TangleSession.bridge_classes(fixed_points=None)` collects
partitions from every cached per-fixed-point trellis (the same scan
`partition_element_for` does), warns about any unpartitioned stable branch, and calls
A.3 on the all-fixed-points trellis. Cached by `workbench.generation`.

Tests: every non-partial bridge lands in exactly one class; on k=10 the anchor bridge
is its own class; on p3 both tangles' bridges are classed and no class mixes fixed
points (there are no heteroclinic crossings in the fixture).

✅ DONE 2026-09-03 (commit 470cc60; suite 485 -> 529 with A.4). Landed: `ElementRef` with `label`/`fixed_point`/`orbit_index`/`branch_index`, `PartitionInterval.ref`, `StablePartitionResult.ref`; `row_of_end`/`rows_of_bridge` from `crossing_sign` (0 mismatches against `_row_at` on 14 k=10 and 60 p3 bridge ends); `BridgeClass`, `bridge_classes`, `element_sort_key`/`class_sort_key` (deterministic class order by fixed-point index, orbit, branch, side, element id); `TangleSession.bridge_classes` + `_gathered_partitions` + `_partition_signature`, cached on (generation, partition signature) so re-partitioning at one generation invalidates. Fixtures `k10_session`/`k10_partitioned`/`p3_partitioned` live in conftest. Deviations: on k=10 the anchor bridge (0,2) shares its class with bridge (1,6) (both ends same row, elements R#0/R#2; the partition has three elements per side), so "the anchor bridge is its own class" is false on the fixture and the test pins the definitional claim (`lo_id is None` on the class's first element) instead; `ElementRef.label` keys on the fixed point's period only (two same-period saddles would collide), so symbols are compared by object, never by label; the session warning fires for a branch unpartitioned on EITHER side.

## Phase B — The dual graph

New `topology/DualGraph.py`, class `DualGraph(arrangement, partitions)`.

**B.1 Arc nodes.** One `ArcNode` per stable `Arc.edge_key`:
`arc`, `left: ElementRef`, `right: ElementRef`, `filled: bool`, `faces: (left_face, right_face)`.
The element on a side is the one whose span OWNS the arc's midpoint cdist
(`owns_cdist`), so a singleton pinched at an endpoint is never the arc's element.

**B.2 Face nodes and the merge rule.** Start from `arrangement.faces` (all of them).
Union-find merges:
- Each component's **outer face** (the one face of that component whose traversal
  polygon has positive signed area) with the innermost containing face of another
  component that geometrically contains it (smallest `|area|` among
  `containing_faces` whose polygon contains one of the component's nodes).
- Outer faces with no containing face merge into one global **unbounded** node.
Every merged face is one `FaceNode` with `faces: list[Region]`, `kind` in
`{"region", "open", "outer"}`, and the arc nodes it touches. An arc node touches
at most two face nodes; `faces` records which side each is on (definition above).
Assert: every arc node has a face on both sides after merging.

**B.3 Fill.** For each stable branch `j` of the all-fixed-points trellis:
`T_j` = last id in `ordered_ids()`; `key_prev = fp.advance_key(key_j, -1)`;
`T_prev` = last crossing of that branch; `c_img` = `stable_cdist(iterate(T_prev, 1))`
if registered, else `stable_cdist(T_prev) * per_step_beta("stable")` (side flip is
irrelevant here). Fill every arc node on branch `j` whose `lo_cdist >= c_img - tol`.
An arc that straddles `c_img` is filled and logged. Record the derived segment per
branch. When the trellis has a strong pip, log at INFO whether the derived segment
equals `(iterate(q0, k_value), q0]`; a mismatch is a warning, not an error.

**B.4 Plot.** `topology/plotting.plot_dual_graph(dual_graph, ax=None)`: hollow markers at
arc midpoints, filled markers for passable nodes, a dot at each face node's point,
edges as thin lines. Face point: `representative_point` for a region; for a merged or
open face, the mean of its arc midpoints, pushed outside the bounding box for the
unbounded node. `Trellis`/`TangleSession` wrappers follow the existing delegation.

Tests (k=10, p3): one face node per merged face, exactly one `outer` node per level
of nesting (p3: the p3 outer face is merged with the p1 containing face, and the p1
outer face is the unbounded node); the combinatorial face side agrees with
`_side_of(stable frame, representative_point - arc midpoint)` on every closed region;
the fill on k=10 is exactly the arcs in `(f(q0), q0]` and on p3 exactly
`(f^3(q0), q0]` on the pip's branch (q0 from `strong_pip_intersection`); every
filled node has a face on both sides.

✅ DONE 2026-09-03 (commit c9df177; suite 552 -> 560). Landed: `topology/DualGraph.py` with `ArcNode` (element per side = owner of the arc's midpoint cdist), `FaceNode` (`kind` region/open/outer, `is_unbounded`), `DualGraph(arrangement, partitions, *, strong_pips=None)` with `arc_nodes`, `face_nodes`, `unbounded`, `fill_segments`, `arc_nodes_on/of`, `element`, `element_at`, `face_of`, `element_images` (C.3), `graph()`, `summary()`; `plotting.plot_dual_graph`/`face_point` (`clip_to_arcs=True` because a region's representative point can sit far outside the tangle). Verified by an independent reviewer: face side agrees with geometry on all 11 k=10 and 41 p3 closed-region arcs (two constructions); p3 has exactly one unbounded node and one merged outer node holding p1's containing face plus p3's outer face with consistent sides; fill is `(f(q0), q0]` on k=10 (4 arcs) and `(f^3(q0), q0]` on the p3 pip's branch (2 arcs). Deviations: each component's outer face is OPEN on both fixtures, so the outer face is found by the sign of the traversal area computed from the arc polylines (Region.area raises on open faces); on p3 the period-1 pip's own branch is filled too, `(f(q0), q0]` (the rule, not a special case — the brief's "empty elsewhere" meant the non-pip branches); `Trellis.iterate` now composes ±1 links (needed for `f^3(q0)`; committed with Phase C); filled iff `hi_cdist > c_img + tol` (straddlers filled and logged); a predecessor branch with no crossing past its anchor fills the whole branch (logged); `ArcNode.midpoint` is Optional; the session wrapper for the plot arrives with D.6.

## Phase C — Images of crossings and elements without the table

**C.1 `Trellis.image_cdist(intersection_id, n=1) -> (branch_key, cdist, from_table: bool)`.**
Table when `iterate(id, n)` exists; else `advance_key(key, n)` and
`cdist * per_step_beta("stable") ** n` (and the unstable analogue, parameterised by
stability). Sign of `n` handled.

**C.2 `image_of_element` fallback.** When an endpoint iterate is missing, use C.1's
scaled cdists instead of returning None; keep None only for an unbounded end. Return
type unchanged. The existing test that pins None for a missing iterate is retargeted
to an unbounded end.

**C.3 `element_images(partitions) -> dict[ElementRef, list[ElementRef]]`** on
`DualGraph` (all elements of all partitions it was built from), with the side carried
through `orientation_preserving`.

Tests: on the p3 fixture C.1 agrees with the table wherever the table is populated
(relative error under `cdist_tol`); `image_of_element` still passes its existing
tests; the fallback gives the same elements as the table on every element whose
iterates are registered.

✅ DONE 2026-09-03 (commit f9607bc; suite 529 -> 552 incl. B/A.4 rework). Landed: `Trellis.image_cdist(id, n, stability)`; `image_of_element(..., use_table=, partitions=)` falling back to scaled cdists (None only for an unbounded end); `Trellis.iterate` composes ±1 table links (the registry table holds only what inference recorded, n = ±1 on the fixtures). C.3 `element_images` was implemented by the Phase B agent on `DualGraph` (it needs the graph's partitions). Deviations: the scaling law matches the table to 1.35e-3 relative error on p3 (1.8e-4 on k=10; pure polyline discretisation, converging with `area_cutoff`), so the C.1 bound is `SCALING_RTOL = 1e-2` rather than `cdist_tol`; a SCALED endpoint is snapped onto a partition boundary within `SNAP_RTOL = 2e-3` because forced scaling otherwise pulled a neighbouring element into 27/56 p3 covers — the snap never fires on the default path on either fixture, so the true fallback is pinned by coverage only (no ground truth exists for an unregistered iterate).

## Phase D — Walks, words, and the transition graph

**D.1 Endpoints of a walk.** For a bridge with endpoints `p, q` and rows `r_p, r_q`:
`side_p'` = `r_p` flipped once per step when not orientation preserving; the start
element is the owner (`owns_cdist`) of `image_cdist(p)` on `(image branch, side_p')`;
likewise the end. Start nodes = arc nodes of the start element; end nodes = arc
nodes of the end element.

**D.2 Walk.** In the bipartite graph: leave a start node into its face on `side_p'`;
from a face, step to any arc node on its boundary; an intermediate arc node must be
`filled` and is crossed to its other face; finish on an end node entered from the
face on `side_q'`. Simple paths only (no arc node twice), DFS with a hard cap (raise
past 10 000 paths). Each path maps to a word: for every face traversal
`(arc_in, face, arc_out)` the letter is `BridgeClass(element of arc_in on the face's
side, element of arc_out on the face's side)`. A zero-face walk (start element is
the end element and the bridge does not leave) is the empty word and is rejected —
every image bridge crosses at least one face.

**D.3 Uniqueness.** Distinct words are collected; exactly one is accepted. Otherwise
`AmbiguousWalkError` listing every distinct word with one path each. A letter that is
not a known class is a **new symbol**: kept, labelled, and logged as a warning (the
trellis is not closed under the map at this extent).

**D.4 Per bridge, then per class.** The word is computed per bridge (D.1 uses the
bridge's own endpoints) and grouped by class. All bridges of one class must yield
the same word; a disagreement is a `ValueError` naming the class and both words — it
means the partition is not fine enough, which is exactly the diagnostic the user
wants surfaced.

**D.5 `SymbolicDynamics`** result in new `topology/SymbolicDynamics.py`:
`symbols: dict[str, BridgeClass]` (a, b, …, then a1, a2 … past 26), `rules:
dict[str, str]` (`"a": "ab"`), `words: dict[BridgeClass, Word]`,
`transition_graph: nx.DiGraph` (edge weight = multiplicity),
`transition_matrix() -> NDArray[int]`, `describe() -> str` printing the symbol table,
the rules and the matrix. Symbol order: by fixed point, then branch, then element ids
of the first endpoint, so two runs of the same trellis name symbols identically.

**D.6 Session entry points.** `TangleSession.dual_graph(fixed_points=None)` and
`TangleSession.symbolic_dynamics(fixed_points=None)`, both cached by generation on
the all-fixed-points trellis, mirroring `arrangement()`.

Tests:
- **Geometric validation (the important one).** Build k=10 and p3 sessions, compute
  the symbolic dynamics, THEN grow the unstable manifold one more step, recompute
  with `preserve_ids=True`, re-create bridges, and for every original bridge with a
  registered image (`image_bridges`) map each image bridge to its class (A.3) and
  compare with the walked word. They must be equal letter for letter.
- Every anchor class's word begins with itself (the image of `[fp, p0]` contains it).
- No class's word is empty; every letter of every word is a known class on the
  fixtures (assert, so a new symbol on a fixture is a test failure, since the fixtures
  are grown far enough).
- The transition graph is strongly connected on the k=10 fixture restricted to the
  classes inside the resonance zone (pin the count of nodes and edges).
- `AmbiguousWalkError` on a hand-built dual graph with two filled nodes on one face.

✅ DONE 2026-09-03 (commit c349571; suite 560 -> 597, incl. D.6). Landed: `DualGraph.walk_bridge` / `Walk` / `AmbiguousWalkError` / `Arrangement.face_at_stub`; `topology/SymbolicDynamics.py` (`symbolic_dynamics`, `word_of_bridge`, `SymbolicDynamics` with `symbols`, `rules`, `words`, `members`, `walks`, `new_symbols`, `transition_graph`, `transition_matrix()`, `case_counts()`, `display_label()`, `describe()`); `TangleSession.dual_graph` / `symbolic_dynamics` / `plot_dual_graph` / `plot_transition_graph` (cached on generation + partition signature); `plotting.plot_transition_graph`. Independently verified: every bridge of both fixtures (k=10: 7, cases i/ii/iii 3/1/3; p3: 30, 26/2/2) yields letter for letter the classes of its image bridges after one more growth step; the k=10 words were hand-read (anchor `d -> dce`, lobes `c -> a`, `e -> dbe`, all forced by the geometry) and refine correctly at ten steps (`f -> feg`); p3 words cycle orbit 0 -> 1 -> 2 -> 0 and the anchor classes begin with themselves. Deviations: (1) the plan's walk presumes an unregistered image, but most bridges in a grown trellis already have theirs registered and it runs ALONG existing arcs, so the word is a registered prefix (image bridges tiling, sliced by index on the image unstable branch from the iterate table) plus a walked suffix that starts from the u+ stub of the last registered crossing (case ii) or from the arc owning the scaled image (case iii); (2) a walk's endpoint arcs are never filled — they are images of crossings the trellis has seen and lie at or below `c_img` — only intermediate arcs must be; the brief's opposite requirement was wrong and was dropped; (3) p3 at 13 unstable steps has four unrealised innermost self-classes (`L#0 -> L#0`, `R#1 -> R#1` on fp3 orbit 1 and on fp1), so "no new symbols on the fixtures" is pinned on k=10 only and the four are pinned on p3; (4) the geometric validation classes the grown image bridges against the ORIGINAL partitions (growth re-partitions, so element ids do not survive) by `row_of_end` + `owns_cdist`; (5) D.4 never fired — the two-member k=10 classes agree across their case-(i) and case-(iii) members; (6) a first version used a 1%-of-cdist window to find the registered span and broke one growth step out (crossings at u ~ 6000 with gaps of 20–40); fixed by slicing by index, with a regression test that grows k=10 one step (cases 7/1/7, 7 symbols, 0 new, coarsened words equal the originals).

## Phase E — Script, docs, memory

- `scripts/henon_symbolic_dynamics.py`: k=10 and nested p3, prints `describe()`, draws
  the tangle with the dual graph overlaid and the transition graph beside it.
- CLAUDE.md table rows for `BridgeClass`, `DualGraph`, `SymbolicDynamics`, the
  `ElementRef` name, and the fill rule; the `Note` on same-side bridge rows gains the
  combinatorial rule.
- Re-run the plan-doc checklist: each phase gets a ✅ DONE line with deviations.

✅ DONE 2026-09-03 (commit 7ff12b4; suite 597 passed / 1 skipped, baseline 485). Landed: `scripts/henon_symbolic_dynamics.py` (both fixtures, `describe()` to stdout, tangle + dual graph + a zoom on the period-3 tangle beside the transition graph, `figures/henon_symbolic_dynamics_{k10,p3}.png`); CLAUDE.md rows and the "Bridge classes, the dual graph, and words" section; `topology/__init__` export symmetry (`rows_of_bridge`, `element_sort_key`/`class_sort_key`, `WalkCase`, `symbol_label`, `plot_dual_graph`/`plot_transition_graph`/`face_point`); new-symbol warnings use `display_label`. The style pass found the branch already clean (ruff -F, Google docstrings, future annotations, logging). Deviations: the CLAUDE.md sentences about endpoint arcs and the row rule live in the new section as expectations, not in the "Invariants the code asserts" paragraph (nothing asserts them). Two Sonnet agents were cut off by a rate limit mid-Phase E and restarted; nothing was lost.

## Files

| File | Change |
|---|---|
| `topology/TopologyResults.py` | `ElementRef`; `PartitionInterval.ref`; `StablePartitionResult.ref` |
| `topology/StablePartition.py` | `row_of_end` (combinatorial) beside `_row_at` |
| `topology/BridgeClass.py` | new: `BridgeClass`, `bridge_classes` |
| `topology/DualGraph.py` | new: `ArcNode`, `FaceNode`, `DualGraph`, walks, `AmbiguousWalkError` |
| `topology/SymbolicDynamics.py` | new: `SymbolicDynamics`, `symbolic_dynamics(dual_graph, classes)` |
| `topology/Trellis.py` | `image_cdist`; `image_of_element` fallback |
| `topology/plotting.py` | `plot_dual_graph`, `plot_transition_graph` |
| `topology/__init__.py` | exports |
| `loom/TangleSession.py` | `bridge_classes`, `dual_graph`, `symbolic_dynamics` + caches |
| `tests/test_bridge_class.py`, `tests/test_dual_graph.py`, `tests/test_symbolic_dynamics.py` | as listed per phase |
| `scripts/henon_symbolic_dynamics.py` | new |

## Order

A -> B -> C -> D -> E. A and C are independent of each other and can run in parallel;
B needs A (element refs on arc nodes); D needs all three.

## Risks and where they are caught

- **Handedness.** The partition's `left` (positive cross in the stable dynamical
  direction) and the arrangement's counter-clockwise rotation must be the same
  orientation of the plane. B's test against `_side_of` on closed regions catches a
  flip; A.2's test against `_row_at` catches it independently.
- **Non-invariant partition.** The Deferred note says the hole set is not forward
  invariant at the computed extent. D.4 turns that into a named error per class
  rather than a silent wrong word.
- **Straddling fill boundary.** When `f(T_{j-1})` is not a registered crossing the
  boundary falls inside an arc; that arc is filled (conservative) and logged.
- **Inversion.** Not validated (four coincident anchor nodes, see the previous plan).
  `advance_key` and `per_step_beta` are used throughout so nothing is period-1-specific,
  but no inversion test is added.
- **Path explosion.** The DFS is capped; the fixtures have tens of faces, not
  thousands.

## Deferred (found during implementation, not fixed)

- Registry ids are not reproducible between builds in one process (three k=10 builds gave three id assignments for the same eight crossings; `Tangle` keys segment edges on `id()`), so no test may hard-code a crossing or bridge id; class labels, element ids, words and `describe()` are id-free and stable. Pre-existing (also in the previous plan's Deferred list).
- The registry's `IterateTable` stores only what inference recorded (n = ±1 on the fixtures); `Trellis.iterate` now composes single-step links, but `iterate_orbit` stops at a repeated stable branch (so it cannot express `f^3(q0)` on p3) and nothing multi-step is ever tabled. Deeper iterate inference would give the C.2 fallback a ground truth; today the true fallback (unregistered iterate) is pinned by coverage only.
- `Trellis.scale_cdist` still uses the unstable beta for a stable distance; correct only because `per_step_beta("stable")` is the reciprocal on an area-preserving map. `image_cdist` deliberately does not route through it.
- `_snap_to_partition_boundary` re-sorts the image branch's boundary set per call; memoise per `StablePartitionResult` and bisect if the dual graph makes element images hot.
- `ElementRef.label` / `BridgeClass.label` key on the fixed point's period only; two distinct saddles of the same period collide. `SymbolicDynamics.display_label` uses the fixed-point index instead; never use labels as identities.
- `bridge_classes` accepts a reversed `BridgeId` silently (rows swapped) and raises a bare `KeyError` for an unknown id; `_element_at_end`'s guard is an `assert` (degrades under `python -O`).
- Nothing in the suite exercises a bridge end where the geometric `_row_at` is unreadable (0 such ends on both fixtures); the claim that `row_of_end` answers there is unverified by a fixture.
- `tests/test_session_trellis_cache.py` defines its own, differently built `k10_session` that shadows the conftest one.
- The p3 fixture builds in ~0.15 s, not the ~10 s the briefs assumed; the `@pytest.mark.slow` marks are cosmetic.
- A merged outer node stands for a non-simply-connected piece of plane (p3: p1's containing face plus p3's outer face); a walk may pass between any two of its walls, which is right for one tangle but over-generous once two tangles are linked heteroclinically. No fixture has a heteroclinic crossing.
- A predecessor stable branch with no crossing past its anchor makes the whole branch filled (logged at INFO); conservative, never exercised.
- `plot_dual_graph`: a region's `representative_point` can sit far outside the tangle (k=10 region 3 at x = 72), hence `clip_to_arcs=True`; the mean-of-midpoints hub of a merged outer node is not guaranteed to lie in the annulus. `plot_dual_graph`/`face_point`/`plot_transition_graph` exports from `topology/__init__` are added by the Phase E style pass.
- The dual-graph and symbolic-dynamics session caches key on (generation, partition signature) and NOT on the strong pips (which only drive the fill log line); `set_strong_pip` without a re-partition returns the cached graph unless `rebuild=True`. `trellis(fp)` and `trellis([fp])` are distinct cache keys, so one pip can be gathered twice (duplicates a log line).
- Walk letters use arc-midpoint element ownership while classes use crossing ownership; they agree because no partition boundary falls strictly inside an arc and because an image landing ON a registered crossing without a table entry is rejected (`_arc_owning`) rather than assigned to one of the two arcs.
- Case (i) validations (a bridge whose image is fully registered) are near-tautological — both sides read the same iterate table; the independent evidence is the case (ii)/(iii) bridges (4 on k=10, 4 on p3) plus D.4 agreement inside the two-member k=10 classes. The plan's "no new symbols on the fixtures" is an extent property, not a fixture property: p3 at 13 unstable steps has 4 unrealised innermost self-classes (pinned), k=10 at 9 steps has none.
- `test_k10_zone_classes_are_strongly_connected` pins a 2-node subgraph (weak, but what the plan asked for).
- Inversion is still unexercised by any of this: `advance_key`/`per_step_beta` are used throughout, the `orientation_preserving` side flip in `image_of_element`/`element_images`/the walk is read-checked only.
- `plot_transition_graph` uses `nx.spring_layout(seed=0)` with no `k`/`iterations` knobs; on p3 the labels `t`/`z` and `b`/`c1` overlap whatever the figure size.
- The `SymbolicDynamics` new-symbol WARNING and `describe()` now both use `display_label`; `BridgeClass.label` (period-keyed) remains for `__str__`.
