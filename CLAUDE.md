# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Project Is

`tanglepack` is a Python library for computing and visualizing **heteroclinic/homoclinic tangles** — the stable and unstable manifolds of saddle fixed points in 2D discrete dynamical systems (area-preserving maps). The core scientific workflow is: define a map → find a fixed point → initialize manifolds → grow them iteratively → detect intersections → extract bridges and resonance zones.

There is currently no GUI frontend: the old Dash dashboard (`tanglepack_webdash`) and PySide6 GUI (`tanglepack_gui`) were removed in July 2026 (recoverable from git history) and a new GUI will be built from scratch later. The library is used through scripts, notebooks, and the `TangleWorkbench`/`TangleSession` APIs.

## Fundamental Invariant — Which Manifolds Can Intersect

This is a hard physical law of the system, not a convention, and it must inform every piece of intersection/bridge logic:

- **Two unstable manifolds can never intersect each other.** Two stable manifolds can never intersect each other. This holds whether the manifolds belong to the same fixed point or to different fixed points.
- **Only an unstable manifold can cross a stable manifold** (homoclinic when they share a fixed point, heteroclinic when they belong to different fixed points).

Why: a point on a stable manifold converges to that manifold's fixed point under forward iteration; a point on an unstable manifold converges to its fixed point under backward iteration. A shared point of two stable manifolds would have to converge forward to *both* fixed points at once (and a shared point of two unstable manifolds would have to converge backward to both) — impossible. The same uniqueness argument forbids a manifold from crossing itself.

**Implication for the code:** a same-stability pair (u×u or s×s) appearing in `_intersecting_segments` is *always* a numerical artifact — two near-parallel polygonal approximations straddling near a tangency, not a real crossing. Such pairs should be filtered out (and are worth logging), but they must never be treated as a legitimate geometric case to model.

## Area Preservation

These maps are **area-preserving**: mapping any region forward (or backward) produces a region of exactly the same area. This is a fundamental property of the system, not an approximation.

A consequence used by the topological layer: for an intersection point, the product of its stable and unstable canonical distances behaves like a preserved area (one forward iterate scales the unstable cdist up by the per-step eigenvalue factor and the stable cdist down by the same factor, so the product is unchanged along an iterate chain). Note, however, that two *different* iterate chains generally have *different* such products — equal products do not imply two points are on the same chain, so do not use the product alone to decide chain membership. Compare the stable and unstable canonical distances individually (a collision on both means the same point) when you need to know whether two intersections are iterates.

## Bridge classes and the dual graph

The topology layer reads the combinatorics of a trellis's arrangement without growing anything further, in three steps:

- **Row (which side a bridge end is on).** The stable manifold is oriented anchorward. At a crossing with sign `s = crossing_sign(p)`, `row(first) = left iff s > 0` and `row(second) = left iff s < 0` — a bridge leaves its first endpoint along `u+` and enters its second along `u-`, and `u+` sits on the left of the dynamical stable direction exactly when `s > 0`. This combinatorial rule needs no geometry, agrees with the same-branch invariant above whenever the two signs differ, and works even at a synthetic anchor.
- **Bridge class.** A `BridgeClass` is the UNORDERED pair `{X, Y}` of `ElementRef`s the bridge's two rows resolve to — the homotopy class of arcs between those two elements — oriented anchor outward (`source` is the anchor-nearer element). A member bridge traversed in the unstable direction runs `source -> target` (direction `+1`, symbol `a`) or back (`-1`, symbol `a^-1`); both belong to the one class. A bridge with both ends in ONE element is a loop: it is the image of a bridge that did connect two elements and it will keep mapping into that element, so it is FOLDED into the class of its first non-loop ancestor along the registry preimage chain (`Trellis.iterate(-1)` + `bridge_between`, never the map), and a class holding a loop is INERT (it never appears in an itinerary). A loop with no registered ancestor stands alone as a `BridgeClass(x, x)` entry and is logged. Inertness is evidence-based: the k=10 exterior class is active at 9 unstable steps and inert at 10. `bridge_classes(trellis, partitions)` returns a `BridgeClassTable` of `BridgeClassEntry`s (members with directions, `letter`, `zone_key`); the session letters ACTIVE classes from its persistent `BridgeAlphabet` (a class seen before keeps its letter, a new one takes the next) and records each class's resonance zone via `classify_bridge`. Refined symbols that subdivide a class are future work and will hang off the entry, not change the class.
- **Face side of an arc.** The arrangement traverses every bounded face CLOCKWISE (the face on the right of each half-edge). A stable `Arc` traversed `hi -> lo` (`reverse=True`) runs in the stable dynamical direction, so that traversal's face is on the arc's `right`; traversed `lo -> hi` the face is on the `left`.
- **Minimal trellis.** The dual graph is built over a REDUCED trellis (`MinimalTrellis.py`): every bridge with a hole attached (`Hole.bounding_ids`), plus, for each hole bridge whose bridge class is ACTIVE, the bridges tiling its first forward image (`image_chain`, iterate table only; a subdivided image contributes every consecutive pair; a pair with no `Bridge` object is skipped with a WARNING). Inert hole bridges are kept but not mapped. The stable manifold is kept whole; the surviving nodes are the kept bridges' endpoints, the anchors and each stable branch's outermost crossing. The reduced snapshot is a synthetic `Trellis` sharing the parent's registry, and its `Arrangement` is built in `sparse=True` mode: unstable arcs come from the bridge list (not consecutive crossings) and an absent unstable ray is an absent wall (the two sectors merge), not a slit.
- **Partition families.** `PartitionFamily.py` holds one partition per `(branch, side)` behind a base class: `HomotopyPartition` wraps the hole-cut partition (`partition_stable_manifold`, unchanged), `IteratedHomotopyPartition` refines it by the minimal trellis's IMAGE bridges. The cut rule is in bridges: each image bridge cuts the side of the stable manifold it lies on (its row), the element under it closed at both crossings `[c1, c2]` and the neighbours open; an existing hole boundary is never re-cut; the other side is untouched at an image endpoint (a crossing shared by two consecutive image pieces is cut on both sides, once per piece). Iterated elements record `parent_element_id` (the homotopy element they refine) and `cut_by` (the innermost image bridge they lie under). `ElementRef` carries no family tag: never mix refs of the two families in one structure.
- **Fundamental segment.** The strong pip `q0` DECLARES the stretch where a walk may cross the stable manifold: on the pip's own stable branch, the edges meeting `(cdist(f^k(q0)), cdist(q0)]` with `k = k_value` (via `Trellis.image_cdist`), and nothing on any other branch. There the two side nodes of a stable edge are UNIFIED into one solid, traversable `StableNode`; everywhere else each edge has two open side nodes that are walls. Choosing a different candidate pip moves the unified set without trimming; no pip unifies nothing and warns. `TangleSession.dual_graph` keys its cache on the gathered pips for that reason. Faces carry NO bridge class (the relation is not one to one); a walk alternates face -> stable node -> face and collects the stable elements it visits.
- **The row invariant.** `row_of_end` reproduces the same-row invariant exactly when the two crossing signs differ — with equal signs the rows come out opposite, which is what a bridge that crosses the stable branch would look like, so the invariant check `check_bridge_rows_consistent` is the geometric guard.

## Commands

**Install (editable):**
```bash
pip install -e .
```

**Run tests:**
```bash
pytest
pytest tests/test_manifold_machine.py          # single test file
pytest tests/test_manifold_machine.py::test_fn # single test
pytest --cov=tanglepack                        # with coverage
```

**Run a scripted example:**
```bash
python scripts/tangle_workbench_test.py
```

## Core Library Architecture (`src/tanglepack/`)

Three subpackages plus examples. **Dependency rule: `numerics` and `topology` never import `loom`** (loom sits on top and weaves the other two together).

- `numerics/` — the numerical engine: maps, fixed points, manifold growth, crossing detection, bridges, the intersection registry.
- `topology/` — combinatorial structures over the registry: the `Trellis` snapshot, strong pips, pseudoneighbors, holes, the stable partition, and the planar `Arrangement` of regions.
- `loom/` — cross-layer algorithms that read topology and act on numerics: `TangleSession`, resonance zones, blasting.
- `examples/henon.py` — the one Hénon map definition (`henon_map(k, b)`, `henon_map_inverse`, `henon_jacobian`, `saddle_guesses`) used by tests and scripts (the notebooks still carry their own copies; see the plan's Deferred list).

The entry point for programmatic use is `TangleSession` (a facade over one `TangleWorkbench`, caching one `Trellis` per fixed-point selection and delegating unknown attributes to the workbench):

```python
import numpy as np, tanglepack
from tanglepack.examples.henon import henon_map, henon_map_inverse, henon_jacobian, saddle_guesses

session = tanglepack.TangleSession(henon_map(10), henon_map_inverse(10), henon_jacobian(10))
fp = session.construct_fixed_point(saddle_guesses(10, 1)["saddle"])
session.orient_eigenvectors(fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])})
session.initialize_both_manifolds(fp)
session.grow_n_times(fp, "unstable", num_iterations=8)
session.grow_until_turnaround(fp, "stable")
session.compute_intersections([fp])          # registers anchors + crossings, infers iterates
session.trim_stable_manifolds(fp)
session.create_bridges(fp)                   # bridges keyed by BridgeId
T = session.trellis(fp)                      # cached by generation
T.classify_strong_pips(); T.compute_pseudoneighbors(); T.punch_holes(); T.partition_stable_manifold()
A = session.arrangement()                    # faces -> Region objects
M = session.minimal_trellis()                # hole bridges + images of the active ones
P = session.iterated_partition()             # homotopy partition cut by M's image bridges
G = session.dual_graph()                     # StableNodes / FaceNodes over M's sparse arrangement
```

**Data model — how a manifold is stored.** Points form two simultaneous doubly-linked lists inside `Point`: the *geometric* list (`forward`/`backward`, order along the curve) and the *iterate* list (`next_iterate`/`prev_iterate`, which point maps to which). `BranchPoint` is the root anchored at the periodic point. `BaseManifold` holds `root`/`tail` (setters bump `version`) and one `_collect` walker behind every getter; bridges memoise their point arrays (read-only).

**Manifold key.** `ManifoldKey = (FixedPoint, Stability, orbit_index, branch_index)`. Every `BaseManifold` and `Bridge` carries its key (required at construction); `workbench.register_manifold(key, manifold)` asserts they agree. `FixedPoint.advance_key(key, n)` is the ONE rule for what n map steps do to a key (orbit index advances; the branch flips only when the orbit index wraps and the point has inversion). `FixedPoint.num_branches` is 2 iff inversion (`k_value == 2 * period`); `per_step_beta(stability)` is the per-map-step canonical-distance factor `|λ|^(1/period)` (a branch return is `k_value` steps and costs `λ^num_branches`); `branch_cycle(stability)` lists the `k_value` keys in forward-map order.

**Key classes and modules**

| Layer | Name | File | Role |
|---|---|---|---|
| numerics | `DynamicalSystem` | `DynamicalSystem.py` | Wraps `map`, `map_inv`, optional `jacobian` (axis-0 batches) |
| numerics | `FixedPoint`, `FixedPointSolver` | `FixedPoint.py`, `FixedPointSolver.py` | Periodic orbit, eigendata, `k_value`, `advance_key`, `per_step_beta`, `branch_cycle`; Newton/fsolve locator with saddle validation and an `orient` hook |
| numerics | `ManifoldInitializer`, `ManifoldMachine` | same names | Fundamental segments ("Kevin way", chain-walked with `advance_key`); growth by one map step with refinement |
| numerics | `BaseManifold`, `Bridge` | `BaseManifold.py`, `Bridge.py` | Linked-list manifold; `Bridge` = unstable arc between two consecutive crossings, `Bridge.id: BridgeId = (first_id, second_id)` ordered by unstable cdist, `None` when `partial` |
| numerics | `Intersection`, `IntersectionRegistry`, `IterateTable` | same names | THE single store of crossings (ids, cdists, keys, `crossing_sign`, `unstable_segment`); `generation` counter; sorted cdist arrays; per-branch `graph(bridges=)`; `iterate_table[id, n]` |
| numerics | `Tangle` | `Tangle.py` | R-tree index of segments (index state only); `resolve_crossings()`; `create_bridges(crossings)` cuts at the exact registry crossings |
| numerics | `TangleWorkbench` | `TangleWorkbench.py` | Orchestrator: fixed points, manifolds, growth (`grow_n_times`, `grow_until_*`, `grow_until(predicate)`), `compute_intersections` (registers one synthetic **anchor** per (unstable branch, stable branch) pair at cdist (0,0)), bridges by id (`bridge`, `bridges_at`), `generation` |
| numerics | `BridgeIterator`, `IterateInference`, `graphviz` | same names | Bridge forward images and derived genealogy (`image_bridges`/`preimage_bridges`); iterate inference (both cdists compared individually, never the product); networkx/matplotlib views |
| numerics | `geometry` | `geometry.py` | `arc_polyline`, `oriented_bridge_polyline`, `signed_polygon_area`, `point_in_polygon`, `polyline_midpoint`, `polygon_interior_point` |
| topology | `Trellis`, `TrellisBranch` | same names | Per-generation snapshot of registry + bridges + per-branch orderings; strong pips, pseudoneighbors, holes, partitions, `element_for`, `image_of_element`, lazy `arrangement`; `image_cdist(id, n, stability)` (table first, else `advance_key`/`per_step_beta`) and `iterate(id, n)` (composes `n = ±1` table links when the direct entry is missing) |
| topology | `StrongPip`, `Pseudoneighbor`, `StablePartition` | same names | The three algorithms (see the PDFs in docs/); `StablePartition` also holds the invariants `check_holes_share_bridge_side` / `check_bridge_rows_consistent`, plus `row_of_end`/`rows_of_bridge` (the combinatorial row, from `crossing_sign` alone) |
| topology | `Hole`, `PseudoneighborPair`, `PartitionInterval`, `StablePartitionResult`, `Arc`, `Region` | `TopologyResults.py` | Result dataclasses; partition elements carry `element_id`; `Region` = closed minimal face with lazy `area`, `representative_point`, `contains` |
| topology | `ElementRef`, `BridgeClass` | `TopologyResults.py`, `BridgeClass.py` | `ElementRef` = the global identity of a partition element, `(branch_key, side, element_id)`; a `BridgeClass` is the unordered, anchor-outward-oriented pair of `ElementRef`s a bridge connects, resolved by the ROW at each end (see "Bridge classes..." above); `bridge_classes(trellis, partitions)` returns a `BridgeClassTable` (`BridgeClassEntry` per class: `BridgeMember`s with direction ±1/0, loops folded into their preimage's class, `inert`, `letter`, `zone_key`, `symbol(bridge_id)`, `describe()`) |
| topology | `Arrangement` | `Arrangement.py` | Half-edge planar arrangement of all arcs: rotation at each node from `crossing_sign` alone, virtual tails for dangling ends, faces by linear traversal; `sparse=True` builds unstable arcs from the bridge list and stubs no absent unstable ray; `regions`, `containing_faces`, `open_faces`, `image_of`/`preimage_of` (area-verified) |
| topology | `MinimalTrellis`, `minimal_trellis`, `image_chain` | `MinimalTrellis.py` | The reduced trellis (hole bridges + images of the active ones): kept/dropped bridge ids, image chains, skipped pairs, the synthetic `sparse` snapshot and its lazy sparse `arrangement` (see "Minimal trellis" above) |
| topology | `PartitionFamily`, `HomotopyPartition`, `IteratedHomotopyPartition`, `Cut` | `PartitionFamily.py` | Partition families over `StablePartitionResult`s: indexed access, `element_at`/`interval_at` (midpoint ownership), `signature()`, `describe()`; the iterated cut with its `Cut` audit records (see "Partition families" above) |
| topology | `DualGraph`, `StableNode`, `FaceNode` | `DualGraph.py` | Dual of a minimal trellis's sparse arrangement: one `FaceNode` per merged face (union-find glues each component's outer face onto a containing face or the global unbounded node), two open `StableNode`s per stable edge (one per side) or one solid unified node on the fundamental segment, each carrying its element, interval, homotopy parent, cutting bridge and image elements; `image_face` (lazy, map-verified), `graph()`, `summary()` |
| topology | `plotting` | `plotting.py` | All topology drawing routines (Trellis `plot_*` delegate here), including `plot_minimal_trellis` (hole / image / dropped bridges), `plot_dual_graph` (open side nodes off each edge, solid unified nodes on it, face points, face-to-node edges), `stable_node_point`, and the two legend-handle helpers; `plot_stable_partition(labels=)` for homotopy and iterated rows on one axes |
| loom | `TangleSession` | `TangleSession.py` | Facade: trellis/arrangement caches keyed by `workbench.generation`, resonance zones, fan-outs over fixed points; `bridge_classes`, `minimal_trellis` and `iterated_partition` (cached on `(generation, partition signature)`), `homotopy_partition` (uncached wrapper of the gathered partitions), and `dual_graph` (also keyed on the gathered strong pips), all gathering partitions, holes and pips from every cached per-fixed-point trellis |
| loom | `BridgeAlphabet` | `BridgeAlphabet.py` | Persistent class letters `a..z, aa, ...`: `letter_for(bridge_class)` reuses or assigns the next unused letter; `class_of(letter)`; `reset()` |
| loom | `ResonanceZone`, `blast_zone` | `ResonanceZone.py`, `Blast.py` | Zone = stable manifold trimmed at a strong pip and its iterates; boundary is a list of `Arc`s; blasting iterates interior bridges |

**Invariants the code asserts (keep them):** cdist non-decreasing along a manifold; every registered crossing is unstable × stable; every bridge is an arc between consecutive crossings on one unstable branch; holes of one origin share `bridge_side` (parity-flipped under an orientation-reversing map); a bridge whose crossings share a stable branch has the same row at both ends; exactly one partition element owns each crossing per side; `image_of(region)` preserves area.

**Generation counter.** `IntersectionRegistry.generation` bumps on every registry mutation; `TangleWorkbench.generation` is a monotone token over its own mutations, the registry generation and every manifold's `version`. A `Trellis` records `_built_generation`; `TangleSession.trellis()` / `arrangement()` rebuild when it differs. There is no manual invalidation (`invalidate_trellises()` is a deprecated no-op).

**Intersection detection.** `compute_intersections` bulk-loads all manifold segments into the rtree (`Tangle.add_manifolds`); `iterate_bridge` registers image bridges query-only (a bridge is unstable, so everything it can cross is stable and already indexed). Crossings are resolved once, registered in the registry, then bridges are cut at those exact crossings. Same-stability pairs are discarded as numerical artifacts; near-parallel pairs are logged and discarded. Registry ids survive `compute_intersections(preserve_ids=True)`.

## Coding Style

- **Docstrings**: Google-style with `Args:`, `Returns:`, `Raises:`, and `Note:` sections. Write them for every public class, `__init__`, and method. One-sentence private helper docstrings are fine.
- **Dev Notes**: Module-level `"""Dev Notes: ..."""` blocks are the preferred place for longer-term design questions and open issues — not inline TODOs scattered through methods.
- **Type hints**: Always annotate function signatures. Use `Literal["unstable", "stable"]` for the stability parameter. Use `Optional[X]` and `NDArray[np.float64]` from `numpy.typing` for array return types.
- **imports**: `from __future__ import annotations` at the top of every file that uses forward references. Group: stdlib, then numpy/scipy, then local `.` imports.
- **Naming**: `PascalCase` for classes, `snake_case` for everything else. `_single_leading_underscore` for private methods.
- **Avoid bare print statements**: Use `logging` (the module already configures a `NullHandler` logger in `ManifoldMachine.py`). Debug output that ends up in production code is a recurring issue — prefer `logger.debug(...)`.
- **Comments**: Sparse inline comments only when the invariant is non-obvious. Do not leave blocks of commented-out old code in committed files — the dev notes pattern or git history is the right place for those.
- **Assertions**: Used to check invariants (e.g., cdist ordering after merge). Keep them; they are the primary correctness guard here.
