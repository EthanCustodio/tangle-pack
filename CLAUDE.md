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
| topology | `Trellis`, `TrellisBranch` | same names | Per-generation snapshot of registry + bridges + per-branch orderings; strong pips, pseudoneighbors, holes, partitions, `element_for`, `image_of_element`, lazy `arrangement` |
| topology | `StrongPip`, `Pseudoneighbor`, `StablePartition` | same names | The three algorithms (see the PDFs in docs/); `StablePartition` also holds the invariants `check_holes_share_bridge_side` / `check_bridge_rows_consistent` |
| topology | `Hole`, `PseudoneighborPair`, `PartitionInterval`, `StablePartitionResult`, `Arc`, `Region` | `TopologyResults.py` | Result dataclasses; partition elements carry `element_id`; `Region` = closed minimal face with lazy `area`, `representative_point`, `contains` |
| topology | `Arrangement` | `Arrangement.py` | Half-edge planar arrangement of all arcs: rotation at each node from `crossing_sign` alone, virtual tails for dangling ends, faces by linear traversal; `regions`, `containing_faces`, `open_faces`, `image_of`/`preimage_of` (area-verified) |
| topology | `plotting` | `plotting.py` | All topology drawing routines (Trellis `plot_*` delegate here) |
| loom | `TangleSession` | `TangleSession.py` | Facade: trellis/arrangement caches keyed by `workbench.generation`, resonance zones, fan-outs over fixed points |
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
