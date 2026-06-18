# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Project Is

`tanglepack` is a Python library for computing and visualizing **heteroclinic/homoclinic tangles** — the stable and unstable manifolds of saddle fixed points in 2D discrete dynamical systems (area-preserving maps). The core scientific workflow is: define a map → find a fixed point → initialize manifolds → grow them iteratively → detect intersections → extract bridges and resonance zones.

The frontend built on top of the core library is **`tanglepack_webdash`**: a browser-based app using Plotly Dash.

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

**Run the web dashboard:**
```bash
tanglepack-dash
# or: python -m tanglepack_webdash.app
```

**Run a scripted example:**
```bash
python scripts/tangle_workbench_test.py
```

## Core Library Architecture (`src/tanglepack/`)

The entry point for all programmatic use is `TangleWorkbench`. It wires together the other components and exposes a high-level API:

```python
import tanglepack, numpy as np

wb = tanglepack.TangleWorkbench(my_map, my_map_inverse)
fp = wb.construct_fixed_point([4, -4])
wb.orient_eigenvectors(fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])})
wb.initialize_both_manifolds(fp)
wb.grow_n_times(fp, "unstable", num_iterations=6)
wb.grow_until_turnaround(fp, "stable")
intersections = wb.compute_intersections(fp)
bridges = wb.create_bridges(fp)
graph = wb.build_intersection_graph(fp)
```

**Data model — how a manifold is stored:**

Points form two simultaneous doubly-linked lists inside `Point`:
1. **Geometric list** (`forward`/`backward`): ordering along the manifold curve.
2. **Iterate list** (`next_iterate`/`prev_iterate`): which point maps to which under the dynamical map.

`BranchPoint` is the root of a manifold (anchored to the fixed point). `BaseManifold` holds a `root` and `tail` pointer into the geometric list and exposes walk/plot helpers. `Bridge(BaseManifold)` is a truncated segment of the unstable manifold between two intersection points.

**Key classes and their roles:**

| Class | File | Role |
|---|---|---|
| `DynamicalSystem` | `DynamicalSystem.py` | Wraps `map`, `map_inv`, optional `jacobian` |
| `FixedPoint` | `FixedPoint.py` | Stores period, eigenstuff, `BranchPoint` list, k-value (doubles for inversion) |
| `FixedPointSolver` | `FixedPointSolver.py` | Newton's method to locate fixed points |
| `ManifoldInitializer` | `ManifoldInitializer.py` | Builds the initial short manifold segment from eigendata ("Kevin way") |
| `ManifoldMachine` | `ManifoldMachine.py` | Grows manifolds by one iteration of the map; inserts new `Point`s geometrically |
| `BaseManifold` | `BaseManifold.py` | Linked-list manifold with walk/plot methods |
| `Tangle` | `Tangle.py` | R-tree spatial index over manifold segments; computes intersections via `_Segment` pairs |
| `Bridge` | `Bridge.py` | Subclass of `BaseManifold` representing one piece of unstable manifold between intersection points |
| `TangleWorkbench` | `TangleWorkbench.py` | Orchestrates all of the above; manifolds keyed by `(FixedPoint, stability, orbit_index, branch_index)` |

**Manifold key convention:**

`TangleWorkbench.manifolds` is a dict keyed by `(FixedPoint, Stability, orbit_index, branch_index)`. `orbit_index` runs over the periodic orbit; `branch_index` is 0 or 1 for fixed points with inversion (negative eigenvalue, `k_value = 2 * period`).

**Intersection detection (`Tangle`):**

The `Tangle` class uses an `rtree` spatial index. Each adjacent pair of `Point`s on a manifold is registered as a `_Segment`. Candidate intersecting segments are found via bounding-box queries, then exact intersection is computed. Results stored in `_intersecting_segments` (set of `(seg_id1, seg_id2)` pairs), `_intersecting_coords`, and `_intersecting_points`. Always call `Tangle.clear_all()` before recomputing to avoid stale references. Per the fundamental invariant above, every legitimate crossing is exactly one unstable + one stable segment; any same-stability pair that slips into `_intersecting_segments` is a numerical artifact and must be filtered before building `Intersection`s.

## Web Dashboard Architecture (`src/tanglepack_webdash/`)

Built with Plotly Dash. Per-session state is stored in `sessions.py` (`WBState` dataclass in an in-memory `_REGISTRY` dict). Each browser session gets its own `TangleWorkbench` instance.

- `layout/` — Dash component trees (sections for fixed point, manifolds, bridges, orientation, etc.)
- `callbacks/` — one module per UI feature; each exports `register(app)` which attaches `@app.callback` handlers
- `assets/` — clientside JS for click handling (`clicked_points.js`, `bridge_selection.js`, `cursor_readout.js`)
- `utils/figures.py` — converts manifold/bridge data into Plotly figure traces
- `maps.py` — string-to-lambda parser for user-supplied map expressions

## Files with `.disabled` Extension

Several files have a `.disabled` extension (e.g., `ManifoldMachine.disabled`, `Tangle.disabled`, `clicked_points.disabled`). These are archived old implementations kept for reference — they are not imported anywhere.

## Coding Style

- **Docstrings**: Google-style with `Args:`, `Returns:`, `Raises:`, and `Note:` sections. Write them for every public class, `__init__`, and method. One-sentence private helper docstrings are fine.
- **Dev Notes**: Module-level `"""Dev Notes: ..."""` blocks are the preferred place for longer-term design questions and open issues — not inline TODOs scattered through methods.
- **Type hints**: Always annotate function signatures. Use `Literal["unstable", "stable"]` for the stability parameter. Use `Optional[X]` and `NDArray[np.float64]` from `numpy.typing` for array return types.
- **imports**: `from __future__ import annotations` at the top of every file that uses forward references. Group: stdlib, then numpy/scipy, then local `.` imports.
- **Naming**: `PascalCase` for classes, `snake_case` for everything else. `_single_leading_underscore` for private methods.
- **Avoid bare print statements**: Use `logging` (the module already configures a `NullHandler` logger in `ManifoldMachine.py`). Debug output that ends up in production code is a recurring issue — prefer `logger.debug(...)`.
- **Comments**: Sparse inline comments only when the invariant is non-obvious. Do not leave blocks of commented-out old code in committed files — the dev notes pattern or git history is the right place for those.
- **Assertions**: Used to check invariants (e.g., cdist ordering after merge). Keep them; they are the primary correctness guard here.
