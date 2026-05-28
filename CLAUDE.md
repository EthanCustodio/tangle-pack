# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Project Is

`tanglepack` is a Python library for computing and visualizing **heteroclinic/homoclinic tangles** — the stable and unstable manifolds of saddle fixed points in 2D discrete dynamical systems (area-preserving maps). The core scientific workflow is: define a map → find a fixed point → initialize manifolds → grow them iteratively → detect intersections → extract bridges and resonance zones.

There are two frontends built on top of the core library:
- **`tanglepack_gui`**: A desktop app using PySide6 + pyqtgraph
- **`tanglepack_webdash`**: A browser-based app using Plotly Dash

The desktop app is depcricated and should not be considered further

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

**Run the desktop GUI:**
```bash
tanglepack-gui
# or: python -m tanglepack_gui.app
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

The `Tangle` class uses an `rtree` spatial index. Each adjacent pair of `Point`s on a manifold is registered as a `_Segment`. Candidate intersecting segments are found via bounding-box queries, then exact intersection is computed. Results stored in `_intersecting_segments` (set of `(seg_id1, seg_id2)` pairs), `_intersecting_coords`, and `_intersecting_points`. Always call `Tangle.clear_all()` before recomputing to avoid stale references.

## Web Dashboard Architecture (`src/tanglepack_webdash/`)

Built with Plotly Dash. Per-session state is stored in `sessions.py` (`WBState` dataclass in an in-memory `_REGISTRY` dict). Each browser session gets its own `TangleWorkbench` instance.

- `layout/` — Dash component trees (sections for fixed point, manifolds, bridges, orientation, etc.)
- `callbacks/` — one module per UI feature; each exports `register(app)` which attaches `@app.callback` handlers
- `assets/` — clientside JS for click handling (`clicked_points.js`, `bridge_selection.js`, `cursor_readout.js`)
- `utils/figures.py` — converts manifold/bridge data into Plotly figure traces
- `maps.py` — string-to-lambda parser for user-supplied map expressions

## Desktop GUI Architecture (`src/tanglepack_gui/`)

Built with PySide6 + pyqtgraph. `MainWindow` owns a `Canvas` (pyqtgraph `PlotWidget`) and a side dock with controls. Map expressions are parsed from text input via `sympy` (`adapters/map_parser.py`). `adapters/workbench_view.py` converts `TangleWorkbench` manifold data into numpy arrays for plotting.

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
