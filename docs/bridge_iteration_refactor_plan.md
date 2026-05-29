# Bridge Iteration Pipeline — Implementation Plan

## What We Are Building

The goal is to make the following script workflow natural and fully public:

```python
wb = tanglepack.TangleWorkbench(f, f_inv)
fp = wb.construct_fixed_point([4, -4])
wb.orient_eigenvectors(fp, {...})
wb.initialize_both_manifolds(fp)
wb.grow_n_times(fp, "unstable", 6)
wb.grow_until_turnaround(fp, "stable")

# Setup: index manifolds, detect crossings, cut into bridges
wb.compute_intersections(fp)
wb.trim_stable_manifolds(fp)
bridges = wb.create_bridges(fp)

# Iterate one bridge — fully automatic
new_bridges = wb.iterate_bridge(bridges[2])

# Or iterate everything that hasn't been iterated yet
more_bridges = wb.iterate_all_bridges()

# Plot it all
wb.plot_tangle(fp, "stable", color="r")
wb.plot_all_bridges()
wb.plot_intersections(fp)
plt.show()
```

**Problems being fixed:**

1. `workbench._man_machine.iterate_bridge(b)` — private internals exposed in user scripts.
2. Iterated bridges are disconnected from the tangle; no automatic re-intersection or bridge splitting.
3. No `Intersection` class — intersection data is scattered across three parallel internal dicts in `Tangle`.
4. `Bridge` has no `iterated` flag, no parent/child links, no manifold-order chain.
5. Bridge registry lives inside `Tangle.bridges` (wrong layer — `Tangle` is a geometry indexer).

---

## Implementation Order

**Follow this order strictly.** Each step only depends on what came before it.

1. `src/tanglepack/Intersection.py` — new file, no deps to update
2. `src/tanglepack/Bridge.py` — add fields, import `Intersection`
3. `src/tanglepack/Tangle.py` — use `Intersection`, add incremental methods
4. `src/tanglepack/TangleWorkbench.py` — add bridge registry and public API
5. `src/tanglepack/__init__.py` — add exports
6. `scripts/tangle_workbench_test.py` — use new API

---

## File 1: `src/tanglepack/Intersection.py` (NEW FILE)

This is a clean, standalone dataclass. No circular imports.

### Purpose

Represents a single crossing between the stable and unstable manifolds. Can be either:
- **Geometric**: backed by two detected crossing segments (has `seg_ids`)
- **Synthetic**: user-created or derived from the fixed point itself (has `seg_ids = None`)

Decoupled from `BranchPoint` — `BranchPoint` is a geometry node in the manifold linked list. `Intersection` is the user-facing, attribute-carrying object.

### Full implementation

```python
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from numpy.typing import NDArray


@dataclass(eq=False)
class Intersection:
    """
    Represents a single crossing between the stable and unstable manifolds.

    Two instances at the same coordinates are distinct objects (eq=False keeps
    identity-based equality and hash, so Intersections can live in sets/dicts).

    Attributes:
        coords: Geometric (x, y) of the crossing.
        unstable_cdist: Position along the unstable manifold at this crossing.
        stable_cdist: Position along the stable manifold at this crossing.
        seg_ids: The pair of R-tree segment IDs that produced this crossing.
            None for synthetic intersections.
        label: Optional human-readable name.
    """

    coords: tuple[float, float]
    unstable_cdist: float
    stable_cdist: float
    seg_ids: Optional[frozenset[int]] = None
    label: Optional[str] = None

    # --- constructors ---

    @classmethod
    def from_segments(
        cls,
        coords: tuple[float, float],
        unstable_cdist: float,
        stable_cdist: float,
        seg1_id: int,
        seg2_id: int,
        label: Optional[str] = None,
    ) -> Intersection:
        """Create an Intersection backed by two R-tree segment IDs."""
        return cls(coords, unstable_cdist, stable_cdist,
                   frozenset({seg1_id, seg2_id}), label)

    @classmethod
    def synthetic(
        cls,
        coords: tuple[float, float],
        unstable_cdist: float,
        stable_cdist: float,
        label: Optional[str] = None,
    ) -> Intersection:
        """
        Create an Intersection not backed by a detected segment crossing.

        Use this for:
        - The fixed point itself
        - Manually specified turning points
        - Any crossing you want to declare programmatically
        """
        return cls(coords, unstable_cdist, stable_cdist, None, label)

    # --- helpers ---

    @property
    def is_synthetic(self) -> bool:
        """True if this intersection was not detected from crossing segments."""
        return self.seg_ids is None

    def get_point(self) -> NDArray[np.float64]:
        """Return coords as a (2,) array, consistent with Point.get_point()."""
        return np.array(self.coords, dtype=np.float64)
```

### Notes

- `@dataclass(eq=False)` — two Intersection objects at the same location are different
  objects (e.g., different orbit iterates hitting the same region). Identity-based equality
  is correct. They can safely be used as dict keys or in sets.
- Do NOT import `BranchPoint` here. `Intersection` is intentionally decoupled from the
  manifold geometry layer.

---

## File 2: `src/tanglepack/Bridge.py` (MODIFY)

### Imports to add

```python
from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from .Intersection import Intersection
```

Use `TYPE_CHECKING` guard to avoid a circular import at runtime (Tangle imports Bridge,
Bridge would import Intersection, which is fine — but being explicit is good practice).
At runtime these are just `Optional[...]` annotations on instance attributes, not enforced.

### `__init__` — add fields after `super().__init__()`

Do NOT change the method signature. All new fields are set as instance attributes after
the super call:

```python
# --- bridge metadata ---
self.iterated: bool = False
self.parent: Optional[Bridge] = None
self.children: list[Bridge] = []          # bridges produced by iterating this one
self.next_bridge: Optional[Bridge] = None  # next bridge along the unstable manifold (by cdist)
self.prev_bridge: Optional[Bridge] = None  # prev bridge along the unstable manifold (by cdist)
self.intersection_in: Optional[Intersection] = None   # left bounding intersection
self.intersection_out: Optional[Intersection] = None  # right bounding intersection
```

**Why `children: list[Bridge] = []` is safe here**: this is set per-instance in `__init__`,
not as a class-level default. No shared-list bug.

### `map_forward()` — update docstring only

Replace the `pass` body with:
```python
def map_forward(self):
    """
    Not implemented here; iteration is handled by TangleWorkbench.iterate_bridge().
    ManifoldMachine owns the map logic and BranchPoint insertion.
    """
```

Leave `pass` as the body. The docstring now explains why it's empty.

---

## File 3: `src/tanglepack/Tangle.py` (MODIFY)

This is the most substantial change. We are adding three things:
1. `Intersection` objects stored alongside the existing `_intersecting_coords` / `_intersecting_points` dicts.
2. A `_processed_pairs` set so `populate_intersection_dict()` becomes incremental (idempotent on already-processed pairs).
3. A new `populate_intersections_for_manifold()` that resolves only crossings involving a specific manifold.
4. A `for_manifold` parameter on `create_bridges()`.

We are **not** breaking any existing call sites. All existing behavior is preserved.

### Imports to add

```python
from .Intersection import Intersection
```

### `__init__` — add new attributes, remove stale ones

**Add:**
```python
self._intersections: list[Intersection] = []
self._intersection_by_seg: dict[int, Intersection] = {}
self._processed_pairs: set[frozenset[int]] = set()
```

**Remove:**
```python
self.bridges = None   # bridge registry moves to TangleWorkbench
self.graph = None     # never used here
```

### `clear_all()` — add clears for new attributes

After the existing `.clear()` calls, add:
```python
self._intersections.clear()
self._intersection_by_seg.clear()
self._processed_pairs.clear()
```

Remove the stale comment `# Note: bridges are kept as they may still be valid`.

### `populate_intersection_dict()` — make incremental, create `Intersection` objects

**Changes:**

1. At the top of the loop, skip already-processed pairs:
   ```python
   if seg_id_pair in self._processed_pairs:
       continue
   ```

2. After computing `branch_point` (which we keep for backward compat), also create an
   `Intersection` object:
   ```python
   intersection = Intersection.from_segments(
       coords=tuple(point),
       unstable_cdist=unstable_cdist,
       stable_cdist=stable_cdist,
       seg1_id=seg1_id,
       seg2_id=seg2_id,
   )
   self._intersections.append(intersection)
   self._intersection_by_seg[seg1_id] = intersection
   self._intersection_by_seg[seg2_id] = intersection
   ```

3. At the bottom of the loop body, mark the pair as processed:
   ```python
   self._processed_pairs.add(seg_id_pair)
   ```

The existing `_intersecting_coords` and `_intersecting_points` population remains unchanged
for backward compatibility.

### New method: `populate_intersections_for_manifold()`

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

        point = self._find_true_intersection(seg_1, seg_2)

        seg_1_cdist = 0.5 * (
            seg_1.p0.get_cdist(seg_1.manifold.stability)
            + seg_1.p0_seg1.get_cdist(seg_1.manifold.stability)
        )
        seg_2_cdist = 0.5 * (
            seg_2.p0.get_cdist(seg_2.manifold.stability)
            + seg_2.p0_seg1.get_cdist(seg_2.manifold.stability)
        )

        unstable_cdist = seg_1_cdist if seg_1.manifold.stability == "unstable" else seg_2_cdist
        stable_cdist = seg_1_cdist if seg_1.manifold.stability == "stable" else seg_2_cdist

        # keep existing BranchPoint creation for backward compat
        branch_point = BranchPoint(2, (unstable_cdist, stable_cdist), point[0], point[1])
        self._intersecting_coords[seg1_id] = point
        self._intersecting_coords[seg2_id] = point
        self._intersecting_points[seg1_id] = branch_point
        self._intersecting_points[seg2_id] = branch_point

        # new Intersection object
        intersection = Intersection.from_segments(
            coords=tuple(point),
            unstable_cdist=unstable_cdist,
            stable_cdist=stable_cdist,
            seg1_id=seg1_id,
            seg2_id=seg2_id,
        )
        self._intersections.append(intersection)
        self._intersection_by_seg[seg1_id] = intersection
        self._intersection_by_seg[seg2_id] = intersection
        new_intersections.append(intersection)

        self._processed_pairs.add(seg_id_pair)

    return new_intersections
```

### `create_bridges()` — add `for_manifold` parameter, wire `Intersection` refs, wire chain

**Signature change:**
```python
def create_bridges(
    self, for_manifold: Optional[BaseManifold] = None
) -> list[Bridge]:
```

**Filtering logic** — add at the start of the pair-iteration loop:
```python
manifold_seg_ids = (
    self._manifold_segs.get(for_manifold, set()) if for_manifold is not None else None
)

for sid_pair in self._intersecting_segments:
    if manifold_seg_ids is not None and not (sid_pair & manifold_seg_ids):
        continue   # skip pairs not involving the target manifold
    # ... rest of existing loop unchanged ...
```

**Wire `intersection_in` / `intersection_out`** — after constructing each bridge in the loop:
```python
bridge.intersection_in = self._intersection_by_seg.get(seg1.id)
bridge.intersection_out = self._intersection_by_seg.get(seg2.id)
```

**Wire `next_bridge` / `prev_bridge` chain** — after the bridge-creation loop,
before the return statement:
```python
bridge_list = list(bridges.values())
for i in range(len(bridge_list) - 1):
    bridge_list[i].next_bridge = bridge_list[i + 1]
    bridge_list[i + 1].prev_bridge = bridge_list[i]
```

**Remove the internal bridge storage** — delete this block at the end:
```python
# REMOVE:
if self.bridges is None:
    self.bridges = bridges
else:
    merged = self.bridges.copy()
    merged.update(bridges)
    self.bridges = merged

return list(self.bridges.values())
```

**Replace with:**
```python
return bridge_list
```

### Summary of all Tangle changes

| What | Type | Notes |
|---|---|---|
| `_intersections` | new attr | `list[Intersection]` |
| `_intersection_by_seg` | new attr | `dict[int, Intersection]` |
| `_processed_pairs` | new attr | `set[frozenset[int]]` |
| `self.bridges` | removed | lives in `TangleWorkbench` |
| `self.graph` | removed | never used |
| `clear_all()` | modified | clear three new attrs |
| `populate_intersection_dict()` | modified | incremental + create Intersection |
| `populate_intersections_for_manifold()` | new method | targeted incremental resolve |
| `create_bridges()` | modified | `for_manifold` param, wire Intersection refs, wire chain, remove self.bridges storage |

---

## File 4: `src/tanglepack/TangleWorkbench.py` (MODIFY)

### Imports to add

```python
from .Bridge import Bridge
from .Intersection import Intersection
```

(Both may already be transitively imported, but be explicit.)

### `__init__` — add bridge registry

```python
self._bridges: list[Bridge] = []
```

### `create_bridges()` — populate registry

```python
def create_bridges(self, fixed_point: FixedPoint) -> list[Bridge]:
    bridges = self.Tangle.create_bridges()   # for_manifold=None → all intersections
    self._bridges.extend(bridges)
    return bridges
```

### New property: `uniiterated_bridges`

```python
@property
def uniiterated_bridges(self) -> list[Bridge]:
    """All bridges that have not yet been iterated forward."""
    return [b for b in self._bridges if not b.iterated]
```

### New method: `iterate_bridge()`

```python
def iterate_bridge(self, bridge: Bridge) -> list[Bridge]:
    """
    Map a bridge forward one iterate, add the result to the tangle, detect new
    intersections with the stable manifold, cut the result into new bridges,
    and return those bridges.

    Marks the original bridge as iterated and wires parent/child links.

    Args:
        bridge: A bridge created by create_bridges() or a previous iterate_bridge().

    Returns:
        List of new Bridge objects from cutting the iterated result.
        If the iterated bridge makes no new crossings, returns a single-element
        list containing the unsplit iterated bridge.

    Raises:
        ValueError: If bridge has already been iterated.
        ValueError: If create_bridges() has not been called yet.
    """
    if bridge.iterated:
        raise ValueError(
            "This bridge has already been iterated. Check bridge.children for the results."
        )
    if not self._bridges:
        raise ValueError(
            "No bridges registered. Call create_bridges() before iterate_bridge()."
        )

    # 1. map forward
    iterated = self._man_machine.iterate_bridge(bridge)

    # 2. add to tangle (stable manifold already indexed from compute_intersections)
    self.Tangle.add_manifold(iterated)

    # 3. resolve only new crossings involving the iterated bridge
    new_intersections = self.Tangle.populate_intersections_for_manifold(iterated)

    # 4. cut at crossings
    if new_intersections:
        new_bridges = self.Tangle.create_bridges(for_manifold=iterated)
    else:
        # no crossings — iterated bridge is already a valid unsplit bridge
        new_bridges = [iterated]

    # 5. wire genealogy
    bridge.iterated = True
    bridge.children = new_bridges
    for nb in new_bridges:
        nb.parent = bridge

    # 6. register
    self._bridges.extend(new_bridges)

    return new_bridges
```

### New method: `iterate_all_bridges()`

```python
def iterate_all_bridges(self) -> list[Bridge]:
    """
    Iterate all bridges that have not yet been mapped forward.

    Returns:
        All new bridges produced across all iterations.
    """
    pending = list(self.uniiterated_bridges)  # snapshot before loop mutates _bridges
    all_new: list[Bridge] = []
    for bridge in pending:
        all_new.extend(self.iterate_bridge(bridge))
    return all_new
```

### `plot_all_bridges()` — make registry-aware

```python
def plot_all_bridges(self, bridges: Optional[list[Bridge]] = None) -> None:
    """
    Plot a list of bridges. If no list is supplied, plots all registered bridges.

    Args:
        bridges: List of bridges to plot. Defaults to self._bridges.
    """
    if bridges is None:
        bridges = self._bridges
    n = len(bridges)
    if n == 0:
        return
    colors = cm.get_cmap("tab20", n)
    for i, bridge in enumerate(bridges):
        bridge.plot(color=colors(i))
```

### Summary of all TangleWorkbench changes

| What | Type | Notes |
|---|---|---|
| `self._bridges` | new attr | `list[Bridge]` |
| `create_bridges()` | modified | extends `_bridges`, unchanged return |
| `uniiterated_bridges` | new property | filters `_bridges` by `not iterated` |
| `iterate_bridge()` | new method | full pipeline: iterate → index → intersect → cut → wire |
| `iterate_all_bridges()` | new method | applies `iterate_bridge` to all pending |
| `plot_all_bridges()` | modified | `bridges` param now optional, defaults to registry |

---

## File 5: `src/tanglepack/__init__.py` (MODIFY)

Add these two exports:

```python
from .Intersection import Intersection
from .Bridge import Bridge
```

The full file should become:

```python
from .DynamicalSystem import DynamicalSystem

from .BasePoint import BasePoint
from .Point import Point
from .BranchPoint import BranchPoint

from .FixedPointSolver import FixedPointSolver
from .FixedPoint import FixedPoint

from .BaseManifold import BaseManifold
from .ManifoldMachine import ManifoldMachine
from .ManifoldInitializer import ManifoldInitializer

from .Intersection import Intersection
from .Bridge import Bridge
from .Tangle import Tangle
from .TangleWorkbench import TangleWorkbench
```

---

## File 6: `scripts/tangle_workbench_test.py` (MODIFY)

Replace the two lines that use private internals:

```python
# BEFORE:
new_bridges = workbench._man_machine.iterate_bridge(bridges[2])
# ...
new_bridges.plot()

# AFTER:
new_bridges = workbench.iterate_bridge(bridges[2])
# ...
workbench.plot_all_bridges(new_bridges)
```

Or to use the registry-aware form:
```python
workbench.plot_all_bridges()   # plots everything including new_bridges
```

---

## Invariants and Gotchas

### The tangle must NOT be cleared after `create_bridges()`

`compute_intersections()` calls `Tangle.clear_all()` and rebuilds. After that call, the
stable manifold segments are in the R-tree. `iterate_bridge()` relies on them being there
to detect new crossings. **Do not call `compute_intersections()` again after the pipeline
has started** — doing so clears the tangle and loses all iterated bridge state.

The correct single-setup pattern is:
```python
wb.compute_intersections(fp)   # clears + rebuilds
wb.trim_stable_manifolds(fp)
bridges = wb.create_bridges(fp)
# from here on: only iterate_bridge() and iterate_all_bridges()
```

### `_processed_pairs` prevents double-resolving

Once a pair of segment IDs has been resolved to an `Intersection`, it's in
`_processed_pairs`. Calling `populate_intersection_dict()` or
`populate_intersections_for_manifold()` again with the same pairs is a no-op for those
pairs. This is what makes incremental updates safe.

### `create_bridges(for_manifold=iterated)` only cuts the iterated bridge

The `for_manifold` filter checks that at least one segment ID in each intersecting pair
belongs to `for_manifold`. Old bridge-manifold segments won't appear in new pairs, so
old bridges are never re-processed.

### `iterate_bridge()` returns `[iterated]` if no crossings

If the iterated bridge doesn't cross the stable manifold (possible in early iterates),
`new_bridges = [iterated]` is returned as a single-element list. The iterated bridge is
still registered in `_bridges`. Subsequent `iterate_bridge(new_bridges[0])` will continue
the chain.

### `children: list[Bridge] = []` is per-instance

Because this default is set inside `__init__`, not at the class level, each Bridge gets its
own list. There is no shared-list class-attribute bug here.

### Stale `BranchPoint` objects in `_intersecting_points`

`_intersecting_points` continues to hold `BranchPoint` markers purely for backward
compatibility (e.g., any existing test or code that accesses it directly). These are NOT
inserted into the manifold linked lists — they are just coordinate/cdist carriers. The new
`Intersection` objects in `_intersection_by_seg` are the preferred way to access this data
going forward.

---

## Complete "After" Script Example

```python
import tanglepack, numpy as np
import matplotlib.pyplot as plt


def henon_map(point):
    k, b = (10, 1)
    x, y = point
    return np.array([y - k + x**2, -b * x])


def henon_map_inverse(point):
    k, b = (10, 1)
    x, y = point
    return np.array([-y / b, x + k - (y**2) / (b**2)])


wb = tanglepack.TangleWorkbench(henon_map, henon_map_inverse)
fp = wb.construct_fixed_point([4, -4])
wb.orient_eigenvectors(fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])})
wb.initialize_both_manifolds(fp)
wb.grow_n_times(fp, "unstable", num_iterations=6)
wb.grow_until_turnaround(fp, "stable")

# One-time setup: index, intersect, cut
wb.compute_intersections(fp)
wb.trim_stable_manifolds(fp)
bridges = wb.create_bridges(fp)
print(f"{len(bridges)} initial bridges")

# Iterate a specific bridge
new_bridges = wb.iterate_bridge(bridges[2])
print(f"bridge[2] produced {len(new_bridges)} child bridges")

# Check genealogy
assert bridges[2].iterated
assert bridges[2].children == new_bridges
assert new_bridges[0].parent is bridges[2]

# Iterate everything else in one call
more_bridges = wb.iterate_all_bridges()
print(f"iterate_all produced {len(more_bridges)} additional bridges")

# Create a synthetic intersection (not from a segment crossing)
synth = tanglepack.Intersection.synthetic(
    coords=(1.5, -0.3),
    unstable_cdist=12.4,
    stable_cdist=3.1,
    label="manually placed",
)

# Plot
plt.figure()
wb.plot_tangle(fp, "stable", color="r")
wb.plot_all_bridges()           # uses internal registry
wb.plot_intersections(fp)
plt.xlim([-15, 15])
plt.ylim([-15, 15])
plt.show()
```
