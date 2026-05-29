# Intersection Registry and Topological Layer — Implementation Plan

## Overview and Goals

After generating enough intersections via bridge iteration, the workflow transitions from
manifold numerics (growing linked lists, managing segments, R-tree geometry) to topology
(the combinatorial structure of how intersection points map to each other, how they're
ordered on each manifold, and what that implies about transport and mixing).

The topological layer must be:
- **Separable** — once the registry is populated, `TangleWorkbench` and `Tangle` are
  no longer needed. The topology can be studied from the registry alone.
- **Connected** — the registry is populated automatically as bridges are iterated.
- **Easy to query** — `registry[id]` → `Intersection`; `table[3, 2]` → id of 2nd
  forward iterate of intersection 3; sorted orderings on W^u and W^s immediately accessible.
- **Extensible with synthetic intersections** — you can add intersections that weren't
  found geometrically (e.g., the fixed point itself, manually placed reference points).

---

## Mathematical Background (the structures we're capturing)

### What intersections are

Every point p in W^u ∩ W^s maps to another intersection under f:
- f(p) is also in W^u ∩ W^s (f preserves both manifolds)
- The iterate map on intersections is well-defined and discrete

### How cdists transform under iteration

For a period-1 saddle fixed point with unstable eigenvalue λ_u > 1:
```
unstable_cdist(f^n(p)) = λ_u^n  ×  unstable_cdist(p)   # grows under forward map
stable_cdist(f^n(p))   = λ_u^{-n} × stable_cdist(p)    # shrinks under forward map
```

This is the **key formula**. Given any intersection p with known cdists, we can predict
exactly where f^n(p) lives in (unstable_cdist, stable_cdist) space. If that predicted
point has already been computed and registered, we can link them in the iterate table
without any additional manifold numerics.

### The two orderings

Intersections live on both manifolds simultaneously. They can be sorted two ways:
- **Unstable order** — sorted by `unstable_cdist`. This is the sequence of crossings
  encountered walking along W^u away from the fixed point.
- **Stable order** — sorted by `stable_cdist`. This is the sequence of crossings
  encountered walking along W^s away from the fixed point.

These two orderings define the combinatorial template of the tangle. The permutation
that takes you from unstable order to stable order encodes the topological type of
the tangle and drives topological entropy calculations.

### What the iterate table captures

If you label intersections 0, 1, 2, ... by their position in unstable order, the
iterate map permutes these labels. The iterate table `F[id, n]` is a record of
which intersection ID maps to which. The stable-order ranking of the image tells
you the "symbolic dynamics" — the fundamental object in the topological analysis.

---

## Architecture Decision: Two Layers

```
┌──────────────────────────────────────────────────────────┐
│  NUMERIC LAYER                                           │
│  TangleWorkbench, ManifoldMachine, Tangle                │
│  - Grows manifolds as linked lists                       │
│  - Detects geometric crossings (R-tree)                  │
│  - Cuts bridges                                          │
│  - Produces Intersection objects with coords + cdists    │
└───────────────────────┬──────────────────────────────────┘
                        │  wb.populate_registry()
                        │  (or auto-populated during iterate_bridge)
┌───────────────────────▼──────────────────────────────────┐
│  TOPOLOGICAL LAYER                                       │
│  IntersectionRegistry, IterateTable                      │
│  - All intersections with unique IDs                     │
│  - IterateTable: registry[id, n] → id                    │
│  - Sorted orderings on W^u and W^s                       │
│  - Lobe structure (future)                               │
│  - Markov matrix (future)                                │
│  - No manifold objects needed from here on               │
└──────────────────────────────────────────────────────────┘
```

The separation point is `IntersectionRegistry`. Everything below it is geometry.
Everything in it and above is topology.

---

## Update: `src/tanglepack/Intersection.py`

One field is added to the existing `Intersection` dataclass: `id`. The registry
assigns this value after insertion. It is `None` until the intersection has been
registered.

### Change: add `id` field

```python
@dataclass(eq=False)
class Intersection:
    coords: tuple[float, float]
    unstable_cdist: float
    stable_cdist: float
    seg_ids: Optional[frozenset[int]] = None
    label: Optional[str] = None
    id: Optional[int] = None          # assigned by IntersectionRegistry.add()
```

No other changes. All existing constructors remain valid. The field is last so no
call sites break.

---

## New File: `src/tanglepack/IterateTable.py`

### Purpose

A 2D lookup structure supporting `table[intersection_id, n]` → `intersection_id`,
where `n` is the iterate depth:
- `n > 0`: forward iterates (under f^n)
- `n < 0`: backward iterates (under f^{-n})
- `n = 0`: identity (returns the same id)

This is **two logical structures in one** (the forward iterate map and the backward
iterate map), accessible through a single unified interface.

### Full implementation

```python
from __future__ import annotations

from typing import Optional
import numpy as np
from numpy.typing import NDArray


class IterateTable:
    """
    2D lookup: (intersection_id, n) → intersection_id.

    Supports:
        table[3, 2]  → ID of f^2(intersection 3)
        table[3, -1] → ID of f^{-1}(intersection 3)
        table[3, 0]  → 3  (identity)

    Entries are None until set. Setting one direction auto-sets the reverse:
        table[3, 2] = 7  also records  table[7, -2] = 3

    Attributes:
        _forward: dict[int, dict[int, int]]  — _forward[id][n] = target_id  (n > 0)
        _backward: dict[int, dict[int, int]] — _backward[id][n] = target_id (n > 0)
    """

    def __init__(self):
        self._forward: dict[int, dict[int, int]] = {}
        self._backward: dict[int, dict[int, int]] = {}

    # --- primary interface ---

    def __getitem__(self, key: tuple[int, int]) -> Optional[int]:
        """
        table[id, n] → ID of n-th iterate of intersection `id`.
        Returns None if the iterate is not yet recorded.
        """
        source_id, n = key
        if n == 0:
            return source_id
        if n > 0:
            return self._forward.get(source_id, {}).get(n)
        else:
            return self._backward.get(source_id, {}).get(-n)

    def __setitem__(self, key: tuple[int, int], target_id: int):
        """
        table[id, n] = target_id  records f^n(id) = target_id.
        Also records the reverse: f^{-n}(target_id) = id.
        """
        source_id, n = key
        if n == 0:
            return  # identity — nothing to store
        if n > 0:
            self._forward.setdefault(source_id, {})[n] = target_id
            self._backward.setdefault(target_id, {})[n] = source_id
        else:
            self._backward.setdefault(source_id, {})[-n] = target_id
            self._forward.setdefault(target_id, {})[-n] = source_id

    def __contains__(self, key: tuple[int, int]) -> bool:
        """Check if table[id, n] is recorded."""
        return self[key] is not None

    # --- bulk operations ---

    def forward_depth(self, source_id: int) -> int:
        """Maximum forward iterate depth recorded for this intersection."""
        return max(self._forward.get(source_id, {}).keys(), default=0)

    def backward_depth(self, source_id: int) -> int:
        """Maximum backward iterate depth recorded for this intersection."""
        return max(self._backward.get(source_id, {}).keys(), default=0)

    def forward_chain(self, source_id: int) -> list[int]:
        """Return [id, f(id), f^2(id), ...] for all recorded forward iterates."""
        chain = [source_id]
        n = 1
        while (nxt := self[source_id, n]) is not None:
            chain.append(nxt)
            n += 1
        return chain

    def backward_chain(self, source_id: int) -> list[int]:
        """Return [id, f^{-1}(id), f^{-2}(id), ...] for all recorded backward iterates."""
        chain = [source_id]
        n = 1
        while (prev := self[source_id, -n]) is not None:
            chain.append(prev)
            n += 1
        return chain

    def all_registered_ids(self) -> set[int]:
        """All intersection IDs that appear in any entry of this table."""
        ids = set(self._forward.keys()) | set(self._backward.keys())
        for sub in self._forward.values():
            ids.update(sub.values())
        for sub in self._backward.values():
            ids.update(sub.values())
        return ids

    # --- dense array export ---

    def as_forward_array(
        self, ids: list[int], max_depth: int
    ) -> NDArray[np.int64]:
        """
        Dense 2D array A where A[i, d-1] = ID of f^d(ids[i]), or -1 if unknown.

        Args:
            ids: Ordered list of intersection IDs (defines row order).
            max_depth: Number of forward iterate columns.

        Returns:
            Array of shape (len(ids), max_depth), dtype int64.
            -1 = iterate not yet recorded.

        Example:
            arr = table.as_forward_array(registry.all_ids(), max_depth=5)
            arr[3, 1]  # ID of f^2(ids[3])  (column 1 = depth 2)
        """
        arr = np.full((len(ids), max_depth), fill_value=-1, dtype=np.int64)
        id_to_row = {iid: i for i, iid in enumerate(ids)}
        for i, iid in enumerate(ids):
            for d in range(1, max_depth + 1):
                target = self[iid, d]
                if target is not None:
                    arr[i, d - 1] = target
        return arr

    def as_backward_array(
        self, ids: list[int], max_depth: int
    ) -> NDArray[np.int64]:
        """
        Same as as_forward_array but for backward iterates (f^{-d}).
        B[i, d-1] = ID of f^{-d}(ids[i]), or -1 if unknown.
        """
        arr = np.full((len(ids), max_depth), fill_value=-1, dtype=np.int64)
        for i, iid in enumerate(ids):
            for d in range(1, max_depth + 1):
                target = self[iid, -d]
                if target is not None:
                    arr[i, d - 1] = target
        return arr

    def register_iterate(self, source_id: int, n: int, target_id: int):
        """
        Explicit named method — delegates to __setitem__ for readability.
        Prefer this in code that populates the table; prefer [] in code that reads.
        """
        self[source_id, n] = target_id
```

### Notes

- **Automatic reverse recording**: setting `table[3, 2] = 7` also sets `table[7, -2] = 3`.
  This means forward and backward tables stay consistent without manual effort.
- **`as_forward_array`** uses depth index `d-1` so that column 0 = first forward iterate,
  matching the natural "n=1, n=2, ..." language. This is different from `table[id, n]`
  where n starts at 1. Be careful of the off-by-one when converting between the two.
  Consider renaming: column `d` in the array = `table[id, d+1]`.
- **No numpy dependency at instantiation** — only needed for the `as_*_array` methods.

---

## New File: `src/tanglepack/IntersectionRegistry.py`

### Purpose

The master store of all intersection points in the tangle. Provides:
- Unique auto-incrementing integer IDs
- Collision detection so duplicate intersections (same cdist pair within tolerance)
  are not registered twice
- Fast lookup by ID: `registry[id]` → `Intersection`
- Two sorted orderings: by unstable cdist, by stable cdist
- An `IterateTable` instance for the iterate map
- A method to infer the iterate map automatically from cdist relationships and a
  known eigenvalue — no need to track iterate links through the bridge code

### Full implementation

```python
from __future__ import annotations

import bisect
from typing import Optional

import numpy as np
from numpy.typing import NDArray

from .Intersection import Intersection
from .IterateTable import IterateTable


class IntersectionRegistry:
    """
    Master store of all intersection points in the tangle.

    Design:
        IDs are contiguous integers 0, 1, 2, ..., N-1 assigned in insertion order.
        If a new intersection collides (same cdists within tolerance) with an existing
        one, the existing ID is returned and no duplicate is stored.

    Primary interface:
        registry.add(intersection)           → int (the assigned ID)
        registry[id]                         → Intersection
        registry.iterate_table[id, n]        → int or None
        registry.by_unstable_cdist           → list[int] (IDs sorted by u-cdist)
        registry.by_stable_cdist             → list[int] (IDs sorted by s-cdist)

    Attributes:
        _store: dict[int, Intersection]
        _next_id: int
        cdist_tol: float — collision tolerance in (unstable_cdist, stable_cdist) space
        iterate_table: IterateTable
        _unstable_order: list[int] — IDs sorted by unstable_cdist (maintained sorted)
        _stable_order: list[int]   — IDs sorted by stable_cdist (maintained sorted)
    """

    def __init__(self, cdist_tol: float = 1e-6):
        self._store: dict[int, Intersection] = {}
        self._next_id: int = 0
        self.cdist_tol = cdist_tol
        self.iterate_table = IterateTable()

        # Sorted lists of IDs maintained in cdist order.
        # Each list stores IDs; lookups sort by the corresponding cdist value.
        self._unstable_order: list[int] = []   # sorted ascending by unstable_cdist
        self._stable_order: list[int] = []     # sorted ascending by stable_cdist

        # Fast collision lookup: (rounded_u, rounded_s) → id
        # This is a secondary index; _store is authoritative.
        self._cdist_index: dict[tuple[float, float], int] = {}

    # --- core insert / lookup ---

    def add(self, intersection: Intersection) -> int:
        """
        Register an intersection, returning its unique ID.

        If a collision is detected (another intersection with cdists within
        self.cdist_tol), the existing ID is returned and intersection is NOT added.
        Otherwise a new ID is assigned, intersection.id is set, and the intersection
        is stored.

        Args:
            intersection: The Intersection to register.

        Returns:
            The integer ID (new or existing on collision).
        """
        existing = self._find_collision(intersection)
        if existing is not None:
            return existing

        new_id = self._next_id
        self._next_id += 1
        intersection.id = new_id
        self._store[new_id] = intersection

        self._insert_into_unstable_order(new_id, intersection.unstable_cdist)
        self._insert_into_stable_order(new_id, intersection.stable_cdist)

        # update secondary collision index
        key = self._cdist_key(intersection)
        self._cdist_index[key] = new_id

        return new_id

    def add_synthetic(
        self,
        coords: tuple[float, float],
        unstable_cdist: float,
        stable_cdist: float,
        label: Optional[str] = None,
    ) -> int:
        """
        Convenience method to add an intersection that is not from a detected
        segment crossing — e.g., a manually placed reference point or the fixed
        point itself.
        """
        return self.add(Intersection.synthetic(coords, unstable_cdist, stable_cdist, label))

    def __getitem__(self, id: int) -> Intersection:
        """registry[id] → Intersection. Raises KeyError if id is unknown."""
        return self._store[id]

    def __len__(self) -> int:
        return len(self._store)

    def __contains__(self, id: int) -> bool:
        return id in self._store

    def __iter__(self):
        """Iterate over all (id, Intersection) pairs in insertion order."""
        return iter(self._store.items())

    # --- ordering views ---

    @property
    def by_unstable_cdist(self) -> list[int]:
        """IDs sorted ascending by unstable_cdist. Updated on every add()."""
        return list(self._unstable_order)

    @property
    def by_stable_cdist(self) -> list[int]:
        """IDs sorted ascending by stable_cdist. Updated on every add()."""
        return list(self._stable_order)

    def unstable_rank(self, id: int) -> int:
        """Position of intersection `id` in the W^u ordering (0-based)."""
        return self._unstable_order.index(id)

    def stable_rank(self, id: int) -> int:
        """Position of intersection `id` in the W^s ordering (0-based)."""
        return self._stable_order.index(id)

    def all_ids(self) -> list[int]:
        """All registered IDs in insertion order."""
        return list(self._store.keys())

    # --- cdist-based lookup ---

    def find_by_cdist(
        self,
        unstable_cdist: float,
        stable_cdist: float,
        tol: Optional[float] = None,
    ) -> Optional[int]:
        """
        Find the ID of an intersection whose cdists are within `tol` of the
        given values. Returns None if no match found.

        Args:
            unstable_cdist: Target unstable canonical distance.
            stable_cdist: Target stable canonical distance.
            tol: Search tolerance. Defaults to self.cdist_tol.

        Returns:
            Matching intersection ID, or None.
        """
        if tol is None:
            tol = self.cdist_tol

        for id, existing in self._store.items():
            if (abs(existing.unstable_cdist - unstable_cdist) < tol
                    and abs(existing.stable_cdist - stable_cdist) < tol):
                return id
        return None

    # --- iterate table population ---

    def infer_iterates(
        self,
        lambda_u: float,
        max_depth: int = 10,
        tol_multiplier: float = 10.0,
    ) -> int:
        """
        Scan all registered intersections and infer iterate relationships by
        predicting the cdists of their forward/backward iterates and checking
        for matches in the registry.

        This is the PRIMARY METHOD for filling in the iterate table. Call it
        after each batch of bridge iterations to link up newly computed
        intersections.

        Under forward iteration for a period-1 saddle with eigenvalue lambda_u:
            unstable_cdist(f^n(p)) = lambda_u^n * unstable_cdist(p)
            stable_cdist(f^n(p))   = lambda_u^{-n} * stable_cdist(p)

        Args:
            lambda_u: Unstable eigenvalue (> 1).
            max_depth: How many iterate levels to search (forward and backward).
            tol_multiplier: Scale factor on self.cdist_tol for cdist matching.
                Use > 1 to account for numerical error accumulation over many iterates.

        Returns:
            Number of new iterate relationships recorded.
        """
        tol = self.cdist_tol * tol_multiplier
        recorded = 0

        for source_id, source in self._store.items():
            for n in range(1, max_depth + 1):
                # forward
                if (source_id, n) not in self.iterate_table:
                    pred_u = (lambda_u ** n) * source.unstable_cdist
                    pred_s = source.stable_cdist / (lambda_u ** n)
                    target_id = self.find_by_cdist(pred_u, pred_s, tol)
                    if target_id is not None:
                        self.iterate_table.register_iterate(source_id, n, target_id)
                        recorded += 1

                # backward (no need to search separately — auto-linked by IterateTable)
                # but check if the backward entry is already there from a forward link
                if (source_id, -n) not in self.iterate_table:
                    pred_u = source.unstable_cdist / (lambda_u ** n)
                    pred_s = source.stable_cdist * (lambda_u ** n)
                    target_id = self.find_by_cdist(pred_u, pred_s, tol)
                    if target_id is not None:
                        self.iterate_table.register_iterate(source_id, -n, target_id)
                        recorded += 1

        return recorded

    def register_iterate(self, source_id: int, n: int, target_id: int):
        """
        Explicitly record that f^n(source) = target.
        Delegates to iterate_table; provided here for convenience so callers
        only need one object.
        """
        self.iterate_table.register_iterate(source_id, n, target_id)

    # --- array exports ---

    def as_forward_array(self, max_depth: int = 5) -> NDArray[np.int64]:
        """
        Dense array A[i, d-1] = ID of f^d(ids[i]) where ids = all_ids() in order.
        Shape: (N, max_depth). -1 = unknown.

        This is the 'matrix' where:
            A[3, 1] = ID of f^2(intersection with index 3 in all_ids())
        """
        return self.iterate_table.as_forward_array(self.all_ids(), max_depth)

    def as_backward_array(self, max_depth: int = 5) -> NDArray[np.int64]:
        """
        Dense array B[i, d-1] = ID of f^{-d}(ids[i]).
        Shape: (N, max_depth). -1 = unknown.
        """
        return self.iterate_table.as_backward_array(self.all_ids(), max_depth)

    def unstable_order_array(self) -> NDArray[np.float64]:
        """
        (N, 3) array where each row is [id, unstable_cdist, stable_cdist],
        sorted by unstable_cdist. Useful for analysis and export.
        """
        rows = []
        for iid in self._unstable_order:
            x = self._store[iid]
            rows.append([iid, x.unstable_cdist, x.stable_cdist])
        return np.array(rows, dtype=np.float64)

    def stable_order_array(self) -> NDArray[np.float64]:
        """
        (N, 3) array sorted by stable_cdist. Same columns as unstable_order_array.
        """
        rows = []
        for iid in self._stable_order:
            x = self._store[iid]
            rows.append([iid, x.unstable_cdist, x.stable_cdist])
        return np.array(rows, dtype=np.float64)

    # --- internal helpers ---

    def _find_collision(self, intersection: Intersection) -> Optional[int]:
        """Linear scan for an existing intersection within cdist_tol."""
        for id, existing in self._store.items():
            if (abs(existing.unstable_cdist - intersection.unstable_cdist) < self.cdist_tol
                    and abs(existing.stable_cdist - intersection.stable_cdist) < self.cdist_tol):
                return id
        return None

    def _cdist_key(self, intersection: Intersection) -> tuple[float, float]:
        """Discretized key for the secondary collision index."""
        digits = max(0, -int(np.floor(np.log10(self.cdist_tol))) - 1)
        return (
            round(intersection.unstable_cdist, digits),
            round(intersection.stable_cdist, digits),
        )

    def _insert_into_unstable_order(self, id: int, unstable_cdist: float):
        """Insert `id` into _unstable_order maintaining sort by unstable_cdist."""
        keys = [self._store[i].unstable_cdist for i in self._unstable_order]
        pos = bisect.bisect_left(keys, unstable_cdist)
        self._unstable_order.insert(pos, id)

    def _insert_into_stable_order(self, id: int, stable_cdist: float):
        """Insert `id` into _stable_order maintaining sort by stable_cdist."""
        keys = [self._store[i].stable_cdist for i in self._stable_order]
        pos = bisect.bisect_left(keys, stable_cdist)
        self._stable_order.insert(pos, id)
```

### Key design decisions

**Contiguous IDs**: IDs are 0-based contiguous integers assigned in insertion order.
On a collision, no new ID is generated. This means `all_ids()[i]` gives the i-th
registered intersection, and `as_forward_array()[3, ...]` is about the 4th intersection
registered — not the one with ID=3 (unless IDs and indices happen to match, which they
do when there are no collisions). 

Note: if you want `A[3, 1]` to mean "intersection with ID=3", use:
```python
arr = registry.as_forward_array()
id_3_row = arr[registry.all_ids().index(3), :]
```
OR: if IDs are dense (no collisions), `arr[3, :]` works directly.

**Linear collision scan**: Fine for tangles with up to a few thousand intersections.
If the registry grows very large (10k+), replace `_find_collision` with a KD-tree
on (unstable_cdist, stable_cdist).

**Sorted lists are rebuilt on every insert**: `_insert_into_unstable_order` is O(N)
due to the list-comprehension to extract keys. For typical tangle sizes (< 1000
intersections) this is negligible. If needed, replace with `sortedcontainers.SortedList`
keyed on cdist.

---

## Updated: `src/tanglepack/TangleWorkbench.py`

### New attribute in `__init__`

```python
self._intersection_registry = IntersectionRegistry()
```

### New property

```python
@property
def intersection_registry(self) -> IntersectionRegistry:
    return self._intersection_registry
```

### Updated `compute_intersections()`

After calling `self.Tangle.populate_intersection_dict()`, populate the registry:

```python
def compute_intersections(self, fixed_point: FixedPoint):
    self.Tangle.clear_all()
    self._intersection_registry = IntersectionRegistry()   # fresh registry on rebuild
    
    self.index_manifolds(fixed_point, "unstable")
    self.index_manifolds(fixed_point, "stable")
    self.Tangle.populate_intersection_dict()

    for intersection in self.Tangle._intersections:
        self._intersection_registry.add(intersection)

    return list(self.Tangle._intersecting_coords.values())
```

### Updated `iterate_bridge()`

After the new intersections are detected, add them to the registry:

```python
# ... existing steps 1-4 ...

# register new intersections
for bridge in new_bridges:
    if bridge.intersection_in is not None:
        self._intersection_registry.add(bridge.intersection_in)
    if bridge.intersection_out is not None:
        self._intersection_registry.add(bridge.intersection_out)
```

### New method: `infer_iterate_table()`

```python
def infer_iterate_table(self, fixed_point: FixedPoint, max_depth: int = 10) -> int:
    """
    Populate the registry's iterate table by predicting cdist values of forward
    and backward iterates and checking for matches in the registry.

    Call this after any batch of iterate_bridge() calls.

    Args:
        fixed_point: The fixed point whose eigenvalue drives the cdist transform.
        max_depth: How many iterate levels to search forward and backward.

    Returns:
        Number of new iterate relationships added.
    """
    lambda_u = abs(fixed_point.unstable_eigenvalues[0])
    return self._intersection_registry.infer_iterates(lambda_u, max_depth)
```

### New method: `populate_registry()` (manual rebuild)

```python
def populate_registry(self, fixed_point: FixedPoint) -> IntersectionRegistry:
    """
    Rebuild the intersection registry from scratch from the current Tangle state.
    Useful if the registry got out of sync or needs to be reset.
    """
    self._intersection_registry = IntersectionRegistry()
    for intersection in self.Tangle._intersections:
        self._intersection_registry.add(intersection)
    return self._intersection_registry
```

---

## Updated: `src/tanglepack/__init__.py`

Add exports:

```python
from .IterateTable import IterateTable
from .IntersectionRegistry import IntersectionRegistry
```

---

## Complete Usage Example

```python
import tanglepack
import numpy as np
import matplotlib.pyplot as plt


def henon_map(point):
    k, b = (10, 1)
    x, y = point
    return np.array([y - k + x**2, -b * x])

def henon_map_inverse(point):
    k, b = (10, 1)
    x, y = point
    return np.array([-y / b, x + k - (y**2) / (b**2)])


# ── Numeric phase ──────────────────────────────────────────────────────────
wb = tanglepack.TangleWorkbench(henon_map, henon_map_inverse)
fp = wb.construct_fixed_point([4, -4])
wb.orient_eigenvectors(fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])})
wb.initialize_both_manifolds(fp)
wb.grow_n_times(fp, "unstable", num_iterations=6)
wb.grow_until_turnaround(fp, "stable")

wb.compute_intersections(fp)        # also populates registry
wb.trim_stable_manifolds(fp)
bridges = wb.create_bridges(fp)

# Iterate everything three rounds
for _ in range(3):
    wb.iterate_all_bridges()

# Infer iterate table from cdist relationships
new_links = wb.infer_iterate_table(fp, max_depth=8)
print(f"Recorded {new_links} iterate relationships")


# ── Topological phase ──────────────────────────────────────────────────────
# From here, TangleWorkbench is no longer needed.
registry = wb.intersection_registry

print(f"Total intersections: {len(registry)}")

# Look up an intersection by ID
p = registry[3]
print(f"Intersection 3: coords={p.coords}, u_cdist={p.unstable_cdist:.4f}")

# Forward and backward iterate chains
fwd_chain = registry.iterate_table.forward_chain(3)
print(f"Forward chain from 3: {fwd_chain}")

# The iterate table: which intersection is the 2nd forward iterate of intersection 3?
iterate_id = registry.iterate_table[3, 2]
print(f"f^2(intersection 3) = intersection {iterate_id}")

# Backward iterates
pre_id = registry.iterate_table[3, -1]
print(f"f^{-1}(intersection 3) = intersection {pre_id}")

# Unstable and stable orderings (the two fundamental orderings)
u_order = registry.by_unstable_cdist    # [id0, id1, id2, ...] sorted by u-cdist
s_order = registry.by_stable_cdist     # same IDs in stable-cdist order
print("Unstable order:", u_order)
print("Stable order:", s_order)

# Dense numpy array for analysis: F[i, d-1] = ID of f^d(ids[i]), -1=unknown
F = registry.as_forward_array(max_depth=5)
B = registry.as_backward_array(max_depth=5)
print("Forward iterate array shape:", F.shape)
print("F[3, 1] =", F[3, 1], "(should match iterate_table[ids[3], 2])")

# Add a synthetic intersection (not from a geometric crossing)
fixed_pt_id = registry.add_synthetic(
    coords=tuple(fp.coordinates[0].ravel()),
    unstable_cdist=0.0,
    stable_cdist=0.0,
    label="fixed_point",
)
print(f"Fixed point registered as intersection {fixed_pt_id}")

# Find any intersection near a given cdist
near_id = registry.find_by_cdist(unstable_cdist=5.2, stable_cdist=1.1, tol=0.1)
```

---

## The Two Structures — Clarified

The two structures the user described:

**Structure 1 — `registry.iterate_table[id, n]`**  
The forward/backward iterate table. `n > 0` is forward, `n < 0` is backward, `n = 0` is identity. One object handles both directions.

**Structure 2 — `registry.by_unstable_cdist` and `registry.by_stable_cdist`**  
The two orderings. These are ordered lists of IDs that define the combinatorial template of the tangle. The permutation taking `by_unstable_cdist` to `by_stable_cdist` is the topological invariant that drives all further analysis (topological entropy, lobe areas, transport coefficients).

Dense array forms: `registry.unstable_order_array()` and `registry.stable_order_array()` give (N, 3) numpy arrays for export and analysis.

If the user's "second structure" was intended to be a separate matrix (not orderings), the most natural candidate is:

```python
# Rank matrix: R[i, d-1] = UNSTABLE RANK of the d-th forward iterate of ids[i]
# tells you: after d steps, where does this intersection land in the W^u ordering?
def rank_matrix(registry, max_depth):
    ids = registry.all_ids()
    R = np.full((len(ids), max_depth), fill_value=-1, dtype=np.int64)
    for i, iid in enumerate(ids):
        for d in range(1, max_depth + 1):
            target = registry.iterate_table[iid, d]
            if target is not None:
                R[i, d - 1] = registry.unstable_rank(target)
    return R
```

This is the "symbolic dynamics matrix" — it tells you not just which intersection
maps to which (the iterate table), but what symbolic position (W^u rank) the image
occupies. This is the building block for topological entropy calculations. Whether to
build it as a method on `IntersectionRegistry` or as a standalone utility depends on
how much topological analysis will be done.

---

## Future Structures (not in this implementation phase)

### Lobe Table

A lobe is the region bounded by two consecutive intersections on W^u and the
corresponding arcs of W^s between them. There are N-1 lobes for N intersections (plus
two half-infinite lobes at the ends). The lobe table maps:
- Lobe index → the two intersection IDs that bound it
- Lobe index → its pre-image lobe (under f^{-1})
- Lobe index → area (computed from the manifold arclength integrals)

### Markov / Transition Matrix

An (N-1) × (N-1) matrix where `M[i, j] = 1` if lobe i maps (under one iterate) to
overlap with lobe j. This is the topological Markov chain. Its largest eigenvalue's
log is a lower bound on topological entropy.

### Symbol Sequences

For each intersection ID, its full forward and backward symbol sequences:
```python
sequence = [registry.unstable_rank(id)
            for id in registry.iterate_table.forward_chain(some_id)]
```
These are the "addresses" of intersection points in the symbolic dynamics.

### The Full Tangle Template

A graph where nodes are intersections (sorted by unstable cdist) and edges connect
adjacent intersections in the unstable ordering. The iterate map on this graph gives
the full topological template (the "tangle matrix" in MacKay–Meiss–Percival theory).

---

## Implementation Order

Follow this order. Each step only depends on what came before.

| Step | File | Type | Depends on |
|------|------|------|-----------|
| 1 | `Intersection.py` | Modify | — |
| 2 | `IterateTable.py` | New | — |
| 3 | `IntersectionRegistry.py` | New | `Intersection`, `IterateTable` |
| 4 | `TangleWorkbench.py` | Modify | `IntersectionRegistry`, bridge plan changes |
| 5 | `__init__.py` | Modify | all above |

The bridge plan changes from `bridge_iteration_refactor_plan.md` should be implemented
**before** step 4 here, since `TangleWorkbench.iterate_bridge()` is what generates the
`Intersection` objects that feed the registry.

---

## Summary of New Files

| File | Purpose |
|------|---------|
| `IterateTable.py` | `table[id, n]` → id; auto-reverse; forward/backward chain; dense array export |
| `IntersectionRegistry.py` | Master store; unique IDs; collision detect; both orderings; `infer_iterates()` |

## Summary of Modified Files

| File | Changes |
|------|---------|
| `Intersection.py` | Add `id: Optional[int] = None` field |
| `TangleWorkbench.py` | Add `_intersection_registry`; update `compute_intersections()` and `iterate_bridge()` to populate it; add `infer_iterate_table()` and `populate_registry()` |
| `__init__.py` | Export `IterateTable` and `IntersectionRegistry` |
