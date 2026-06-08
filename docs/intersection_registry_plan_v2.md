# Intersection Registry and Topological Layer — Implementation Plan (v2)

> **What changed from v1:**  
> This revision builds on the *existing* `Intersection` class (a plain class, not a dataclass).
> It eliminates stored "fictitious" intersections in favour of live query methods.  
> It generalises `Tangle` to detect crossings between *any* manifold pair and records
> manifold provenance on every intersection.  
> It adds a live NetworkX graph maintained automatically inside `IntersectionRegistry`.  
> The filter interface is designed so adding a new query method is a one-liner.

---

## Overview and Goals

After generating enough intersections via bridge iteration, the workflow transitions from
manifold numerics to topology — the combinatorial structure of how intersection points
order on each manifold, how they permute under the iterate map, and what that implies
about transport and mixing.

The topological layer must be:

- **Separable** — once the registry is populated, `TangleWorkbench` and `Tangle` are no
  longer needed. All topology lives in the registry and its graph.
- **Connected** — the registry is populated automatically as bridges are iterated.
- **Provenance-aware** — every intersection records which manifolds (fixed point, stability,
  orbit index, branch index) were involved in the crossing.
- **General** — the same infrastructure handles homoclinic intersections (W^u ∩ W^s of
  one fixed point), heteroclinic crossings (between two different fixed points), and any
  other manifold pair.
- **Queryable** — `registry[id]` → `Intersection`; `registry.on_interval(lo, hi)` →
  source intersections whose forward iterates land in [lo, hi]; sorted orderings on W^u
  and W^s immediately accessible.
- **Graphed** — a live `nx.MultiDiGraph` is maintained inside the registry. New
  intersections are added as nodes automatically on registration; iterate edges are wired
  as the iterate table is populated.
- **Extensible** — adding a new query method is a 5-line function delegating to
  `_filter_intersections()`.

---

## Mathematical Background

### What intersections are

Every point p in W^u ∩ W^s maps to another intersection under f:
- f(p) is also in W^u ∩ W^s (f preserves both manifolds)
- The iterate map on intersections is well-defined and discrete

### How cdists transform under iteration

For a period-1 saddle fixed point with unstable eigenvalue λ_u > 1:
```
unstable_cdist( f^n(p) ) = λ_u^n  × unstable_cdist(p)     # grows under forward map
stable_cdist(   f^n(p) ) = λ_u^{-n} × stable_cdist(p)     # shrinks under forward map
```

This is the **key formula**. Given any intersection p with known cdists, we can predict
exactly where f^n(p) lives in (unstable_cdist, stable_cdist) space without any further
manifold numerics. In the heteroclinic case (p on W^u of FP_1 and W^s of FP_2), λ_u is
taken from the manifold that owns the cdist being scaled — each side uses its own fixed
point's eigenvalue.

### The two orderings

Intersections live on both manifolds simultaneously. They can be sorted two ways:
- **Unstable order** — sorted by `unstable_cdist`: the sequence of crossings encountered
  walking along W^u away from the fixed point.
- **Stable order** — sorted by `stable_cdist`: the sequence of crossings encountered
  walking along W^s away from the fixed point.

The permutation taking unstable order to stable order encodes the topological type of the
tangle and drives topological-entropy calculations.

### The iterate table

If intersections are labelled by their position in unstable order, the iterate map permutes
these labels. `IterateTable[id, n]` records which intersection ID maps to which under f^n.

---

## Architecture

```
┌──────────────────────────────────────────────────────────┐
│  NUMERIC LAYER                                           │
│  TangleWorkbench, ManifoldMachine, Tangle                │
│  - Grows manifolds as linked lists                       │
│  - Detects geometric crossings between ANY pair of       │
│    manifolds (R-tree + segment intersection tests)       │
│  - Cuts bridges                                          │
│  - Produces Intersection objects with coords, cdists,    │
│    and ManifoldKey provenance for each side              │
└───────────────────────┬──────────────────────────────────┘
                        │  populate / auto-populate
                        │
┌───────────────────────▼──────────────────────────────────┐
│  TOPOLOGICAL LAYER                                       │
│  IntersectionRegistry, IterateTable                      │
│  - Unique IDs; collision detection                       │
│  - ManifoldKey provenance per intersection               │
│  - IterateTable: registry[id, n] → id                    │
│  - Sorted orderings on W^u and W^s                       │
│  - Live nx.MultiDiGraph (nodes = IDs, edges = adjacency  │
│    on manifold + iterate map)                            │
│  - Query interface: on_interval(), on_cdist_range(),     │
│    from_fixed_point(), from_branch(), filter()           │
│  - No manifold objects needed from here on               │
└──────────────────────────────────────────────────────────┘
```

The separation point is `IntersectionRegistry`. Everything below is geometry; everything
in it and above is topology.

---

## Type Definition: `ManifoldKey`

`ManifoldKey` is a 4-tuple that uniquely identifies a manifold within `TangleWorkbench`.
It uses **exactly the same format** as the keys of `TangleWorkbench.manifolds`, so a
key retrieved from an intersection can look up the manifold object directly:

```python
manifold = wb.manifolds[intersection.manifold_a_key]
```

### Definition

```python
# In src/tanglepack/Intersection.py
from typing import Literal, TYPE_CHECKING

if TYPE_CHECKING:
    from .FixedPoint import FixedPoint

# ManifoldKey = (fixed_point, stability, orbit_index, branch_index)
# Identical to the key type used in TangleWorkbench.manifolds.
ManifoldKey = tuple["FixedPoint", Literal["unstable", "stable"], int, int]
```

- `[0]` `FixedPoint` — the anchor point (saddle) the manifold emanates from
- `[1]` `stability` — `"unstable"` or `"stable"`
- `[2]` `orbit_index` — index into the periodic orbit (0 for period-1 saddles)
- `[3]` `branch_index` — 0 or 1 for fixed points with inversion (negative eigenvalue),
  always 0 otherwise

Import `ManifoldKey` from `Intersection` anywhere it is needed:
```python
from .Intersection import ManifoldKey
```

---

## Updated: `src/tanglepack/Intersection.py`

### What changes

1. Add `ManifoldKey` TypeAlias (see above).
2. Add `label: Optional[str]` to `__init__` — this field was accepted by the class methods
   but silently dropped because `__init__` had no `label` parameter.
3. Add `manifold_a_key: Optional[ManifoldKey]` and `manifold_b_key: Optional[ManifoldKey]`.
   **Convention**: `manifold_a_key` owns `unstable_cdist`; `manifold_b_key` owns
   `stable_cdist`. For a standard homoclinic intersection, `manifold_a_key` will point to
   the unstable manifold and `manifold_b_key` to the stable manifold of the same fixed
   point. For a heteroclinic intersection, they will point to different fixed points.
4. Update `from_segments()` to accept and forward the new keys.
5. Add a `fixed_points` property for quick access to the involved fixed points.
6. Keep `is_synthetic` and `get_point()` unchanged.
7. The class remains a **plain class** (no `@dataclass`).

### Updated `__init__`

```python
def __init__(
    self,
    coords: tuple[float, float] = None,
    unstable_cdist: float = None,
    stable_cdist: float = None,
    seg_ids: Optional[frozenset[int]] = None,
    id: Optional[int] = None,
    label: Optional[str] = None,
    manifold_a_key: Optional[ManifoldKey] = None,
    manifold_b_key: Optional[ManifoldKey] = None,
):
    self.coords = coords
    self.unstable_cdist = unstable_cdist
    self.stable_cdist = stable_cdist
    self.seg_ids = seg_ids
    self.id = id
    self.label = label
    self.manifold_a_key = manifold_a_key
    self.manifold_b_key = manifold_b_key
```

### Updated `from_segments` classmethod

```python
@classmethod
def from_segments(
    cls,
    coords: tuple[float, float],
    unstable_cdist: float,
    stable_cdist: float,
    seg1_id: int,
    seg2_id: int,
    manifold_a_key: Optional[ManifoldKey] = None,
    manifold_b_key: Optional[ManifoldKey] = None,
    label: Optional[str] = None,
) -> Intersection:
    """Create an Intersection backed by two R-tree segment IDs."""
    return cls(
        coords=coords,
        unstable_cdist=unstable_cdist,
        stable_cdist=stable_cdist,
        seg_ids=frozenset({seg1_id, seg2_id}),
        label=label,
        manifold_a_key=manifold_a_key,
        manifold_b_key=manifold_b_key,
    )
```

### New `fixed_points` property

```python
@property
def fixed_points(self) -> tuple:
    """Return the distinct FixedPoint objects involved in this intersection."""
    fps = []
    if self.manifold_a_key is not None:
        fps.append(self.manifold_a_key[0])
    if self.manifold_b_key is not None and self.manifold_b_key[0] is not fps[0]:
        fps.append(self.manifold_b_key[0])
    return tuple(fps)
```

### Note on non-standard manifold pairs

When both manifolds have the same stability (e.g., two unstable branches from different
fixed points), the field names `unstable_cdist` and `stable_cdist` are misnomers — the
manifold keys clarify the true meaning. A future refactor could rename them to
`cdist_a`/`cdist_b`, but the current names are kept for backward compatibility.

---

## Updated: `src/tanglepack/BaseManifold.py`

Add `manifold_key` to `__init__`. `TangleWorkbench` sets this immediately after creating
each manifold so that `Tangle` can read it when building intersection objects.

### Change: add `manifold_key` parameter

```python
def __init__(
    self,
    root: Point | BranchPoint,
    stability: Literal["stable", "unstable"],
    stretch_param: float,
    fixed_point: FixedPoint,
    name: str = "unnamed",
    tail: Optional[Point | BranchPoint] = None,
    branch_index: Optional[int] = None,
    manifold_key: Optional[ManifoldKey] = None,   # ← new
):
    ...
    self.manifold_key = manifold_key
```

`manifold_key` defaults to `None` and is filled in by `TangleWorkbench.initialize_manifold()`
right after the manifold is constructed. `Bridge` inherits `BaseManifold`, so iterated
bridges carry their key automatically (via the Bridge constructor calling `super().__init__`
with the parent manifold's key).

---

## Updated: `src/tanglepack/Tangle.py`

### Goal

Generalise intersection detection so it works between **any two manifold objects**,
including:
- W^u(FP_1) ∩ W^s(FP_1) — standard homoclinic
- W^u(FP_1) ∩ W^s(FP_2) — heteroclinic
- W^u(FP_1) ∩ W^u(FP_2) — same-stability crossing between different fixed points

The geometric detection machinery (R-tree, segment intersection test) already handles all
these cases — `_insert_segment` only skips segments from the exact same Python manifold
object. The only change needed is in how `Intersection` objects are built: manifold keys
must be extracted from each segment's manifold and attached to the intersection.

### Change: helper to extract manifold key from a segment

Add a private helper to `Tangle`:

```python
@staticmethod
def _key_of(seg: _Segment) -> Optional[ManifoldKey]:
    """Read the ManifoldKey stored on the segment's manifold, if set."""
    return getattr(seg.manifold, "manifold_key", None)
```

### Change: generalised intersection building in `populate_intersection_dict`

Replace the hardcoded `unstable`/`stable` cdist assignment with the generalised form:

```python
# Determine which segment's manifold "owns" each cdist field.
# Convention: the unstable-stability manifold becomes manifold_a;
# if neither (or both) are unstable, seg_1 becomes manifold_a by default.
if seg_1.manifold.stability == "unstable" or seg_2.manifold.stability != "unstable":
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
```

Apply **the same change** in `populate_intersections_for_manifold` — same logic, same
placement after `point = self._find_true_intersection(...)`.

No other changes to `Tangle`. The R-tree and segment-intersection machinery is already
general.

---

## New File: `src/tanglepack/IterateTable.py`

Unchanged from v1. Full implementation reproduced below for completeness.

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

    Setting one direction auto-records the reverse:
        table[3, 2] = 7  also records  table[7, -2] = 3

    Attributes:
        _forward: dict[int, dict[int, int]]  — _forward[id][n] = target_id  (n > 0)
        _backward: dict[int, dict[int, int]] — _backward[id][n] = target_id (n > 0)
    """

    def __init__(self):
        self._forward: dict[int, dict[int, int]] = {}
        self._backward: dict[int, dict[int, int]] = {}

    def __getitem__(self, key: tuple[int, int]) -> Optional[int]:
        source_id, n = key
        if n == 0:
            return source_id
        if n > 0:
            return self._forward.get(source_id, {}).get(n)
        else:
            return self._backward.get(source_id, {}).get(-n)

    def __setitem__(self, key: tuple[int, int], target_id: int):
        source_id, n = key
        if n == 0:
            return
        if n > 0:
            self._forward.setdefault(source_id, {})[n] = target_id
            self._backward.setdefault(target_id, {})[n] = source_id
        else:
            self._backward.setdefault(source_id, {})[-n] = target_id
            self._forward.setdefault(target_id, {})[-n] = source_id

    def __contains__(self, key: tuple[int, int]) -> bool:
        return self[key] is not None

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

    def as_forward_array(self, ids: list[int], max_depth: int) -> NDArray[np.int64]:
        """
        Dense 2D array A where A[i, d-1] = ID of f^d(ids[i]), or -1 if unknown.

        Args:
            ids: Ordered list of intersection IDs (defines row order).
            max_depth: Number of forward iterate columns.

        Returns:
            Array of shape (len(ids), max_depth), dtype int64.
        """
        arr = np.full((len(ids), max_depth), fill_value=-1, dtype=np.int64)
        for i, iid in enumerate(ids):
            for d in range(1, max_depth + 1):
                target = self[iid, d]
                if target is not None:
                    arr[i, d - 1] = target
        return arr

    def as_backward_array(self, ids: list[int], max_depth: int) -> NDArray[np.int64]:
        """Dense array B[i, d-1] = ID of f^{-d}(ids[i]), or -1 if unknown."""
        arr = np.full((len(ids), max_depth), fill_value=-1, dtype=np.int64)
        for i, iid in enumerate(ids):
            for d in range(1, max_depth + 1):
                target = self[iid, -d]
                if target is not None:
                    arr[i, d - 1] = target
        return arr

    def register_iterate(self, source_id: int, n: int, target_id: int):
        """Explicit named method — delegates to __setitem__."""
        self[source_id, n] = target_id
```

---

## New File: `src/tanglepack/IntersectionRegistry.py`

This is the most heavily revised component relative to v1.

### Key additions relative to v1

1. **`_graph: nx.MultiDiGraph`** — maintained live; nodes added on `add()`, iterate
   edges added on `register_iterate()`, adjacency edges rebuilt lazily.
2. **`_filter_intersections(predicate)`** — the extensibility core. Every named query
   method is a thin wrapper around this.
3. **`on_interval(lo, hi, stability, fixed_point, branch_index)`** — returns intersections
   whose forward iterate's cdist falls in [lo, hi] on the given manifold. Handles
   multiple fixed points and branches transparently.
4. **`on_cdist_range(lo, hi, stability)`** — direct cdist filter (no iteration applied).
5. **`from_fixed_point(fp)`** and **`from_branch(branch_index)`** — provenance filters.
6. **`infer_iterates()`** — no longer takes `lambda_u` as a parameter; reads it per
   intersection from `manifold_a_key`. Handles multi-FP registries automatically.
7. **`add_synthetic()`** — retained as a low-level utility but no longer the primary
   workflow. Prefer computing real intersections and querying them.

### Full implementation

```python
from __future__ import annotations

import bisect
from typing import Callable, Literal, Optional, TYPE_CHECKING

import networkx as nx
import numpy as np
from numpy.typing import NDArray

from .Intersection import Intersection, ManifoldKey
from .IterateTable import IterateTable

if TYPE_CHECKING:
    from .FixedPoint import FixedPoint


class IntersectionRegistry:
    """
    Master store of all intersection points in the tangle.

    IDs are contiguous integers (0, 1, 2, …) assigned in insertion order.
    Duplicate intersections (same cdists within tolerance) are detected and
    deduped — the existing ID is returned rather than creating a second entry.

    Primary interface:
        registry.add(intersection)               → int (assigned ID)
        registry[id]                             → Intersection
        registry.iterate_table[id, n]            → int or None
        registry.by_unstable_cdist               → list[int]  (sorted by u-cdist)
        registry.by_stable_cdist                 → list[int]  (sorted by s-cdist)
        registry.graph                           → nx.MultiDiGraph (live)

    Query interface (all return list[Intersection]):
        registry.on_interval(lo, hi)             → pre-images that map into [lo, hi]
        registry.on_cdist_range(lo, hi)          → intersections with cdist in [lo, hi]
        registry.from_fixed_point(fp)            → intersections involving fp
        registry.from_branch(branch_index)       → intersections on given branch
        registry.filter(predicate)               → arbitrary predicate

    Attributes:
        _store: dict[int, Intersection]
        _next_id: int
        cdist_tol: float
        iterate_table: IterateTable
        _unstable_order: list[int]  — IDs sorted ascending by unstable_cdist
        _stable_order: list[int]    — IDs sorted ascending by stable_cdist
        _cdist_index: dict[tuple[float, float], int]  — secondary collision index
        _graph: nx.MultiDiGraph
        _graph_adjacency_dirty: bool
    """

    def __init__(self, cdist_tol: float = 1e-6):
        self._store: dict[int, Intersection] = {}
        self._next_id: int = 0
        self.cdist_tol = cdist_tol
        self.iterate_table = IterateTable()

        self._unstable_order: list[int] = []
        self._stable_order: list[int] = []
        self._cdist_index: dict[tuple[float, float], int] = {}

        self._graph: nx.MultiDiGraph = nx.MultiDiGraph()
        self._graph_adjacency_dirty: bool = False

    # ── core insert / lookup ───────────────────────────────────────────────

    def add(self, intersection: Intersection) -> int:
        """
        Register an intersection, returning its unique ID.

        If a collision is detected (another intersection with cdists within
        self.cdist_tol), the existing ID is returned and no duplicate is stored.
        On a new insertion, the node is added to self.graph immediately.

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

        key = self._cdist_key(intersection)
        self._cdist_index[key] = new_id

        # Add node to the live graph.
        self._graph.add_node(
            new_id,
            coords=intersection.coords,
            unstable_cdist=intersection.unstable_cdist,
            stable_cdist=intersection.stable_cdist,
            manifold_a_key=intersection.manifold_a_key,
            manifold_b_key=intersection.manifold_b_key,
            label=intersection.label,
        )
        self._graph_adjacency_dirty = True

        return new_id

    def add_synthetic(
        self,
        coords: tuple[float, float],
        unstable_cdist: float,
        stable_cdist: float,
        label: Optional[str] = None,
        manifold_a_key: Optional[ManifoldKey] = None,
        manifold_b_key: Optional[ManifoldKey] = None,
    ) -> int:
        """
        Low-level utility: add an intersection not backed by a detected segment
        crossing (e.g., a manually placed reference point).

        For most workflows, prefer computing real intersections and querying them
        with on_interval() and related methods.
        """
        return self.add(
            Intersection(
                coords=coords,
                unstable_cdist=unstable_cdist,
                stable_cdist=stable_cdist,
                label=label,
                manifold_a_key=manifold_a_key,
                manifold_b_key=manifold_b_key,
            )
        )

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

    # ── ordering views ─────────────────────────────────────────────────────

    @property
    def by_unstable_cdist(self) -> list[int]:
        """IDs sorted ascending by unstable_cdist."""
        return list(self._unstable_order)

    @property
    def by_stable_cdist(self) -> list[int]:
        """IDs sorted ascending by stable_cdist."""
        return list(self._stable_order)

    def unstable_rank(self, id: int) -> int:
        """0-based position of intersection `id` in the W^u ordering."""
        return self._unstable_order.index(id)

    def stable_rank(self, id: int) -> int:
        """0-based position of intersection `id` in the W^s ordering."""
        return self._stable_order.index(id)

    def all_ids(self) -> list[int]:
        """All registered IDs in insertion order."""
        return list(self._store.keys())

    # ── live graph ─────────────────────────────────────────────────────────

    @property
    def graph(self) -> nx.MultiDiGraph:
        """
        The live intersection graph.

        Nodes are intersection IDs. Two edge types are maintained:
          - type="adjacency", stability="unstable"/"stable":
              connects consecutive intersections in the respective sorted ordering.
              These are rebuilt lazily whenever the sorted order changes.
          - type="iterate", n=<int>:
              directed edge (p → f^n(p)) added by register_iterate().

        The graph is always up to date w.r.t. nodes and iterate edges.
        Adjacency edges are rebuilt on first access after any add().
        """
        if self._graph_adjacency_dirty:
            self._rebuild_adjacency_edges()
        return self._graph

    def _rebuild_adjacency_edges(self):
        """Remove and rebuild all adjacency-type edges from the current sorted lists."""
        stale = [
            (u, v, k)
            for u, v, k, d in self._graph.edges(keys=True, data=True)
            if d.get("type") == "adjacency"
        ]
        for u, v, k in stale:
            self._graph.remove_edge(u, v, k)

        for i in range(len(self._unstable_order) - 1):
            u, v = self._unstable_order[i], self._unstable_order[i + 1]
            self._graph.add_edge(
                u, v,
                key=f"adj_u_{i}",
                type="adjacency",
                stability="unstable",
            )

        for i in range(len(self._stable_order) - 1):
            u, v = self._stable_order[i], self._stable_order[i + 1]
            self._graph.add_edge(
                u, v,
                key=f"adj_s_{i}",
                type="adjacency",
                stability="stable",
            )

        self._graph_adjacency_dirty = False

    # ── iterate table ──────────────────────────────────────────────────────

    def register_iterate(self, source_id: int, n: int, target_id: int):
        """
        Record f^n(source) = target.

        Delegates to iterate_table and wires the directed edge into the graph.
        """
        self.iterate_table.register_iterate(source_id, n, target_id)
        if source_id in self._graph and target_id in self._graph:
            self._graph.add_edge(
                source_id, target_id,
                key=f"iter_{source_id}_{n}",
                type="iterate",
                n=n,
            )

    def infer_iterates(
        self,
        max_depth: int = 10,
        tol_multiplier: float = 10.0,
    ) -> int:
        """
        Scan all registered intersections and infer iterate relationships by
        predicting cdists and checking for matches.

        Reads lambda_u per intersection from its manifold_a_key, so this works
        correctly with multiple fixed points in the same registry.

        Unlike v1, this method takes no lambda_u parameter.

        Args:
            max_depth: How many iterate levels to search forward and backward.
            tol_multiplier: Scale factor on cdist_tol for matching.

        Returns:
            Number of new iterate relationships recorded.
        """
        tol = self.cdist_tol * tol_multiplier
        recorded = 0

        for source_id, source in self._store.items():
            lambda_u = self._get_lambda_u(source)
            if lambda_u is None:
                continue

            for n in range(1, max_depth + 1):
                if (source_id, n) not in self.iterate_table:
                    pred_u = (lambda_u ** n) * source.unstable_cdist
                    pred_s = source.stable_cdist / (lambda_u ** n)
                    target_id = self.find_by_cdist(pred_u, pred_s, tol)
                    if target_id is not None:
                        self.register_iterate(source_id, n, target_id)
                        recorded += 1

                if (source_id, -n) not in self.iterate_table:
                    pred_u = source.unstable_cdist / (lambda_u ** n)
                    pred_s = source.stable_cdist * (lambda_u ** n)
                    target_id = self.find_by_cdist(pred_u, pred_s, tol)
                    if target_id is not None:
                        self.register_iterate(source_id, -n, target_id)
                        recorded += 1

        return recorded

    # ── query interface ────────────────────────────────────────────────────

    def filter(self, predicate: Callable[[Intersection], bool]) -> list[Intersection]:
        """
        Return all registered intersections that satisfy predicate.

        This is the extensibility core. All named query methods below are thin
        wrappers that build a predicate and delegate here.

        Args:
            predicate: A callable that takes an Intersection and returns bool.

        Returns:
            List of matching Intersection objects (no guaranteed ordering).

        Example:
            registry.filter(lambda ix: ix.unstable_cdist > 10.0)
        """
        return [ix for ix in self._store.values() if predicate(ix)]

    def on_interval(
        self,
        lo: float,
        hi: float,
        stability: Literal["unstable", "stable"] = "unstable",
        fixed_point: Optional[FixedPoint] = None,
        branch_index: Optional[int] = None,
    ) -> list[Intersection]:
        """
        Return all intersections p such that f(p) has cdist on `stability` in [lo, hi].

        Uses the eigenvalue formula — no additional manifold numerics required:
            f(p).unstable_cdist = lambda_u × p.unstable_cdist
            f(p).stable_cdist   = p.stable_cdist / lambda_u

        lambda_u is read per intersection from manifold_a_key (for stability="unstable")
        or manifold_b_key (for stability="stable"), so a registry containing intersections
        from multiple fixed points with different eigenvalues is handled correctly.

        Args:
            lo: Lower cdist bound for f(p) on the given manifold.
            hi: Upper cdist bound for f(p) on the given manifold.
            stability: Which manifold's cdist to apply the interval on.
            fixed_point: If given, only consider intersections involving this FP.
            branch_index: If given, only consider intersections on this branch.

        Returns:
            List of Intersection objects (the sources, not their images).
        """
        results = []
        for ix in self._store.values():
            # Optional provenance filters
            if fixed_point is not None or branch_index is not None:
                key = ix.manifold_a_key if stability == "unstable" else ix.manifold_b_key
                if key is None:
                    continue
                if fixed_point is not None and key[0] is not fixed_point:
                    continue
                if branch_index is not None and key[3] != branch_index:
                    continue

            lambda_u = self._get_lambda_u(ix)
            if lambda_u is None:
                continue

            if stability == "unstable":
                image_cdist = lambda_u * ix.unstable_cdist
            else:
                image_cdist = ix.stable_cdist / lambda_u

            if lo <= image_cdist <= hi:
                results.append(ix)

        return results

    def on_cdist_range(
        self,
        lo: float,
        hi: float,
        stability: Literal["unstable", "stable"] = "unstable",
    ) -> list[Intersection]:
        """
        Return all intersections whose CURRENT cdist on `stability` is in [lo, hi].

        Unlike on_interval(), no iteration is applied — this filters by where the
        intersection already sits, not where it maps to.

        Args:
            lo: Lower bound of the cdist range.
            hi: Upper bound of the cdist range.
            stability: Which manifold cdist to filter on.

        Returns:
            List of matching Intersection objects sorted by the chosen cdist.
        """
        attr = "unstable_cdist" if stability == "unstable" else "stable_cdist"
        return self.filter(lambda ix: lo <= getattr(ix, attr) <= hi)

    def from_fixed_point(
        self,
        fp: FixedPoint,
        stability: Optional[Literal["unstable", "stable"]] = None,
    ) -> list[Intersection]:
        """
        Return all intersections that involve the given fixed point.

        Checks both manifold_a_key and manifold_b_key. Optionally restrict to
        intersections where fp appears on the specified stability side.

        Args:
            fp: The FixedPoint to filter by.
            stability: If given, only match on that side (manifold_a for "unstable",
                manifold_b for "stable").

        Returns:
            List of matching Intersection objects.
        """
        def pred(ix: Intersection) -> bool:
            a_match = ix.manifold_a_key is not None and ix.manifold_a_key[0] is fp
            b_match = ix.manifold_b_key is not None and ix.manifold_b_key[0] is fp
            if stability == "unstable":
                return a_match
            if stability == "stable":
                return b_match
            return a_match or b_match

        return self.filter(pred)

    def from_branch(
        self,
        branch_index: int,
        stability: Optional[Literal["unstable", "stable"]] = None,
    ) -> list[Intersection]:
        """
        Return all intersections on the given branch.

        Args:
            branch_index: 0 or 1.
            stability: If given, only check the corresponding manifold key side.

        Returns:
            List of matching Intersection objects.
        """
        def pred(ix: Intersection) -> bool:
            a_match = ix.manifold_a_key is not None and ix.manifold_a_key[3] == branch_index
            b_match = ix.manifold_b_key is not None and ix.manifold_b_key[3] == branch_index
            if stability == "unstable":
                return a_match
            if stability == "stable":
                return b_match
            return a_match or b_match

        return self.filter(pred)

    # ── cdist-based lookup ─────────────────────────────────────────────────

    def find_by_cdist(
        self,
        unstable_cdist: float,
        stable_cdist: float,
        tol: Optional[float] = None,
    ) -> Optional[int]:
        """
        Find the ID of an intersection whose cdists are within tol of the given values.

        Args:
            unstable_cdist: Target unstable cdist.
            stable_cdist: Target stable cdist.
            tol: Search tolerance. Defaults to self.cdist_tol.

        Returns:
            Matching ID, or None.
        """
        if tol is None:
            tol = self.cdist_tol
        for id, existing in self._store.items():
            if (abs(existing.unstable_cdist - unstable_cdist) < tol
                    and abs(existing.stable_cdist - stable_cdist) < tol):
                return id
        return None

    # ── array exports ──────────────────────────────────────────────────────

    def as_forward_array(self, max_depth: int = 5) -> NDArray[np.int64]:
        """
        Dense array A[i, d-1] = ID of f^d(ids[i]) where ids = all_ids() in order.
        Shape: (N, max_depth). -1 = unknown.
        """
        return self.iterate_table.as_forward_array(self.all_ids(), max_depth)

    def as_backward_array(self, max_depth: int = 5) -> NDArray[np.int64]:
        """Dense array B[i, d-1] = ID of f^{-d}(ids[i]). Shape (N, max_depth). -1 = unknown."""
        return self.iterate_table.as_backward_array(self.all_ids(), max_depth)

    def unstable_order_array(self) -> NDArray[np.float64]:
        """(N, 3) array [id, unstable_cdist, stable_cdist] sorted by unstable_cdist."""
        rows = [[iid, self._store[iid].unstable_cdist, self._store[iid].stable_cdist]
                for iid in self._unstable_order]
        return np.array(rows, dtype=np.float64)

    def stable_order_array(self) -> NDArray[np.float64]:
        """(N, 3) array [id, unstable_cdist, stable_cdist] sorted by stable_cdist."""
        rows = [[iid, self._store[iid].unstable_cdist, self._store[iid].stable_cdist]
                for iid in self._stable_order]
        return np.array(rows, dtype=np.float64)

    # ── internal helpers ───────────────────────────────────────────────────

    def _get_lambda_u(self, intersection: Intersection) -> Optional[float]:
        """
        Read the unstable eigenvalue magnitude from the intersection's manifold keys.

        For stability="unstable", uses manifold_a_key (the unstable manifold side).
        For stability="stable", uses manifold_b_key.
        If neither key is set, returns None and the intersection is skipped by callers.
        """
        key = intersection.manifold_a_key or intersection.manifold_b_key
        if key is None:
            return None
        fp = key[0]
        if not hasattr(fp, "unstable_eigenvalues") or not fp.unstable_eigenvalues:
            return None
        return abs(fp.unstable_eigenvalues[0])

    def _find_collision(self, intersection: Intersection) -> Optional[int]:
        """Linear scan for an existing intersection within cdist_tol."""
        for id, existing in self._store.items():
            if (abs(existing.unstable_cdist - intersection.unstable_cdist) < self.cdist_tol
                    and abs(existing.stable_cdist - intersection.stable_cdist) < self.cdist_tol):
                return id
        return None

    def _cdist_key(self, intersection: Intersection) -> tuple[float, float]:
        digits = max(0, -int(np.floor(np.log10(self.cdist_tol))) - 1)
        return (
            round(intersection.unstable_cdist, digits),
            round(intersection.stable_cdist, digits),
        )

    def _insert_into_unstable_order(self, id: int, unstable_cdist: float):
        keys = [self._store[i].unstable_cdist for i in self._unstable_order]
        pos = bisect.bisect_left(keys, unstable_cdist)
        self._unstable_order.insert(pos, id)

    def _insert_into_stable_order(self, id: int, stable_cdist: float):
        keys = [self._store[i].stable_cdist for i in self._stable_order]
        pos = bisect.bisect_left(keys, stable_cdist)
        self._stable_order.insert(pos, id)
```

---

## Updated: `src/tanglepack/TangleWorkbench.py`

### New attribute in `__init__`

```python
self._intersection_registry = IntersectionRegistry()
self._bridges: list[Bridge] = []     # was None; initialise as empty list
```

### New property

```python
@property
def intersection_registry(self) -> IntersectionRegistry:
    return self._intersection_registry
```

### Updated `initialize_manifold()` — set `manifold_key` on each new manifold

```python
for (orbit_index, branch_index), manifold in initial_segments.items():
    key = (fixed_point, stability, orbit_index, branch_index)
    manifold.manifold_key = key                              # ← new
    self.manifolds[key] = manifold
```

This ensures every `BaseManifold` object carries a `manifold_key` that `Tangle` can read.

### Updated `compute_intersections()`

```python
def compute_intersections(self, fixed_point: FixedPoint):
    self.Tangle.clear_all()
    self._intersection_registry = IntersectionRegistry()    # fresh on rebuild

    self.index_manifolds(fixed_point, "unstable")
    self.index_manifolds(fixed_point, "stable")
    self.Tangle.populate_intersection_dict()

    for intersection in self.Tangle._intersections:
        self._intersection_registry.add(intersection)

    return list(self.Tangle._intersecting_coords.values())
```

### Updated `iterate_bridge()`

After step 3 (resolve new crossings), register them:

```python
new_intersections = self.Tangle.populate_intersections_for_manifold(iterated)

for ix in new_intersections:
    self._intersection_registry.add(ix)
```

### New method: `infer_iterate_table()`

No longer needs `fixed_point` or `lambda_u` — the registry reads eigenvalues from each
intersection's stored manifold key:

```python
def infer_iterate_table(self, max_depth: int = 10) -> int:
    """
    Fill in the registry's iterate table by predicting cdist values of forward
    and backward iterates and checking for matches.

    Args:
        max_depth: How many iterate levels to search forward and backward.

    Returns:
        Number of new iterate relationships added.
    """
    return self._intersection_registry.infer_iterates(max_depth)
```

### New method: `populate_registry()`

```python
def populate_registry(self) -> IntersectionRegistry:
    """Rebuild the intersection registry from the current Tangle state."""
    self._intersection_registry = IntersectionRegistry()
    for intersection in self.Tangle._intersections:
        self._intersection_registry.add(intersection)
    return self._intersection_registry
```

### New method: `build_intersection_graph()`

```python
def build_intersection_graph(self) -> nx.MultiDiGraph:
    """
    Return the live intersection graph from the registry.

    Accessing this property triggers a rebuild of adjacency edges if any
    intersections have been added since the last access.
    """
    return self._intersection_registry.graph
```

---

## Updated: `src/tanglepack/__init__.py`

```python
from .IterateTable import IterateTable
from .IntersectionRegistry import IntersectionRegistry
from .Intersection import Intersection, ManifoldKey
```

---

## Extensibility Guide: Adding New Query Methods

All query methods follow the same pattern: build a predicate, delegate to `filter()`.
Adding a new one is typically 5-10 lines. Examples:

```python
# Filter intersections near a physical coordinate (within radius r)
def near_coords(self, x: float, y: float, radius: float) -> list[Intersection]:
    return self.filter(
        lambda ix: np.hypot(ix.coords[0] - x, ix.coords[1] - y) <= radius
    )

# Filter intersections whose second forward iterate is registered
def with_second_iterate(self) -> list[Intersection]:
    return self.filter(lambda ix: (ix.id, 2) in self.iterate_table)

# Filter heteroclinic intersections only
def heteroclinic(self) -> list[Intersection]:
    def pred(ix: Intersection) -> bool:
        if ix.manifold_a_key is None or ix.manifold_b_key is None:
            return False
        return ix.manifold_a_key[0] is not ix.manifold_b_key[0]
    return self.filter(pred)

# Filter by minimum cdist (intersections deep enough into the tangle)
def beyond_cdist(self, threshold: float, stability: str = "unstable") -> list[Intersection]:
    attr = "unstable_cdist" if stability == "unstable" else "stable_cdist"
    return self.filter(lambda ix: getattr(ix, attr) >= threshold)
```

Rules for new methods:
1. Name it by what it returns, not how it works.
2. Use `self.filter(predicate)` as the implementation body.
3. Accept `Optional[FixedPoint]` and `Optional[int]` (for `branch_index`) when the filter
   is parametrised by manifold provenance, so callers can narrow results to one FP/branch.
4. Document the predicate in the docstring so callers can replicate or combine it.

---

## Complete Usage Example

```python
import tanglepack, numpy as np


def henon_map(point):
    k, b = 10, 1
    x, y = point
    return np.array([y - k + x**2, -b * x])

def henon_map_inverse(point):
    k, b = 10, 1
    x, y = point
    return np.array([-y / b, x + k - (y**2) / (b**2)])


# ── Numeric phase ──────────────────────────────────────────────────────────
wb = tanglepack.TangleWorkbench(henon_map, henon_map_inverse)
fp = wb.construct_fixed_point([4, -4])
wb.orient_eigenvectors(fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])})
wb.initialize_both_manifolds(fp)
wb.grow_n_times(fp, "unstable", num_iterations=6)
wb.grow_until_turnaround(fp, "stable")

wb.compute_intersections(fp)          # also populates registry
wb.trim_stable_manifolds(fp)
bridges = wb.create_bridges(fp)

for _ in range(3):
    wb.iterate_all_bridges()

new_links = wb.infer_iterate_table(max_depth=8)
print(f"Recorded {new_links} iterate relationships")


# ── Topological phase ──────────────────────────────────────────────────────
registry = wb.intersection_registry
print(f"Total intersections: {len(registry)}")

# Lookup by ID
p = registry[3]
print(f"Intersection 3: coords={p.coords}, u_cdist={p.unstable_cdist:.4f}")
print(f"  Fixed point: {p.manifold_a_key[0].coordinates[0].ravel()}")
print(f"  Branch: {p.manifold_a_key[3]}")

# Iterate chains
fwd = registry.iterate_table.forward_chain(3)
print(f"Forward chain from 3: {fwd}")

# Orderings
u_order = registry.by_unstable_cdist
s_order = registry.by_stable_cdist

# --- Query interface ---

# Which intersections will map into the cdist range [5, 10] on the unstable manifold?
sources = registry.on_interval(5.0, 10.0, stability="unstable")
print(f"{len(sources)} intersections map into u-cdist [5, 10]")

# All intersections currently sitting in a cdist range (no iteration)
in_range = registry.on_cdist_range(2.0, 8.0, stability="unstable")
print(f"{len(in_range)} intersections have u-cdist in [2, 8]")

# All intersections involving the specific fixed point
from_fp = registry.from_fixed_point(fp)

# Only branch-0 intersections on the stable side
branch0_stable = registry.from_branch(0, stability="stable")

# Custom filter: only intersections whose iterate is also registered
has_iterate = registry.filter(lambda ix: (ix.id, 1) in registry.iterate_table)

# --- Live graph ---
G = wb.build_intersection_graph()
print(f"Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

# Adjacency edges along W^u
adj_u = [(u, v) for u, v, d in G.edges(data=True) if d.get("stability") == "unstable"]
# Iterate edges
iter_edges = [(u, v, d["n"]) for u, v, d in G.edges(data=True) if d.get("type") == "iterate"]

# Visualise (existing method still works)
wb.visualize_intersection_graph(G)

# --- Dense array exports ---
F = registry.as_forward_array(max_depth=5)   # shape (N, 5), -1 = unknown
B = registry.as_backward_array(max_depth=5)
```

---

## The Two Fundamental Structures

**Structure 1 — `registry.iterate_table[id, n]`**  
The forward/backward iterate map. `n > 0` forward, `n < 0` backward, `n = 0` identity.

**Structure 2 — `registry.by_unstable_cdist` and `registry.by_stable_cdist`**  
The two orderings that define the combinatorial template. The permutation mapping
`by_unstable_cdist` → `by_stable_cdist` drives all topological analysis.

**Structure 3 (new) — `registry.graph`**  
The live `nx.MultiDiGraph` with adjacency edges (reconstructed from the sorted lists)
and iterate edges (wired by `register_iterate`). Nodes carry coordinates and manifold
keys as attributes. Both `TangleWorkbench.visualize_intersection_graph` and any networkx
algorithms operate directly on this graph.

---

## Future Structures (not in this implementation phase)

### Lobe Table

Lobes are bounded regions between consecutive intersections on W^u and the corresponding
W^s arcs. The lobe table maps: lobe index → bounding intersection IDs, pre-image lobe,
and area. With `registry.by_unstable_cdist` giving the ordered boundary points and the
iterate table giving pre-images, lobe-table construction is a straightforward next step.

### Markov / Transition Matrix

An (N-1) × (N-1) matrix `M[i, j] = 1` if lobe i overlaps with lobe j under one iterate.
Its leading eigenvalue's log lower-bounds topological entropy. Follows directly from lobe
table + iterate table.

### Symbol Sequences and Rank Matrix

```python
# Rank matrix: R[i, d-1] = unstable_rank of f^d(ids[i])
def rank_matrix(registry: IntersectionRegistry, max_depth: int) -> NDArray[np.int64]:
    ids = registry.all_ids()
    R = np.full((len(ids), max_depth), fill_value=-1, dtype=np.int64)
    for i, iid in enumerate(ids):
        for d in range(1, max_depth + 1):
            target = registry.iterate_table[iid, d]
            if target is not None:
                R[i, d - 1] = registry.unstable_rank(target)
    return R
```

Symbol sequences drive topological entropy calculations and are the natural output of
the two structures + iterate table.

---

## Implementation Order

Each step depends only on what precedes it.

| Step | File | Type | Key change |
|------|------|------|-----------|
| 1 | `Intersection.py` | Modify | Add `ManifoldKey`, `label`, `manifold_a_key`, `manifold_b_key`; fix `from_segments` |
| 2 | `BaseManifold.py` | Modify | Add `manifold_key` attribute |
| 3 | `Tangle.py` | Modify | Extract and attach `ManifoldKey` in both populate methods |
| 4 | `IterateTable.py` | New | No dependencies |
| 5 | `IntersectionRegistry.py` | New | Depends on `Intersection`, `IterateTable` |
| 6 | `TangleWorkbench.py` | Modify | Set `manifold_key` on manifolds; wire registry; update compute/iterate methods |
| 7 | `__init__.py` | Modify | Export `IterateTable`, `IntersectionRegistry`, `ManifoldKey` |

Steps 4 and 1–3 can be done in parallel since `IterateTable` has no dependencies on the
rest of the library. Step 6 should come after all of steps 1–5 are complete.

---

## Summary of Changes

### New files

| File | Purpose |
|------|---------|
| `IterateTable.py` | `table[id, n]` → id; auto-reverse; chain methods; dense array export |
| `IntersectionRegistry.py` | Master store with IDs, collision dedup, sorted orderings, live graph, extensible filter interface |

### Modified files

| File | Changes |
|------|---------|
| `Intersection.py` | `ManifoldKey` TypeAlias; add `label`, `manifold_a_key`, `manifold_b_key` to `__init__`; fix `from_segments` signature; add `fixed_points` property |
| `BaseManifold.py` | Add `manifold_key: Optional[ManifoldKey]` attribute |
| `Tangle.py` | Generalise both populate methods to extract and attach `ManifoldKey`; no longer hardcode stable/unstable labelling |
| `TangleWorkbench.py` | Set `manifold.manifold_key` at manifold creation; add `_intersection_registry`; update `compute_intersections()` and `iterate_bridge()`; add `infer_iterate_table()`, `populate_registry()`, `build_intersection_graph()` |
| `__init__.py` | Export `IterateTable`, `IntersectionRegistry`, `ManifoldKey` |
