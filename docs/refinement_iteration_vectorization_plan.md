# Refinement & Iteration Vectorization Plan

## Goal

Speed up the two hot paths in `ManifoldMachine` — **manifold refinement** (profiled at
~80–90% of total runtime) and **iteration** — without changing the data structure or the
mathematical results. Three concrete changes:

1. **Closed-form curvature area** — replace the per-pair `np.linalg.inv` of a 3×3
   Vandermonde matrix in `_parabolic_fit`/`_curvature_area` with a handful of scalar ops
   derived in closed form. Mathematically identical, no LAPACK dispatch.
2. **Vectorized refinement sweep** — restructure refinement from depth-first scalar
   recursion into breadth-first *waves*: one vectorized sweep computes every edge's
   curvature area at once (decision sweep), a second vectorized sweep maps all required
   midpoints at once (iteration sweep), then insertions are applied to the linked list.
3. **Vectorized iteration** — collapse the per-point `[map_fwd(p) for p in coords]`
   comprehensions in `iterate_manifold` / `_iterate_without_refine` into single batched
   map calls. (This path *looked* vectorized but is actually a Python loop.)

### Non-goals

- **No GPU.** Everything here is NumPy on CPU. GPU is explicitly out of scope; see the
  note at the bottom for how this design makes a future GPU port trivial.
- **No change to the data structure.** The doubly-linked `Point` lists stay. O(1)
  insertion is the right tool for an insertion-dominated workload; insertions remain a
  short Python loop (cheap) — only the map evaluations and curvature math get vectorized.
- **No change to the refinement criterion.** The area-cutoff test produces the same
  accept/reject decisions (up to floating point) it does today.

---

## Background — where the time actually goes

Per refinement step (`refine_two_points` → `_curvature_area` + `_get_refined_point`):

- `_curvature_area` does **two `np.linalg.inv` calls on 3×3 matrices** (`_parabolic_fit`,
  once for the left stencil and once for the right) plus several `np.vstack` allocations.
  It calls **no map** — it is pure geometry on coordinates already in the list.
- `_get_refined_point` does **exactly one `map_fwd`** on a single point, plus a `Point`
  allocation.

So the decision (neighbor-coupled, map-free) and the new-point computation (map-heavy,
edge-local) are cleanly separable. That separation is what makes both vectorization and
the closed-form rewrite possible. For analytic maps the `np.linalg.inv` overhead is a
large share of the 80–90%; for expensive maps the batched map sweep is the win. This plan
captures both.

Relevant locations (line numbers approximate, `src/tanglepack/ManifoldMachine.py`):

| Method | ~Line | Role |
|---|---|---|
| `iterate_manifold` | 220 | per-point map loop at the `np.vstack([... for p in ...])` |
| `_iterate_without_refine` | 477 | same loop; also two stray `print` statements |
| `refine_manifold` | 716 | outer DFS walk over adjacent pairs |
| `refine_two_points` | 779 | per-pair adaptive bisection (the stack) |
| `_get_refined_point` | 862 | midpoint preiterate → map → new `Point` |
| `_curvature_area` | 1051 | left+right parabola areas, returns max |
| `_parabolic_fit` | 1029 | 3×3 Vandermonde + `np.linalg.inv` |
| `_compute_single_area` | 1118 | analytic ∫(parabola − line) |
| `_linear_fit` | 1008 | chord through the two endpoints |

---

## Part A — Closed-form curvature area

### The math

For an edge with endpoints `p0 = (x0, y0)`, `p1 = (x1, y1)` and a third stencil point
`p2 = (x2, y2)` (the left neighbor for the left area, the right neighbor for the right
area), the current code:

1. fits a parabola `P(x) = a x² + b x + c` through the three points (Vandermonde inverse),
2. fits the chord line `L(x) = m x + d` through `p0, p1`,
3. integrates `P(x) − L(x)` over `[x0, x1]`.

We never need the coefficients. Because `P` and `L` agree at `x0` and `x1`, their
difference is a quadratic vanishing at both, so it factors exactly:

```
P(x) − L(x) = a · (x − x0)(x − x1)
```

where `a` is the parabola's leading coefficient, which equals the **second divided
difference** of the three points:

```
a = f[x0, x1, x2]
  = ( (y2 − y1)/(x2 − x1) − (y1 − y0)/(x1 − x0) ) / (x2 − x0)
```

Integrating the factored form over the segment (with `h = x1 − x0`) uses the standard
identity `∫_{x0}^{x1} (x − x0)(x − x1) dx = −h³ / 6`, so:

```
area = |a| · h³ / 6
```

That is the entire curvature area — no matrix, no LAPACK, exactly equal (up to floating
point) to today's `_parabolic_fit` + `_compute_single_area` pipeline. The left area uses
`p2 = left neighbor`, the right area uses `p2 = right neighbor`, and `_curvature_area`
returns the larger, exactly as now.

### Scalar implementation (drop-in, keeps the existing control flow)

Add a helper and rewrite `_curvature_area` to call it; `_parabolic_fit`,
`_compute_single_area`, and `_linear_fit` become unused by refinement and can be removed
(or kept only if referenced elsewhere — grep first).

```python
@staticmethod
def _segment_curvature_area(p0, p1, p2) -> float:
    """Area between the chord p0-p1 and the parabola through p0, p1, p2,
    integrated over [x0, x1]. Closed form via the second divided difference."""
    x0, y0 = p0
    x1, y1 = p1
    x2, y2 = p2
    h = x1 - x0
    denom = (x1 - x0) * (x2 - x1) * (x2 - x0)
    if denom == 0.0:               # near-vertical / coincident x: matches current skip
        return 0.0
    a = ((y2 - y1) * (x1 - x0) - (y1 - y0) * (x2 - x1)) / denom
    return abs(a) * abs(h) ** 3 / 6.0
```

### Equivalence & edge cases

- **Identical values.** The result equals the current `max(left_area, right_area)` up to
  floating point. A regression test asserts `np.isclose` against the old implementation on
  a battery of random triples *before* the old code is deleted.
- **Near-vertical segments.** Today `_parabolic_fit` raises `LinAlgError` on a singular
  Vandermonde and `refine_two_points` catches it and `continue`s (skips). The closed form
  reproduces this exactly: the denominator `(x1−x0)(x2−x1)(x2−x0)` is the Vandermonde
  determinant; a zero (or near-zero) denominator → return `0.0` → below cutoff → skip.
  This preserves the existing (known) limitation that vertical segments are
  under-refined; fixing that — e.g. an arc-length parametrization — is out of scope.
- **Boundary edges.** When `p0` is the root (no left neighbor) the left area is `0`; when
  `p1` is the tail (no right neighbor) the right area is `0`. Same as today.

---

## Part B — Vectorized map contract (the linchpin)

Both the iteration sweep (Part C) and the refinement map sweep (Part D) require the map to
accept a batch. Today `ManifoldView.map_fwd` is `system.map` (or `system.map_inv`) and is
called on a single `(2,)` point.

**Decision: standardize the internal map contract on `(N, 2) → (N, 2)`.** A single point is
just `N = 1`.

Changes:

1. **`DynamicalSystem`** — document the batched contract. Add a thin adapter
   `map_batch(points)` / `map_inv_batch(points)` that guarantees 2-D in/out:
   `np.atleast_2d(points)` on the way in, and a fallback that applies a scalar map
   row-wise (`np.stack([f(p) for p in pts])`) when a user map is not yet vectorized, so
   nothing breaks during migration. Vectorized maps hit the fast path; legacy scalar maps
   still work (just no speedup).
2. **`tanglepack_webdash/maps.py`** — emit vectorized lambdas. Write expressions against
   the last axis (`p[..., 0]`, `p[..., 1]`) and return `np.stack([...], axis=-1)` so they
   broadcast over `(N, 2)` unchanged. Existing single-point call sites keep working
   because `(2,)` is `(N=1, 2)` under `atleast_2d`.
3. **Single-point callers** (`FixedPointSolver`, `ManifoldInitializer`,
   `_get_iterate`) — route through the adapter and index `[0]` when they need a scalar
   result. No behavior change.

`ManifoldView` exposes `map_fwd`/`map_back` bound to the batched versions, keeping the
stability-direction wiring it already has.

---

## Part C — Vectorize the iteration sweep

In `iterate_manifold` (and the twin loop in `_iterate_without_refine`) `non_iterated_coords`
is *already* an `(N, 2)` array (`BaseManifold.get_non_iterated_point_array` returns
`np.vstack(...)`), but it is consumed by a Python loop:

```python
# before
iterated_points = np.vstack([viewer.map_fwd(p) for p in non_iterated_coords])
# after
iterated_points = viewer.map_fwd(non_iterated_coords)   # one batched call
```

That single line is the whole iteration speedup once Part B is in place. Also in this pass:

- Replace the two stray `print(f"Num Incorrectly labeled points: ...")` statements in
  `_iterate_without_refine` with `logger.debug(...)` (per the project's no-bare-print rule).
- Guard the empty case: `viewer.map_fwd` on a `(0, 2)` array must return `(0, 2)`; the
  adapter handles this so the existing `if len(non_iterated_coords):` branch is unchanged.

No structural change — the surrounding `Point` construction, iterate-link insertion, and
merge/refine calls stay exactly as they are.

---

## Part D — Vectorized wave-based refinement

Restructure refinement from depth-first scalar recursion into breadth-first **waves**.
Each wave is two vectorized sweeps plus one insertion loop:

1. **Decision sweep (vectorized, map-free).** Materialize the current working set of edges
   into arrays and compute every edge's curvature area at once using the Part A closed
   form. Produce a boolean mask of edges over `area_cutoff`.
2. **Iteration sweep (vectorized, batched map).** For the masked edges, build the midpoint
   pre-iterate coordinates `0.5 * (pre0 + pre1)` as one `(M, 2)` array and map them all in
   a single `viewer.map_fwd` call.
3. **Insertion (Python loop, O(1) each).** Create the new `Point`s, cache their
   preiterates, and splice each into the list between its `p0` and `p1`. Collect the child
   edges `(p0, new)` and `(new, p1)` as the next wave's working set.

Repeat until a sweep finds no edge over cutoff.

### Why insertions are safe within a wave

- New-point coordinates depend only on each edge's two endpoints' preiterates (read in the
  sweep, before any write), so there is no read-after-write hazard — the batch is
  order-independent.
- Adjacent edges `(a,b)` and `(b,c)` share node `b`, but inserting into `(a,b)` writes
  `b.backward` while inserting into `(b,c)` writes `b.forward` — different fields. So the
  full frontier can be refined in one wave; no independent-set selection is needed.
- The only neighbor coupling is in the **decision** (the left/right stencil), and that is
  resolved by computing all decisions against the frozen pre-wave snapshot.

### Sketch

```python
def refine_manifold_waves(self, manifold, *, full_rescan=False):
    viewer = ManifoldView(manifold, self.system)
    edges = self._materialize_edges(manifold)   # working set: list of (p0, p1) nodes

    while edges:
        # ---- 1. DECISION SWEEP (vectorized, no map) ----
        P0  = np.array([e[0].get_point() for e in edges])      # (E, 2)
        P1  = np.array([e[1].get_point() for e in edges])
        L   = self._neighbor_coords(edges, side="left")        # (E, 2), NaN if none
        R   = self._neighbor_coords(edges, side="right")
        PRE0 = self._preiterate_coords(edges, 0, manifold.stability)   # NaN if missing
        PRE1 = self._preiterate_coords(edges, 1, manifold.stability)

        area = np.maximum(
            self._curvature_area_vec(P0, P1, L),   # left stencil
            self._curvature_area_vec(P0, P1, R),   # right stencil
        )                                           # NaN stencils -> 0 inside helper

        too_close = np.einsum("ij,ij->i", P1 - P0, P1 - P0) < 1e-16
        pre_ok    = np.isfinite(PRE0).all(1) & np.isfinite(PRE1).all(1)
        split     = (area >= self.area_cutoff) & ~too_close & pre_ok
        if not split.any():
            break

        # ---- 2. ITERATION SWEEP (vectorized, one batched map) ----
        mid_pre = 0.5 * (PRE0[split] + PRE1[split])            # (M, 2)
        new_xy  = viewer.map_fwd(mid_pre)                      # single call
        cdist   = 0.5 * (
            np.array([e[0].cdist for e in edges])[split]
            + np.array([e[1].cdist for e in edges])[split]
        )

        # ---- 3. INSERTION (O(1) splices) ----
        next_edges = []
        for (p0, p1), pre, xy, cd in zip(
            [e for e, s in zip(edges, split) if s], mid_pre, new_xy, cdist
        ):
            new = Point(xy[0], xy[1], float(cd), stretch_param=p0.stretch_param)
            self._cache_preiterate(new, pre, manifold.stability)
            self._insert_point_geometrically(p0, new, manifold)
            next_edges += [(p0, new), (new, p1)]

        edges = self._materialize_edges(manifold) if full_rescan else next_edges
```

`_curvature_area_vec(P0, P1, P2)` is the array form of Part A: compute the second divided
difference across all rows, return `|a| * |h|³ / 6`, with rows whose denominator is ~0 or
whose `P2` is NaN forced to `0.0` (the boundary / vertical skip cases).

### `full_rescan` flag

- **`full_rescan=False` (default):** each wave processes only the children of the previous
  wave's splits. Fast (work ∝ points actually added). Evaluation order differs from
  today's DFS, so the exact final point set differs slightly — benign for an adaptive
  heuristic, since both converge to "every edge under cutoff."
- **`full_rescan=True`:** every wave re-scans the entire edge list. Guarantees that *no*
  final edge sits over cutoff (catches an edge pushed over by a neighbor's new point).
  Costs extra decision sweeps — but those are map-free and cheap, especially after Part A.

Keep the existing `refine_manifold` as the scalar fallback during migration; switch
callers (`grow_manifold`, `iterate_manifold`, `new_grow_manifold`) to
`refine_manifold_waves` behind a flag, then make it the default once parity is confirmed.

### Helpers to add

- `_materialize_edges(manifold)` — one walk producing the ordered list of adjacent
  `(p0, p1)` node pairs (reuses the existing walk logic).
- `_neighbor_coords(edges, side)` — left/right stencil coords aligned to `edges`, `NaN`
  where the neighbor is absent (root/tail).
- `_preiterate_coords(edges, which, stability)` — gather `prev_iterate`/`next_iterate`
  coords (stability-aware, via `_get_preiterate`), `NaN` where missing — this is the
  vectorized form of the existing "skip if preiterate missing" guard.
- `_curvature_area_vec(...)` — vectorized Part A.

---

## Correctness & testing

1. **Closed-form equivalence (Part A).** Random-triple test: assert
   `np.isclose(_segment_curvature_area(...), old _curvature_area(...))` across many
   stencils, including near-vertical (both should skip) and boundary (area 0) cases.
   Land this *before* deleting `_parabolic_fit` / `_compute_single_area`.
2. **Iteration equivalence (Part C).** A vectorized `map_fwd` over `(N,2)` must equal the
   row-wise loop element-for-element; assert on a real map from a test script.
3. **Refinement parity (Part D).** On a fixed map + fixed `area_cutoff`, compare
   `refine_manifold` vs `refine_manifold_waves(full_rescan=True)`: every final edge under
   cutoff, and point counts/coords within tolerance. `full_rescan=False` is checked for
   the weaker invariant (all *frontier* edges under cutoff) plus a visual/area sanity pass.
4. **End-to-end.** Re-run `scripts/tangle_workbench_test.py` and the existing
   `tests/` suite; intersection counts and bridge structure must be unchanged.
5. **Profile.** Re-profile after Part A alone, then after Parts C+D, to confirm the
   refinement share drops as expected and to see the new dominant cost.

---

## Rollout order

1. **Part A** — closed-form `_curvature_area` (pure win, no contract change, smallest
   blast radius). Ship behind nothing; it's a drop-in with a regression test.
2. **Part B** — batched map contract + `maps.py` vectorization + adapter fallback.
3. **Part C** — one-line iteration vectorization + print→logger cleanup.
4. **Part D** — `refine_manifold_waves` behind a flag; flip to default after parity tests.

Each step is independently valuable and independently testable.

---

## Future note — GPU

This design makes a GPU port mechanical, with **no further restructuring**:

- The only heavy numerical kernels are now (a) the batched `viewer.map_fwd` over `(M, 2)`
  in the iteration sweep, and (b) the elementwise `_curvature_area_vec` over `(E, 2)` in
  the decision sweep. Both are pure array math with no Python-level per-element work.
- To run on GPU, swap the array backend (CuPy or JAX) inside `DynamicalSystem.map_batch`
  and the two `*_vec` helpers — `np.*` → `xp.*` — and keep the host-side linked-list
  insertion loop on the CPU. The wave structure already batches all device work, so the
  host↔device boundary is exactly the wave boundary: upload the wave's coordinate arrays,
  run the kernel, download results, splice.
- The win is largest for expensive maps (each `map_fwd` heavy) and for long manifolds
  (large `M`/`E` per wave amortizing launch overhead). With JAX you additionally get
  `jit`+`vmap` on the map and free Jacobians for `FixedPointSolver`.

Until then, everything above is plain NumPy on CPU and already a substantial speedup.
