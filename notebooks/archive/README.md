# Archived notebooks

These notebooks are kept for historical reference only. **None of them run against
the current library** — every one calls an API that was deleted long before the
`regions-and-identity` refactor. Do not use them as examples; start from
`notebooks/henon_nested_period_3.ipynb` or the `scripts/*.ipynb` walkthroughs instead.

| Notebook | Why it was archived |
|---|---|
| `curvature_area_book.ipynb` | Uses the pre-`TangleWorkbench` API: `tanglepack.InitialManifold`, `tanglepack.Manifold.linear_fit/parabolic_fit`, and `FixedPoint(map, guess, ...)` as a solver. All removed. |
| `initial_manifold_book.ipynb` | Same removed API (`InitialManifold.iterate_manifold`, `manifold.points.points`); superseded by `ManifoldInitializer` + `ManifoldMachine`. |
| `intersection_sandbox.ipynb` | Uses `tanglepack.FixedPoint2` / `tanglepack.Manifold2` (never-shipped prototypes) and hand-rolls an rtree/shapely intersection search that `Tangle` now owns. |
| `period_three.ipynb` | Stale signatures (`construct_fixed_point(guess, 2)`, `grow_x_times(manifold, fp, n)`) and duplicates `scripts/henon_period_3_manifolds.py`, which is maintained. |
