# Request: Blasting nested tangles through to symbolic dynamics

*From the author, 2026-10-02. This is the input for one `tangle_refactor` run.*

## The question

**Can a blasted tangle go through the whole pipeline (blast → trellis → holes and partition → minimal trellis → iterated partition → dual graph → symbolic dynamics) and come out correct? Does that still hold when the tangle is nested inside another one and both zones are blasted?**

I'm worried about errors that appear only after blasting, and only after blasting *different* zones of one session. The nested period 1 + period 3 Hénon case (`HENON_P3 = (2, 1)`) is the test bed:

1. **Each tangle on its own can take many blasts.** The period-3 orbit alone, and the period-1 saddle alone, should each survive many blasts (not just one or two) and still reach symbolic dynamics, with every invariant intact.
2. **The constructed nested case can do the same.** Blast the inner zone, the outer zone, or both, in any order and any reasonable number of times. Every stage should still run, and the answers should agree with each tangle studied on its own wherever they ought to.

This is a **study first** and a **refactor second**. Answer the question with evidence. Then, along the path the pipeline actually walks, look for places that could be faster, clearer, less complex and easier to read. Recommend the smallest change that fixes what you find.

## Where to start

| What | Where |
|---|---|
| Case recipes | `src/tanglepack/examples/henon_cases.py`: `build_period3(blasts=…)`, `build_nested(outer_blasts=…, inner_blasts=…)`, `TangleBuild.blast_sizes` |
| Blasting | `src/tanglepack/loom/Blast.py` (`blast_zone`, the frontier, `seen`, `min_separation`, `_resolve_zone`) and `TangleSession.blast_zone` |
| Zones | `src/tanglepack/loom/ResonanceZone.py` (`define_resonance_zone`, `_trim_stable_at`, `capture_boundary`, `restore`, `resolve_boundary_intersection_id`), `TangleSession.resonance_zone`, `add_resonance_zones`, `cartoon_zones`, `cartoon_enclosures` |
| Downstream | `Trellis` (pips, pseudoneighbors, holes, partition), `MinimalTrellis.py`, `PartitionFamily.py`, `DualGraph.py`, `DualWalk.py`, `SymbolicDynamics.py`, `BridgeClass.py` |
| Scripts and tests | `scripts/henon_symbolic_itineraries_period3.py` (`--blasts`), `scripts/henon_symbolic_itineraries_nested.py` (`--outer-blasts`, `--inner-blasts`), `scripts/symbolic_figures.py`, `tests/test_higher_period_cartoon.py`, `tests/test_loom_blast_restore.py`, `tests/test_stable_partition_period3.py` |

## What I already know (from `henon_cases.py` Dev Notes and `CLAUDE.md`)

- **Period 3 at 13 unstable steps is CLOSED.** The words are `a -> b -> c -> a u^-1 w^-1`, there are no image bridges, and blasting "changes nothing topologically". Verify this claim at high blast counts rather than taking it on trust. `figures/henon_symbolic_itineraries/p3_blast4/` and `p3_blast5/` exist from my own attempts.
- **Period 3 at 16 unstable steps trips `check_holes_share_bridge_side`.** Holes of one origin disagree on bridge side between iterate 0 and iterate −3. This is open and has not been investigated. If many blasts reach the same failure, say whether it is the same bug.
- **Nested needs two outer blasts for every class to resolve.** With none, the outer class `e` is a singleton landing with no registered chain. With two, the outer classes reproduce the k=2.8 two-blast words.
- **`build_nested` runs the inner blasts first**, "so the outer blasts do not consume its bridges". Both blasts pass `fixed_point=fixed_points` (both saddles). So an outer blast may also iterate the period-3 bridges lying inside the outer zone. Is that intended, harmless or the source of errors? Does the outcome depend on blast order?
- **`build_nested` uses one slug (`nested_p1_p3`) whatever the blast counts.** Figure runs at different counts overwrite each other. This is minor; note it if it matters.

## What "works" means

For each run, report pass/fail and the first failure with its traceback or warning:

- **Every stage runs.** Nothing raises from `blast_zone` through `symbolic_dynamics()`.
- **Every asserted invariant holds.** That covers the `CLAUDE.md` invariants list, `check_holes_share_bridge_side` and `check_bridge_rows_consistent`, and area-verified `image_of`.
- **No warnings that point to a bug.** Watch for same-stability (u×u or s×s) crossings logged and dropped, landings that straddle a cut (`contained=False`), `no empty stretch`, image pairs skipped with no `Bridge`, `virtual` (`newN`) classes, `ambiguous` walks, unresolved classes, `unmatched_members`, and `is_reliable == False`. Say which ones are expected and which are not.
- **The symbolic dynamics are stable under more blasting.** Once a class resolves, its word and refined classes should not change as blasts are added. If they do, explain why. More blasts may *refine* the answer (more registered images, more classes resolved), but should never contradict it.
- **The nested case agrees with each tangle alone.** The inner tangle's words in the nested session match `build_period3` at the same growth and blast count. The outer tangle's words match the period-1-alone run. Any connecting (heteroclinic) classes are explained.
- **Cost stays sensible.** Track trellis size, `blast_sizes`, wall time per stage, and where the time goes as blasts increase.

Suggested matrix (adjust if a cell is pointless or too slow, and say why):

- period 3 alone at 0, 1, 2, 4, 8 blasts
- period 1 alone at 0, 1, 2, 4, 8 blasts (the alone recipe may need writing as a probe script, not in the repo)
- nested at outer ∈ {0, 2, 4, 8} × inner ∈ {0, 2, 4, 8}, plus the reversed order (outer first) for at least one cell

## Refactor lens

While tracing the path above, look for code that is harder to read than the mathematics it computes (see `house_style`). Candidates I suspect, not conclusions:

- **The blast frontier and filters in `Blast.py`.** Containment by test point, the fixed-point filter normalisation, the `seen` guard and the proximity point cloud.
- **The `_repin` / `_partition` dance in `henon_cases.py`** that every blast has to be followed by. Should blasting leave the session in a consistent state by itself?
- **Zone lookup and identity across blasts.** `ResonanceZone` boundary capture and restore, `_anchor_id` and `_live_id`, and resolving boundary ids after the registry changes.
- **`TangleSession` (≈2000 lines).** Which parts of the blast and zone path could live closer to what they do, or be shorter?
- **Repeated work.** Rebuilds or recomputation that blasting triggers (cache keys on generation, pips and partition signature) that could be avoided.

Prioritise: correctness errors first, then clarity and complexity, then speed. Keep the recommendation to the **smallest** change that fixes what you found. Put everything else under *out of scope* for later runs.

## Constraints

- `CLAUDE.md` invariants are physics. Code may assert them and must never model their violation. In particular, a same-stability crossing is always a numerical artifact.
- Topology code resolves iterates by table lookup (`Trellis.iterate`, registered cdists) and never calls the map.
- Holes propagate backward only. The +1..+(k−1) exemption is provisional, but stays as-is for this run unless it is shown to cause a nested-blast error.
- Pip-branch-only unification for period-k stays (author, 2026-09-30).
- `numerics` and `topology` never import `loom`.
- The working tree has uncommitted `.claude/` changes and untracked figure folders. They are unrelated to this run, and the change branches from `HEAD`.
