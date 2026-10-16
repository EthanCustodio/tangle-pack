# Plan — B9 `session-three-concerns` (planner, 2026-10-15, base `517aa01`)

## In short

`TangleSession` (1,744 lines) does four jobs: it is the facade, it forwards plotting calls, it computes the cartoon's geometry, and it fans stage calls out over fixed points. B9 keeps the facade and removes the other three.

**What leaves the session:**
- **Plot delegates.** All thirteen `plot_*` methods go. Twelve only forward a call; callers will call `topology.plotting` or `Trellis.plot_*` directly. The thirteenth, `plot_resonance_zones`, has its own body, so it moves.
- **Cartoon geometry.**
  - `cartoon_zones` moves unchanged to a new loom module, `loom/cartoon.py`. It needs loom's `ResonanceZone` and `ring_order`, so it cannot live in topology.
  - `cartoon_enclosures` moves to `topology.plotting`. It reads only the trellis, the manifolds and the fixed points.
- **Fan-outs.** `_fanout_call`, `_fanout_plot`, the manual stage calls (`compute_pseudoneighbors`, `punch_holes`, `partition_stable_manifold`), `partition_element_for`, and the describe fan-outs are deleted. Their callers read `session.topology(fp)` or the trellis.

**What stays on the session:** the snapshots, `topology()` and its stage readers, the pip memory, zones, `classify_bridge`, `blast_zone`, and the delegation to the workbench.

**Nothing computed changes.** The analyst's probe found that the four "element at a bridge end" readings agree on all 63 corpus items (0 disagreements). Building the zones from pips would change one circle cartoon, so the zones stay as they are stored (Q2).

**The split.** One branch would touch about 30 files. I propose two slices:
- **B9a (this run):** the session's three concerns. 24 files, almost all of them call-site rewrites.
- **B9b (next run):** fold the four element-at-end implementations into `owner_at_end`. This touches the class rule and an import cycle, so it needs a different kind of review. Skeleton at the end.

**What you decide:** Q0 (the split) and Q1–Q6. Each has a default, so "approved" accepts all of them.

---

## 1. The answer

Three of the four concerns can leave with no change to anything computed or drawn:
- About 860 lines leave the session: 438 of plot delegates, 180 of cartoon geometry, 204 of fan-outs and 39 of describe fan-outs.
- Neither cartoon function reads the pip memory. `cartoon_zones` reads the stored zones and the manifolds; `cartoon_enclosures` reads a trellis, the manifolds and the fixed points.
- The four element-at-end readings agree on all 7,090 row ends, all 14,180 (end, side) pairs and all 1,776 iterated ends across the 63 corpus items.
- Computing the zones from `Trellis.pips` matches the stored zones on all 41 cases that store a zone. It would change `inversion_period2`'s circle cartoon, though, and it would pull `ResonanceZone.py` into the slice. So that part of the plan row is declined (Q2 (a)).

## 2. Recommendation

Change it, as B9a now, with B9b next.

## 3. Scope (B9a)

**Slug:** `session-three-concerns`.
- Branch: `refactor/session-three-concerns`.
- Worktree: `../tangle-pack-session-three-concerns`, created off `517aa01`.
- Base probe worktree: `../tangle-pack-session-three-concerns-base`, detached at `517aa01`.
- Suite inside a worktree: `MPLBACKEND=Agg PYTHONPATH=src ../tangle-pack/env/bin/python -m pytest -q -rxX -p no:cacheprovider`.

**File list.** The diff may touch these files and no others. There are 24, but 20 of them are call-site rewrites of a few lines each.

| Group | Files |
|---|---|
| src (4) | `src/tanglepack/loom/TangleSession.py`; `src/tanglepack/loom/cartoon.py` (**new**); `src/tanglepack/loom/__init__.py`; `src/tanglepack/topology/plotting.py` |
| tests (8) | `tests/plotting/test_plot_smoke.py`; `tests/plotting/test_plot_topology.py`; `tests/facade/test_session_fanouts.py`; `tests/facade/test_session_bridge_classes.py`; `tests/facade/test_session_caches.py`; `tests/regression/test_nested_own_blast.py`; `tests/unit/topology/test_arrangement.py`; `tests/unit/topology/test_partition_elements.py` |
| scripts (11 + 1) | `scripts/symbolic_figures.py`, `henon_pseudoneighbors.py`, `henon_pseudoneighbors_nested.py`, `henon_strong_pip_nested.py`, `henon_minimal_dual_graph.py`, `henon_dual_graph.py`, `henon_k10_dual_graph.py`, `henon_bridge_classes.py`, `henon_resonance_zone.py`, `henon_bridge_classification.py`, `henon_blast_period_3.py`; `scripts/tangle_workbench_test.ipynb` (one source line only, outputs untouched) |

**Never touched:**
- `ResonanceZone.py`, `Topology.py`, `Blast.py`.
- `topology/__init__.py`: `cartoon_enclosures` is reached as `plotting.cartoon_enclosures`, and `__all__` stays at 48 names.
- `tanglepack/__init__.py`.
- `tests/helpers/`, `tests/cases.py`, `tests/invariants/`, `tests/conftest.py`.
- `CLAUDE.md`, `docs/`, `.refactor/`, `artifacts/`, `figures/`, and the two untracked sweep scripts.

**Size:**
- `TangleSession.py` goes from 1,744 to about 848 lines: the analyst's 856, minus about 8 for unused imports. **Cap: ≤ 870.** If you answer Q4 "keep", the cap becomes ≤ 900.
- `loom/cartoon.py` is new, at about 175 lines.
- `plotting.py` grows by about 75 lines (cap +85), only for what moved.

**Net for B9a:**

| Part | Lines |
|---|---|
| src | about −650 (−896 deleted, +250 moved) |
| tests | about −150 |
| scripts | about ±20 |
| **Total** | **about −800**, about 250 of them moved |

## 4. Steps — three commits, one per concern, each green

### C1 "Cartoon geometry and the zones plotter leave the session"

`loom/cartoon.py` (new) holds two functions. Its module docstring says why it is in loom: it takes `ResonanceZone`s and `ring_order`, and topology never imports loom.

```python
def cartoon_zones(zones: Iterable[ResonanceZone], workbench) -> tuple[dict, dict, dict]:
    # body of TangleSession.cartoon_zones :781-889, verbatim except:
    #   self.workbench -> workbench
    #   sorted(self.resonance_zones.values(), …) -> sorted(zones, …)

def plot_resonance_zones(zones: Iterable[ResonanceZone], ax=None, colors=None,
                         alpha: float = 0.3, **fill_kwargs) -> list:
    # body of :1688-1742; the boundary_polygon(workbench) fallback goes
    # (define_resonance_zone always captures, ResonanceZone.py:547)
```

- `loom/__init__.py` exports `cartoon_zones` and `plot_resonance_zones`, and its docstring's entry-point list gains one line.
- `topology/plotting.py` gains `cartoon_enclosures(trellis, bridge_ids, parents, fixed_points)`.
  - The body is `:709-779`, with `workbench.manifolds` replaced by `trellis.manifolds` (the same live dict, `Trellis.py:203`).
  - `workbench.fixed_points` is replaced by the `fixed_points` argument. `test_plot_topology.py:419`'s sub-selection needs the child fixed point, which its trellis does not hold.
  - The docstrings at `plotting.py:6`, `:2519` and `:2896` stop citing the session.
- Deleted from `TangleSession`: `plot_dual_graph_cartoon` (:663-707), `cartoon_enclosures`, `cartoon_zones` and `plot_resonance_zones`.
- A circle cartoon is now drawn by passing the three zone values and the enclosures explicitly. Compute the enclosures **only when `parents` is non-empty**, as `:703` did:

```python
interior, parents, order = cartoon_zones(session.resonance_zones.values(), session.workbench)
extra = {"enclosures": plotting.cartoon_enclosures(dual.trellis, dual.minimal.kept_bridge_ids,
                                                   parents, session.fixed_points)} if parents else {}
plotting.plot_dual_graph_cartoon(dual, dynamics, ax=ax, shape="circle", interior_side=interior,
                                 parents=parents, arc_order=order, **extra, **style)
```

- Callers:
  - `scripts/symbolic_figures.py:225`. The circle branch is as above; a `shape="line"` call passes no zone values.
  - `tests/plotting/test_plot_topology.py:391`, `:400`, `:412`, `:415`, `:419`, `:422`.
  - `tests/plotting/test_plot_smoke.py:311-320`.
  - `plot_resonance_zones(session.resonance_zones.values(), alpha=…)` at `henon_resonance_zone.py:92`, `henon_bridge_classification.py:145`, `henon_blast_period_3.py:136` and the notebook cell (`tangle_workbench_test.ipynb`, json line 500).
  - `henon_blast_period_3.py` is still dead at `:95`. Rewrite that one line, report it as still dead, and do not revive or delete the script.

### C2 "Plot delegates leave the session"

Delete these from `TangleSession`:
- `plot_transition_graph`, `plot_walk`, `plot_bridges_by_class`, `plot_itinerary_table`, `plot_dual_graph`, `plot_minimal_trellis`.
- `plot_strong_pip_candidates`, `plot_strong_pip`, `plot_pseudoneighbors`, `plot_holes`, `plot_stable_partition`.
- `_fanout_plot` (:999-1039).

Rewrite the callers as the analyst's table says:
- `symbolic_figures.py:181` → `session.trellis(fp).plot_strong_pip(ax=ax, s=120, …)`.
- `henon_pseudoneighbors.py:95-97` → `trellis.plot_*` (`trellis` is already a local at :89).
- `henon_pseudoneighbors_nested.py:97-98` → `for fp in (fp3, fp1): session.trellis(fp).plot_pseudoneighbors(s=40); session.trellis(fp).plot_holes()`.
- `henon_strong_pip_nested.py:81-82` → the same loop. Also update the docstring at `:6` and the comment at `:80`.
- `henon_minimal_dual_graph.py`:
  - `:118-119` → `session.trellis(fp).plot_holes` / `.plot_strong_pip`.
  - `:149` → `plotting.plot_dual_graph(dual, ax=…, clip_to_arcs=False)`.
- `henon_dual_graph.py:244` → `plotting.plot_dual_graph(…)`.
- `henon_k10_dual_graph.py`:
  - `:219` → `plotting.plot_dual_graph`.
  - `:225`, `:226`, `:230`, `:239` → `session.trellis(fp).plot_*`.
- `henon_bridge_classes.py:253`, `:267`, `:332` → `session.trellis(fp).plot_strong_pip` / `.plot_stable_partition`.

Tests:
- Delete `test_plot_smoke.py:59-71`, `:74-80` and `:277-308`.
- Delete `test_session_bridge_classes.py:151-168` (Q5).

### C3 "Fan-outs and `partition_element_for` deleted"

Delete these from `TangleSession`:
- `_fanout_call` (:967-997).
- `compute_pseudoneighbors`, `punch_holes`, `partition_stable_manifold`.
- `partition_element_for` (:1362-1412).
- `describe_holes`, `describe_stable_partitions`, `_describe` (Q4).
- The imports nothing reads any more: `Callable`, `endpoint_index`, `ring_order`, `CartoonLayout`, `Walk`, `Table`, `Line2D`, `networkx`, and `plotting` if it is unused.

Keep `_resolve_fixed_points`, because `Blast.py:309` reads it.

Scripts:
- `henon_pseudoneighbors.py:83` → `trellis.compute_pseudoneighbors(verbose=True)`.
- `henon_dual_graph.py`:
  - `:83-85` → `session.topology()`.
  - `:134-136` → `for fp in fps: session.topology(fp)`. A per-fixed-point call fills each trellis the way the fan-out did; a whole-selection `topology()` would leave them empty.
- `henon_k10_dual_graph.py`:
  - `:186-188` → `session.topology(fp)`.
  - `:196` → `session.trellis(fp).describe_stable_partitions()`.

Tests:
- `test_session_fanouts.py`:
  - Keep only the `classify_strong_pips` halves of `:17-44` and `:47-65`.
  - Rename the two tests to `test_classify_strong_pips_covers_every_fixed_point` and `test_classify_strong_pips_of_one_fixed_point_is_its_list`.
  - Delete `:68-83`.
  - Rewrite the module docstring. The file name stays.
- `test_session_bridge_classes.py:36-38`: delete the three lines from `_repartition_at`.
- `test_session_caches.py:135-137` → `trellis.compute_pseudoneighbors(); trellis.punch_holes(); trellis.partition_stable_manifold()` on the local `trellis` (:127).
- `test_nested_own_blast.py:78-80`: delete.
- `test_arrangement.py:316-318` → `session.topology(fp)`.
- `test_partition_elements.py`:
  - Drop line `:162`, the `partition_element_for` assert. `:161`'s `element_for` assert stays until B9b.
  - Delete `test_session_closes_the_cross_trellis_gap` (:165-191).

## 5. Behaviour

None. No word, matrix, count, `is_reliable`, figure or script stdout changes.

The only reachable difference is that a session method someone calls is gone, which raises `AttributeError` through `__getattr__`'s workbench delegation.

`plot_resonance_zones` loses a fallback that no `define_resonance_zone` zone can reach.

## 6. Not doing

- **Zones from pips** (Q2 (b)). It would change `inversion_period2`'s circle cartoon and move the boundary builder and `ring_cycles` / `ring_order` into topology.
- **B9b's element-at-end fold.** `Trellis.element_for`, `_elements_at_bridge` / `elements_at_bridge` and `BridgeClass._element_at_end` stay as they are.
- **A one-call circle-cartoon helper.** Three callers spell out the inputs instead; ask if you want one in `loom/cartoon.py`.
- **`Trellis.plot_*` thin delegates** (SIMP-11 / B12), and the plotting split.
- **Stale CLAUDE.md rows** (B12): the `TangleSession` row lists the delegates and `cartoon_zones`.
- **The four stale committed PNGs** (`p3`, `nested`, `heteroclinic`, `inversion`).
- **Reviving `henon_blast_period_3.py`.**
- **`classify_strong_pips`'s dict-or-list return.** It is pip memory and stays.

## 7. Rubric

Each item is checked by running it, never by reading a report. "Base" is `../tangle-pack-session-three-concerns-base` at `517aa01`; "branch" is the worktree's tip.

1. **The suite passes on the branch.** `-rxX` lists exactly one xfail (C3) and no XPASS.
2. **Exactly the planned tests are gone.** The collected ids differ from base by these deletions only:
   - `test_plot_smoke.py::test_session_plot_fanouts_draw_every_fixed_point`
   - `::test_session_plot_pseudoneighbors_computes_on_demand`
   - `::test_session_plot_delegates_return_their_types`
   - `test_session_fanouts.py::test_describe_fanouts_mention_every_fixed_point`
   - `test_session_bridge_classes.py::test_a_pip_trimmed_away_is_not_replaced_by_plotting[plot_strong_pip]` and `[plot_strong_pip_candidates]`
   - `test_partition_elements.py::test_session_closes_the_cross_trellis_gap`

   …and by these renames only:
   - the two `test_session_fanouts` tests (§4 C3);
   - `test_session_cartoon_delegate_honours_its_toggles` → `test_cartoon_honours_its_toggles`.

   The expected count is **1208 passed / 47 skipped / 5 deselected / 1 xfailed**. Any other difference fails this item.
3. **Every rewritten test that uses only base API passes on base too:** the `test_session_fanouts` pair, `test_session_caches`, `test_arrangement::test_stable_arcs_report_their_partition_elements`, `test_nested_own_blast`, `test_session_bridge_classes`, `test_partition_elements`, and the toggles test once it is written against `plotting.plot_dual_graph_cartoon`. The two `test_plot_topology` circle tests need `loom.cartoon` and are exempt; rubric 9 pins their figures instead.
4. **Teeth: each rewrite keeps its assertions.** Read every rewritten test against base:
   - **`test_cartoon_honours_its_toggles`** keeps all three assertions, now on `plotting.plot_dual_graph_cartoon`: the result is a `CartoonLayout`, `bridges_drawn == walks_drawn == 0`, and `not ax.patches`.
   - **`test_p3_circle_follows_the_zone_boundary`** keeps every assertion:
     - `parents == {}`;
     - the interior keys equal the branch cycle;
     - the arc order is a permutation of the branch cycle;
     - the cyclic anchor order;
     - `_assert_interior_inside`.
   - **`test_nested_inner_circle_sits_inside_the_outer_one`** keeps every assertion:
     - `parents`;
     - the enclosures are non-empty and every one equals the inner fixed point;
     - the sub-selection `alone` gives non-empty enclosures;
     - the circle-in-keep-out inequalities;
     - the node-inside loop.
   - **The `classify_strong_pips` halves** are verbatim.
   - **`test_session_caches`** still repartitions at the alternate pip on the local trellis; its later assertions are unchanged.
   - **`test_arrangement`'s** element loop still checks at least one arc. The engineer checks it is not vacuous on the branch.
5. **The facade is clean.** `grep -n "def plot_\|def cartoon_\|_fanout\|def partition_element_for\|def describe_holes\|def describe_stable_partitions\|def _describe" src/tanglepack/loom/TangleSession.py` prints nothing.
6. **The session is ≤ 870 lines.** `wc -l src/tanglepack/loom/TangleSession.py` reports ≤ 870 (≤ 900 if Q4 keeps the describe fan-outs).
7. **No one calls a deleted method.** `grep -rn "session\.\(plot_\|cartoon_\|compute_pseudoneighbors\|punch_holes\|partition_stable_manifold\|partition_element_for\|describe_holes\|describe_stable_partitions\)" src tests scripts stress` prints nothing. `classify_strong_pips` and the stage readers are not matched.
8. **The scope holds.**
   - `git diff --stat 517aa01..` lists only the 24 files of §3.
   - The notebook's diff is one source line.
   - `cartoon_enclosures` lives in `topology/plotting.py`.
   - `cartoon_zones` and `plot_resonance_zones` live in `loom/cartoon.py` and are importable from `tanglepack.loom`.
   - `grep -n "loom" src/tanglepack/topology/plotting.py src/tanglepack/topology/*.py` shows no import of loom.
   - `plotting.py` has grown by ≤ 85 lines.
   - The net line count is stated per group (src, tests, scripts) and is about −800.
9. **The cartoon, the zones plotter and every rewritten plot site draw byte-identical canvases.** A scratch script (on the model of `scratchpad/b9/probe/zones_fig.py`) draws the following on base, through the session methods, and on the branch, through the moved functions, and compares the sha256 of `fig.canvas.buffer_rgba()` plus an artist dump:
   - the circle cartoon on p3, nested and k10;
   - the line cartoon with every toggle off;
   - `plot_resonance_zones` on p3 and nested;
   - `plot_strong_pip`, `plot_holes`, `plot_pseudoneighbors` and `plot_stable_partition` per fixed point on the nested fixture;
   - `plot_dual_graph` on k10.

   Every pair must be equal.
10. **The 42 symbolic-itinerary PNGs are sha256-equal.**
    - Run in `git archive` exports of base and branch, **outside the repo**:
      - `henon_symbolic_itineraries.py`;
      - `_period3`;
      - `_nested`;
      - `_nested --heteroclinic`;
      - `_inversion`.
    - All of them run in full; each takes ≤ 14 s.
    - Compare base against branch, never against the committed PNGs (four are stale).
11. **The scripts behave the same.** Each of the 11 rewritten scripts is run in full in both exports:
    - Exit codes match.
    - stdout is identical.
    - Every PNG the script writes is sha256-equal base vs branch.
    - `henon_blast_period_3.py` fails on both with the same last traceback line (`:95`, "Intersection 12 is not a strong-pip candidate"). It is reported as still dead.
12. **No number changes on the corpus.** The 63-item probe (`scratchpad/b10/probe/run.sh`, re-pointed at the two worktrees) gives `all_base.jsonl` and `all_branch.jsonl` with 0 differences.
13. **Each commit is green.**
    - The branch has three commits, C1 → C2 → C3, in that order.
    - The suite passes at each commit (run the full suite at each, or the touched files plus `tests/facade tests/plotting` at C1 and C2 and the full suite at C3).
    - Each commit message names its concern.

## 8. Questions

Each has a default; "approved" accepts all of them.

- **Q0. Split B9 into B9a (this plan) and B9b (skeleton below)?** [**Yes: run B9a now and write B9b's request next.**]
- **Q1. Do thin plot delegates stay on the session?** [**None.** Delegates with docstrings would add about 230 lines and break the cap.]
- **Q2. Where do the cartoon's zones come from?** [**(a)**: `cartoon_zones` moves verbatim to `loom/cartoon.py` and reads the stored zones. No figure changes and `ResonanceZone.py` is untouched. (b) would compute them from pips, which changes `inversion_period2`'s circle and moves the boundary builder into topology.]
- **Q3. Where does `plot_resonance_zones` go?** [**`loom/cartoon.py`**, taking zones.]
- **Q4. Do the describe fan-outs go too?** [**Yes.** That is `describe_holes`, `describe_stable_partitions` and `_describe`, plus the test at `:68-83`.]
- **Q5. What happens to `test_a_pip_trimmed_away_is_not_replaced_by_plotting`?** [**Delete it.** `Trellis.plot_*` never touches the pip memory: `Trellis.py:622-637` classifies without choosing.]
- **Q6. Should "element on any side" become API?** [**No new API.** `owner_of_intersection` plus `stable_key_at` are enough. This belongs to B9b.]

## B9b skeleton: `one-element-at-end`

**Goal:** fold three implementations into `PartitionFamily.owner_at_end`, which reads the row through `row_of_end` and the branch through `stable_key_at` (D6, C7):
- `Trellis.element_for` (:1381) is deleted. After B9a its only readers are tests.
- `StablePartition._elements_at_bridge` (:1059) is deleted, along with:
  - the `elements_at_bridge` field (`TopologyResults.py:415-424`, `:442`);
  - its two writers (`StablePartition.py:536`, `PartitionFamily.py:898`);
  - the import at `PartitionFamily.py:140`.

  The field is written but never read in src.
- `BridgeClass._element_at_end` (:763, read at :597-598 and through `ends` at :740) becomes one call into a shared core.
  - `PartitionFamily.py:136` imports `BridgeClass`, so that core lives beside `row_of_end` in `StablePartition.py` (or behind a local import).
  - The planner decides where, and states the `AssertionError` → `ValueError` difference on an unowned crossing.

**Files:** `Trellis.py`, `StablePartition.py`, `PartitionFamily.py`, `BridgeClass.py`, `TopologyResults.py`, `tests/unit/topology/test_partition_elements.py`.

**Rewriting `test_partition_elements.py`:**
- The both-side checks (:45-78, :99-119, :122-161) are rewritten with `family.owner_of_intersection(id, side, trellis.stable_key_at(id, trellis.unstable_key_of(bid)))`, so they lose no strength.
- :77, :81-96 and :158 go with the field.
- `test_dual_walk.py:877` still binds `owner_at_end` onto a fake, so `owner_at_end` must keep calling `self.owner_of_intersection`.

**Proof:** the element-at-end agreement probe (`scratchpad/b9/probe/`, re-pointed from `partition_element_for` and `element_for` to the survivor) must give 7,090 / 14,180 / 1,776 agreements on the branch, and the 63-item probe 0 differences. The law tier reads none of the four, and it stays untouched.
