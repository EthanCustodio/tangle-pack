# Change: B4a `holes-algorithm-vs-diagnostics` (coder, iteration 1)

## Where

- **Base commit:** `ae1885947cfc682a08b730293471ebacab4464b3` (`ae18859`, the author's `main` HEAD). Diff everything against it.
- **Branch:** `refactor/holes-algorithm-vs-diagnostics`, in the worktree `/home/dezhu/Development/TanglePack/tangle-pack-holes-algorithm-vs-diagnostics`. Tip: `fc9f9a4`.
- **Base worktree:** `/home/dezhu/Development/TanglePack/tangle-pack-holes-algorithm-vs-diagnostics-base`, detached at `ae18859`.
- **Python:** `/home/dezhu/Development/TanglePack/tangle-pack/env/bin/python`. The editable install points at the main repo's `src`, but `pyproject.toml` sets `pythonpath = ["src", "tests"]`, so pytest run inside a worktree imports that worktree's code. Scripts run outside pytest need `PYTHONPATH=<worktree>/src:<worktree>/tests`.
- I did not touch the author's `main` or their uncommitted changes.

## Commits

| # | sha | what |
|---|---|---|
| 1 | `95ea780` | **Move.** The 16 names are cut out of `StablePartition.py` byte for byte, by AST span, and pasted in the same order into the new `topology/diagnostics.py`. That module has its own `logger` and imports only `row_of_end, rows_of_bridge` from `StablePartition`.<br>`check_bridge_rows_consistent` and `check_holes_share_bridge_side` stay in `StablePartition` and import `bridge_row_violation` / `bridge_side_violations` inside their bodies. `_punch_in_bridge` reaches `_bridge_side_of` through a lazy import.<br>`Trellis.punch_holes` takes `_warn_geometry_flags` from `.diagnostics` and keeps the call-time `StablePartition` imports, so the spies still work. `__init__` re-exports from both modules.<br>The test imports and logger names in files 6-11 follow. **Full suite green at this commit** (1171 passed, 38 skipped, 5 deselected, 1 xfailed, 211 s). |
| 2 | `728dcbb` | **Measured side on demand.** `Hole` drops `measured_side` and `side_source`, along with their docstring paragraphs. `_punch_in_bridge` no longer measures anything; `carried` still places the hole and `_map_backward` stays. `punch_holes` drops `side_source="signs"`.<br>`diagnostics.py` gains `side_source(hole)` and `measured_sides(trellis, holes)`, and `side_diagnoses(..., measured=None)`. `geometry_flags` is the one caller that builds `measured_sides(...)`. `_warn_geometry_flags` becomes `warn_geometry_flags`, and `Trellis` calls it.<br>The Dev Notes and docstrings no longer name the fields. The `SideSource` alias moves from `TopologyResults` to `diagnostics`, next to its only reader. **The suite is red at this commit alone**: 12 tests that build or read the removed fields fail until commit 3. This follows the plan's split (src in 2, tests in 3); the plan requires green only for commit 1. |
| 3 | `fc9f9a4` | **Tests.**<br>`laws.check_holes_are_classified` drops the `side_source` clause.<br>In `test_topology_invariants`, `_hole(origin, iterate, side, *, direct=False)` builds a direct hole with a `PseudoneighborPair` and a propagated one without. The three measured-side tests pass `measured={(origin, iterate): side}`. The `geometry_flags` test stubs `diagnostics.measured_sides`. The unmeasured-INFO test also asserts that the orbit check alone (no `measured`) logs nothing.<br>In `test_production_checks`, "carry their orbit side" reads `side_source()` and `measured_sides()`, with `assert measured` so it cannot pass vacuously. "Flags a contradicting measured side" monkeypatches `diagnostics.measured_sides` to flip one entry and expects exactly one WARNING on `tanglepack.topology.diagnostics`. |

## `git diff --stat ae18859`

```
 src/tanglepack/topology/StablePartition.py         | 657 +------------------
 src/tanglepack/topology/TopologyResults.py         |  19 +-
 src/tanglepack/topology/Trellis.py                 |   8 +-
 src/tanglepack/topology/__init__.py                |  18 +-
 src/tanglepack/topology/diagnostics.py             | 712 +++++++++++++++++++++
 tests/helpers/laws.py                              |   4 +-
 tests/plotting/test_plot_topology.py               |   2 +-
 .../regression/test_coarse_blasts_rows_by_signs.py |   4 +-
 tests/unit/topology/test_production_checks.py      |  29 +-
 tests/unit/topology/test_stable_partition.py       |   7 +-
 tests/unit/topology/test_topology_invariants.py    |  59 +-
 11 files changed, 815 insertions(+), 704 deletions(-)
```

## Suite result

`MPLBACKEND=Agg env/bin/python -m pytest -q -rxX`, run in the branch worktree on the code of the tip:

```
XFAIL tests/invariants/test_law_case_sanity.py::test_orientation_reversing_case_builds[orientation_reversing] - C3: INV-11: det J < 0 is out of scope (C3); FixedPoint.set_k_value rejects a saddle with exactly one negative eigenvalue
1171 passed, 38 skipped, 5 deselected, 1 xfailed in 201.29s (0:03:21)
```

- **Which tip.** The run was on `85eff04`. Afterwards I rebuilt commits 2-3 to correct one docstring Note in `measured_sides`, which wrongly claimed it was "the one place the topology calls the map" while `_map_backward` still places holes. The new tip `fc9f9a4` differs from `85eff04` by that 4-line docstring only.
- **Re-run on `fc9f9a4`.** `tests/unit/topology tests/plotting/test_plot_topology.py tests/regression/test_coarse_blasts_rows_by_signs.py tests/invariants`: 495 passed, 38 skipped, 1 xfailed.
- **The only xfail** is `cases.KNOWN_ISSUES[("orientation_reversing_case_builds", "orientation_reversing")]` (`tests/cases.py:1323`).
- **Commit 1 alone** has the same totals: 1171 / 38 / 5 / 1.

## Behaviour changes

There is no result change. Every hole, measured side, row diagnosis, flag and fingerprint is identical on base and branch (rubric 4 below). What does change:

1. **Logger name.** The geometry-flag WARNINGs and the "Hole-side check skips ..." INFO now come from `tanglepack.topology.diagnostics`, with the same text (plan section 5).
2. **Hand-built hole description.** A hand-built `Hole` with no `pair` now describes itself as `carried` in `SideDiagnosis.describe()`, where it used to say `None` (plan section 5).
3. **The I1 path no longer counts unmeasured holes.** `check_holes_share_bridge_side` → `side_diagnoses` without `measured` never counts "unmeasured", so it never logs that INFO for missing measurements, and it never calls the map. The INFO is still logged on the I1 path when holes are unplaced, with the unmeasured count 0.
   - On base the INFO could appear twice per `punch_holes(diagnose=True)` (I1 + flags). On the branch it appears at most once for unmeasured holes, on the flag path.
   - On every build probed no propagated hole was unmeasured, so this is not observable there.
4. **Not observed on any build: orbits started from a fallback.** `measured_sides` replays from each orbit's REFERENCE hole only, as the plan specifies. If an orbit's backward chain starts from the fallback (reference without a spanning bridge, `StablePartition.py` propagation fallback), its propagated holes have no entry. They are then counted "unmeasured" (INFO) where base measured them.
   - The analyst found the fallback never runs on any probed build, and B4b removes it.
5. **API.**
   - `Hole` no longer accepts `measured_side=` / `side_source=`.
   - `SideSource` is no longer importable from `TopologyResults`; it is in `diagnostics`. Nothing else imported it, and it is not in `__all__`.
   - `_warn_geometry_flags` is now `warn_geometry_flags`.

**Out-of-scope breakage, not fixed (outside the file list, as the plan's section 6 says):**
- `scripts/refactor_proof/hole-side-and-own-blast/hole_sides.py` imports moved names from `StablePartition` and will break.
- The docstring at `src/tanglepack/numerics/geometry.py:7` still names `topology.StablePartition._build_oriented_bridge_polyline`.

## Rubric, my reading

1. **Pass.** Full suite on the branch: 1171 passed, 38 skipped, 5 deselected, 1 xfailed. `-rxX` lists only the `orientation_reversing` KNOWN_ISSUE.
2. **Pass.** On base `ae18859`:
   - `test_topology_invariants`, `test_production_checks`, `test_stable_partition`, `test_plot_topology` and `test_coarse_blasts_rows_by_signs`: 94 passed.
   - `tests/invariants`, which exercises `helpers/laws.py`: 133 passed, 38 skipped, 1 xfailed.

   On the branch, the same files pass: the 5 files plus `test_exact_ownership` give 100 passed, and `tests/invariants` gives 133 / 38 / 1. Files 6-11 are edited on the branch (imports, helper, logger names), so "the same tests" means the branch versions on the branch and the base versions on base.
3. **Pass.** `git diff --stat ae18859..` lists exactly the eleven files of plan section 3 (above).
4. **Pass.**
   - **Script.** The comparison script is below. It was run in both worktrees over the 11 `cases.LAW_CASES` plus `p3_coarse_blasts`, on `session.topology().trellis` and on `session.topology([fp]).trellis` for each fixed point.
   - **Normalisation.** `FixedPoint` reprs in the fingerprint carry object addresses, so ` at 0x...>` is stripped with `sed` on both sides. The law-tier `fingerprint` is dumped minus its `products` entry, which holds `id()`s that are not comparable across processes.
   - **Result.** 1258 lines on each side and `diff` is EMPTY (md5 `74385b0439e6dec2619241cef8c9ca2c` both). The dump holds 435 hole lines, 382 measured sides and 240 row diagnoses.
   - **Non-empty where it can fail.** On `p3_coarse_blasts` it reports `geometry_flags: 2 rows, 20 sides` and 20 side diagnoses on both the topology and the per-fixed-point trellis, with 80 measured sides each. On base the measured sides come from `hole.measured_side`; on the branch they come from `measured_sides(trellis, trellis.holes)`.
   - **How `side_diagnoses` is called.** On the branch it gets `measured=measured_sides(...)`, so it diagnoses the same thing base did with the field.
   - **Per-case counts** (rows flagged, sides flagged, measured), one triple per trellis:
     - k10 (0,0,2)×2
     - k28_one_blast (0,0,2)×2; k28_two_blasts (0,0,3)×2
     - p3 (0,0,6)×2
     - nested and nested_inner_blasts (0,0,9),(0,0,3),(0,0,6)
     - inversion (0,0,4)×2; inversion_one_blast (0,0,6)×2; inversion_two_blasts (0,0,6)×2; inversion_period2 (0,0,8)×2
     - heteroclinic (0,0,69),(0,0,3),(0,0,40)
     - p3_coarse_blasts (2,20,80)×2

   This matches the analyst's §7-8 numbers.
5. **Pass.** At commit `95ea780`, for each of the 16 names, `ast.get_source_segment` in `diagnostics.py` equals its segment in `StablePartition.py` at `ae18859` (16/16 identical), and none of the 16 is still defined in `StablePartition.py`. `_punch_in_bridge`'s lazy `_bridge_side_of` import exists only in commit 1, as planned.
6. **Pass.** `wc -l src/tanglepack/topology/StablePartition.py` = 1354 (≤ 1370; base 1947).
7. **Pass**, reading the plan's rename as intended.
   - **Names defined.** `diagnostics.py` defines 18 top-level names: the 16 with `_warn_geometry_flags` renamed `warn_geometry_flags` (plan step 2), plus `side_source` and `measured_sides`. It also has the `SideSource` type alias moved from `TopologyResults`.
   - **Imports.** Its only import from `StablePartition` is the top-level `from .StablePartition import bridge_for_pair, row_of_end, rows_of_bridge`. `StablePartition.py` imports from `.diagnostics` only inside `check_bridge_rows_consistent` (`:713`) and `check_holes_share_bridge_side` (`:738`).
   - **Fresh interpreter.** `python -c "import tanglepack.topology.diagnostics"` with `PYTHONPATH=<branch>/src`, run from `/tmp`, works.
8. **Pass.** `grep -rn "measured_side\|side_source" src` hits only `topology/diagnostics.py` (17 lines). The same grep on `TopologyResults.py` is empty.
   - The `StablePartition.py` Dev Notes were reworded so they do not even contain those substrings ("the MEASURED side", "the measured-side flag"), because the grep is substring-based.
9. **Pass.** With `diagnostics.measured_sides` replaced by a function that raises, on `p3_coarse_blasts`:
   - `trellis.punch_holes(diagnose=False)` succeeds (84 holes);
   - `trellis.punch_holes(diagnose=True)` raises `RuntimeError("the map was called")`.

   This was run as a script, not added as a test; the test engineer may want to pin it.
10. **Pass.** `grep -n "map_inv\|dynamical_system" src/tanglepack/topology/StablePartition.py` hits only `:955` and `:958`, both inside `_map_backward` (`:951`), which still exists (default answer to question 1).
11. **Pass.** `sorted(tanglepack.topology.__all__)` is identical on base and branch (both printed, `uniq -c` = 2).
12. **Pass.** `test_coarse_blasts_rows_by_signs.py` passes with `_LOGGER = "tanglepack.topology.diagnostics"`. With that line reverted to `"tanglepack.topology.StablePartition"` it FAILS on the branch: `assert (0 - 0) > (2 * 0)` at `_warnings(caplog) - before > 2 * quiet`. The file was restored afterwards and the worktree is clean.
13. **Pass.** "its orbit expects", "on its image sub-arc" and "approaches its two defining crossings" each appear exactly once in `diagnostics.py`, verbatim (moved byte for byte; the message code was not edited).

## Small things beyond the letter of the plan (all inside the file list)

- **`SideSource` moved** from `TopologyResults` to `diagnostics`: its only reader is there, and the `Hole` docstring no longer refers to it.
- **`Trellis.py` docstrings.** The two references to `topology.StablePartition.geometry_flags` / `.bridge_side_violations` now say `topology.diagnostics.*`.
- **`test_coarse_blasts_rows_by_signs._warnings` docstring** says "diagnostics logger" instead of "partition logger".
- **Two test assertions added**, both only to keep existing tests non-vacuous: the unmeasured-INFO test also checks that no INFO is logged without `measured`, and the carried-side test asserts `measured_sides` is non-empty on p3.

## Rubric-4 comparison script (verbatim; run as `PYTHONPATH=<wt>/src:<wt>/tests MPLBACKEND=Agg python b4a_compare.py out.txt`, then `sed -E 's/ at 0x[0-9a-f]+>/>/g'` on both outputs and `diff`)

```python
"""B4a rubric 4: dump holes, measured sides, row diagnoses and flags, base or branch alike."""
from __future__ import annotations

import sys

from cases import BUILDERS, LAW_CASES
from helpers.law_tier import fingerprint
from tanglepack.topology import bridge_row_diagnosis, geometry_flags, side_diagnoses

try:  # branch: the measured side is computed on demand
    from tanglepack.topology.diagnostics import measured_sides
except ImportError:  # base: it is a Hole field
    measured_sides = None


def measured(trellis):
    if measured_sides is not None:
        return measured_sides(trellis, trellis.holes)
    return {(h.origin, h.iterate): h.measured_side for h in trellis.holes if h.measured_side is not None}


def end_fields(end):
    return (end.intersection_id, end.crossing_sign, end.combinatorial_row, end.geometric_row,
            end.approach_angle, end.displacement, end.unstable_cdist, end.stable_cdist)


def dump(label, trellis, out):
    holes = trellis.holes
    out.append(f"== {label}: {len(holes)} holes")
    for h in holes:
        out.append(repr(("hole", h.origin, h.iterate, h.bounding_ids, h.near_intersection_id,
                         h.bridge_side, [(i, w, r, k[1:]) for i, w, r, k in (h.openings or [])],
                         h.coords)))
    sides = measured(trellis)
    for key in sorted(sides, key=repr):
        out.append(repr(("measured", key, sides[key])))
    punched = {frozenset(h.bounding_ids) for h in holes if h.bounding_ids is not None}
    bridges = [b for b in trellis.bridges
               if frozenset((b.first_intersection, b.second_intersection)) in punched]
    for b in bridges:
        d = bridge_row_diagnosis(trellis, b)
        out.append(repr(("rows", b.id, None if d is None else (
            d.bridge_id, None if d.branch_key is None else d.branch_key[1:], d.n_nodes, d.chord,
            [end_fields(e) for e in d.ends], d.signs_alternate, d.geometry_disagrees, d.verdict))))
    rows, flagged = geometry_flags(trellis, bridges, holes)
    out.append(f"geometry_flags: {len(rows)} rows, {len(flagged)} sides")
    out.extend(d.describe() for d in rows)
    out.extend(d.describe() for d in flagged)
    kwargs = {"orientation_preserving": trellis.orientation_preserving}
    if measured_sides is not None:
        kwargs["measured"] = sides
    diagnoses = side_diagnoses(holes, **kwargs)
    out.append(f"side_diagnoses: {len(diagnoses)}")
    out.extend(d.describe() for d in diagnoses)
    return len(rows), len(flagged), len(sides)


def main(path):
    out = []
    for name in (*LAW_CASES, "p3_coarse_blasts"):
        case = BUILDERS[name]()
        fp = fingerprint(case)
        fp.pop("products")  # object ids, not comparable across processes
        out.append(f"#### {name}")
        out.append(repr(fp))
        counts = [dump("topology", case.session.topology().trellis, out)]
        for i, point in enumerate(case.fixed_points):
            counts.append(dump(f"per-fixed-point {i}", case.session.topology([point]).trellis, out))
        print(name, counts, flush=True)
    with open(path, "w") as fh:
        fh.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main(sys.argv[1])
```

**The diff (base vs branch, after address normalisation):** empty. `diff` exit status 0, and both files have md5 `74385b0439e6dec2619241cef8c9ca2c`, 1258 lines. Session copies are in the scratchpad (`b4a_compare.py`, `base.norm.txt`, `branch.norm.txt`, `rubric4.diff`); the script above is the complete source for re-running it.
