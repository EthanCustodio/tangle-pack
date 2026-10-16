# Plan: B0 `oracle-laws` (tests only)

*Planner, 2026-10-12, traversal 1 of `planner_human_2`. No analyst stage ran for this run, so the
findings are slice B0 of `docs/alignment_plan_2026-10-09.md` (with A.4, A.5, A.6 and A.8d) and section 5 of
`docs/code_fidelity_audit_2026-10-09.md`. I checked their claims against the code and ran probes on the
base commit `a3bd8f0` (scratch scripts only; no file changed). Three probe results change the B0 text, and
they are marked **[changed]** below.*

## 1. The answer

Yes. Four new laws check the definitions by computations that never call the code they check. On the base
commit, three of them already pass on all 11 law cases. Measured with sketches of checks 1 and 3:

- **Check 1 (the pip box).** No crossing lands inside the box on any case. Between 5 and 108 crossings are
  tested per case.
- **Check 2 (the pair intervals).** The orchestrator probed this on 2026-10-09: zero punctures.
- **Check 4 (k=10 is the full 2-shift).** On k=10 the build gives `a -> a u^-1 a^-1`, with matrix `[[2]]`
  and refined matrix `[[1,1],[1,1]]`.

**[changed]** Check 3, the hole containment, fails as B0 wrote it ("strictly between") on `k28_one_blast`
and `inversion`, and on `inversion_period2` at today's 9 stable steps.

- **What happens.** The reference pair's far end maps back exactly onto the container's own end crossing.
  The iterate table links them: `iterate(5,-1) = 3` on `k28_one_blast`, and `iterate(3,-2) = 1` on
  `inversion`. The real inverse map then overshoots that crossing by noise (relative 1e-7 to 7e-4).
- **What it means.** This is physics, not a bug. When a reference end's preimage is a crossing of N, it
  is an end of the container.
- **The fix.** The law is closed at such a shared end and strict elsewhere. With that change it passes on
  all 11 cases.
- **A wrong claim in the alignment plan.** A.0/A.6 say "the table never holds these points". That is
  false on three law cases. It does not matter here, but B4b should know.

**[changed]** Step 5 (`stable_steps=10`) passes every layer on `inversion_period2` except one. The law
`every_active_class_walks` now checks nothing there, because every active class lands on a singleton.
So it needs a not-applicable entry.

## 2. Recommendation

Make the change. It touches tests only, and the library behaves exactly as before.

## 3. Scope

- **Slug:** `oracle-laws`.
- **Branch and worktree:** branch `refactor/oracle-laws` in the worktree `../tangle-pack-oracle-laws`,
  made from the author's `HEAD` `a3bd8f0` on `main`.
  - The base is a detached worktree at `a3bd8f0`: `../tangle-pack-oracle-laws-base`.
  - The uncommitted documentation in the author's tree is never touched or committed: `.refactor/decisions.md`,
    `artifacts/request_origin_analyst.md`, and the two `docs/*_2026-10-09.md` files.

**File list.** Nothing else gets edited.

1. `tests/helpers/laws.py`: the four checks, registered in `LAYERS`.
2. `tests/cases.py`: the `NOT_APPLICABLE` entries, any `KNOWN_ISSUES` entries from check 2, and the
   `build_inversion_period2` default.
3. `tests/invariants/test_law_partition.py`: one sentence in the module docstring.
4. `tests/invariants/test_law_symbolic.py`: one sentence in the module docstring.
5. `tests/conftest.py`: one sentence in the docstring's no-pin paragraph recording the Q11 exception.
   This file is not in B0; see Question 1. Answer "no" to keep four files.

**What I verified for the file list:**

- **Registration.** Checks are NOT found by prefix. `LAYERS` in `laws.py` (line 2885) is an explicit
  dict, so a new `check_*` must be added to its layer's list there. Nothing else needs it:
  - `law_tier.layer_tests` reads `LAYERS`.
  - The layer test's own `__doc__` is built from each check's first docstring line.
  - The law id is the name minus `check_`.
- **Set-aside skips.** `NOT_APPLICABLE[(law_id, case)]` takes the pair out of `run_layer`. The pair then
  gets a `test_<layer>_set_aside` parameter that runs the check, asserts `count == 0`, and calls
  `pytest.skip("not applicable: …")`. So a not-applicable check must RETURN 0 on those cases, not
  assert. `KNOWN_ISSUES` gives `xfail(strict=True)` through `issue_marks`.
- **The two layer modules.** Their docstrings do not list laws by name. They describe every law of
  their layer in prose, and the closing sentences ("reference pairs are consecutive on both branches…",
  "the matrix counts the word tokens") play that role. So they stay in scope, one sentence each. They
  are documentation only: the test functions come from `layer_tests`.

**Size.** About 200 changed lines: four checks of about 35 lines each, about 20 lines in `cases.py`, and a few
docstring lines. One slice.

## 4. Steps

**Step 1: `check_chosen_pip_is_strong` (layer `partition`).** This is A.4's box test, re-derived from
the registry and `FixedPoint` alone. It imports nothing from `StrongPip`, uses no `image_cdist` and no
iterate table. It reads the whole registry, because classification un-chooses first.

```python
def check_chosen_pip_is_strong(case: "Case") -> int:
    """A.4: no crossing of the pip's fixed point, carried onto the pip's stable branch by the
    fewest backward steps, lies in the open box (0, c_s(q0)) x (0, c_u(q0)); q0's own orbit
    is skipped by a 2-D collision within SCALING_RTOL, the one identity assumption."""
    trellis, count = _staged(case), 0
    for fp in case.fixed_points:
        q0 = case.registry[trellis.pips[fp]]
        cycle = fp.branch_cycle("stable")
        s0, u0 = float(q0.stable_cdist), float(q0.unstable_cdist)
        for _rid, r in case.registry:
            if r.is_anchor or r.manifold_a_key[0] is not fp or r.manifold_b_key not in cycle:
                continue  # homoclinic crossings of fp only: a heteroclinic u-cdist is another saddle's
            d = (cycle.index(r.manifold_b_key) - cycle.index(q0.manifold_b_key)) % len(cycle)
            s = float(r.stable_cdist) / fp.per_step_beta("stable") ** d      # backward: stable grows
            u = float(r.unstable_cdist) / fp.per_step_beta("unstable") ** d  # unstable shrinks
            count += 1
            if abs(s - s0) <= SCALING_RTOL * s0 and abs(u - u0) <= SCALING_RTOL * u0:
                continue
            assert not (0 < s < s0 and 0 < u < u0), (_rid, s, s0, u, u0)
    return count
```

**Step 2: `check_reference_pairs_clear_of_X` (layer `partition`).** For each reference pair, test both
open intervals the same way, by scaling only.

1. **The unstable interval.** Use the cycle and β of the unstable side's fixed point (D7: A's).
2. **The stable interval.** Use the stable side's cycle and β.
3. **The crossings tested.** Every non-anchor crossing `r` of the cut N (`trellis.own_intersection_ids`)
   whose key of that stability lies on that cycle.
4. **The landings.** Every `n ≡ (pos(key) - pos(r_key)) mod k` whose scaled cdist `c(r)·β^n` falls
   strictly inside `(lo, hi)`. There are finitely many because `β^k > 1`: solve for the lap index
   between two logarithms.
5. **Allowed landings.** A landing the registry's iterate table links to `x0` or `x1` is allowed
   (`registry.iterate_table`, a lookup of numerics data, never `_unstable_interval_is_clear`). So is a
   scaled landing that collides with `x0` or `x1` on both cdists within `SCALING_RTOL`, where its other
   cdist is scaled the same way, on `advance_key`.
6. **Anything else fails the law.**

```python
    for pair in trellis.reference_pseudoneighbors:
        x0, x1 = (case.registry[i] for i in pair.as_tuple())
        for stability, key in (("unstable", x0.manifold_a_key), ("stable", x0.manifold_b_key)):
            fp = key[0]; cycle = fp.branch_cycle(stability); beta = fp.per_step_beta(stability)
            lo, hi = sorted((_cdist(x0, stability), _cdist(x1, stability)))
            for rid, r in crossings:
                for n in _laps(...):            # n ≡ n0 (mod k), lo < c(r) beta**n < hi
                    if case.registry.iterate_table... in pair.as_tuple(): continue
                    if _collides(...): continue
                    raise AssertionError((pair.as_tuple(), stability, rid, n))
            count += 1
```

- **Small helpers.** `_cdist` and the lap solver can be two-line local helpers; the step uses each
  twice.
- **Run it on base first.** The orchestrator's probe found no case failing. If the base run does find
  one, each failing case gets
  `KNOWN_ISSUES[("reference_pairs_clear_of_X", name)] = KnownIssue("NUM-18: …", wp="B2")`.
- **Correction to B0.** B0's Behaviour line ("two strict xfails, p3 and inversion") contradicts its own
  step 2 and the probe. The expected number of entries is zero.

**Step 3: `check_propagated_hole_inside_container` (layer `partition`).** **[changed]** The law is now
closed at a shared end. For each hole with `iterate = n < 0`:

- **The container.** It is `trellis.bridge_between(*hole.bounding_ids)`.
- **The reference bridge's ends.** They are `hole.origin`, and are never read through `bridge_for_pair`.
- **Carrying the ends back.** Map each end's `coords` back `|n|` times with the REAL
  `case.workbench.dynamical_system.map_inv`. This is the test side, so calling the map is allowed.
- **Projecting.** Project each carried point onto the container's polyline: its nodes and their
  `get_cdist("unstable")`, nearest segment, cdist interpolated along it.
- **The containment test.**
  - If the iterate table links the end to `n`, that linked crossing must be one of the container's two
    crossings. It is a shared end, and the real map's overshoot (1e-7 to 7e-4 relative) is noise.
  - Otherwise the projected cdist lies strictly inside the container's two crossing cdists.
- **What it never uses.** `image_cdist` and `_containing_bridge`.
- **Non-vacuous everywhere.** `propagated_holes_land_on_predicted_branch` already runs on every law case,
  so every case has a propagated hole.

**Step 4: `check_k10_is_the_full_shift` (layer `symbolic`; Q11).** The check returns 0 unless
`case.name == "k10"`, which is what the set-aside test needs. On k10 it is letter-free:

```python
    dyn = case.session.symbolic_dynamics()
    active = [cls for cls, cd in dyn.classes.items() if cd.kind == "active"]
    assert len(active) == 1
    tokens = [s for s in dyn.classes[active[0]].symbols if s.kind == "active"]
    assert [s.bridge_class for s in tokens] == active * 2
    assert sorted(s.direction for s in tokens) == [-1, 1]
    assert len({s.occurrence for s in tokens}) == 2              # two P1 pairs -> two children
    assert dyn.transition_matrix(refined=False)[1].tolist() == [[2]]
    assert dyn.transition_matrix()[1].tolist() == [[1, 1], [1, 1]]   # its in-splitting (A.8d)
```

In `cases.py`, add the entry below for the other ten law cases:

```python
**{("k10_is_the_full_shift", name): "the textbook answer is known for k=10 only"
   for name in LAW_CASES if name != "k10"},
```

**Step 5: `build_inversion_period2(stable_steps=10)` (Q12).**

1. Change the default.
2. Rewrite the docstring: it currently says "9 unstable and 9 stable steps", and "10 runs the stable
   manifold past the pip orbit" goes (moot since D10).
3. **[changed]** Add `"inversion_period2"` to the existing `every_active_class_walks` not-applicable
   tuple (`k28_one_blast`, `p3`, `inversion`). The probe shows it now checks nothing there.
4. The other four `inversion_period2` not-applicable entries stay: each still checks 0 in the probe.

Unaffected: `inversion_period2_10` and `test_pip_orbit_cut.py` pass both step counts explicitly, at
(10, 10).

**Step 6: docstrings.**

- `test_law_partition.py` gains: "the chosen pip passes the box test, every reference pair's two open
  intervals hold no scaled point of X, and every propagated hole's container holds the reference
  bridge's ends carried back by the real inverse map."
- `test_law_symbolic.py` gains: "on k=10, and only there, the dynamics is the full 2-shift (author
  Q11, the one textbook pin)."
- `tests/conftest.py` gains the same Q11 exception in its no-pin paragraph (Question 1).

## 5. Behaviour

None in the library: no `src/` file changes. The suite gains:

- four laws inside the existing layer tests;
- 11 visible set-aside skips: the k10 law on 10 cases, and `every_active_class_walks` on
  `inversion_period2`;
- possibly strict xfails from check 2 (expected none).

The law case `inversion_period2` becomes the 9-unstable / 10-stable build: 6 classes, 10 holes, class `e`
resolved, `is_reliable` True (probe).

## 6. Not doing

- **No deletions.** No tautological law is deleted; B12 revisits them.
- **No deep-series oracle.** Scaling-only is wrong at depth by design (A.5), so the laws run on
  `LAW_CASES` only.
- **No committed teeth tests.** The break-it checks below monkeypatch private library functions, which
  the test policy forbids pinning. The reviewer runs them as scratch scripts.
- **No library fix.** The A.0 / A.6 remark ("the table never holds these points") is not corrected
  here; it is noted for B4b.

## 7. Rubric

1. On the branch, `MPLBACKEND=Agg env/bin/python -m pytest -q -rxX` passes. Its xfails are exactly the
   base's (`orientation_reversing_case_builds`) plus check-2 entries, if any.
2. Every pinning test passes on the base commit and on the branch. Here the pinning tests are the four new
   laws, and the library under them is unchanged: `git diff a3bd8f0 refactor/oracle-laws -- src` is empty.
   So `pytest tests/invariants -q` is the same run on both. The reviewer also runs it in the base
   worktree with the branch's `tests/` copied in.
3. `git diff --stat a3bd8f0 refactor/oracle-laws` lists only the five files of section 3, or four if
   Question 1 is answered no.
4. No behaviour change beyond section 5. The only build that changes is `inversion_period2`; on the
   branch, `build_inversion_period2()` gives 6 classes, 10 holes and `is_reliable` True. The reviewer
   reads these numbers; no test pins them.
5. Independence: grep the bodies of the four new functions. Every one of `StrongPip`, `is_strong_pip`,
   `Pseudoneighbor`, `_unstable_interval_is_clear`, `StablePartition`, `_containing_bridge`,
   `bridge_for_pair` and `image_cdist` has zero hits. Check 1 also has no `iterate`. The module-level
   imports already in `laws.py` are not the new code.
6. Each new check is listed in its `LAYERS` entry: three under `"partition"` and one under
   `"symbolic"`. The `-v` ids of `test_partition_laws` / `test_symbolic_laws` are unchanged.
7. Teeth for check 1: patch `tanglepack.topology.StrongPip.is_strong_pip` to accept every crossing,
   build `k10`, `set_strong_pip` to a crossing the unpatched function rejects, then run `topology()`.
   `check_chosen_pip_is_strong` must raise. **[changed]** B0 said "the second candidate", which is
   really strong on k10, so the test would not bite. In the probe, crossings 5 and 6 are the rejected
   ones; pip 5 puts 2, 3 and 4 in the box.
8. Teeth for check 2: patch `Pseudoneighbor._unstable_interval_is_clear` to return True and build
   `k28_two_blasts`. `check_reference_pairs_clear_of_X` must raise. Not probed; if it does not bite,
   the coder reports it and does not weaken the law.
9. Teeth for check 3: patch `StablePartition._containing_bridge` to return the next bridge on the branch,
   and `StablePartition._map_backward` to return None (without it the build crashes in
   `_bridge_side_of`). Build `k10`; `check_propagated_hole_inside_container` must raise. Probed: all 4
   carried ends fall outside.
10. Check 3 with the strict-only rule (no shared-end closure) fails on `k28_one_blast` and `inversion` on
    base. With the closure it passes on all 11. The reviewer confirms both by toggling the clause.
11. Check 2 on base: the `("reference_pairs_clear_of_X", …)` entries in `KNOWN_ISSUES` are exactly the
    failing cases, each citing NUM-18 with `wp="B2"` (expected: none).
12. `check_k10_is_the_full_shift` passes on `k10`. `-rs` lists exactly 10 skips
    `not applicable: k10_is_the_full_shift on <case>`, plus one for
    `every_active_class_walks on inversion_period2`.
13. `build_inversion_period2`'s signature reads `stable_steps: int = 10`, and its docstring no longer
    contains "past the pip orbit".

## 8. Questions

1. **Should `tests/conftest.py` record the Q11 exception?** Its docstring says "No test pins a
   case-specific fact: no word, … transition matrix …", and check 4 is exactly such a pin.
   **Default: yes**, one sentence citing Q11, which makes it the fifth file. Answer "no" to keep B0's
   four files; then the policy text stays out of date until B12.
