# Findings: B7b `base-chords-in-stable-order` (analyst, iteration 1)

*(Saved verbatim by the orchestrator: the analyst's Write of this file was refused.)*

**Working tree is dirty.** The changes are uncommitted docs and run records: `.refactor/decisions.md`, `.refactor/knowledge/*.md`, `CLAUDE.md`, `docs/alignment_plan_2026-10-09.md`, `docs/code_fidelity_audit_2026-10-09.md`, a staged rename of the old `artifacts/request_origin_analyst.md`, untracked `.refactor/runs/*`, `artifacts/`, `figures/henon_period3_*` and two `scripts/henon_period3_*.py`. The change branches from `HEAD` = `a6a4f68`. None of these files is in the slice's file list.

**To the planner:** after you read this, your plan goes to the person running the session. They are asked directly and grade it "approved", "revise" or "answered". Two of the request's numbers are wrong, and the rubric has to be re-derived from section 2 below (see section 5).

## The answer, up front

- **The question.** Which stretches of a lobe image's own stable branch are empty? Measured on every lobe the corpus has, the stretches inside `f(H)` on the base branch are exactly the crossings paired in **STABLE order**. That is C19's Jordan rule, and it holds on the base branch as on a foreign one. The stop rule does not fire:
  - 47 lobe readings (stretched holes with 2 or more interior crossings on their base branch, over 16 builds and 47 series steps) all alternate out/in from the anchor, with the base counted as inside.
  - On every one of them, the covered gaps are the stable-order pairs.
  - Chain order is never the covered set where the two orders differ. *Verified* (section 3).
- **Where the orders differ.** Only one lobe in the corpus carries 4 or more base-branch interior crossings: hole `(6, 7)` on `k32_series`. It does so at **steps 2-10, not only 6-7**. The orders differ at all nine steps; T fails only at 6 and 7.
- **The B7 patch experiment was wrong.** Its scratch patch (`b7/patched_empty.py`) lost the line `empty.setdefault(b, []).append((a, kind))`, so every empty stretch was indexed only at its smaller id. That buggy patch reproduces B7's "split inert classes −1 on 15 builds and steps" exactly, list and all. *Verified:* I ran the buggy variant and got heteroclinic, p3_16, p3_6_blasts, nested inner×6, p3_15, p3_series 6-9, k32 2-5 and 8-10, each 1→0 or 3→2.
- **What the correct patch does.** The correct patch drops only the `branch_key != base_branch` test. It changes **no split count anywhere**. What it does change, all on `k32_series` only:
  - **Iterated elements** change at steps 2-10.
  - **`unmatched_members`** goes 1 → 0 at 6-7.
  - **`refined_rules`** change at 8-10: one INERT class loses a child (4 → 3), so one inert token's refined index changes. `M_1` does not change.
- **Rubric items 2 and 3 are false as written.** Item 2 ("`refined_rules` identical") and item 3 ("−1 split inert on the 15") must be re-stated (section 5).

## 1. The code at HEAD (verified by reading `a6a4f68`)

| Item | Request says | HEAD |
|---|---|---|
| `_empty_stretches` | `:927` | `:927` (Note `:946-956`; chord comment `:985-988`; called from `from_minimal` `:778`) |
| the base-branch test | `:993` | `:993` `if branch_key != base_branch:` (block `:985-996`) |
| `_incomplete_footprints` | `:1007` | `:1007` (called `:912`); `_hole_stretches` `:1070` |
| module docstring | `:30-46` | the bullets run `:31-47`; the chain/stable clause is `:39-41`, the odd-count clause `:43-47` |
| Dev Notes | `:88-92`, `:101-118` | History `:88-92`; Foreign chords `:101-118` (the B7b/known-issue sentences `:109-116`) |
| SymbolicDynamics docstring | `:28-31` | "laws on every build" `:27-28`; the k32/B7b sentence `:29-31` |
| KNOWN_ISSUES entry | `:1354-1361` | `:1354-1361` (dict `:1346`) |
| chain-order unit pin | `test_iterated_cut.py:92-110` | `:92-110`. The module docstring `:5-7` also states the chain-order rule, and must change too |
| sanity guard | — | `_STANDALONE_ISSUES` `test_law_case_sanity.py:36-43`; unknown-id check `:89-95`; wp regex `:100` `WP\d+|C\d|DEC-\d+|B\d+[a-z]?` |

**Stale example in the module docstring.** `PartitionFamily.py:49-58` uses pre-WP8 ids: "hole `(6, 5)` has base `(10, 7)` and chord `(8, 9)` … `L_1^2 = [3, 8]`". At HEAD the k28_two_blasts lobe is hole `(6,7)`, chain `8,9,10,11`, base `(8,11)`, chord `(9,10)`, as CLAUDE.md says. *Verified* by the probe. The slice edits this docstring, so the planner may fix the example or leave it. It is outside the request.

## 2. The patch run at HEAD (verified; 16 builds + 47 series steps = 63)

**Method.** `scratchpad/b7b/probe.py` takes HEAD's `_empty_stretches` source, replaces `if branch_key != base_branch:` with `if True:`, and execs that, so nothing else changes. I ran base and patch in separate processes.

- **Builds:** the 11 law cases plus `p3_16`, `p3_6_blasts`, `build_k28(blasts=8)` and `build_nested(inner_blasts=6)`. That is B7's 15. I added `p3_15` as a 16th.
- **Series steps** replicate `cases.blast_series`: blast, probe, then `_recomputed(workbench)`.
- **T and in-splitting** are measured by the test helpers themselves: `helpers.laws.p1_itinerary_disagreements` and `check_refined_matrix_is_an_in_splitting`.

**Identical on all 63** (base vs patch):
- letter-free `rules` (class label, direction)
- `M` (with class labels) and `M_1` (names and entries)
- class count, active count and inert count
- the number of split ACTIVE and split INERT classes
- `is_reliable`

**Theorem T.**
- HEAD: 1 disagreement on k32 at step 6 and at step 7 (class `k`, `A:p1@0.0/R#7 <-> A:p1@0.0/R#9`, member `(6,7)` reads `v(L_(0.0;1)^4, L_(0.0;3)^2) l^-1(…, R_(0.0;19)^4)`; siblings `(54,55)` at step 6 and `(54,55)`,`(64,65)`,`(96,97)` at step 7 read `L_(0.0;1)^5, L_(0.0;3)^1` and `R_(0.0;19)^3`).
- Patch: **0 everywhere**.
- `p1_tested` is equal under both. k32 reads `[0,0,0,2,3,4,9,13,13,13,13]`.

**In-splitting:** "ok" on all 63, under both.

**What differs, and only on these steps:**

| step | iterated elements | other |
|---|---|---|
| k32_series 2-5 | 23 → 24. Right row: the crossing 9 is pinched to a singleton `[9,9]` between the hole `(10,9)` and the new chord `(9,12)`, so `R_(0.0;12)^2..^4` are renamed | none |
| k32_series 6-7 | 40 → 39. The singleton `L_(0.0;1)^4 = [10,10]` goes; `(94,9)` and `(9,12)` re-bounded on both rows | `unmatched_members` 1 → 0 (member `(6,7)`) |
| k32_series 8-10 | 113 → 111. The left row loses the singletons `[10,10]` and `[9,9]` | `refined_rules` differ: inert class `v = {L_(0.0;3), L_(0.0;5)}` (positional `L#2<->L#4`) has children `v_1..v_4` at HEAD and `v_1..v_3` patched. Member `(9,10)` (alone in `v_3 = (L^5, L_5^2)` at HEAD) joins `(93,94)`… in `v_3 = (L^5, L_5^1)`. One active word (class `R#23<->R#27`) has its token `v_4` become `v_3`. `refined_rules` key count 55 → 54 |

There are no other differences on any build or step. In particular nothing changes on heteroclinic, p3_15/16, p3_6_blasts, nested inner×6 or p3_series 6-9, contrary to the request's list.

**Pytest under the patch.** I ran the patch as a scratch pytest plugin (`-p b7b_patch_plugin`, `PYTHONPATH=scratchpad/b7b`) over `test_iterated_cut.py`, `test_foreign_chords.py`, `test_partition_family.py`, `test_iterated_partition.py` and `test_class_count_under_blasting.py`. The only failure is `test_members_share_the_P1_itinerary_under_blasting[k32_series]`, which goes XPASS(strict). At HEAD the same selection gives 43 passed and the two expected xfails. *Verified.*

## 3. Base-branch lobes and the stop rule (verified)

**Method.** For every hole bridge with a stretch, I took its `_image_chain` (consecutive), its base branch and its interior crossings on that branch.

- **The polygon.** Where there are 2 or more such crossings, I built the lobe image polygon as `_lobe_image` does: each chain arc from its `Bridge` (all arcs had one), closed by the base along the stable manifold.
- **The gaps.** I walked the base branch from the anchor (root) to the workbench manifold's tail. The gaps are the stretches between consecutive chain crossings on the branch.
- **The reading.** Each gap is closed at its crossing coordinates, and its arc-length midpoint is read with `point_in_polygon(tol=0)`. The base gap is not read, because it lies on the boundary.
- **Robustness.** The midpoint's distance to the ring, over the gap length, is at least 0.0497 everywhere (inversion_period2) and about 0.46 typically.

**Counts:**
- **Readings with 2 or more base-branch interior crossings: 47.** 38 have 2; 9 have 4. No count is odd, and none is 6 or more.
- **Builds/steps with a 4+ base branch: `k32_series` steps 2, 3, 4, 5, 6, 7, 8, 9, 10.** Every time it is the same lobe:
  - hole `(6,7)`, chain `8→9→10→11→12→13`, base `(8,13)`;
  - interior stable order `11, 10, 9, 12`;
  - cdists 13: 0.036221, 8: 0.176023, 11: 6.571306, 10: 6.664764, 9: 6.785874, 12: 7.417581, tail 7.575983.
- **Readings on that lobe (identical at steps 2-10):**
  - (anchor, 13) out
  - (13, 8) base
  - (8, 11) out
  - **(11, 10) IN**
  - (10, 9) out
  - **(9, 12) IN**
  - (12, tail) out
- **Covered set** = `{(10,11), (9,12)}` = the STABLE pairs. The chain pairs `(9,10), (11,12)` are not covered: `(10,9)` lies outside and `(11,12)` spans three gaps. Verdict "stable".
- **The other 38 readings** have one chord, so the two orders coincide. Each alternates, and the one covered gap is the chord. On every one of them that chord is already a hole stretch (kind `"hole"`).
- **Lobes starting AT the anchor** (inversion `(0,13)`, inversion_period2 `(0,5)`, inversion_series@0) read base, out, in, out. That is still alternation: the anchor is a corner of `f(H)`, on its boundary, not inside it. A test that insists the first gap is "out" would wrongly fail there; it must count the base as inside.

**Stop rule: not triggered.** No lobe is covered in chain order and not in stable order, and every lobe alternates.

**One clause of step 2 is false as written.** "no gap outside is [an empty stretch]" fails, on the branch as well as on base:
- inversion `(3,2)`, inversion_period2 `(15,14)` and inversion_series@0 `(3,2)` are outside their lobe but are another hole's stretch;
- the same holds for k32 steps 2-5 at `(10,9)`.

The statement that holds on all 47 is that every inside gap is empty and no outside gap is a base or chord of THIS lobe.

**Which form fails on base.**
- The inside-gap clause fails on base at k32 steps 2-7: `(9,12)` is inside and not empty under chain order.
- At steps 8-10 `(9,12)` is itself a hole stretch, so on base only the "not this lobe's chord" clause fails there: chain order's `(9,10)` is outside.
- k32 step 6 fails both, so it satisfies the request's "FAILS on base on k32 step 6".

**Parametrisation.** A representative set covering both kinds is k32_series step 6 (4 crossings, orders differ) and any 2-crossing case. Law cases with a 2-crossing base lobe: k10, k28_two_blasts, nested, nested_inner_blasts, inversion, inversion_one_blast, inversion_two_blasts, inversion_period2 and heteroclinic. The anchor-cornered inversion lobe is the shape most worth including. Builds with none: k28_one_blast, p3, p3_15, p3_16 and p3_6_blasts (their lobes have 0 base-branch interior crossings). These are the NOT_APPLICABLE candidates.

**Midpoints in `test_foreign_chords.py`.** It takes `stretch[len//2]` of `arc_polyline`, which holds real nodes only, and skips a node-less stretch. On k32 the gaps `(11,10)` (2 nodes) and `(79,80)` / `(81,78)` (0 nodes) show that a base-branch reading must close each gap at its crossings, as my probe does.

## 4. Folding `_incomplete_footprints`; the monkeypatch sites; the guard (verified by reading and by running)

**Callers outside the module: six sites in four files** (the request says five):
1. `test_iterated_cut.py:85`: the harness. It returns the dict, which `:106`, `:120` and `:129` index.
2. `test_iterated_cut.py:158-170`: the WRAPPER `with_foreign_far_ends`, which calls the original and rebuilds a dict from `.items()`.
3. `test_iterated_partition.py:70`: a direct call that indexes `empty[cut.intersection_id]`. This one is not in the request's list.
4. `test_iterated_partition.py:90-100`: the WRAPPER `without_bases`, which rebuilds from `.items()`.
5. `test_foreign_chords.py:60`: `.get`.
6. `test_partition_family.py:169-170`: calls both helpers.

**If `_empty_stretches` returns a pair or a result object,** all six must change:
- The two wrappers must pass the second value through.
- The four direct calls must unpack.
- `from_minimal` must stop calling `_incomplete_footprints` at `:912`.

**If `_empty_stretches` keeps its dict and both helpers become thin wrappers** over one private chain walk, the six sites need nothing. But `from_minimal` must keep calling `cls._empty_stretches`, because the two wrappers monkeypatch exactly that name and expect `from_minimal` to read it. A dict subclass carrying the footprints would be lost through both wrappers, which build plain dicts. That is harmless to those two tests, but fragile.

**The odd-count meaning.** Folding does not change it. `_incomplete_footprints` counts per branch, so it does not depend on order. Under stable order, though, the unpaired crossing of an odd count is the OUTERMOST along the branch, not "its last crossing" in chain order. That fits "the chain left the computed stable manifold there" better. Two wordings mention "last": the module docstring `:44-45` and the `_empty_stretches` Note. The odd case occurs nowhere on the corpus.

**Re-pinning the unit tests.** The HEAD pin `test_foreign_crossings_pair_in_stable_order` and `test_the_base_is_read_at_a_non_anchor_end` (`:112-123`) give every OWN crossing `stable_cdist=1.0`. A stable sort keeps chain order on ties, so both pass under the patch: **the current pin has no teeth for B7b**. The k32 shape built by hand (hole `(1,2)`; cdists 8: 0.176, 9: 6.786, 10: 6.665, 11: 6.571, 12: 7.418, 13: 0.036; chain `8..13`) gives these chords. *Verified:* `scratchpad/b7b/handbuilt.py`.
- HEAD: `(9,10), (11,12)`.
- Patched: `(10,11), (9,12)`.

`test_partition_family.py:150` (`stable_cdist=float(iid)`) is indifferent to the order and passes under both.

**The guard.** `_STANDALONE_ISSUES` is an allow-list, and four of its six ids already have no KNOWN_ISSUES entry. Dropping `series_members_share_the_P1_itinerary` from it is safe, and keeping it is equally safe. Either way `:89-95` still rejects any `(law, case)` pair whose id is in neither the layered laws nor the allow-list. *Verified:* the guard test passes at HEAD, and the logic is a pure set check.

## 5. What differs from the request

1. **"−1 split inert on 15 builds and steps" is an artefact.** It comes from B7's buggy scratch patch (section 2). The correct patch changes **no** split count. The only class-structure change is one inert child fewer (k32 8-10).
2. **"letter-free `rules` and `refined_rules` identical" is false for `refined_rules` on k32 8-10.** One inert token's index changes and one inert child vanishes. `rules`, `M`, `M_1`, counts and `is_reliable` are identical. Rubric 2 should compare `refined_rules` restricted to active refined symbols, or should state the k32 8-10 inert change.
3. **"the one place the two orders differ … k32_series at 6 and 7 blasts."** It is k32_series **2-10**: one lobe, nine steps. Iterated names and intervals change on all nine. T fails only at 6-7. Rubric 3's "names differ only on the measured list" holds with that list.
4. **"pinches crossing 10 into a singleton" is verified** at steps 6-7 (`L_(0.0;1)^4 = [10,10]`) and at 8-10. At steps 2-5 it is the PATCH that pinches 9 (between hole `(10,9)` and chord `(9,12)`).
5. **The step-2 clause "no gap outside is [empty]" is false on the branch too** (section 3). The clause should read "no outside gap is this lobe's base or chord".
6. **The base-commit failure.** On base, the geometric test fails on k32 steps 2-7 through an inside gap, and at 8-10 only through the outside-chord clause. Step 6 fails both ways.
7. **Call sites.** There are six `_empty_stretches` call sites, not five (`test_iterated_partition.py:70` is missing from the list). `test_iterated_cut.py`'s module docstring `:5-7` also states the chain-order rule.
8. **Line numbers:** module bullets `:31-47` (not `:30-46`); SymbolicDynamics `:27-31`. The rest match.
9. **New behaviour change not in the request:** `unmatched_members` 1 → 0 on k32 6-7 (member `(6,7)` now matches its child).

## Recommendation

The smallest change that answers the question is the request's own one-line change: delete the `if branch_key != base_branch:` condition at `PartitionFamily.py:993` and sort on every branch. The evidence above supports it on every lobe the corpus has, and the stop rule is not triggered. The planner owns the rest (docstrings, tests, the fold); I do not re-plan it.

## Out of scope (noticed)

- The stale k=2.8 ids in the `PartitionFamily` module docstring `:49-58` (section 1).
- The CLAUDE.md "Partition families" bullet still states chain order. It is B12, and the tree's CLAUDE.md carries uncommitted edits anyway.
- `test_foreign_chords.py` skips node-less stretches (`if not len(stretch): continue`), so a foreign gap with no node is never read.

## Scratch (rerun, never trust)

Everything is in `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/b7b/`:
- `probe.py` (modes `base` / `patch` / `buggy`), `run_all.sh`
- `compare.py`, `lobes_summary.py`, `dig.py`
- `b7b_patch_plugin.py` (pytest `-p`), `handbuilt.py`
- raw `out/*.jsonl`

The knowledge files are updated: `.refactor/knowledge/topology.md` (B7's inert-split claim corrected, a new B7b analyst section) and `.refactor/knowledge/testing.md` (the toothless pins, the guard, the midpoints).
