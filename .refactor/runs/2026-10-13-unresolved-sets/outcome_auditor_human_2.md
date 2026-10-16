# Outcome — B8 `unresolved-sets` (auditor, 2026-10-13)

Branch `refactor/unresolved-sets` @ `ae18859`, base `8a1cdc3`, worktree `../tangle-pack-unresolved-sets`. The review passed in round 1. Not merged: the merge is yours.

Your main tree is dirty, as it was before this run: uncommitted `docs/` and `.refactor/` edits, the staged rename of `artifacts/request_origin_analyst.md`, and untracked run folders. This run added one entry to `.refactor/decisions.md` and edits to three knowledge files. The branch touches none of these.

## The two problems, and what happens to them now

1. **The iterate table broke ties silently.** Sometimes two crossings both looked like the preimage of the same crossing `j`. The code linked whichever came first and logged a WARNING. The map is a bijection, so at most one of them can be right, and nothing in the definition picks one. Now neither is linked. Both go into a new set, `registry.unresolved_links`, as `{source: (j,)}`. If `j` already had a link, that link is removed.
2. **A stale pip could not be forgotten.** `resonance_zone(id)` adopted the zone's pip without checking it, and a pip that stopped being a candidate made `topology()` raise for good. Deleting it from the session was not enough, because the cached trellis re-supplied it. Now `TangleSession.forget_strong_pip(fp)` drops it and uncuts the trellis, and the next `topology()` takes the default. `resonance_zone` keeps a pip only if it is a candidate on the trimmed trellis. Otherwise it logs one WARNING, stores the zone anyway and keeps the previous choice.

**The answer to your question now:** at neither site does the code decide something the definition leaves open and then carry on. What the brackets cannot decide is stored as data:
- a source whose brackets meet several images is recorded with all of them;
- an image claimed more than once is recorded with every claimant.

The law `iterate_table_complete` keeps that set apart from the linked sources. The pip path no longer adopts a non-candidate, and the error a stale pip raises now names the fix. On real data none of this fires on any law case or series step, so the table rule is a safety net. The pip API is the part you can see.

## The most important hunk (`src/tanglepack/numerics/IterateInference.py`, `infer_iterates`)

Before. The first claimant wins, decided by the order of passes and of the dict:
```python
for i, j in _map_gated(workbench, registry, matches).items():
    if (j, -1) in table:
        logger.warning("crossing %d: its image %d already has the preimage %d; not linked", ...)
        continue
    registry.register_iterate(i, 1, j)
    recorded += 1
```
After. The claimants are grouped by image before anything is linked:
```python
claimed_by: dict[int, list[int]] = defaultdict(list)
for i, j in _map_gated(workbench, registry, matches).items():
    claimed_by[j].append(i)
for j, claimants in claimed_by.items():
    held = table[j, -1]
    earlier = [s for s, images in registry.unresolved_links.items() if images == (j,)]
    contested = [*claimants, *([] if held is None else [held]), *earlier]
    if len(contested) == 1:
        registry.register_iterate(contested[0], 1, j)
        recorded += 1
    else:
        logger.warning("image %d is claimed by crossings %s; none linked", j, contested)
        if held is not None:
            registry.remove_iterate(held, 1)
        registry.unresolved_links.update(dict.fromkeys(contested, (j,)))
if registry.unresolved_links != before:
    registry._generation += 1
```
How the rest of the pass works:
- Each pass first drops its own sources' entries, so it decides them afresh.
- The `earlier` term stops a later partial pass, such as one `iterate_bridge`, from relinking an image an earlier pass left contested.
- The old early `return` when nothing is detected is gone, so the generation check at the end runs on every path.

## What got shorter or clearer

Nothing got shorter. Source is +107/−28 lines over five files:

| file | lines |
|---|---|
| `IterateInference.py` | +43/−22 |
| `TangleSession.py` | +31/−4 |
| `IntersectionRegistry.py` | +20 |
| `IterateTable.py` | +9 |
| `Topology.py` | +4/−2 |

Tests are +194/−3. What is clearer:
- The collision loop now says the rule plainly ("contested by more than one: link none").
- `IterateTable` has a real `remove` instead of having no way to undo a link.
- The `infer_iterates` docstring states what is recorded and why.

## Every behaviour change

- **Law cases and series: no state change.** Verified by the test engineer (`te/pin.py`) and re-read by the reviewer. On all 11 `LAW_CASES` these are identical on base and branch:
  - crossings, bridge ids, holes, pips, partition signature, alphabet, registry generation;
  - the full `+1` link set.

  The `+1` link sets are also identical at every step of all five `SERIES`, both after the blast and after the recompute. `unresolved_links == {}` everywhere.
- **`workbench.generation` reads higher** on the seven cases that define a zone: 2→3 on k28 ×1/×2 and inversion ×1/×2, and 2→4 on nested, nested_inner_blasts and heteroclinic. The candidacy test builds one more trellis, and this counter counts reads that see a change. The test engineer reproduced the same bump on base with nothing but an added trellis read. No state differs.
- **`near_tangency` now records about 25 `unresolved_links` entries**, for example `25, 27 -> (61, 59)`. Nothing reads them, and it is not a law case.
- **A collision, if one ever fired, links no claimant.** Before, the first one was linked.
- **`resonance_zone` at a non-candidate** now logs a WARNING and keeps the previous pip, and `topology()` runs. Before, it adopted the pip and `topology()` raised `ValueError`. Shown on inversion ×2, crossing 33.
- **New API:** `forget_strong_pip(fp)`. The stale-pip WARNING and `ValueError` now end with "…or drop it with forget_strong_pip()".
- **Two departures inside `resonance_zone`** (`TangleSession.py:1520-1536`), both inside the eight files and both needed:
  - (a) **It pops the previous pip before building the trellis**, and puts it back only on the non-candidate path. Taken literally, the plan's step would have rebuilt the trellis with the old pip restored. Where the zone trims the old pip away, that logs a false "no longer a candidate" WARNING: 1 WARNING on the inversion ×2 build, against 0 on base and 0 on the branch. With no previous pip, the session is left without one, so the next `topology()` takes the default.
  - (b) **It sets `trellis.strong_pip = found` on the candidate path.** Otherwise `Topology._choose_pip` would re-adopt whatever the cached trellis already held.

## What was verified, and how

- **Suite, rerun by me on the branch:** `MPLBACKEND=Agg env/bin/python -m pytest -q -rxX` gives 1171 passed, 38 skipped, 5 deselected, 1 xfailed (198 s). The only XFAIL is `orientation_reversing` (C3), the one `KNOWN_ISSUES` entry. Base has 1162; the difference is the 9 new tests.
- **Scope, checked by me:** `git diff --stat 8a1cdc3` lists exactly the eight planned files, and the worktree is clean.
- **Diff read in full by me.** I checked that `contested` cannot hold duplicates:
  - `held` is linked, so it is not among this pass's pending sources;
  - this pass's own entries are popped before `earlier` is read;
  - a linked `held` cannot be an `earlier` entry, which the law guarantees.

  A full pass after a collision re-records the same pair without bumping the generation. `_map_gated({})` returns early, so dropping the early return costs no map call.
- **Add-zones path, probed by me** (`scratchpad/aud/batch.py`), plan question 1. On nested and heteroclinic, which use `add_resonance_zones` and so test candidacy before the batch recompute, the branch matches base exactly:
  - chosen pips A=4 and B=7 on both session and trellis;
  - no loom WARNING;
  - `unresolved_links` empty.
- **New tests fail on base** (test engineer): base 9 failed / 11 passed, branch 20 passed. 15 of 17 break-its turn a test red. Break 2 is an equivalent mutant: an in-order loop that keeps the live `earlier` term reaches the same result. The law's three new clauses were shown to reject bad data by a one-off probe only.
- **Not verified by anyone:** a real collision on a built case, because none exists in the corpus.

## What did not go to plan

- The two `resonance_zone` departures above; the reviewer judged both correct.
- The coder renamed the plan's `i := sources[0]`, which would have shadowed the `sources` parameter.
- The literal rubric 4 (fingerprints equal) fails on the raw observation counter only. It cannot hold together with the step-4 candidacy test you approved, and the reviewer passed it on state.
- The test engineer added 3 tests and strengthened 1. Without them, three breaks went undetected: dropping the `earlier` term, dropping the per-pass pop, and dropping the generation bump on recording alone.

## Points I am least sure were done the way you meant

1. **Under `add_resonance_zones` the candidacy test runs before the batch's single recompute.** It reads trimmed manifolds against the not-yet-recomputed registry. This was your accepted default, and every recipe's pip passes there; I checked nested and heteroclinic myself. But a pip that is a candidate only after the recompute, or only before it, would be judged on the wrong trellis.
2. **The law's new clauses are guarded only by a review probe.** The set is empty on every law case, so the committed law callers pass even with each clause deleted. If you want them to have teeth in the suite, that would take a unit test feeding the law a hand-built registry. The plan chose not to add one.
3. **Candidates that make `topology()` assert are left for later:** k10 2, k28 ×1/×2 4, and inversion 2 and 14. `forget_strong_pip` does not help there, because the pip is a legitimate candidate and the failure is the I2 row assertion on a bridge joined across the cut. They are reachable through the public `set_strong_pip`, and this is the next real defect on the pip path.

Also noted, not graded:
- `IterateInference` writes `registry._generation += 1` into another module's private state (the plan asked for it; a public bump method would be cleaner).
- `check_iterate_table_complete` still scans the whole registry rather than N.

## Files

- Branch worktree: `/home/dezhu/Development/TanglePack/tangle-pack-unresolved-sets`. The key changes:
  - `src/tanglepack/numerics/IterateInference.py` (~:118-156)
  - `src/tanglepack/loom/TangleSession.py` (:1119 `forget_strong_pip`; :1520-1536 `resonance_zone`)
  - `tests/helpers/laws.py` (~:770)
- Decisions entry: `/home/dezhu/Development/TanglePack/tangle-pack/.refactor/decisions.md` (2026-10-13 B8 outcome).
- Knowledge files updated:
  - `.refactor/knowledge/numerics.md` (`unresolved_links`, the collision rule, `remove`)
  - `.refactor/knowledge/loom.md` (pip memory, `forget_strong_pip`, the zone candidacy test)
  - `.refactor/knowledge/testing.md` (new tests, law clauses, the fingerprint caveat)
