# Plan — B8 `unresolved-sets`

## The two problems, plainly

1. **The iterate table breaks ties silently.** Sometimes two crossings both look like the preimage of the same crossing `j`. The map is a bijection, so at most one of them is right, and the brackets cannot tell which. Today the code links whichever came first, logs a WARNING and moves on. From then on the table holds a link that was decided by order, not by the definition. A second, milder case is one source whose brackets meet several images. That source is already left unlinked, but nothing records that it was ambiguous.
2. **A stale pip cannot be forgotten.** `resonance_zone(id)` adopts the zone's pip as the chosen pip without checking that it is a strong-pip candidate. Once a chosen pip stops being a candidate, `topology()` raises, and no API exists to drop it. Deleting it from the session dict is not enough either, because the cached trellis still holds it and `topology()` picks it up again.

The fix records both as data instead of guessing. `registry.unresolved_links` holds the links the brackets cannot decide. `TangleSession.forget_strong_pip(fp)` drops a pip, and `resonance_zone` stops adopting a pip that is not a candidate.

## 1. The answer

Both sites are where alignment-plan section B8 puts them, a few lines lower at `HEAD` `8a1cdc3`:
- the collision: `IterateInference.py:128-136`;
- the multi-match: `:232-240`;
- the pip adoption: `TangleSession.py:1509`;
- the raise: `Topology.py:121-125`.

The analyst ran a probe after B1 over every builder, every step of all five series and k28 x8, each with a recompute:
- **No collision fires anywhere.** An order-free census agrees: no image is claimed twice, even with the unresolved detections counted.
- **The multi-match fires only on `near_tangency`**, which no law reads.
- **Every recipe's zone pip is a candidate**, before and after its zone.

So on the law cases the collision rule is a safety net the definition keeps honest, and the visible change is to the pip API. The plan's example collision on k28 x8 (sources 125 and 251 claiming 213) uses pre-B1 ids and no longer exists.

## 2. Recommendation

Change it. Keep slice B8 as decided, with the analyst's corrections:
- three more files;
- an order-free collision rule;
- `forget_strong_pip` also resets the trellis;
- the candidacy test happens after the zone is defined;
- the generation bump compares before and after;
- fixtures chosen by rule.

## 3. Scope

Slug **`unresolved-sets`**. Branch `refactor/unresolved-sets`, worktree `../tangle-pack-unresolved-sets` from `HEAD` `8a1cdc3`. The uncommitted `docs/` and `.refactor/` edits and the staged rename stay untouched.

Files (eight):
1. `src/tanglepack/numerics/IterateTable.py`: a `remove` method. The table has none, and its docstring forbids touching `_forward`/`_backward` from outside.
2. `src/tanglepack/numerics/IntersectionRegistry.py`: `unresolved_links` and `remove_iterate`.
3. `src/tanglepack/numerics/IterateInference.py`: recording, and the order-free collision rule.
4. `src/tanglepack/loom/TangleSession.py`: `forget_strong_pip`, the candidacy test in `resonance_zone`, and the warning text.
5. `src/tanglepack/loom/Topology.py`: the `_choose_pip` error text and docstring.
6. `tests/helpers/laws.py`: the step-5 law, inside `check_iterate_table_complete`.
7. `tests/unit/numerics/test_iterate_inference.py`
8. `tests/facade/test_session_topology.py`

Estimate: about 50 changed source lines and about 120 test lines.

## 4. Steps

**Step 1. The registry gets a set and a removal.**
```python
# IterateTable
def remove(self, source_id: int, n: int) -> None:
    """Forget f^n(source) and its backward twin."""
    target = self._forward[source_id].pop(n)
    del self._backward[target][n]

# IntersectionRegistry
self.unresolved_links: dict[int, tuple[int, ...]] = {}   # source -> candidate images, none linked

def remove_iterate(self, source_id: int, n: int) -> None:
    self.iterate_table.remove(source_id, n)          # drops the graph edge iter_{source}_{n} too
    ...; self._generation += 1
```
`compute_intersections` builds a fresh registry, so the set starts empty after every recompute.

**Step 2. Multi-match.** `_bracketed_image` records `registry.unresolved_links[src] = tuple(images)` and keeps its INFO line.

**Step 3. Collision, decided by no order.** At the start of the pass, `infer_iterates` takes a snapshot `before = dict(registry.unresolved_links)` and drops the entries of its own sources. This happens before the `if not detected` early return. After gating, it groups the matches by image:
```python
claimants = defaultdict(list)
for i, j in _map_gated(workbench, registry, matches).items():
    claimants[j].append(i)
for j, sources in claimants.items():
    held = table[j, -1]
    others = [s for s, c in registry.unresolved_links.items() if c == (j,)]  # contested in an earlier pass
    contested = [*sources, *([held] if held is not None else []), *others]
    if len(contested) == 1:
        registry.register_iterate(i := sources[0], 1, j)
        recorded += 1
    else:
        logger.warning("image %d is claimed by crossings %s; none linked", j, contested)
        if held is not None:
            registry.remove_iterate(held, 1)
        registry.unresolved_links.update(dict.fromkeys(contested, (j,)))
if registry.unresolved_links != before:
    registry._generation += 1
```
Grouping removes the third-claimant hole. The `others` term means a partial pass, such as a bare `iterate_bridge`, cannot relink a source that an earlier pass left contested.

**Step 4. The pip.**
- `TangleSession.forget_strong_pip(fp)`: `self._strong_pips.pop(fp, None)` and `self.trellis(fp).strong_pip = None`. The second line restores the uncut trellis, so the next `topology()` classifies and picks the default.
- `resonance_zone`: after the zone is defined, run the same test `_restore_strong_pip` uses: `t = self.trellis(rz.key); found = t.registry.find(pip); found is not None and t._is_candidate(found)`. If it passes, remember the pip. Otherwise log one WARNING; the zone is still stored and the previous choice stays.
- `Topology._choose_pip`'s `ValueError` and `_restore_strong_pip`'s WARNING keep the phrase "no longer a strong-pip candidate", which existing tests match. Both end with "choose another with set_strong_pip() or drop it with forget_strong_pip()".

**Step 5. The law** (layer `crossings`, inside `check_iterate_table_complete`):
- skip any source in `registry.unresolved_links`, the same way a linked one is skipped;
- assert that no unresolved source has a `+1` link;
- assert that every entry has at least one candidate;
- assert that no entry `(j,)` coexists with a `(j, -1)` link.

## 5. Behaviour

**None on any law case**: no collision or multi-match fires there, and no recipe adopts a non-candidate pip. Outside the law cases:
- `near_tangency` records about 25 `unresolved_links` entries, which nothing reads.
- If a collision ever fires, neither claimant is linked (before, the first was).
- `resonance_zone()` at a non-candidate logs a WARNING and keeps the previous pip. Before, it adopted the pip and `topology()` later raised.
- The new `forget_strong_pip` exists.

## 6. Not doing

- Odd foreign-chord counts and I2 on image bridges (`PartitionFamily.py:760, 960-967`): that is B7. Marker disagreement was B6a.
- **Legitimate candidates that make `topology()` assert** (equal signs on a bridge joined across the pip cut: k10 2, k28 x1/x2 4, inversion 2 and 14). These are reachable through the public `set_strong_pip`. This is a real defect for a later run.
- **A zone at a non-candidate inward of the chosen pip** trims the pip away, and `topology()` still raises. `forget_strong_pip`, named in the message, is the remedy.
- `check_iterate_table_complete` scans the whole registry instead of N. That is harmless while the two are equal on the law cases.
- `CLAUDE.md` wording: that is B12.

## 7. Rubric

1. `MPLBACKEND=Agg env/bin/python -m pytest -q` passes on the branch, and `-rxX` lists only `cases.KNOWN_ISSUES`.
2. Every new test passes on the branch and fails on base `8a1cdc3`. Every existing test touching these files passes on both.
3. `git diff --stat 8a1cdc3` lists exactly the eight files of section 3.
4. **Fingerprints match:** for every case in `cases.LAW_CASES`, `tests/helpers/law_tier.py::fingerprint(case)` on the base worktree equals the one on the branch. Diff the two dumps.
5. **The set is empty where it should be:** on every `LAW_CASE` and every step of every `cases.SERIES`, `registry.unresolved_links == {}` on the branch, and the `+1` link count equals the base's.
6. **The old warning is gone:** `grep -n "already has the preimage" src` finds nothing. `grep -n "none linked" src/tanglepack/numerics/IterateInference.py` finds exactly one WARNING, and its branch calls `remove_iterate` and writes `unresolved_links`. No `continue` sits in that loop.
7. **`remove` exists:** `IterateTable.remove` and `IntersectionRegistry.remove_iterate` exist. After `register_iterate(a, 1, b)` then `remove_iterate(a, 1)`, the table holds neither `(a, 1)` nor `(b, -1)`, and the graph has no edge `iter_{a}_1`.
8. **Unit test, collision across passes,** on the module's `henon_tangle_with_bridges` fixture, chosen by rule:
   - take a linked pair `x -> y` and another detected crossing `z`;
   - empty the table and hand-link `z -> y`;
   - run `infer_iterates(workbench, [x])`.

   Expected: neither `x` nor `z` is linked to `y`, both map to `(y,)` in `unresolved_links`, and one WARNING names `y`. On base, `z` keeps `y`.
9. **Unit test, collision within one pass:** monkeypatch `_bracketed_image` and `_map_gated` so that three sources name one image. None is linked and all three are recorded. A version of step 3 that walks the matches in order would link the third.
10. **Unit test, multi-match:** a monkeypatched multi-match records the source with its candidates. Running `infer_iterates` twice leaves `registry.generation` unchanged on the second pass.
11. **Facade test, non-candidate zone:** on `inversion_two_blasts`, pick by rule a non-candidate whose stable branch differs from the chosen pip's. Then `resonance_zone(id)` must:
    - leave `chosen_strong_pip(fp)` unchanged;
    - log exactly one WARNING;
    - let `topology()` run without raising.

    On base, `topology()` raises `ValueError`.
12. **Facade test, forget:** on k10, `set_strong_pip(fp, c)` with `c` any candidate other than the default, then `forget_strong_pip(fp)`, then `topology()`. Expected: `chosen_strong_pip(fp)` and `trellis(fp).strong_pip` are the default pip. No topology runs at `c`, so candidates that assert are not an issue.
13. **The error names the remedy:** a stale pip (the existing `test_session_bridge_classes.py` stale-pip set-up) raises a `ValueError` matching both "no longer a strong-pip candidate" and "forget_strong_pip".
14. **The law rejects bad data:** `check_iterate_table_complete` raises on a hand-built registry where an unresolved source `(j,)` coexists with a `(j, -1)` link. Shown by a one-off probe in the review, not a committed test.

## 8. Questions

1. Under `add_resonance_zones`, `resonance_zone` runs with `recompute=False`, so the candidacy test of step 4 reads the trellis before the batch's single recompute. The analyst found every recipe's pip passes there, including the nested build. **Default: test there anyway, and keep `add_resonance_zones` unchanged.** The alternative is to defer the test to after the batch recompute, which costs about 6 more lines in `add_resonance_zones`.
