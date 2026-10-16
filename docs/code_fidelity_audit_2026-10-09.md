# How true is tanglepack to its mathematics? An audit (2026-10-09)

An independent reading of the library against the three algorithm notes in `docs/`
(strong pip, pseudoneighbors, stable-manifold partition) and `docs/pipeline_chalkboard.md`.
Audited state: `main` at `69371fd` plus the uncommitted D15 working tree (forward-iterate
holes removed). Method: every layer was read file by file with probes run against the
case builders in `tests/cases.py`; the full suite was run once (1122 passed, 26 skipped,
1 xfail, 5 perf deselected, 156 s). Line numbers refer to that working tree. `CLAUDE.md`,
the knowledge files and `open_issues.md` were treated as claims to verify, not as truth.

## 1. Verdict

**Stages 0-6 (manifolds through the stable partition) are faithful to the written
mathematics, with two definitional gaps and one wrong heuristic. Stages 7-9 (classes,
words, transition matrix) have no written mathematics to be faithful to, and that is
where the complexity and the piecemeal rules live.**

1. The strong-pip test and the partition are the PDF algorithms. The pseudoneighbor test
   checks only half the definition (the unstable interval against all iterates, the stable
   interval against registered crossings only). The backward hole chain is β-scaling on
   every call in practice, while the documentation says "table first".
2. One numerics heuristic is wrong by the project's own law: a crossing-collapse rule runs
   only on the blast path, so a plain recompute after deep blasting changes the crossing
   set. Reproduced here at k=2.8, 8 blasts: 248 crossings become 252.
3. The word of a class is a dual-graph search with a sound filter (the marker rule) and
   heuristic tie-breaks (shortest, walls, detour cap, judgement budget, first in order).
   Measured on every law case and three deep runs, it reproduces the registered image chain
   wherever one exists. What it adds, net of that, is words for classes no computed bridge
   spans. Every refined child inherits the parent's word verbatim; that is correct iff all
   members of a class share one iterated-level itinerary, a theorem the code never stated
   or checked (revised, see stage 9).
4. The code is three to five times the size of a faithful implementation. 31.8k library
   lines; the session facade (1.7k) contains no mathematics; stage 5+6 needs ~150 lines
   and has ~1,950; 303 bug-ticket ids are cited inside `src/`.
5. Tests are honest about not pinning arbitrary facts, but about a third of the 98 laws
   can only fail if the code reports its own failure, no law checks a surviving
   pseudoneighbor pair or the chosen pip against the PDF definition on a real build, and
   the one case with a textbook answer (k=10, the full 2-shift) is never checked end to end.

## 2. The mathematics, as written

Three one-page notes in `docs/`:

- **Strong pip.** q0 is a strong pip iff no crossing, mapped back onto q0's stable branch
  by cdist scaling, lands in the open box (0, c_s(q0)) x (0, c_u(q0)).
- **Pseudoneighbor pair.** Consecutive crossings x, x' are a pair iff the open stable
  interval AND the open unstable interval between them contain no point of X, the set of
  all iterates of all crossings. Reference pairs are found in W^S(r_n, r_{n+p}).
- **Partition.** Each reference bridge is mapped backward until it maps into a bridge
  already visited; every landing bridge gets a hole; intervals bounding a hole are open.

Everything after that (bridge classes, inertness, the minimal trellis, the iterated
partition, the dual graph, the word, refinement, the transition matrix) exists only in the
chalkboard and `CLAUDE.md`, written after the code as a record of decisions D1-D15 and
C1-C22. **Half the pipeline has no specification independent of its implementation.**
That is the first and largest finding, because it is what makes the later rules
unfalsifiable.

## 3. Fidelity, stage by stage

### Stage 0, manifold growth

Faithful in shape: fundamental segment `[f^-k(x0), x0]`, mapped with refinement, cdist
monotone by construction (`ManifoldMachine.py:355` sets `cdist = stretch * cdist`).

Two deviations from "a polyline within `area_cutoff` of the true curve":

- A refined point is `f(midpoint of the two preiterates' coordinates)`, cdist the mean of
  its neighbours (`ManifoldMachine.py:759-792`). Sub-refinements are chords of chords, so
  the polyline is `f(chord)` and the error compounds by the local stretch each lap. The
  area test bounds the smoothness of the image polyline, not its distance to W^u. No
  ground truth exists for this in the repo.
- Refinement re-walks the whole merged curve every step (`:386`), which is why every
  bridge memo is invalidated on growth.

Heuristics, none derived from the math: `_MIN_SEED_STEP = 5e-6` (the "dynamic" seed is
~1e-5 everywhere, so it is a constant), a 1e-8 spacing guard, the merge tie rule "keep the
iterate-linked point" (17 drops per k=2.8 build), `area_cutoff = 1e-7` (chosen empirically),
the 2M point budget.

The data model (two doubly-linked lists inside `Point`, `BranchPoint` slots, walk caches,
about 1 KB per point per `budget.py:21`) serves four uses that are each a per-node
`(preimage_xy, mapped_flag)`. `edist` and the stretch-param slots are written and never
read. Arrays would do.

### Stage 1, crossings

One rtree per stability, unstable segments queried against the stable index
(`Tangle.py:526-538`), so "no same-stability pair" is structural. Clean.

The anchor is a synthetic `Intersection` at (0,0) threaded through 13 guard sites in
numerics plus special rules in the registry, `create_bridges` and `_register_bridge`. The
D4 model (one point on every branch) is right; the implementation is a special case, not a
node of the curves.

**Wrong.** `Tangle._collapse_noise_crossings` (`Tangle.py:569-627`) collapses a run of
near-coincident detections to its median if the run is odd and drops it entirely if even.
It is reached only through `populate_intersections_for_manifold`, which only
`BridgeIterator` calls (`BridgeIterator.py:363`). `compute_intersections` never runs it.
Reproduced in this session on k=2.8 at 8 blasts:

```
build:            crossings=248  +1 links=145  bridges=175   (2 collapses logged)
plain recompute:  crossings=252  +1 links=144  bridges=179
                  + 3 "multi-match (unlinked)", 2 "image already has preimage"
```

Chalkboard stage 10 ("grow, blast and recompute commute", NUM-12/C20) is false exactly
where the heuristic fires. `tests/regression/test_growth_after_a_blast.py` covers only
k28x2 and nested, where it never fires. The rule also votes away u×s detections by parity,
which is not the author's artifact law (that law is about same-stability pairs).

Follow-up probes (planning round, same day) sharpened this. The four gained detections are
the outer members of two ODD runs of three, `(75, 1, 2)` and `(125, 250, 251)`, spaced
1e-13 to 1e-12 relative along the unstable curve; the build kept each median. Their
`resolution_ratio` is 624 to 4642, so the resolution flag cannot see them, and the kept
median 125 (ratio 12) is the least resolved of its run. The heuristic is load-bearing:
built with it disabled, k28x8 has 40 holes and 23 classes instead of 24 and 15, with
different words, and both results report reliable. It fired on no law case and on no other
blast series. The pipeline has no rule for a registered under-resolved crossing; the
alignment plan's Q13 proposes one (`docs/alignment_plan_2026-10-09.md`, A.2).

### Stage 2, iterate table

The bracket rule (`IterateInference.py:199-277`) is a sound necessary condition plus
uniqueness plus a map gate. "Exact" overclaims (cdist is linearly interpolated between
nodes), but no fallback registers a guessed link. The partial bijection is asserted at
registration and **pre-empted by a warning** one level up: "image already has a preimage"
logs and the first source in registry order wins (`:125-132`). Order-dependent; it fired in
the probe above.

`FixedPoint` is genuinely the one rule for branch bookkeeping (`advance_key`,
`per_step_beta`, `branch_cycle`); topology reads no eigenvalue. Two local duplicates
remain (`TangleWorkbench._stretch_per_round`, `DualGraph.py:1171`).

### Stage 3, strong pip

`StrongPip.is_strong_pip` is the PDF algorithm: every crossing of the fixed point is
mapped to its minimal representative on q0's branch (table first, else scaling) and tested
against the open box. One representative suffices because the box lies inside the computed
trellis, so any orbit point in the box is itself registered. One tolerance: an unlinked
representative within `collision_rtol = 1e-2` of q0 on both cdists is taken to be q0's own
orbit. With the exact table this path rarely fires, but it is a 1% identity assumption.

The pip-trellis cut (D10/D12, `Trellis.cut_at_pip_orbit`, `Trellis.py:256-309`) is a real
idea and the right place to decide N; the code is clean. But its depth m is "the last
forward iterate of q0 the table holds", and that depth is set by blast policy:
`Blast.min_separation` (`Blast.py:397-412`) drops a child "too close" and never feeds it
again, which is what stalls the pip orbit (`open_issues.md:473`: f^10(q0) dropped at k32
blast 9). The zone-interior filter (`Blast.py:271-277`) is why BCG-8's reference bridge is
never iterated. **A resolution guard and a zone polygon decide the definition of N.** D12's
"class count never decreases" is therefore a measured law, not a theorem.

### Stage 4, pseudoneighbors

`compute_pseudoneighbors` walks the single reference window from the pip toward
M^k(pip) and keeps a consecutive pair iff `_unstable_interval_is_clear`
(`Pseudoneighbor.py:329-470`). That function does implement X over all n in Z, placing
unlinked iterates by scaling and recognising the pair's own orbit by a 2-D collision at
`collision_rtol = 1e-2`. Faithful for the unstable interval.

**The stable interval is never tested against X.** Consecutiveness on the computed branch
only excludes registered crossings. NUM-18 documents orphan window pairs on p3 and inversion
whose open stable interval holds the predicted images of tip crossings; they "pass only
because the stable-side check reads registered crossings alone". The author recorded this
as a known asymmetry with no change wanted (memory, 2026-09-15). It remains the clearest
case of a definition implemented by half, and the reference switches at depth (BCG-4,
BCG-8) have not been re-examined under the full definition.

Other deltas: pairs touching the anchor are never pairs (`lo <= 0` returns False). D15 (today,
uncommitted) removed the provisional `+1..+(k-1)` forward-hole rule, the right direction
(net 184 lines deleted, p3 9 to 5 classes, the four lost were inert).

### Stage 5, holes

`propagate_reference_holes` (`StablePartition.py:397-546`) is the PDF chain: move both
bridge ends back by `image_cdist`, bisect the branch's crossing order for the container,
stop on a visited bridge. Faithful. But the documented "lookup first, scaling as fallback"
is a fiction here: a reference crossing's preimage lies beyond the trimmed stable branch
and is never registered, so `from_table` was False on 100% of calls on k10, k28, inversion
and heteroclinic and on most of nested and p3, and the flag is discarded at `:513` with no
log. Containment is decided to scaling accuracy, silently. This is also the mechanism
behind the NUM-5 drops at coarse cutoff.

The sign arithmetic is correct. `row_of_end` (`:784-835`) and `_hole_openings`
(`:1333-1378`) follow from sign = sgn(u+ × s+) and the bridge leaving its first end along u+;
the carried side (`side_after`) is sound and does not depend on INT-2. The measured side
disagreed zero times on every law case. The anchor-branch rule (`:1356-1368`) is an inline
second copy of `Trellis.stable_key_at` plus an inversion flip stated nowhere else.

Two paths punch one hole (direct backward pairs, then propagation) with two dedups; no
case produces a direct backward pair.

### Stage 6, partition

`_hole_marks_on_branch` + `intervals_from_marks` (`:1657-1793`): cut at hole-bridge
endpoints, open at a hole, singletons falling out of "both flanks open". About 75 lines,
faithful, clean. The bare-branch element (C1) and the synthetic anchor are small honest
additions.

### Stage 7, bridge classes and inertness

The class as the unordered element pair by rows, forward folding of loops, orphans: as the
chalkboard says (`BridgeClass.py:582-722, 852-873`).

Inertness (`:826-849`) is "inert iff it has image evidence and every image is a loop or
inert; no evidence keeps it active". That is monotone in evidence, not a property of the
dynamics: growth can flip a class either way (the Dev Notes say so: "grow until the images
are registered"). It is load-bearing: forcing every class active adds 1-2 image bridges and
3-6 iterated elements on six of eleven cases, because inert hole bridges' chains are long
runs of loops whose crossings would cut the partition. Whether dropping those cuts is right
is a modelling question the record does not answer; the chalkboard gives the definition but
not the reason inert images may be omitted from the minimal trellis.

### Stage 8, iterated partition, dual graph, word

**Iterated partition.** `IteratedHomotopyPartition.from_minimal` (`PartitionFamily.py:668-889`)
opens, at each image-bridge crossing that is an END of an empty stretch (a hole, a lobe
base, a chord), the flank facing it. It is an endpoint rule, not a set. A set-theoretic
candidate ("a stretch is empty if a hole or a hole's image covers it") reproduces the code
on all eleven shallow cases and on k28 at 4 and 6 blasts, and diverges at depth (k28 8
blasts, p3 6 blasts, inversion2, heteroclinic), always by cutting more. The Dev Notes
record three rewrites against fixture expectations ("R_5^3 must hold the a_1 bridge"). The
odd foreign-chord count (C19, `:960-967`) is a parity violation of the Jordan argument
treated as a WARNING. **No mathematical object is named that this stage computes.**

**Arrangement and dual graph.** The rotation-from-sign arrangement (`Arrangement.py`) is a
correct combinatorial-map construction and sparse mode is elegant. Geometry is used for
the outer face per component and one containment test each. Clean; the complexity is not
here.

**The word.** Two readings exist. (a) `trellis_itinerary` reads the registered image chain
of a member bridge, which numerics already computed. (b) `shortest_walks` searches the
dual graph between landing spans under the marker rule, walls, `MAX_DETOUR = 4` and a
judgement budget `max_walks = 64`, and D3 declares (b) the answer and (a) evidence.

The marker rule is a sound necessary condition (f is a homeomorphism, so f(D_b) holds
exactly the images of D_b's markers). It is not a complete invariant of the homotopy class:
two non-homotopic arcs in a punctured disk can enclose the same punctures. So (b) is a
search with a sound filter and heuristic tie-breaks. Measured:

| Probe | Result |
|---|---|
| Walk vs registered chain, every case where a chain gives a verdict | 100% agreement, except inversion at 6 blasts (BCG-8, already flagged unreliable) |
| Share of classes already read off the trellis path | p3 4/5, nested 4/9, heteroclinic 23/75, p3 6 blasts 20/43, k28 8 blasts 9/23 |
| Remove the marker rule | 2 classes change, both on heteroclinic, both with chains that already give the answer |
| Remove walls | nothing changes; one tie reopens (`bn`, the member-less connecting class) |
| Reverse seed order | nothing changes |
| `max_walks = 8` | nothing changes |

Net of (a), the walk machinery decides the classes no computed bridge spans (D11) and the
few spanned classes with no registered image. Fifteen distinct rules are needed to predict
a walked word; one is needed for the trellis reading.

### Stage 9, refinement and the transition matrix

`refine` splits a class by occurrence; `match_members` is ~90 clean lines. Every
refined child receives the parent's word verbatim (`SymbolicDynamics.py:762`,
`for child in children: refined_rules[child.name] = list(refined_word)`); nothing computes
f(a_i) separately.

Revised after the planning round (same day). The inheritance is correct exactly when every
member of a class has the same itinerary at the iterated level (the alignment plan's
"Theorem T", with a planar argument: a hole-bridge image separating two members' images
would pull back to a hole inside the hole-free strip between them). Under T the refined
graph is the in-splitting of the class graph: identical sibling rows, children's columns
summing to the parent's, conjugate and of the same entropy, and k=10's `[[1,1],[1,1]]` is
derived from `[[2]]`. T holds on every class where it can be tested today, which is only 8
classes across the corpus (k10 1, inversion2 2, heteroclinic 2, k28 8 blasts 3); no other
build has two members with registered chains. So the defect is not the arithmetic but that
T was unstated and unchecked, and the matrix had not been derived. The plan's slice B7 states
T and adds the law.

### Stage 10, recompute and the session

Zones and trims feed no topology output (probed: un-trimming and clearing zones on k28x2
left holes, classes, letters and words identical). Pure-function holds for the
mathematics. The pip choice is session state with four writers; `resonance_zone()` adopts
its pip without a candidacy test (`TangleSession.py:1509` vs `:1113`), and a pip that stops
being a candidate makes `topology()` raise with no API to forget it. Letters are
history-dependent labels by design.

## 4. Complexity

| Measure | Value |
|---|---|
| Library lines | 31,810 (plotting 3,207; graphviz 358; gpu 183; numerics diagnostics 437) |
| Test lines | 24,937 (`helpers/laws.py` 3,043; 98 laws) |
| Functions, numerics / topology / loom | 307 / 537 / 103 |
| `logger.warning` + `info`, numerics / topology / loom | 20 / 77 / 8 |
| Named tolerances and constants (non-plotting) | 17; none derived from the math, ~8 argued from numerics, the rest tuned (18.0, 3, 1e-7, 5e-6, `MAX_DETOUR` 4) |
| Decision/WP ids cited in `src/` | 303 occurrences; 128 distinct issue ids in `open_issues.md` |
| Commits | 368 total; 246 since 2026-10-01 |
| Churn, `StablePartition.py` | 35 commits, 3,725 added / 1,766 deleted, 1,947 now |
| Churn, `TangleWorkbench.py` | 59 commits, 5,665 / 3,422, 2,239 now |
| Churn, `DualGraph.py` | 20 commits, 2,918 / 1,468, 1,450 now |

**Algorithm versus everything else, per file** (code lines, docstrings excluded):

| File | code | algorithm | diagnostics, placement, data, plumbing |
|---|---|---|---|
| `StablePartition.py` | ~990 | ~410 | ~340 geometric diagnostics, ~170 plot placement |
| `PartitionFamily.py` | ~740 | ~325 | ~370 accessors, `describe`, `Cut` records |
| `BridgeClass.py` | ~500 | ~235 | ~260 dataclasses, table, `describe` |
| `DualWalk.py` + `SymbolicDynamics.py` | ~1,700 | ~? | 8 walk statuses, 3-valued logic in three places, 21 warn-and-continue sites |
| `TangleSession.py` | 1,735 total | 0 | caches, pip memory, 243 plot delegates, 182 cartoon geometry, 264 of fan-outs that "never feed a product" |
| `ResonanceZone.py` | 574 total | ~50 (`ring_cycles`) | a trim that is now a view, stale `intersection_ids`, cartoon stitching |
| `TangleWorkbench.py` | 2,239 total | ~300 | 617 of nine growth drivers, 181 resolution reporting, 149 plotting, id epochs |
| `IntersectionRegistry.py` | 860 total | ~400 | ~230 networkx graph and queries with zero topology/loom callers |

Given a `Trellis` with signs and `image_cdist`, the PDF's stages 5+6 need about 120-150
lines. The file has 1,947.

Dead or test-only API: `_get_refined_point`/`_get_iterate` ("kept as the pinned reference"
for one unit test), `grow_x_times`, `refine_manifold(final_node)`, `release_manifold`,
`preimage_bridges`, `trim_stable_at_intersection` (exported, zero users), `classify_bridges`,
`partition_element_for`, `terminated_early`, `StableNode.other_face`, `Arrangement.preimage_of`,
`euler_characteristic`, `DualGraph.graph()`. `topology.__all__` has 101 names, about half
plotting, legend handles, diagnosis dataclasses and sort keys. Four implementations of
"element at a bridge end" exist (`Trellis.element_for`, `BridgeClass._element_at_end`,
`PartitionFamily.owner_at_end`, `StablePartition._elements_at_bridge`).

**A faithful implementation, sized.**

- Numerics, ~2.5k lines: per branch, arrays `xy, c, pre_xy, mapped`; a growth step maps the
  unmapped rows, scales c, merges and refines in vectorised layers (the same kernel as
  today); one rtree per stability; a registry of `(u_key, c_u, s_key, c_s, sign, brackets)`;
  the bracket rule with collisions and multi-matches surfaced as an unresolved set; one
  recompute path and no collapse heuristic. Keep `FixedPoint` and `geometry.py`.
- Stages 3-6, ~400 lines: the three PDF algorithms over the cut trellis, the pseudoneighbor
  test applied to both intervals, signs deciding rows and openings, one WARNING when a
  scaled span lies within `SCALING_RTOL` of a crossing, no geometry in the algorithm module.
- Stage 7, ~150 lines, once inertness is either a definition or dropped.
- Stages 8-9: cannot be sized until the word and the iterated partition are defined
  independently of the search. If the word is the registered image chain of one member,
  about 1,500 lines of walk machinery go.
- Loom, ~400 lines: `workbench + Topology(generation, selection, pips) + BridgeAlphabet`,
  with `set_strong_pip`, `forget_strong_pip`, and `blast(fp, n)` = iterate every un-iterated
  bridge whose midpoint lies in the fixed point's pip-orbit lobe. Zones computed from
  `Trellis.pips` on demand; cartoon geometry in `plotting`.

## 5. Cheats, hacks and questionable tests

### In the library

| # | Finding | Where | Severity |
|---|---|---|---|
| 1 | Crossing collapse only on the blast path; recompute does not commute at depth; the only rule the pipeline has for under-resolution, and it decides k28x8's topology | `Tangle.py:569-627`, `BridgeIterator.py:363` | wrong |
| 2 | Refined children inherit the parent's word; correct iff Theorem T, which was unstated and unchecked (8 classes testable today) | `SymbolicDynamics.py:762` | smell (was: wrong); see stage 9 |
| 3 | `min_separation` and the zone polygon set the depth of N | `Blast.py:271-277, 397-412` | wrong by the chalkboard's own definition |
| 4 | Stable interval of a pseudoneighbor pair never tested against X | `Pseudoneighbor.py:149-230` | definitional gap, author-accepted |
| 5 | Backward hole chain is scaling on every call; "table first" is documentation only; flag discarded | `StablePartition.py:513` | smell now, risk at depth |
| 6 | First-source-wins on an iterate-table collision, warn and continue | `IterateInference.py:125-132` | smell |
| 7 | `is_reliable` ignores `truncated`; `max_walks=1` yields "unique" and reliable on heteroclinic | `SymbolicDynamics.py:839-848` | latent wrong |
| 8 | Row invariant I2 asserted for hole bridges, warned for image bridges | `PartitionFamily.py:760` vs `StablePartition.py:1073-1088` | smell |
| 9 | Zone's pip adopted without a candidacy test; lost pip has no recovery | `TangleSession.py:1509, 1082`, `Topology.py:122-126` | trap |
| 10 | Member disagreement on markers: warn and read the first | `DualWalk.py:817-828` | smell |
| 11 | Missing crossing sign: "assuming the positive rotation" | `Arrangement.py:288` | wrong if hit, never hit |
| 12 | `cdist_tol = 1e-6` absolute decides tiling of an image span (10% near the anchor, meaningless at depth); 11 topology/loom sites still read it | `BridgeIterator.py:182-189` | smell |
| 13 | Dead code kept for one test | `ManifoldMachine.py:815-850` | smell |
| 14 | Thresholds set so the fixtures stay silent (`resolution_ratio_threshold = 3`, `MAX_DETOUR = 4` "the detours found so far are two crossings long") | `TangleWorkbench.py:231`, `DualWalk.py:86` | cosmetic |
| 15 | Recipe constants with no stated reason: zones at f(q0) not q0 for k28/k32/inversion; heteroclinic 16/10 chosen because a since-removed law failed at 16/9; `outer_blasts=2` because otherwise a class is unresolved | `tests/cases.py:313-328`, `henon_cases.py:42-65` | smell |

The `is_reliable` exclusion list deserves its own line. It counts: no itinerary, ambiguous,
contradicted by evidence, straddling landing. It excludes: "stable end" and "order flip"
evidence (no verdict), member-less classes, cross-side pairs, truncation, the span
fallback, marker disagreement between members, unplaced marker images, unmatched members,
every partition-cut reason. The "order flip" rule is a valid noise detector, but its
consequence is one-sided: the noisy cdists are discarded as evidence while the same cdists
built the hole-bridge ends, the landing spans and the partition that produced the walk.
`CLAUDE.md` states the motive: k=2.8 at 9-13 blasts is reliable "since an order flip gives
no verdict". No code branch is keyed on a fixture; the shaping is in which conditions count.

### In the tests

**Independent oracles, genuinely good.** Area preservation by finite differences, the
same-stability sweep (kept independent on purpose), crossing signs and brackets rebuilt
from nodes, the iterate-table laws gated by the real map, `backward_holes_keep_bridge_side`
(real `map_inv` carries the hole back), `direct_hole_side_is_coordinate_side`,
`row_of_end_is_geometry`, `lobe_enclosure_is_geometry`, `arrangement_image_of`, the D1
series laws, letters recomputed from the partition.

**Tautological or status reads (about a third of the 98).** `bridge_rows_consistent`
calls the library's own violation function. `loops_fold_forward` and
`minimal_trellis_bridges` build their oracle from `image_chain`, the function under test.
`itinerary_pairs_are_classes` and `class_orientation_anchor_outward` use `oriented_class`
and `anchor_outward_key`, the same functions the table uses. `openings_follow_signs` and
`partition_covers_branch` restate D6 with the same formula; the opened half at a non-anchor
bound is never validated geometrically. `every_active_class_walks`, `landings_contained`,
`alphabet_closed_under_words` and a dozen more read a status field. `assert_iterate_relation`
checks `c_iter = stretch * c`, the machine's own assignment. `walks_agree_with_registered_images`
is real on k10 and k28 two blasts and circular on p3, k28 one blast and inversion, where
the word and the evidence come from the same chain.

**No test checks the PDF definitions on a real build.** The strong pip has two synthetic
unit tests; no law re-checks the chosen pip's open box. The pseudoneighbor X-condition has
synthetic tests; surviving pairs on a real build are never verified. Hole propagation is
checked for termination and landing branch, not that the carried bridge lies inside the
container's lobe.

**No test checks a known result.** `a -> a u^-1 a^-1` and `[[1,1],[1,1]]` are pinned only
on a hand-built layout that is handed the k=10 itinerary
(`tests/unit/topology/test_symbolic_dynamics.py:39-130, 291-327`). The real k=10 build is
never compared with the full 2-shift. The no-fixture-fact policy is right for arbitrary
counts; a textbook result is an oracle.

**Weakened or codified demotions.** `test_coarse_blasts_rows_by_signs.py` builds a case
whose polyline rows contradict the signs, then runs only the sign-based laws and asserts
on log counts. The policy says `is_reliable` is never pinned; three regression tests pin it
True (`test_open_p3_deep_runs.py:74`, `test_pip_orbit_cut.py:96`, `test_evidence_diagnosis.py:122`).
Of the 26 not-applicable skips, `every_active_class_walks` on k28 one blast, p3 and
inversion is literally true (the only active class is singleton to singleton) but means
the central walk is exercised by six real builds and the walls law by four.
`arrangement_image_of` is not applicable on the horseshoe.

**Brittle and plumbing.** 280 `caplog` lines, warning counts pinned, private state read
in a dozen tests, 98 facade tests of cache freshness. `helpers/fakes.py` supplies every
pipeline invariant (cut, anchors, ownership, faces) explicitly, so kernel tests never
derive them.

## 6. What to do, in the author's terms

Clarify the algorithm; do not add rules. In this order:

1. **Write the two missing definitions.** The word of a class (the homotopy class of f(b)
   relative to the holes, as a sequence of elements) and the iterated partition (what set
   it is). Until they exist, every walk heuristic and every cut rewrite is unfalsifiable.
   The honest options for the word are: a class needs one computed representative and the
   word is its registered chain (D14 already says a class exists only where bridges exist,
   which contradicts D11); or the walk stays and is labelled a conjecture with its flags.
   The data say both agree wherever both exist.
2. **Finish the pseudoneighbor definition.** Test the open stable interval against X as the
   unstable one is tested. Then re-measure BCG-4 and BCG-8 before deciding DEC-25.
3. **Make N's depth part of the definition.** The pip trellis is cut at the deepest
   registered orbit point; the blast must register it. No `min_separation`, no zone-interior
   filter deciding it, or m becomes an explicit parameter.
4. **One crossing path.** Delete `_collapse_noise_crossings`; an even run is a resolution
   failure reported through the existing flag on every path.
5. **Turn warn-and-continue into unresolved sets.** Table collisions, multi-matches, lost
   pips, marker disagreements, odd chord counts, I2 on image bridges: surface them as data
   with the same consequence everywhere, not as log lines the pipeline walks past.
6. **Compute the refined matrix or stop calling it one.** f(a_i) needs the next partition
   level; the horseshoe's answer is a coincidence of the full shift.
7. **Decide inertness.** Either a definition that says why inert images may be dropped from
   the minimal trellis, or keep and map every hole bridge.
8. **Then shrink.** Arrays for manifolds; the session to its three concerns; zones from
   `Trellis.pips`; plotting and diagnostics out of the algorithm modules; one "element at a
   bridge end"; the test-only API deleted; the 303 decision ids moved from the code into
   one design note.
9. **Add the two oracles the tests lack.** The real k=10 build must give the full 2-shift,
   and every surviving pseudoneighbor pair and chosen pip on a real build must satisfy the
   PDF condition by an independent computation.

## 7. Not verified

- Whether the refined polyline is within `area_cutoff` of W^u at depth (no ground truth).
- The effect of the k=2.8 8-blast recompute drift on words (only crossing counts diffed).
- Whether the endpoint rule or a set-cover rule gives correct words at depth.
- Whether a wrong early wall ever propagates (never observed; walls bite once, on one case).
- Orientation-reversing side carrying; no built case exercises it.
- Completeness of the marker rule, argued by counterexample reasoning, not exhibited.
