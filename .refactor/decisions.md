# Decisions

What the author has decided. Binding on every later run. Newest last.

## 2026-10-02 — blasting nested tangles (slug `hole-side-and-own-blast`, base `2315204`)

- **Said** (plan verdict, `plan_ok_planner_human.json`): approved with no notes; the plan's
  defaults stand.
- **Consequence 1 — direct-hole side is combinatorial.** A directly punched hole's
  `Hole.bridge_side` comes from the crossing signs, in the style of `row_of_end`:
  `first` = lower unstable cdist end, `outward = +1 if second.stable_cdist > first.stable_cdist
  else -1`, `left iff first.crossing_sign * outward > 0`. Not the lobe-polygon area, not a
  nearest-vertex measurement.
- **Consequence 2 — a propagated hole's side is read on its image sub-arc** (the `span` inside
  the containing bridge), never against the whole bridge polyline.
- **Consequence 3 — each nested zone blasts only its own fixed point**
  (`fixed_point=[zone.fixed_point]` in `build_nested`); blast order no longer matters.
- **Consequence 4 — CLAUDE.md** carries one sentence: the period-3 side trip at 15-16 steps and
  6+ blasts is fixed, but those runs are unreliable (an unreachable class, a virtual `new1`;
  open). No other CLAUDE.md change.
- Left open by the plan's "Not doing" list (not decided, not to be assumed): p3 >= 6 blasts
  unreliable; `ManifoldMachine.iterate_manifold -> merge_manifolds` non-contiguous crash;
  `_hole_openings` whole-polyline nearest vertex; pip choice lost on generation bump (`_repin`);
  `Blast.py` tidy-up; constant nested slug; noisy inert-word WARNING; blast speed.

## 2026-10-02 — outcome verdict on `hole-side-and-own-blast` (first outcome)

- **Said** (`debrief_ok_auditor_human.json`): "figures", no further notes.
- **Consequence:** proof figures were made on the branch (commit `f8a6d73`, three scripts under
  `scripts/refactor_proof/hole-side-and-own-blast/`; PNGs untracked under
  `figures/refactor_proof/hole-side-and-own-blast/` in the branch worktree) and the auditor
  re-audited (`outcome_auditor_human_2.md`, second pass). The proof scripts are outside the
  plan's five-file list by instruction; rubric 5 is judged against `3b33baa`.

## 2026-10-07 — open-issues round (catalog `docs/open_issues.md`, base `684d70b`)

- **Said** (planning session answers, 2026-10-07), recorded as D1-D8 in `docs/open_issues.md`:
  - **D1 full recompute.** "We only start doing the topology after all the numerics is done …
    any growth that happens requires entirely recomputing the topology." Planned loop: blast,
    compute the topology, blast again, recompute, look for topological changes. The class
    count "may increase multiple times, strictly increasing, never decreasing. Linear growth
    over many iterates is a clear bug though."
  - **D2** the alphabet is the minimal trellis's bridges only.
  - **D3** "the shortest walk on the dual graph gives the word, that is the algorithm" and a
    walk must always exist; the virtual `new1` class is "a very clear bug"; f(bridge) makes the
    word for singleton bridges that cannot walk. The period-3 path "seems buggy".
  - **D4** "Each anchor point should have two stable branches and two unstable branches …
    Each of those individual branches of course have a left side and right side." The left of
    one inversion branch is the right of the other; zones break the space up and must be kept
    in mind; cross-side is a flag only for the SAME branch, never for a bridge between two
    branches.
  - **D5** keep the period-k +1..+(k-1) iterate rule, fix the gap where table links are missing.
  - **D6** signs decide hole openings and rows; geometry only checks.
  - **D7** heteroclinic pairs give holes in B's partition.
  - **D8** simplify freely ("make as many choices as you'd like to make the design simpler …
    robust and professional"); public-API alignment is a later, separate plan.
- **Consequence:** WP0-WP12 in `docs/open_issues.md`, each run analyst -> planner ->
  test-engineer -> coder -> code-reviewer on its own branch. Claude's calls C1-C6 under D8
  (anchor-only branch = one element, element identity by bounding holes, det<0 out of scope,
  innermost-zone blasting, new-image-only same-stability count, per-call `max_points`) are
  revisable by the author.

## 2026-10-08 — DEC-19, the period-k growth stagger (raised by WP5)

- **Said:** "Accept the asymmetry" (offered: trim the trellis to complete orbits; warn and grow more; accept).
- **Consequence (D9):** a stable branch whose reference pair's +j images lie beyond the computed manifold has no fundamental-segment hole; that is correct for the finite trellis, not a bug. WP5's WARNING naming the branch stays. Laws that assume every branch has a window hole (e.g. the inversion anchor class leading its tangle) must be scoped to branches that have one, not left as known issues.

## 2026-10-08 — DEC-20, DEC-21, DEC-22 (raised by WP7 and WP8)

- **DEC-20, said:** "Cut S at the pip orbit." **Consequence (D10):** the partition and the dual graph use the stable manifold cut at the strong pip's orbit (chalkboard stage 3, as pseudoneighbors already do); this replaces the 2026-09-14 "the pip moves the unified set without trimming" rule.
- **DEC-21, said:** "I'm not sure this is an issue? We do not need computed bridges to find words, we strictly use the dual graph unless it is a singleton. That allows us to potentially find words even though we haven't computed the actual iterates. A word can include the class itself and any other classes." Confirmed on follow-up: "Yes, but mark them." **Consequence (D11):** a walk pair no computed bridge spans is a legitimate class: lettered, walked from its elements' landings, NOT a reliability failure, and marked "no computed bridge yet" (informational). Supersedes C9's "a pair no bridge spans is never a class".
- **DEC-22, said:** "Balanced trellis." **Consequence (D12):** cut the crossing set at an unstable depth matched to the stable cut at the pip, so blasting cannot reshape the fundamental segment at depth; the class count then never drops under blasting.

## 2026-10-08 — WP10 (numerics and performance leftovers, slice 1), Claude's calls under D3 / D8

Revisable by the author; recorded so later runs know they are calls, not the author's words.

- **C11 (under D3) — incomplete evidence gives no verdict.** A registered image chain whose
  word parts from the walk because the image passes the end of the computed stable branch
  (diagnosis `"stable end"`: two consecutive chain crossings on one stable branch with one
  sign) misses its crossings with the trimmed-away part, so it is incomplete: `word_agrees`
  is None, the diagnosis is kept and logged at INFO, and it never makes `verified` False or
  the dynamics unreliable. A `"walk"` disagreement stays a WARNING and a reliability failure.
  Effect: k28 at 6-7 blasts is reliable; at 8+ blasts the DEC-23 ties keep it unreliable.
- **C12 (under D8) — the law tier builds each case once per RUN.** `law_case` is session
  scoped (was module scoped: every one of ~12 law modules rebuilt all 11 cases, 40 s of
  setup); the fingerprint guard still fails a test that mutates the shared build.

## 2026-10-09 — DEC-23, equal shortest walks (raised by WP9)

- **Said:** chose "Non-crossing constraint" (offered: non-crossing constraint; registered image decides; stay ambiguous).
- **Consequence (D13):** when equal-length shortest walks read different words, the walk is enriched with the constraint that f(bridge) cannot cross the images of the unstable arcs bounding the lobe; the side is picked combinatorially, with no registered image needed. Applies to the heteroclinic lobe tie and the same-kind ties on k28 at >= 8 blasts.

## 2026-10-09 — WP9b (the non-crossing rule, D13), Claude's calls

Revisable by the author; recorded so later runs know they are calls, not the author's words.

- **C13 (under D13) — the marker rule applies to every walk with a member lobe, not only to
  ties.** D13 is worded for equal walks, but its physics (f is a homeomorphism, so the image of a
  member's lobe holds exactly the images of the markers the lobe holds) binds every walk. On the
  law case `heteroclinic` the class `bh` has a UNIQUE shortest walk (trivial, empty word) that
  holds the wrong markers and that its registered image contradicts; a rule restricted to ties
  would leave that contradiction (and `is_reliable` False) on a law case. Effect: `bh` reads
  `bi bj bi^-1`, its registered word, through a two-crossing detour.
- **C14 (under D3 / D13) — the walk runs between the landing SPANS.** D3's "from the landing of
  one end to the landing of the other" is read exactly: the walk leaves from and arrives at the
  stable edges the image span `f(X)` meets, not every edge of the iterated element (the hole
  bridges bounding `X` map onto the kept image bridges at the ends of `f(X)`, which `f(b)` cannot
  cross). This alone removes every tie on k28 at 8-13 blasts and k32 at 8-10 blasts, and is a
  no-op on every law case.
- **DEC-24 raised** (open): a member-less connecting class whose equal walks differ by a loop
  around another tangle, with no heteroclinic iterate link in N (`heteroclinic`'s `bn`). See
  `docs/open_issues.md`.

## 2026-10-09 — WP11 (simplification sweep, slice 1), Claude's calls under D8 / C3

Revisable by the author; recorded so later runs know they are calls, not the author's words.

- **C15 (under D8) — the topology reads no absolute cdist slack.** Anchors are told by
  identity (`Intersection.is_anchor`: synthetic at exactly (0, 0)), partition ownership,
  containment, the pip-orbit cut, the fundamental-segment overlap and `interval_at` are
  exact, and the only slack left is the RELATIVE error a scaled cdist carries
  (`SCALING_RTOL` at a scaled landing end, `collision_rtol` in `Trellis.same_point`).
  Measured before the change: forcing the `cdist_tol` that topology and loom read to 0
  (and to 1e-12) moved no output on any law case, deep run or blast series, and the full
  suite passed. `cdist_tol` stays for geometry only: degenerate bridge cuts in numerics
  and clipping a polyline to its arc. The cartoon identifies partition boundaries exactly
  too (they are registry floats).
- **C16 (under C3) — orientation is a FixedPoint fact.** `FixedPoint.preserves_orientation`
  is the sign of the solver's single-step Jacobian (`partial_jacobians[0]`, any period;
  None before a solve), and `Trellis.orientation_preserving` reads it; the topology never
  evaluates the dynamical system's Jacobian. The six side flips are one rule,
  `topology.side_after(side, steps, preserving)`, kept although det J < 0 stays out of
  scope (C3).

## 2026-10-09 — DEC-24, a member-less connecting class (raised by WP9b)

- **Said:** "Take whatever path is the most logical with our foundational setup. It is not that tricky. 'One class connects the two tangles and has no computed bridge' is a bit nonsense. If there is no computed bridge that connects the tangles then there is no bridge class which connects the tangles either. Treat every bridge the same. They just connect two intersection points. Those two intersections lie in stable partition elements. Side matters for those partition elements. A bridge connecting two partition elements on the same tangle is no different than a bridge connecting two partition elements on different tangles."
- **Consequence (D14):** no special case for connecting classes anywhere. A member-less connecting class must be traced to whatever produced it (the walk or landing that named the pair) and fixed there; read together with D11, a walk-named pair is a class like any other, but a pair across tangles that no bridge connects signals a construction error to root out, not a class to model.

## 2026-10-09 — WP9c (every bridge is the same, D14), Claude's calls under D13 / D14

Revisable by the author; recorded so later runs know they are calls, not the author's words.

- **Flag for the author: D14's premise does not hold on the fixture.** `bn` = {A:R_1, B:L_(2.0;1)}
  is no construction error: `named_by` gives `bn <- bl <- bi`, and `bi` has two COMPUTED members,
  `(91, 92)` and `(93, 94)` (arcs of W^u(p3 2.0) from p3's stable branch to the heteroclinic
  crossings on p1's window). `bl = f(bi)` and `bn = f^2(bi)` are member-less only because the
  images of 92 and 93 lie past p3's unstable cut (D12): D11's situation, and at p3 17 unstable steps
  they are computed. D14's rule ("treat every bridge the same") is what fixes it: the gap was that
  D13 had a rule only for classes with a member LOBE.
- **C17 (under D13 / D14) — walls.** D13's own wording, "f(bridge) cannot cross the images of the
  unstable arcs", for every class: a class's arcs lie on its members' unstable branches, or one map
  step past its namer's (D11); arcs on different unstable branches are disjoint, and so are their
  images, so the decided walk (`unique` / `trivial`) of a class walked EARLIER whose branches are
  disjoint from this class's is a wall, and a tie loses the walks whose chords interleave with a
  wall's on a face ring (`shortest_walks(walls=)`, `walk_chords`, `DualGraph.chords_cross`). A tie
  filter only (a unique walk is never second-guessed); when every walk crosses a wall the tie
  stands. Only a walk of multiplicity 1 is a wall or is judged by one: a walk standing for
  several face paths names none of them (review fix). Effect: `heteroclinic`'s `bn` resolved with the word base printed first, the dynamics
  reliable; at 17u the connecting class `bm` resolved to its registered word; nothing else moves
  on any law case, deep run or series step.
- **C18 (under D14) — what is tested.** D14 is tested as "a class connects two tangles only where
  a computed bridge in the trellis connects them" (the clause in `alphabet_closed_under_words`),
  not as "every connecting class has a computed member", which contradicts D11 (`bl`, `bn`).

## 2026-10-09 — WP10b (leftovers), Claude's calls

Revisable by the author; recorded so later runs know they are calls, not the author's words.

- **C19 (under D8; PK-4) — foreign chords pair in STABLE order.** A hole's image chain can
  leave its base branch for an excursion across another stable branch (another orbit
  branch or fixed point). The lobe H (bridge plus the stable segment between its ends)
  holds no point of any stable branch of the trellis, so f(H) holds no anchor; walking any
  OTHER branch from its anchor, the branch enters and leaves f(H) at successive crossings
  with the image chain (Jordan alternation), so the empty stretches there are its crossings
  paired in stable order. The base branch keeps the author's chain-order pairing
  (2026-09-23). The base is read at a non-anchor end (an anchor lies on every branch, D4),
  and the two "cannot happen" branches (a base across two branches; a non-anchor crossing
  abutting a stretch on another branch) are assertions. Checked by point-in-polygon on
  p3_6_blasts, p3_15, p3_16 and heteroclinic; changes refinement labels only (class tables,
  holes, words, refined rules and matrices identical).
- **C20 (under D1; NUM-12) — blast image curves are computed manifold.** The workbench keeps
  every image curve a blast computes and re-indexes it (query-only) in every recompute until
  growth has merged it into a branch; `create_bridges` cuts it like any curve, and one
  crossed fewer than twice is held whole as a partial piece. A recompute, a growth round or
  a zone (re)definition after a blast keeps every crossing, bridge, partial piece and
  `iterated` flag, and blast / recompute / blast equals blasting straight.

## 2026-10-09 — WP12 (re-verify and document), Claude's calls

Revisable by the author; recorded so later runs know they are calls, not the author's words.

- **C21 (under D3) — an order flip at a member end is noise, no verdict.** The twin of C11. An
  evidence (a registered member image) whose word parts from the walk, where a member end and a
  bound of that end's own element have registered `+1` images in the REVERSE stable order on
  one branch, is diagnosed `"order flip"`: the map preserves the order along a branch, so the
  registered cdists contradict the dynamics there (k=2.8 at 9 blasts: end 8.5e-6 outward of its
  bound, images 1.2e-6 the other way, both crossings resolution ratio ~1e4: a deep lambda-lemma
  copy hugging the hole-bridge end). `word_agrees` None, INFO, never a reliability failure.
  Lookup only. Effect: k=2.8 at 9-13 blasts reliable; k=3.2 from 9 blasts still unreliable
  (INT-11 kind 2: no registered preimage to prove a flip, the resolution floor, NUM-3/4).
  Law cases unchanged.
- **C22 (under D1 / D8) — growth un-trims a zone; restore never shortens.** A stable trim is a
  view: `grow_n_times` re-finds every branch's physical end (keeping a trim across growth would
  hide every crossing the growth was meant to discover, D1), and the topology never depends on
  the trim (the trellis is cut at the pip orbit, D10). `ResonanceZone.restore` leaves a branch
  that growth has run past alone (it used to cut it back: k=2.8 two blasts 362.9 -> 42.25). A
  zone that re-trims itself after growth is D8 API work.
- **Reporting only — the unmarked lobe.** A crossing on the image lobe of a hole with no
  stretch of its own marks nothing by the author's 2026-09-23 rule; it is now skipped as
  `"unmarked lobe"` at INFO instead of the "grow or blast" WARNING (the fold IS registered).
  No partition, class or word changes.
- **DEC-25 raised (BCG-8, open): the inversion series loses classes at blast 6.** The
  references `(23,24), (26,27)` are replaced by `(56,57)` alone, a clean witness (0.16 / 0.26 in
  stable cdist), no noise; its `+1` pair is never spanned because the reference bridge lies
  outside the zone (NUM-21), so the second stable branch has no window hole (D9's stagger):
  holes 9 -> 5, seed 9 -> 6, two classes ambiguous and contradicted. D12 does not stop it (the
  pip-orbit laps keep deepening). It contradicts D1's "a decrease is a bug" while following the
  pseudoneighbor definition. Options: (a) accept (the count converges as N deepens; D1 holds up
  to reference switches that are pure discoveries); (b) compute the reference arcs' images even
  outside the zone (needs growth, the image is partial); (c) refine old reference pairs rather
  than replace them (DEC-1's open sub-question). Also for the author: does D1 bind the CLOSED
  alphabet across a switch? (k=3.2 at blast 9: closed 45 -> 39 while the seed stays at 34 and
  the holes grow 56 -> 58.)
- **DEC-14 recommendation (the author's call): the default blasting cutoff 1e-7 defines the
  symbolic result.** At 1e-5 the seed decreases on k28, k32, p3 and nested_inner (E12); the k28
  drop is a lambda-lemma copy 5.6e-5 from a pair end that the coarse polyline orders on the
  wrong side, a pair-ordering resolution limit the per-crossing ratio (183) does not see. Coarse
  blasting cutoffs are for exploration only.

## 2026-10-09 — D15, holes propagate backward only on every orbit (reverses D5)

- **Author:** looking at the period-3 itinerary figures, "it seems like we are iterating the
  holes forward for some reason. Whereas in all other cases we only iterate the original holes
  backwards. That is the same for higher period cases."
- **Decision:** no hole at any forward iterate, the `+1 .. +(k_value-1)` appearances on the
  other branches of a period-k (or inversion) orbit included. The other branches get their
  holes from the backward chain, which walks the whole unstable branch cycle.
  `_is_forward_beyond_fundamental` -> `_is_forward_iterate` (`iterate > 0`); the PK-3 gap
  WARNING and `_missing_appearance` are deleted (a missing `+j` costs nothing now); law
  `no_direct_hole_beyond_fundamental` + `window_holes_on_every_branch` -> `no_forward_hole`;
  regression `test_period_k_window_gaps.py` deleted; `("walks_avoid_disjoint_images", "p3")`
  NOT_APPLICABLE (one decided walk left).
- **Measured:** p3 0-5 blasts 12 -> 8 holes, 9 -> 5 classes (the 4 lost were inert; words
  `a -> b -> c -> a u^-1 v^-1`); p3 8 blasts and p3_16 25/27 -> 23 seed, 47/49 -> 43 closed;
  p3_15, p3_6_blasts, nested_inner6 unchanged (their `+j` were never registered); inversion 0
  blasts 8 -> 6 holes, 6 -> 4 classes, every stable branch still holed. Suite 1117 pass.
- **Open (probed, not changed):** past 5 blasts the p3 class count is driven by the reference
  switch (BCG-4): the references become blast-born crossings at unstable cdist ~2600-3100 whose
  backward chain runs ~20 iterates (42 holes). DEC-25's "no window hole on the second branch"
  is no longer a gap under D15.

## 2026-10-09 — the alignment plan (`docs/alignment_plan_2026-10-09.md`, section E, Q1-Q13)

- **Said** (via the orchestrator, verbatim): "Do not drop the refined matrix and iterated partition.
  They are very important. The rest is fine."
- **Consequence — Q5 (the refined level) is a GOAL.** The audit's finding (refined children inherit the
  parent's word; the matrix is "assumed") is to be SOLVED, not removed: plan A.8d defines the iterated
  partition `P1 = P ∨ f(∂P)` with the author's empty-stretch rule (2026-09-23), the refined symbols as
  occurrences, and the refined word; states Theorem T (every member of a class has the same
  `P1`-level freely reduced itinerary), under which the inherited word is the IN-SPLITTING of the class
  graph (conjugate; identical sibling rows; children's columns sum to the parent's); slice B7 makes T
  and the in-splitting identities laws on every build and derives k=10's `[[1,1],[1,1]]` from `[[2]]`.
  Nothing of the iterated partition or refinement is deleted.
- **Consequence — every other question took the plan's default:**
  - Q1 (S-a): the pseudoneighbor test checks BOTH open intervals against X (the stable one too),
    reversing the 2026-09-15 "known asymmetry, no change" (slice B2; the series re-measured before
    DEC-25).
  - Q2 (N-a): the blast never drops a child holding a pip-orbit point; `m` stays "deepest registered"
    (B3). N-b (an explicit `depth`) not chosen.
  - Q3 (I-a): the minimal trellis maps EVERY hole bridge; inertness = transience in the class
    transition graph, read after the words; the evidence fixed point goes (B5). Measured: no
    class-level word moves on 14 builds; the iterated partition gets finer on 6.
  - Q4 (W-a): D3 stands (the walk is the word, the registered chain an oracle); search bounds become
    honest — a walk found under a truncated budget, a span meeting no node, members disagreeing on
    markers are `unresolved` (B6a); the evidence comparison moves out of the algorithm (B6b); the word
    is FREELY reduced (loops dropped, adjacent `c c^-1` cancelled).
  - Q6: the spectral radius stays unshown (2026-09-23 rule).
  - Q7: the manifold data model (arrays) is a later plan (B11).
  - Q8: D1 restated — at fixed depth `m` the hole set and class count are non-decreasing under
    blasting; a deeper `m` is a discovery; D1 does not bind the closed alphabet.
  - Q9: a class exists iff a computed bridge spans it or the word of an existing class names it
    (closure from the spanned classes; D11 and D14 read together).
  - Q10: recipe constants kept, their reasons written (B12).
  - Q11: the suite may pin ONE textbook result, the k=10 full 2-shift (the single exception to the
    2026-10-05 no-fixture-fact policy; B0).
  - Q12: `build_inversion_period2` takes `stable_steps=10` (class `e` resolves, same 6 classes /
    10 holes; B0).
  - Q13 (a): a detection the polylines cannot resolve — `resolution_ratio` below the threshold, or a
    neighbour of the same two curves closer along the unstable polyline than the chord deviation (a
    near-coincident run is one unresolved stretch, whatever its parity) — is registered and reported
    with the remedy "refine" and is NOT a crossing of N on any path; the blast-path parity vote
    (`Tangle._collapse_noise_crossings`) goes. The margin stays 3 and the corpus minimum (heteroclinic
    3.04) is recorded in the rubric (B1). Reverses "flags never discard" for this one case; NUM-4 is
    thereby decided.
- **Claude's defaults pending the author's objection (plan Q14-Q16):** an incomplete footprint (odd
  chain crossings on a branch) is data, never a reliability failure; a build contradicting Theorem T
  stops slice B7 and comes back to the author; no second refinement level is built.

## 2026-10-12 — B0 `oracle-laws` (alignment plan slice B0; plan `artifacts/plan_planner_human_2.md`)

- **Said** (via the orchestrator, verbatim): "Go ahead and run B0 through the planner and coder
  agents", then "Do you need to do more? If you do just do it".
- **Consequence — B0 approved as planned** (`plan_ok_planner_human.json`, iteration 1), with the
  planner's corrections standing: the hole-containment law is CLOSED at an end the iterate table
  links to a container crossing and strict elsewhere; no check-2 `KNOWN_ISSUES` (none failed on
  base); `inversion_period2` joins `every_active_class_walks`' NOT_APPLICABLE.
- **Plan question 1 took its default (yes):** the Q11 exception (the k=10 full-2-shift pin) is
  recorded in one sentence of `tests/conftest.py`'s no-pin paragraph.
- **Scope widened by one file** (orchestrator addendum, under the second instruction):
  `tests/regression/test_open_p3_deep_runs.py`, because its loop reads `LAYERS` and failed the
  k10 pin as vacuous on the deep builds. Rule adopted: `laws.LAW_CASES_ONLY =
  ("reference_pairs_clear_of_X", "k10_is_the_full_shift")`, skipped by the deep-run loop, so the
  scaling-only pair oracle never runs at depth (A.5). Six files in all, none under `src/`.
- **Test-engineer and code-reviewer stages ran on the orchestrator's judgement** (the second
  instruction); review passed in one round (`cycle_ok_test_engineer_code_reviewer.json`).
- **Not changed, left for the author:** check 3 bounds no distance between a carried end and the
  container polyline (reviewer's note; measured, see `artifacts/outcome_auditor_human_2.md`).

## 2026-10-12 — B1 `resolved-crossings`, review round 1 (author decision at the revise)

- **Said** (plan assay): approved the plan as written; the three questions take their defaults.
- **Said** (after the reviewer's `revise` on rubric 1-2): "Rule stands; widen scope to fix the four tests."
- **Consequence.** Q13 (a) holds everywhere: `N` is the resolved set, the trellis of the pip is cut at the
  deepest RESOLVED point of the pip's orbit, and an unresolved orbit point (k28 x8's f^8(q0) at ratio
  2.991; k32 x9 at 1.48; nested_outer x9 at 1.78) shortens `T(q0)` by one iterate. The remedy is
  "refine". The deep-run numbers that follow (k28 x8 24/15/23 -> 4/4/4; k32 x9-10 holes 58 -> 56,
  closed 39 -> 45, reliable) are reported, never pinned. The plan's file list is widened by the four
  files the failures live in: `tests/helpers/laws.py` (the pip-orbit-cut oracle reads `N`),
  `tests/regression/test_coarse_blasts_rows_by_signs.py`, `tests/regression/test_evidence_diagnosis.py`,
  `tests/regression/test_pip_orbit_cut.py` (rework their preconditions or builds, e.g. a finer
  `area_cutoff` so the detection resolves; never delete the defect they guard). Rubric 1's "none of the
  tests named in section 6 is edited" is amended to these four; rubric 3 lists 18 paths.

## 2026-10-12 — B1 `resolved-crossings`, outcome (review passed in round 2; branch `refactor/resolved-crossings` @ `76d0710`, base `1029fcc`)

- **Said:** nothing new since the revise-round decision above (the rule stands; scope 18 paths).
  The review passed in round 2 (`cycle_ok_test_engineer_code_reviewer.json`, iteration 2); no
  debrief verdict yet. Not merged: the merge is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): N = `resolved_ids(threshold)`
  on seven readers; law cases unchanged; every deep-run change is one mechanism -- the pip orbit's
  deepest detection is unresolved (k28 x8 `f^9(q0)` id 245, ratio 2.991 -- the entry above calls it
  f^8, counted from the pip it is f^9; k32 x9 `f^10(q0)` 523, 1.48; nested_outer x9 `f^10(q0)` 447,
  1.78), is never linked, and `T(q0)` stays at the previous blast's trellis under every further
  blast until refined (k32 x9-10 and nested_outer x9-11 equal x8; k28 x8 equals x7).
- **Open for the author (from the audit):** (1) freeze vs placing an unresolved orbit point by
  `image_cdist`; (2) the two `threshold = None` regression builds (keep / NOT_APPLICABLE / unit
  re-guard); (3) `CLAUDE.md`'s "flags never discard" and k=2.8 x8 fixture facts stay stale until
  B12; `Trellis.py:204` default for the `MinimalTrellis` sparse snapshot (fix at `MinimalTrellis.py:411`).

## 2026-10-12 — B6a `walk-bounds-honest`, outcome (review passed in round 1; branch `refactor/walk-bounds-honest` @ `8a1cdc3`, base `76d0710`)

- **Said** (plan assay): asked for a plain description of the issue first, then approved the plan
  as written; question 1 takes its default (`lobe_image_markers` raises `MembersDisagree`, caught
  in `_resolve_itinerary` before any walk). No debrief verdict yet. Not merged: the merge is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): a walk is a word only when
  established. A truncated search (any status), a landing span overlapping no stable node (new status
  `span_meets_no_node`) and member lobes disagreeing on a marker each leave the class unresolved with
  its reason; an unresolved class lends no wall. `_nodes_meeting` keeps a second tier (overlap > 0) so a
  span narrower than its own scaling error keeps the one edge it lies in (k32 blasts 6/8). No rule,
  reliability or closed-size change on any law case, deep run or series step; heteroclinic at
  `max_walks=2` becomes unreliable (was reliable on a truncated `unique`).
- **Left open:** export `MembersDisagree` from `tanglepack.topology` (out of scope); `CLAUDE.md` still
  describes the old "falls back to all of them" span rule; whether `SCALING_RTOL * |c|` is too coarse
  for narrow landing spans (numerics, later).

## 2026-10-13 — B8 `unresolved-sets`, outcome (review passed in round 1; branch `refactor/unresolved-sets` @ `ae18859`, base `8a1cdc3`)

- **Said** (plan assay, 2026-10-13): approved the plan as written -- eight files, the order-free
  collision rule, `forget_strong_pip` resets the cached trellis, the candidacy test runs after the
  zone is defined; question 1 takes its default (test candidacy inside `resonance_zone` as is,
  `add_resonance_zones` unchanged). No debrief verdict yet. Not merged: the merge is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): an image claimed by two or
  more crossings links none of them (held link removed by `IntersectionRegistry.remove_iterate`), and
  a multi-match is recorded; both go to `registry.unresolved_links`, which the law
  `iterate_table_complete` keeps disjoint from the `+1` sources. `TangleSession.forget_strong_pip`
  exists; `resonance_zone` adopts its pip only when it is a candidate on the trimmed trellis (else one
  WARNING, previous choice kept). No state change on any law case or series step (`unresolved_links`
  empty everywhere but `near_tangency`); the workbench generation counter reads one higher per zone.
- **Left open:** candidates that make `topology()` assert once set (k10 2, k28 x1/x2 4, inversion
  2/14); `registry._generation += 1` written from `IterateInference`; the law's new clauses guarded
  only by a review probe; `check_iterate_table_complete` scans the registry, not N.

## 2026-10-13 — B4a `holes-algorithm-vs-diagnostics`, outcome (review passed in round 1; branch `refactor/holes-algorithm-vs-diagnostics` @ `f166573`, base `ae18859`)

- **Said** (plan assay, 2026-10-13): approved the plan as written -- eleven files, three commits
  (verbatim move, measured side on demand, tests), no behaviour change; question 1 takes its default:
  the inverse-map placement of propagated holes (`_map_backward`) stays in `StablePartition` until
  B4b, and the "no map in `StablePartition`" rubric moves to B4b. No debrief verdict yet. Not merged:
  the merge is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): the 16 geometric readings
  live in `topology/diagnostics.py` (byte-identical move); `Hole` has no `measured_side` /
  `side_source`; the measured side is `diagnostics.measured_sides(trellis, holes)`, built only on the
  flag path (`punch_holes(diagnose=True)`), and `side_source(hole)` is a function. I1/I2 stay
  assertions in `StablePartition` as wrappers over cores in `diagnostics`. Every hole, measured side,
  row diagnosis, flag and fingerprint is identical on the 11 law cases and `p3_coarse_blasts`.
  StablePartition.py 1,947 -> 1,354 lines.
- **Left open:** `_map_backward` (B4b); stale references outside the file list (`numerics/geometry.py:7`,
  `scripts/refactor_proof/hole-side-and-own-blast/hole_sides.py`, `CLAUDE.md` holes bullet and `Hole`
  row); `measured_sides` cannot replay an orbit started from the propagation fallback (B4b removes it).

## 2026-10-14 — B4b `holes-are-combinatorial`, outcome (review passed in round 1; branch `refactor/holes-are-combinatorial` @ `ba23a1f`, base `f166573`)

- **Said** (plan assay, 2026-10-13): approved the plan as written -- nine files plus one
  `TangleSession` docstring line; a hole is its five combinatorial fields; placement moves to
  `plotting.hole_point`; question 1 takes its default (a): `hole_point` carries the reference point
  back with `trellis.dynamical_system.map_inv`, so figures are identical and `StablePartition` calls
  no map; the new WARNING watches scaled span ends only. No debrief verdict yet. Not merged: the merge
  is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): `Hole` is
  `(near_intersection_id, pair, bridge_side, bounding_ids, openings, iterate, origin)`; `punch_holes`
  punches references only and asserts no backward pair; propagation starts at the reference, no
  dedups, no fallback; `_hole_openings` names the anchor branch by `stable_key_at` plus the one-curve
  identification; `StablePartition.py` 1,091 lines, no numpy, no map; drawing is
  `plotting.hole_point` (which calls `map_inv`), also read by `measured_sides` and three laws. Every
  hole's five fields, fingerprint, word and measured side identical on the 11 law cases and
  `p3_coarse_blasts`; drawn points identical (0.0). Logs: `near_tangency`'s 14 "cannot propagate"
  WARNINGs gone; one new NUM-5 WARNING (fires on `k28_deep_blast_sliver`, no law case).
- **Left open:** `hole_point` raises a bare `StopIteration` on a propagated hole with no reference
  (unreachable through `Trellis.punch_holes`); stale "seeds already-punched orbits and regions"
  comments (`Trellis.py:1250`, `StablePartition.py:306`); I1/I2 cores in diagnostics (B4a);
  `CLAUDE.md` holes bullet and `Hole` row (B12).

## 2026-10-14 — B2 `pseudoneighbor-both-intervals`, outcome (review passed in round 1; branch `refactor/pseudoneighbor-both-intervals` @ `9cd29d5`, base `ba23a1f`)

- **Said** (plan assay, 2026-10-14): approved the plan as written -- one
  `_interval_is_clear(..., stability, ...)` with a min/max m-range, `per_step_beta(stability)`,
  candidates `own_intersection_ids` with the D7 guard mirrored; a pair kept only when both open
  intervals are clear; three unit tests; e13 over the seven E12 cases at both cutoffs (14 cells); the
  one question takes its default: e13's results jsonl is committed as the fourth file. No debrief
  verdict yet. Not merged: the merge is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): the pseudoneighbor test now
  matches the PDF definition (both open intervals against X, one function). No behaviour change on
  any law case or any of e13's 14 cells (0 stable fires in 182 steps). Seven unit-test fixtures whose
  data contradicted the definition were changed in data only (probe crossings moved to the other
  stable branch, or nudged); the reviewer accepted this against rubric 2's literal wording. DEC-25
  gets no new input: the inversion series still drops 2 -> 1 references at blast 6, with the stable
  clause playing no part.
- **Left open:** the m-range `min/max` is guarded only by an uncommitted hand check; two older tests
  (`walk_stops_at_reference_window_end`, `heteroclinic_candidate_ignored`) have no teeth for their
  named rule; `CLAUDE.md` "a crossing past a cut is a point of X" wording (B12); the old function name
  in `docs/` and `docs/open_issues.md`.

## 2026-10-14 — B3 `depth-of-n`, outcome (review passed in round 1; branch `refactor/depth-of-n` @ `7d3a294`, base `9cd29d5`)

- **Said** (plan assay, 2026-10-14): approved the plan as written, question 1 = (a) -- N-a as
  decided and restated: a child ending on the chosen pip's composed `+1` orbit is exempt from
  `min_separation` (one INFO), `BlastResult.pip_depth` composed (never `forward_depth`), the
  `blast_zone` docstring, `SeriesStep.pip_depth` measured test-side, the series law
  "non-decreasing" (question 2's default), the proximity-guard test restated, a new unit-test module
  on a multi-iteration call. Five files. The series stalls -- (ii) unresolved `f^{m+1}(q0)`, (iii)
  undetected, (iv) the zone filter -- left for a later decision. No debrief verdict yet. Not merged:
  the merge is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): `m` = the registered `+1`
  chain from the chosen pip; `min_separation` never decides it. No change on any law case or on any
  step of the five series (pip_depth k28 `[1..8, 8]`, k32 `[1..9, 9, 9]`, nested_outer
  `[1..9, 9, 9, 9]`, p3 `[14..23]`, inversion `[5..9]`). Inside one multi-iteration call a
  pip-bearing child is kept: nested_outer 11 iterations ends at m = 9 (was 7), and that call now
  fires the existing 361 same-stability and 1 unresolved-detection WARNINGs (the exemption keeps
  near-coincident curve). One INFO per exempted child (one per single-iteration blast with the
  recipes' 1e-5). The child filter runs after every parent of the iteration is mapped (no effect of
  its own, instrumented by the reviewer). Suite 1193 pass / 1 xfail (KNOWN_ISSUES).
- **Left open:** the stalls (ii) on k28 / k32 / nested_outer / inversion, (iii) on p3 /
  heteroclinic, (iv) the zone filter -- options: refine (B1), place by `image_cdist` (PK-6-like,
  reverses B1 for the pip orbit, drifts toward N-b), a finer image near the anchor for (iii); whether
  the guard should win over the exemption when the curve is truly under-resolved; the stale B3
  section and decision N in `docs/alignment_plan_2026-10-09.md` (`:149-154`, `:903-943`,
  `:1166-1167`); `scripts/henon_blast_period_3.py` dead on `main` (pinned pip 12 not a candidate);
  the test-side orbit composition written three times.

## 2026-10-14 — B5 `inertness-is-transience`, plan round 1 sent back (author note at the revise)

- **Said** (`artifacts/plan_ok_planner_human.json`, status `revise`): "I believe you can determine if
  a bridge is inert before constructing the dual graph. Every bridge that has an initial hole punched
  in it should be inert. The zeroth holes."
- **Consequence:** the analyst tested the rule on all 15 builds (`artifacts/findings_analyst_planner.t2.md`).
  It holds everywhere: every class holding a zeroth-hole (reference-hole) bridge has an EMPTY word,
  and it is a theorem when no selected orbit sits inside the zeroth lobe and no other hole sits beside
  the bridge (the lobe holds one marker, its image holds none). It is sufficient, not the definition:
  heteroclinic `u` is a hole-free inert class. Re-planned: the minimal trellis maps every hole bridge
  EXCEPT one carrying a zeroth hole, decided from `holes` before the dual graph (no class table read);
  transience over the words (Q3, I-a) stays the verdict for every class; a new law
  `zeroth_hole_classes_have_empty_words` checks the theorem. Minimal trellis, `P1`, image bridges,
  refined classes and words stay as at `HEAD` on all 15 builds (round 1's "P1 finer on 8 builds" and
  its three new D9 WARNINGs are gone); only five verdicts change.

## 2026-10-14 — B5 `inertness-is-transience`, outcome (plan approved at round 2; review passed in round 1 of coding; branch `refactor/inertness-is-transience` @ `3fb310a`, base `7d3a294`)

- **Said** (plan assay, round 2, 2026-10-14): approved the round-2 plan as written -- the author's
  zeroth-hole rule (a bridge carrying a reference hole, iterate 0, is inert: its lobe holds one
  marker whose image is no marker, so its word is empty) decides, before the dual graph, which hole
  bridges the minimal trellis maps (every hole bridge except the zeroth-hole ones, read from
  `holes`); transience over the words is the verdict for every class (Q3, I-a); new law
  `zeroth_hole_classes_have_empty_words`; nine files (question 1's default). No debrief verdict yet.
  Not merged: the merge is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): the registered-evidence
  inertness is deleted (`_resolve_inertness`, `_collect_image_evidence`, `has_image_evidence`, the
  `image_classes` / `image_loops` / `unresolved` slots); `entry.inert` is `None` until
  `symbolic_dynamics` sets it from `inert_classes(words)`, then letters (active via `letter_for`,
  inert series); `Topology.classes` runs the dynamics; PK-8's WARNING gone. Words, minimal
  trellis, `P1`, refined classes and `is_reliable` unchanged on all 15 builds; five verdicts active
  -> inert (inversion_two_blasts `g`, heteroclinic `d` / `bj`, p3_6_blasts `ab`, nested inner 6
  `ae`); `unmapped` shrinks by those bridges. Suite 1193 pass / 1 xfail (KNOWN_ISSUES).
- **Left open:** the zeroth-hole rule is a sufficient condition, not the definition (heteroclinic
  `u` is a hole-free inert class); the unblasted nested fixture `p3_partitioned` has a zeroth-hole
  class (bridge `(4,5)`) that is UNRESOLVED and therefore active while its bridge is not mapped
  (outside the theorem's hypothesis); `inert_classes` not re-exported from `tanglepack.topology`;
  `CLAUDE.md` inertness text stale (B12).

## 2026-10-14 — B7 `refined-level-defined`, outcome (plan approved round 1 with question 1's default, option (ii); review passed in round 1 of coding; branch `refactor/refined-level-defined` @ `a6a4f68`, base `3fb310a`)

- **Said** (`artifacts/plan_ok_planner_human.json`, approved): "Question 1 takes the default, option
  (ii): keep the base-branch chain-order chord pairing in B7, register
  (members_share_the_P1_itinerary, k32_series) as a strict-xfail KNOWN_ISSUE naming B7b, and plan B7b
  base-chords-in-stable-order as the next slice. Seven files; no number changes." No debrief verdict
  yet. Not merged: the merge is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): `P1 = P ∨ f(∂P)` and the
  empty-stretch rule (base chords chain order 2026-09-23, foreign chords stable order C19) stated in
  the `PartitionFamily` module docstring; the odd-fold-count WARNING is gone, replaced by
  `IteratedHomotopyPartition.incomplete_footprints` (+ one INFO per build, a `describe()` line),
  collected by a sibling `_incomplete_footprints` because `_empty_stretches`' signature is
  monkeypatched by out-of-scope tests; Theorem T and the in-splitting stated in `SymbolicDynamics`
  (docstrings only, AST-equal); laws `members_share_the_P1_itinerary` (law cases only) and
  `refined_matrix_is_an_in_splitting`; the series run of T under its own id
  `series_members_share_the_P1_itinerary`, strict-xfail on `k32_series` (wp `B7b`). No number changes
  on the 15 builds; suite 1200 passed / 46 skipped / 5 deselected / 2 xfailed.
- **Deviations the author did not approve in advance:** an EIGHTH file
  (`tests/invariants/test_law_case_sanity.py`, +2/-1: the series id in `_STANDALONE_ISSUES`, the wp
  regex accepts `B\d+[a-z]?`), ordered by the orchestrator because the seven-file plan could not carry
  the approved known issue past that guard; the known-issue key carries the `series_` prefix, not the
  verdict's literal `members_share_the_P1_itinerary`.
- **Left open (for B7b or the author):** the k32 class is `{R_(0.0;8), R_(0.0;10)}` in ElementNaming
  (the KnownIssue text's "R_7 <-> R_9" is the positional `R#7`/`R#9`); the walked word agrees with the
  siblings and only the hole-bridge member `(6, 7)` disagrees; the module docstring's "laws on every
  build" overstates T (law cases + series only); B7b reverses the author's 2026-09-23 base-branch
  chain-order rule; commit trailers say `Claude Opus 5.5` though the orchestrator asked for another name.

## 2026-10-15 — B7b `base-chords-in-stable-order`, outcome (plan approved round 1 with both defaults; review passed in round 1 of coding; branch `refactor/base-chords-in-stable-order` @ `63477e3`, base `a6a4f68`)

- **Said** (`artifacts/plan_ok_planner_human.json`, approved 2026-10-15): "Approved as written. Both
  questions take their defaults: question 1 no fold of _incomplete_footprints into _empty_stretches;
  question 2 yes, the geometric reading goes in tests/helpers/laws.py with a series test in the
  regression tier. Seven files. The base-branch chords pair in stable order, reversing the 2026-09-23
  chain-order rule on the base branch." No debrief verdict yet. Not merged: the merge is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): ONE chord rule on every stable
  branch -- a lobe image's interior crossings pair in STABLE order, base branch included (C19's Jordan
  argument; `_empty_stretches` drops the `branch_key != base_branch` test and the dead
  `(base_branch,) = ends`). The 2026-09-23 chain-order clause is REVERSED. Evidence: on all 47 base-lobe
  readings of the corpus (16 builds + 47 series steps) the gaps alternate out/in from the anchor (base
  inside) and the covered gaps are the stable-order pairs; the one lobe where the orders differ is
  `k32_series` hole `(6,7)` at steps 2-10; the stop rule did not fire. Behaviour changes on `k32_series`
  only: iterated elements (2-5 23->24, 6-7 40->39, 8-10 113->111), `unmatched_members` 1->0 and Theorem T
  1->0 disagreements at 6-7, `refined_rules` at 8-10 (inert `v` 4 children -> 3, one token `v_4 -> v_3`).
  Letter-free words, `M`, `M_1`, all class counts, split counts and `is_reliable` identical on all 63.
  The request's "split inert -1 on 15 builds" was an artefact of B7's buggy scratch patch (it lost the
  `empty.setdefault(b, ...)` line). KNOWN_ISSUES shrinks to C3 alone. New test reading
  `helpers.laws.base_chord_misreadings` (point-in-polygon, clauses a-d) on 4 builds + 5 series. Suite
  1211 passed / 47 skipped / 5 deselected / 1 xfailed.
- **Left open:** the CLAUDE.md "Partition families" bullet still states base chords in chain order
  (B12); `_incomplete_footprints`' INFO log, `describe()` line and the class attribute doc still say
  "last" crossing (now the outermost); the new build test reads 4 of the 11 build lobes (the other 7
  read clean by hand); outside the k32 lobe the base tests cannot catch a DROPPED chord (the one inside
  gap is also a hole stretch), only a wrong one; commit trailers say `Claude Fable 5.1`.

## 2026-10-15 — B6b `word-finished`, outcome (plan approved round 1 with all six defaults; review passed in round 1 of coding; branch `refactor/word-finished` @ `65edec4`, base `63477e3`)

- **Said** (`artifacts/plan_ok_planner_human.json`, approved 2026-10-15): "Approved as written. All six
  questions take their defaults: (1) WalkSearch.status has three values resolved/ambiguous/unresolved
  and 'trivial' belongs only to a class (empty word); (2) a search cut short with candidates becomes
  unresolved with reason 'truncated', walks kept, the truncated field goes; (3) class-level free
  reduction per A.8c as one free_reduce applied to cd.symbols and to the closure/inertness words,
  translate_itinerary stays unreduced so the T law keeps its same-occurrence reading, the registered
  evidence word is reduced by the same function before comparison; (4) members_disagree stays on the
  class only (six walk reasons); (5) MembersDisagree exported from tanglepack.topology; (6) one slice,
  ten files." No debrief verdict yet. Not merged: the merge is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`):
  - "trivial" now means only an empty class word. `WalkSearch.status` takes the values `resolved`,
    `ambiguous` and `unresolved`. `WalkSearch.reason` holds one of six failure values (`unreachable`,
    `no_start_node`, `no_goal_node`, `no_admissible_walk`, `truncated`, `span_meets_no_node`). The
    `WalkSearch.truncated` flag is gone, and a search cut short by the budget keeps its walks but is
    `unresolved`.
  - `SymbolicDynamics.free_reduce` cancels at class level. It has three readers: `cd.symbols`, the
    closure/inertness words and the registered evidence word. `translate_itinerary` and the T law
    are unreduced.
  - The evidence comparison now lives in `topology/diagnostics.py`, behind
    `verify_words(records, trellis, naming, table, letters)`. `_compare_evidence` was folded into
    `verify_words`. `_member_pair` stays in SD.
  - No number changes on the 63 corpus items. Free reduction fires 0 times there.
  - Suite: 1223 passed / 47 skipped / 5 deselected / 1 xfailed (C3).
  - Line counts: SD 1537 -> 1312, DualWalk 1370 -> 1383, diagnostics 713 -> 990 (exactly at its cap).
- **Deviations the author did not approve in advance:**
  - `WalkReason` and `free_reduce` were also exported.
  - The per-class DEBUG line `class %s: %s, word %r, verified %s` was dropped.
  - A capped search with several walks logs the max_walks WARNING (now ending "the search decides
    nothing") and no longer also the ambiguous WARNING.
  - `plotting._class_status` lost its `-` fallback and now returns `cd.status`.
  - The loops unit test still pins the unreduced `translate_itinerary` and adds `free_reduce == []`.
  - The plan's rubric 7 said the capped tie reads `ambiguous` on base. It actually reads `unique`
    with `truncated=True`.
- **Left open:**
  - The `CLAUDE.md` symbolic-dynamics bullet and the `DualWalk` table row still list the old
    statuses and `WalkSearch.truncated` (B12).
  - `diagnostics.py` has zero headroom under its 990 cap.
  - `free_reduce` runs twice per class (closure, then `cd.symbols`), so a cancellation's DEBUG line
    would print twice.
  - `plotting._class_status` is untested.
  - The ambiguous, failure and truncated rows of the status mapping are exercised by unit tests
    only, because the corpus has none.
  - Commit trailers say `Claude Fable 5.1`.

## 2026-10-15 — B10 `dead-api`, outcome (plan approved round 1 with all eight defaults; review passed in round 1 of coding; branch `refactor/dead-api` @ `517aa01`, base `65edec4`)

- **Said** (`artifacts/plan_ok_planner_human.json`, approved 2026-10-15): "Author approved the B10
  dead-api plan 2026-10-15 with all eight defaults": (Q1) rewrite the two `*_new_intersection_structures.py`
  demo scripts with `filter`, dropping their graph and array sections; (Q2) delete
  `TangleSession.classify_bridges` and rewrite its script and test as a comprehension over
  `classify_bridge`; (Q3 (a)) change only `henon_blast_period_3.py:123` and report the script dead;
  (Q4) delete the `from_branch` test with the method; (Q5) drop the `"_get_iterate"` string in
  `e02_gpu_end_to_end.py`; (Q6) leave the notebooks; (Q7) delete `IterateTable.as_forward_array` /
  `as_backward_array` in C5; (Q8) trim `topology.__all__` to 48 names and keep every import line.
  No debrief verdict yet. Not merged: the merge is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`):
  - Deleted: `grow_x_times`, `_get_refined_point`, `_get_iterate`, `trim_stable_at_intersection`,
    the registry's networkx graph, `graphviz.py`, `build_intersection_graph`,
    `visualize_intersection_graph`, `on_cdist_range`, `from_branch`, `nearest_by_unstable_cdist`,
    both `as_*_array` pairs, `classify_bridges`, `BlastResult.terminated_early`. Kept:
    `release_manifold` (live), every law-tier oracle, `refine_manifold`, `filter` /
    `from_fixed_point` / `iterate_orbit`.
  - `import tanglepack` no longer loads networkx. `topology.__all__` 104 -> 48.
  - No number changes on the 63 corpus items (probe byte-identical). Suite 1215 passed / 47 skipped /
    5 deselected / 1 xfailed (C3); 8 test ids gone as planned, 1 renamed.
  - Net −1013 lines over src + tests (−1034 with scripts and stress).
- **Deviations the author did not approve in advance:** an unused `TYPE_CHECKING` import dropped in
  C3; the `Iterable` import went in C4, not C5; `remove_iterate` no longer reads `target_id`; a few
  extra docstrings reworded so none mentions a registry graph; the comment in
  `henon_bridge_classification.py` moved (2 lines, not 1).
- **Left open:** CLAUDE.md rows `:133` (`graph(bridges=)`), `:136` (`graphviz`), `:154`
  (`terminated_early`) and `docs/source/tanglepack.numerics.rst:183-189` (automodule of the deleted
  `graphviz`) are stale (B12); the 56 F401s in `topology/__init__.py` (B12's import block);
  `scripts/manifold_cutting_test.ipynb` and `notebooks/archive/period_three.ipynb` still call
  `grow_x_times` (Q6); `henon_blast_period_3.py` still dead at `:95`;
  `test_map_step_and_graph.py` keeps "graph" in its name; commit trailers say `Claude Fable 5.1`.

## 2026-10-15 — B9a `session-three-concerns`, outcome (plan approved round 1 with every default Q0-Q6; review passed in round 1 of coding; branch `refactor/session-three-concerns` @ `5f0db5f`, base `517aa01`)

- **Said** (`artifacts/plan_ok_planner_human.json`, approved 2026-10-15): "Author approved the B9
  plan 2026-10-15 with every default": (Q0) split B9 into B9a `session-three-concerns` now and B9b
  `one-element-at-end` next; (Q1) no thin plot delegate stays on the session; (Q2 (a))
  `cartoon_zones` moves verbatim to `loom/cartoon.py` and reads the STORED zones -- zones are NOT
  computed from `Trellis.pips` (that would change `inversion_period2`'s circle cartoon and pull the
  zone boundary builder and `ring_cycles` / `ring_order` into topology); (Q3) `plot_resonance_zones`
  goes to the same module, taking zones; (Q4) the describe fan-outs (`describe_holes`,
  `describe_stable_partitions`, `_describe`) are deleted; (Q5) the
  `test_a_pip_trimmed_away_is_not_replaced_by_plotting` test is deleted (`Trellis.plot_*` never
  touches the pip memory); (Q6) no new "element on any side" API. No debrief verdict yet. Not
  merged: the merge is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`):
  - `TangleSession.py` 1744 -> 843 lines; it holds no `plot_*`, no cartoon geometry, no fan-out,
    no `partition_element_for`, no describe fan-out. `cartoon_zones` / `plot_resonance_zones` live
    in `loom/cartoon.py` (exported from `tanglepack.loom`); `cartoon_enclosures(trellis, bridge_ids,
    parents, fixed_points)` in `topology/plotting.py` (not in `topology.__all__`, still 48 names).
  - Callers draw a circle cartoon in three explicit steps (`cartoon_zones` -> `cartoon_enclosures`
    when `parents` -> `plotting.plot_dual_graph_cartoon`); per-fixed-point plots go through
    `session.trellis(fp).plot_*`; manual stage calls through `session.topology(fp)` or the trellis.
  - No number, figure or script stdout changes (63-item probe byte-identical, 32 canvases, 42
    itinerary PNGs, 11 scripts). Suite 1208 passed / 47 skipped / 5 deselected / 1 xfailed (C3);
    7 test ids deleted, 3 renamed, as planned. Net -749 lines (src -613, tests -143, scripts +7).
- **Deviations the author did not approve in advance:** imports pruned in the commit that stops
  using them; the `plotting.py:6` docstring fix in C2, not C1; small docstring edits in listed
  scripts and `test_partition_elements.py`; a test-engineer commit strengthening
  `test_session_caches._repartition_by_hand` (pre-existing toothless assertion).
- **Left open:** B9b (the element-at-end fold); stale docstring citations of
  `TangleSession.partition_element_for` at `Trellis.py:1413`, `TopologyResults.py:427` (B9b) and of
  "the nested-zone classifier in TangleSession" at `ResonanceZone.py:295`; CLAUDE.md `TangleSession`
  row (B12); `tanglepack.loom.plot_resonance_zones` untested; `henon_blast_period_3.py` still dead
  at `:95`; commit trailers say `Claude Fable 5.1`.

## 2026-10-15 — B9b `one-element-at-end`, outcome (plan approved round 1 with every default Q1-Q5; review passed in round 1 of coding; branch `refactor/one-element-at-end` @ `94fca6b`, base `5f0db5f`)

- **Said** (`artifacts/plan_ok_planner_human.json`, approved 2026-10-15): "Author approved the B9b
  one-element-at-end plan 2026-10-15 with all five defaults": (Q1 (B)) a function-level import of
  `HomotopyPartition` in `bridge_classes`, which calls `owner_at_end`; `_element_at_end` deleted
  (not option A's two shared helpers beside `row_of_end`); (Q2) the class rule's unowned-crossing
  AssertionError becomes `owner_of_intersection`'s ValueError on that unreachable path; (Q3) the
  three tests citing `element_for` are renamed; (Q4) no `plot_resonance_zones` smoke test here
  (B12); (Q5) `ResonanceZone.py` untouched (its two `TangleSession` citations are accurate). No
  debrief verdict yet. Not merged: the merge is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`):
  - `PartitionFamily.owner_at_end` is the one reader of "which element owns this bridge end on its
    row", on every computed path; the class rule's `ends` calls it. `Trellis.element_for`,
    `StablePartition._elements_at_bridge`, the `StablePartitionResult.elements_at_bridge` field
    (with its two writers) and `BridgeClass._element_at_end` are deleted. The any-side reading is
    `owner_of_intersection` on the `stable_key_at` branch; no "any side" API added (B9a Q6).
  - No number changes: agreement probe 7,090 / 14,180 / 1,776 with readings byte-identical to
    base; 63-item outcome probe byte-identical. Suite 1207 passed / 47 skipped / 5 deselected / 1
    xfailed (C3). One test deleted, three renamed. Net -168 lines (src -149, tests -19).
- **Deviations the author did not approve in advance:** docstring wording beyond the analyst's
  candidate (the `bridge_classes` Note, the singular "it" in the `StablePartitionResult` Note);
  the missing-partition ValueError text is now `owner_of_intersection`'s (type unchanged, still
  names the branch); commit trailers say `Claude Opus 5.5` (the spawn asked for `Fable 5.1`).
- **Left open:** CLAUDE.md `Trellis` row still lists `element_for` (`:139` in the working tree,
  B12); `plot_resonance_zones` untested (B12); the unused `Intersection` import at
  `StablePartition.py:238`; `henon_blast_period_3.py` still dead at `:95`; four stale committed
  higher-period PNGs; the rewritten `test_partition_elements.py` does not catch a broken C7 anchor
  rule or a wrong row in `owner_at_end` (the law tier does).

## 2026-10-15 — B12a `one-design-note`, outcome (plan approved round 1 with the B12a/B12b split and every default Q1-Q12; review revise in round 1, pass in round 2; branch `refactor/one-design-note` @ `b39c7e2`, base `94fca6b`)

- **Said** (`artifacts/plan_ok_planner_human.json`, approved 2026-10-15): "Author approved the B12
  plan 2026-10-15 with the split (Q11: B12a now = commit 0 + design note + src de-citation +
  CLAUDE.md; B12b next = laws, smoke test, __init__ trim, rst lists, chalkboard, Q10 reasons, PNG
  regeneration, plan Outcome) and all twelve defaults": (Q1) the 649 MB `figures/henon_period3_*`
  sweeps are NOT committed, their two scripts and CLAUDE.md lines are; (Q2) B12b deletes the 9
  TAUTOLOGY laws, keeps `dual_bipartite_degree`, every MIXED and ORACLE law (-4 ids); (Q3) B12b
  regenerates only the 3 content-stale PNG folders (`heteroclinic_p1_p3`, `nested_p1_p3_inner6`,
  `inversion_blast2`) in their own commit; (Q4) the note is `docs/design_note.md`, sections
  A.0-A.11 (+ H, S, O, I); (Q5) CLAUDE.md's physics sections stay byte-identical, their ids
  (`WP10`, `NUM-20`, `WP10b`, a date) exempt; (Q6) the 8 log strings (7 ids + 1 date) stay in
  `src/`; (Q7) B12b: `topology.__all__` -> 51 with `MembersDisagree`, `WalkReason`, `free_reduce`,
  the other 53 unused imports dropped; (Q8) the other id families (`B*`, `Q*`, `E*`, `F7`, `N-a`)
  and dates leave `src/` too, "design note §A.x" stays; (Q9) Sphinx: a list check only (+12
  modules, -graphviz, `conf.py` path `../../src`), no `-W` build; (Q10) the chalkboard's 17 rows
  closed in place in B12b, the 9 audit divergences closed without a row listed in the note's H now;
  (Q11) two branches, B12a first; (Q12) the recipe reason for zones at `f(q0)` is the measured
  equivalence (zone at `f(q0)` after n blasts = zone at `q0` after n-1, k28/k32, n <= 2), written
  in B12b. No debrief verdict yet. Not merged: the merge is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`):
  - `docs/design_note.md` (722 lines) is the design: one section per stage, every one of the 86
    ids `src/` cited at base anchored, plus H (slices, WPs, superseded rules), S, O, I.
    `src/` keeps only the 7 log-string ids and the 1 log-string date; 15 "design note §A.x"
    pointers. Stripped AST equal on 51 files; src -60 lines (+441/-501), 35 files (34 planned +
    `FixedPointSolver.py`, a recorded deviation the reviewer asked for).
  - CLAUDE.md: "Bridge classes and the dual graph" replaced by "The pipeline, stage by stage" (one
    paragraph per stage, each naming its functions and ending "design note §A.x"); table at HEAD;
    80,246 -> 62,737 bytes (173 -> 185 lines); `:1-27` and Coding Style byte-identical.
  - Suite 1207 passed / 47 skipped / 5 deselected / 1 xfailed (C3), collected ids identical.
- **Round-1 defects** (caught by reading, no script could): the note said `image_of` maps the
  boundary with the real map, listed `from_table` as deleted by B4b (B4b added it), said the
  in-splitting runs on the series (only T does), and that the minimal trellis ignores the class
  table; plus history sentences in `FixedPoint.py` / `FixedPointSolver.py`. All fixed in round 2.
- **Left open:** B12b (laws, smoke test, `__init__` trim, rst, chalkboard, Q10 reasons, PNGs, the
  plan's Outcome section); old-code narration outside the rubric's word list (`BridgeIterator.py`
  "plan row 1.1/2.2 removed", `BaseManifold.py:16-17`, `ManifoldMachine.py:632,700-701`,
  `Tangle.py:593,1017-1018`, `DualGraph.py:50-51`, `geometry.py:9-11`); the note and CLAUDE.md say
  "no fan-outs" though `TangleSession.classify_strong_pips` still fans out over fixed points (kept
  by B9a for pip memory); note cosmetics (`:506` unwrapped, `:377` broken line, `:695` "990-line
  cap" vs 988); commit trailers say `Claude Opus 5.5`. These notes were appended in the MAIN tree
  after commit 0 copied it: the branch's `.refactor/` is one entry behind main's working copy.
