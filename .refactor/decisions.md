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
