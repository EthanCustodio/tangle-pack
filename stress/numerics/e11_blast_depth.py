"""Deep blasting: many blasts, several cutoffs and separations, topology after each.

e07's ``blasts`` study blasted three cases up to 10 times and looked at the
topology once, at the end. This one blasts up to 20 times and snapshots the
whole topology stack after EVERY blast, so the series shows what each blast
uncovered.

Author's clarification (2026-10-07): blasting exists to uncover new
information. New crossings, bridges, bridge classes, holes, split elements,
new letters and changed words or transition matrices are EXPECTED and
common; they are recorded as discoveries (the ``*_changed`` / ``*_added``
flags), never as failures. Only genuine errors fail: the blast itself raising
(the config outcome), a topology exception or invariant ``AssertionError``
(recorded per blast in ``topology_outcome``, the blasting goes on), the
watchdog (timeout / OOM; every blast is checkpointed), or the numerics
violating a law (``same_stability`` crossings, ``cdist_decreasing`` along a
new bridge). The lobe-area ratio is a sanity probe of the polygonal
approximation, not a law check: a lobe and its preimage lobe are compared as
polygons, which a fold can make self-intersecting.

Cases (``case``)
----------------
* ``k28`` -- the k=2.8 saddle, zone at ``f(q0)`` (``cases.build_k28``).
* ``k32`` -- the k=3.2 saddle by the same recipe: a second single saddle whose
  first blasts already change the words (probe 2026-10-07; k=2.2 has no
  registered ``f(q0)`` at 10 unstable steps and k >= 4 has a single crossing).
* ``p3`` -- the period-3 orbit of the nested map alone, its zone at the
  default pip (``henon_cases.build_period3``).
* ``nested_outer`` -- period 1 + period 3, the outer (period-1) zone blasted.
* ``nested_inner`` -- the same, the outer zone first blasted twice at
  ``min_separation = 1e-4`` (the ``build_nested`` default, so every outer
  class resolves), then the inner (period-3) zone blasted.

Every case is BUILT at ``area_cutoff = 1e-7`` (the canonical recipes), so all
cells of one case start from the same tangle; ``area_cutoff`` is the
refinement threshold while BLASTING (the bridge images are refined by the
manifold machine). ``num_iterations`` is passed to each ``blast_zone`` call.

Grid (full mode, 53 configs; calibrated 2026-10-07)
---------------------------------------------------
For each case ``min_separation`` in {None, 1e-6, 1e-5, 1e-4, 1e-3} at
``area_cutoff = 1e-7`` (plus 1e-2 for ``k28``/``k32``: through 1e-3 the
separation never bites there, see below); ``area_cutoff`` in {1e-5, 1e-6,
1e-8} at the case's canonical separation; ``num_iterations = 2`` at every
canonical cell; repeats 1 and 2 of the canonical cell of ``k28``, ``p3`` and
``nested_inner`` (determinism of the topology sequence: compare
``topology_hash`` series). Dropped: the full cutoff x separation product
(4 x 5 x 5 = 100 cells; every heavy cell runs to OOM or its budget, ~4-15
min, so the product is ~10 h of worker time) and repeats of the heavy
``k32``/``nested_outer`` cells (the k28 repeats and the identical
separations below already show determinism).

A config stops early only on an error, the watchdog (``TIMEOUT``,
``RSS_CAP_GB``; every blast is checkpointed) or its own budgets, which are
reached first: a blast predicted to overrun the wall budget (``_budget`` =
timeout - 120 s, cost ~2x the last blast) or the memory budget
(``_rss_budget_gb`` = 0.9 x ``RSS_CAP_GB``, peak RSS x the last blast's
growth factor) is not started (``stop_reason`` ``budget`` /
``memory_budget``), so the config ends ``ok`` with its last topology kept.

Calibration (2026-10-07, 11 cells, 5-6 workers)
------------------------------------------------
* k28 / k32 / nested_outer at ``area_cutoff = 1e-7``: bridge points grow
  ~1.6x per blast; the 16th blast crosses 12 GB (the watchdog killed these
  calibration cells after 15 blasts, 230-470 s; at 1e-8 after 13 blasts,
  228 s; the memory budget now stops them one blast earlier, cleanly). At 1e-5 k28 is
  budget-stopped after 18 blasts (708 s, 5.8 GB): the longest cell.
* p3 and nested_inner: 20 blasts in 3-5 s.
* Estimate: 31 heavy cells, ~210 worker-minutes, 30-40 min wall at 8
  workers (critical path one budget-bound cell, 13-15 min), peak RAM
  8 x 12 GB = 96 GB.
* Discoveries: k28 at min_separation None, 1e-5 and 1e-3 gives IDENTICAL
  crossing counts and an identical ``topology_hash`` sequence (the separation
  never drops a child); ``num_iterations = 2`` at blast n tracks
  ``num_iterations = 1`` at blast 2n for the first 8 single blasts, then
  differs. The topology changes often: k28's three classes
  become 20 at blast 8 and 57 at blast 15, letters run past ``bw``; the
  canonical words (letters renamed by table position, ``structure_changed``)
  separate a relabelling of a class (its elements renumbered) from a new
  symbolic structure. A cutoff of 1e-5 already changes what blast 7
  uncovers (12 classes vs 4 at 1e-7). nested_outer's outer classes drop
  from 13 to 11 at blast 10. Genuine topology errors, recorded per blast:
  ``iterated element ... is not inside its parent`` on a degenerate
  ``[7.37881, 7.37881]`` element (k28 / k32, blasts 9-14), ``Holes of origin
  ... disagree on bridge_side`` (k28 at 1e-5, blasts 16-18), ``Bridge ...
  approaches its two defining crossings from opposite sides`` (nested_inner,
  ``num_iterations = 2``, blasts 16-20). No law violation anywhere.

Record schema
-------------
metrics: ``case``, ``area_cutoff``, ``min_separation``, ``num_iterations``,
``build_s``, ``blasts_done``, ``stop_reason`` (``max_blasts`` | ``budget`` |
``memory_budget``),
``topology_errors``, ``law_violations``, ``topology_changes`` (blasts with any
change flag), ``structure_changes`` (blasts whose canonical words changed),
``first_change_blast``, ``letters_final``, ``n_classes_final``.
series (index 0 = before any blast, one entry per blast after):
numerics ``blast``, ``blast_s``, ``crossings``,
``new_crossings``, ``bridges``, ``new_bridges``, ``interior_bridges``,
``unstable_points``, ``stable_points``, ``bridge_points``, ``max_abs``,
``rss_mb``; laws ``same_stability``, ``cdist_decreasing``,
``cdist_worst_drop``, ``lobe_n``, ``lobe_ratio_max``, ``lobe_ratio_p50``,
``lobe_sign_flips``; topology scalars ``repin_s`` (``classify_strong_pips`` +
the pip restore, inside the topology guard; None when it raised),
``topology_s``, ``topology_outcome``
(ok | exception | invariant_assert), ``n_classes``, ``n_active``,
``n_inert``, ``n_letters``, ``is_reliable``, ``n_unresolved``,
``n_virtual``, ``n_ambiguous``, ``n_holes``, ``n_homotopy_elements``,
``n_iterated_elements``, ``minimal_kept``, ``minimal_image``,
``topology_hash``, ``canonical_hash``; change flags vs the last successful
snapshot (None at blast 0 and when this snapshot failed) ``classes_added``,
``classes_removed``, ``words_changed`` (classes present in both whose full
word differs, compared per class since inert letters are re-dealt every
snapshot), ``new_letters`` (counts),
``matrix_changed``, ``reliability_changed``, ``structure_changed``
(booleans); ``topology`` (one dict per blast: ``error`` (the traceback only at
a message's first occurrence), ``pips``, ``default_pips`` (what
``classify_strong_pips`` chose before the restore), ``classes``
``[label, letter, inert]``, ``rules``, ``refined_rules``, ``canonical_rules``,
``matrix`` / ``refined_matrix`` ``{names, n, nnz, rows, hash}`` (rows
omitted past 40 symbols), ``unresolved`` ``[letter, reason]``, ``virtual``,
``ambiguous``, and the change lists ``classes_added``, ``classes_removed``,
``words_changed``, ``new_letters``; lists capped at 60, words at 300 chars).
"""

from __future__ import annotations

import hashlib
import json
import time
import traceback

import numpy as np

import harness
from maps import curve_arrays, lobe_polygon, max_abs_coord, saddle_session, total_points

EXPERIMENT = "e11_blast_depth"
TIMEOUT = 900
WORKERS = 8
RSS_CAP_GB = 12

BUILD_CUTOFF = 1e-7
#: The separation each recipe blasts with (build_k28, build_period3, build_nested).
CANONICAL_SEPARATION = {"k28": 1e-5, "k32": 1e-5, "p3": 1e-5, "nested_outer": 1e-4,
                        "nested_inner": 1e-4}
CAP = 60
WORD_CHARS = 300
_SEEN_ERRORS: set[str] = set()


def configs(quick: bool) -> list[dict]:
    if quick:
        return [
            {"case": "k28", "max_blasts": 3, "area_cutoff": 1e-7, "min_separation": 1e-5,
             "num_iterations": 1},
            {"case": "p3", "max_blasts": 2, "area_cutoff": 1e-7, "min_separation": 1e-5,
             "num_iterations": 1},
        ]

    def cell(case, *, area_cutoff=1e-7, min_separation="canonical", num_iterations=1,
             repeat=0):
        separation = CANONICAL_SEPARATION[case] if min_separation == "canonical" else min_separation
        return {"case": case, "max_blasts": 20, "area_cutoff": area_cutoff,
                "min_separation": separation, "num_iterations": num_iterations,
                "repeat": repeat, "_timeout": TIMEOUT, "_budget": TIMEOUT - 120,
                "_rss_budget_gb": 0.9 * RSS_CAP_GB}

    grid = []
    for case in CANONICAL_SEPARATION:
        separations = (None, 1e-6, 1e-5, 1e-4, 1e-3) + ((1e-2,) if case in ("k28", "k32") else ())
        grid += [cell(case, min_separation=s) for s in separations]
        grid += [cell(case, area_cutoff=c) for c in (1e-5, 1e-6, 1e-8)]
        grid += [cell(case, num_iterations=2)]
    grid += [cell(case, repeat=r) for case in ("k28", "p3", "nested_inner") for r in (1, 2)]
    # Longest first (the budget-bound coarse cutoffs, then the OOM-bound cases),
    # so the many second-long cells fill in behind them.
    return sorted(grid, key=lambda c: (c["case"] in ("p3", "nested_inner"),
                                       c["area_cutoff"] != 1e-5))


# --------------------------------------------------------------------------- #
# Building the cases
# --------------------------------------------------------------------------- #
def build(case: str):
    """Session, fixed points to repin, their pips, and the zone to blast."""
    from tanglepack.examples import henon_cases as hc

    if case in ("k28", "k32"):
        session, fp = saddle_session(2.8 if case == "k28" else 3.2, 1, area_cutoff=BUILD_CUTOFF)
        session.grow_n_times(fp, "unstable", num_iterations=10)
        session.grow_until_turnaround(fp, "stable")
        hc._intersect_and_bridge(session, [fp])
        session.classify_strong_pips()
        trellis = session.trellis(fp)
        pip = trellis.iterate(trellis.strong_pip, 1)
        zone = session.resonance_zone(pip)
        hc._repin(session, [fp], [pip])
        return session, [fp], [pip], zone

    session = hc._session()
    fp3 = hc._grow_period3(session, 13, 9, BUILD_CUTOFF)
    if case == "p3":
        hc._intersect_and_bridge(session, [fp3])
        session.classify_strong_pips()
        pip = session.trellis(fp3).strong_pip
        zone = session.resonance_zone(pip)
        hc._repin(session, [fp3], [pip])
        return session, [fp3], [pip], zone

    fp1 = hc._grow_period1(session, 11, BUILD_CUTOFF)
    fixed_points = [fp1, fp3]
    hc._intersect_and_bridge(session, [fp3, fp1])
    session.classify_strong_pips()
    pips = [session.trellis(fp).strong_pip for fp in fixed_points]
    session.add_resonance_zones(pips)
    zones = sorted(session.resonance_zones.values(), key=lambda zone: zone.area)
    inner, outer = zones[0], zones[-1]
    hc._repin(session, fixed_points, pips)
    if case == "nested_outer":
        return session, fixed_points, pips, outer
    for _ in range(2):
        session.blast_zone(outer, num_iterations=1, fixed_point=[outer.fixed_point],
                           min_separation=1e-4)
        hc._repin(session, fixed_points, pips)
    return session, fixed_points, pips, inner


def repin(session, fixed_points, pips) -> list:
    """Classify, restore the chosen pips; return the default pips the classification chose."""
    session.classify_strong_pips()
    defaults = [session.strong_pip(fp) for fp in fixed_points]
    for fp, pip in zip(fixed_points, pips):
        session.set_strong_pip(fp, pip)
    return defaults


# --------------------------------------------------------------------------- #
# Numerics and the laws
# --------------------------------------------------------------------------- #
def numerics_snapshot(session, new_bridges: list) -> dict:
    wb = session.workbench
    return {
        "crossings": len(wb.intersection_registry),
        "bridges": len(wb.bridges),
        "new_bridges": len(new_bridges),
        "bridge_points": sum(len(b.get_point_array()) for b in wb.bridges),
        "unstable_points": total_points(session, "unstable"),
        "stable_points": total_points(session, "stable"),
    }


def law_checks(session, new_bridges: list, curves: dict, sample: int = 100) -> dict:
    """Same-stability crossings, cdist monotonicity of the new bridges, lobe areas.

    A new bridge ``(a, b)`` whose ends have registered preimages ``(a', b')``
    bounds the image of the lobe ``(a', b')``: the two polygons must have the
    same area (area preservation), up to the polygonal approximation.
    """
    wb = session.workbench
    registry, table = wb.intersection_registry, wb.intersection_registry.iterate_table
    same_stability = sum(
        ix.manifold_a_key[1] == ix.manifold_b_key[1] for _id, ix in registry
    )
    drops = [np.diff(np.asarray(b.get_cdist_array(), dtype=float).ravel()) for b in new_bridges]
    drops = np.concatenate(drops) if drops else np.zeros(0)
    ratios, flips = [], 0
    for bridge in [b for b in new_bridges if b.id is not None][:sample]:
        a, b = bridge.id
        pa, pb = table[a, -1], table[b, -1]
        if pa is None or pb is None:
            continue
        lobe, preimage = lobe_polygon(wb, curves, a, b), lobe_polygon(wb, curves, pa, pb)
        if lobe is None or preimage is None:
            continue
        area, pre_area = _signed_area(lobe), _signed_area(preimage)
        if pre_area:
            ratios.append(abs(area / pre_area - 1))
            flips += int(np.sign(area) != np.sign(pre_area))
    return {
        "same_stability": int(same_stability),
        "cdist_decreasing": int((drops < 0).sum()),
        "cdist_worst_drop": float(-drops.min()) if len(drops) and drops.min() < 0 else 0.0,
        "lobe_n": len(ratios),
        "lobe_ratio_max": max(ratios) if ratios else None,
        "lobe_ratio_p50": float(np.median(ratios)) if ratios else None,
        "lobe_sign_flips": flips,
    }


def _signed_area(polygon: np.ndarray) -> float:
    from tanglepack.numerics.geometry import signed_polygon_area

    return float(signed_polygon_area(polygon))


# --------------------------------------------------------------------------- #
# The topology snapshot
# --------------------------------------------------------------------------- #
def _words(rules: dict) -> dict:
    return {name: " ".join(s.text for s in word)[:WORD_CHARS]
            for name, word in list(rules.items())[:CAP]}


def _canonical_rules(dynamics, classes: list) -> dict:
    """The unrefined rules with every letter renamed by its class's table position.

    The session alphabet never reuses a letter, so a class whose elements are
    renumbered by a blast reappears under a new letter; renaming active classes
    ``A0, A1, ...`` and inert ones ``I0, I1, ...`` in table order separates a
    relabelling from a change of the symbolic structure itself.
    """
    canonical, counts = {}, {True: 0, False: 0}
    for _label, letter, inert in classes:
        if letter and letter not in canonical:
            canonical[letter] = f"{'I' if inert else 'A'}{counts[inert]}"
            counts[inert] += 1

    def rename(token: str) -> str:
        base, _, power = token.partition("^")
        return canonical.get(base, base) + (f"^{power}" if power else "")

    return {canonical.get(letter, letter): " ".join(rename(s.text) for s in word)
            for letter, word in dynamics.rules.items()}


def _matrix(dynamics, refined: bool) -> dict:
    names, matrix = dynamics.transition_matrix(refined=refined)
    out = {"names": names[:CAP], "n": len(names), "nnz": int(np.count_nonzero(matrix))}
    if len(names) <= 40:
        out["rows"] = matrix.tolist()
    out["hash"] = hashlib.sha1(json.dumps([names, matrix.tolist()]).encode()).hexdigest()[:12]
    return out


def topology_snapshot(session, fixed_points) -> dict:
    """The topology stack on the current state, as small JSON."""
    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()
    table = session.bridge_classes()
    dynamics = session.symbolic_dynamics()
    minimal = session.minimal_trellis()
    iterated = session.iterated_partition()
    trellises = [session.trellis(fp) for fp in fixed_points]
    letters = {cls: cd.letter for cls, cd in dynamics.classes.items()}
    classes = [[e.bridge_class.label, letters.get(e.bridge_class), e.inert] for e in table]
    rules, refined_rules = _words(dynamics.rules), _words(dynamics.refined_rules)
    # Compared in full (the stored words are capped): per class, because an inert
    # letter is re-dealt (u, v, ...) at every snapshot while an active one persists.
    words = {cls.label: [cd.letter, hashlib.sha1(cd.word.encode()).hexdigest()[:12]]
             for cls, cd in dynamics.classes.items()}
    canonical = _canonical_rules(dynamics, classes)
    matrix, refined_matrix = _matrix(dynamics, False), _matrix(dynamics, True)
    return {
        "classes": classes[:CAP],
        "class_labels": [label for label, _l, _i in classes],
        "words": words,
        "letters": sorted({letter for _c, letter, _i in classes if letter}),
        "rules": rules,
        "canonical_hash": hashlib.sha1(json.dumps(sorted(canonical.items())).encode()
                                       ).hexdigest()[:12],
        "canonical_rules": dict(list(canonical.items())[:CAP]),
        "refined_rules": refined_rules,
        "matrix": matrix,
        "refined_matrix": refined_matrix,
        "is_reliable": bool(dynamics.is_reliable),
        "unresolved": [[cd.letter, (cd.unresolved_reason or "")[:120]]
                       for cd in dynamics.unresolved][:CAP],
        "virtual": list(dynamics.virtual_classes.values())[:CAP],
        "ambiguous": [cd.letter for cd in dynamics.classes.values() if cd.ambiguous][:CAP],
        "n_active": len(table.active),
        "n_inert": len(table.inert),
        "n_holes": sum(len(t.holes) for t in trellises),
        "n_homotopy_elements": sum(len(r.intervals) for t in trellises
                                   for r in t.stable_partitions),
        "n_iterated_elements": sum(len(r.intervals) for r in iterated.as_list()),
        "minimal_kept": len(minimal.kept_bridge_ids),
        "minimal_image": len(minimal.image_bridge_ids),
        "hash": hashlib.sha1(json.dumps(
            [sorted(classes, key=str), sorted(words.items()), matrix["hash"],
             refined_matrix["hash"]]
        ).encode()).hexdigest()[:12],
    }


def changes(previous: dict | None, current: dict | None) -> dict:
    """What this blast uncovered relative to the last successful snapshot."""
    if previous is None or current is None:
        return {}
    before, after = set(previous["class_labels"]), set(current["class_labels"])
    words_before, words_after = previous["words"], current["words"]
    return {
        "classes_added": sorted(after - before)[:CAP],
        "classes_removed": sorted(before - after)[:CAP],
        "words_changed": sorted(
            words_after[label][0] for label in set(words_before) & set(words_after)
            if words_before[label][1] != words_after[label][1]
        )[:CAP],
        "matrix_changed": (previous["matrix"]["hash"] != current["matrix"]["hash"]
                           or previous["refined_matrix"]["hash"]
                           != current["refined_matrix"]["hash"]),
        "new_letters": sorted(set(current["letters"]) - set(previous["letters"])),
        "reliability_changed": previous["is_reliable"] != current["is_reliable"],
        "structure_changed": previous["canonical_hash"] != current["canonical_hash"],
    }


def record_topology(rec, session, fixed_points, pips, previous):
    """Repin and snapshot the topology into the series; a failure is recorded, never raised.

    The repin (``classify_strong_pips``) is topology too, so it sits inside the
    guard: a strong-pip failure is a topology error of this blast, not of the run.
    """
    t0 = time.perf_counter()
    snapshot, outcome, error, defaults = None, "ok", None, None
    try:
        defaults = repin(session, fixed_points, pips)
        rec.append("repin_s", time.perf_counter() - t0)
        t0 = time.perf_counter()
        snapshot = topology_snapshot(session, fixed_points)
    except Exception as exc:  # noqa: BLE001 - the failure is the datum
        outcome = "invariant_assert" if isinstance(exc, AssertionError) else "exception"
        error = f"{type(exc).__name__}: {exc}"[:500]
        if error not in _SEEN_ERRORS:  # the traceback once per distinct message
            _SEEN_ERRORS.add(error)
            error += "\n" + traceback.format_exc()[-1500:]
    if defaults is None:
        rec.append("repin_s", None)
    rec.append("topology_s", time.perf_counter() - t0)
    rec.append("topology_outcome", outcome)
    delta = changes(previous, snapshot)
    snap = snapshot or {}
    for name in ("n_active", "n_inert", "is_reliable", "n_holes", "n_homotopy_elements",
                 "n_iterated_elements", "minimal_kept", "minimal_image"):
        rec.append(name, snap.get(name))
    rec.append("n_classes", len(snap["class_labels"]) if snapshot else None)
    rec.append("n_letters", len(snap["letters"]) if snapshot else None)
    for name in ("unresolved", "virtual", "ambiguous"):
        rec.append(f"n_{name}", len(snap[name]) if snapshot else None)
    rec.append("topology_hash", snap.get("hash"))
    for name in ("classes_added", "classes_removed", "words_changed", "new_letters"):
        rec.append(name, len(delta[name]) if delta else None)
    for name in ("matrix_changed", "reliability_changed", "structure_changed"):
        rec.append(name, delta.get(name))
    rec.append("canonical_hash", snap.get("canonical_hash"))
    detail = {k: v for k, v in snap.items()
              if k not in ("class_labels", "words", "hash", "canonical_hash")}
    rec.append("topology", {"error": error, "pips": pips, "default_pips": defaults,
                            **detail, **{k: v for k, v in delta.items()
                                         if isinstance(v, list)}})
    return snapshot, delta


# --------------------------------------------------------------------------- #
# The run
# --------------------------------------------------------------------------- #
def run(config: dict, rec: harness.Recorder) -> None:
    t_start = time.perf_counter()
    for name in ("case", "area_cutoff", "min_separation", "num_iterations"):
        rec.metric(name, config[name])
    with rec.stage("build"):
        session, fixed_points, pips, zone = build(config["case"])
    rec.metric("build_s", rec.stages["build"])
    session.workbench._man_machine.area_cutoff = config["area_cutoff"]
    wb = session.workbench
    budget = config.get("_budget", float("inf"))

    known = {id(b) for b in wb.bridges}
    curves, curve_size = curve_arrays(session), total_points(session)
    max_abs = max_abs_coord(session)
    previous, _ = record_topology(rec, session, fixed_points, pips, None)
    blast_step = dict(blast=0, blast_s=0.0, new_crossings=0, interior_bridges=0)
    totals = dict(topology_errors=0, law_violations=0, topology_changes=0, structure_changes=0)
    first_change, stop_reason, last_cost = None, "max_blasts", 0.0

    def append_numerics(new_bridges: list) -> None:
        for name, value in {**blast_step, **numerics_snapshot(session, new_bridges),
                            **law_checks(session, new_bridges, curves)}.items():
            rec.append(name, value)
        rec.append("max_abs", max_abs)
        rec.append("rss_mb", rec.peak_rss_mb)

    append_numerics([])
    rec.checkpoint()
    rss_cap_mb = config.get("_rss_budget_gb", float("inf")) * 1024
    for blast in range(1, config["max_blasts"] + 1):
        # Graceful timeout: stop before a blast the budget cannot hold (cost grows ~2x).
        if time.perf_counter() - t_start + 2 * last_cost > budget:
            stop_reason = "budget"
            break
        # Likewise before a blast that would cross the memory cap: the peak grows by
        # the last blast's factor (~1.6x), and a watchdog kill is not a finding.
        peaks = rec.series["rss_mb"]
        growth = max(peaks[-1] / peaks[-2], 1.0) if len(peaks) > 1 and peaks[-2] else 1.0
        if peaks[-1] * growth > rss_cap_mb:
            stop_reason = "memory_budget"
            break
        t_blast = time.perf_counter()
        crossings_before = len(wb.intersection_registry)
        result = session.blast_zone(zone, num_iterations=config["num_iterations"],
                                    fixed_point=[zone.fixed_point],
                                    min_separation=config["min_separation"])
        t1 = time.perf_counter()
        new_bridges = [b for b in wb.bridges if id(b) not in known]
        known |= {id(b) for b in new_bridges}
        if total_points(session) != curve_size:
            curves, curve_size = curve_arrays(session), total_points(session)
        for b in new_bridges:
            points = b.get_point_array()
            if len(points):
                max_abs = max(max_abs, float(np.nanmax(np.abs(points))))
        blast_step = dict(blast=blast, blast_s=t1 - t_blast,
                          new_crossings=len(wb.intersection_registry) - crossings_before,
                          interior_bridges=len(result.all_interior_bridges()))
        append_numerics(new_bridges)
        snapshot, delta = record_topology(rec, session, fixed_points, pips, previous)
        changed = any(bool(v) for v in delta.values())
        previous = snapshot or previous
        totals["topology_errors"] += snapshot is None
        totals["law_violations"] += rec.series["same_stability"][-1] > 0
        totals["law_violations"] += rec.series["cdist_decreasing"][-1] > 0
        totals["topology_changes"] += changed
        totals["structure_changes"] += bool(delta.get("structure_changed"))
        if changed and first_change is None:
            first_change = blast
        rec.metric("blasts_done", blast)
        rec.metrics.update(totals)
        rec.metric("first_change_blast", first_change)
        rec.checkpoint()
        last_cost = time.perf_counter() - t_blast

    rec.metric("stop_reason", stop_reason)
    rec.metrics.update(totals)
    if previous is not None:
        rec.metric("letters_final", previous["letters"])
        rec.metric("n_classes_final", len(previous["class_labels"]))


if __name__ == "__main__":
    harness.main(globals())
