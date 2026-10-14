"""Both open intervals of a pseudoneighbor pair (B2): what does the stable one kill?

The definition (``Pseudoneighbor_Algorithm.pdf``) keeps a consecutive window
pair only when NEITHER open interval, the unstable one nor the stable one,
holds a point of X. Until B2 the code tested the unstable interval alone and
took the stable one as empty by construction; ``_interval_is_clear(...,
stability, ...)`` now tests both. This experiment re-measures the E12 series
with the one function called twice per window pair:

* ``before`` -- the pairs whose UNSTABLE interval is clear (the old rule);
* ``after``  -- the pairs whose both intervals are clear (the library's
  references, checked against ``trellis.reference_pseudoneighbors``).

Every pair in ``before`` but not in ``after`` is a reference the stable clause
lost; it is recorded with the stable reading that fired (``table``, ``own
frame`` or ``scaled``, read off the function's debug line), the witness
crossing and its step count, and ``reference_witness`` of the law helpers
(None for a pair killed by scaling alone). The clause can only reject, never
add, so a loss is a finding to report, not a failure.

Cases and depths are E12's (``e12_reverify.DEPTH`` / ``SEPARATION``): k28 and
k32 to 13 blasts, p3 to 20, nested_outer and nested_inner to 12, inversion to
10, heteroclinic to 8, each at 1e-7 and 1e-5, one iteration per call: 14 cells.
E12's two-iteration cells and depth study are not rerun. Each step also
records the seed / closed / hole counts and N (``n_trellis_crossings``): once
N stops growing the later blasts measure the same N again (B3).

Report only: not wired into ``stress/numerics/report/``.
"""

from __future__ import annotations

import inspect
import logging
import sys
import time

import harness
import e11_blast_depth as e11
from e12_reverify import DEPTH, MAX_POINTS, SEPARATION

sys.path.insert(0, str(harness.REPO / "tests"))

from helpers.laws import reference_witness  # noqa: E402
from tanglepack.topology import Pseudoneighbor  # noqa: E402
from tanglepack.topology.Pseudoneighbor import _interval_is_clear, _reference_branch  # noqa: E402

EXPERIMENT = "e13_both_intervals"
TIMEOUT = 900
WORKERS = 8
RSS_CAP_GB = 12
COLLISION_RTOL = inspect.signature(
    Pseudoneighbor.compute_pseudoneighbors).parameters["collision_rtol"].default


def configs(quick: bool) -> list[dict]:
    budgets = {"_timeout": TIMEOUT, "_budget": TIMEOUT - 120, "_rss_budget_gb": 0.9 * RSS_CAP_GB}
    if quick:
        return [{"case": "k28", "max_blasts": 3, "area_cutoff": 1e-7},
                {"case": "p3", "max_blasts": 2, "area_cutoff": 1e-7}]
    cells = [{"case": case, "max_blasts": depth, "area_cutoff": cutoff, **budgets}
             for case, depth in DEPTH.items() for cutoff in (1e-7, 1e-5)]
    # The heavy period-1 series first, so the short cells fill in behind them.
    return sorted(cells, key=lambda c: c["case"] not in ("k28", "k32", "nested_outer"))


class _Puncture(logging.Handler):
    """Keeps the arguments of ``_interval_is_clear``'s puncture debug line:
    ``(x0, x1, stability, r_id, steps, reading, value, lo, hi)``."""

    def __init__(self) -> None:
        super().__init__(logging.DEBUG)
        self.args: tuple | None = None

    def emit(self, record: logging.LogRecord) -> None:
        if "interval punctured" in str(record.msg):
            self.args = record.args


def window_pairs(trellis):
    """Every consecutive pair of each fixed point's reference window, as
    ``compute_pseudoneighbors`` walks it."""
    for fixed_point in trellis.fixed_points:
        branch, r_n = _reference_branch(trellis, fixed_point)
        if branch is None or len(branch) < 2:
            continue
        r_end = trellis.iterate(r_n, branch.return_period)
        if r_end is None or trellis.intersection(r_end).manifold_b_key != branch.key:
            continue
        x0, x1 = r_n, trellis.next_intersection(r_n, "stable", toward_anchor=True)
        while x0 != r_end and x1 is not None:
            yield x0, x1
            x0, x1 = x1, trellis.next_intersection(x1, "stable", toward_anchor=True)


def measure(trellis, puncture: _Puncture) -> dict:
    """Before / after reference sets of one trellis and the pairs the stable clause lost."""
    candidates = trellis.own_intersection_ids
    before, after, lost, stable_fires = [], [], [], 0
    for x0, x1 in window_pairs(trellis):
        unstable_clear = _interval_is_clear(trellis, x0, x1, "unstable", candidates, COLLISION_RTOL)
        puncture.args = None
        stable_clear = _interval_is_clear(trellis, x0, x1, "stable", candidates, COLLISION_RTOL)
        pair = tuple(sorted((x0, x1)))
        stable_fires += not stable_clear
        if not unstable_clear:
            continue
        before.append(pair)
        if stable_clear:
            after.append(pair)
            continue
        _, _, _, witness, steps, reading, value, lo, hi = puncture.args or (None,) * 9
        lost.append({"pair": pair, "reading": reading, "witness": witness, "steps": steps,
                     "value": value, "interval": [lo, hi],
                     "reference_witness": reference_witness(trellis, pair)})
    library = sorted(p.as_tuple() for p in trellis.reference_pseudoneighbors)
    return {"before": sorted(before), "after": sorted(after), "lost": lost,
            "n_stable_fires": stable_fires, "matches_library": sorted(after) == library}


def record_step(rec: harness.Recorder, session, puncture: _Puncture) -> None:
    """One recompute (D1): the pair sets and the counts, a failure recorded, never raised."""
    t0 = time.perf_counter()
    try:
        topology = session.topology()
        trellis = topology.trellis
        values = {**measure(trellis, puncture),
                  "n_classes": len(topology.classes),
                  "n_closed": len(topology.dynamics.classes),
                  "n_holes": len(topology.holes),
                  "n_trellis_crossings": len(trellis.own_intersection_ids),
                  "topology_error": None}
    except Exception as exc:  # noqa: BLE001 - the failure is the datum
        values = {"topology_error": f"{type(exc).__name__}: {exc}"[:500]}
    values["topology_s"] = time.perf_counter() - t0
    for name in ("before", "after", "lost", "n_stable_fires", "matches_library", "n_classes",
                 "n_closed", "n_holes", "n_trellis_crossings", "topology_error", "topology_s"):
        rec.append(name, values.get(name))


def run(config: dict, rec: harness.Recorder) -> None:
    for name, value in config.items():
        if not name.startswith("_"):
            rec.metric(name, value)
    t_start = time.perf_counter()
    with rec.stage("build"):
        session, _fixed_points, _pips, zone = e11.build(config["case"])
    session.area_cutoff = config["area_cutoff"]
    puncture = _Puncture()
    logger = logging.getLogger(Pseudoneighbor.__name__)
    logger.addHandler(puncture)
    logger.setLevel(logging.DEBUG)

    rec.append("blast", 0)
    rec.append("blast_stop", None)
    record_step(rec, session, puncture)
    rec.checkpoint()
    stop_reason, last_cost = "max_blasts", 0.0
    budget = config.get("_budget", float("inf"))
    for blast in range(1, config["max_blasts"] + 1):
        if time.perf_counter() - t_start + 2 * last_cost > budget:
            stop_reason = "budget"
            break
        t_blast = time.perf_counter()
        result = session.blast_zone(zone, num_iterations=1, fixed_point=[zone.fixed_point],
                                    min_separation=SEPARATION[config["case"]],
                                    diagnose=False, max_points=MAX_POINTS)
        rec.append("blast", blast)
        rec.append("blast_stop", result.stop_reason)
        record_step(rec, session, puncture)
        rec.metric("blasts_done", blast)
        rec.checkpoint()
        last_cost = time.perf_counter() - t_blast
    rec.metric("stop_reason", stop_reason)


def summarise(records: list[dict]) -> None:
    """Print one line per cell, and every lost pair of the cells with any."""
    for record in records:
        config, series = record["config"], record["series"]
        cell = f"{config['case']:13s} {config['area_cutoff']:.0e}"
        lost = [(blast, entry) for blast, entries in zip(series.get("blast", []),
                                                         series.get("lost", []))
                for entry in entries or ()]
        errors = [e for e in series.get("topology_error", []) if e]
        mismatches = series.get("matches_library", []).count(False)
        print(f"{cell} {record['outcome']:9s} blasts {record['metrics'].get('blasts_done', 0):2d}"
              f" ({record['metrics'].get('stop_reason')}) lost {len(lost)}"
              f" stable fires {sum(n or 0 for n in series.get('n_stable_fires', []))}"
              f" library mismatches {mismatches} topology errors {len(errors)}")
        for blast, entry in lost:
            print(f"    blast {blast}: {entry}")


if __name__ == "__main__":
    harness.main(globals())
    if "--worker" not in sys.argv and "--list" not in sys.argv:
        summarise(harness.load(EXPERIMENT, quick="--quick" in sys.argv))
