"""Fixed-point solver robustness: random seeds, every period 1..6, four maps.

One config per ``(k, b, period)``; each solves ``n_seeds`` random seeds in two
variants:

* ``orbit`` -- one point uniform in the box and its first ``period - 1`` images
  under the map (a guess that already obeys the orbit's ordering; for longer
  periods it often escapes to huge or non-finite coordinates);
* ``independent`` -- ``period`` independent uniform points.

Each guess goes through ``FixedPointSolver.compute_fixed_point`` (the timed
fsolve shooting solve), and a converged orbit then through
``construct_fixed_point`` (eigendata, saddle validation, ``set_k_value``).
Splitting the two keeps the converged orbit of an elliptic or otherwise
rejected solution, so distinct orbits of every kind can be counted. The
residual ``max |f(x_{i-1}) - x_i|`` and the minimal period are measured here,
independently of the solver.

For ``k = 10`` the map is a complete binary horseshoe, so the number of orbits
of minimal period ``p`` is ``(1/p) sum_{d|p} mu(p/d) 2^d`` = 2, 1, 2, 3, 6, 9
(the ``horseshoe_orbits`` reference).
"""

from __future__ import annotations

import time

import numpy as np

import harness
from maps import henon_fixed_points
from tanglepack.examples.henon import henon_jacobian, henon_map, henon_map_inverse
from tanglepack.numerics.DynamicalSystem import DynamicalSystem
from tanglepack.numerics.FixedPointSolver import FixedPointSolver

EXPERIMENT = "e08_fixed_point_solver"
TIMEOUT = 900
WORKERS = 8

MAPS = [(10, 1), (2.8, 1), (2, 1), (2.1, 1)]
HORSESHOE_ORBITS = {1: 2, 2: 1, 3: 2, 4: 3, 5: 6, 6: 9}


def configs(quick: bool) -> list[dict]:
    n_seeds = 5 if quick else 200
    return [
        {"k": k, "b": b, "period": p, "n_seeds": n_seeds, "rng": 1000 * i + p}
        for i, (k, b) in enumerate(MAPS)
        for p in range(1, 7)
    ]


def box_half_width(k: float, b: float) -> float:
    """3 for the small maps, scaled with the outer fixed point for large k."""
    x_plus = henon_fixed_points(k, b)["saddle"][0]
    return max(3.0, 1.5 * abs(x_plus))


def orbit_guess(f, seed: np.ndarray, period: int) -> np.ndarray:
    points = [seed]
    for _ in range(period - 1):
        points.append(f(points[-1]))
    return np.array(points)


def residual(f, orbit: np.ndarray) -> float:
    images = np.array([f(x) for x in orbit])
    return float(np.max(np.abs(np.roll(images, 1, axis=0) - orbit)))


def minimal_period(orbit: np.ndarray, tol: float = 1e-8) -> int:
    period = len(orbit)
    for d in range(1, period + 1):
        if period % d == 0 and np.allclose(np.roll(orbit, -d, axis=0), orbit, atol=tol):
            return d
    return period


def canonical(orbit: np.ndarray, digits: int = 8) -> list:
    """The orbit's points of one minimal cycle, rotated to start at the smallest."""
    cycle = np.round(orbit[: minimal_period(orbit)], digits) + 0.0
    start = min(range(len(cycle)), key=lambda i: tuple(cycle[i]))
    return np.roll(cycle, -start, axis=0).tolist()


def failure_class(exc: Exception) -> str:
    message = str(exc)
    if "fsolve did not converge" in message:
        return "no_convergence_ier" + message.split("ier=")[1].split(")")[0]
    if "are complex" in message:
        return "elliptic"
    if "exactly one modulus" in message:
        return "not_hyperbolic"
    if "disagreeing" in message:
        return "sign_disagreement"
    return type(exc).__name__


def solve(solver: FixedPointSolver, f, guess: np.ndarray) -> dict:
    out: dict = {"guess_finite": bool(np.all(np.isfinite(guess)))}
    t0 = time.perf_counter()
    try:
        orbit = solver.compute_fixed_point(guess)
    except Exception as exc:  # noqa: BLE001 - the failure class is the datum
        out.update(solve_s=time.perf_counter() - t0, outcome=failure_class(exc),
                   error=f"{type(exc).__name__}: {exc}"[:200])
        return out
    out["solve_s"] = time.perf_counter() - t0
    out["residual"] = residual(f, orbit)
    out["minimal_period"] = minimal_period(orbit)
    out["orbit"] = canonical(orbit)
    t0 = time.perf_counter()
    try:
        fixed_point = solver.construct_fixed_point(orbit)
    except Exception as exc:  # noqa: BLE001
        out.update(construct_s=time.perf_counter() - t0, outcome=failure_class(exc),
                   error=f"{type(exc).__name__}: {exc}"[:200])
        return out
    lam_u = float(np.ravel(fixed_point.unstable_eigenvalues[0])[0])
    lam_s = float(np.ravel(fixed_point.stable_eigenvalues[0])[0])
    out.update(construct_s=time.perf_counter() - t0, outcome="saddle", lambda_u=lam_u,
               lambda_s=lam_s, det_defect=abs(lam_u * lam_s - 1.0),
               k_value=fixed_point.k_value, accuracy=float(fixed_point.accuracy))
    return out


def run(config: dict, rec: harness.Recorder) -> None:
    k, b, period = config["k"], config["b"], config["period"]
    f = henon_map(k, b)
    solver = FixedPointSolver(DynamicalSystem(f, henon_map_inverse(k, b), henon_jacobian(k, b)))
    rng = np.random.default_rng(config["rng"])
    half = box_half_width(k, b)
    rec.metric("box_half_width", half)

    with np.errstate(all="ignore"):
        for i in range(config["n_seeds"]):
            seed = rng.uniform(-half, half, size=2)
            guesses = {
                "orbit": orbit_guess(f, seed, period),
                "independent": rng.uniform(-half, half, size=(period, 2)),
            }
            for variant, guess in guesses.items():
                with rec.stage(f"solve_{variant}", sync=False):
                    result = solve(solver, f, guess)
                rec.append("solves", {"seed": i, "variant": variant, **result})
            if i % 20 == 19:
                rec.checkpoint()

    solves = rec.series.get("solves", [])
    outcomes: dict[str, int] = {}
    for s in solves:
        outcomes[f"{s['variant']}:{s['outcome']}"] = outcomes.get(f"{s['variant']}:{s['outcome']}", 0) + 1
    rec.metric("outcomes", outcomes)

    distinct: dict[str, set] = {}
    for s in solves:
        if "orbit" in s:
            label = f"p{s['minimal_period']}:{'saddle' if s['outcome'] == 'saddle' else s['outcome']}"
            distinct.setdefault(label, set()).add(str(s["orbit"]))
    rec.metric("distinct_orbits", {label: len(orbits) for label, orbits in distinct.items()})
    rec.metric("distinct_orbit_list", {label: sorted(orbits) for label, orbits in distinct.items()})
    converged = [s["residual"] for s in solves if "residual" in s]
    rec.metric("max_residual", max(converged) if converged else None)
    if k == 10:
        rec.metric("horseshoe_orbits", HORSESHOE_ORBITS[period])


if __name__ == "__main__":
    harness.main(globals())
