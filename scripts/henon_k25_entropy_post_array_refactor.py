"""
The k=2.5 entropy figures rerun after the numerics arrays campaign (main 0d890b8 on).

Every k=2.5 entropy script, unchanged, pointed at
``figures/henon_k25_entropy_post_array_refactor/``: the arrays campaign took
the 1e-7 series from 1914 to 232 MiB held, so the sweeps run to finer cutoffs
and further blasts than ``figures/henon_k25_entropy/``. A pruned sweep
(``min_separation`` 1e-5) runs until its frontier is exhausted (the first
blast that keeps no new interior bridge in any zone), or until its time or
memory budget runs out.

Subcommands (one process each; ``PYTHONPATH=src env/bin/python scripts/henon_k25_entropy_post_array_refactor.py``):

* ``sweep --cases p1 --cutoffs 1e-10`` -- the pruned entropy sweep (``data/``);
* ``pool --cases p1 p3 --cutoffs 1e-9 1e-10 --processes 4 --workers 4`` -- one
  ``sweep`` job per case and cutoff in ``scripts/sweep_pool.py`` (``pool/``:
  status and logs), merged into ``data/sweeps.csv``; ``--max-rss-gb`` and
  ``--max-seconds`` are the sweep's own soft stop (before a blast), and the
  pool's hard cap (kills the job) is the same memory and ``HARD_TIME_FACTOR``
  times the time, so the soft time stop fires first; a rerun skips every
  finished job (one whose JSON records a stop or its last blast) and every
  capped one unless ``--retry``;
* ``unpruned --cases p1 --cutoff 1e-9`` -- the unpruned series of
  ``henon_k25_components_unpruned.py`` (``components/unpruned/data/``);
* ``length ...`` -- ``henon_k25_length_growth.py`` (``length_growth/``), its arguments passed on;
* ``compare ...`` -- ``henon_k25_length_compare.py``, its arguments passed on;
* ``panels ...`` -- ``henon_k25_blast_panels.py`` (``blast_panels/``), its arguments passed on;
* ``plot --cutoff 1e-9`` -- every sweep figure, ``components/`` at that cutoff,
  ``components/finest/`` with each case at its finest exhausted sweep and
  ``components/unpruned/`` from the data on disk.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))

import henon_k25_blast_panels as blast_panels  # noqa: E402
import henon_k25_components as components  # noqa: E402
import henon_k25_components_unpruned as unpruned  # noqa: E402
import henon_k25_entropy as entropy  # noqa: E402
import henon_k25_length_compare as length_compare  # noqa: E402
import henon_k25_length_growth as length_growth  # noqa: E402
from sweep_pool import Job, run_pool  # noqa: E402
from symbolic_figures import FIGURES_DIR, configure_logging  # noqa: E402

OUT_DIR = FIGURES_DIR / "henon_k25_entropy_post_array_refactor"
#: The pool's hard time cap over ``--max-seconds``: the sweep only checks its clock between
#: blasts, and the last blast of a series can take about as long as all before it.
HARD_TIME_FACTOR = 2.0


def install() -> None:
    """Point every script at this run's folder and stop sweeps on an exhausted frontier."""
    entropy.OUT_DIR = OUT_DIR
    entropy.DATA_DIR = OUT_DIR / "data"
    entropy.STOP_WHEN_EXHAUSTED = True
    components.OUT_DIR = OUT_DIR / "components"
    unpruned.OUT_DIR = components.OUT_DIR / "unpruned"
    blast_panels.OUT_DIR = OUT_DIR / "blast_panels"
    length_growth.OUT_DIR = OUT_DIR / "length_growth"
    length_growth.DATA_DIR = length_growth.OUT_DIR / "data"


def finest_exhausted() -> dict[str, float]:
    """Per case, the finest cutoff whose pruned sweep ran until its frontier was exhausted."""
    best: dict[str, float] = {}
    for (case, cutoff), record in entropy.load_records().items():
        if (record["stopped"] or "").startswith("frontier exhausted") and cutoff < best.get(case, math.inf):
            best[case] = cutoff
    return best


class _Label(float):
    """A float that formats as a fixed text: the components figures print their cutoff with ``:g``."""

    def __new__(cls, text: str) -> "_Label":
        obj = super().__new__(cls, 0.0)
        obj.text = text
        return obj

    def __format__(self, spec: str) -> str:
        return f"finest per case: {self.text}"


def components_at_finest() -> None:
    """``components/finest/``: the components figures with each case read at its finest exhausted sweep."""
    finest = {case: cutoff for case, cutoff in sorted(finest_exhausted().items(), key=lambda kv: entropy.CASES.index(kv[0]))}
    label = _Label(", ".join(f"{case} {cutoff:g}" for case, cutoff in finest.items()))
    results = {}
    for case, cutoff in finest.items():
        cache: dict = {}
        results[case] = [components.analyse(row, cache) for row in entropy.load_records()[(case, cutoff)]["rows"]]
    components.OUT_DIR = OUT_DIR / "components" / "finest"
    components.OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(components.write_polynomials(results, label).read_text())
    for draw in (components.figure_components, components.figure_polynomials, components.figure_spectra,
                 components.figure_convergence):
        draw(results, label)


def sweep_done(case: str, cutoff: float, max_blasts: int) -> bool:
    """A sweep is finished when its JSON records a stop or a row at its last blast (an unreadable one is not)."""
    path = entropy.data_path(case, cutoff)
    try:
        record = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return False
    return bool(record["stopped"] or (record["rows"] and record["rows"][-1]["blasts"] >= max_blasts))


def pool(args: argparse.Namespace) -> None:
    """One ``sweep`` job per (case, cutoff) in a :func:`sweep_pool.run_pool`, merged into ``data/sweeps.csv``."""
    limits = [*(["--max-seconds", f"{args.max_seconds:g}"] if args.max_seconds is not None else []),
              "--max-rss-gb", f"{args.max_rss_gb:g}"]
    jobs = [
        Job(f"{case}_a{cutoff:.0e}",
            [sys.executable, str(Path(__file__).resolve()), "sweep", "--cases", case, "--cutoffs", f"{cutoff:g}",
             "--max-blasts", str(args.max_blasts), *limits],
            lambda case=case, cutoff=cutoff: sweep_done(case, cutoff, args.max_blasts),
            csv=entropy.data_path(case, cutoff).with_suffix(".csv"))
        for cutoff in args.cutoffs for case in args.cases
    ]
    run_pool(jobs, processes=args.processes, workers=args.workers, out=OUT_DIR / "pool",
             max_rss_gb=args.max_rss_gb,
             max_seconds=None if args.max_seconds is None else HARD_TIME_FACTOR * args.max_seconds,
             merged_csv=entropy.DATA_DIR / "sweeps.csv", retry=args.retry)


def main() -> None:
    command, rest = (sys.argv[1], sys.argv[2:]) if len(sys.argv) > 1 else ("plot", [])
    configure_logging()
    plt.switch_backend("Agg")
    install()
    if command in ("length", "compare", "panels"):
        module = {"length": length_growth, "compare": length_compare, "panels": blast_panels}[command]
        sys.argv = [sys.argv[0], *rest]
        module.main()
        return
    parser = argparse.ArgumentParser(description="The k=2.5 entropy figures, post arrays campaign.")
    parser.add_argument("--cases", nargs="+", choices=entropy.CASES, default=list(entropy.CASES))
    parser.add_argument("--cutoffs", nargs="+", type=float, default=[1e-7])
    parser.add_argument("--cutoff", type=float, default=1e-9)
    parser.add_argument("--max-blasts", type=int, default=400)
    parser.add_argument("--max-seconds", type=float, default=None)
    parser.add_argument("--max-rss-gb", type=float, default=64.0)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--processes", type=int, default=1)
    parser.add_argument("--retry", action="store_true")
    args = parser.parse_args(rest)
    if command == "sweep":
        for cutoff in args.cutoffs:
            for case in args.cases:
                entropy.sweep(case, cutoff, args.max_blasts, args.max_seconds, args.max_rss_gb, args.workers)
    elif command == "pool":
        pool(args)
    elif command == "unpruned":
        unpruned.install()
        for case in args.cases:
            entropy.sweep(case, args.cutoff, args.max_blasts, args.max_seconds, args.max_rss_gb, args.workers)
    elif command == "plot":
        entropy.draw_figures()
        sys.argv = [sys.argv[0], "--cutoff", f"{args.cutoff:g}"]
        components.main()
        components_at_finest()
        unpruned.install()
        entropy.OUT_DIR = OUT_DIR
        sys.argv = [sys.argv[0], "--cutoff", f"{args.cutoff:g}"]
        components.main()
    else:
        raise SystemExit(f"unknown command {command!r}")


if __name__ == "__main__":
    main()
