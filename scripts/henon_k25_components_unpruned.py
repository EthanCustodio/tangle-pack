"""
The k=2.5 components figures from an UNPRUNED blast series.

The entropy sweep (``henon_k25_entropy.py``) prunes every blast with
``min_separation=1e-5``, which stops its frontier from about blast 8 on. This
driver reruns the same sweep with ``min_separation=None`` (and without the
blast's same-stability count, ``BLAST_DIAGNOSE=False``, the cost that grows
fastest), then draws every figure of ``henon_k25_components.py`` from it.

Writes ``figures/henon_k25_entropy/components/unpruned/``: the components,
polynomials (+ ``polynomials.txt``), spectra_by_blast and convergence figures,
with the sweep's data in ``data/``.

Run one sweep per process, then the figures::

    PYTHONPATH=src env/bin/python scripts/henon_k25_components_unpruned.py --cases p1 --max-rss-gb 60
    ...
    PYTHONPATH=src env/bin/python scripts/henon_k25_components_unpruned.py --plot-only

(``--cutoff 1e-9``, ``--max-blasts 15``.)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))

import henon_k25_components as components  # noqa: E402
import henon_k25_entropy as entropy  # noqa: E402
from symbolic_figures import configure_logging  # noqa: E402

OUT_DIR = components.OUT_DIR / "unpruned"


def install() -> None:
    """Unpruned blasts, no same-stability count, every path in the subfolder."""
    entropy.MIN_SEPARATION = None
    entropy.BLAST_DIAGNOSE = False
    entropy.DATA_DIR = OUT_DIR / "data"
    components.OUT_DIR = OUT_DIR


def main() -> None:
    parser = argparse.ArgumentParser(description="The k=2.5 components figures from an unpruned blast series.")
    parser.add_argument("--cases", nargs="+", choices=entropy.CASES, default=list(entropy.CASES))
    parser.add_argument("--cutoff", type=float, default=1e-9)
    parser.add_argument("--max-blasts", type=int, default=15)
    parser.add_argument("--max-rss-gb", type=float, default=None)
    parser.add_argument("--plot-only", action="store_true")
    args = parser.parse_args()
    configure_logging()
    plt.switch_backend("Agg")
    install()
    if not args.plot_only:
        for case in args.cases:
            entropy.sweep(case, args.cutoff, args.max_blasts, None, args.max_rss_gb)
        return
    sys.argv = [sys.argv[0], "--cutoff", f"{args.cutoff:g}", "--cases", *args.cases]
    components.main()


if __name__ == "__main__":
    main()
