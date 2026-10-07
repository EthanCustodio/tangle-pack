"""Run the whole numerics stress campaign.

Phase A (timing, serial, GPU never shared): e01, e02, e03, e04, e10.
Phase B (robustness sweeps, each script's own worker pool): e05..e09.

    python stress/numerics/run_campaign.py [--quick] [--only e01,e07] [--skip e03]

Each experiment is run as its own CLI (``eNN_*.py [--quick]``) so a broken
script never stops the campaign; progress lines go to results/progress.log.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

import harness

ROOT = Path(__file__).resolve().parent

TIMING = ["e01", "e02", "e03", "e04", "e10"]
SWEEPS = ["e05", "e06", "e07", "e08", "e09"]


def _script(prefix: str) -> Path:
    matches = sorted(ROOT.glob(f"{prefix}_*.py"))
    if not matches:
        raise FileNotFoundError(f"no script for {prefix}")
    return matches[0]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--only", type=str, default=None)
    parser.add_argument("--skip", type=str, default="")
    args = parser.parse_args()

    order = TIMING + SWEEPS
    if args.only:
        order = [p for p in order if p in args.only.split(",")]
    order = [p for p in order if p not in args.skip.split(",")]

    t0 = time.time()
    harness.log_progress(f"=== campaign start quick={args.quick} order={order}")
    for prefix in order:
        script = _script(prefix)
        cmd = [sys.executable, str(script)] + (["--quick"] if args.quick else [])
        t1 = time.time()
        code = subprocess.call(cmd, cwd=harness.REPO)
        harness.log_progress(
            f"=== {script.name} exit={code} in {(time.time() - t1) / 60:.1f} min"
        )
    harness.log_progress(f"=== campaign done in {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
