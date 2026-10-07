"""Shared harness for the numerics stress-test campaign.

Every experiment script is a thin module that defines

    EXPERIMENT = "e03_scaling"
    def configs(quick: bool) -> list[dict]: ...
    def run(config: dict, rec: Recorder) -> None: ...

and ends with ``if __name__ == "__main__": harness.main(globals())``.

``main`` gives every script the same CLI:

    python stress/numerics/eXX_*.py [--quick] [--workers N] [--timeout S]
                                    [--rss-cap-gb G] [--out PATH] [--only I,J]

The PARENT process launches one fresh subprocess per config (``--worker``), so
a crash, a hang or a memory blow-up of one config never takes the sweep down.
A psutil watchdog kills a worker whose resident set exceeds the cap (we cannot
use RLIMIT_AS: CUDA reserves huge virtual ranges). Each finished config is one
JSON line in ``results/<experiment>.jsonl``:

    {experiment, index, config, outcome, error, traceback_tail, wall_s,
     stages: {name: seconds}, metrics: {...}, series: {...},
     warnings: {"logger | message template": count}, peak_rss_mb,
     gpu_peak_mb, git_sha, host, started}

``outcome`` is one of ``ok | timeout | oom | exception | invariant_assert |
crash``. Failures are data: the sweep always continues.
"""

from __future__ import annotations

import argparse
import contextlib
import gc
import json
import logging
import os
import platform
import subprocess
import sys
import tempfile
import threading
import time
import traceback
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Callable, Iterator, Optional

import numpy as np
import psutil

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent.parent
RESULTS = ROOT / "results"
PROGRESS_LOG = RESULTS / "progress.log"
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(ROOT))


# --------------------------------------------------------------------------- #
# In-worker recording
# --------------------------------------------------------------------------- #
class _CountingHandler(logging.Handler):
    """Counts WARNING+ records by (logger, unformatted message template)."""

    def __init__(self) -> None:
        super().__init__(level=logging.WARNING)
        self.counts: dict[str, int] = {}

    def emit(self, record: logging.LogRecord) -> None:
        key = f"{record.name} | {str(record.msg)[:120]}"
        self.counts[key] = self.counts.get(key, 0) + 1


def _jsonable(value: Any) -> Any:
    """Convert numpy scalars/arrays (recursively) into JSON-friendly values."""
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return _jsonable(value.tolist())
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        f = float(value)
        return f if np.isfinite(f) else repr(f)
    if isinstance(value, float) and not np.isfinite(value):
        return repr(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def gpu_sync() -> None:
    """Block until the GPU is idle (no-op without CuPy)."""
    try:
        import cupy as cp  # type: ignore

        cp.cuda.runtime.deviceSynchronize()
    except Exception:
        pass


class Recorder:
    """Collects timings, metrics and per-step series for one config.

    Use ``with rec.stage("grow"):`` for a timed stage (stages with the same
    name accumulate), ``rec.metric(name, value)`` for scalars and
    ``rec.append(series, value)`` for per-step series (e.g. points per
    iteration). Timed stages synchronise the GPU at both ends.
    """

    def __init__(self) -> None:
        self.stages: dict[str, float] = {}
        self.metrics: dict[str, Any] = {}
        self.series: dict[str, list] = {}
        self.handler = _CountingHandler()
        self._peak_rss = 0
        self._stop = threading.Event()
        self._proc = psutil.Process()
        self._thread = threading.Thread(target=self._sample, daemon=True)

    def _sample(self) -> None:
        while not self._stop.is_set():
            with contextlib.suppress(Exception):
                self._peak_rss = max(self._peak_rss, self._proc.memory_info().rss)
            self._stop.wait(0.05)

    def start(self) -> None:
        """Attach the warning counter and start the RSS sampler."""
        logging.getLogger().addHandler(self.handler)
        logging.getLogger().setLevel(logging.WARNING)
        self._thread.start()

    def stop(self) -> None:
        """Stop the sampler and detach the warning counter."""
        self._stop.set()
        self._thread.join(timeout=1)
        logging.getLogger().removeHandler(self.handler)

    @contextlib.contextmanager
    def stage(self, name: str, sync: bool = True) -> Iterator[None]:
        """Time a block; accumulates into ``stages[name]``."""
        if sync:
            gpu_sync()
        t0 = time.perf_counter()
        try:
            yield
        finally:
            if sync:
                gpu_sync()
            self.stages[name] = self.stages.get(name, 0.0) + time.perf_counter() - t0

    def metric(self, name: str, value: Any) -> None:
        """Record one scalar (or small JSON-able) metric."""
        self.metrics[name] = value

    def append(self, series: str, value: Any) -> None:
        """Append one value to a per-step series."""
        self.series.setdefault(series, []).append(value)

    @property
    def peak_rss_mb(self) -> float:
        """Peak resident set observed so far, in MB."""
        with contextlib.suppress(Exception):
            self._peak_rss = max(self._peak_rss, self._proc.memory_info().rss)
        return self._peak_rss / 2**20


def _gpu_peak_mb() -> Optional[float]:
    if "cupy" not in sys.modules:
        return None
    try:
        import cupy as cp  # type: ignore

        return cp.get_default_memory_pool().total_bytes() / 2**20
    except Exception:
        return None


def _worker(module_globals: dict, config: dict, out_path: str) -> None:
    """Run one config in this process and dump the partial record to out_path."""
    rec = Recorder()
    rec.start()
    gc.collect()
    outcome, error, tb = "ok", None, None
    t0 = time.perf_counter()

    def dump() -> None:
        record = {
            "outcome": outcome,
            "error": error,
            "traceback_tail": tb,
            "wall_s": time.perf_counter() - t0,
            "stages": rec.stages,
            "metrics": rec.metrics,
            "series": rec.series,
            "warnings": rec.handler.counts,
            "peak_rss_mb": rec.peak_rss_mb,
            "gpu_peak_mb": _gpu_peak_mb(),
        }
        tmp = out_path + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(_jsonable(record), fh)
        os.replace(tmp, out_path)

    # Partial results survive a kill: the parent reads the last checkpoint.
    rec.checkpoint = dump  # type: ignore[attr-defined]
    try:
        module_globals["run"](config, rec)
    except AssertionError as exc:
        outcome, error = "invariant_assert", f"AssertionError: {exc}"[:500]
        tb = traceback.format_exc()[-3000:]
    except MemoryError as exc:
        outcome, error = "oom", f"MemoryError: {exc}"[:500]
        tb = traceback.format_exc()[-3000:]
    except BaseException as exc:  # noqa: BLE001 - failures are data
        outcome, error = "exception", f"{type(exc).__name__}: {exc}"[:500]
        tb = traceback.format_exc()[-3000:]
    rec.stop()
    dump()


# --------------------------------------------------------------------------- #
# Parent: sweep driver
# --------------------------------------------------------------------------- #
def _git_sha() -> str:
    with contextlib.suppress(Exception):
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=REPO, text=True
        ).strip()
    return "unknown"


def log_progress(line: str) -> None:
    """Append a timestamped line to results/progress.log and echo it."""
    RESULTS.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%H:%M:%S")
    with open(PROGRESS_LOG, "a") as fh:
        fh.write(f"{stamp} {line}\n")
    print(f"{stamp} {line}", flush=True)


def _run_one(
    script: str,
    experiment: str,
    index: int,
    config: dict,
    timeout: float,
    rss_cap_gb: float,
) -> dict:
    """Launch one worker subprocess, police it, and return its full record."""
    with tempfile.TemporaryDirectory() as tmp:
        out_path = os.path.join(tmp, "record.json")
        cmd = [
            sys.executable,
            script,
            "--worker",
            "--config",
            json.dumps(config),
            "--record",
            out_path,
        ]
        started = time.strftime("%Y-%m-%dT%H:%M:%S")
        t0 = time.perf_counter()
        proc = subprocess.Popen(
            cmd,
            cwd=REPO,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
            env={**os.environ, "PYTHONPATH": str(REPO / "src")},
        )
        stderr_chunks: list[str] = []
        reader = threading.Thread(
            target=lambda: stderr_chunks.append(proc.stderr.read()), daemon=True
        )
        reader.start()
        ps = psutil.Process(proc.pid)
        killed: Optional[str] = None
        peak = 0
        while proc.poll() is None:
            with contextlib.suppress(psutil.Error):
                peak = max(peak, ps.memory_info().rss)
                if ps.memory_info().rss > rss_cap_gb * 2**30:
                    killed = "oom"
            if time.perf_counter() - t0 > timeout:
                killed = "timeout"
            if killed:
                proc.kill()
                proc.wait()
                break
            time.sleep(0.05)
        reader.join(timeout=5)
        wall = time.perf_counter() - t0
        record: dict = {}
        if os.path.exists(out_path):
            with open(out_path) as fh:
                record = json.load(fh)
        if killed:
            record["outcome"] = killed
            record["error"] = (
                f"killed by watchdog after {wall:.1f}s"
                if killed == "timeout"
                else f"RSS exceeded {rss_cap_gb} GB"
            )
        elif not record:
            record = {
                "outcome": "crash",
                "error": f"worker exited {proc.returncode} without a record",
                "traceback_tail": "".join(stderr_chunks)[-3000:],
            }
        record.setdefault("peak_rss_mb", peak / 2**20)
        record["peak_rss_mb"] = max(record.get("peak_rss_mb") or 0, peak / 2**20)
        record.update(
            {
                "experiment": experiment,
                "index": index,
                "config": config,
                "process_wall_s": wall,
                "git_sha": _git_sha(),
                "host": platform.node(),
                "started": started,
            }
        )
        return record


def sweep(
    script: str,
    experiment: str,
    configs: list[dict],
    *,
    out: Path,
    workers: int = 1,
    timeout: float = 600,
    rss_cap_gb: float = 24,
) -> list[dict]:
    """Run every config in its own worker and append records to ``out``."""
    out.parent.mkdir(parents=True, exist_ok=True)
    lock = threading.Lock()
    records: list[dict] = []
    log_progress(f"[{experiment}] start: {len(configs)} configs, workers={workers}")

    def task(item: tuple[int, dict]) -> None:
        index, config = item
        cfg_timeout = float(config.get("_timeout", timeout))
        record = _run_one(script, experiment, index, config, cfg_timeout, rss_cap_gb)
        with lock:
            with open(out, "a") as fh:
                fh.write(json.dumps(record) + "\n")
            records.append(record)
            log_progress(
                f"[{experiment}] {len(records)}/{len(configs)} #{index} "
                f"{record['outcome']} {record.get('process_wall_s', 0):.1f}s "
                f"rss={record.get('peak_rss_mb', 0):.0f}MB "
                f"{json.dumps({k: v for k, v in config.items() if not k.startswith('_')})[:140]}"
                + (f" ERR {record.get('error')}"[:200] if record["outcome"] != "ok" else "")
            )

    with ThreadPoolExecutor(max_workers=workers) as pool:
        list(pool.map(task, list(enumerate(configs))))
    outcomes: dict[str, int] = {}
    for r in records:
        outcomes[r["outcome"]] = outcomes.get(r["outcome"], 0) + 1
    log_progress(f"[{experiment}] done: {outcomes}")
    return records


def main(module_globals: dict) -> None:
    """The shared CLI of every experiment script (see module docstring)."""
    parser = argparse.ArgumentParser(description=module_globals.get("__doc__"))
    parser.add_argument("--quick", action="store_true", help="tiny grid, seconds")
    parser.add_argument("--workers", type=int, default=module_globals.get("WORKERS", 1))
    parser.add_argument("--timeout", type=float, default=module_globals.get("TIMEOUT", 600))
    parser.add_argument("--rss-cap-gb", type=float, default=module_globals.get("RSS_CAP_GB", 24))
    parser.add_argument("--out", type=str, default=None)
    parser.add_argument("--only", type=str, default=None, help="comma list of config indices")
    parser.add_argument("--list", action="store_true", help="print configs and exit")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--config", type=str, help=argparse.SUPPRESS)
    parser.add_argument("--record", type=str, help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.worker:
        _worker(module_globals, json.loads(args.config), args.record)
        return

    experiment = module_globals["EXPERIMENT"]
    configs = module_globals["configs"](args.quick)
    if args.only:
        keep = {int(i) for i in args.only.split(",")}
        configs = [c for i, c in enumerate(configs) if i in keep]
    if args.list:
        for i, c in enumerate(configs):
            print(i, json.dumps(c))
        return
    suffix = "_quick" if args.quick else ""
    out = Path(args.out) if args.out else RESULTS / f"{experiment}{suffix}.jsonl"
    if out.exists():
        out.rename(out.with_suffix(f".jsonl.bak{int(time.time())}"))
    sweep(
        str(Path(module_globals["__file__"]).resolve()),
        experiment,
        configs,
        out=out,
        workers=args.workers,
        timeout=args.timeout,
        rss_cap_gb=args.rss_cap_gb,
    )


def load(experiment: str, quick: bool = False) -> list[dict]:
    """Read every record of one experiment's results file."""
    path = RESULTS / f"{experiment}{'_quick' if quick else ''}.jsonl"
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
