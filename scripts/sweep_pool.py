"""
A small process pool for independent sweep jobs: memory and time caps, resume, one merged CSV.

The author runs many sweeps (k x cutoff x case) side by side, each one process
of a script. :func:`run_pool` launches them ``processes`` at a time, appends
``--workers W`` to every command (so ``processes x workers`` cores are busy,
never more than ``cores``), samples each job's whole process tree once a
second and kills a tree past ``max_rss_gb`` or ``max_seconds``. Every job ends
with one row of ``out/pool_status.csv`` (``done``, ``failed``, ``rss_cap``,
``time_cap`` or ``skipped``), written as it finishes, so a rerun resumes:

* a job whose ``done()`` is True is ``skipped`` without running;
* a job the status file records as ended by a cap keeps that row and is not
  relaunched (a job that hits a memory cap would hit it again) unless
  ``retry`` is set;
* every other job runs, a ``failed`` one included.

Resume is per finished job: a killed job starts over (its script may resume
inside, as the zone sweep does per ``k``). The pool lives in ``scripts/`` on
purpose: the library has no scheduler (D6).

Dev Notes:
    * A job's command is run as given; its stdout and stderr go to
      ``out/logs/<name>.log``, overwritten on every run.
    * Ctrl-C kills every running tree, records it ``failed`` and re-raises.
"""

from __future__ import annotations

import csv
import os
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Callable, Optional

import psutil

#: Seconds between two samples of the running jobs.
POLL_SECONDS = 1.0
STATUS_COLUMNS = ["name", "status", "returncode", "seconds", "peak_rss_gb"]
CAPS = ("rss_cap", "time_cap")


@dataclass(frozen=True)
class Job:
    """
    One independent command of a sweep.

    Attributes:
        name: Unique; the key of the job's status row and its log file.
        argv: The command; the pool appends ``["--workers", str(workers)]``.
        done: True when the job's output is complete on disk (resume).
        csv: The rows this job writes, merged into the pool's one CSV.
    """

    name: str
    argv: list[str]
    done: Callable[[], bool]
    csv: Optional[Path] = None


def read_status(out: Path) -> dict[str, dict[str, str]]:
    """The rows of ``out/pool_status.csv`` by job name (empty without the file)."""
    path = out / "pool_status.csv"
    if not path.exists():
        return {}
    with path.open(newline="") as handle:
        return {row["name"]: row for row in csv.DictReader(handle)}


def write_status(out: Path, rows: dict[str, dict[str, str]]) -> None:
    """Rewrite ``out/pool_status.csv`` (written whole, after every finished job)."""
    with (out / "pool_status.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=STATUS_COLUMNS)
        writer.writeheader()
        writer.writerows(rows.values())


def tree_rss_gb(process: psutil.Process) -> float:
    """The resident memory of a process and all its descendants, in GB."""
    total = 0
    for member in [process, *process.children(recursive=True)]:
        try:
            total += member.memory_info().rss
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
    return total / 1e9


def kill_tree(process: psutil.Process) -> None:
    """Kill a process and every descendant, and reap the process."""
    try:
        members = [*process.children(recursive=True), process]
    except psutil.NoSuchProcess:
        members = [process]
    for member in members:
        try:
            member.kill()
        except psutil.NoSuchProcess:
            pass
    psutil.wait_procs(members, timeout=5)


def merge_csvs(paths: list[Path], merged: Path) -> int:
    """Concatenate the existing CSVs in order over the union of their columns; returns the row count."""
    columns: dict[str, None] = {}
    rows: list[dict[str, str]] = []
    for path in paths:
        if not path.exists():
            continue
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            columns.update(dict.fromkeys(reader.fieldnames or []))
            rows.extend(reader)
    merged.parent.mkdir(parents=True, exist_ok=True)
    with merged.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), restval="")
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


@dataclass
class _Running:
    """A launched job and what the pool samples of it."""

    job: Job
    popen: subprocess.Popen
    process: psutil.Process
    started: float
    log: IO[str]
    peak_gb: float = 0.0


def run_pool(
    jobs: list[Job],
    *,
    processes: int,
    workers: int = 1,
    out: Path,
    max_rss_gb: Optional[float] = None,
    max_seconds: Optional[float] = None,
    merged_csv: Optional[Path] = None,
    retry: bool = False,
    cores: int = os.cpu_count() or 1,
) -> list[dict[str, str]]:
    """
    Run the jobs ``processes`` at a time under the caps, resuming from ``out``.

    Args:
        jobs: The jobs, in the order they are launched and merged.
        processes: Jobs running at once.
        workers: Appended to every command as ``--workers``.
        out: The folder of ``pool_status.csv`` and ``logs/``.
        max_rss_gb: Kill a job whose process tree holds more resident memory.
        max_seconds: Kill a job running longer.
        merged_csv: Where every existing ``job.csv`` is concatenated.
        retry: Relaunch the jobs a previous run ended by a cap.
        cores: The cores ``processes * workers`` may use.

    Returns:
        The status row of every job, in job order.

    Raises:
        ValueError: ``processes * workers`` exceeds ``cores``, or two jobs
            share a name.
    """
    if processes * workers > cores:
        raise ValueError(f"{processes} processes x {workers} workers exceed the {cores} cores")
    if len({job.name for job in jobs}) != len(jobs):
        raise ValueError("job names must be unique")
    (out / "logs").mkdir(parents=True, exist_ok=True)
    status = read_status(out)

    def record(job: Job, outcome: str, returncode: object = "", seconds: float = 0.0, peak: float = 0.0) -> None:
        status[job.name] = {"name": job.name, "status": outcome, "returncode": str(returncode),
                            "seconds": f"{seconds:.1f}", "peak_rss_gb": f"{peak:.3f}"}
        write_status(out, status)
        print(f"[pool] {job.name}: {outcome}" + (f" ({seconds:.0f} s, {peak:.2f} GB)" if seconds else ""), flush=True)

    pending: list[Job] = []
    for job in jobs:
        if job.done():
            record(job, "skipped")
        elif status.get(job.name, {}).get("status") in CAPS and not retry:
            print(f"[pool] {job.name}: {status[job.name]['status']} before, not relaunched (retry relaunches)", flush=True)
        else:
            pending.append(job)

    running: list[_Running] = []
    try:
        while pending or running:
            while pending and len(running) < processes:
                job = pending.pop(0)
                log = (out / "logs" / f"{job.name}.log").open("w")
                popen = subprocess.Popen([*job.argv, "--workers", str(workers)], stdout=log, stderr=subprocess.STDOUT)
                running.append(_Running(job, popen, psutil.Process(popen.pid), time.monotonic(), log))
            time.sleep(POLL_SECONDS)
            for run in list(running):
                seconds = time.monotonic() - run.started
                returncode = run.popen.poll()
                if returncode is None:
                    run.peak_gb = max(run.peak_gb, tree_rss_gb(run.process))
                    cap = ("rss_cap" if max_rss_gb is not None and run.peak_gb > max_rss_gb
                           else "time_cap" if max_seconds is not None and seconds > max_seconds else None)
                    if cap is None:
                        continue
                    kill_tree(run.process)
                    returncode, outcome = run.popen.wait(), cap
                else:
                    outcome = "done" if returncode == 0 else "failed"
                run.log.close()
                running.remove(run)
                record(run.job, outcome, returncode, seconds, run.peak_gb)
    except KeyboardInterrupt:
        for run in running:
            kill_tree(run.process)
            run.log.close()
            record(run.job, "failed", run.popen.wait(), time.monotonic() - run.started, run.peak_gb)
        raise

    if merged_csv is not None:
        rows = merge_csvs([job.csv for job in jobs if job.csv is not None], merged_csv)
        print(f"[pool] {rows} rows merged into {merged_csv}", flush=True)
    return [status[job.name] for job in jobs if job.name in status]
