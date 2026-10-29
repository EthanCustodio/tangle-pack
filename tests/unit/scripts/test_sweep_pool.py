"""The sweep pool of ``scripts/sweep_pool.py``: caps, resume, merge, on tiny ``python -c`` jobs."""

from __future__ import annotations

import csv
import signal
import sys
from pathlib import Path

import psutil
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts"))

import sweep_pool  # noqa: E402
from sweep_pool import Job, read_status, run_pool  # noqa: E402

#: Allocates about 300 MB and holds it.
ALLOCATE = "x = bytearray(300_000_000); import time; time.sleep(30)"


@pytest.fixture(autouse=True)
def fast_sampling(monkeypatch):
    """Sample five times a second instead of once, so the suite stays quick."""
    monkeypatch.setattr(sweep_pool, "POLL_SECONDS", 0.2)


def python(code: str) -> list[str]:
    return [sys.executable, "-c", code]


def writer(path: Path, rows: list[dict[str, str]]) -> str:
    """Code that writes ``rows`` to ``path`` as a CSV."""
    columns = list(dict.fromkeys(key for row in rows for key in row))
    return (f"import csv; w = csv.DictWriter(open({str(path)!r}, 'w', newline=''), fieldnames={columns!r}); "
            f"w.writeheader(); w.writerows({rows!r})")


def marker_job(tmp_path: Path, name: str, code: str = "", **kwargs) -> Job:
    """A job that runs ``code`` and then leaves ``<name>.ran``; done once it exists."""
    marker = tmp_path / f"{name}.ran"
    return Job(name, python(f"{code}\nopen({str(marker)!r}, 'w').close()"), marker.exists, **kwargs)


def statuses(rows: list[dict[str, str]]) -> dict[str, str]:
    return {row["name"]: row["status"] for row in rows}


def test_merged_csv_holds_rows_in_job_order_over_the_union_of_columns(tmp_path):
    first, second = tmp_path / "first.csv", tmp_path / "second.csv"
    jobs = [
        Job("second", python(writer(second, [{"a": "3", "c": "x"}])), lambda: False, csv=second),
        Job("first", python(writer(first, [{"a": "1", "b": "2"}, {"a": "4", "b": "5"}])), lambda: False, csv=first),
    ]
    run_pool(jobs, processes=2, out=tmp_path, merged_csv=tmp_path / "merged.csv")
    with (tmp_path / "merged.csv").open(newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    assert reader.fieldnames == ["a", "c", "b"]
    assert rows == [{"a": "3", "c": "x", "b": ""}, {"a": "1", "c": "", "b": "2"}, {"a": "4", "c": "", "b": "5"}]


def test_done_jobs_are_skipped_and_a_rerun_runs_only_the_missing_ones(tmp_path):
    (tmp_path / "old.ran").touch()
    jobs = [marker_job(tmp_path, name) for name in ("old", "new")]
    rows = run_pool(jobs, processes=2, out=tmp_path)
    assert statuses(rows) == {"old": "skipped", "new": "done"}
    assert not (tmp_path / "logs" / "old.log").exists()

    (tmp_path / "new.ran").unlink()
    rows = run_pool(jobs, processes=2, out=tmp_path)
    assert statuses(rows) == {"old": "skipped", "new": "done"}
    assert statuses(read_status(tmp_path).values()) == {"old": "skipped", "new": "done"}


def test_a_job_past_the_memory_cap_is_killed_with_its_whole_tree(tmp_path):
    pid_file = tmp_path / "child.pid"
    parent = (f"import os, subprocess, sys, time\n"
              f"child = subprocess.Popen([sys.executable, '-c', {ALLOCATE!r}])\n"
              f"open({str(pid_file)!r}, 'w').write(str(child.pid))\n"
              f"time.sleep(30)")
    rows = run_pool([Job("tree", python(parent), lambda: False)], processes=1, out=tmp_path, max_rss_gb=0.2)
    assert statuses(rows) == {"tree": "rss_cap"}
    assert float(rows[0]["peak_rss_gb"]) > 0.2
    child = int(pid_file.read_text())
    try:
        assert psutil.Process(child).status() == psutil.STATUS_ZOMBIE
    except psutil.NoSuchProcess:
        pass


def test_a_job_past_the_time_cap_ends_as_time_cap(tmp_path):
    rows = run_pool([Job("sleeper", python("import time; time.sleep(30)"), lambda: False)],
                    processes=1, out=tmp_path, max_seconds=1.5)
    assert statuses(rows) == {"sleeper": "time_cap"}
    assert rows[0]["returncode"] == str(-signal.SIGKILL)
    assert float(rows[0]["seconds"]) < 10


def test_a_capped_job_is_not_relaunched_unless_retry(tmp_path):
    job = marker_job(tmp_path, "capped", "import time; time.sleep(30)")
    run_pool([job], processes=1, out=tmp_path, max_seconds=1.5)
    rows = run_pool([job], processes=1, out=tmp_path, max_seconds=1.5)
    assert statuses(rows) == {"capped": "time_cap"}
    assert not (tmp_path / "capped.ran").exists()

    quick = marker_job(tmp_path, "capped")
    rows = run_pool([quick], processes=1, out=tmp_path, retry=True)
    assert statuses(rows) == {"capped": "done"}
    assert (tmp_path / "capped.ran").exists()


def test_more_processes_times_workers_than_cores_is_refused(tmp_path):
    with pytest.raises(ValueError, match="cores"):
        run_pool([], processes=3, workers=2, out=tmp_path, cores=5)


def test_a_failing_command_is_recorded_and_the_others_still_run(tmp_path):
    jobs = [Job("broken", python("raise SystemExit(3)"), lambda: False), marker_job(tmp_path, "fine")]
    rows = run_pool(jobs, processes=1, out=tmp_path)
    assert statuses(rows) == {"broken": "failed", "fine": "done"}
    assert rows[0]["returncode"] == "3"


def test_the_pool_appends_the_workers(tmp_path):
    seen = tmp_path / "argv.txt"
    job = Job("argv", python(f"import sys; open({str(seen)!r}, 'w').write(' '.join(sys.argv[1:]))"), lambda: False)
    run_pool([job], processes=1, workers=2, out=tmp_path, cores=2)
    assert seen.read_text() == "--workers 2"
