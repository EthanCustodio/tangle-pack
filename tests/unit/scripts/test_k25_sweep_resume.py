"""The k25 entropy sweep's resume rule (``sweep_done``), its pool's jobs and its summary row, on hand-written records."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts"))

import henon_k25_entropy as entropy  # noqa: E402
import henon_k25_entropy_post_array_refactor as post  # noqa: E402
from henon_k25_entropy_post_array_refactor import sweep_done  # noqa: E402


@pytest.fixture(autouse=True)
def data_dir(tmp_path, monkeypatch):
    monkeypatch.setattr(entropy, "DATA_DIR", tmp_path)


def record(rows: list[dict], stopped=None) -> dict:
    return {"case": "p1", "area_cutoff": 1e-5, "rows": rows, "stopped": stopped}


def write(rec: dict) -> None:
    entropy.data_path(rec["case"], rec["area_cutoff"]).write_text(json.dumps(rec))


@pytest.mark.parametrize("rec, done", [
    (record([]), False),
    (record([{"blasts": 0}, {"blasts": 1}]), False),
    (record([{"blasts": 0}, {"blasts": 1}, {"blasts": 2}]), True),
    (record([], stopped="RuntimeError: before the first step"), True),
    (record([{"blasts": 0}], stopped="time limit (60 s) before blast 1"), True),
])
def test_a_sweep_is_done_on_a_stop_or_a_row_at_its_last_blast(rec, done):
    write(rec)
    assert sweep_done("p1", 1e-5, max_blasts=2) is done


def test_a_sweep_without_its_json_is_not_done():
    assert not sweep_done("p1", 1e-5, max_blasts=2)


def test_a_truncated_json_is_not_done():
    entropy.data_path("p1", 1e-5).write_text(json.dumps(record([], stopped="x"))[:-7])
    assert not sweep_done("p1", 1e-5, max_blasts=2)


def pool_jobs(monkeypatch, **options) -> tuple[list, dict]:
    """The jobs and keywords :func:`post.pool` hands :func:`sweep_pool.run_pool` (stubbed: nothing runs)."""
    seen: dict = {}
    monkeypatch.setattr(post, "run_pool", lambda jobs, **kwargs: seen.update(jobs=jobs, kwargs=kwargs))
    args = argparse.Namespace(cases=["p1", "p3"], cutoffs=[1e-5, 1e-6], max_blasts=2, max_seconds=None,
                              max_rss_gb=4.0, processes=2, workers=3, retry=False)
    vars(args).update(options)
    post.pool(args)
    return seen["jobs"], seen["kwargs"]


def test_the_pool_has_one_sweep_job_per_case_and_cutoff_each_reading_its_own_files(monkeypatch):
    jobs, kwargs = pool_jobs(monkeypatch)
    pairs = [("p1", 1e-5), ("p3", 1e-5), ("p1", 1e-6), ("p3", 1e-6)]
    assert len({job.name for job in jobs}) == len(jobs) == len(pairs)
    for job, (case, cutoff) in zip(jobs, pairs):
        command = job.argv[job.argv.index("sweep"):]
        assert command[:5] == ["sweep", "--cases", case, "--cutoffs", f"{cutoff:g}"]
        assert "--workers" not in command
        assert job.csv == entropy.data_path(case, cutoff).with_suffix(".csv")
    assert (kwargs["processes"], kwargs["workers"], kwargs["merged_csv"]) == (2, 3, entropy.DATA_DIR / "sweeps.csv")

    write({**record([], stopped="done"), "case": "p3", "area_cutoff": 1e-6})
    assert [job.done() for job in jobs] == [False, False, False, True]


def test_the_pool_time_cap_is_behind_the_sweep_soft_stop(monkeypatch):
    jobs, kwargs = pool_jobs(monkeypatch, max_seconds=60.0)
    assert all(job.argv[job.argv.index("--max-seconds") + 1] == "60" for job in jobs)
    assert kwargs["max_seconds"] == post.HARD_TIME_FACTOR * 60.0 > 60.0
    assert kwargs["max_rss_gb"] == 4.0
    assert pool_jobs(monkeypatch)[1]["max_seconds"] is None


def read_row(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        return reader.fieldnames, list(reader)


def test_the_summary_row_is_the_last_step(tmp_path):
    rows = [{"blasts": 0, "entropy": 0.1, "points": 10}, {"blasts": 1, "entropy": 0.2, "points": 20, "extra": 9}]
    entropy.write_summary_row(record(rows, stopped="frontier exhausted at blast 1"), tmp_path / "s.csv")
    columns, [row] = read_row(tmp_path / "s.csv")
    assert columns == entropy.SUMMARY_COLUMNS
    assert (row["case"], row["area_cutoff"], row["blasts"], row["entropy"], row["points"]) == ("p1", "1e-05", "1", "0.2", "20")
    assert row["stopped"] == "frontier exhausted at blast 1"
    assert row["reliable"] == ""


def test_the_summary_row_before_the_first_step_is_empty_but_its_keys(tmp_path):
    entropy.write_summary_row(record([]), tmp_path / "s.csv")
    _, [row] = read_row(tmp_path / "s.csv")
    assert row["case"] == "p1" and row["area_cutoff"] == "1e-05"
    assert all(row[key] == "" for key in entropy.SUMMARY_COLUMNS if key not in ("case", "area_cutoff"))
