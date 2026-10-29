"""The zone sweep's parallel path (``--processes`` / ``--job``) with a stubbed pool and frame builder: no tangle is built."""

from __future__ import annotations

import csv
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts"))

import henon_period3_zone_sweep as zone  # noqa: E402

#: Two runs from the seed k = 2: one k below it, three above.
RANGE = ["--k-min", "1.9", "--k-max", "2.3", "--k-step", "0.1"]
KS = [1.9, 2.0, 2.1, 2.2, 2.3]


def run_main(monkeypatch, argv: list[str]) -> None:
    monkeypatch.setattr(sys, "argv", ["henon_period3_zone_sweep.py", *argv])
    zone.main()


def pool_jobs(monkeypatch, tmp_path: Path, *options: str) -> tuple[list, dict]:
    """The jobs and keywords ``--processes`` hands the pool, which runs nothing."""
    seen: dict = {}
    monkeypatch.setattr(zone, "run_pool", lambda jobs, **kwargs: seen.update(jobs=jobs, kwargs=kwargs))
    monkeypatch.setattr(zone, "draw_summary", lambda rows, path: None)
    for k in KS:  # rows on disk, so the merge has something to read
        write_row(tmp_path, k)
    run_main(monkeypatch, [*RANGE, "--out", str(tmp_path), *options])
    return seen["jobs"], seen["kwargs"]


def write_row(out: Path, k: float) -> None:
    zone.row_path(out, k).parent.mkdir(exist_ok=True)
    zone.write_csv([{"k": f"{k:.4f}"}], zone.row_path(out, k))


def k_range(job) -> tuple[float, float]:
    return float(job.argv[job.argv.index("--k-min") + 1]), float(job.argv[job.argv.index("--k-max") + 1])


@pytest.mark.parametrize("processes", ["1", "2", "4"])
def test_every_k_and_the_seed_is_in_exactly_one_job(monkeypatch, tmp_path, processes):
    jobs, _ = pool_jobs(monkeypatch, tmp_path, "--processes", processes)
    for k in KS:
        assert sum(low <= k <= high for low, high in map(k_range, jobs)) == 1, k
    for job in jobs:  # no chunk straddles the seed: each lies in one run of sweep_order
        low, high = k_range(job)
        assert high < zone.SEED_K or low >= zone.SEED_K


def test_a_job_forwards_the_users_options_but_not_the_pools(monkeypatch, tmp_path):
    jobs, kwargs = pool_jobs(monkeypatch, tmp_path, "--processes", "2", "--workers", "3", "--max-seconds", "9",
                             "--max-rss-gb", "2", "--retry", "--blasts", "2", "--no-fill")
    for job in jobs:
        assert job.argv.count("--job") == 1
        assert job.argv[job.argv.index("--blasts") + 1] == "2" and "--no-fill" in job.argv
        assert job.argv[job.argv.index("--out") + 1] == str(tmp_path)
        for flag in ("--processes", "--workers", "--max-seconds", "--max-rss-gb", "--retry"):
            assert flag not in job.argv, flag
    assert kwargs == dict(processes=2, workers=3, out=tmp_path, max_rss_gb=2.0, max_seconds=9.0, retry=True)


def test_a_job_is_done_once_each_k_has_a_row_or_a_skip_record(monkeypatch, tmp_path):
    jobs, _ = pool_jobs(monkeypatch, tmp_path, "--processes", "1")
    for k in KS:
        zone.row_path(tmp_path, k).unlink()
    assert not any(job.done() for job in jobs)
    for k in KS[:-1]:
        write_row(tmp_path, k)
    zone.skip_path(tmp_path, KS[-1]).write_text("ValueError: no orbit\n")
    assert all(job.done() for job in jobs)


def test_the_merge_sorts_the_rows_by_k(monkeypatch, tmp_path):
    pool_jobs(monkeypatch, tmp_path, "--processes", "2")
    with (tmp_path / "zones.csv").open(newline="") as handle:
        assert [row["k"] for row in csv.DictReader(handle)] == [f"{k:.4f}" for k in KS]


class FakeBuilder:
    """Stands in for ``build_zone_frame``: records ``(k, blasts)`` and fails at ``fail_at``."""

    def __init__(self, fail_at: float = None):
        self.calls: list[tuple[float, int]] = []
        self.fail_at = fail_at

    def __call__(self, k, guess, orientation, wanted, **growth):
        self.calls.append((round(k, 10), growth["blasts"]))
        if k == self.fail_at:
            raise ValueError("no orbit")
        pip = SimpleNamespace(label=lambda: "q", point=(0.0, 0.0))
        frame = SimpleNamespace(k=k, pip=pip, pip_jump=0.0, zone_area=1.0, inside=[], mean_area=0.5,
                                lobe_over_zone=0.5, blast_sizes=[], blast_error=None, unstable_steps=1, stable_steps=1)
        return frame, guess, orientation


def run_job(monkeypatch, tmp_path: Path, builder: FakeBuilder) -> list[float]:
    """One ``--job`` run over 2.0..2.3 at two blasts; returns the ``k`` whose outputs it wrote."""
    written: list[float] = []

    def write_job_frame(frame, args):
        written.append(frame.k)
        write_row(args.out, frame.k)

    monkeypatch.setattr(zone, "build_zone_frame", builder)
    monkeypatch.setattr(zone, "write_job_frame", write_job_frame)
    run_main(monkeypatch, ["--job", "--k-min", "2.0", "--k-max", "2.3", "--blasts", "2", "--out", str(tmp_path)])
    return written


def test_a_job_writes_the_seed_row_and_blasts_every_k(monkeypatch, tmp_path):
    builder = FakeBuilder()
    assert run_job(monkeypatch, tmp_path, builder) == [2.0, 2.1, 2.2, 2.3]
    assert builder.calls == [(2.0, 2), (2.1, 2), (2.2, 2), (2.3, 2)]


def test_a_resumed_job_warm_steps_a_finished_k_and_blasts_the_missing_ones(monkeypatch, tmp_path):
    for k in (2.0, 2.1):
        write_row(tmp_path, k)
    builder = FakeBuilder()
    assert run_job(monkeypatch, tmp_path, builder) == [2.2, 2.3]
    assert builder.calls[1:] == [(2.1, 0), (2.2, 2), (2.3, 2)]


def test_a_failed_k_is_recorded_and_not_retried(monkeypatch, tmp_path):
    assert run_job(monkeypatch, tmp_path, FakeBuilder(fail_at=2.2)) == [2.0, 2.1, 2.3]
    assert zone.skip_path(tmp_path, 2.2).read_text().startswith("ValueError: no orbit")

    zone.row_path(tmp_path, 2.3).unlink()
    builder = FakeBuilder(fail_at=2.2)
    assert run_job(monkeypatch, tmp_path, builder) == [2.3]
    assert builder.calls[1:] == [(2.1, 0), (2.3, 2)]
