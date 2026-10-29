"""The k25 entropy sweep's resume rule (``sweep_done``) and its summary row, on hand-written records."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts"))

import henon_k25_entropy as entropy  # noqa: E402
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
