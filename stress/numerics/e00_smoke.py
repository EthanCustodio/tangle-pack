"""Harness smoke test: one ok k=10 growth, one forced exception, one timeout."""
import time
import harness
from maps import saddle_session, manifold_sizes, total_points

EXPERIMENT = "e00_smoke"


def configs(quick):
    return [{"mode": "ok", "steps": 5, "gpu": False}, {"mode": "ok", "steps": 5, "gpu": True},
            {"mode": "raise"}, {"mode": "hang", "_timeout": 3}]


def run(config, rec):
    if config["mode"] == "raise":
        raise RuntimeError("boom")
    if config["mode"] == "hang":
        time.sleep(60)
    with rec.stage("setup"):
        s, fp = saddle_session(10, gpu=config["gpu"])
    for i in range(config["steps"]):
        with rec.stage("grow"):
            s.grow_n_times(fp, "unstable", num_iterations=1)
        rec.append("points", total_points(s, "unstable"))
    rec.metric("sizes", manifold_sizes(s))


if __name__ == "__main__":
    harness.main(globals())
