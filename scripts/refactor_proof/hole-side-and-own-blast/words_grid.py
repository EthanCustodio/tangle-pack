"""
Proof of ``hole-side-and-own-blast``: the nested words against each tangle alone.

Two passes. The first runs inside one source tree (``--tree``, the branch by
default, the base worktree for the old code) and records every case as JSON:

* period 3 alone at 0..8 blasts and period 1 alone at the outer blast counts
  (the references, ``min_separation = 1e-4`` like the nested recipe);
* the nested tangle at every (outer, inner) cell of the grid, inner blasts
  first (the recipe), and at a few cells outer first too.

A case's words are spelled in element names (orbit code and short name, no
fixed-point letter) with directions, because letters differ between
sessions. The second pass (``--plot``) reads the base and branch JSONs and
draws ``nested_grid.png`` (each nested cell: does each tangle read exactly as
it does alone, is it reliable, does the blast order matter) and
``p3_blast_sweep.png`` (period 3 alone by blast count: holes, classes,
reliability, and whether the words still equal those at 0 blasts).

Run from the branch worktree::

    MPLBACKEND=Agg ../tangle-pack/env/bin/python scripts/refactor_proof/hole-side-and-own-blast/words_grid.py --tree ../tangle-pack-hole-side-and-own-blast-base --label base --blast-both
    MPLBACKEND=Agg ../tangle-pack/env/bin/python scripts/refactor_proof/hole-side-and-own-blast/words_grid.py --label branch
    MPLBACKEND=Agg ../tangle-pack/env/bin/python scripts/refactor_proof/hole-side-and-own-blast/words_grid.py --plot

``--blast-both`` makes the outer-first variant blast the inner zone the way
the base recipe did (every fixed point at once); the branch blasts a zone's
own fixed point only.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
BRANCH = HERE.parents[2]
OUT = BRANCH / "figures" / "refactor_proof" / "hole-side-and-own-blast"

OUTER = (0, 2, 4, 6, 8)
INNER = (0, 2, 4, 6)
ORDER_CELLS = ((2, 2), (4, 4), (8, 4))
P3_BLASTS = range(9)
MIN_SEPARATION = 1e-4


class _Warnings(logging.Handler):
    """Collects the library's WARNING messages of one case."""

    def __init__(self) -> None:
        super().__init__(logging.WARNING)
        self.messages: list[str] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.messages.append(record.getMessage())


def tangle_record(build, fixed_point) -> dict:
    """One tangle's classes, words (spelled) and reliability."""
    dynamics = build.session.symbolic_dynamics()

    def spelled(bridge_class) -> str:
        names = (dynamics.naming.homotopy_name(ref) for ref in (bridge_class.source, bridge_class.target))
        return "{" + ", ".join(f"{n.orbit_code}:{n.short_text}" for n in names) + "}"

    mine = [cd for cd in dynamics.classes.values() if cd.bridge_class.source.branch_key[0] is fixed_point]
    words = {
        spelled(cd.bridge_class): " ".join(
            spelled(s.bridge_class) + ("^-1" if s.direction < 0 else "") for s in cd.symbols
        )
        for cd in mine
    }
    reliable = all(
        cd.itinerary is not None and not cd.ambiguous and cd.verified is not False for cd in mine
    )
    trellis = build.session.trellis(fixed_point)
    return {
        "words": words,
        "letters": {cd.letter: cd.word for cd in mine if cd.kind == "active"},
        "reliable": reliable,
        "classes": len(mine),
        "active": sum(cd.kind == "active" for cd in mine),
        "holes": len(trellis.holes),
        "bridges": len(trellis.bridges),
    }


def run_case(name: str, build_fn, results: dict) -> None:
    """Build one case and record each of its tangles, or the exception it raised."""
    handler = _Warnings()
    logging.getLogger().addHandler(handler)
    start = time.time()
    try:
        build = build_fn()
        record = {
            "ok": True,
            "reliable": build.session.symbolic_dynamics().is_reliable,
            "tangles": [tangle_record(build, fp) for fp in build.fixed_points],
            "blast_sizes": build.blast_sizes,
        }
    except Exception as error:  # the old code raises; that is the evidence
        record = {"ok": False, "error": f"{type(error).__name__}: {error}"[:400]}
        traceback.print_exc(limit=1)
    finally:
        logging.getLogger().removeHandler(handler)
    record["seconds"] = round(time.time() - start, 2)
    record["skipped_bridges"] = sum("skipping bridge" in m for m in handler.messages)
    record["warnings"] = sorted({m[:90] for m in handler.messages})
    results[name] = record
    print(name, "ok" if record["ok"] else record["error"][:120], f"{record['seconds']}s", flush=True)


def collect(tree: Path, label: str, blast_both: bool) -> None:
    """Pass one: every case inside ``tree``, written to ``grid_<label>.json``."""
    sys.path.insert(0, str(tree / "src"))
    import tanglepack
    from tanglepack.examples import henon_cases as hc

    assert Path(tanglepack.__file__).is_relative_to(tree), tanglepack.__file__
    logging.basicConfig(level=logging.ERROR)
    logging.getLogger().setLevel(logging.WARNING)

    def period1(blasts: int):
        session = hc._session()
        fp1 = hc._grow_period1(session, 11, 1e-7)
        hc._intersect_and_bridge(session, [fp1])
        session.classify_strong_pips()
        pip = session.trellis(fp1).strong_pip
        zone = session.resonance_zone(pip)
        hc._repin(session, [fp1], [pip])
        for _ in range(blasts):
            session.blast_zone(zone, num_iterations=1, fixed_point=[fp1], min_separation=MIN_SEPARATION)
            hc._repin(session, [fp1], [pip])
        hc._partition(session)
        return hc.TangleBuild(session, [fp1], [pip], None, "period 1", "p1")

    def outer_first(outer: int, inner: int):
        build = hc.build_nested(outer_blasts=outer)
        session = build.session
        zone = min(session.resonance_zones.values(), key=lambda z: z.area)
        fixed_point = build.fixed_points if blast_both else [zone.fixed_point]
        for _ in range(inner):
            session.blast_zone(zone, 1, fixed_point=fixed_point, min_separation=MIN_SEPARATION)
            hc._repin(session, build.fixed_points, build.pips)
        hc._partition(session)
        return build

    results: dict = {}
    for blasts in P3_BLASTS:
        run_case(f"p3_{blasts}", lambda: hc.build_period3(blasts=blasts, min_separation=MIN_SEPARATION), results)
    for outer in OUTER:
        run_case(f"p1_{outer}", lambda: period1(outer), results)
    for outer in OUTER:
        for inner in INNER:
            run_case(f"nested_{outer}_{inner}", lambda: hc.build_nested(outer_blasts=outer, inner_blasts=inner), results)
    for outer, inner in ORDER_CELLS:
        run_case(f"outer_first_{outer}_{inner}", lambda: outer_first(outer, inner), results)

    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"grid_{label}.json"
    path.write_text(json.dumps({"tree": str(tree), "label": label, "results": results}, indent=1))
    print("wrote", path)


# --------------------------------------------------------------------------- #
# Pass two: the figures
# --------------------------------------------------------------------------- #
GREEN, PALE, RED, GREY = "#7fc97f", "#fdf0a6", "#f08080", "#9e9e9e"


def cell_verdict(results: dict, outer: int, inner: int) -> tuple[str, str]:
    """A nested cell's colour and text: each tangle against itself alone, reliability, order."""
    cell = results[f"nested_{outer}_{inner}"]
    if not cell["ok"]:
        return GREY, "crash\n" + cell["error"].split(":")[0]
    outer_tangle, inner_tangle = cell["tangles"]
    lines, matches, reliable = [], [], True
    for tag, tangle, alone in (("P1", outer_tangle, f"p1_{outer}"), ("P3", inner_tangle, f"p3_{inner}")):
        reference = results[alone]
        if not reference["ok"]:
            lines.append(f"{tag}: alone crashes")
            matches.append(False)
            continue
        same = tangle["words"] == reference["tangles"][0]["words"]
        matches.append(same)
        reliable &= tangle["reliable"]
        lines.append(f"{tag} {'=' if same else '≠'} alone" + ("" if tangle["reliable"] else " (unrel.)"))
    order = results.get(f"outer_first_{outer}_{inner}")
    if order is not None:
        same_order = order["ok"] and [t["words"] for t in order["tangles"]] == [
            t["words"] for t in cell["tangles"]
        ]
        lines.append("order: " + ("same" if same_order else "differs" if order["ok"] else "crash"))
        matches.append(same_order)
    if cell["skipped_bridges"]:
        lines.append(f"{cell['skipped_bridges']} skipped bridge(s)")
    colour = RED if not all(matches) else GREEN if reliable else PALE
    return colour, "\n".join(lines)


def draw_grid(ax, data: dict) -> None:
    """One tree's nested grid: rows = outer blasts, columns = inner blasts."""
    from matplotlib.patches import Rectangle

    results = data["results"]
    for row, outer in enumerate(OUTER):
        for col, inner in enumerate(INNER):
            colour, text = cell_verdict(results, outer, inner)
            ax.add_patch(Rectangle((col, row), 1, 1, facecolor=colour, edgecolor="white", linewidth=2))
            ax.text(col + 0.5, row + 0.5, text, ha="center", va="center", fontsize=8.5)
    ax.set_xlim(0, len(INNER))
    ax.set_ylim(len(OUTER), 0)
    ax.set_xticks([c + 0.5 for c in range(len(INNER))], [str(i) for i in INNER])
    ax.set_yticks([r + 0.5 for r in range(len(OUTER))], [str(o) for o in OUTER])
    ax.set_xlabel("inner (period-3) blasts")
    ax.set_ylabel("outer (period-1) blasts")
    ax.set_title(f"{data['label']}: nested words vs each tangle alone", fontsize=12)
    ax.set_aspect("equal")


def draw_sweep(axes, datasets: list[dict]) -> None:
    """Period 3 alone by blast count: hole and class counts, and word stability."""
    for data, marker in zip(datasets, ("s", "o")):
        results = data["results"]
        counts = {key: [] for key in ("holes", "classes", "active")}
        for blasts in P3_BLASTS:
            case = results[f"p3_{blasts}"]
            for key in counts:
                counts[key].append(case["tangles"][0][key] if case["ok"] else None)
        for key, colour in zip(counts, ("#1f77b4", "#ff7f0e", "#2ca02c")):
            xs = [b for b, v in zip(P3_BLASTS, counts[key]) if v is not None]
            ys = [v for v in counts[key] if v is not None]
            axes[0].plot(xs, ys, marker=marker, color=colour, linestyle="-" if data["label"] == "branch" else ":",
                         label=f"{key} ({data['label']})")
        for blasts in P3_BLASTS:
            if not results[f"p3_{blasts}"]["ok"]:
                axes[0].annotate(f"{data['label']}\nraises", (blasts, 0), ha="center", va="bottom", color="dimgray", fontsize=8)
    axes[0].set_ylim(bottom=-1)
    axes[0].set_xlabel("blasts of the period-3 zone")
    axes[0].set_ylabel("count")
    axes[0].set_title("period 3 alone: holes and classes", fontsize=12)
    axes[0].legend(fontsize=8, ncol=2)

    branch = next(d for d in datasets if d["label"] == "branch")["results"]
    ax = axes[1]
    ax.set_axis_off()
    first = branch["p3_0"]["tangles"][0]["words"]
    rows = []
    for blasts in P3_BLASTS:
        tangle = branch[f"p3_{blasts}"]["tangles"][0]
        same = "same as 0 blasts" if tangle["words"] == first else "DIFFERENT words"
        letters = ", ".join(f"{k}->{v or '1'}" for k, v in sorted(tangle["letters"].items(), key=lambda kv: (len(kv[0]), kv[0])))
        if len(letters) > 70:
            letters = letters[:67] + "..."
        rows.append([str(blasts), same, "yes" if tangle["reliable"] else "NO", letters])
    table = ax.table(
        cellText=rows, colLabels=["blasts", "spelled words", "reliable", "active words (letters)"],
        loc="center", cellLoc="left", colWidths=[0.07, 0.17, 0.08, 0.68],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.5)
    for (row, _col), cell in table.get_celld().items():
        if row and rows[row - 1][1] != "same as 0 blasts":
            cell.set_facecolor(PALE if rows[row - 1][2] == "NO" else "white")
    ax.set_title("branch: period 3 alone, words by blast count", fontsize=12)


def plot() -> None:
    """Pass two: both figures from the two JSONs."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    datasets = [json.loads((OUT / f"grid_{label}.json").read_text()) for label in ("base", "branch")]

    fig, axes = plt.subplots(1, 2, figsize=(17, 9.5))
    for ax, data in zip(axes, datasets):
        draw_grid(ax, data)
    fig.legend(
        handles=[
            Patch(color=GREEN, label="both tangles read exactly as alone, reliable"),
            Patch(color=PALE, label="both read as alone, but not reliable (expected: outer 0, inner >= 6)"),
            Patch(color=RED, label="a tangle reads differently from alone, or the blast order matters"),
            Patch(color=GREY, label="raises"),
        ],
        loc="lower center", ncol=2, fontsize=10,
    )
    fig.suptitle("Nested Hénon k=2: period 1 + period 3, before (2315204) and after the change", fontsize=14)
    fig.tight_layout(rect=(0, 0.07, 1, 0.96))
    fig.savefig(OUT / "nested_grid.png", dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(20, 6.5), gridspec_kw={"width_ratios": [1, 1.9]})
    draw_sweep(axes, datasets)
    fig.tight_layout()
    fig.savefig(OUT / "p3_blast_sweep.png", dpi=150)
    plt.close(fig)
    print("wrote", OUT / "nested_grid.png", OUT / "p3_blast_sweep.png")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--tree", type=Path, default=BRANCH, help="source tree to run")
    parser.add_argument("--label", default="branch", help="JSON label (base / branch)")
    parser.add_argument("--blast-both", action="store_true", help="outer-first blasts every fixed point (base recipe)")
    parser.add_argument("--plot", action="store_true", help="draw the figures from both JSONs")
    args = parser.parse_args()
    if args.plot:
        plot()
    else:
        collect(args.tree.resolve(), args.label, args.blast_both)


if __name__ == "__main__":
    main()
