"""
The recurrent core of the k=2.5 transition matrices, blast by blast.

Reads the entropy sweep's data (``henon_k25_entropy.py``,
``figures/henon_k25_entropy/data/<case>_a<cutoff>.json``) at ONE cutoff and,
for each case and blast, splits the refined transition matrix ``M_1`` into the
strongly connected components of its graph. Only a component that carries a
cycle (a RECURRENT component) contributes nonzero eigenvalues: the matrix is
block triangular in the component order, so

    char(M_1)(x) = x^Z * prod over recurrent components C of char(M_1[C])(x),

``Z`` the number of transient symbols. Each recurrent block's characteristic
polynomial is factored over the integers (sympy), so the factor carrying the
spectral radius rho, and whether it changes from blast to blast, is read off
exactly.

Figures, in ``figures/henon_k25_entropy/components/``:

* ``components`` -- per case, against blasts: active classes, refined symbols,
  symbols in recurrent components, symbols in the leading component (the one
  holding rho), and the number of recurrent components;
* ``polynomials`` -- per case, the factored characteristic polynomial at each
  blast where it changes (runs of unchanged blasts merged), the factor holding
  rho marked; also written as ``polynomials.txt``;
* ``spectra_by_blast`` -- per case, the nonzero eigenvalues (roots of the
  factors, with multiplicity) at each such run, in the complex plane.
* ``convergence`` -- where h settles: p1 and nested zoomed to their plateau
  (each change named by the factor holding rho), the gap h_nested - h_p1, p3
  zoomed, and the size of every change of h on a log axis.

Run with ``env/bin/python scripts/henon_k25_components.py`` (``PYTHONPATH=src``;
``--cutoff 1e-9``, ``--cases p1 p3 nested``). No computation of the tangle: it
only reads the JSON.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import networkx
import numpy as np
import sympy

sys.path.insert(0, str(Path(__file__).resolve().parent))

import henon_k25_entropy as entropy  # noqa: E402
from symbolic_figures import configure_logging, save_panel  # noqa: E402

logger = logging.getLogger(__name__)

OUT_DIR = entropy.OUT_DIR / "components"
X = sympy.Symbol("x")
#: An irreducible factor of higher degree is shown by a name, ``P_d`` (``P_d'``, ... for a
#: second factor of the same degree), and written in full in polynomials.txt.
MAX_SHOWN_DEGREE = 8
#: The name of each abbreviated factor, assigned in order of first appearance.
_NAMES: dict = {}


def factor_name(factor) -> tuple[int, int]:
    """``(degree, index)`` naming an abbreviated factor: the index-th distinct one of its degree."""
    if factor not in _NAMES:
        degree = sympy.degree(factor, X)
        _NAMES[factor] = (degree, sum(1 for d, _ in _NAMES.values() if d == degree))
    return _NAMES[factor]


@dataclass
class Step:
    """One blast's recurrent core."""

    blasts: int
    active: int
    symbols: int
    recurrent_symbols: int
    leading_symbols: int
    recurrent_components: int
    zeros: int  # Z: transient symbols plus zero roots inside recurrent blocks
    factors: tuple  # ((sympy factor, multiplicity), ...) of the nonzero part, canonical order
    rho: float
    entropy: float


def factor_block(block: np.ndarray) -> tuple[int, list[tuple]]:
    """``(power of x, [(irreducible factor, multiplicity), ...])`` of a block's characteristic polynomial."""
    poly = sympy.Matrix(block.tolist()).charpoly(X).as_expr()
    _, factors = sympy.factor_list(poly, X)
    zero_power = sum(m for f, m in factors if f == X)
    return zero_power, [(f, m) for f, m in factors if f != X]


def largest_root(factor) -> float:
    """The largest root modulus of one factor."""
    return max(abs(complex(r)) for r in sympy.Poly(factor, X).nroots(n=15))


def analyse(row: dict, cache: dict) -> Step:
    """Split one row's refined matrix into components and factor the recurrent ones."""
    matrix = np.array(row["matrix"], dtype=np.int64).reshape(len(row["symbols"]), -1)
    n = len(matrix)
    graph = networkx.DiGraph()
    graph.add_nodes_from(range(n))
    graph.add_edges_from(zip(*np.nonzero(matrix)))
    zeros, factors, sizes = 0, {}, []
    leading = (0.0, 0)
    for component in networkx.strongly_connected_components(graph):
        index = sorted(component)
        block = matrix[np.ix_(index, index)]
        if len(index) == 1 and block[0, 0] == 0:
            zeros += 1
            continue
        key = block.tobytes() + bytes(block.shape)
        if key not in cache:
            cache[key] = factor_block(block)
        zero_power, block_factors = cache[key]
        zeros += zero_power
        sizes.append(len(index))
        block_rho = max((largest_root(f) for f, _ in block_factors), default=0.0)
        if block_rho > leading[0]:
            leading = (block_rho, len(index))
        for f, m in block_factors:
            factors[f] = factors.get(f, 0) + m
    ordered = tuple(sorted(factors.items(), key=lambda item: (-largest_root(item[0]), sympy.degree(item[0], X),
                                                               str(item[0]))))
    rho = leading[0]
    return Step(
        blasts=row["blasts"], active=row["active"], symbols=n, recurrent_symbols=sum(sizes),
        leading_symbols=leading[1], recurrent_components=len(sizes), zeros=zeros, factors=ordered,
        rho=rho, entropy=entropy.entropy_of(rho),
    )


def runs(steps: list[Step]) -> list[tuple[int, int, Step]]:
    """``(first blast, last blast, step)`` per run of blasts with the same factored polynomial."""
    out: list[tuple[int, int, Step]] = []
    for step in steps:
        if out and out[-1][2].factors == step.factors and out[-1][2].zeros == step.zeros:
            out[-1] = (out[-1][0], step.blasts, out[-1][2])
        else:
            out.append((step.blasts, step.blasts, step))
    return out


def _factor_text(factor, multiplicity: int, rho: float, *, mathtext: bool) -> str:
    """One factor, abbreviated past :data:`MAX_SHOWN_DEGREE`, starred when it holds rho."""
    degree = sympy.degree(factor, X)
    holds_rho = abs(largest_root(factor) - rho) < 1e-9 * max(rho, 1.0)
    if degree <= MAX_SHOWN_DEGREE:
        body = sympy.latex(factor) if mathtext else str(factor).replace("**", "^")
    else:
        degree, index = factor_name(factor)
        primes = "'" * index
        body = (f"P_{{{degree}}}{primes}(x)" if mathtext else f"P{degree}{primes}(x)")
    text = f"({body})" + (f"^{{{multiplicity}}}" if mathtext and multiplicity > 1 else
                          f"^{multiplicity}" if multiplicity > 1 else "")
    if holds_rho:
        text = (r"\mathbf{[}" + text + r"\mathbf{]}") if mathtext else f"[{text}]"
    return text


def polynomial_text(step: Step, *, mathtext: bool) -> str:
    """The factored polynomial ``x^Z * ...``, the factor holding rho in brackets."""
    parts = []
    if step.zeros:
        parts.append(f"x^{{{step.zeros}}}" if mathtext else f"x^{step.zeros}")
    parts += [_factor_text(f, m, step.rho, mathtext=mathtext) for f, m in step.factors]
    joined = r" \cdot ".join(parts) if mathtext else " * ".join(parts)
    return f"${joined}$" if mathtext else joined


def big_factor_notes(steps: list[Step]) -> list[str]:
    """Full text of every abbreviated factor (for the .txt), once each."""
    seen, notes = set(), []
    for step in steps:
        for f, _ in step.factors:
            if sympy.degree(f, X) > MAX_SHOWN_DEGREE and f not in seen:
                seen.add(f)
                degree, index = factor_name(f)
                notes.append(f"  P{degree}{chr(39) * index}(x), largest root {largest_root(f):.10f}: "
                             f"{str(sympy.expand(f)).replace('**', '^')}")
    return notes


def roots_of(step: Step) -> np.ndarray:
    """The nonzero eigenvalues, with multiplicity, from the factors."""
    roots = []
    for f, m in step.factors:
        roots += [complex(r) for r in sympy.Poly(f, X).nroots(n=15)] * m
    return np.array(roots, dtype=complex)


# ── figures ─────────────────────────────────────────────────────────────────


def figure_components(results: dict, cutoff: float) -> Path:
    """Alphabet against recurrent-core sizes, per case."""
    fig, axes = plt.subplots(len(results), 1, figsize=(14, 3.6 * len(results)), squeeze=False)
    for ax, (case, steps) in zip(axes[:, 0], results.items()):
        x = [s.blasts for s in steps]
        ax.plot(x, [s.active for s in steps], color="0.6", marker=".", label="active classes")
        ax.plot(x, [s.symbols for s in steps], color="k", marker=".", label=r"refined symbols (size of $M_1$)")
        ax.plot(x, [s.recurrent_symbols for s in steps], color="tab:orange", marker=".",
                label="symbols in recurrent components")
        ax.plot(x, [s.leading_symbols for s in steps], color=entropy.CASE_COLORS[case], marker="o",
                label=r"symbols in the leading component (holds $\rho$)")
        ax.step(x, [s.recurrent_components for s in steps], where="mid", color="tab:purple", linestyle="--",
                label="recurrent components (count)")
        twin = ax.twinx()
        twin.plot(x, [s.entropy for s in steps], color="tab:green", linewidth=0.9, alpha=0.6)
        twin.set_ylabel(r"$h = \ln\rho$", color="tab:green")
        twin.tick_params(axis="y", colors="tab:green")
        ax.set_ylabel("symbols")
        ax.grid(alpha=0.3)
        ax.set_title(f"{entropy.CASE_TITLES[case]}", fontsize=10)
        # Outside the axes, past the h axis: every corner inside holds data.
        ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.09, 1.0), borderaxespad=0)
    axes[-1, 0].set_xlabel("blasts")
    fig.suptitle(f"Hénon k={entropy.K:g}: the recurrent core of $M_1$ against blasts (cutoff {cutoff:g})")
    fig.tight_layout()
    return save_panel(fig, OUT_DIR, "components")


def figure_polynomials(results: dict, cutoff: float) -> Path:
    """The factored polynomial at each change, per case, as a text table."""
    lines: list[tuple[str, str, str]] = []
    for case, steps in results.items():
        lines.append((entropy.CASE_TITLES[case], "", "header"))
        for first, last, step in runs(steps):
            span = f"{first}" if first == last else f"{first}–{last}"
            lines.append((f"blasts {span}  (h = {step.entropy:.6f})", polynomial_text(step, mathtext=True), "row"))
    fig, ax = plt.subplots(figsize=(16, 0.36 * len(lines) + 1.2))
    ax.axis("off")
    for i, (label, poly, kind) in enumerate(lines):
        y = 1 - (i + 0.5) / len(lines)
        if kind == "header":
            ax.text(0.0, y, label, fontsize=11, fontweight="bold", va="center", transform=ax.transAxes)
            continue
        ax.text(0.01, y, label, fontsize=9, va="center", family="monospace", transform=ax.transAxes)
        ax.text(0.27, y, poly, fontsize=10, va="center", transform=ax.transAxes)
    fig.suptitle(rf"Hénon k={entropy.K:g}: characteristic polynomial of $M_1$, factored over $\mathbb{{Z}}$ "
                 rf"(cutoff {cutoff:g}; the factor holding $\rho$ in bold brackets; $P_d$ = an irreducible "
                 rf"factor of degree $d$, in full in polynomials.txt)", fontsize=10)
    fig.tight_layout()
    return save_panel(fig, OUT_DIR, "polynomials")


def write_polynomials(results: dict, cutoff: float) -> Path:
    """The same table as plain text, with every abbreviated factor in full."""
    out = [f"Henon k={entropy.K:g}, cutoff {cutoff:g}: char(M_1) factored over Z; [ ] = the factor holding rho",
           ""]
    for case, steps in results.items():
        out.append(entropy.CASE_TITLES[case])
        for first, last, step in runs(steps):
            span = f"{first}" if first == last else f"{first}-{last}"
            out.append(f"  blasts {span:>7}  h={step.entropy:.10f}  size {step.symbols:3d}  "
                       f"recurrent {step.recurrent_symbols:3d}  leading {step.leading_symbols:3d}: "
                       f"{polynomial_text(step, mathtext=False)}")
        out += big_factor_notes(steps)
        out.append("")
    path = OUT_DIR / "polynomials.txt"
    path.write_text("\n".join(out))
    return path


def figure_spectra(results: dict, cutoff: float) -> Path:
    """The nonzero eigenvalues at each change of the polynomial, one row per case."""
    case_runs = {case: runs(steps) for case, steps in results.items()}
    columns = max(len(r) for r in case_runs.values())
    fig, axes = plt.subplots(len(case_runs), columns, figsize=(3.2 * columns, 3.5 * len(case_runs)),
                             squeeze=False)
    theta = np.linspace(0, 2 * np.pi, 300)
    limit = 1.15 * max(step.rho for r in case_runs.values() for *_, step in r)
    for row, (case, case_run) in enumerate(case_runs.items()):
        for col in range(columns):
            ax = axes[row, col]
            if col >= len(case_run):
                ax.axis("off")
                continue
            first, last, step = case_run[col]
            z = roots_of(step)
            ax.plot(np.cos(theta), np.sin(theta), color="0.7", linewidth=0.7, linestyle=":")
            if step.rho > 0:
                ax.plot(step.rho * np.cos(theta), step.rho * np.sin(theta), color=entropy.CASE_COLORS[case],
                        linewidth=0.8, alpha=0.6)
            if len(z):
                ax.scatter(z.real, z.imag, s=14, color=entropy.CASE_COLORS[case], alpha=0.8, zorder=3)
            ax.scatter([0], [0], s=20, marker="x", color="k", zorder=3)
            ax.set_xlim(-limit, limit)
            ax.set_ylim(-limit, limit)
            ax.set_aspect("equal")
            ax.tick_params(labelsize=6)
            ax.grid(alpha=0.25)
            span = f"{first}" if first == last else f"{first}–{last}"
            ax.set_title(f"blasts {span}: rho = {step.rho:.5f}\n{len(z)} nonzero, x at 0: {step.zeros}",
                         fontsize=8)
        axes[row, 0].set_ylabel(entropy.CASE_TITLES[case], fontsize=9)
    fig.suptitle(rf"Hénon k={entropy.K:g}: nonzero eigenvalues of $M_1$ at each change of its polynomial "
                 rf"(cutoff {cutoff:g}; circle: $\rho$; dotted: unit circle)")
    fig.tight_layout()
    return save_panel(fig, OUT_DIR, "spectra_by_blast")


def _zoom(ax, values: list[float], pad: float = 0.15, within: float = 0.1) -> None:
    """Limit the y axis to the plateau: the values within ``within`` (a fraction) of the largest."""
    top = max(values)
    near = [v for v in values if v >= (1 - within) * top]
    lo, hi = min(near), max(near)
    span = max(hi - lo, 1e-5 * max(abs(hi), 1.0))
    ax.set_ylim(lo - pad * span, hi + pad * span)


def _label_plateaus(ax, steps: list[Step], color: str, above: bool) -> None:
    """Name the factor holding rho at the start of each run of an unchanged polynomial, on the plateau."""
    lo, hi = ax.get_ylim()
    shown = None
    for first, _, step in runs(steps):
        if not step.factors or step.entropy < lo or step.entropy > hi:
            continue
        lead = _factor_text(step.factors[0][0], 1, step.rho, mathtext=True).replace(r"\mathbf{[}", "").replace(
            r"\mathbf{]}", "")
        if lead == shown:
            continue
        shown = lead
        ax.annotate(f"${lead}$\nh = {step.entropy:.6f}", (first, step.entropy),
                    xytext=(0, 9 if above else -11), textcoords="offset points", fontsize=7, color=color,
                    ha="left", va="bottom" if above else "top")


def figure_convergence(results: dict, cutoff: float) -> Path:
    """
    Where h settles: the plateaus on a zoomed axis, the nested-minus-p1 gap, and the size of every change.

    (a) p1 and nested, the y axis zoomed to the plateau, each change labelled with
    the factor holding rho; (b) h_nested - h_p1 per blast; (c) p3, zoomed;
    (d) |h_n - h_{n-1}| at every blast where h changed, on a log axis.
    """
    fig, axes = plt.subplots(2, 2, figsize=(15, 9.5))
    (a, b), (c, d) = axes
    for case, above in (("p1", False), ("nested", True)):
        if case not in results:
            continue
        steps = results[case]
        a.plot([s.blasts for s in steps], [s.entropy for s in steps], color=entropy.CASE_COLORS[case],
               marker="o", markersize=4, drawstyle="steps-post", label=entropy.CASE_TITLES[case])
    plateau = [s.entropy for case in ("p1", "nested") if case in results for s in results[case]]
    if plateau:
        _zoom(a, plateau, pad=0.35)
        for case, above in (("p1", False), ("nested", True)):
            if case in results:
                _label_plateaus(a, results[case], entropy.CASE_COLORS[case], above)
    a.set_title("(a) period 1 alone and nested, zoomed to the plateau", fontsize=10)
    a.set_ylabel(r"$h = \ln\rho(M_1)$")
    a.legend(fontsize=8, loc="lower right")

    if "p1" in results and "nested" in results:
        p1 = {s.blasts: s.entropy for s in results["p1"]}
        pairs = [(s.blasts, s.entropy - p1[s.blasts]) for s in results["nested"]
                 if s.blasts in p1 and p1[s.blasts] > 0 and s.entropy > 0]
        x, gap = zip(*pairs)
        b.plot(x, gap, color=entropy.CASE_COLORS["nested"], marker="o", markersize=4, drawstyle="steps-post")
        b.axhline(0, color="0.5", linewidth=0.8)
        b.annotate(f"last blast: {gap[-1]:+.2e}", (x[-1], gap[-1]), xytext=(-5, 8), textcoords="offset points",
                   ha="right", fontsize=8)
        top = max(max(gap), 1e-12)
        b.set_ylim(min(min(gap), 0.0) - 0.15 * top, 1.35 * top)
    b.set_title(r"(b) $h_{\mathrm{nested}} - h_{\mathrm{p1}}$ at the same blast count (both nonzero)", fontsize=10)
    b.set_ylabel(r"$\Delta h$")

    if "p3" in results:
        steps = results["p3"]
        c.plot([s.blasts for s in steps], [s.entropy for s in steps], color=entropy.CASE_COLORS["p3"], marker="o",
               markersize=4, drawstyle="steps-post", label=entropy.CASE_TITLES["p3"])
        _zoom(c, [s.entropy for s in steps if s.entropy > 0] or [0.0], pad=0.35, within=0.25)
        _label_plateaus(c, steps, entropy.CASE_COLORS["p3"], above=False)
        c.legend(fontsize=8, loc="lower right")
    c.set_title("(c) period 3 alone, zoomed to the plateau", fontsize=10)
    c.set_ylabel(r"$h = \ln\rho(M_1)$")

    offsets = {"p1": -0.15, "p3": 0.0, "nested": 0.15}
    for case, steps in results.items():
        jumps = [(s.blasts, abs(s.entropy - prev.entropy)) for prev, s in zip(steps, steps[1:])
                 if abs(s.entropy - prev.entropy) > 1e-12 and prev.entropy > 0]
        if jumps:
            x, size = zip(*jumps)
            d.scatter(np.array(x) + offsets.get(case, 0.0), size, color=entropy.CASE_COLORS[case], s=28,
                      label=entropy.CASE_TITLES[case], zorder=3)
            d.vlines(np.array(x) + offsets.get(case, 0.0), 1e-12, size, color=entropy.CASE_COLORS[case],
                     linewidth=0.8, alpha=0.6)
    quiet = [entropy.CASE_TITLES[case] for case, steps in results.items()
             if not any(abs(s.entropy - prev.entropy) > 1e-12 and prev.entropy > 0 for prev, s in zip(steps, steps[1:]))]
    if quiet:
        d.text(0.98, 0.02, "no change after the first plateau: " + ", ".join(quiet), transform=d.transAxes,
               ha="right", va="bottom", fontsize=8, color="0.3")
    d.set_yscale("log")
    d.set_ylim(bottom=1e-7)
    d.set_title(r"(d) size of every change of $h$ after it first leaves 0, $|h_n - h_{n-1}|$", fontsize=10)
    d.set_ylabel(r"$|\Delta h|$ per blast")
    d.legend(fontsize=8, loc="upper left")
    last = max(s.blasts for steps in results.values() for s in steps)
    for ax in axes.flat:
        ax.grid(alpha=0.3)
        ax.set_xlabel("blasts")
        ax.set_xlim(-0.5, last + 0.5)
    fig.suptitle(rf"Hénon k={entropy.K:g}: where the entropy settles (cutoff {cutoff:g})")
    fig.tight_layout()
    return save_panel(fig, OUT_DIR, "convergence")


def main() -> None:
    parser = argparse.ArgumentParser(description="The recurrent core of the k=2.5 transition matrices.")
    parser.add_argument("--cases", nargs="+", choices=entropy.CASES, default=list(entropy.CASES))
    parser.add_argument("--cutoff", type=float, default=1e-9)
    args = parser.parse_args()
    configure_logging()
    plt.switch_backend("Agg")

    results: dict[str, list[Step]] = {}
    for case in args.cases:
        record = json.loads(entropy.data_path(case, args.cutoff).read_text())
        cache: dict = {}
        results[case] = [analyse(row, cache) for row in record["rows"]]
        for step, row in zip(results[case], record["rows"]):
            assert abs(step.rho - row["radius"]) < 1e-6 * max(1.0, row["radius"]), (case, step.blasts)
        logger.info("%s: %d blasts, %d distinct polynomials", case, len(results[case]), len(runs(results[case])))
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(write_polynomials(results, args.cutoff).read_text())
    figure_components(results, args.cutoff)
    figure_polynomials(results, args.cutoff)
    figure_spectra(results, args.cutoff)
    figure_convergence(results, args.cutoff)


if __name__ == "__main__":
    main()
