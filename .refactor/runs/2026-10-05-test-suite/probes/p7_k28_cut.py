"""P7: the k=2.8 two-blast cut facts of CLAUDE.md, with ids and without.

CLAUDE.md writes ``R_1^1=[0,10] R_1^2=(10,7) R_1^3=[7,6]``,
``R_5^1=[3,8] R_5^2=(8,9) R_5^3=[9,4]``, ``L_1^1=[0,3) L_1^2=[3,8]`` and "the
hole (6,5)'s lobe has base (10,7) and chord (8,9)". This probe builds the case
several times, prints every iterated element (short name, bracket, lo/hi ids,
cdists) and the holes' bounding ids, then evaluates the id-free relations the
golden tier will pin.
"""

from __future__ import annotations

import logging
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from _builds import build_k28  # noqa: E402

logging.basicConfig(level=logging.ERROR)


def bracket(iv) -> str:
    """The closedness of an interval as '[ ]', '( )', '[ )' or '( ]'."""
    return ("[" if iv.closed_lo else "(") + " " + ("]" if iv.closed_hi else ")")


def one_build(tag: str) -> tuple[dict, object]:
    """Build, print and return ({short name: interval}, trellis)."""
    session, fp = build_k28(2)
    T = session.trellis(fp)
    P = session.iterated_partition()
    naming = session.dual_graph().naming
    rows = {}
    print(f"\n--- build {tag}: holes {[ (h.iterate, h.bounding_ids) for h in T.holes]}")
    for ref, name in naming.items():
        iv = P.element(ref)
        rows[name.short_text] = iv
        print(f"   {name.short_text:8s} {name.text:14s} {bracket(iv)}  lo={iv.lo_id} hi={iv.hi_id}  "
              f"[{iv.lo_cdist:.4g}, {iv.hi_cdist if iv.hi_cdist is None else round(iv.hi_cdist, 4)}]")
    return rows, T


def relations(rows: dict, T) -> list[tuple[str, bool]]:
    """The id-free relations (plan §D, k28 item 6)."""
    g = rows.get
    out = [
        ("R_1^1 begins at the anchor (lo cdist 0, closed)", g("R_1^1").lo_cdist == 0.0 and g("R_1^1").closed_lo),
        ("R_1^1.hi == R_1^2.lo", g("R_1^1").hi_id == g("R_1^2").lo_id),
        ("R_1^2.hi == R_1^3.lo", g("R_1^2").hi_id == g("R_1^3").lo_id),
        ("R_1^2 is open both ends (lobe base)", bracket(g("R_1^2")) == "( )"),
        ("R_1^1, R_1^3 closed both ends", bracket(g("R_1^1")) == "[ ]" and bracket(g("R_1^3")) == "[ ]"),
        ("R_5^1.hi == R_5^2.lo", g("R_5^1").hi_id == g("R_5^2").lo_id),
        ("R_5^2.hi == R_5^3.lo", g("R_5^2").hi_id == g("R_5^3").lo_id),
        ("R_5^2 is open both ends (chord)", bracket(g("R_5^2")) == "( )"),
        ("R_5^1, R_5^3 closed both ends", bracket(g("R_5^1")) == "[ ]" and bracket(g("R_5^3")) == "[ ]"),
        ("L_1^2 has the same two ends as R_5^1", (g("L_1^2").lo_id, g("L_1^2").hi_id) == (g("R_5^1").lo_id, g("R_5^1").hi_id)),
        ("L_1^1/L_1^2 boundary is R_5^1.lo", g("L_1^1").hi_id == g("L_1^2").lo_id == g("R_5^1").lo_id),
        ("L_1^1 = [ ) and L_1^2 = [ ]", bracket(g("L_1^1")) == "[ )" and bracket(g("L_1^2")) == "[ ]"),
    ]
    holes = {h.iterate: set(h.bounding_ids) for h in T.holes}
    base = {T.iterate(i, 1) for i in holes[-1]}
    out += [
        ("holes at iterates exactly {0,-1,-2,-3}", set(holes) == {0, -1, -2, -3}),
        ("R_2 ends == hole(-1) bounds (the hole whose lobe is cut)", {g("R_2").lo_id, g("R_2").hi_id} == holes[-1]),
        ("R_1^2 ends == images of hole(-1) bounds (lobe base)", {g("R_1^2").lo_id, g("R_1^2").hi_id} == base),
        ("R_5^2 ends == hole(0) bounds (the chord is the reference hole's stretch)", {g("R_5^2").lo_id, g("R_5^2").hi_id} == holes[0]),
        ("L_2 ends == hole(0) bounds", {g("L_2").lo_id, g("L_2").hi_id} == holes[0]),
        ("L_1^1.hi in hole(-2) bounds (anchor-lobe fold abuts that hole)", g("L_1^1").hi_id in holes[-2]),
    ]
    return out


if __name__ == "__main__":
    for tag in ("1", "2", "3"):
        rows, T = one_build(tag)
        for name, ok in relations(rows, T):
            print(f"   {'OK ' if ok else 'NO '} {name}")
