# tanglepack_webdash/parser.py
from __future__ import annotations
import numpy as np
from sympy import symbols, sympify, lambdify, Tuple

x, y = symbols("x y")


def _split_vector_expr(txt: str):
    txt = txt.strip()
    if not txt:
        raise ValueError("Empty function text.")
    expr = sympify(f"({txt})", convert_xor=True)
    if isinstance(expr, Tuple):
        parts = tuple(expr)
    elif isinstance(expr, (tuple, list)):
        parts = tuple(expr)
    else:
        raise ValueError(
            "Expected two expressions separated by a comma, e.g. 'x + y**2, -x'."
        )
    if len(parts) != 2:
        raise ValueError(f"Expected exactly 2 expressions, got {len(parts)}.")
    return parts[0], parts[1]


def _point_only_any_axis(fxy):
    """
    Public callable F accepts:
      - F(point) with shape (2,)
      - F(arr) where exactly one axis has length 2 (coords), anywhere in the shape.
    We move that axis to the end, evaluate, then move it back to preserve shape.
    """

    def F(v):
        a = np.asarray(v)

        if a.ndim == 0:
            raise ValueError("Expected a 2D point; got scalar.")

        # Simple vector
        if a.ndim == 1:
            if a.shape[0] != 2:
                raise ValueError(f"Expected a 2D point of length 2; got {a.shape}.")
            fx, fy = fxy(a[0], a[1])
            return np.array([fx, fy], float)

        # Batched: find a single coordinate axis of length 2
        axes_len2 = [ax for ax, n in enumerate(a.shape) if n == 2]
        if not axes_len2:
            raise ValueError(f"No axis of length 2 found in shape {a.shape}.")
        # Heuristic: prefer the LAST axis of length 2
        coord_ax = axes_len2[-1]

        moved = np.moveaxis(a, coord_ax, -1)  # (..., 2)
        X, Y = moved[..., 0], moved[..., 1]
        fx, fy = fxy(X, Y)
        out = np.stack([np.asarray(fx), np.asarray(fy)], axis=-1)  # (..., 2)
        out = np.moveaxis(out, -1, coord_ax)  # restore original axis position
        return out

    return F


def parse_map_text(txt: str):
    fx_expr, fy_expr = _split_vector_expr(txt)
    fxy = lambdify((x, y), (fx_expr, fy_expr), modules="numpy")
    return _point_only_any_axis(fxy)
