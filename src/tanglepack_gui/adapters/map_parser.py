# src/tanglepack_gui/adapters/map_parser.py
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
            "Expected two expressions separated by a comma, e.g. 'x + y**2, 0.3*x - 0.2*y'."
        )
    if len(parts) != 2:
        raise ValueError(f"Expected exactly 2 expressions, got {len(parts)}.")
    return parts[0], parts[1]


def _vectorized_point_only(fxy):
    """
    Public callable F only accepts:
      - F(point) with shape (2,)
      - F(points) with any shape that includes ONE 2-length axis for coordinates,
        not necessarily the last (e.g. (N,2), (2,N), (2,2,8), etc.).
    Returns an array with the SAME batch axes layout; the coordinate axis stays
    where it was in the input.
    """

    def F(v):
        arr = np.asarray(v)

        # Scalar -> error
        if arr.ndim == 0:
            raise ValueError("Expected a 2D point; got scalar.")

        # Simple vector (2,)
        if arr.ndim == 1:
            if arr.shape[0] != 2:
                raise ValueError(
                    f"Expected a 2D point of length 2; got shape {arr.shape}."
                )
            fx, fy = fxy(arr[0], arr[1])
            return np.array([fx, fy], dtype=float)

        # Batched: find a coordinate axis of length 2
        axes_len_2 = [ax for ax, n in enumerate(arr.shape) if n == 2]
        if not axes_len_2:
            raise ValueError(
                f"No coordinate axis of length 2 found in shape {arr.shape}."
            )

        # Heuristic: prefer the LAST axis with length 2 (works for (2,2,8) -> use axis=1)
        coord_ax = axes_len_2[-1]

        # Move that coord axis to the end → (..., 2)
        moved = np.moveaxis(arr, coord_ax, -1)

        # Now split coords and evaluate
        x_ = moved[..., 0]
        y_ = moved[..., 1]
        fx, fy = fxy(x_, y_)

        fx = np.asarray(fx)
        fy = np.asarray(fy)
        out = np.stack([fx, fy], axis=-1)  # (..., 2)

        # Move coord axis back where it was in the input
        out = np.moveaxis(out, -1, coord_ax)
        return out

    return F


def parse_map_text(txt: str):
    fx_expr, fy_expr = _split_vector_expr(txt)
    free = fx_expr.free_symbols | fy_expr.free_symbols
    bad = free - {x, y}
    if bad:
        raise ValueError(
            f"Unknown symbol(s): {', '.join(sorted(str(s) for s in bad))}. Only 'x' and 'y' are allowed."
        )
    fxy = lambdify((x, y), (fx_expr, fy_expr), modules="numpy")
    return _vectorized_point_only(fxy)
