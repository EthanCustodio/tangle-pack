"""P4: the b = -1 (orientation-reversing, det J = -1) placeholder case.

Hénon ``(x, y) -> (y - k + x**2, x)`` at b = -1 has fixed points x = y = +-sqrt(k)
and Jacobian [[2x, 1], [1, 0]] with eigenvalues x +- sqrt(x**2 + 1): always a
real saddle, product -1, so exactly one eigenvalue is negative. A single
``k_value`` cannot represent that, and ``FixedPoint.set_k_value`` must reject it
with ValueError. Checks k = 4 (fixed points (2, 2) and (-2, -2)) and k = 10.
"""

from __future__ import annotations

import logging
import traceback

import numpy as np

from tanglepack import TangleSession
from tanglepack.examples import henon_jacobian, henon_map, henon_map_inverse

logging.basicConfig(level=logging.ERROR)

for k in (4.0, 10.0):
    for sign in (+1, -1):
        guess = [sign * np.sqrt(k), sign * np.sqrt(k)]
        x = guess[0]
        eig = np.linalg.eigvals(np.array([[2 * x, 1.0], [1.0, 0.0]]))
        session = TangleSession(henon_map(k, -1), henon_map_inverse(k, -1), henon_jacobian(k, -1))
        try:
            fp = session.construct_fixed_point(guess)
            outcome = f"NO ERROR: fp={np.ravel(fp.coordinates[0])} k_value={fp.k_value}"
        except Exception as exc:  # noqa: BLE001
            frame = [f for f in traceback.extract_tb(exc.__traceback__)][-1]
            outcome = f"raises {type(exc).__name__} from {frame.name} ({frame.filename.rsplit('/', 1)[-1]}:{frame.lineno})"
        print(f"k={k:g}, b=-1, guess={np.round(guess, 4).tolist()}, eigenvalues={np.round(sorted(eig), 4).tolist()} "
              f"(det={np.prod(eig):+.3f}): {outcome}")
