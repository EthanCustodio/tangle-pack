"""The hydrogen atom in parallel electric and magnetic fields: its Poincare map.

A Python port of the MATLAB code in ``../hydrogen`` (``hdxds.m``, ``hpotential.m``,
``x_full.m``, ``map/hmap.m``, ``map/hmapinv.m``, ``vpvidentify.m``,
``coord_transforms/qp2vpv.m``, ``hfindx.m``).

In semiparabolic coordinates ``(u, v)`` with regularized time ``s`` (``dt/ds =
u^2 + v^2``) and electric-field scaling (``F = 1``), the motion is that of the
separable Hamiltonian

    K = (pu^2 + pv^2) / 2 + V(u, v) = 2,
    V = -E (u^2 + v^2) + B^2 (u^4 v^2 + u^2 v^4) / 8 + F (u^4 - v^4) / 2.

The surface of section is ``u = 0`` (``pu`` fixed by ``K = 2``) in the
coordinates ``(v, pv)``; the map sends a point to the next crossing of ``u = 0``
in either direction (the potential is even in ``u``, so a crossing with
``pu < 0`` is the mirror of one with ``pu > 0``), exactly as ``hmap.m``.

Two symmetries are used:

* ``R(v, pv) = (v, -pv)`` (time reversal composed with ``u -> -u``) is a
  reversor, ``M^-1 = R M R``, so the inverse map is the forward integrator.
* ``-I`` commutes with ``M`` (the potential is even in ``v``). The two period-1
  saddles ``(+-v*, 0)`` are swapped by it; identifying ``z ~ -z`` with the
  area-preserving ``z -> z^2 / (sqrt(2) |z|)`` (``vpvidentify.m``) turns the
  heteroclinic pair into ONE saddle with a homoclinic tangle (``hmap_id.m``).

The integrator is a batch over points (axis 1 of a ``(2, N)`` array, the
tanglepack convention): Yoshida's 6th-order composition of Stormer-Verlet with
a step size chosen per point and per step from the local stiffness, and the
crossing of ``u = 0`` located by Newton on the step length.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

# Yoshida (1990) 6th order, solution A: the seven leapfrog weights.
_W1, _W2, _W3 = -1.17767998417887, 0.235573213359357, 0.784513610477560
_W0 = 1.0 - 2.0 * (_W1 + _W2 + _W3)
_WEIGHTS = (_W3, _W2, _W1, _W0, _W1, _W2, _W3)


class HydrogenSystem:
    """
    The hydrogen Poincare map for one set of scaled parameters.

    Args:
        b: The scaled magnetic field ``B`` (``HB`` in the MATLAB code).
        e: The scaled energy ``E`` (``HE``).
        f: The electric field (``HF``); 1 is electric-field scaling.
        h0: The base step in ``s``; the step used is ``h0 / stiffness``.
        s_max: A trajectory not back on the section by this ``s`` is an error.
    """

    def __init__(self, b: float = 3.5, e: float = 1.0, f: float = 1.0,
                 h0: float = 0.01, s_max: float = 200.0):
        self.b, self.e, self.f = float(b), float(e), float(f)
        self.h0 = float(h0)
        self.s_max = float(s_max)

    # ── the flow ─────────────────────────────────────────────────────────────

    def potential(self, u, v):
        """``hpotential.m``."""
        b, e, f = self.b, self.e, self.f
        return (-e * (u**2 + v**2) + b**2 * (u**4 * v**2 + u**2 * v**4) / 8.0
                + f * (u**4 - v**4) / 2.0)

    def force(self, u, v):
        """``-grad V`` (the ``pu``, ``pv`` rows of ``hdxds.m``)."""
        b, e, f = self.b, self.e, self.f
        u2, v2 = u * u, v * v
        fu = 2.0 * e * u - b**2 * (2.0 * u2 * u * v2 + u * v2 * v2) / 4.0 - 2.0 * f * u2 * u
        fv = 2.0 * e * v - b**2 * (2.0 * v2 * v * u2 + v * u2 * u2) / 4.0 + 2.0 * f * v2 * v
        return fu, fv

    def _step(self, state, h):
        """One Yoshida-6 step of length ``h`` (scalar or per point)."""
        u, v, pu, pv = state
        for w in _WEIGHTS:
            hw = w * h
            fu, fv = self.force(u, v)
            pu = pu + 0.5 * hw * fu
            pv = pv + 0.5 * hw * fv
            u = u + hw * pu
            v = v + hw * pv
            fu, fv = self.force(u, v)
            pu = pu + 0.5 * hw * fu
            pv = pv + 0.5 * hw * fv
        return u, v, pu, pv

    def _step_size(self, u, v):
        """``h0`` over the largest local frequency (cyclotron ~ B r^2 / 2, Stark ~ |v|)."""
        r2 = u * u + v * v
        return self.h0 / (1.0 + 0.5 * self.b * r2 + np.sqrt(6.0 * self.f) * np.abs(v))

    def lift(self, v, pv):
        """``x_full.m``: the full state on ``u = 0`` with ``pu >= 0``."""
        pu2 = 4.0 - 2.0 * self.potential(0.0, v) - pv**2
        if np.any(pu2 < 0):
            raise ValueError("a point lies outside the energetically allowed region")
        return np.zeros_like(v), v, np.sqrt(pu2), pv

    def flow_to_section(self, v0, pv0):
        """
        Integrate every point from ``u = 0`` (``pu > 0``) to its next ``u = 0``.

        Returns:
            ``(v, pv, s)``: the crossing and the regularized time to reach it.
        """
        u, v, pu, pv = self.lift(np.asarray(v0, float), np.asarray(pv0, float))
        n = u.size
        out_v, out_pv, out_s = np.empty(n), np.empty(n), np.empty(n)
        idx = np.arange(n)
        s = np.zeros(n)
        first = True
        while idx.size:
            h = self._step_size(u, v)
            nu, nv, npu, npv = self._step((u, v, pu, pv), h)
            # The start is ON the section (u = 0, pu > 0): the first step is never a crossing.
            crossed = np.zeros(idx.size, bool) if first else (np.sign(nu) != np.sign(u)) & (u != 0)
            crossed |= nu == 0
            first = False
            if crossed.any():
                c = crossed
                # Newton on the step length from the last state before the crossing.
                st = (u[c], v[c], pu[c], pv[c])
                d = -st[0] / st[2]  # linear guess
                for _ in range(6):
                    cu, cv, cpu, cpv = self._step(st, d)
                    d = d - cu / cpu
                cu, cv, cpu, cpv = self._step(st, d)
                k = idx[c]
                out_v[k], out_pv[k], out_s[k] = cv, cpv, s[c] + d
            keep = ~crossed
            idx = idx[keep]
            u, v, pu, pv = nu[keep], nv[keep], npu[keep], npv[keep]
            s = s[keep] + h[keep]
            if idx.size and s.max() > self.s_max:
                raise RuntimeError("a trajectory did not return to u = 0 by s_max")
        return out_v, out_pv, out_s

    # ── the maps on the (v, pv) plane: tanglepack's (2,) / (2, N) convention ──

    def map(self, point: NDArray) -> NDArray:
        """``hmap.m``: one iterate on the ``(v, pv)`` plane."""
        point = np.asarray(point, dtype=float)
        flat = point.reshape(2, -1)
        v, pv, _ = self.flow_to_section(flat[0], flat[1])
        return np.stack([v, pv]).reshape(point.shape)

    def map_inverse(self, point: NDArray) -> NDArray:
        """``hmapinv.m``, as ``R M R`` with the reversor ``R(v, pv) = (v, -pv)``."""
        point = np.asarray(point, dtype=float)
        flat = point.reshape(2, -1)
        v, pv, _ = self.flow_to_section(flat[0], -flat[1])
        return np.stack([v, -pv]).reshape(point.shape)

    # ── the identified plane (homoclinic form) ───────────────────────────────

    @staticmethod
    def identify(point: NDArray) -> NDArray:
        """``vpvidentify.m``: ``z -> z^2 / (sqrt(2)|z|)``, area-preserving, ``z ~ -z``."""
        point = np.asarray(point, dtype=float)
        v, pv = point[0], point[1]
        r = np.sqrt(2.0 * (v**2 + pv**2))
        with np.errstate(invalid="ignore", divide="ignore"):
            q = np.where(r != 0, (v**2 - pv**2) / r, 0.0)
            p = np.where(r != 0, 2.0 * v * pv / r, 0.0)
        return np.stack([q, p])

    @staticmethod
    def unidentify(point: NDArray) -> NDArray:
        """``qp2vpv.m``: the preimage of ``identify`` with ``v >= 0`` (either one works)."""
        point = np.asarray(point, dtype=float)
        q, p = point[0], point[1]
        # z = sqrt(2 |w| w); on the cut (p = 0, q < 0) take pv > 0.
        R = np.sqrt(q**2 + p**2)
        v = np.sqrt(np.maximum(R**2 + R * q, 0.0))
        pv = np.where(p >= 0, 1.0, -1.0) * np.sqrt(np.maximum(R**2 - R * q, 0.0))
        return np.stack([v, pv])

    def map_identified(self, point: NDArray) -> NDArray:
        """``hmap_id.m``: the map on the identified plane."""
        point = np.asarray(point, dtype=float)
        flat = point.reshape(2, -1)
        return self.identify(self.map(self.unidentify(flat))).reshape(point.shape)

    def map_identified_inverse(self, point: NDArray) -> NDArray:
        """``hmapinv_id.m``: the inverse map on the identified plane."""
        point = np.asarray(point, dtype=float)
        flat = point.reshape(2, -1)
        return self.identify(self.map_inverse(self.unidentify(flat))).reshape(point.shape)

    @staticmethod
    def jacobian_of(fn, eps: float = 1e-5):
        """A central-difference Jacobian of ``fn`` (error ~1e-9 at ``eps = 1e-5``)."""

        def jacobian(point: NDArray) -> NDArray:
            point = np.asarray(point, dtype=float)
            shifts = np.array([[eps, -eps, 0.0, 0.0], [0.0, 0.0, eps, -eps]])
            images = fn(point[:, None] + shifts)
            return np.column_stack([(images[:, 0] - images[:, 1]) / (2 * eps),
                                    (images[:, 2] - images[:, 3]) / (2 * eps)])

        return jacobian

    def jacobian(self, point: NDArray) -> NDArray:
        """The Jacobian of :meth:`map` at one point."""
        return self.jacobian_of(self.map)(point)

    def jacobian_identified(self, point: NDArray) -> NDArray:
        """The Jacobian of :meth:`map_identified` at one point."""
        return self.jacobian_of(self.map_identified)(point)

    # ── the period-1 saddle (``hfindx.m``) ───────────────────────────────────

    def fixed_point_guess(self) -> float:
        """``hfindx.m``'s analytic starting value for ``v*``."""
        if self.b >= 1:
            return np.sqrt(2) * (1 - (3 / 8) * (2 + self.e) / (0.25 * self.b**2 + 1))
        return np.sqrt(-self.e + 0.25 * self.b**2)

    def find_fixed_point(self, v_start: float | None = None) -> NDArray:
        """The saddle ``(v*, 0)``: the root of ``M(v, 0)_v - v`` (``to_solve.m``)."""
        from scipy.optimize import newton

        v_start = self.fixed_point_guess() if v_start is None else v_start
        g = lambda x: self.map(np.array([x, 0.0]))[0] - x
        v_star = newton(g, v_start, tol=1e-14, maxiter=50)
        return np.array([v_star, 0.0])
