"""
The blasted k=2.8 session, as a fixture module.

Mirrors the ``k28_partitioned`` fixture the concurrent
``fix-backward-hole-propagation`` work adds to ``tests/conftest.py`` (same
name, same recipe: the k=2.8 case of ``scripts/henon_bridge_classes.py``).
It lives in its own module so the two branches do not both edit conftest;
after the merge delete this file and the conftest fixture takes over with no
test edits (tests import it with ``from k28_fixture import k28_partitioned``).

Session-scoped (about a minute to build); the tests that use it must not
mutate the session. Registry ids are not reproducible between builds, so
locate crossings by canonical-distance order, never by id.
"""

from __future__ import annotations

import numpy as np
import pytest

from tanglepack import TangleSession
from tanglepack.examples import (
    henon_jacobian as _henon_jacobian_factory,
    henon_map as _henon_map_factory,
    henon_map_inverse as _henon_map_inverse_factory,
)

HENON_K28 = (2.8, 1)
_k28_map = _henon_map_factory(*HENON_K28)
_k28_map_inverse = _henon_map_inverse_factory(*HENON_K28)
_k28_jacobian = _henon_jacobian_factory(*HENON_K28)


@pytest.fixture(scope="session")
def k28_partitioned():
    """``(session, fp)``: ten unstable steps at ``area_cutoff = 1e-7``, the stable
    manifold trimmed at the pip f(q0) (a resonance zone), ONE single-iteration
    blast of that zone, then classify / pseudoneighbors / holes / partition."""
    session = TangleSession(_k28_map, _k28_map_inverse, _k28_jacobian)
    session.workbench._man_machine.area_cutoff = 1e-7
    fp = session.construct_fixed_point([4, -4])
    session.orient_eigenvectors(
        fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}
    )
    session.initialize_both_manifolds(fp)
    session.grow_n_times(fp, "unstable", num_iterations=10)
    session.grow_until_turnaround(fp, "stable")
    session.compute_intersections([fp])
    session.trim_stable_manifolds(fp)
    session.create_bridges(fp)
    session.infer_iterate_table()

    session.classify_strong_pips()
    trellis = session.trellis(fp)
    pip = trellis.iterate(trellis.strong_pip, 1)
    assert pip is not None, "the default strong pip must have a registered image"

    zone = session.resonance_zone(pip)
    session.classify_strong_pips()
    session.set_strong_pip(fp, pip)
    session.blast_zone(zone, num_iterations=1, fixed_point=[fp], min_separation=1e-5)
    session.classify_strong_pips()
    session.set_strong_pip(fp, pip)

    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()
    assert session.trellis(fp).stable_partitions
    return session, fp
