"""Session-level pseudoneighbor, hole and partition plot helpers.

The compute / punch / partition / describe fan-out shapes are
``tests/facade/test_session_fanouts.py``.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # headless: exercise the plot helpers without a display
import matplotlib.pyplot as plt
import pytest


def test_plot_helpers_draw_pairs_and_holes(k10_session):
    """plot_pseudoneighbors computes on demand; plot_holes draws punched holes."""
    session, fp = k10_session

    plt.figure()
    try:
        pair_handles = session.plot_pseudoneighbors()
        assert len(pair_handles) == 1

        session.trellis(fp).punch_holes()
        hole_handles = session.plot_holes()
        assert hole_handles
    finally:
        plt.close()


def test_plot_stable_partition_fans_out(k10_session):
    """The session plot helper draws one row set per partitioned trellis."""
    session, fp = k10_session
    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()

    plt.figure()
    try:
        handles = session.plot_stable_partition()
        assert len(handles) == 1
    finally:
        plt.close()
