"""Session-level strong-pip plot helpers on the nested session.

The classify fan-out shapes are ``tests/facade/test_session_fanouts.py``; the
trellis staleness guard (a recompute swaps the registry) is a row of
``tests/facade/test_session_caches.py``.
"""

import matplotlib

matplotlib.use("Agg")  # headless: exercise the plot helpers without a display
import matplotlib.pyplot as plt


def test_plot_helpers_cover_every_tangle(henon_p3_session):
    """One plot call draws candidates/strong pips for both inner and outer tangles."""
    session, _fp3, _fp1, _zone = henon_p3_session

    plt.figure()
    try:
        candidate_handles = session.plot_strong_pip_candidates()
        strong_pip_handles = session.plot_strong_pip()
    finally:
        plt.close()

    # One handle per fixed point (the nested session has two).
    assert len(candidate_handles) == 2
    assert len(strong_pip_handles) == 2
