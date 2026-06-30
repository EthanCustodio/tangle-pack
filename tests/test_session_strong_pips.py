"""Session-level strong-pip facade and trellis staleness guard.

These cover the one-call API that fans classification/plotting across every
fixed point of a nested session, and the guarantee that a Trellis handed back by
:meth:`TangleSession.trellis` is never bound to a registry the workbench has
since replaced (e.g. by a resonance-zone recompute).
"""

import matplotlib

matplotlib.use("Agg")  # headless: exercise the plot helpers without a display
import matplotlib.pyplot as plt


def test_classify_strong_pips_all_fixed_points(henon_p3_session):
    """No-arg classify returns a per-fixed-point dict, each tangle classified."""
    session, fp3, fp1, _zone = henon_p3_session

    result = session.classify_strong_pips()

    assert set(result) == {fp3, fp1}
    for fp in (fp3, fp1):
        assert result[fp], f"{fp} should have strong-pip candidates"
        # Default selection also picks the one strong pip per tangle.
        assert session.strong_pip(fp) in result[fp]


def test_classify_single_fixed_point_returns_list(henon_p3_session):
    """Passing one fixed point returns just its candidate list (not a dict)."""
    session, _fp3, fp1, _zone = henon_p3_session

    candidates = session.classify_strong_pips(fp1)

    assert isinstance(candidates, list)
    assert candidates == session.strong_pip_candidates(fp1)


def test_trellis_auto_rebuilds_after_registry_change(henon_p3_session):
    """A cached Trellis bound to a replaced registry is transparently rebuilt."""
    session, fp3, _fp1, _zone = henon_p3_session

    stale = session.trellis(fp3)
    # Recompute swaps in a brand-new registry, stranding the cached snapshot.
    session.compute_intersections([fp3, _fp1])
    assert stale.registry is not session.workbench.intersection_registry

    fresh = session.trellis(fp3)
    assert fresh is not stale
    assert fresh.registry is session.workbench.intersection_registry


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
