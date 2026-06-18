import pytest
from tanglepack import Point
from tanglepack import BranchPoint
from tanglepack import BaseManifold
import numpy as np
import matplotlib
matplotlib.use("Agg") #suppresses graphics
import matplotlib.pyplot as plt

def make_point_chain(xys):
    """Creates a forward/backward-linked list of Points from (x, y) tuples."""

    nodes = [Point(x, y) for x, y in xys]
    for i in range(len(nodes) - 1):
        nodes[i].forward = nodes[i + 1]
        nodes[i + 1].backward = nodes[i]
    return nodes


def make_branch_point_pair():
    """Creates a BranchPoint with two incoming and two outgoing branches, all linked."""

    bp = BranchPoint(num_branches=2)

    a, b = Point(0, 0), Point(1, 0)  # incoming branches
    c, d = Point(0, 1), Point(1, 1)  # outgoing branches

    # Link branches to branchpoint
    bp.backward_branches[0] = a
    bp.backward_branches[1] = d
    bp.forward_branches[0] = c
    bp.forward_branches[1] = b

    # Link nodes to branchpoint
    a.forward = bp
    d.forward = bp
    c.backward = bp
    b.backward = bp

    return a, b, bp, c, d


def test_walk_fwd_back_unstable_linear():
    """Tests walking along an unstable manifold"""

    pts = make_point_chain([(0, 0), (1, 1), (2, 2)])
    manifold = BaseManifold(root=pts[0], stability="unstable", stretch_param=1.0)

    # Forward walk
    assert manifold.walk_fwd(None, pts[0]) is pts[1]
    assert manifold.walk_fwd(pts[0], pts[1]) is pts[2]

    # Backward walk
    assert manifold.walk_back(pts[2], pts[1]) is pts[0]
    assert manifold.walk_back(pts[1], pts[0]) is None


def test_walk_fwd_back_stable_linear():
    """Tests walking along a stable manifold"""

    pts = make_point_chain([(0, 0), (1, 1), (2, 2)])
    manifold = BaseManifold(root=pts[-1], stability="stable", stretch_param=1.0)

    # Forward walk goes backward
    assert manifold.walk_fwd(None, pts[2]) is pts[1]
    assert manifold.walk_fwd(pts[2], pts[1]) is pts[0]

    # Backward walk goes forward
    assert manifold.walk_back(pts[0], pts[1]) is pts[2]
    assert manifold.walk_back(pts[1], pts[2]) is None


def test_make_point_list():
    """tests if a list of points is extracted correctly"""

    coords = np.array([(0, 0), (1, 1), (2, 2)])
    points = make_point_chain(coords)

    manifold = BaseManifold(root=points[0], stability="unstable", stretch_param=1.0)

    recovered_coords = manifold.get_point_array()
    recovered_points = manifold.get_point_array(return_nodes=True)

    np.testing.assert_array_equal(coords, recovered_coords)

    for i in range(len(recovered_points)):
        assert points[i] is recovered_points[i]


def test_make_empty_list():
    """tests if the list produced is empty if there are no nodes"""

    manifold = BaseManifold(root=None, stability="unstable", stretch_param=1.0)

    np.testing.assert_array_equal(manifold.get_point_array(), [])


def test_branch_forward_unstable_toggles_branches():
    """Tests walking through a branch point along an unstable manifold"""

    a, b, bp, c, d = make_branch_point_pair()
    manifold = BaseManifold(root=a, stability="unstable", stretch_param=1.0)

    # Entering bp from a should exit on d
    assert manifold.walk_fwd(None, a) is bp
    assert manifold.walk_fwd(a, bp) is c
    # Entering from b should exit on c
    assert manifold.walk_back(None, c) is bp
    assert manifold.walk_back(c, bp) is a


def test_branch_forward_stable_toggles_branches():
    """Tests walking through a branch point along a stable manifold"""

    a, b, bp, c, d = make_branch_point_pair()
    manifold = BaseManifold(root=b, stability="stable", stretch_param=1.0)

    # Entering from c should exit on b
    assert manifold.walk_fwd(None, b) is bp
    assert manifold.walk_fwd(b, bp) is d
    # Entering from d should exit on a
    assert manifold.walk_back(None, d) is bp
    assert manifold.walk_back(d, bp) is b


def test_branch_forward_raises_on_invalid_entry():
    """checks if an error is raised when the previous node is not connected"""

    a, b, bp, c, d = make_branch_point_pair()
    fake = Point(9, 9)
    manifold = BaseManifold(root=bp, stability="unstable", stretch_param=1.0)

    with pytest.raises(ValueError):
        manifold.walk_fwd(fake, bp)


def test_branch_backward_raises_on_invalid_entry():
    """checks if an error is raised when the next node is not connected"""

    a, b, bp, c, d = make_branch_point_pair()
    fake = Point(9, 9)
    manifold = BaseManifold(root=bp, stability="unstable", stretch_param=1.0)

    with pytest.raises(ValueError):
        manifold.walk_back(fake, bp)


def test_walk_back_from_branch_point_error_raised():
    """checks if an error is raised when not supplying previous point from a branch"""

    a, b, bp, c, d = make_branch_point_pair()

    manifold = BaseManifold(root=bp, stability="unstable", stretch_param=1.0)

    with pytest.raises(ValueError):
        manifold.walk_back(None, bp)


def test_walk_fwd_from_branch_point_error_raised():
    """checks if an error is raised when not supplying previous point from a branch"""

    a, b, bp, c, d = make_branch_point_pair()

    manifold = BaseManifold(root=bp, stability="unstable", stretch_param=1.0)

    with pytest.raises(ValueError):
        manifold.walk_fwd(None, bp)


def test_walk_fwd_with_only_branch_index():
    """walks forward from a branch point with an index"""

    a, b, bp, c, d = make_branch_point_pair()

    manifold = BaseManifold(root=bp, stability="unstable", stretch_param=1.0)

    assert manifold.walk_fwd(None, bp, branch_index=0) is c

    manifold = BaseManifold(root=bp, stability="stable", stretch_param=1.0)

    assert manifold.walk_fwd(None, bp, branch_index=1) is d


def test_walk_back_with_only_branch_index():
    """walks forward from a branch point with an index"""

    a, b, bp, c, d = make_branch_point_pair()

    manifold = BaseManifold(root=bp, stability="unstable", stretch_param=1.0)

    assert manifold.walk_back(None, bp, branch_index=0) is a

    manifold = BaseManifold(root=bp, stability="stable", stretch_param=1.0)

    assert manifold.walk_back(None, bp, branch_index=1) is b


def test_manifold_plot(monkeypatch):
    """Tests plotting. Made by chat GPT"""

    # Arrange – build a tiny manifold with known coords
    pts = np.array([[0, 0], [1, 1], [2, 0]])
    m   = BaseManifold(root=pts[0], stability="unstable", stretch_param=1.0)
    # mock get_point_array so we don't need the full linked list
    monkeypatch.setattr(m, "get_point_array", lambda branch_index=None: pts)

    # Act
    fig = plt.figure()          # create a figure you'll inspect
    ax  = fig.add_subplot(111)
    with monkeypatch.context() as mctx:
        # patch plt.* inside the method so they draw on *this* axes
        mctx.setattr(plt, "plot",  ax.plot)
        mctx.setattr(plt, "scatter", ax.scatter)
        m.plot(color="red", show_points=True)

    # Assert – check the Line2D data and scatter offsets
    lines = ax.get_lines()
    assert len(lines) == 1
    x, y = lines[0].get_data()
    np.testing.assert_array_equal(x, pts[:, 0])
    np.testing.assert_array_equal(y, pts[:, 1])

    col = ax.collections[0]       # first PathCollection (scatter)
    np.testing.assert_array_equal(col.get_offsets(), pts)

    plt.close(fig)   # cleanup


def test_empty_manifold_plot():

    m   = BaseManifold(None, stability="unstable", stretch_param=1.0)

    with pytest.raises(ValueError):
        m.plot()


