import pytest
from tanglepack.Point import Point
from tanglepack.BranchPoint import BranchPoint
from tanglepack.BaseManifold import BaseManifold


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
    bp.backward_branches[1] = b
    bp.forward_branches[0] = c
    bp.forward_branches[1] = d

    # Link nodes to branchpoint
    a.forward = bp
    b.forward = bp
    c.backward = bp
    d.backward = bp

    return a, b, bp, c, d


def test_walk_fwd_back_unstable_linear():
    pts = make_point_chain([(0, 0), (1, 1), (2, 2)])
    manifold = BaseManifold(root=pts[0], stability="unstable", stretch_param=1.0)

    # Forward walk
    assert manifold.walk_fwd(None, pts[0]) is pts[1]
    assert manifold.walk_fwd(pts[0], pts[1]) is pts[2]

    # Backward walk
    assert manifold.walk_back(pts[2], pts[1]) is pts[0]
    assert manifold.walk_back(pts[1], pts[0]) is None


def test_walk_fwd_back_stable_linear():
    pts = make_point_chain([(0, 0), (1, 1), (2, 2)])
    manifold = BaseManifold(root=pts[-1], stability="stable", stretch_param=1.0)

    # Forward walk goes backward
    assert manifold.walk_fwd(None, pts[2]) is pts[1]
    assert manifold.walk_fwd(pts[2], pts[1]) is pts[0]

    # Backward walk goes forward
    assert manifold.walk_back(pts[0], pts[1]) is pts[2]
    assert manifold.walk_back(pts[1], pts[2]) is None


def test_branch_forward_unstable_toggles_branches():
    a, b, bp, c, d = make_branch_point_pair()
    manifold = BaseManifold(root=bp, stability="unstable", stretch_param=1.0)

    # Entering bp from a should exit on d
    assert manifold.walk_fwd(a, bp) is d
    # Entering from b should exit on c
    assert manifold.walk_fwd(b, bp) is c


def test_branch_forward_stable_toggles_branches():
    a, b, bp, c, d = make_branch_point_pair()
    manifold = BaseManifold(root=bp, stability="stable", stretch_param=1.0)

    # Entering from c should exit on b
    assert manifold.walk_fwd(c, bp) is b
    # Entering from d should exit on a
    assert manifold.walk_fwd(d, bp) is a


def test_branch_forward_raises_on_invalid_entry():
    a, b, bp, c, d = make_branch_point_pair()
    fake = Point(9, 9)
    manifold = BaseManifold(root=bp, stability="unstable", stretch_param=1.0)

    with pytest.raises(ValueError):
        manifold.walk_fwd(fake, bp)


def test_branch_backward_unstable():
    a, b, bp, c, d = make_branch_point_pair()
    manifold = BaseManifold(root=bp, stability="unstable", stretch_param=1.0)

    # Coming *from* c (which is forward branch), should exit back to b
    assert manifold.walk_back(c, bp) is b
    # Coming *from* d, should exit to a
    assert manifold.walk_back(d, bp) is a


def test_branch_backward_stable():
    a, b, bp, c, d = make_branch_point_pair()
    manifold = BaseManifold(root=bp, stability="stable", stretch_param=1.0)

    # Coming *from* a (which is forward), should exit to d
    assert manifold.walk_back(a, bp) is d
    # Coming *from* b, should exit to c
    assert manifold.walk_back(b, bp) is c


def test_branch_backward_raises_on_invalid_entry():
    a, b, bp, c, d = make_branch_point_pair()
    fake = Point(99, 99)
    manifold = BaseManifold(root=bp, stability="unstable", stretch_param=1.0)

    with pytest.raises(ValueError):
        manifold.walk_back(fake, bp)

