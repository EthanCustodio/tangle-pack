import pytest
import numpy as np
from tanglepack.BranchPoint import BranchPoint
from tanglepack.Point import Point

def test_point_creation():

    num_branches = 2
    p = BranchPoint(num_branches, 1.0, 2.0)

    assert p.x == 1.0
    assert p.y == 2.0
    assert p.next_iterate == None
    assert p.prev_iterate == None

    assert len(p.cdists) == num_branches
    assert len(p.forward_branches) == num_branches


def test_insert_point_forward():

    num_branches = 2
    p = BranchPoint(num_branches)

    p1 = Point()
    p2 = Point()

    p.insert_point_forward(p1, 0)
    p.insert_point_forward(p2, 1)

    assert p.forward_branches[0] == p1
    assert p.forward_branches[1] == p2

    assert p1.backward == p
    assert p2.backward == p


def test_insert_point_backward():

    num_branches = 2
    p = BranchPoint(num_branches)

    p1 = Point()
    p2 = Point()

    p.insert_point_backward(p1, 0)
    p.insert_point_backward(p2, 1)

    assert p.backward_branches[0] == p1
    assert p.backward_branches[1] == p2
    assert p.forward_branches[0] == None
    assert p.forward_branches[1] == None

    assert p1.forward == p
    assert p2.forward == p
    assert p1.backward == None
    assert p2.backward == None


def test_insert_point_forward_connected():
    """Tests inserting into an already made linked list"""

    num_branches = 2
    p = BranchPoint(num_branches)

    p1 = Point()
    p2 = Point()

    p.insert_point_forward(p2, 0)
    p.insert_point_forward(p1, 0)

    # check forward connections
    assert p.forward_branches[0] == p1
    assert p1.forward == p2
    assert p2.forward == None

    # check backward connections
    assert p2.backward == p1
    assert p1.backward == p
    assert p.backward_branches[0] == None


def test_insert_point_backward_connected():
    """Tests inserting into an already made linked list"""

    num_branches = 2
    p = BranchPoint(num_branches)

    p1 = Point()
    p2 = Point()

    p.insert_point_backward(p2, 0)
    p.insert_point_backward(p1, 0)

    # check forward connections
    assert p.backward_branches[0] == p1
    assert p1.backward == p2
    assert p2.backward == None

    # check backward connections
    assert p2.forward == p1
    assert p1.forward == p
    assert p.forward_branches[0] == None

