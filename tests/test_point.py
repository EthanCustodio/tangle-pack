import pytest
import numpy as np
from tanglepack.Point import Point

def test_point_creation():

    p = Point(1.0, 2.0)

    assert p.x == 1.0
    assert p.y == 2.0
    assert p.forward == None
    assert p.backward == None


def test_get_point():

    p = Point(1.0, 2.0)

    test_point = np.array([1.0, 2.0])

    assert p.get_point().all() == test_point.all()


def test_insert_point_forward():

    p1 = Point()
    p2 = Point()

    p1.insert_point_forward(p2)

    # check if they are connected
    assert p1.forward == p2
    assert p2.backward == p1

    # check if their ends are not
    assert p1.backward == None
    assert p2.forward == None


def test_insert_point_backward():

    p1 = Point()
    p2 = Point()

    p1.insert_point_backward(p2)

    # check if they are connected
    assert p1.backward == p2
    assert p2.forward == p1

    # check if their ends are not
    assert p1.forward == None
    assert p2.backward == None


def test_insert_point_forward_connected():
    """Tests inserting into an already made linked list"""

    p1 = Point()
    p2 = Point()
    p3 = Point()

    p1.insert_point_forward(p3)
    p1.insert_point_forward(p2)

    # check forward connections
    assert p1.forward == p2
    assert p2.forward == p3
    assert p3.forward == None

    # check backward connections
    assert p3.backward == p2
    assert p2.backward == p1
    assert p1.backward == None


def test_insert_point_backward_connected():
    """Tests inserting into an already made linked list"""

    p1 = Point()
    p2 = Point()
    p3 = Point()

    p3.insert_point_backward(p1)
    p3.insert_point_backward(p2)

    # check forward connections
    assert p1.forward == p2
    assert p2.forward == p3
    assert p3.forward == None

    # check backward connections
    assert p3.backward == p2
    assert p2.backward == p1
    assert p1.backward == None


def test_insert_next_iterate():

    p1 = Point()
    p2 = Point()

    p1.insert_next_iterate(p2)

    # check if they are connected
    assert p1.next_iterate == p2
    assert p2.prev_iterate == p1

    # check if their ends are not
    assert p1.prev_iterate == None
    assert p2.next_iterate == None


def test_insert_prev_iterate():

    p1 = Point()
    p2 = Point()

    p1.insert_prev_iterate(p2)

    # check if they are connected
    assert p1.prev_iterate == p2
    assert p2.next_iterate == p1

    # check if their ends are not
    assert p1.next_iterate == None
    assert p2.prev_iterate == None

