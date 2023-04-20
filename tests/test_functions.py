import numpy as np
from mylib.functions import unit_vector, angle_between


def test_unit_vector():
    a = np.empty(3)
    a.fill(1.)
    ua = unit_vector(a)
    b = np.empty(3)
    b.fill(1 / np.sqrt(3.))
    assert np.array_equal(ua, b)


def test_angle_between():
    ab = angle_between([1, 0, 0], [-1, 0, 0])
    assert ab == 3.141592653589793
