from mylib.potentials import morse, lennard_jones


def test_morse():
    assert morse(1, 2, 3, 1) == 0


def test_lj():
    assert lennard_jones(1, 2, 1) == 0
