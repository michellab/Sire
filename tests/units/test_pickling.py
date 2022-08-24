

def _get_dims(u):
    return [u.angle(), u.charge(), u.length(), u.mass(),
            u.quantity(), u.temperature(), u.time()]


def test_pickling():
    import pickle
    from sire.units import kcal_per_mol, nanometer2

    a = kcal_per_mol / nanometer2

    a_dims = _get_dims(a)

    s = pickle.dumps(a)
    b = pickle.loads(s)

    b_dims = _get_dims(b)

    assert a_dims == b_dims

    c = b / a

    assert _get_dims(c) == [0, 0, 0, 0, 0, 0, 0]


if __name__ == "__main__":
    test_pickling()

