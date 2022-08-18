
__all__ = ["measure"]


def measure(item0, item1=None,
            item2=None, item3=None,
            improper_angle: bool = False):
    """
    Measure and return the distance, angle, torsion angle or improper
    angle for the passed items. The items can be points in space,
    atoms, residues or any molecule view.

    If one item is passed, then it should be Bond, Angle, Dihedral
    or Improper. In this case, the ``measure()`` function of that
    item will be returned.

    If two items are passed, then the distance between them
    is returned.

    If three items are passed, then the angle between then
    is returned.

    If four items are passed, then the torsion angle between
    them is returned if ``improper_angle`` is ``False`` (the default).
    If ``improper_angle`` is ``True``, then the improper angle
    is returned.

    Args:
        item0: Any sire object that be converted to coordinates. This
               is either sire.maths.Vector, or if this is a molecule view,
               then it is the result of calling view.coordinates()
        item1: Any sire object that be converted to coordinates. This
               is either sire.maths.Vector, or if this is a molecule view,
               then it is the result of calling view.coordinates()
        item2: Any sire object that be converted to coordinates. This
               is either sire.maths.Vector, or if this is a molecule view,
               then it is the result of calling view.coordinates()
        item3: Any sire object that be converted to coordinates. This
               is either sire.maths.Vector, or if this is a molecule view,
               then it is the result of calling view.coordinates()
        improper_angle: Whether or not the improper angle should
               be returned (default False, as the torsion angle
               is calculated by default between four items)

      Returns:
        measurement : Either a distance or an angle depending on the
                      number of items passed.
    """
    if item0 is None:
        # They are measuring nothing...
        return 0

    if item1 is None:
        # this must be an object with a `.measure()` function
        try:
            return item0.measure()
        except AttributeError:
            pass

        raise AttributeError(
            "You can only call `measure` with a single item if that item "
            "is a Bond, Angle, Dihedral or Improper. Asking to measure "
            f"a single {item0} is not supported."
        )


    def _to_coords(item):
        if item is None:
            return item
        else:
            try:
                return item.coordinates()
            except Exception:
                return item

    item0 = _to_coords(item0)
    item1 = _to_coords(item1)
    item2 = _to_coords(item2)
    item3 = _to_coords(item3)

    if item3 is None:
        if item2 is None:
            from .maths import Vector
            from .units import angstrom
            return Vector.distance(item0, item1) * angstrom
        else:
            from .maths import Triangle
            return Triangle(item0, item1, item2).angle()

    else:
        from .maths import Torsion

        t = Torsion(item0, item1, item2, item3)

        if improper_angle:
            return t.improper_angle()
        else:
            return t.angle()
