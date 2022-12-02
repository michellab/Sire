
__all__ = ["measure"]


def measure(item0, item1=None,
            item2=None, item3=None,
            improper_angle: bool = False,
            ignore_space: bool = False,
            map=None):
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

    Note that this will take into account any periodic boundary
    conditions. The first-found space will be used to map all
    points into a minimum-image-convention box, and the measurements
    will be made from these. Set 'ignore_space=True' if you want
    to ignore the space and perform measurements without
    periodic boundary conditions.

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
        improper_angle: bool. Whether or not the improper angle should
               be returned (default False, as the torsion angle
               is calculated by default between four items)
        ignore_space: bool. Whether or not to ignore any space found
               in the passed sire objects that instead to perform
               the measurements in an infinite cartesian space.

      Returns:
        measurement : Either a distance or an angle depending on the
                      number of items passed.
    """
    if item0 is None:
        # They are measuring nothing...
        return 0

    from .base import create_map

    map = create_map(map)

    if item1 is None:
        # this must be an object with a `.measure()` function
        if hasattr(item0, "measure"):
            return item0.measure(map=map)

        # or it could be a list of objects
        try:
            nvals = len(item0)
        except Exception:
            nvals = 0

        if nvals < 2 or nvals > 4:
            raise AttributeError(
                "You can only call `measure` with a single item if that item "
                "is a Bond, Angle, Dihedral or Improper. Asking to measure "
                f"a single {item0} is not supported."
            )

        items = [i for i in item0]
        return measure(*items, improper_angle=improper_angle, map=map)

    from .maths import Vector

    def _to_coords(item, map):
        if item is None:
            return item
        else:
            try:
                return item.coordinates(map=map)
            except Exception:
                pass

            return Vector.to_vector(item)


    def _get_space(items, map):
        """Return the first space property that we find from the
           passed items
        """
        try:
            space_property = map["space"]
        except Exception:
            space_property = "space"

        for item in items:
            try:
                # does this naturally have a space property?
                return item.property(space_property)
            except Exception:
                pass

            try:
                # does the molecule container of this view have a space
                # property?
                return item.molecule().property(space_property)
            except Exception:
                pass

            try:
                # maybe this was a multi-molecule container?
                # In which case, ask the first item for its space property
                return item[0].molecule().property(space_property)
            except Exception:
                pass

        return None

    # try to find a space from these objects
    if ignore_space:
        space = None
    else:
        space = _get_space([item0, item1, item2, item3], map=map)

    item0 = _to_coords(item0, map=map)
    item1 = _to_coords(item1, map=map)
    item2 = _to_coords(item2, map=map)
    item3 = _to_coords(item3, map=map)

    if space is not None:
        # map each point into the same space as item0 - this should
        # mean that the minimum image distance / angle will be calculated
        if item1 is not None:
            item1 = space.get_minimum_image(item1, item0)

        if item2 is not None:
            item2 = space.get_minimum_image(item2, item0)

        if item3 is not None:
            item3 = space.get_minimum_image(item3, item0)


    if item3 is None:
        if item2 is None:
            return Vector.distance(item0, item1)
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
