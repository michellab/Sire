from ..legacy.Maths import Sphere as _Sphere


def _fix_sphere():
    def radius(obj):
        from ..units import angstrom

        return obj.__old__radius() * angstrom

    _Sphere.__old__radius = _Sphere.radius
    _Sphere.radius = radius

    def volume(obj):
        from ..units import angstrom3

        return obj.__old__volume() * angstrom3

    _Sphere.__old__volume = _Sphere.volume
    _Sphere.volume = volume

    def surface_area(obj):
        from ..units import angstrom2

        return obj.__old__surface_area() * angstrom2

    try:
        _Sphere.__old__surface_area = _Sphere.surfaceArea
        delattr(_Sphere, "surfaceArea")
    except AttributeError:
        _Sphere.__old__surface_area = _Sphere.surface_area

    _Sphere.surface_area = surface_area

    def translate(obj, delta):
        from ._vector import Vector

        return obj.__old__translate(Vector(delta))

    _Sphere.__old__translate = _Sphere.translate
    _Sphere.translate = translate

    def set_radius(obj, radius):
        from ..units import angstrom

        # get the default unit of length
        length_unit = angstrom.get_default()

        def _is_number(v):
            return isinstance(v, int) or isinstance(v, float)

        if _is_number(radius):
            radius = radius * length_unit

        return obj.__old__set_radius(radius.to(angstrom))

    try:
        _Sphere.__old__set_radius = _Sphere.setRadius
        delattr(_Sphere, "setRadius")
    except AttributeError:
        _Sphere.__old__set_radius = _Sphere.set_radius

    _Sphere.set_radius = set_radius

    def intersection_volume(obj, other):
        from ..units import angstrom3

        return obj.__old__intersection_volume(other) * angstrom3

    try:
        _Sphere.__old__intersection_volume = _Sphere.intersectionVolume
        delattr(_Sphere, "intersectionVolume")
    except AttributeError:
        _Sphere.__old__intersection_volume = _Sphere.intersection_volume

    _Sphere.intersection_volume = intersection_volume

    @staticmethod
    def combined_volume(spheres):
        from ..units import angstrom3

        if not hasattr(spheres, "__len__"):
            spheres = [spheres]

        return _Sphere.__old__combined_volume(spheres) * angstrom3

    try:
        _Sphere.__old__combined_volume = _Sphere.combinedVolume
        delattr(_Sphere, "combinedVolume")
    except AttributeError:
        _Sphere.__old__combined_volume = _Sphere.combined_volume

    _Sphere.combined_volume = combined_volume

    @staticmethod
    def combined_volume_mc(spheres, nsamples=-1):
        from ..units import angstrom3

        if not hasattr(spheres, "__len__"):
            spheres = [spheres]

        return _Sphere.__old__combined_volume_mc(spheres, nsamples) * angstrom3

    try:
        _Sphere.__old__combined_volume_mc = _Sphere.combinedVolumeMC
        delattr(_Sphere, "combinedVolumeMC")
    except AttributeError:
        _Sphere.__old__combined_volume_mc = _Sphere.combined_volume_mc

    _Sphere.combined_volume_mc = combined_volume_mc

    def __str__(obj):
        return f"Sphere( center={str(obj.center())} radius={obj.radius()} )"

    _Sphere.__str__ = __str__
    _Sphere.__repr__ = __str__


class Sphere(_Sphere):
    """
    A sphere in 3D space
    """

    def __init__(self, *args, **kwargs):
        from ..units import angstrom

        # get the default unit of length
        length_unit = angstrom.get_default()

        def _is_number(v):
            return isinstance(v, int) or isinstance(v, float)

        def _is_vector(v):
            try:
                return len(v) == 3
            except Exception:
                return False

        # mix of doubles and lengths?
        new_args = []

        for i in range(0, len(args)):
            if _is_number(args[i]):
                new_args.append((args[i] * length_unit).to(angstrom))
            elif _is_vector(args[i]):
                from ._vector import Vector

                new_args.append(Vector(args[i]))
            elif hasattr(args[i], "to"):
                new_args.append(args[i].to(angstrom))
            else:
                new_args.append(args[i])

        for key in kwargs.keys():
            if _is_number(kwargs[key]):
                kwargs[key] = (kwargs[key] * length_unit).to(angstrom)
            elif _is_vector(kwargs[key]):
                from ._vector import Vector

                kwargs[key] = Vector(args[i])
            elif hasattr(kwargs[key], "to"):
                kwargs[key] = kwargs[key].to(angstrom)

        super().__init__(*new_args, **kwargs)


if not hasattr(_Sphere, "__old__radius"):
    _fix_sphere()
