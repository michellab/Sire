"""
.. currentmodule:: sire.mm

"""

from ..legacy import MM as _MM

from .. import use_new_api as _use_new_api
_use_new_api()


def calculate_energy(*args, **kwargs):
    from ..mol import _to_molecules

    new_args = []
    new_kwargs = {}

    for arg in args:
        try:
            new_args.append(_to_molecules(arg))
        except Exception:
            new_args.append(arg)

    for key, value in kwargs.items():
        try:
            new_kwargs[key] = _to_molecules(value)
        except Exception:
            new_kwargs[key] = value

    return _MM.calculate_energy(*new_args, **new_kwargs)


def create_forcefield(*args, map=None, **kwargs):
    from ..mol import _to_molecules

    new_args = []
    new_kwargs = {}

    for arg in args:
        try:
            new_args.append(_to_molecules(arg))
        except Exception:
            new_args.append(arg)

    for key, value in kwargs.items():
        try:
            new_kwargs[key] = _to_molecules(value)
        except Exception:
            new_kwargs[key] = value

    if map is None:
        map = {}

    new_kwargs["map"] = map

    return _MM.create_forcefield(*new_args, **new_kwargs)


Bond = _MM.Bond
SelectorBond = _MM.SelectorBond
SelectorMBond = _MM.SelectorMBond

Angle = _MM.Angle
SelectorAngle = _MM.SelectorAngle
SelectorMAngle = _MM.SelectorMAngle

Dihedral = _MM.Dihedral
SelectorDihedral = _MM.SelectorDihedral
SelectorMDihedral = _MM.SelectorMDihedral

Improper = _MM.Improper
SelectorImproper = _MM.SelectorImproper
SelectorMImproper = _MM.SelectorMImproper

try:
    Bond.__len__ = Bond.nAtoms
except AttributeError:
    Bond.__len__ = Bond.num_atoms

try:
    Angle.__len__ = Angle.nAtoms
except AttributeError:
    Angle.__len__ = Angle.num_atoms

try:
    Dihedral.__len__ = Dihedral.nAtoms
except AttributeError:
    Dihedral.__len__ = Dihedral.num_atoms

try:
    Improper.__len__ = Improper.nAtoms
except AttributeError:
    Improper.__len__ = Improper.num_atoms

from ..mol import __fix_getitem, _add_evals, _add_property_func, \
    _add_apply_func

for C in [Bond, SelectorBond, SelectorMBond,
          Angle, SelectorAngle, SelectorMAngle,
          Dihedral, SelectorDihedral, SelectorMDihedral,
          Improper, SelectorImproper, SelectorMImproper]:
    __fix_getitem(C)
    _add_evals(C)
    _add_property_func(C)
    _add_apply_func(C)

from ..mol import _cursor, _cursors, _cursorsm, _trajectory, _viewfunc

Bond.cursor = _cursor
SelectorBond.cursor = _cursors
Angle.cursor = _cursor
SelectorAngle.cursor = _cursors
Dihedral.cursor = _cursor
SelectorDihedral.cursor = _cursors
Improper.cursor = _cursor
SelectorImproper.cursor = _cursors

SelectorMBond.trajectory = _trajectory
SelectorMAngle.trajectory = _trajectory
SelectorMDihedral.trajectory = _trajectory
SelectorMImproper.trajectory = _trajectory

SelectorBond.view = _viewfunc
SelectorAngle.view = _viewfunc
SelectorDihedral.view = _viewfunc
SelectorImproper.view = _viewfunc

SelectorMBond.view = _viewfunc
SelectorMAngle.view = _viewfunc
SelectorMDihedral.view = _viewfunc
SelectorMImproper.view = _viewfunc

SelectorMBond.cursor = _cursorsm
SelectorMAngle.cursor = _cursorsm
SelectorMDihedral.cursor = _cursorsm
SelectorMImproper.cursor = _cursorsm
