"""
.. currentmodule:: sire.mm

"""

from ..legacy import MM as _MM

from .. import use_new_api as _use_new_api
_use_new_api()

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

from sire.mol import _cursor, _cursors

Bond.cursor = _cursor
SelectorBond.cursor = _cursors
Angle.cursor = _cursor
SelectorAngle.cursor = _cursors
Dihedral.cursor = _cursor
SelectorDihedral.cursor = _cursors
Improper.cursor = _cursor
SelectorImproper.cursor = _cursors

from ..mol import _trajectory

SelectorMBond.trajectory = _trajectory
SelectorMAngle.trajectory = _trajectory
SelectorMDihedral.trajectory = _trajectory
SelectorMImproper.trajectory = _trajectory
