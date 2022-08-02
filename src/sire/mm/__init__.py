"""
.. currentmodule:: sire.mm

"""

from ..legacy import MM as _MM

from .. import use_new_api as _use_new_api
_use_new_api()

Bond = _MM.Bond
SelectorBond = _MM.SelectorBond
SelectorMBond = _MM.SelectorMBond

try:
    Bond.__len__ = Bond.nAtoms
except AttributeError:
    Bond.__len__ = Bond.num_atoms

from ..mol import __fix_getitem, _add_evals, _add_property_func, \
    _add_apply_func

for C in [Bond, SelectorBond, SelectorMBond]:
    __fix_getitem(C)
    _add_evals(C)
    _add_property_func(C)
    _add_apply_func(C)
