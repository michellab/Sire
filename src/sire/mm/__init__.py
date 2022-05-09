"""
.. currentmodule:: sire.mm

"""

from ..legacy import MM as _MM

from .. import use_new_api as _use_new_api
_use_new_api()

Bond = _MM.Bond
SelectorBond = _MM.SelectorBond

from ..mol import __fix_getitem

__fix_getitem(SelectorBond)
