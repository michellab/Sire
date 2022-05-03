"""
.. currentmodule:: sire.io

"""

from .. import system as _system
from .. import mol as _mol

_mol_version = _mol.__file__

from ..legacy import IO as _IO

from .. import use_new_api as _use_new_api
_use_new_api()


def load_molecules(*args, **kwargs):
    from ..legacy.IO import load_molecules as _load_molecules
    from ..system import System
    return System(_load_molecules(*args, **kwargs))
