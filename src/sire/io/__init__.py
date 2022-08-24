"""
.. currentmodule:: sire.io

"""

from ..legacy import IO as _IO
from .. import use_new_api as _use_new_api
_use_new_api()

from .. import system as _system

def load_molecules(*args, **kwargs):
    from ..legacy.IO import load_molecules as _load_molecules
    from ..system import System
    return System(_load_molecules(*args, **kwargs))
