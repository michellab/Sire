"""
.. currentmodule:: sire.system

"""

from ..legacy import System as _System
from .. import use_new_api as _use_new_api
_use_new_api()

from .. import mol as _mol

from ._system import *
