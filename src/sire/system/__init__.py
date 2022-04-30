"""
.. currentmodule:: sire.system

"""

from .. import mol as _mol
from ..legacy import System as _System

from ..utils import pythonize_module as _pythonize_module
_pythonize_module(_System)

from ._system import *
