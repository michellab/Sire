
__all__ = ["System"]

from ..legacy import System as _System
from .. import use_new_api as _use_new_api
_use_new_api()

from ._system import *
