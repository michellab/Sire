"""
.. currentmodule:: sire.base

"""

from ..legacy import Base as _Base

from .. import use_new_api as _use_new_api
_use_new_api()

Property = _Base.Property
Properties = _Base.Properties

