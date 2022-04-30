"""
.. currentmodule:: sire.legacy.Units

Classes
=======

.. autosummary::
    :toctree: generated/

    Celsius
    Fahrenheit
    GeneralUnit
    Unit

Functions
=========

.. autosummary::
    :toctree: generated/

    acute
    convert
    convertFrom
    convertTo

"""
from ._Units import *

from .. import Base as _Base

wrap = _Base._add_wrap_function(wrap)
