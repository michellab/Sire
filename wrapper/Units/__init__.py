"""
.. currentmodule:: Sire.Units

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
from  Sire.Units._Units import *

__all__ = [ "Celsius", "Fahrenheit", "GeneralUnit", "Unit",
            "acute", "convert", "convertFrom", "convertTo" ]

import Sire.Base

wrap = Sire.Base._add_wrap_function(wrap)

