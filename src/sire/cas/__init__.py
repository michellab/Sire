"""
.. currentmodule:: sire.cas

"""

from ..legacy import CAS as _CAS

from .. import use_new_api as _use_new_api
_use_new_api()

Symbol = _CAS.Symbol
Expression = _CAS.Expression
