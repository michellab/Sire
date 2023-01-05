__all__ = ["Expression", "lam", "Symbol", "x", "y"]

from ..legacy import CAS as _CAS

from .. import use_new_api as _use_new_api

_use_new_api()

Symbol = _CAS.Symbol
Expression = _CAS.Expression

lam = Symbol("lambda")
x = Symbol("x")
y = Symbol("y")
