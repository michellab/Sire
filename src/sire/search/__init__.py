"""
.. currentmodule:: sire.search

"""

from ..legacy import Search as _Search

from .. import use_new_api as _use_new_api
_use_new_api()

approx_equal = _Search.approx_equal
set_approx_epsilon = _Search.set_approx_epsilon
get_approx_epsilon = _Search.get_approx_epsilon
