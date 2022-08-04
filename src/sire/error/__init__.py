"""
.. currentmodule:: sire.error

"""

from ..legacy import Error as _Error

from .. import use_new_api as _use_new_api
_use_new_api()

from ..legacy.Error import get_last_error_details, \
    disable_backtrace_exceptions \
    enable_backtrace_exceptions


