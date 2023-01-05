__all__ = [
    "disable_backtrace_exceptions",
    "enable_backtrace_exceptions",
    "get_last_error_details",
]

from ..legacy import Error as _Error

from .. import use_new_api as _use_new_api

_use_new_api()

get_last_error_details = _Error.get_last_error_details
disable_backtrace_exceptions = _Error.disable_backtrace_exceptions
enable_backtrace_exceptions = _Error.enable_backtrace_exceptions
