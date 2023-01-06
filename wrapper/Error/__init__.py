
from .. import Qt as _Qt

from ._Error import *

__old_printError = printError

def printError(e):
    __old_printError( str(e) )

def get_back_trace():
    """Print the current backtrace (including the C++ part)"""
    return getBackTrace()

def print_error(e):
    """Print the passed error"""
    printError(e)

# we are going to start by disabling the exception backtraces.
# These are only really needed by developers, and slow some things
# down. If you are a developer, than call
# enable_backtrace_exceptions() anywhere to re-enable them
disable_backtrace_exceptions()
