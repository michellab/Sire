"""
.. currentmodule:: Sire.Error

This module is used to map between Sire's C++ exceptions and
standard Python exceptions.


Functions
=========

.. autosummary::
    :toctree: generated/

    get_back_trace
    print_error
    get_last_error_details

"""

import Sire.Qt

from Sire.Error._Error import *

__all__ = ["get_back_trace", "print_error", "get_last_error_details"]

__old_printError = printError

def printError(e):
    __old_printError( str(e) )

def get_back_trace():
    """Print the current backtrace (including the C++ part)"""
    return getBackTrace()

def print_error(e):
    """Print the passed error"""
    printError(e)
