==========
Sire.Error
==========

This module is used to convert between Sire's C++ and Python exceptions.
There are some key functions that are useful if you want more detail
about an error or issue that has occurred.

:func:`~sire.error.get_last_error_details`
    Get more information about the last C++ Sire exception that was
    raised. This includes information like the C++ backtrace from
    where the exception was raised, and the line of C++ code that
    triggered the issue.

:func:`~sire.error.get_back_trace`
    Get the full current backtrace.

.. toctree::
   :maxdepth: 3

   index_api_Sire_Error
