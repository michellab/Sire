===========
Sire.Stream
===========

Nearly all Sire objects can be streamed (marshalled) to and from a portable, versioned
and compact binary format. This is interfaced to
`Python's pickle module <https://docs.python.org/3/library/pickle.html>`__,
meaning that nearly all Sire objects can be safely pickled, e.g.

.. code-block:: python

    >>> s = pickle.dumps(mol)
    >>> new_mol = pickle.loads(s)

You wouldn't normally need to use any of the functionality in this module
directly.

.. toctree::
   :maxdepth: 3

   index_api_Sire_Stream
