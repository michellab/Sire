=======
Sire.IO
=======

This module provides all of the parsers which are used to load and
save molecular information from files. Sire comes with parsers for
many popular molecule file formats. The Sire parsers are written
using a simple design philosophy:

1. The parser should read as much information as possible from the file,
   ideally discarding nothing

2. The parser should be symmetric, meaning that it should be capable
   of writing any information that it can read.

These two principles mean that Sire can be used to convert between
molecule file formats. You can also use Sire to load a molecular
system from a file, make some edits, and then save it back to
the same file format.

.. note::
    Note that you rarely need to use any of the functions or classes
    from this module directly. The :func:`sire.load` and :func:`sire.save`
    functions will do nearly everything you need.

The package centers around a few core classes

:class:`~sire.io.MoleculeParser`
    Load one or more molecules from the specified file(s) or URL(s).

.. toctree::
   :maxdepth: 3

   index_api_Sire_IO
