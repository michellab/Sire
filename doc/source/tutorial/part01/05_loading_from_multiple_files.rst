===========================
Loading from multiple files
===========================

It is often the case that molecular information needs to be read from
multiple files, e.g. a separate topology and coordinate file.

You load from multiple files simply by passing multiple filenames and/or
URLs to :func:`sire.load`.

>>> mols = sr.load("https://siremol.org/m/ala.top",
...                "https://siremol.org/m/ala.crd")
Downloading from 'https://siremol.org/m/ala.top'...
Unzipping './ala.top.bz2'...
Downloading from 'https://siremol.org/m/ala.crd'...
Unzipping './ala.crd.bz2'...

>>> print(mols)
System( name=ACE num_molecules=631 num_residues=633 num_atoms=1912 )

You can pass in the filenames as multiple arguments or as a list,
whichever you find easiest.

>>> mols = sr.load(["https://siremol.org/m/ala.top",
...                 "https://siremol.org/m/ala.crd"])

If the files or URLs have a common base, then you can save some typing
by using :func:`sire.expand`, e.g.

>>> mols = sr.load(sr.expand("https://siremol.org/m",
...                          "ala.top", "ala.crd"))

or

>>> mols = sr.load(sr.expand(sr.tutorial_url,
...                          ["ala.top", "ala.crd"]))

.. note::

   ``sr.tutorial_url`` expands to the base URL for tutorial files
   (https://siremol.org/m). It is worth using this variable for
   the tutorial as it auto-completes and will reduce errors.

If you are loading files, you can also make use of glob expressions
(wildcard expansions), e.g.

>>> mols = sr.load("ala.???")

.. note::

   This line loads the ``ala.top`` and ``ala.crd`` files that
   were downloaded by the above lines. This is because Sire downloads
   files at URLs to the current directory. You can tell it to use
   a different directory by passing that in via the ``directory``
   argument, e.g. ``sr.load(sr.expand(sr.tutorial_url,"cholesterol.sdf"), directory="tmp")``.
   The directory will be created automatically if it doesn't exist.

.. note::

   We couldn't use ``ala.*`` because the directory contains the compressed
   input files, ``ala.crd.bz2`` and ``ala.top.bz2``, which were downloaded
   by ``sr.load``. If you remove both these files, then you could
   use ``sr.load("ala.*")``.
