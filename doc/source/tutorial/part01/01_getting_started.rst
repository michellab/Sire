===============
Getting Started
===============

Sire is Python library that is designed to make it easier for you
to build and manipulate molecular systems.

This tutorial assumes that you have installed ``Sire`` following
:doc:`the instructions here. <../../install>`, and have activated
the Anaconda / miniconda Python environment into which Sire was
installed.

.. warning::

  This tutorial is written for ``metawards`` version |SireVersion| or
  higher. If you are using an older version then please upgrade.

Importing Sire
--------------

You import ``Sire`` by typing

```{python}
import Sire
```

As a convention, we will import ``Sire`` under the alias `sr` to reduce
the amount of typing. We will use `sr` throughout this tutorial to
mean `Sire`.

```{python}
import Sire as sr
```

Loading a molecule
------------------

We load molecules using the :fun:`Sire.load` function. This accepts either
a filename, a URL, or a `PDB code<pdb_code>`.

For example, let's load the XXX.

```{python}
mol = sr.load("https://path/to/something.sdf")
```

Supports lots of file formats.

Saving a molecule
-----------------


