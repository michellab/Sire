==========
Sire.Units
==========

This module implements the dimensionally-consistent units system
that is used throughout Sire. Most of Sire takes dimensioned values.
These support type checking (e.g. making sure that we don't pass
in a length to a function that expects a time), plus automatic
unit conversion between a range of commonly-used units.

For example;

.. code-block:: python

    >>> from Sire.Units import meter, second, kilogram, joule, kcal, mole
    >>> energy = 5000 * kilogram * (meter / second)**2 / mole
    >>> print(energy.to(joule / mole))
    5000.000000000001

    >>> print(energy.to(kcal / mole))
    1.1950286806883368

.. toctree::
   :maxdepth: 3

   index_api_Sire_Units
