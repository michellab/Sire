==============
Indexing Atoms
==============

You can get an atom from any container by indexing it, as you would
a Python list.

>>> print(mol[0])
Atom( N:1     [ -54.07,   11.27,   41.93] )

would get the first atom in the container, while

>>> print(mol[-1])
Atom( O:11732 [   4.73,    5.03,   46.82] )

gets the last atom.

Passing in an index of an atom that doesn't exist results in an
``IndexError`` exception being raised, e.g.

>>> print(mol[100000])
---------------------------------------------------------------------------
IndexError                                Traceback (most recent call last)
Input In [32], in <cell line: 1>()
----> 1 print(mol[100000])
<BLANKLINE>
File ~/sire.app/lib/python3.8/site-packages/Sire/Mol/__init__.py:414, in __fixed__getitem__(obj, key)
    412         return obj.residue(key)
    413     else:
--> 414         return obj.atom(key)
    415 else:
    416     if __is_chain_class(obj):
<BLANKLINE>
IndexError: SireError::invalid_index: No item at index 100000. Index range is from -11728 to 11727. (call Sire.Error.get_last_error_details() for more info)

You can request multiple atoms by using a slice, e.g.

>>> print(mol[0:10])
Selector<SireMol::Atom>( size=10
0:  Atom( N:1     [ -54.07,   11.27,   41.93] )
1:  Atom( CA:2    [ -55.43,   11.35,   42.54] )
2:  Atom( C:3     [ -56.06,    9.95,   42.55] )
3:  Atom( O:4     [ -57.04,    9.73,   41.82] )
4:  Atom( CB:5    [ -56.32,   12.33,   41.76] )
5:  Atom( CG1:6   [ -55.68,   13.72,   41.72] )
6:  Atom( CG2:7   [ -57.70,   12.40,   42.39] )
7:  Atom( CD1:8   [ -55.42,   14.31,   43.09] )
8:  Atom( N:9     [ -55.50,    9.04,   43.36] )
9:  Atom( CA:10   [ -56.02,    7.64,   43.47] )
)

gives the first 10 atoms, while

>>> print(mol[9::-1])
Selector<SireMol::Atom>( size=10
0:  Atom( CA:10   [ -56.02,    7.64,   43.47] )
1:  Atom( N:9     [ -55.50,    9.04,   43.36] )
2:  Atom( CD1:8   [ -55.42,   14.31,   43.09] )
3:  Atom( CG2:7   [ -57.70,   12.40,   42.39] )
4:  Atom( CG1:6   [ -55.68,   13.72,   41.72] )
5:  Atom( CB:5    [ -56.32,   12.33,   41.76] )
6:  Atom( O:4     [ -57.04,    9.73,   41.82] )
7:  Atom( C:3     [ -56.06,    9.95,   42.55] )
8:  Atom( CA:2    [ -55.43,   11.35,   42.54] )
9:  Atom( N:1     [ -54.07,   11.27,   41.93] )
)

gives the first 10 atoms in reverse order.

The returned :class:`~Sire.Mol.Selector_Atom_` is itself a
molecular container, which has all of the same indexing and
other functions as :class:`~Sire.Mol.Atom` et al.

>>> print(mol[0:10][-1:4:-1])
Selector<SireMol::Atom>( size=5
0:  Atom( CA:10   [ -56.02,    7.64,   43.47] )
1:  Atom( N:9     [ -55.50,    9.04,   43.36] )
2:  Atom( CD1:8   [ -55.42,   14.31,   43.09] )
3:  Atom( CG2:7   [ -57.70,   12.40,   42.39] )
4:  Atom( CG1:6   [ -55.68,   13.72,   41.72] )
)

gets the 5th to 9th atoms in reverse order.

You can also pass in a list of indicies of atoms.

>>> print(mol[[0, 2, 4, 6, 8]])
Selector<SireMol::Atom>( size=5
0:  Atom( N:1     [ -54.07,   11.27,   41.93] )
1:  Atom( C:3     [ -56.06,    9.95,   42.55] )
2:  Atom( CB:5    [ -56.32,   12.33,   41.76] )
3:  Atom( CG2:7   [ -57.70,   12.40,   42.39] )
4:  Atom( N:9     [ -55.50,    9.04,   43.36] )
)

This gives the atoms at indicies 0, 2, 4, 6, and 8.

This indexing is a shorthand for calling the :func:`~Sire.Mol.Atom.atom`
(for single) and :func:`~Sire.Mol.Atom.atoms` (for multiple) atom functions.
You can call these directly, e.g.

>>> print(mol.atom(0))
Atom( N:1     [ -54.07,   11.27,   41.93] )

>>> print(mol.atoms([0, 2, 4, 6, 8]))
Selector<SireMol::Atom>( size=5
0:  Atom( N:1     [ -54.07,   11.27,   41.93] )
1:  Atom( C:3     [ -56.06,    9.95,   42.55] )
2:  Atom( CB:5    [ -56.32,   12.33,   41.76] )
3:  Atom( CG2:7   [ -57.70,   12.40,   42.39] )
4:  Atom( N:9     [ -55.50,    9.04,   43.36] )
)

The :func:`~Sire.Mol.Atom.atom` function is guaranteed to always return
either a single atom, or raise an exception if this is not possible.

The :func:`~Sire.Mol.Atom.atoms` function will return either a single
atom or a selector of atoms, or raise an exception if this is not possible.

Accessing by name
-----------------

So far we have been accessing the molecular containers by index, as if
they were Python lists. We can also treat the molecular containers
like Python dictionaries, and get atoms by their name.

For example, to get the atoms called "C" we would use

>>> print(mol["C"])
Selector<SireMol::Atom>( size=1494
0:  Atom( C:3     [ -56.06,    9.95,   42.55] )
1:  Atom( C:11    [ -56.14,    7.05,   42.06] )
2:  Atom( C:18    [ -54.61,    4.90,   40.03] )
3:  Atom( C:26    [ -54.80,    2.14,   38.44] )
4:  Atom( C:35    [ -53.58,   -0.34,   36.88] )
...
1489:  Atom( C:11625 [  23.84,   -6.28,   31.63] )
1490:  Atom( C:11630 [  26.57,   -5.11,   30.71] )
1491:  Atom( C:11637 [  28.49,   -2.95,   31.83] )
1492:  Atom( C:11646 [  30.02,   -0.22,   31.11] )
1493:  Atom( C:11654 [  32.09,   -0.82,   34.12] )
)

Note that there are multiple atoms in this molecule called "C", hence
several are returned. This would raise an exception if you called
the shorthand :func:`~Sire.Mol.Atom.atom` function with this name,

>>> print(mol.atom("C"))
---------------------------------------------------------------------------
KeyError                                  Traceback (most recent call last)
Input In [10], in <cell line: 1>()
----> 1 mol.atom("C")
<BLANKLINE>
KeyError: "SireMol::duplicate_atom: More than one atom matches the ID
AtomName('C') (number of matches is 1494).
(call Sire.Error.get_last_error_details() for more info)"

A ``KeyError`` exception has been raised because there are multiple
atoms in this protein that are called "C" and Sire does not know which
one you want.

In this case, you would have to use the shorthand
:func:`~Sire.Mol.Atom.atoms` function.

>>> print(mol.atoms("C"))
Selector<SireMol::Atom>( size=1494
0:  Atom( C:3     [ -56.06,    9.95,   42.55] )
1:  Atom( C:11    [ -56.14,    7.05,   42.06] )
2:  Atom( C:18    [ -54.61,    4.90,   40.03] )
3:  Atom( C:26    [ -54.80,    2.14,   38.44] )
4:  Atom( C:35    [ -53.58,   -0.34,   36.88] )
...
1489:  Atom( C:11625 [  23.84,   -6.28,   31.63] )
1490:  Atom( C:11630 [  26.57,   -5.11,   30.71] )
1491:  Atom( C:11637 [  28.49,   -2.95,   31.83] )
1492:  Atom( C:11646 [  30.02,   -0.22,   31.11] )
1493:  Atom( C:11654 [  32.09,   -0.82,   34.12] )
)

.. note::

    Using the index operator (``mol["C"]``) is easiest, as it will always
    do the right thing. Use the :func:`~Sire.Mol.Atom.atom` and
    :func:`~Sire.Mol.Atom.atoms` functions only when you want to
    ensure that the container will return atoms.

As before, the returned :class:`~Sire.Mol.Selector_Atom_` is itself a container,
and so also has its own ``.atom()``, ``.atoms()`` and indexing functions, e.g.

>>> print(mol["C"][0])
Atom( C:3     [ -56.06,    9.95,   42.55] )

gives the atom at the index 0 in the container of atoms that are called "C",
and

>>> print(mol("C")[-1])
Atom( C:11654 [  32.09,   -0.82,   34.12] )

gives the last atom in the container of atoms that are called "C".

Asking for an atom that doesn't exist will result in a ``KeyError``
exception being raised.

>>> print(mol["X"])
---------------------------------------------------------------------------
KeyError                                  Traceback (most recent call last)
Input In [24], in <cell line: 1>()
----> 1 print(mol["X"])
<BLANKLINE>
File ~/sire.app/lib/python3.8/site-packages/Sire/Mol/__init__.py:419, in __fixed__getitem__(obj, key)
    417     return obj.residues(key)
    418 else:
--> 419     return obj.atoms(key)
<BLANKLINE>
File ~/sire.app/lib/python3.8/site-packages/Sire/Mol/__init__.py:428, in __fixed__atoms__(obj, idx)
    426     return obj.__orig__atoms(list(idx))
    427 else:
--> 428     return obj.__orig__atoms(idx)
<BLANKLINE>
KeyError: 'SireMol::missing_atom: There is no atom called "X" in the layout "{c4d51f89-f4f7-4e0c-854d-da27affe1baf}". (call Sire.Error.get_last_error_details() for more info)'

Searching by Name and Number
---------------------------

So far, we have identified atoms by their index or by their name. There is
a third way to identify atoms. This is by their atom number. The atom number
is an identifying number that is assigned to atoms, typically via an
atom number column in the input file. We can look up atoms by number
using a search string.

>>> print(mol["atomnum 1"])
Atom( N:1     [ -54.07,   11.27,   41.93] )

Here, ``atomnum 1`` is a search string that looks for atoms with number 1.

.. note::

   Note how the name ``N`` and number ``1`` of the atom are printed
   out, along with its coordinates.

This search string is very powerful. For example, you can search for
atoms that have a number that is less than 10.

>>> print(mol["atomnum < 10"])
Selector<SireMol::Atom>( size=9
0:  Atom( N:1     [ -54.07,   11.27,   41.93] )
1:  Atom( CA:2    [ -55.43,   11.35,   42.54] )
2:  Atom( C:3     [ -56.06,    9.95,   42.55] )
3:  Atom( O:4     [ -57.04,    9.73,   41.82] )
4:  Atom( CB:5    [ -56.32,   12.33,   41.76] )
5:  Atom( CG1:6   [ -55.68,   13.72,   41.72] )
6:  Atom( CG2:7   [ -57.70,   12.40,   42.39] )
7:  Atom( CD1:8   [ -55.42,   14.31,   43.09] )
8:  Atom( N:9     [ -55.50,    9.04,   43.36] )
)

You can combine search strings with logical operators (e.g. ``and``,
``or`` and ``not``).

>>> print(mol["atomnum >= 5 and atomnum < 10"])
Selector<SireMol::Atom>( size=5
0:  Atom( CB:5    [ -56.32,   12.33,   41.76] )
1:  Atom( CG1:6   [ -55.68,   13.72,   41.72] )
2:  Atom( CG2:7   [ -57.70,   12.40,   42.39] )
3:  Atom( CD1:8   [ -55.42,   14.31,   43.09] )
4:  Atom( N:9     [ -55.50,    9.04,   43.36] )
)

and can also search for multiple atom numbers

>>> print(mol["atomnum 1, 3, 5, 7"])
Selector<SireMol::Atom>( size=4
0:  Atom( N:1     [ -54.07,   11.27,   41.93] )
1:  Atom( C:3     [ -56.06,    9.95,   42.55] )
2:  Atom( CB:5    [ -56.32,   12.33,   41.76] )
3:  Atom( CG2:7   [ -57.70,   12.40,   42.39] )
)

You can include atom names (``atomname``) in the search string.

>>> print(mol["atomname CA"])
Selector<SireMol::Atom>( size=1494
0:  Atom( CA:2    [ -55.43,   11.35,   42.54] )
1:  Atom( CA:10   [ -56.02,    7.64,   43.47] )
2:  Atom( CA:17   [ -54.99,    6.39,   39.98] )
3:  Atom( CA:25   [ -55.33,    2.58,   39.80] )
4:  Atom( CA:34   [ -52.97,    1.03,   37.19] )
...
1489:  Atom( CA:11624 [  22.43,   -6.30,   32.21] )
1490:  Atom( CA:11629 [  25.36,   -5.51,   29.89] )
1491:  Atom( CA:11636 [  27.51,   -3.84,   32.59] )
1492:  Atom( CA:11645 [  28.74,   -0.85,   30.58] )
1493:  Atom( CA:11653 [  31.65,   -0.00,   32.91] )
)

can search for multiple atom names

>>> print(mol["atomname CA, C, N"])
Selector<SireMol::Atom>( size=4482
0:  Atom( N:1     [ -54.07,   11.27,   41.93] )
1:  Atom( CA:2    [ -55.43,   11.35,   42.54] )
2:  Atom( C:3     [ -56.06,    9.95,   42.55] )
3:  Atom( N:9     [ -55.50,    9.04,   43.36] )
4:  Atom( CA:10   [ -56.02,    7.64,   43.47] )
...
4477:  Atom( CA:11645 [  28.74,   -0.85,   30.58] )
4478:  Atom( C:11646 [  30.02,   -0.22,   31.11] )
4479:  Atom( N:11652 [  30.40,   -0.51,   32.35] )
4480:  Atom( CA:11653 [  31.65,   -0.00,   32.91] )
4481:  Atom( C:11654 [  32.09,   -0.82,   34.12] )
)

and can use mixtures of any identifiers, e.g.

>>> print(mol["atomnum > 1000 and atomname N"])
Selector<SireMol::Atom>( size=1369
0:  Atom( N:1005  [ -30.70,  -11.17,   28.73] )
1:  Atom( N:1014  [ -32.11,  -10.01,   25.69] )
2:  Atom( N:1023  [ -32.43,  -11.41,   22.41] )
3:  Atom( N:1027  [ -32.33,  -11.69,   18.84] )
4:  Atom( N:1038  [ -35.13,  -10.12,   17.49] )
...
1364:  Atom( N:11623 [  22.09,   -7.64,   32.65] )
1365:  Atom( N:11628 [  24.07,   -5.44,   30.61] )
1366:  Atom( N:11635 [  26.39,   -4.38,   31.81] )
1367:  Atom( N:11644 [  28.07,   -1.72,   31.54] )
1368:  Atom( N:11652 [  30.40,   -0.51,   32.35] )
)

Uniquely identifying atoms
--------------------------

It is often useful to have a unique identifier for an atom in a molecule.
The name and number can't do this, as it is possible that many atoms
could have the same name, and/or the same number.

To solve this, Sire uses the index of the atom in the molecule as its
unique identifier. This index is called ``atomidx``, and can also be
used in searches.

>>> print(mol["atomidx 0"])
Atom( N:1     [ -54.07,   11.27,   41.93] )

This has printed the first atom in the molecule.

You can get an atom's index using the :func:`~Sire.Mol.Atom.index` function.

>>> print(mol["atomidx 0"].index())
AtomIdx(0)

.. warning::

    Be careful searching with ``atomidx``. This is the unique
    index of the atom within its parent molecule, not the index
    of the atom in a container. So ``mol[5:10]["atomidx 0"]`` would
    raise a KeyError as the first atom in the molecule is not
    contained in the slice of atoms 5 to 9.

Iterating over atoms
--------------------

The :class:`~Sire.Mol.Selector_Atom_` class is iterable, meaning that
it can be used in loops.

>>> for atom in mol["atomnum < 10"]:
...     print(atom)
Atom( N:1     [ -54.07,   11.27,   41.93] )
Atom( CA:2    [ -55.43,   11.35,   42.54] )
Atom( C:3     [ -56.06,    9.95,   42.55] )
Atom( O:4     [ -57.04,    9.73,   41.82] )
Atom( CB:5    [ -56.32,   12.33,   41.76] )
Atom( CG1:6   [ -55.68,   13.72,   41.72] )
Atom( CG2:7   [ -57.70,   12.40,   42.39] )
Atom( CD1:8   [ -55.42,   14.31,   43.09] )
Atom( N:9     [ -55.50,    9.04,   43.36] )
