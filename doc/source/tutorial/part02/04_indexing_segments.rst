=================
Indexing Segments
=================

Segments are collections of atoms, sometimes non-contiguously collections
of atoms. They typically represent a user-defined segment within a molecule
or protein. An atom can only belong to one segment at a time. Atoms do
not need to be assigned to segments. Segments are implemented via the
:class:`~sire.mol.Segment` class, which is itself a molecular container
for :class:`~sire.mol.Atom` objects.

You can access the segments in a molecule container using the
:func:`~sire.mol.Segment.segment` and :func:`~sire.mol.Segment.segments` functions,
which are available on all of the molecular container types.

Not many molecules have named segments, so first let's load a molecule
that does.

>>> mols = sr.load(sr.expand(sr.tutorial_url, "alanin.psf"))
Downloading from 'https://siremol.org/m/alanin.psf'...
>>> mol = mols[0]

This molecule contains only a single segment called "MAIN".

>>> print(mol.segments())
Selector<SireMol::Segment>( size=1
0:  Segment( MAIN num_atoms=66 )
)
>>> print(mol.segment(0))
Segment( MAIN num_atoms=66 )
>>> print(mol.segment("MAIN"))
Segment( MAIN num_atoms=66 )

Search for segments
-------------------

You can search for segments using their name (``segname``) or their
index (``segidx``).

>>> print(mol.segments("segname MAIN"))
Selector<SireMol::Segment>( size=1
0:  Segment( MAIN num_atoms=66 )
)
>>> print(mol.segment("segidx 0"))
Segment( MAIN num_atoms=66 )

.. note::

   Unlike atoms and residues, segments do not have a number. They
   are identified only by their index in their parent molecule, or
   their name

You can do a segment search via the containers index operator too!

>>> print(mol["segname MAIN"])
Molecule( 2.137 : num_atoms=66, num_residues=12 )

.. note::

    Sire will automatically convert a result from a search string
    called via the index operator to the largest matching view.
    In this case, the single segment contains all of the atoms
    of the whole molecule. So Sire has converted the result up
    to the whole molecule view.

You can combine the search string with chain, residue and/or atom search
terms too.

>>> print(mol["segname MAIN and atomname C"])
Selector<SireMol::Atom>( size=11
0:  Atom( C:2 )
1:  Atom( C:8 )
2:  Atom( C:14 )
3:  Atom( C:20 )
4:  Atom( C:26 )
...
6:  Atom( C:38 )
7:  Atom( C:44 )
8:  Atom( C:50 )
9:  Atom( C:56 )
10:  Atom( C:62 )
)

>>> print(mol["segname MAIN and resname ACE"])
Residue( ACE:1   num_atoms=3 )

As for other types, you can search for multiple segment names using
a comma, and can do wildcard (glob) searching too!

>>> print(mol.segment("segname /M*/"))
Segment( MAIN num_atoms=66 )

Finding the atoms in a segment
------------------------------

Because both :class:`~sire.mol.Segment` and :class:`~sire.mol.Selector_Segment_`
are molecular containers, they also have their own
:func:`~sire.mol.Segment.atom` and :func:`~sire.mol.Segment.atoms` functions,
which behave as you would expect.

>>> print(mol["segname MAIN"].atoms("C"))
Selector<SireMol::Atom>( size=11
0:  Atom( C:2 )
1:  Atom( C:8 )
2:  Atom( C:14 )
3:  Atom( C:20 )
4:  Atom( C:26 )
...
6:  Atom( C:38 )
7:  Atom( C:44 )
8:  Atom( C:50 )
9:  Atom( C:56 )
10:  Atom( C:62 )
)

You can also use ``atoms in``, ``chains in`` or ``residues in`` to get the
atoms, residues or chains in a segment.

>>> print(mol["residues in segname MAIN"])
Selector<SireMol::Residue>( size=12
0:  Residue( ACE:1   num_atoms=3 )
1:  Residue( ALA:2   num_atoms=6 )
2:  Residue( ALA:3   num_atoms=6 )
3:  Residue( ALA:4   num_atoms=6 )
4:  Residue( ALA:5   num_atoms=6 )
...
7:  Residue( ALA:8   num_atoms=6 )
8:  Residue( ALA:9   num_atoms=6 )
9:  Residue( ALA:10  num_atoms=6 )
10:  Residue( ALA:11  num_atoms=6 )
11:  Residue( CBX:12  num_atoms=3 )
)

>>> print(mol["atoms in segname MAIN"])
Selector<SireMol::Atom>( size=66
0:  Atom( CA:1 )
1:  Atom( C:2 )
2:  Atom( O:3 )
3:  Atom( N:4 )
4:  Atom( H:5 )
...
61:  Atom( C:62 )
62:  Atom( O:63 )
63:  Atom( N:64 )
64:  Atom( H:65 )
65:  Atom( CA:66 )
)

A ``KeyError`` will be raised if there are no residues or chains within
a segment, e.g.

>>> print(mol["chains within segname MAIN"])
---------------------------------------------------------------------------
KeyError                                  Traceback (most recent call last)
Input In [24], in <cell line: 1>()
----> 1 print(mol["chains in segname MAIN"])
<BLANKLINE>
File ~/sire.app/lib/python3.8/site-packages/Sire/Mol/__init__.py:462, in __fixed__getitem__(obj, key)
    458 elif type(key) is str:
    459     # is this a search object - if so, then return whatever is
    460     # most relevant from the search
    461     try:
--> 462         return __from_select_result(obj.search(key))
    463     except SyntaxError:
    464         pass
<BLANKLINE>
KeyError: 'SireMol::missing_chain: This view does not contain any chains. (call Sire.Error.get_last_error_details() for more info)'

You can go to segments from atoms or residues using ``segments with``, e.g.

>>> print(mol["segments with atomname C"])
Molecule( 2.137 : num_atoms=66, num_residues=12 )

Finding the atoms, residues or chains in a segment
--------------------------------------------------

Like all molecular containers, you can find the contained atoms,
residues or chains by calling the appropriate functions;

>>> print(mol["segname MAIN"].atoms())
Selector<SireMol::Atom>( size=66
0:  Atom( CA:1 )
1:  Atom( C:2 )
2:  Atom( O:3 )
3:  Atom( N:4 )
4:  Atom( H:5 )
...
61:  Atom( C:62 )
62:  Atom( O:63 )
63:  Atom( N:64 )
64:  Atom( H:65 )
65:  Atom( CA:66 )
)

>>> print(mol["segidx 0"].residues())
Selector<SireMol::Residue>( size=12
0:  Residue( ACE:1   num_atoms=3 )
1:  Residue( ALA:2   num_atoms=6 )
2:  Residue( ALA:3   num_atoms=6 )
3:  Residue( ALA:4   num_atoms=6 )
4:  Residue( ALA:5   num_atoms=6 )
...
7:  Residue( ALA:8   num_atoms=6 )
8:  Residue( ALA:9   num_atoms=6 )
9:  Residue( ALA:10  num_atoms=6 )
10:  Residue( ALA:11  num_atoms=6 )
11:  Residue( CBX:12  num_atoms=3 )
)

Uniquely identifying a segment
------------------------------

You uniquely identify a segment in a molecule using its segment index
(``segidx``). You can get the index of a segment in a molecule by
calling its :func:`~sire.mol.Segment.index` function.

>>> print(mol.segment(0).index())
SegIdx(0)

.. warning::

    Be careful indexing by segment index. This is the index of the segment
    that uniquely identifies it within its parent molecule. It is not the
    index of the segment in an arbitrary molecular container.

Segment identifying types
-------------------------

Another way to index segments is to use the segment identifying types, i.e.
:class:`~sire.mol.SegName` and :class:`~sire.mol.SegIdx`. The
easiest way to create these is by using the function
:func:`sire.segid`.

Use strings to create :class:`~sire.mol.SegName` objects,

>>> print(sr.segid("MAIN"))
SegName('MAIN')
>>> print(mol[sr.segid("MAIN")])
Segment( MAIN num_atoms=66 )

and integers to create :class:`~sire.mol.SegIdx` objects.

>>> print(sr.segid(0))
SegIdx(0)
>>> print(mol[sr.segid(0)])
Segment( MAIN num_atoms=66 )

You can set both a name and an index by passing in two arguments.

>>> print(mol[sr.segid("MAIN", 0)])
Segment( MAIN num_atoms=66 )
>>> print(mol[sr.segid(name="MAIN", idx=0)])
Segment( MAIN num_atoms=66 )

.. note::

    Sire will return the Segment from an index operator if a segment
    identifying type is used as the index. This is slightly different
    behaviour to how the search string operates. In practice though,
    all molecular container classes behave in the same way, so you will
    often not notice or need to know which molecular container class
    has been returned.

Iterating over segments
-----------------------

The :class:`~sire.mol.Selector_Segment_` class is iterable, meaning that
it can be used in loops.

>>> for segment in mol.segments():
...     print(segment)
Segment( MAIN num_atoms=66 )

This is particularly helpful when combined with loops over the atoms in
a segment.

>>> for segment in mol.segments():
...    for atom in segment.atoms("element carbon"):
...        print(segment, atom.residue(), atom)
Segment( MAIN num_atoms=66 ) Residue( ACE:1   num_atoms=3 ) Atom( CA:1 )
Segment( MAIN num_atoms=66 ) Residue( ACE:1   num_atoms=3 ) Atom( C:2 )
Segment( MAIN num_atoms=66 ) Residue( ALA:2   num_atoms=6 ) Atom( CA:6 )
Segment( MAIN num_atoms=66 ) Residue( ALA:2   num_atoms=6 ) Atom( CB:7 )
Segment( MAIN num_atoms=66 ) Residue( ALA:2   num_atoms=6 ) Atom( C:8 )
Segment( MAIN num_atoms=66 ) Residue( ALA:3   num_atoms=6 ) Atom( CA:12 )
Segment( MAIN num_atoms=66 ) Residue( ALA:3   num_atoms=6 ) Atom( CB:13 )
...
Segment( MAIN num_atoms=66 ) Residue( ALA:11  num_atoms=6 ) Atom( C:62 )
Segment( MAIN num_atoms=66 ) Residue( CBX:12  num_atoms=3 ) Atom( CA:66 )

Finding all segment names
-------------------------

You can find the names of all segments using the :class:`~sire.mol.Select_Segment_.names`
function.

>>> print(mol.segments().names())
[SegName('MAIN')]
