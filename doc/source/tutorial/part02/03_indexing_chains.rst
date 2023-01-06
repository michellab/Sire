===============
Indexing Chains
===============

Chains are collections of residues. They typically represent a single
chain in a protein. Chains are implemented via the :class:`~sire.mol.Chain`
class, which itself is a molecular container for
:class:`~sire.mol.Residue` objects. A residue can only belong to one
chain at a time (and residues do not need to be assigned to chains).

You can access the chains in a molecule container using the
:func:`~sire.mol.Chain.chain` and :func:`~sire.mol.Chain.chains` functions,
which are available on all of the molecular container types.

>>> print(mol.chain(0))
Chain( A num_residues=123 num_atoms=985)

gives the chain at index 0, while

>>> print(mol.chain("B"))
Chain( B num_residues=638 num_atoms=4881)

gives the chain called "B", and

>>> print(mol.chains())
Selector<SireMol::Chain>( size=4
0:  Chain( A num_residues=123 num_atoms=985)
1:  Chain( B num_residues=638 num_atoms=4881)
2:  Chain( C num_residues=126 num_atoms=1000)
3:  Chain( D num_residues=631 num_atoms=4862)
)

returns all of the chains.

Like all other containers, you can slice the function call or result,
so both

>>> print(mol.chains(range(3,-1,-1)))
Selector<SireMol::Chain>( size=4
0:  Chain( D num_residues=631 num_atoms=4862)
1:  Chain( C num_residues=126 num_atoms=1000)
2:  Chain( B num_residues=638 num_atoms=4881)
3:  Chain( A num_residues=123 num_atoms=985)
)

>>> print(mol.chains()[3:-1:-1])
Selector<SireMol::Chain>( size=4
0:  Chain( D num_residues=631 num_atoms=4862)
1:  Chain( C num_residues=126 num_atoms=1000)
2:  Chain( B num_residues=638 num_atoms=4881)
3:  Chain( A num_residues=123 num_atoms=985)
)

slice to get the four chains in reverse order.

Search for chains
-----------------

You can also search for chains, using their name (``chainname``),
and/or their index in their parent molecule (``chainidx``).

>>> print(mol.chain("chainname A"))
Chain( A num_residues=123 num_atoms=985)

>>> print(mol.chain("chainidx 1"))
Chain( B num_residues=638 num_atoms=4881)

.. note::

   Unlike atoms and residues, chains do not have a number. They
   are identified only by their index in their parent molecule, or
   their name

You can use the chain search string in a molecular container's index
operator too!

>>> print(mol["chainname A"])
Chain( A num_residues=123 num_atoms=985)

and you can combine it with residue and/or atom identifiers, e.g.

>>> print(mol["chainname A and resname ALA"])
Selector<SireMol::Residue>( size=11
0:  Residue( ALA:23  num_atoms=5 )
1:  Residue( ALA:30  num_atoms=5 )
2:  Residue( ALA:53  num_atoms=5 )
3:  Residue( ALA:65  num_atoms=5 )
4:  Residue( ALA:85  num_atoms=5 )
...
6:  Residue( ALA:96  num_atoms=5 )
7:  Residue( ALA:104 num_atoms=5 )
8:  Residue( ALA:105 num_atoms=5 )
9:  Residue( ALA:122 num_atoms=5 )
10:  Residue( ALA:158 num_atoms=5 )
)

>>> print(mol["chainname A and element O"])
Selector<SireMol::Atom>( size=197
0:  Atom( O:4     [ -57.04,    9.73,   41.82] )
1:  Atom( O:12    [ -57.27,    6.71,   41.66] )
2:  Atom( O:19    [ -53.43,    4.61,   40.31] )
3:  Atom( O:27    [ -55.45,    2.35,   37.41] )
4:  Atom( O:36    [ -54.71,   -0.60,   37.34] )
...
192:  Atom( OE1:977 [ -30.32,   17.02,    0.48] )
193:  Atom( OE2:978 [ -31.39,   15.27,   -0.31] )
194:  Atom( O1:11662 [ -40.43,   -8.92,    3.99] )
195:  Atom( O2:11664 [ -40.18,   -8.25,    0.56] )
196:  Atom( O4:11667 [ -42.40,   -9.39,   -2.03] )
)

You can also search for multiple chain names

>>> print(mol["chainname A, B"])
Selector<SireMol::Chain>( size=2
0:  Chain( A num_residues=123 num_atoms=985)
1:  Chain( B num_residues=638 num_atoms=4881)
)

Wildcard (glob) searching is also supported for chain names

>>> print(mol["chainname /[cd]/i"])
Selector<SireMol::Chain>( size=2
0:  Chain( C num_residues=126 num_atoms=1000)
1:  Chain( D num_residues=631 num_atoms=4862)
)

Finding the residues in a chain
-------------------------------

Because both :class:`~sire.mol.Chain` and :class:`~sire.mol.Selector_Chain_`
are molecular containers, they also have their own
:func:`~sire.mol.Residue.residue` and :func:`~sire.mol.Residue.residues` functions,
which behave as you would expect.

>>> print(mol["chainname A"].residue(sr.resid("ALA", 23)))
Residue( ALA:23  num_atoms=5 )

You can get all of the residues in a chain by calling the
:func:`~sire.mol.Chain.residues` function without any arguments,

>>> print(mol["chainname A"].residues())
Selector<SireMol::Residue>( size=123
0:  Residue( ILE:6   num_atoms=8 )
1:  Residue( VAL:7   num_atoms=7 )
2:  Residue( LEU:8   num_atoms=8 )
3:  Residue( LYS:9   num_atoms=9 )
4:  Residue( SER:10  num_atoms=6 )
...
118:  Residue( TRP:157 num_atoms=14 )
119:  Residue( ALA:158 num_atoms=5 )
120:  Residue( PHE:159 num_atoms=11 )
121:  Residue( GLU:160 num_atoms=9 )
122:  Residue( PEG:801 num_atoms=7 )
)

In addition, the index operator for chains searches by default for residues,
not for atoms. Thus

>>> print(mol["chainname A"][0])
Residue( ILE:6   num_atoms=8 )

gives the first *residue* in chain "A", not the first atom. Similarly

>>> print(mol["chainname A"]["ALA"])
Selector<SireMol::Residue>( size=11
0:  Residue( ALA:23  num_atoms=5 )
1:  Residue( ALA:30  num_atoms=5 )
2:  Residue( ALA:53  num_atoms=5 )
3:  Residue( ALA:65  num_atoms=5 )
4:  Residue( ALA:85  num_atoms=5 )
...
6:  Residue( ALA:96  num_atoms=5 )
7:  Residue( ALA:104 num_atoms=5 )
8:  Residue( ALA:105 num_atoms=5 )
9:  Residue( ALA:122 num_atoms=5 )
10:  Residue( ALA:158 num_atoms=5 )
)

searches for the *residues* called "ALA", not the atoms called "ALA".
This is because chains are containers for residues, not containers for atoms.

Another route is to use ``residues in`` in the search string.

>>> print(mol["residues in chainname B"])
Selector<SireMol::Residue>( size=123
0:  Residue( ILE:6   num_atoms=8 )
1:  Residue( VAL:7   num_atoms=7 )
2:  Residue( LEU:8   num_atoms=8 )
3:  Residue( LYS:9   num_atoms=9 )
4:  Residue( SER:10  num_atoms=6 )
...
118:  Residue( TRP:157 num_atoms=14 )
119:  Residue( ALA:158 num_atoms=5 )
120:  Residue( PHE:159 num_atoms=11 )
121:  Residue( GLU:160 num_atoms=9 )
122:  Residue( PEG:801 num_atoms=7 )
)

or, go from residues to chains using ``chains with``

>>> print(mol["chains with resname ALA"])
Selector<SireMol::Chain>( size=4
0:  Chain( A num_residues=123 num_atoms=985)
1:  Chain( B num_residues=638 num_atoms=4881)
2:  Chain( C num_residues=126 num_atoms=1000)
3:  Chain( D num_residues=631 num_atoms=4862)
)

Finding the atoms in a chain
----------------------------

You can still get the atoms in a chain by calling the
:func:`~sire.mol.Chain.atom` and :func:`~sire.mol.Chain.atoms` functions.

>>> print(mol["chainidx 0"].atoms("CA"))
Selector<SireMol::Atom>( size=122
0:  Atom( CA:2    [ -55.43,   11.35,   42.54] )
1:  Atom( CA:10   [ -56.02,    7.64,   43.47] )
2:  Atom( CA:17   [ -54.99,    6.39,   39.98] )
3:  Atom( CA:25   [ -55.33,    2.58,   39.80] )
4:  Atom( CA:34   [ -52.97,    1.03,   37.19] )
...
117:  Atom( CA:932  [ -35.65,   13.78,   -7.21] )
118:  Atom( CA:941  [ -38.64,   16.18,   -7.49] )
119:  Atom( CA:955  [ -38.15,   18.05,   -4.20] )
120:  Atom( CA:960  [ -35.38,   20.38,   -3.07] )
121:  Atom( CA:971  [ -33.67,   19.61,    0.23] )
)

>>> print(mol["chainname B"].atom(0))
Atom( N:980   [ -31.52,  -13.85,   36.51] )

Calling the :func:`~sire.mol.Chain.atoms` function without any arguments
returns all of the atoms in the chain.

>>> print(mol["chainname C"].atoms())
Selector<SireMol::Atom>( size=1000
0:  Atom( N:5824  [ -29.76,   20.54,   63.08] )
1:  Atom( CA:5825 [ -30.98,   21.16,   63.60] )
2:  Atom( C:5826  [ -31.03,   22.64,   63.24] )
3:  Atom( O:5827  [ -30.95,   23.01,   62.07] )
4:  Atom( CB:5828 [ -31.07,   20.98,   65.13] )
...
995:  Atom( O2:11698 [   8.19,   41.99,   41.66] )
996:  Atom( C3:11699 [   7.27,   42.05,   42.71] )
997:  Atom( C4:11700 [   5.85,   42.13,   42.15] )
998:  Atom( O4:11701 [   4.93,   42.21,   43.20] )
999:  Atom( O:11726 [  -2.64,   40.59,   37.49] )
)

Another route is to use the ``atoms in`` phrase in the search string,

>>> print(mol["atoms in chainname B"])
Selector<SireMol::Atom>( size=4881
0:  Atom( N:980   [ -31.52,  -13.85,   36.51] )
1:  Atom( CA:981  [ -32.98,  -13.85,   36.21] )
2:  Atom( C:982   [ -33.24,  -13.01,   34.96] )
3:  Atom( O:983   [ -33.85,  -13.54,   34.01] )
4:  Atom( CB:984  [ -33.76,  -13.31,   37.41] )
...
4876:  Atom( O:11721 [ -33.17,   25.70,    2.16] )
4877:  Atom( O:11722 [ -44.37,   53.10,  -15.15] )
4878:  Atom( O:11723 [ -39.25,   59.56,   -5.85] )
4879:  Atom( O:11724 [ -52.78,   63.80,   12.56] )
4880:  Atom( O:11725 [ -36.82,   45.49,    3.24] )
)

and to use ``chains with`` to go from atoms to chains.

>>> print(mol["chains with atomname CA"])
Selector<SireMol::Chain>( size=4
0:  Chain( A num_residues=123 num_atoms=985)
1:  Chain( B num_residues=638 num_atoms=4881)
2:  Chain( C num_residues=126 num_atoms=1000)
3:  Chain( D num_residues=631 num_atoms=4862)
)

Uniquely identifying a chain
----------------------------

You uniquely identify a chain in a molecule using its chain index
(``chainidx``). You can get the index of a chain in a molecule by
calling its :func:`~sire.mol.Chain.index` function.

>>> print(mol.chain(0).index())
ChainIdx(0)

.. warning::

    Be careful indexing by chain index. This is the index of the chain
    that uniquely identifies it within its parent molecule. It is not the
    index of the chain in an arbitrary molecular container.

Chain identifying types
-----------------------

Another way to index chains is to use the chain identifying types, i.e.
:class:`~sire.mol.ChainName` and :class:`~sire.mol.ChainIdx`. The
easiest way to create these is by using the function
:func:`sire.chainid`.

Use strings to create :class:`~sire.mol.ChainName` objects,

>>> print(sr.chainid("A"))
ChainName('A')
>>> print(mol[sr.chainid("A")])
Chain( A num_residues=123 num_atoms=985)

and integers to create :class:`~sire.mol.ChainIdx` objects.

>>> print(sr.chainid(0))
ChainIdx(0)
>>> print(mol[sr.chainid(0)])
Chain( A num_residues=123 num_atoms=985)

You can set both a name and an index by passing in two arguments.

>>> print(mol[sr.chainid("A", 0)])
Chain( A num_residues=123 num_atoms=985)
>>> print(mol[sr.chainid(name="A", idx=0)])
Chain( A num_residues=123 num_atoms=985)

Iterating over chains
---------------------

The :class:`~sire.mol.Selector_Chain_` class is iterable, meaning that
it can be used in loops.

>>> for chain in mol.chains():
...     print(chain)
Chain( A num_residues=123 num_atoms=985)
Chain( B num_residues=638 num_atoms=4881)
Chain( C num_residues=126 num_atoms=1000)
Chain( D num_residues=631 num_atoms=4862)

This is particulary useful when combined with looping over the
residues and/or atoms in the residues.

>>> for chain in mol.chains():
...     for residue in chain.residues():
...         print(chain, residue)
Chain( A num_residues=123 num_atoms=985) Residue( ILE:6   num_atoms=8 )
Chain( A num_residues=123 num_atoms=985) Residue( VAL:7   num_atoms=7 )
Chain( A num_residues=123 num_atoms=985) Residue( LEU:8   num_atoms=8 )
Chain( A num_residues=123 num_atoms=985) Residue( LYS:9   num_atoms=9 )
Chain( A num_residues=123 num_atoms=985) Residue( SER:10  num_atoms=6 )
Chain( A num_residues=123 num_atoms=985) Residue( SER:11  num_atoms=6 )
Chain( A num_residues=123 num_atoms=985) Residue( ASP:12  num_atoms=8 )
Chain( A num_residues=123 num_atoms=985) Residue( GLY:13  num_atoms=4 )
...
Chain( D num_residues=631 num_atoms=4862) Residue( HOH:905 num_atoms=1 )
Chain( D num_residues=631 num_atoms=4862) Residue( HOH:906 num_atoms=1 )

>>> for chain in mol["chainname A, B"]:
...     for atom in chain["element C"]:
...         print(chain, atom.residue(), atom)
Chain( A num_residues=123 num_atoms=985) Residue( ILE:6   num_atoms=8 ) Atom( CA:2    [ -55.43,   11.35,   42.54] )
Chain( A num_residues=123 num_atoms=985) Residue( ILE:6   num_atoms=8 ) Atom( C:3     [ -56.06,    9.95,   42.55] )
Chain( A num_residues=123 num_atoms=985) Residue( ILE:6   num_atoms=8 ) Atom( CB:5    [ -56.32,   12.33,   41.76] )
Chain( A num_residues=123 num_atoms=985) Residue( ILE:6   num_atoms=8 ) Atom( CG1:6   [ -55.68,   13.72,   41.72] )
Chain( A num_residues=123 num_atoms=985) Residue( ILE:6   num_atoms=8 ) Atom( CG2:7   [ -57.70,   12.40,   42.39] )
Chain( A num_residues=123 num_atoms=985) Residue( ILE:6   num_atoms=8 ) Atom( CD1:8   [ -55.42,   14.31,   43.09] )
...
Chain( B num_residues=638 num_atoms=4881) Residue( CIT:803 num_atoms=13 ) Atom( C4:11688 [ -28.14,   10.72,   -2.13] )
Chain( B num_residues=638 num_atoms=4881) Residue( CIT:803 num_atoms=13 ) Atom( C5:11689 [ -28.94,   10.62,   -3.43] )
Chain( B num_residues=638 num_atoms=4881) Residue( CIT:803 num_atoms=13 ) Atom( C6:11692 [ -27.91,    9.78,    0.15] )

Finding all chain names
-----------------------

You can find the names of all chains using the :class:`~sire.mol.Select_Chain_.names`
function.

>>> print(mol.chains().names())
[ChainName('A'), ChainName('B'), ChainName('C'), ChainName('D')]
