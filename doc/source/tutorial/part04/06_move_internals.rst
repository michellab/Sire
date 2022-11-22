==================================
Moving Bonds, Angles and Dihedrals
==================================

Moving internals (e.g. bonds, angles and dihedrals) within molecules
is also very simple. You do this via a cursor that edits the internal
(or internals).

Moving Bonds
============

For example, here we will get a cursor that edits all of the
carbon-hydrogen bonds in the first molecule.

>>> mols = sr.load(sr.expand(sr.tutorial_url, ["ala.top", "ala.crd"]))
>>> cursor = mols[0].bonds("element C", "element H").cursor()

Cursors that edit bonds have two additional key functions;

* ``set_length`` / ``set_lengths`` and
* ``change_length`` / ``change_lengths``

These functions either set the bond length to a specified value,
or change the bond length by the specified delta.

For example, we can set the lengths of all of the carbon-hydrogen
bonds to 1 Å using ``set_length``

>>> cursor.set_length(1.0)
>>> print(cursor.lengths())
[1 Å, 1 Å, 1 Å, 1 Å, 1 Å, 1 Å, 1 Å, 1 Å, 1 Å, 1 Å]

.. note::

    Note that you can specify units of length if you want, e.g.
    ``set_length(1 * sr.units.angstrom)``. The default length units
    are used if units aren't specified.

.. note::

    Note also that ``set_length`` sets the length of all bonds edited
    by the cursor to the passed length. Use ``set_lengths`` and pass
    in a list of lengths if you want to set the bonds to
    different lengths.

We can change the length of the third carbon-hydrogen bond by 0.5 Å using
``change_length``, e.g.

>>> cursor[2].change_length(0.5)
>>> print(cursor[2].length())
1.5 Å
>>> print(cursor.lengths())
[1 Å, 1 Å, 1.5 Å, 1 Å, 1 Å, 1 Å, 1 Å, 1 Å, 1 Å, 1 Å]

Tracking movement with trajectories
===================================

You can visualise and analyse the movements you are performing by
saving trajectory frames. For example, lets gradually stretch the
carbon-carbon bond in the first molecule, saving trajectory
frames as we go.

>>> mol = mols[0]
>>> cursor = mol.cursor()
>>> cursor.save_frame()
>>> bond_cursor = cursor.bond("atomname CA", "resname ALA and atomname C")
>>> for i in range(0, 10):
...     bond_cursor.change_length(0.05)
...     bond_cursor.save_frame()
>>> mol = cursor.commit()

You can view a movie of this movement using

>>> mol.view()

.. image:: images/04_06_01.jpg
   :alt: Still from a movie of the aladip molecule being stretched.

You can get the energy for each frame of the trajectory using

>>> mol.trajectory().energy()
    frame  time    1-4_LJ  1-4_coulomb     angle       bond  dihedral  improper  intra_LJ  intra_coulomb       total
0       0   0.0  2.920343    44.880519  6.791894   4.548210  9.532594  0.285078 -1.102960     -46.370186   21.485492
1       1   0.0  2.607594    43.652433  6.791894   6.645588  9.532594  0.285078 -1.076241     -45.394975   23.043966
2       2   0.0  2.366888    42.448164  6.791894  10.327967  9.532594  0.285078 -1.046408     -44.442886   26.263291
3       3   0.0  2.182746    41.267230  6.791894  15.595345  9.532594  0.285078 -1.014466     -43.513222   31.127199
4       4   0.0  2.043079    40.109161  6.791894  22.447723  9.532594  0.285078 -0.981224     -42.605298   37.623007
5       5   0.0  1.938407    38.973497  6.791894  30.885102  9.532594  0.285078 -0.947316     -41.718473   45.740783
6       6   0.0  1.861278    37.859784  6.791894  40.907480  9.532594  0.285078 -0.913235     -40.852145   55.472728
7       7   0.0  1.805813    36.767577  6.791894  52.514858  9.532594  0.285078 -0.879364     -40.005697   66.812753
8       8   0.0  1.767366    35.696436  6.791894  65.707237  9.532594  0.285078 -0.845992     -39.178556   79.756056
9       9   0.0  1.742254    34.645928  6.791894  80.484615  9.532594  0.285078 -0.813341     -38.370160   94.298863
10     10   0.0  1.727556    33.615627  6.791894  96.846993  9.532594  0.285078 -0.781567     -37.579977  110.438198

Moving Angles
=============


Moving Dihedrals
================

Moving all dihedrals


Moving individual dihedrals


Aligning, Anchoring and Weighting
=================================

