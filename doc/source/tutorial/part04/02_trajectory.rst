======================
Measuring Trajectories
======================

Sire has in-built support for trajectories. This allows you to measure
bonds, angles, dihedrals and other properties for across a range of
frames of a trajectory. This can be useful to see how these measurements
changed during, e.g. a molecular dynamics simulation.

For example, let's load a trajectory contained in a DCD file.

>>> import sire as sr
>>> mols = sr.load(sr.expand(sr.tutorial_url, ["h7n9.pdb", "h7n9.dcd"]))
>>> mol = mols[0]
>>> print(mol.num_frames())
501

You load trajectories in the same way as you load any other molecular file
in Sire. Just add the trajectory file (here ``h7n9.dcd``) to the list of
files to load. The file type will be determined automatically, and the
trajectory frames loaded in order (e.g. if you pass in multiple trajectory
files to load).

Sire recognises two main types of molecular data;

1. Topology data - this is information about the structure or topology
   of the molecules, e.g. how the atoms divide into residues, and how
   these divide across molecules.

2. Frame data - this is information that is specific to an individual
   trajectory frame, e.g. coordinates, velocities, space and time
   information.

Individual molecular input files contain topology data and/or frame data.
In this case, ``h7n9.pdb`` contains both topology and frame data, while
``h7n9.dcd`` contains just frame data.

Sire will load the topology information from the first file in the list
of files that contains topology information (e.g. ``h7n9.pdb``). It will
then load frames in order from all files that contain frame data
(so both ``h7n9.pdb`` and ``h7n9.dcd``).

In this case, 501 frames have been loaded. These are the one frame from
``h7n9.pdb``, which is the first frame. And then the next 500 frames
are the five hundred loaded from ``h7n9.dcd``.

In many cases you would only want to analyse the trajectory data contained
in the DCD file, and would not be interested in the frame data from the
PDB file. We can remove this frame data using the ``delete_frame`` function,
e.g.

>>> mol.delete_frame(0)
>>> print(mol.num_frames())
500

Iterating over trajectories
===========================

You iterate over the frames in a trajectory using a
:class:`sire.mol.TrajectoryIterator` that is returned by the
``.trajectory()`` function.

For example, here we will iterate over all 500 frames and just print
out the first atom of the molecule.

>>> for frame in mol.trajectory():
...     print(frame[0])
Atom( N:1     [  15.62,   42.18,   15.18] )
Atom( N:1     [  15.24,   43.24,   15.65] )
Atom( N:1     [  15.05,   42.39,   15.22] )
Atom( N:1     [  16.16,   41.74,   14.25] )
Atom( N:1     [  16.62,   41.16,   15.16] )
Atom( N:1     [  14.37,   42.47,   15.84] )
Atom( N:1     [  16.63,   42.29,   13.91] )
... etc.

You can call ``.trajectory()`` to iterate over any view of a molecule
(or views of molecules). For example, we could loop over the first
atom directly,

>>> for atom_frame in mol[0].trajectory():
...     print(atom_frame)
Atom( N:1     [  15.62,   42.18,   15.18] )
Atom( N:1     [  15.24,   43.24,   15.65] )
Atom( N:1     [  15.05,   42.39,   15.22] )
Atom( N:1     [  16.16,   41.74,   14.25] )
Atom( N:1     [  16.62,   41.16,   15.16] )
Atom( N:1     [  14.37,   42.47,   15.84] )
Atom( N:1     [  16.63,   42.29,   13.91] )
... etc.

or we could loop over the frames of the first residue, print out
its center of mass coordinates;

>>> for res_frame in mol.residue(0).trajectory():
...     print(res_frame, res_frame.coordinates())
Residue( ARG:1   num_atoms=26 ) ( 14.7759 Å, 42.2632 Å, 18.201 Å )
Residue( ARG:1   num_atoms=26 ) ( 14.368 Å, 43.886 Å, 18.7757 Å )
Residue( ARG:1   num_atoms=26 ) ( 14.3506 Å, 42.9594 Å, 18.359 Å )
Residue( ARG:1   num_atoms=26 ) ( 15.4561 Å, 41.4458 Å, 17.4479 Å )
Residue( ARG:1   num_atoms=26 ) ( 15.3717 Å, 41.1966 Å, 18.1695 Å )
Residue( ARG:1   num_atoms=26 ) ( 13.4568 Å, 43.847 Å, 18.778 Å )
Residue( ARG:1   num_atoms=26 ) ( 15.1335 Å, 43.0555 Å, 17.4046 Å )
etc.

The :class:`~sire.mol.TrajectoryIterator` is itself indexable. This allows
you to slice the frames that are iterated over. Here, we will iterate
over just the first 5 frames of the trajectory

>>> for res_frame in mol.residue(0).trajectory()[0:5]:
...     print(res_frame, res_frame.coordinates())
Residue( ARG:1   num_atoms=26 ) ( 14.7759 Å, 42.2632 Å, 18.201 Å )
Residue( ARG:1   num_atoms=26 ) ( 14.368 Å, 43.886 Å, 18.7757 Å )
Residue( ARG:1   num_atoms=26 ) ( 14.3506 Å, 42.9594 Å, 18.359 Å )
Residue( ARG:1   num_atoms=26 ) ( 15.4561 Å, 41.4458 Å, 17.4479 Å )
Residue( ARG:1   num_atoms=26 ) ( 15.3717 Å, 41.1966 Å, 18.1695 Å )

Measuring bond lengths over a trajectory
========================================

Bonds, Angles, Dihedrals and Impropers are also molecular views,
and so you can iterate over their trajectories in the same way.

For example, the ``H7N9`` protein contains a number of disulfide bonds.

We can find these bonds using

>>> print(mol.bonds("element S", "element S"))
SelectorBond( size=9
0: Bond( SG:177 => SG:5098 )
1: Bond( SG:672 => SG:735 )
2: Bond( SG:1479 => SG:1736 )
3: Bond( SG:1586 => SG:2317 )
4: Bond( SG:2343 => SG:2407 )
5: Bond( SG:3041 => SG:3217 )
6: Bond( SG:3062 => SG:3193 )
7: Bond( SG:3641 => SG:3905 )
8: Bond( SG:5163 => SG:5600 )
)

We can iterate over the trajectory frames for these bonds, and then
measure them. For example, lets iterate over the first ten frames of the
first ``S-S`` bond, and print it's length.

>>> for bond_frame in mol.bonds("element S", "element S")[0].trajectory()[0:10]:
...     print(bond_frame, bond_frame.measure())
Bond( SG:177 => SG:5098 ) 1.97428 Å
Bond( SG:177 => SG:5098 ) 2.12148 Å
Bond( SG:177 => SG:5098 ) 2.01916 Å
Bond( SG:177 => SG:5098 ) 2.03011 Å
Bond( SG:177 => SG:5098 ) 1.92681 Å
Bond( SG:177 => SG:5098 ) 2.05971 Å
Bond( SG:177 => SG:5098 ) 2.00386 Å
Bond( SG:177 => SG:5098 ) 1.9727 Å
Bond( SG:177 => SG:5098 ) 2.06847 Å
Bond( SG:177 => SG:5098 ) 1.98437 Å

We could have done this for all of the bonds, using...

>>> for bonds_frame in mol.bonds("element S", "element S").trajectory()[0:10]:
...     print(bonds_frame, bonds_frame.measures())
SelectorBond( size=9
0: Bond( SG:177 => SG:5098 )
1: Bond( SG:672 => SG:735 )
2: Bond( SG:1479 => SG:1736 )
3: Bond( SG:1586 => SG:2317 )
4: Bond( SG:2343 => SG:2407 )
5: Bond( SG:3041 => SG:3217 )
6: Bond( SG:3062 => SG:3193 )
7: Bond( SG:3641 => SG:3905 )
8: Bond( SG:5163 => SG:5600 )
) [1.97428 Å, 2.03112 Å, 2.07281 Å, 2.0191 Å, 2.04427 Å, 2.06217 Å, 2.06375 Å, 1.98086 Å, 2.00846 Å]
SelectorBond( size=9
0: Bond( SG:177 => SG:5098 )
1: Bond( SG:672 => SG:735 )
2: Bond( SG:1479 => SG:1736 )
3: Bond( SG:1586 => SG:2317 )
4: Bond( SG:2343 => SG:2407 )
5: Bond( SG:3041 => SG:3217 )
6: Bond( SG:3062 => SG:3193 )
7: Bond( SG:3641 => SG:3905 )
8: Bond( SG:5163 => SG:5600 )
) [2.12148 Å, 2.02085 Å, 2.01314 Å, 2.02394 Å, 2.03679 Å, 2.05127 Å, 2.13314 Å, 2.09479 Å, 2.01281 Å]
SelectorBond( size=9
0: Bond( SG:177 => SG:5098 )
1: Bond( SG:672 => SG:735 )
2: Bond( SG:1479 => SG:1736 )
3: Bond( SG:1586 => SG:2317 )
4: Bond( SG:2343 => SG:2407 )
5: Bond( SG:3041 => SG:3217 )
6: Bond( SG:3062 => SG:3193 )
7: Bond( SG:3641 => SG:3905 )
8: Bond( SG:5163 => SG:5600 )
) [2.01916 Å, 2.07407 Å, 2.13044 Å, 2.05 Å, 1.94306 Å, 2.02388 Å, 1.99157 Å, 2.0498 Å, 2.11982 Å]
etc...

...but you can see that we quickly reach the limit of what can sensibly
be printed to the screen.

To make things easier, the ``.measures()`` function on the
:class:`~sire.mol.TrajectoryIterator` will calculate all of the measures
across all of its frames, and will return the result as a pandas dataframe.

>>> df = mol.bonds("element S", "element S").trajectory().measures()
>>> print(df)


