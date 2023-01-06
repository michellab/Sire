==================
Measuring Energies
==================

Sire implements a complete molecular mechanics (MM) engine, so you can
easily calculate MM energies for any parameterised molecules.

For example, load up the ``aladip`` system again...

>>> import sire as sr
>>> mols = sr.load(sr.expand(sr.tutorial_url, ["ala.top", "ala.crd"]))
>>> mol = mols[0]

You can get the energy of this molecule using the
:func:`~sire.mol.MoleculeView.energy` function;

>>> print(mol.energy())
21.4855 kcal mol-1

The energy is returned in sire's default units of energy, kilocalories
per mole. You can change this using, e.g.

>>> sr.units.set_energy_unit(sr.units.joule)
>>> print(mol.energy())
89895.3 joule mol-1

or can switch to sire's default SI units with

>>> sr.units.set_si_units()
>>> print(mol.energy())
89.8953 kJ mol-1

or can switch back to default internal units using

>>> sr.units.set_internal_units()
>>> print(mol.energy())
21.4855 kcal mol-1

Getting energy components
=========================

The returned energy is the total of the various MM component energies, e.g.
the bond energy, angle energy etc.

You can get all of these components using

>>> nrg = mol.energy()
>>> print(nrg.components())
{'intra_LJ': -1.10296 kcal mol-1, 'improper': 0.285078 kcal mol-1,
 'angle': 6.79189 kcal mol-1, 'intra_coulomb': -46.3702 kcal mol-1,
 'bond': 4.54821 kcal mol-1, '1-4_coulomb': 44.8805 kcal mol-1,
 '1-4_LJ': 2.92034 kcal mol-1, 'dihedral': 9.53259 kcal mol-1}

The components are indexed by their name, and can be returned directly by
treating the energy as a dictionary;

>>> print(nrg["bond"])
4.54821 kcal mol-1
>>> print(nrg["angle"])
6.79189 kcal mol-1

Note that these are the same values obtained via the :func:`~sire.mm.Bond.energy`
function on the Bond and Angle classes, e.g.

>>> print(mol.bonds().energy())
4.54821 kcal mol-1
>>> print(mol.angles().energy())
6.79189 kcal mol-1

Calculating energies of selections
==================================

You can calculate the energy of any molecular container. So, to calculate the
energy of all of the water molecules you could use

>>> print(mols["water"].energy())
-5907.76 kcal mol-1

or the energy of all carbon atoms

>>> print(mols["element C"].energy())
11.5485 kcal mol-1

or the components of the energy of all carbon atoms,

>>> print(mols["element C"].energy().components())
{'intra_LJ': -0.186179 kcal mol-1, 'angle': 0.0143363 kcal mol-1,
 'intra_coulomb': -4.6869 kcal mol-1, 'bond': 0.905869 kcal mol-1,
 '1-4_coulomb': 15.5136 kcal mol-1, '1-4_LJ': -0.012299 kcal mol-1}

Calculating interaction energies between selections
===================================================

You can also calculate interaction energies between molecular containers.
For example, to calculate the interaction energy between the first molecule
and all water molecules you would use

>>> print(mols[0].energy(mols["water"]))
-41.4064 kcal mol-1

or to calculate the interaction energy between the first two residues
of the first molecule use

>>> print(mols[0].residues()[0].energy(mols[0].residues()[1]))
-21.3407 kcal mol-1

Again, you can use the ``components()`` function to decompose this
into individual energy components.

>>> print(mols[0].residues()[0].energy(mols[0].residues()[1]).components())
{'intra_LJ': -0.59062 kcal mol-1, 'improper': 0.106471 kcal mol-1,
 'angle': 1.12332 kcal mol-1, 'intra_coulomb': -7.90428 kcal mol-1,
 'bond': 0.00857782 kcal mol-1, '1-4_coulomb': -19.1969 kcal mol-1,
 '1-4_LJ': 0.877224 kcal mol-1, 'dihedral': 4.23544 kcal mol-1}

Again, the values are the same as you would have got calling ``energy`` on
the corresponding bonds, angles, dihedrals etc, e.g.

>>> print(mols[0].bonds("residx 0", "residx 1").energy())
0.00857782 kcal mol-1

Decomposing into individual energies
====================================

The result of calling ``mol.energy()`` and ``mol.atoms().energy()``
is the same,

>>> print(mol.energy())
21.4855 kcal mol-1
>>> print(mol.atoms().energy())
21.4855 kcal mol-1

This is because the ``.energy()`` function returns the sum of the
energies of all views within its molecular container. The total
energy of all of the atoms in a molecule must be equal to the total
energy of the molecule.

Often, you want to see the individual energies of the views. You
may think that you could do this just be looping over the views
in the container, e.g.

>>> for atom in mol.atoms():
...     print(atom.energy())
0
0
0
0
0
..
0
0
0
0
0

but, as you can see above, this is not the case. This is because the
``energy()`` function returns the energy of the view alone, i.e. not
including interactions with any other views. An individual atom has
no energy on its own, hence why we got the zero values above.

We can demonstrate this further by looking at decomposing a molecule's
energy into residue-based components. The total molecular energy...

>>> print(mol.energy())
21.4855 kcal mol-1

is equal to the sum of the energies of its three constituent residues...

>>> total = 0 * sr.units.kcal_per_mol
>>> for residue in mol.residues():
...     print(residue.energy())
...     total += residue.energy()
-13.8435 kcal mol-1
32.8847 kcal mol-1
12.8614 kcal mol-1

...plus the energy of interaction between each pair of residues...

>>> for i in range(0, 2):
...    for j in range(i+1, 3):
...        print(i, j, mol.residues()[i].energy(mol.residues()[j]))
...        total += mol.residues()[i].energy(mol.residues()[j])
0 1 -21.3407 kcal mol-1
0 2 -0.709678 kcal mol-1
1 2 11.6333 kcal mol-1
>>> print(total)
21.4855 kcal mol-1

Decomposing into interaction energies
=====================================

It is much easier to decompose interaction energies. For example,
the interaction energy between the first molecule and all water
molecules is equal to the sum of the interaction energy between
the first molecule's atoms and the water molecules.

>>> print(mol.energy(mols["water"]))
-41.4064 kcal mol-1
>>> total = 0 * sr.units.kcal_per_mol
>>> for atom in mol.atoms():
...     print(atom.energy(mols["water"]))
...     total += atom.energy(mols["water"])
0.609576 kcal mol-1
-1.72327 kcal mol-1
-0.0995291 kcal mol-1
-0.863389 kcal mol-1
3.62748 kcal mol-1
-16.1531 kcal mol-1
-1.379 kcal mol-1
-2.86089 kcal mol-1
-1.39983 kcal mol-1
-0.0325111 kcal mol-1
-4.24037 kcal mol-1
0.112378 kcal mol-1
0.714725 kcal mol-1
0.64597 kcal mol-1
3.88331 kcal mol-1
-16.7354 kcal mol-1
-2.67356 kcal mol-1
-1.99815 kcal mol-1
-2.38936 kcal mol-1
0.477344 kcal mol-1
0.120511 kcal mol-1
0.950679 kcal mol-1
>>> print(total)
-41.4064 kcal mol-1

As well as using a loop, you could use the ``apply()`` function
to call ``energy`` on each view in a container, e.g.

>>> print(mol.apply("energy", mols["water"]))
[0.609576 kcal mol-1, -1.72327 kcal mol-1, -0.0995291 kcal mol-1,
-0.863389 kcal mol-1, 3.62748 kcal mol-1, -16.1531 kcal mol-1,
-1.379 kcal mol-1, -2.86089 kcal mol-1, -1.39983 kcal mol-1,
-0.0325111 kcal mol-1, -4.24037 kcal mol-1, 0.112378 kcal mol-1,
 0.714725 kcal mol-1, 0.64597 kcal mol-1, 3.88331 kcal mol-1,
-16.7354 kcal mol-1, -2.67356 kcal mol-1, -1.99815 kcal mol-1,
-2.38936 kcal mol-1, 0.477344 kcal mol-1, 0.120511 kcal mol-1,
 0.950679 kcal mol-1]

and can calculate the sum automatically using ``apply_reduce()``, e.g.

>>> print(mol.apply_reduce(lambda atom: atom.energy(mols["water"])))
-41.4064 kcal mol-1

Because this is such a common thing that you may want to do, sire provides
the ``.energies()`` function that does this automatically, e.g.

>>> print(mol.atoms().energies(mols["water"]))
[0.609576 kcal mol-1, -1.72327 kcal mol-1, -0.0995291 kcal mol-1,
-0.863389 kcal mol-1, 3.62748 kcal mol-1, -16.1531 kcal mol-1,
-1.379 kcal mol-1, -2.86089 kcal mol-1, -1.39983 kcal mol-1,
-0.0325111 kcal mol-1, -4.24037 kcal mol-1, 0.112378 kcal mol-1,
 0.714725 kcal mol-1, 0.64597 kcal mol-1, 3.88331 kcal mol-1,
-16.7354 kcal mol-1, -2.67356 kcal mol-1, -1.99815 kcal mol-1,
-2.38936 kcal mol-1, 0.477344 kcal mol-1, 0.120511 kcal mol-1,
 0.950679 kcal mol-1]

calculates the energy of each atom in the solute with each water
molecule, while

>>> print(mols[1:].energies(mols[0]))
[-0.0489166 kcal mol-1, -0.0414449 kcal mol-1, -0.0270237 kcal mol-1,
  0.00493017 kcal mol-1, -0.15308 kcal mol-1, -0.805673 kcal mol-1,
...
  0.113608 kcal mol-1, 0.147355 kcal mol-1, -0.0893536 kcal mol-1,
 -0.0101503 kcal mol-1, -0.432695 kcal mol-1, -0.0120518 kcal mol-1]

calculates the energies between the first molecule and every other
molecule in the system.
