==================
Indexing Molecules
==================

Molecules are collections of atoms that (at least notionally) should all
be bonded to each other and represent the concept of a molecule in a
chemical system.

A molecule is a molecular container for atoms, residues, chains and
segments, and is implemented via the :class:`sire.mol.Molecule` class.
At a minimum, a :class:`sire.mol.Molecule` contains at least one
:class:`sire.mol.Atom`. There
is no requirement for this atom to be placed into a residue, chain
or segment, and so it is possible for a molecule to have zero
residues, chains or segments.

There are several classes in ``sire`` that can act as molecular containers
for molecule objects. You have already encountered one, which is the
:class:`~sire.system.System` class, into which the molecules were first
loaded.

You can get a