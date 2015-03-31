
from Sire.Mol import *
from Sire.Maths import *
from Sire.IO import *
from Sire.Units import *
from Sire.Base import *

mol = Molecule("linear")

mol = mol.edit().add( CGName("1") ) \
                .add( AtomName("H1") ).cutGroup() \
                .add( AtomName("C1") ).cutGroup() \
                .add( AtomName("C2") ).cutGroup() \
                .add( AtomName("H2") ).cutGroup() \
                .molecule().commit()

print(mol, mol.nAtoms())

mol = mol.edit() \
         .atom(AtomName("H1")).setProperty("coordinates", Vector(0,0,0)) \
                              .setProperty("element", Element("H")).molecule() \
         .atom(AtomName("C1")).setProperty("coordinates", Vector(0,0,1)) \
                              .setProperty("element", Element("C")).molecule() \
         .atom(AtomName("C2")).setProperty("coordinates", Vector(0,0,2)) \
                              .setProperty("element", Element("C")).molecule() \
         .atom(AtomName("H2")).setProperty("coordinates", Vector(0,0,3)) \
                              .setProperty("element", Element("H")).molecule() \
         .commit()

connectivity = Connectivity(mol)
connectivity = connectivity.edit() \
                           .disconnectAll() \
                           .connect(AtomName("H1"), AtomName("C1")) \
                           .connect(AtomName("C1"), AtomName("C2")) \
                           .connect(AtomName("C2"), AtomName("H2")) \
                           .commit()

mol = mol.edit().setProperty("connectivity", connectivity).commit()

PDB().write(mol, "test0000.pdb")

angle = AngleID(AtomName("H1"),AtomName("C1"),AtomName("C2"))
dihedral = DihedralID(AtomName("H1"),AtomName("C1"),AtomName("C2"),AtomName("H2"))

mol = mol.move().change(angle, 30*degrees).commit()
PDB().write(mol, "test0001.pdb")

# This second angle move moves it back again. The direction for
# positive and negative angle changes is opposite for linear molecules,
# so +30 is actually equal to -30...
mol = mol.move().change(angle, 30*degrees).commit()
PDB().write(mol, "test0002.pdb")

mol = mol.move().change(dihedral, 45*degrees).commit()
PDB().write(mol, "test0003.pdb")
mol = mol.move().change(dihedral, -45*degrees).commit()
PDB().write(mol, "test0004.pdb")

mol = mol.move().change(angle, 30*degrees).commit()
PDB().write(mol, "test0005.pdb")

mol = mol.move().change(dihedral, 45*degrees).commit()
PDB().write(mol, "test0006.pdb")
mol = mol.move().change(dihedral, -45*degrees).commit()
PDB().write(mol, "test0007.pdb")

mol = mol.move().change(angle, 30*degrees).commit()
PDB().write(mol, "test0008.pdb")
