
from Sire.IO import *
from Sire.Mol import *
from Sire.Maths import *
from Sire.Units import *

ethane = PDB().readMolecule("test/io/ethane.pdb")

C01 = AtomName("C01")
C05 = AtomName("C05")

ethane = ethane.edit().setProperty("connectivity", Connectivity(ethane)).commit()

# stretch the C-C bond moving each half by the same amount
for i in range(6, 12):
    r = i / 10.0

    ethane = ethane.move().set( BondID(C01, C05), r*angstrom ).commit()

    PDB().write(ethane, "test_ethane_%003d.pdb" % i)

# now put an anchor on the H02 atom
anchors = ethane.select(AtomName("H02")).selection()
ethane = ethane.edit().setProperty("anchors", anchors).commit()

for i in range(6, 12):
    r = i / 10.0

    ethane = ethane.move().set( BondID(C01, C05), r*angstrom,
                                {"anchors":"anchors"} ).commit()

    PDB().write(ethane, "test_anchor_ethane_%003d.pdb" % i)

