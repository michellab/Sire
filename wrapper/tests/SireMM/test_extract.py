
from Sire.Mol import *
from Sire.MM import *
from Sire.IO import *

from nose.tools import assert_almost_equal

(molecules, space) = Amber().readCrdTop("../io/proteinbox.crd", "../io/proteinbox.top")

molnums = molecules.molNums()

for molnum in molnums:
    protein = molecules[molnum].molecule()

    if protein.nAtoms() > 200:
        break

residues = protein.selectAll( ResIdx(1) )

for i in range(2,11):
    residues = residues.add( ResIdx(i) )

partial = PartialMolecule(residues)
extract = partial.extract()

def test_internal_different():
    internal1 = InternalFF("1")
    internal2 = InternalFF("2")
    internal1.add(protein)
    internal2.add(extract)
    assert( internal1.energy() != internal2.energy() )

def test_internal_same():
    internal1 = InternalFF("1")
    internal2 = InternalFF("2")
    internal1.add(partial)
    internal1.setStrict(True)
    internal2.add(extract)
    assert_almost_equal( internal1.energy().value(), internal2.energy().value(), 3 )

def test_intra_different():
    intra1 = IntraCLJFF("1")
    intra2 = IntraCLJFF("2")
    intra1.add(protein)
    intra2.add(extract)
    assert( intra1.energy() != intra2.energy() )

def test_intra_same():
    intra1 = IntraCLJFF("1")
    intra2 = IntraCLJFF("2")
    intra1.add(partial)
    intra2.add(extract)
    assert_almost_equal( intra1.energy().value(), intra2.energy().value(), 3 )
