
from Sire.IO import *
from Sire.MM import *
from Sire.Maths import *
from Sire.Mol import *

(mols,space) = Amber().readCrdTop("../io/waterbox.crd", "../io/waterbox.top")

cljboxes = CLJBoxes()
idxs = []

natoms = 0

for i in range(0,mols.nMolecules()):
    idxs.append( cljboxes.add( CLJAtoms(mols[MolIdx(i)]) ) )
    natoms += mols[MolIdx(i)].molecule().nAtoms()

def test_cljdelta(verbose = False):
    old_water = mols[MolIdx(0)].molecule()
    old_cljatoms = CLJAtoms(old_water)
    test_cljatoms = cljboxes.atoms(idxs[0])

    if verbose:
        print("\nMaking sure the box has the right number of atoms %s vs. %s",
               cljboxes.nAtoms(), natoms)

    assert( cljboxes.nAtoms() == natoms )

    if verbose:
        print("\nChecking I can get the old atoms back out again...")
        print("OLD:\n%s" % old_cljatoms)
        print("TEST:\n%s" % test_cljatoms)

    assert(old_cljatoms == test_cljatoms)

    new_water = old_water.move().translate( Vector(1) ).commit()
    new_cljatoms = CLJAtoms(new_water)

    cljdelta = CLJDelta(1, cljboxes, idxs[0], new_water)

    test_cljatoms = cljdelta.newAtoms()

    if verbose:
        print("\nChecking changed atoms looks correct")
        print("CHANGED:\n%s" % cljdelta.changedAtoms())
        print("BOX: %s (%s,%s,%s) : %s" % (cljdelta.boxIndex(), \
                 cljdelta.nBoxX(), cljdelta.nBoxY(), cljdelta.nBoxZ(), \
                 cljdelta.isSingleBox()))

    if cljdelta.isSingleBox():
        assert( cljdelta.nBoxX() == 1 )
        assert( cljdelta.nBoxY() == 1 )
        assert( cljdelta.nBoxZ() == 1 )
    else:
        assert( cljdelta.nBoxX() > 1 or cljdelta.nBoxY() > 1 or cljdelta.nBoxZ() > 1 )
    
    if verbose:
        print("\nComparing new atoms are correctly in the delta")
        print("NEW:\n%s" % new_cljatoms)
        print("TEST:\n%s" % test_cljatoms)

    assert(new_cljatoms == test_cljatoms)

    test_cljatoms = cljdelta.oldAtoms()

    if verbose:
        print("\nComparing old atoms are correctly in the delta")
        print("OLD:\n%s" % old_cljatoms)
        print("TEST:\n%s" % test_cljatoms)

    assert(old_cljatoms == test_cljatoms)

    if verbose:
        print("\nTesting that the old indicies are correctly stored in the delta")
        print("OLD:\n%s" % idxs[0])
        print("NEW:\n%s" % cljdelta.oldIndicies())

    assert( idxs[0] == cljdelta.oldIndicies() )

    # now apply the delta on a copy of the boxes
    new_boxes = CLJBoxes(cljboxes)

    new_idxs = new_boxes.apply(cljdelta)

    if verbose:
        print("\nChecking...")
        print(new_idxs)

    assert( new_idxs == idxs[0] )

    test_atoms = new_boxes.atoms(new_idxs)

    if verbose:
        print("\nSeeing if the new atoms are in the box")
        print("NEW:\n%s" % new_cljatoms)
        print("TEST:\n%s" % test_atoms)

    assert(new_cljatoms == test_atoms)

    nold = cljboxes.nAtoms()
    nnew = new_boxes.nAtoms()

    assert( nold == nnew )

if __name__ == "__main__":
    test_cljdelta(True)

