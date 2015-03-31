
from Sire.Mol import *
from Sire.IO import *

# The Sire.Units library contains definitions of all of the physical
# units, e.g. lengths, times, energies, etc.
from Sire.Units import *

# The Sire.MM library contains definitions of molecular mechanics
# parameters and molecular mechanics forcefields
from Sire.MM import *

# The Sire.Stream library contains functions used to save (stream)
# Sire objects to and from a disk
import Sire.Stream

water_dimer = PDB().read("input/water_dimer.pdb")

first_water = water_dimer[ MolIdx(0) ]
second_water = water_dimer[ MolIdx(1) ]

# Lets give TIP4P LJ parameters to the oxygen of each water molecule.
# We can do that by using the ".edit()" function, which allows
# you to edit the atom
oxygen = first_water.atom( AtomName("O00") )
oxygen = oxygen.edit().setProperty("LJ", LJParameter( 3.15363*angstrom,  \
                                                      0.1550*kcal_per_mol) )

# "oxygen" is a copy of the oxygen atom in "first_water". To copy this
# edit back to "first_water", we need to return to the molecule and
# "commit" the edit;
first_water = oxygen.molecule().commit()

# This has added a new "LJ" property to the list of properties of the molecule
print("Available properties of the first water molecule:")
print(first_water.propertyKeys())
print("The LJ property of the \"O00\" atom of the first water molecule:")
print(first_water.atom( AtomName("O00") ).property("LJ"))

# Do the same to the second molecule - in this case, we put everything
# together onto a single line
second_water = second_water.atom( AtomName("O00") ) \
                           .edit().setProperty("LJ", LJParameter( 3.15363*angstrom, \
                                                                  0.1550*kcal_per_mol) ) \
                           .molecule().commit()

# Now we use the same process to add TIP4P charges to the atoms...
first_water = first_water.edit() \
                         .atom( AtomName("M03") ) \
                         .setProperty("charge", -1.04*mod_electron) \
                         .molecule().atom( AtomName("H01") ) \
                         .setProperty("charge", 0.52*mod_electron) \
                         .molecule().atom( AtomName("H02") ) \
                         .setProperty("charge", 0.52*mod_electron) \
                         .molecule().commit()

# We can copy properties from one molecule to another. Here we copy the
# charges from "first_water" to "second_water"
second_water = second_water.edit() \
                           .setProperty("charge", first_water.property("charge") ) \
                           .commit()

# Note that the molecules must have the same number of atoms and residues etc.
# to allow you to copy properties from one to another!

print("\nParameters of the atoms in the second water molecule:")
for i in range(0, second_water.nAtoms()):
    atom = second_water.atom( AtomIdx(i) )
    print("%s : %s, %s" % (atom.name(), atom.property("charge"), atom.property("LJ")))

# Save the parameterised water molecules to the disk
print("\nSaving the parameterised water molecules to disk...")
Sire.Stream.save( (first_water,second_water), "water_dimer.s3" )

