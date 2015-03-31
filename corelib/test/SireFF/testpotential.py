
from Sire.FF import *
from Sire.MM import *
from Sire.Mol import *
from Sire.IO import *
from Sire.Vol import *
from Sire.Units import *
from Sire.Maths import *

na = Molecule("sodium").edit() \
                       .add(CGName("1")) \
                       .add(AtomName("sodium")).commit()
cl = Molecule("chloride").edit() \
                         .add(CGName("1")) \
                         .add(AtomName("chloride")).commit()

na = na.edit() \
       .setProperty("element", Element("Na")) \
       .setProperty("coordinates", Vector(0,0,0)) \
       .setProperty("charge", 1*mod_electron) \
       .setProperty("LJ", LJParameter.dummy()) \
       .molecule() \
       .commit()

cl = cl.edit() \
       .setProperty("element", Element("Cl")) \
       .setProperty("coordinates", Vector(1,0,0)) \
       .setProperty("charge", -1*mod_electron) \
       .setProperty("LJ", LJParameter.dummy()) \
       .molecule() \
       .commit()

cljff = InterCLJFF()
cljff.add(na)
cljff.add(cl)

print(cljff.energy())

grid = RegularGrid(na.evaluate().center(), 50, 0.1*angstrom)
print(grid)

naff = InterCLJFF()
naff.add(na)
potentials = PotentialTable(naff[MGIdx(0)], grid)
naff.potential(potentials)

Cube().write(potentials, naff[MGIdx(0)], "test.cube")
