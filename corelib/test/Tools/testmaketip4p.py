
from Sire.Tools import AmberLoader
from Sire.Tools import Parameter

from Sire.IO import *

params = {}
params["water model"] = "tip4p"

Parameter.push(params)

watersys = AmberLoader.createSystem("test/io/waterbox.top", "test/io/waterbox.crd")

PDB().write(watersys.molecules(), "test_tip4p.pdb")
