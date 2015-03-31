
from Sire.Mol import *
from Sire.Maths import *
from Sire.IO import *

#load up the TIP4P solvent box
tip4pwaters = PDB().read("test/geometry/tip4pbox.pdb")

#change each TIP4P in the box into an SPC water
for water in tip4pwaters:
    #remove the M atom
    water[0].remove("M03")
    
    #set the H-O bond length to 1.0 (SPC bond length)
    water[0].setBond("H01","O00",1.0)
    water[0].setBond("H02","O00",1.0)
    
    #set the angle to 104.5 degrees (SPC H-O-H angle)
    water[0].setAngle("H01","O00","H02",Angle.degrees(104.5))
   
#finally, write out the converted solvent box
PDB().write(tip4pwaters,"spcbox.pdb")

print("Done!")
