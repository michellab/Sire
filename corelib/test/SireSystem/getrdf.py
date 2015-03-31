
import sys

if len(sys.argv) < 2:
   print("USAGE: getrdf.py file.pdb (file.xsc)", file=sys.stderr)
   sys.exit(-1)

from Sire.Mol import *
from Sire.IO import *
from Sire.System import *
from Sire.MM import *
from Sire.Vol import *
from Sire.Maths import *
from Sire.Units import *

#get the name of the PDB file containing the coordinates
pdbfile = sys.argv[1]

#get the name of the xscfile containing the system space
if len(sys.argv) >= 3:
    xscfile = sys.argv[2]
else:
    xscfile = None

#define the RDF to be collected - in this case;
#  RDF extends from 0 A to 10 A with a bin size of 0.05 A
#  RDF is collected between the atom called "O00" in the 
#  first molecule in the PDB file (molecule with index - MolIdx - 0),
#  and the atom called "H01" is all other molecules in the file
rdf = RDFMonitor( 0*angstrom, 10*angstrom, 0.05*angstrom )
rdf.add( MolIdx(0) + AtomName("O00"), AtomName("H01") )

def addRDF(rdf, pdbfile, xscfile=None):
    """Add to the RDF 'rdf' the distances calculated using the coordinates 
       from the PDB file 'pdbfile', using the xscfile 'xscfile' to get the 
       dimensions of the periodic box. If 'xscfile' is None, then 
       no periodic boundary conditions are used."""

    #first get the space in which to calculate intermolecular distances
    space = Cartesian()

    if xscfile:
        lines = open(xscfile,"r").readlines()

        words = lines[0].split()
        mincoords = Vector( float(words[0]), float(words[1]), float(words[2]) )
        maxcoords = Vector( float(words[3]), float(words[4]), float(words[5]) )

        space = PeriodicBox(mincoords, maxcoords)

    #now load all of the molecules
    mols = PDB().read(pdbfile)
                                                
    #create a system to hold the molecules, and add them
    system = System()
    system.add( MoleculeGroup("molecules", mols) )

    #give the space to the system
    system.add( InterCLJFF() )  # bug! need to add InterCLJFF so 
                                # that the system has a space property. This
                                # is fixed in new version of Sire, but is broken
                                # in your version

    system.setProperty("space", space)

    #add the RDF - this calculates the RDF for this PDB file and adds it to 'rdf'
    rdf.monitor(system)

#Now add the RDF for this file
addRDF(rdf, pdbfile, xscfile)

#You could have a loop here looping over lots of PDB files, e.g.
#if all of the PDB files were in 'pdbfiles' and all of the 
#xscfiles were in 'xscfiles' then you could write

#for i in range(0,nfiles):
#    addRDF(rdf, pdbfile[i], xscfile[i])

#Now print out the RDF
print("#RDF")

for i in range(0,rdf.nBins()):
    print(rdf[i].middle().to(angstrom), rdf[i].value())

