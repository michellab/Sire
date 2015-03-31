from Sire.IO import *
from Sire.Mol import *
from Sire.CAS import *
from Sire.Move import *
from Sire.FF import *
from Sire.MM import *
from Sire.Units import *
from Sire.Vol import *
from Sire.Maths import *
from Sire.Base import *
from Sire.System import *

import os,sys
import shutil

#
# Compute the probability distribution of a torsion with a null potential. The distribution should 
# be even. 

temperature = 25*celsius
nmoves = 10000
nblocks = 10

solute_file = "test/io/ABCD.pdb"
solute_name = "ABCD"
solute_params = "test/io/ABCDmover.ff"

protoms_dir = os.getenv("PROTOMSDIR")

ff_parameters = [ "%s/parameter/amber99.ff" % protoms_dir,
                  "%s/parameter/gaff.ff" % protoms_dir ]

combining_rules = "arithmetic"

compress = "none"

def createProtoMS():
    """Create and return a ProtoMS parameteriser object"""
    protoms = ProtoMS("%s/protoms2" % protoms_dir)
    for params in ff_parameters:
         protoms.addParameterFile(params)
    return protoms

def histoAnalysis(nbins=1, values=None, output="BOND.dat", mode="bond", avgnrg=0.0):
    """Histogram analysis"""
    if mode == "angle":
        unit = degrees
    else:
        unit = angstrom
    values.sort()
    minval = convertTo(values[0], unit)
    maxval = convertTo(values[-1], unit)
    histo = Histogram(minval,maxval,nbins=nbins)

    for val in values:
        histo.accumulate( convertTo(val, unit) )

    count = 0.0
    for x in range(0,nbins-1):
        bin = histo[x]
        count += bin.value()

    outfile = output
    stream =  open(outfile,'w')
    stream.write("# distribution Average energy %8.5f kcal.mol-1\n" % avgnrg)
    sumbin = 0.0
    for x in range(0,nbins-1):
        bin = histo[x]
        sumbin += bin.value()/count
        stream.write("%8.5f %8.5f %8.5f\n" % (bin.middle(),bin.value()/count,sumbin))   

    

# Load and parameterise
solute = PDB().readMolecule(solute_file)
solute = solute.edit().rename(solute_name).commit()

protoms = createProtoMS()
protoms.addParameterFile(solute_params)
solute = protoms.parameterise(solute, ProtoMS.SOLUTE)

# Putting the solute in the system
system = System()
solutemolecules = MoleculeGroup("solute", solute)

system.add(solutemolecules)

# Setting a forcefield 
solute_intraff = InternalFF("solute_intraff")
solute_intraff.add(solutemolecules)
solute_intraclj = IntraCLJFF("solute_intraclj")
solute_intraclj.add(solute)

forcefields = [ solute_intraff, solute_intraclj ]

for forcefield in forcefields:
    system.add(forcefield)

# Setup the Moves 

#solute_intra_moves = ZMatMove(solutemolecules)
#solute_intra_moves.setTemperature(temperature)

solute_mover_moves = MoverMove(solutemolecules)
solute_mover_moves.setDihedrals(solute.property('connectivity').getDihedrals())

moves = SameMoves(solute_mover_moves)

# Create monitors 

#r12 = Symbol("r12")
#system.add( DistanceComponent( r12, solute.atom(AtomName("A1")), solute.atom(AtomName("A2")) ) )
#system.add( "bond12", MonitorComponent(r12, RecordValues()), 5)

#r23 = Symbol("r23")
#system.add( DistanceComponent( r23, solute.atom(AtomName("A2")), solute.atom(AtomName("A3")) ) )
#system.add( "bond23", MonitorComponent(r23, RecordValues()), 5)

#r34 = Symbol("r34")
#system.add( DistanceComponent( r34, solute.atom(AtomName("A3")), solute.atom(AtomName("A4")) ) )
#system.add( "bond34", MonitorComponent(r34, RecordValues()), 5)

#theta123 = Symbol("theta123")
#system.add( AngleComponent( theta123, solute.atom(AtomName("A1")), solute.atom(AtomName("A2")), solute.atom(AtomName("A3")) ) )
#system.add( "theta123", MonitorComponent(theta123, RecordValues()), 5)

#theta234 = Symbol("theta234")
#system.add( AngleComponent( theta234, solute.atom(AtomName("A2")), solute.atom(AtomName("A3")), solute.atom(AtomName("A4")) ) )
#system.add( "theta234", MonitorComponent(theta234, RecordValues()), 5)

phi1234 = Symbol("phi1234")
system.add( DihedralComponent( phi1234, solute.atom(AtomName("A1")), solute.atom(AtomName("A2")), solute.atom(AtomName("A3")), solute.atom(AtomName("A4")) ) )
system.add( "phi1234", MonitorComponent(phi1234, RecordValues()), 5)

e_total = system.totalComponent()
system.add( "total_energy", MonitorComponent(e_total, Average()) )

#system.add( "trajectory", TrajectoryMonitor(system[MGIdx(0)]), 1 )

# Run the simulation
for i in range(0,nblocks):
    print("Running nmove...")
    system = moves.move(system, nmoves, True)

print("Analysis...")
# make an histogram...
nbins=31
#histoAnalysis(nbins=nbins,
#              avgnrg = system.monitor( MonitorName("total_energy") ).accumulator().average(),
#              values=system.monitor( MonitorName("bond12") ).accumulator().values(),
#              mode="bond",
#              output="BOND-A1A2.dat")
#histoAnalysis(nbins=nbins,
#              avgnrg = system.monitor( MonitorName("total_energy") ).accumulator().average(),
#              values=system.monitor( MonitorName("bond23") ).accumulator().values(),
#              mode="bond",
#              output="BOND-A2A3.dat")
#histoAnalysis(nbins=nbins,
#              avgnrg = system.monitor( MonitorName("total_energy") ).accumulator().average(),
#              values=system.monitor( MonitorName("bond34") ).accumulator().values(),
#              mode="bond",
#              output="BOND-A3A4.dat")
#histoAnalysis(nbins=nbins,
#              avgnrg = system.monitor( MonitorName("total_energy") ).accumulator().average(),
#              values=system.monitor( MonitorName("theta123") ).accumulator().values(),
#              mode="angle",
#              output="ANGLE-A1A2A3.dat")
#histoAnalysis(nbins=nbins,
#              avgnrg = system.monitor( MonitorName("total_energy") ).accumulator().average(),
#              values=system.monitor( MonitorName("theta234") ).accumulator().values(),
#              mode="angle",
#              output="ANGLE-A2A3A4.dat")
histoAnalysis(nbins=nbins,
              avgnrg = system.monitor( MonitorName("total_energy") ).accumulator().average(),
              values=list(system.monitor( MonitorName("phi1234") ).accumulator().values()),
              mode="angle",
              output="DIHEDRAL-A1A2A3A4.dat")


#**  Output the trajectory
#system.monitor( MonitorName("trajectory") ).writeToDisk("outputXXXXXX.pdb")
