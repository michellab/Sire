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
from Sire.Qt import *

import os,sys
import shutil

temperature = 25*celsius
nmoves = 100000
nblocks = 10

solute_file = "test/io/triatomic.pdb"
solute_name = "triatomic"
solute_params = "test/io/triatomic.ff"

protoms_dir = "%s/Work/Protoms" % os.getenv("HOME")

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

def histoAnalysis(nbins=1, values=None, output="BOND.dat",
mode="bond"):
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
   stream.write("# distribution\n")
   for x in range(0,nbins-1):
       bin = histo[x]
       stream.write("%8.3f %8.3f\n" % (bin.middle(),bin.value()/
count))



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

solute_intra_moves = ZMatMove(solutemolecules)
solute_intra_moves.setTemperature(temperature)
moves = SameMoves(solute_intra_moves)

#solute_mover_moves = MoverMove(solutemolecules)
#solute_mover_moves.setTemperature(temperature)
#moves = SameMoves(solute_mover_moves)

# Create monitors

r12 = Symbol("r12")
system.add( DistanceComponent( r12, solute.atom(AtomName("A1")),
                                    solute.atom(AtomName("A2")) ) )
system.add( "bond12", MonitorComponent(r12, RecordValues()), 5)

r23 = Symbol("r23")
system.add( DistanceComponent( r23, solute.atom(AtomName("A2")), solute.atom(AtomName("A3")) ) )
system.add( "bond23", MonitorComponent(r23, RecordValues()), 5)

theta123 = Symbol("theta123")
system.add( AngleComponent( theta123, solute.atom(AtomName("A1")),
                                      solute.atom(AtomName("A2")),
                                      solute.atom(AtomName("A3")) ) )

system.add( "angle123", MonitorComponent(theta123, RecordValues()), 5)

#system.add( "trajectory", TrajectoryMonitor(system[MGIdx(0)]), 1 )

t = QTime()
t.start()

# Run the simulation
for i in range(0,nblocks):
   print("Running nmove...")
   t.restart()
   system = moves.move(system, nmoves, True)
   print("...took %d ms" % t.elapsed())

print("Analysis...")
# make an histogram...
nbins=31
histoAnalysis(nbins=nbins,
              values=list(system.monitor( MonitorName("bond12") ).accumulator().values()),
              mode="bond",
              output="BOND-A1A2.dat")

histoAnalysis(nbins=nbins,
              values=list(system.monitor( MonitorName("bond23") ).accumulator().values()),
              mode="bond",
              output="BOND-A2A3.dat")

histoAnalysis(nbins=nbins,
              values=list(system.monitor( MonitorName("angle123") ).accumulator().values()),
              mode="angle",
              output="ANGLE-A1A2A3.dat")

#**  Output the trajectory
#system.monitor( MonitorName("trajectory") ).writeToDisk("outputXXXXXX.pdb")
