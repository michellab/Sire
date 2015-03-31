
from Sire.IO import *
from Sire.MM import *
from Sire.Mol import *
from Sire.Move import *
from Sire.MM import *
from Sire.System import *
from Sire.CAS import *
from Sire.Vol import *
from Sire.Units import *

import Sire.Stream

protodir = "/Users/chris/Work/ProtoMS/"

print("Parameterising the oscillator...")

oscillator = PDB().readMolecule("test/io/oscillator.pdb")
oscillator = oscillator.edit().rename("harm1t2").commit()

protoms = ProtoMS("%s/protoms2" % protodir)

protoms.addParameterFile( "test/io/oscillators.ff" )

oscillator = protoms.parameterise(oscillator, ProtoMS.SOLUTE)

print("...parameterisation complete!")

internalff = InternalFF("InternalFF")

internalff.add( oscillator )

system = System()
system.add( internalff )

lam = Symbol("lambda")

system.setComponent(lam, 0.01)

system.setComponent(system.totalComponent(), lam * internalff.components().total())

system.add( "average energy", MonitorComponent(system.totalComponent()) )

lambda_values = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 ]

replicas = Replicas(system, len(lambda_values))

replicas.setLambdaComponent(lam)

for i in range(0,len(lambda_values)):
    replicas.setLambdaValue(i, lambda_values[i])

zmatmove = ZMatMove( internalff.groups()[0] )
zmatmove.setTemperature( 298 * kelvin )

nsubmoves = 1000

replicas.setSubMoves( SameMoves(zmatmove) )
replicas.setNSubMoves(nsubmoves)

# Average energy should be 1/2 kT
theo_nrg = 0.5 * gasr * 298

print("Running a simulation - initial energy = %f kcal mol-1" % system.energy().to(kcal_per_mol))

repexmove = RepExMove()

lambda_trajectory = []

i = -1

def printInfo(replicas):

    lamtraj = replicas.lambdaTrajectory()
    lambda_trajectory.append(lamtraj)

    ids = replicas.replicaIDs()

    for j in range(0,replicas.count()):
        replica = replicas[j]
        system = replica.subSystem()
        zmatmove = replica.subMoves()[0]

        print("Replica %d: lambda = %f: ID = %d" % (j, replica.lambdaValue(), ids[j]))
        print("%d : Energy = %f kcal mol-1" % ( (i+1)*replica.nSubMoves(), \
                                            system.energy().to(kcal_per_mol) ))

        avg_nrg = system[MonitorName("average energy")].accumulator().average()

        print("      Average energy = %f kcal mol-1 : error = %f kcal mol-1" % (avg_nrg, theo_nrg-avg_nrg))
        print("      Acceptance ratio = %f %%" % (100 * zmatmove.acceptanceRatio()))

printInfo(replicas)

#Sire.Stream.save( (replicas, repexmove), "test/SireMove/reti_oscillator.s3" )

for i in range(0,10):
    sim = SupraSim.run( replicas, repexmove, 1 )
    
    replicas = sim.system()
    repexmove = sim.moves()[0]

    printInfo(replicas)

    print("      Replica exchange acceptance ratio: %f %%" % (100*repexmove.acceptanceRatio()))

print("\nReplica trajectory")
for i in range(0, len(lambda_trajectory)):
    print("%d " % i, end=' ')

    for lamval in lambda_trajectory[i]:
        print("%f " % lamval, end=' ')

    print("\n", end=' ') 
