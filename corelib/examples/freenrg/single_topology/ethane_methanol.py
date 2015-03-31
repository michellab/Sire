
from Sire.IO import *
from Sire.Mol import *
from Sire.CAS import *
from Sire.System import *
from Sire.Move import *
from Sire.MM import *
from Sire.Units import *
from Sire.Vol import *
from Sire.Maths import *
from Sire.Base import *

import os

########################################
## SIMULATION PARAMETERS ###############
########################################

simulation_title = "Single topology relative hydration free energy of ethane and methanol"

solute_file = "ethane2methanol.pdb"
solute_name = "ethane2methanol"
solute_params = "ethane2methanol.par"
solvent_file = "boxT4P.pdb"
solvent_name = "T4P"

temperature = 25 * celsius
pressure = 1 * atm

lambda_values = [ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 ]

delta_lambda = 0.001

coulomb_cutoff = 10 * angstrom
coulomb_feather = 9.5 * angstrom
lj_cutoff = 10 * angstrom
lj_feather = 9.5 * angstrom

combining_rules = "arithmetic"

nsubmoves = 5000

random_seed = 140154

compress = "bzip2 -f"

########################################

def printComponents(energies):
    components = energies.symbols()

    components.sort( key=str )

    for component in components:
        print "%s == %f kcal mol-1" % (component, energies[component])

    print "\n",

def createSystem():
    protomsdir = "%s/Work/ProtoMS" % os.getenv("HOME")

    protoms = ProtoMS( "%s/protoms2" % protomsdir )

    protoms.addParameterFile( "%s/parameter/amber99.ff" % protomsdir )
    protoms.addParameterFile( "%s/parameter/solvents.ff" % protomsdir )
    protoms.addParameterFile( "%s/parameter/gaff.ff" % protomsdir )
    protoms.addParameterFile( solute_params )

    solute = PDB().readMolecule(solute_file)
    solute = solute.edit().rename(solute_name).commit()
    solute = protoms.parameterise(solute, ProtoMS.SOLUTE)

    perturbation = solute.property("perturbations")

    lam = Symbol("lambda")
    lam_fwd = Symbol("lambda_{fwd}")
    lam_bwd = Symbol("lambda_{bwd}")

    initial = Perturbation.symbols().initial()
    final = Perturbation.symbols().final()

    solute = solute.edit().setProperty("perturbations",
                perturbation.recreate( (1-lam)*initial + lam*final ) ).commit()

    solute_fwd = solute.edit().renumber().setProperty("perturbations",
                perturbation.substitute( lam, lam_fwd ) ).commit()
    solute_bwd = solute.edit().renumber().setProperty("perturbations",
                perturbation.substitute( lam, lam_bwd ) ).commit()

    solvent = PDB().read(solvent_file)

    tip4p = solvent.moleculeAt(0).molecule()
    tip4p = tip4p.edit().rename(solvent_name).commit()
    tip4p = protoms.parameterise(tip4p, ProtoMS.SOLVENT)

    tip4p_chgs = tip4p.property("charge")
    tip4p_ljs = tip4p.property("LJ")

    for i in range(0,solvent.nMolecules()):
        tip4p = solvent.moleculeAt(i).molecule()
        tip4p = tip4p.edit().rename(solvent_name) \
                            .setProperty("charge", tip4p_chgs) \
                            .setProperty("LJ", tip4p_ljs) \
                            .commit()

        solvent.update(tip4p)

    system = System()

    solutes = MoleculeGroup("solutes")
    solutes.add(solute)
    solutes.add(solute_fwd)
    solutes.add(solute_bwd)

    solvent = MoleculeGroup("solvent", solvent)

    all = MoleculeGroup("all")
    all.add(solutes)
    all.add(solvent)

    system.add(solutes)
    system.add(solvent)
    system.add(all)

    solventff = InterCLJFF("solvent:solvent")
    solventff.add(solvent)

    solute_intraff = InternalFF("solute_intraff")
    solute_intraff.add(solute)

    solute_fwd_intraff = InternalFF("solute_fwd_intraff")
    solute_fwd_intraff.add(solute_fwd)

    solute_bwd_intraff = InternalFF("solute_bwd_intraff")
    solute_bwd_intraff.add(solute_bwd)

    solute_intraclj = IntraCLJFF("solute_intraclj")
    solute_intraclj.add(solute)

    solute_fwd_intraclj = IntraCLJFF("solute_fwd_intraclj")
    solute_fwd_intraclj.add(solute_fwd)

    solute_bwd_intraclj = IntraCLJFF("solute_bwd_intraclj")
    solute_bwd_intraclj.add(solute_bwd)

    solute_solventff = InterGroupCLJFF("solute:solvent")
    solute_solventff.add(solute, MGIdx(0))
    solute_solventff.add(solvent, MGIdx(1))

    solute_fwd_solventff = InterGroupCLJFF("solute_fwd:solvent")
    solute_fwd_solventff.add(solute_fwd, MGIdx(0))
    solute_fwd_solventff.add(solvent, MGIdx(1))

    solute_bwd_solventff = InterGroupCLJFF("solute_bwd:solvent")
    solute_bwd_solventff.add(solute_bwd, MGIdx(0))
    solute_bwd_solventff.add(solvent, MGIdx(1))

    forcefields = [ solventff, solute_intraff, solute_intraclj, solute_solventff,
                               solute_fwd_intraff, solute_fwd_intraclj, solute_fwd_solventff,
                               solute_bwd_intraff, solute_bwd_intraclj, solute_bwd_solventff ]

    for forcefield in forcefields:
        system.add(forcefield)

    xsc_line = open(solvent_file, "r").readlines()[0]
    words = xsc_line.split()

    #HEADER box  -12.5  -12.5  -12.5   12.5   12.5   12.5
    space = PeriodicBox( Vector( float(words[2]), float(words[3]), float(words[4]) ),
                         Vector( float(words[5]), float(words[6]), float(words[7]) ) )

    system.setProperty( "space", space )
    system.setProperty( "switchingFunction", 
                        HarmonicSwitchingFunction(coulomb_cutoff, coulomb_feather,
                                              lj_cutoff, lj_feather) )
    system.setProperty( "combiningRules", VariantProperty(combining_rules) )

    e_total = system.totalComponent()
    e_fwd = Symbol("E_{fwd}")
    e_bwd = Symbol("E_{bwd}")

    total_nrg = solventff.components().total() + \
                solute_intraclj.components().total() + solute_intraff.components().total() + \
                solute_solventff.components().total()

    fwd_nrg = solventff.components().total() + \
              solute_fwd_intraclj.components().total() + solute_fwd_intraff.components().total() + \
              solute_fwd_solventff.components().total()

    bwd_nrg = solventff.components().total() + \
              solute_bwd_intraclj.components().total() + solute_bwd_intraff.components().total() + \
              solute_bwd_solventff.components().total()

    system.setComponent( e_total, total_nrg )
    system.setComponent( e_fwd, fwd_nrg )
    system.setComponent( e_bwd, bwd_nrg )

    system.setConstant(lam, 0.0)
    system.setConstant(lam_fwd, 0.0)
    system.setConstant(lam_bwd, 0.0)

    system.add( SpaceWrapper(Vector(0,0,0), all) )

    system.add( PerturbationConstraint(solutes) )

    system.add( ComponentConstraint( lam_fwd, Min( lam + delta_lambda, 1 ) ) )
    system.add( ComponentConstraint( lam_bwd, Max( lam - delta_lambda, 0 ) ) )

    de_fwd = Symbol("de_fwd")
    de_bwd = Symbol("de_bwd")

    system.setComponent( de_fwd, fwd_nrg - total_nrg )
    system.setComponent( de_bwd, total_nrg - bwd_nrg )

    system.add( "total_energy", MonitorComponent(e_total, Average()) )
    system.add( "de_fwd", MonitorComponent(de_fwd, FreeEnergyAverage(temperature)) )
    system.add( "de_bwd", MonitorComponent(de_bwd, FreeEnergyAverage(temperature)) )

    system.setComponent(lam, 0.0)
    print "LAMBDA=0   : Energy = %f kcal mol-1" % system.energy().to(kcal_per_mol)
    print "             (%f, %f)" % (system.energy(e_fwd).to(kcal_per_mol),
                                   system.energy(e_bwd).to(kcal_per_mol))

    system.setComponent(lam, 0.5)
    print "LAMBDA=0.5 : Energy = %f kcal mol-1" % system.energy().to(kcal_per_mol)
    print "             (%f, %f)" % (system.energy(e_fwd).to(kcal_per_mol),
                                   system.energy(e_bwd).to(kcal_per_mol))

    system.setComponent(lam, 1.0)
    print "LAMBDA=1.0 : Energy = %f kcal mol-1" % system.energy().to(kcal_per_mol)
    print "             (%f, %f)" % (system.energy(e_fwd).to(kcal_per_mol),
                                   system.energy(e_bwd).to(kcal_per_mol))
    return system


def createMoves(system):
    solutes = system[ MGName("solutes") ]
    solvent = system[ MGName("solvent") ]

    solute_mc_weight = 19
    solvent_mc_weight = 980
    volume_mc_weight = 1

    pref_constant = 200 * angstrom2

    max_solute_translation = 0 * angstrom
    max_solute_rotation = 5 * degrees

    max_solvent_translation = 0.15 * angstrom
    max_solvent_rotation = 15 * degrees

    max_volume_change = 0.1 * solvent.nMolecules() * angstrom3

    solvent_moves = RigidBodyMC( PrefSampler(solutes[MolIdx(0)], 
                                             solvent, pref_constant) )

    solvent_moves.setMaximumTranslation(max_solvent_translation)
    solvent_moves.setMaximumRotation(max_solvent_rotation)

    solute_moves = RigidBodyMC( solutes )
    solute_moves.setMaximumTranslation(max_solute_translation)
    solute_moves.setMaximumRotation(max_solute_rotation)
    solute_moves.setSynchronisedTranslation(True)
    solute_moves.setSynchronisedRotation(True)
    solute_moves.setSharedRotationCenter(True)

    solute_intra_moves = ZMatMove( solutes )
    solute_intra_moves.setSynchronisedMotion(True)

    volume_moves = VolumeMove()
    volume_moves.setMaximumVolumeChange(max_volume_change)

    weighted_moves = WeightedMoves()

    weighted_moves.add( solute_moves, solute_mc_weight / 2 )
    weighted_moves.add( solute_intra_moves, solute_mc_weight / 2 )
    weighted_moves.add( solvent_moves, solvent_mc_weight )
    weighted_moves.add( volume_moves, volume_mc_weight )

    weighted_moves.setTemperature(temperature)
    weighted_moves.setPressure(pressure)
    weighted_moves.setGenerator( RanGenerator(random_seed+3) )

    return weighted_moves


def createReplicas():
    system = createSystem()
    moves = createMoves(system)

    replicas = Replicas( len(lambda_values) )

    replicas.setSubSystem(system)
    replicas.setSubMoves(moves)
    replicas.setNSubMoves(nsubmoves)
    replicas.setRecordAllStatistics(True)
    replicas.setLambdaComponent( Symbol("lambda") )

    replicas.setGenerator( RanGenerator(random_seed+7) )
    replicas.packToDisk()

    for i in range(0, len(lambda_values)):
        replicas.setLambdaValue(i, lambda_values[i])

    replicas.add( "total_energy", MonitorMonitor(MonitorName("total_energy"), True) )
    replicas.add( "dg_fwd", MonitorMonitor(MonitorName("de_fwd"), True) )
    replicas.add( "dg_bwd", MonitorMonitor(MonitorName("de_bwd"), True) )
    replicas.add( "trajectory", TrajectoryMonitor(MGName("all")) )

    return replicas


def writeEnergies(energies, filename):
    """This function writes the energies in 'energies' to the file 'filename'"""

    lamvals = energies.keys()
    lamvals.sort()

    FILE = open(filename, "w")

    for lamval in lamvals:
        print >>FILE,"%f   %f"  % (lamval, energies[lamval])


def calculatePMF(gradients):
    """This function calculates and return the PMF given the passed series
       of lambda values and gradients"""
    
    pmf = {}

    lamvals = gradients.keys()
    lamvals.sort()
    
    if lamvals[0] != 0:
        #we need to start from 0
        gradients[0] = gradients[lamvals[0]]
        lamvals.insert(0, 0)
    
    if lamvals[-1] != 1:  
        #we need to end with 1
        gradients[1] = gradients[lamvals[-1]]
        lamvals.append(1)
    
    #start at 0
    pmf[ lamvals[0] ] = 0.0

    for i in range(1,len(lamvals)):
        last_lam = lamvals[i-1]
        this_lam = lamvals[i]

        delta_lam = this_lam - last_lam

        pmf[this_lam] = pmf[last_lam] + (delta_lam * 0.5 * (gradients[this_lam] + \
                                                            gradients[last_lam]))

    return pmf


def writeReplicaData(replicas, outdir, nmoves):
    """This function writes the contents of each replica to disk"""

    nreplicas = replicas.nReplicas()

    pdbs = {}
    total_energy = {}
    dg_f = {}
    dg_b = {}

    replica_ids = {}

    ids = replicas.replicaIDs()

    for i in range(0, nreplicas):
        replica = replicas[i]
        monitors = replica.monitors()
        lamval = replica.lambdaValue()

        pdbs[lamval] = monitors[MonitorName("trajectory")]

        total_energy[lamval] = monitors[MonitorName("total_energy")]
        replica_ids[lamval] = ids[i]

        dg_f[lamval] = monitors[MonitorName("dg_fwd")]
        dg_b[lamval] = monitors[MonitorName("dg_bwd")]

    #now clear all supra-level statistics
    replicas.clearStatistics()

    nrg_avg = {}

    dg_f_avg = {}
    dg_b_avg = {}

    for lamval in lambda_values:
        #create a directory for this lambda value
        lamdir = "%s/%f" % (outdir, lamval)
        id = replica_ids[lamval]

        if os.path.exists(lamdir):
            shutil.rmtree(lamdir)

        # Create this output directory
        os.makedirs(lamdir)         

        # Now write the PDB(s) to this directory
        pdbs[lamval].writeToDisk("%s/output_%003d_%0004d.pdb" % (lamdir,id,nmoves))

        dg_f_avg[lamval] = dg_f[lamval][-1].accumulator().average() / delta_lambda
        dg_b_avg[lamval] = dg_b[lamval][-1].accumulator().average() / delta_lambda

        nrg_avg[lamval] = total_energy[lamval][-1].accumulator().average()

    # write the average energies, gradients and PMF for each
    #Êlambda value to a file in the output directory
    writeEnergies(dg_f_avg, "%s/average_dg_f.txt" % outdir)
    writeEnergies(dg_b_avg, "%s/average_dg_b.txt" % outdir)

    writeEnergies(nrg_avg, "%s/average_energy.txt" % outdir)

    writeEnergies(calculatePMF(dg_f_avg), "%s/pmf_f.txt" % outdir)
    writeEnergies(calculatePMF(dg_b_avg), "%s/pmf_b.txt" % outdir)


def writeReplicaTrajectory(lamvals, replica_ids, trajfile):
    """This function appends the passed replica IDs with the passed
       lambda values to the passed replica trajectory file"""

    if os.path.exists(trajfile):
        FILE = open(trajfile, "a")
    else:
        FILE = open(trajfile, "w")
        print >>FILE,"#move  ",

        for i in range(0,len(lamvals)):
            print >>FILE,"%d  "  % i,

        print >>FILE,"\n",

    print >>FILE,"%d   " % moves.nMoves(),

    for id in replicas.replicaIDs():
        print >>FILE,"%f  " % lamvals[id],

    print >>FILE,"\n",

    FILE.close()


#########################################################
#########################################################
##
##
## ...now the actual work performed by the script...
##
##

import Sire.Stream
import shutil

try:
    # Get the name of the output directory
    outdir = "output"

    restart_file = "repex_restart.s3"

    # Is there an existing restart file (repex_restart.s3)
    if os.path.exists(restart_file):
        (replicas, moves) = Sire.Stream.load(restart_file)

        print "Restart file contains %d replicas, on which we have performed %d moves." % \
                    (replicas.nReplicas(), moves.nMoves())

    else:
        #There isn't - create the replicas using the
        # 'setupSimulation()' function (above)
        print "\n******** %s ********\n" % simulation_title

        print "Setting up the simulation..."
        replicas = createReplicas()

        moves = RepExMove()
        moves.setGenerator( RanGenerator(random_seed+1) )

        #Save these replicas now so that we don't have to set 
        #them up again
        print "Saving and reloading the system..."
        Sire.Stream.save( [replicas,moves], restart_file )

        #Make sure that we can read this file
        replicas = None
        moves = None
        (replicas, moves) = Sire.Stream.load( restart_file )

        #Now remove the output directory, as we will be creating it from scratch
        shutil.rmtree( outdir, True )

    # Get the name of the directory to place all of the output for this move
    outdir_i = "%s/iteration_%0004d" % (outdir, moves.nMoves() + 1)

    # If this output directory exists, then remove it (it is an old version)
    if os.path.exists(outdir_i):
        shutil.rmtree(outdir_i)

    # Create this output directory
    os.makedirs(outdir_i)

    # Now we have the replicas, perform a block of sampling
    # using a replica exchange move
    print "Running replica exchange block %d..." % (moves.nMoves()+1)
    sim = SupraSim.run( replicas, moves, 1, True )
    sim.wait()

    replicas = sim.system()
    moves = sim.moves()

    print "Complete! - aggregate acceptance ratio = %f %%\n" % \
                             (100 * moves[0].acceptanceRatio())

    # Now write out data about each individual replica
    print "Saving output to disk..."
    writeReplicaData(replicas, outdir_i, moves.nMoves())

    # Save the new replica and replica move data
    print "Saving the restart file..."
    if os.path.exists(restart_file):
        shutil.move(restart_file, "old_%s" % restart_file)

    Sire.Stream.save( (replicas, moves), restart_file )

    # Append the replica trajectory to the file
    print "Writing replica trajectory info..."
    trajfile = "%s/replica_trajectory.txt" % outdir
    writeReplicaTrajectory(lambda_values, replicas.replicaIDs(), trajfile)

    # Now compress everything in the output directory to save disk space
    for lamval in lambda_values:
        os.system( "%s %s/%f/*" % (compress, outdir_i, lamval) )

except Exception, e:
    #Something went wrong!!!
    import Sire.Error
    Sire.Error.printError(e)
