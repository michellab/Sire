#
# Evaluates free energy difference between two potential energy functions by
# use of the Zwanzig equation
#
# Based on Shirts et al. Accurate and Efficient Corrections for Missing Dispersion Interactions in Molecular Simulations
# J. Phys. Chem. B, 2007, 111 (45), pp 13052–13063


import os,sys, random
import math
from Sire.Tools.OpenMMMD import *
from Sire.Tools import Parameter, resolveParameters

# Python dependencies
#
try:
    mdtraj = Sire.try_import("mdtraj")
except:
    pass

try:
    numpy = Sire.try_import("numpy")
except:
    pass

import numpy as np

cutoff_type = Parameter("cutoff type", "cutoffperiodic", """The cutoff method to use during the simulation.""")

disable_crf = Parameter("disable crf", False,"""Whether to disable the reaction field crf term.""")

use_solute_inter_nrg = Parameter("use solute intermolecular energy",False,"""Whether to use the solute intermolecular energy instead of the total electrostatic potential energy.""")

trajfile = Parameter("trajfile", "traj000000001.dcd",
                    """File name of the trajectory to process.""")
                    
stepframe = Parameter("step_frame",1,"""The number of frames to step to between two succcessive evaluations.""")

def setupRFFF(system, space, cutoff=10* angstrom):

    print ("Creating force fields... ")

    solutes = system[MGName("solutes")]
    solute = system[MGName("solute_ref")]
    solute_hard = system[MGName("solute_ref_hard")]
    solute_todummy = system[MGName("solute_ref_todummy")]
    solute_fromdummy = system[MGName("solute_ref_fromdummy")]

    solvent = system[MGName("solvent")]

    #Solvent intramolecular CLJ energy
    solvent_intraclj = IntraCLJFF("solvent_intralj")
    solvent_intraclj.add(solvent)
    if (cutoff_type.val != "nocutoff"):
        solvent_intraclj.setUseReactionField(True)
        solvent_intraclj.setReactionFieldDielectric(rf_dielectric.val)
        solvent_intraclj.setDisableReactionFieldShift(disable_crf.val)

    #import pdb; pdb.set_trace()

    # Solvent-solvent LJ energy
    solventff = InterCLJFF("solvent:solvent")
    solventff.add(solvent)
    if (cutoff_type.val != "nocutoff"):
        solventff.setUseReactionField(True)
        solventff.setReactionFieldDielectric(rf_dielectric.val)
        solventff.setDisableReactionFieldShift(disable_crf.val)

    # Solute intramolecular LJ energy
    solute_hard_intraclj = IntraCLJFF("solute_hard_intralj")
    solute_hard_intraclj.add(solute_hard)
    if (cutoff_type.val != "nocutoff"):
        solute_hard_intraclj.setUseReactionField(True)
        solute_hard_intraclj.setReactionFieldDielectric(rf_dielectric.val)
        solute_hard_intraclj.setDisableReactionFieldShift(disable_crf.val)

    solute_todummy_intraclj = IntraSoftCLJFF("solute_todummy_intralj")
    #solute_todummy_intraclj.setShiftDelta(shift_delta.val)
    #solute_todummy_intraclj.setCoulombPower(coulomb_power.val)
    solute_todummy_intraclj.add(solute_todummy)
    if (cutoff_type.val != "nocutoff"):
        solute_todummy_intraclj.setUseReactionField(True)
        solute_todummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
        solute_todummy_intraclj.setDisableReactionFieldShift(disable_crf.val)

    solute_fromdummy_intraclj = IntraSoftCLJFF("solute_fromdummy_intralj")
    #solute_fromdummy_intraclj.setShiftDelta(shift_delta.val)
    #solute_fromdummy_intraclj.setCoulombPower(coulomb_power.val)
    solute_fromdummy_intraclj.add(solute_fromdummy)
    if (cutoff_type.val != "nocutoff"):
        solute_fromdummy_intraclj.setUseReactionField(True)
        solute_fromdummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
        solute_fromdummy_intraclj.setDisableReactionFieldShift(disable_crf.val)


    solute_hard_todummy_intraclj = IntraGroupSoftCLJFF("solute_hard:todummy_intralj")
    #solute_hard_todummy_intraclj.setShiftDelta(shift_delta.val)
    #solute_hard_todummy_intraclj.setCoulombPower(coulomb_power.val)
    solute_hard_todummy_intraclj.add(solute_hard, MGIdx(0))
    solute_hard_todummy_intraclj.add(solute_todummy, MGIdx(1))
    if (cutoff_type.val != "nocutoff"):
        solute_hard_todummy_intraclj.setUseReactionField(True)
        solute_hard_todummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
        solute_hard_todummy_intraclj.setDisableReactionFieldShift(disable_crf.val)

    solute_hard_fromdummy_intraclj = IntraGroupSoftCLJFF("solute_hard:fromdummy_intralj")
    #solute_hard_fromdummy_intraclj.setShiftDelta(shift_delta.val)
    #solute_hard_fromdummy_intraclj.setCoulombPower(coulomb_power.val)
    solute_hard_fromdummy_intraclj.add(solute_hard, MGIdx(0))
    solute_hard_fromdummy_intraclj.add(solute_fromdummy, MGIdx(1))
    if (cutoff_type.val != "nocutoff"):
        solute_hard_fromdummy_intraclj.setUseReactionField(True)
        solute_hard_fromdummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
        solute_hard_fromdummy_intraclj.setDisableReactionFieldShift(disable_crf.val)

    solute_todummy_fromdummy_intraclj = IntraGroupSoftCLJFF("solute_todummy:fromdummy_intralj")
    #solute_todummy_fromdummy_intraclj.setShiftDelta(shift_delta.val)
    #solute_todummy_fromdummy_intraclj.setCoulombPower(coulomb_power.val)
    solute_todummy_fromdummy_intraclj.add(solute_todummy, MGIdx(0))
    solute_todummy_fromdummy_intraclj.add(solute_fromdummy, MGIdx(1))
    if (cutoff_type.val != "nocutoff"):
        solute_todummy_fromdummy_intraclj.setUseReactionField(True)
        solute_todummy_fromdummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
        solute_todummy_fromdummy_intraclj.setDisableReactionFieldShift(disable_crf.val)

    # Solute-solvent LJ energy
    solute_hard_solventff = InterGroupCLJFF("solute_hard:solvent")
    solute_hard_solventff.add(solute_hard, MGIdx(0))
    solute_hard_solventff.add(solvent, MGIdx(1))
    if (cutoff_type.val != "nocutoff"):
        solute_hard_solventff.setUseReactionField(True)
        solute_hard_solventff.setReactionFieldDielectric(rf_dielectric.val)
        solute_hard_solventff.setDisableReactionFieldShift(disable_crf.val)

    solute_todummy_solventff = InterGroupSoftCLJFF("solute_todummy:solvent")
    solute_todummy_solventff.add(solute_todummy, MGIdx(0))
    solute_todummy_solventff.add(solvent, MGIdx(1))
    if (cutoff_type.val != "nocutoff"):
        solute_todummy_solventff.setUseReactionField(True)
        solute_todummy_solventff.setReactionFieldDielectric(rf_dielectric.val)
        solute_todummy_solventff.setDisableReactionFieldShift(disable_crf.val)

    solute_fromdummy_solventff = InterGroupSoftCLJFF("solute_fromdummy:solvent")
    solute_fromdummy_solventff.add(solute_fromdummy, MGIdx(0))
    solute_fromdummy_solventff.add(solvent, MGIdx(1))
    if (cutoff_type.val != "nocutoff"):
        solute_fromdummy_solventff.setUseReactionField(True)
        solute_fromdummy_solventff.setReactionFieldDielectric(rf_dielectric.val)
        solute_fromdummy_solventff.setDisableReactionFieldShift(disable_crf.val)

    # TOTAL
    forcefields = [solute_hard_intraclj, solute_todummy_intraclj, solute_fromdummy_intraclj,
                   solute_hard_todummy_intraclj, solute_hard_fromdummy_intraclj,
                   solute_todummy_fromdummy_intraclj,
                   solventff, solvent_intraclj,
                   solute_hard_solventff, solute_todummy_solventff, solute_fromdummy_solventff]

    for forcefield in forcefields:
        system.add(forcefield)

    system.setProperty("space", space)
    system.setProperty("switchingFunction", CHARMMSwitchingFunction(cutoff))
    system.setProperty("combiningRules", VariantProperty(combining_rules.val))
    system.setProperty("coulombPower", VariantProperty(coulomb_power.val))
    system.setProperty("shiftDelta", VariantProperty(shift_delta.val))


    
    solinter_nrg = solute_hard_solventff.components().coulomb() + \
        solute_todummy_solventff.components().coulomb(0) + \
        solute_fromdummy_solventff.components().coulomb(0)

    total_nrg = solute_hard_intraclj.components().coulomb() + \
                solute_todummy_intraclj.components().coulomb(0) + solute_fromdummy_intraclj.components().coulomb(0) + \
                solute_hard_todummy_intraclj.components().coulomb(0) + solute_hard_fromdummy_intraclj.components().coulomb(0) + \
                solute_todummy_fromdummy_intraclj.components().coulomb(0) + \
                solventff.components().coulomb() + \
                solvent_intraclj.components().coulomb() + \
                solute_hard_solventff.components().coulomb() + \
                solute_todummy_solventff.components().coulomb(0) + \
                solute_fromdummy_solventff.components().coulomb(0)



    e_total = system.totalComponent()

    lam = Symbol("lambda")

    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    print(use_solute_inter_nrg.val)

    if  (use_solute_inter_nrg.val):
        print("Using only solute intermolecular interactions...")
        system.setComponent(e_total,solinter_nrg)
        
    else:
        print("Using the total electrostatic energy...")
        system.setComponent(e_total, total_nrg)


    system.setConstant(lam, 0.0)

    system.add(PerturbationConstraint(solutes))

    # NON BONDED Alpha constraints for the soft force fields
    system.add(PropertyConstraint("alpha0", FFName("solute_todummy_intralj"), lam))
    system.add(PropertyConstraint("alpha0", FFName("solute_fromdummy_intralj"), 1 - lam))
    system.add(PropertyConstraint("alpha0", FFName("solute_hard:todummy_intralj"), lam))
    system.add(PropertyConstraint("alpha0", FFName("solute_hard:fromdummy_intralj"), 1 - lam))
    system.add(PropertyConstraint("alpha0", FFName("solute_todummy:fromdummy_intralj"), Max(lam, 1 - lam)))
    system.add(PropertyConstraint("alpha0", FFName("solute_todummy:solvent"), lam))
    system.add(PropertyConstraint("alpha0", FFName("solute_fromdummy:solvent"), 1 - lam))

    system.setComponent(lam, lambda_val.val)

    return system


def updateSystemfromTraj(system, frame_xyz, cell_lengths, cell_angles):
    print("Here we are processing traj")
    traj_coordinates = frame_xyz[0]

    traj_box_x = cell_lengths[0][0].tolist()
    traj_box_y = cell_lengths[0][1].tolist()
    traj_box_z = cell_lengths[0][2].tolist()

    traj_natoms = len(traj_coordinates)

    # Sire does not support non rectangular boxes
    newmols_coords = {}

    traj_index = 0
    mol_index = 0

    molnums = system.molNums()
    molnums.sort()
    for molnum in molnums:
        #print(molnum)
        mol = system.molecule(molnum)[0].molecule() #-1 to take all the solute atoms
        #print(mol)
        molatoms = mol.atoms()
        molnatoms = mol.nAtoms()
        # Create an empty coord group using molecule so we get the correct layout
        newmol_coords = AtomCoords( mol.property("coordinates") )
        for x in range(0,molnatoms):
            tmparray = traj_coordinates[traj_index]
            atom_coord = Vector( tmparray[0].tolist() , tmparray[1].tolist() , tmparray[2].tolist() )
            atom = molatoms[x]
            cgatomidx = atom.cgAtomIdx()
            newmol_coords.set( cgatomidx, atom_coord)
            traj_index += 1

        newmols_coords[molnum] = newmol_coords
        mol_index += 1

    if traj_natoms != traj_index:
        print ("The number of atoms in the system is not equal to the number of atoms in the trajectory file ! Aborting.")
        sys.exit(-1)
        
    changedmols = MoleculeGroup("changedmols")
    mol_index = 0
    for molnum in molnums:
        mol = system.molecule(molnum)[0].molecule()
        newmol_coords = newmols_coords[molnum]
        mol = mol.edit().setProperty("coordinates", newmol_coords).commit()
        changedmols.add(mol)
        
    system.update(changedmols)
    space = PeriodicBox(Vector( traj_box_x, traj_box_y, traj_box_z ) )
    system.setProperty("space",space)

    return system

def getFreeEnergy(delta_nrgs):
    free_nrg = FreeEnergyAverage(temperature.val)
    for nrg in delta_nrgs:
        free_nrg.accumulate(nrg.value())
    deltaG = free_nrg.average() * kcal_per_mol
    return deltaG

def resample(values):
    nvals = len(values)
    new_values = []
    for x in range(0,nvals):
        i = random.randint(0,nvals-1)
        new_values.append(values[i])
    return new_values

@resolveParameters
def runLambda():
    try:
        host = os.environ['HOSTNAME']
    except KeyError:
        host = "unknown"
    print("### Running RF tail correction calculation on %s ###" % host)
    if True:#verbose.val:
        print("###================= Simulation Parameters=====================###")
        Parameter.printAll()
        print ("###===========================================================###\n")
    print("lambda is %s" % lambda_val.val)

    if os.path.exists(s3file.val):
        (molecules, space) = Sire.Stream.load(s3file.val)
    else:
        amber = Amber()
        (molecules, space) = amber.readCrdTop(crdfile.val, topfile.val)
        Sire.Stream.save((molecules, space), s3file.val)

    system = createSystemFreeEnergy(molecules)

    # THIS IS THE ONE WITH SHORT CUTOFF
    system_shortc = System()
    system_shortc.copy(system)
    #import pdb; pdb.set_trace()
    # Construct a forcefield to compute electrostatic energies
    system_shortc = setupRFFF(system_shortc, space, \
                              cutoff=cutoff_dist.val)
    #import pdb; pdb.set_trace()

    #Compute the greatest box dimension and multiply it by sqrt(3), so to encompass
    #the box into a sphere
    #load and read the trajectory
    start_frame = 1
    end_frame = 1000000000
    step_frame = stepframe.val

    mdtraj_trajfile = mdtraj.open(trajfile.val,'r')
    nframes = len(mdtraj_trajfile)
    if end_frame > (nframes - 1):
        end_frame = nframes - 1
    mdtraj_trajfile.seek(start_frame)
    current_frame = start_frame
    #set a maximum initial value
    maximum = 0.0
    for framenumber in range(0,nframes):
        #for each frame extract the box length (cell_lengths)
        current, cell_lengths, angles = mdtraj_trajfile.read(n_frames=1)
        try :
            box_lengths = cell_lengths[0]
            for length in box_lengths:
                if length> maximum: 
                    maximum = length
        except :
            pass

    #reset the trajectory to the start frame
    mdtraj_trajfile.seek(start_frame)
    #now create a list of radius of cutoff
    maximum = maximum*math.sqrt(3) 
    long_cutoff = np.linspace(cutoff_dist.val.value(),maximum,20)
    print("studying these cutoffs")
    print(long_cutoff)
    ofile = open("electrostatics.dat","w")

    for l_cutoff in long_cutoff:
        l_angstrom = l_cutoff*angstrom
        print("Cutoff for long system:")
        print(l_angstrom)
        #create the "long" System
        system_longc = System()
        system_longc.copy(system)
        # Update this to setupRFFF
        system_longc = setupRFFF(system_longc, space, \
                                cutoff=l_angstrom)
        # TODO) Add code to compute Langevin Dipole term
        #
        #   
        #

        # Now loop over snapshots in dcd and accumulate energies
        delta_nrgs = []

        while (current_frame <= end_frame):
            print ("Processing frame %s " % current_frame)
            print ("CURRENT POSITION %s " % mdtraj_trajfile.tell() )
            frames_xyz, cell_lengths, cell_angles = mdtraj_trajfile.read(n_frames=1)
            system_shortc = updateSystemfromTraj(system_shortc, frames_xyz, cell_lengths, cell_angles)
            system_longc = updateSystemfromTraj(system_longc, frames_xyz, cell_lengths, cell_angles)
            print (system_shortc.energy())
            print (system_longc.energy())
            # TODO) add LD term to longc system
            #delta_nrg = (system_longc.energy()+E_lrc_full - system_shortc.energy())
            delta_nrg = (system_longc.energy() - system_shortc.energy())
            delta_nrgs.append(delta_nrg)
            current_frame += step_frame
            mdtraj_trajfile.seek(current_frame)
            #import pdb; pdb.set_trace()
        #print (delta_nrgs)
        # Now compute free energy change
        deltaG = getFreeEnergy(delta_nrgs)
        #print (deltaG)
        nbootstrap = 100
        deltaG_bootstrap = np.zeros(nbootstrap)
        for x in range(0,nbootstrap):
            resampled_nrgs = resample(delta_nrgs)
            dG = getFreeEnergy(resampled_nrgs)
            deltaG_bootstrap[x] = dG.value()
        dev = deltaG_bootstrap.std()
        #save teh value onto a file
        ofile.write("%.2f, %.4f,%.4f\n" % (l_cutoff,deltaG.value(),dev))
        print ("Rc: %.2f, DG_ELEC = %8.5f +/- %8.5f kcal/mol (1 sigma) " % (l_cutoff,deltaG.value(), dev))
        mdtraj_trajfile.seek(start_frame)
        current_frame=0
