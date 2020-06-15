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
#try:
#    import mdtraj
#except ImportError:
#    print ("LJcutoff.py depends on a working install of the python module mdtraj. Please install mdtraj in your sire python.")
#    sys.exit(-1)
#
#try:
#    import numpy as np
#except ImportError:
#    print ("LJcutoff.py depends on a working install of the python module mdtraj. Please install mdtraj in your sire python.")
#    sys.exit(-1)


bulk_rho = Parameter("bulk_rho", 1.0 * gram/(centimeter*centimeter*centimeter)\
                     ,"""The density of buk solvent.""")

trajfile = Parameter("trajfile", "traj000000001.dcd",
                    """File name of the trajectory to process.""")
stepframe = Parameter("step_frame",1,"""The number of frames to step to between two succcessive evaluations.""")

def setupLJFF(system, space, cutoff=10* angstrom):

    print ("Creating force fields... ")

    solutes = system[MGName("solutes")]
    solute = system[MGName("solute_ref")]
    solute_hard = system[MGName("solute_ref_hard")]
    solute_todummy = system[MGName("solute_ref_todummy")]
    solute_fromdummy = system[MGName("solute_ref_fromdummy")]

    solvent = system[MGName("solvent")]

    #Solvent intramolecular CLJ energy
    solvent_intralj = IntraCLJFF("solvent_intralj")
    solvent_intralj.add(solvent)
    # Solvent-solvent LJ energy
    solventff = InterCLJFF("solvent:solvent")
    solventff.add(solvent)
    # Solute intramolecular LJ energy
    solute_hard_intralj = IntraCLJFF("solute_hard_intralj")
    solute_hard_intralj.add(solute_hard)

    solute_todummy_intralj = IntraSoftCLJFF("solute_todummy_intralj")
    solute_todummy_intralj.setShiftDelta(shift_delta.val)
    solute_todummy_intralj.setCoulombPower(coulomb_power.val)
    solute_todummy_intralj.add(solute_todummy)

    solute_fromdummy_intralj = IntraSoftCLJFF("solute_fromdummy_intralj")
    solute_fromdummy_intralj.setShiftDelta(shift_delta.val)
    solute_fromdummy_intralj.setCoulombPower(coulomb_power.val)
    solute_fromdummy_intralj.add(solute_fromdummy)

    solute_hard_todummy_intralj = IntraGroupSoftCLJFF("solute_hard:todummy_intralj")
    solute_hard_todummy_intralj.setShiftDelta(shift_delta.val)
    solute_hard_todummy_intralj.setCoulombPower(coulomb_power.val)
    solute_hard_todummy_intralj.add(solute_hard, MGIdx(0))
    solute_hard_todummy_intralj.add(solute_todummy, MGIdx(1))

    solute_hard_fromdummy_intralj = IntraGroupSoftCLJFF("solute_hard:fromdummy_intralj")
    solute_hard_fromdummy_intralj.setShiftDelta(shift_delta.val)
    solute_hard_fromdummy_intralj.setCoulombPower(coulomb_power.val)
    solute_hard_fromdummy_intralj.add(solute_hard, MGIdx(0))
    solute_hard_fromdummy_intralj.add(solute_fromdummy, MGIdx(1))

    solute_todummy_fromdummy_intralj = IntraGroupSoftCLJFF("solute_todummy:fromdummy_intralj")
    solute_todummy_fromdummy_intralj.setShiftDelta(shift_delta.val)
    solute_todummy_fromdummy_intralj.setCoulombPower(coulomb_power.val)
    solute_todummy_fromdummy_intralj.add(solute_todummy, MGIdx(0))
    solute_todummy_fromdummy_intralj.add(solute_fromdummy, MGIdx(1))

    # Solute-solvent LJ energy
    solute_hard_solventff = InterGroupCLJFF("solute_hard:solvent")
    solute_hard_solventff.add(solute_hard, MGIdx(0))
    solute_hard_solventff.add(solvent, MGIdx(1))

    solute_todummy_solventff = InterGroupSoftCLJFF("solute_todummy:solvent")
    solute_todummy_solventff.add(solute_todummy, MGIdx(0))
    solute_todummy_solventff.add(solvent, MGIdx(1))

    solute_fromdummy_solventff = InterGroupSoftCLJFF("solute_fromdummy:solvent")
    solute_fromdummy_solventff.add(solute_fromdummy, MGIdx(0))
    solute_fromdummy_solventff.add(solvent, MGIdx(1))

    # TOTAL
    forcefields = [solute_hard_intralj, solute_todummy_intralj, solute_fromdummy_intralj,
                   solute_hard_todummy_intralj, solute_hard_fromdummy_intralj,
                   solute_todummy_fromdummy_intralj,
                   solventff, solvent_intralj,
                   solute_hard_solventff, solute_todummy_solventff, solute_fromdummy_solventff]

    for forcefield in forcefields:
        system.add(forcefield)

    system.setProperty("space", space)
    system.setProperty("switchingFunction", CHARMMSwitchingFunction(cutoff))
    system.setProperty("combiningRules", VariantProperty(combining_rules.val))
    system.setProperty("coulombPower", VariantProperty(coulomb_power.val))
    system.setProperty("shiftDelta", VariantProperty(shift_delta.val))

    # TOTAL
    #total_nrg = solute_hard_intralj.components().total() + \
    #            solute_todummy_intralj.components().total(0) + solute_fromdummy_intralj.components().total(0) + \
    #            solute_hard_todummy_intralj.components().total(0) + solute_hard_fromdummy_intralj.components().total(0) + \
    #            solute_todummy_fromdummy_intralj.components().total(0) + \
    #            solventff.components().total() + \
    #            solvent_intralj.components().total() + \
    #            solute_hard_solventff.components().total() + \
    #            solute_todummy_solventff.components().total(0) + \
    #            solute_fromdummy_solventff.components().total(0)

    total_nrg = solute_hard_intralj.components().lj() + \
                solute_todummy_intralj.components().lj(0) + solute_fromdummy_intralj.components().lj(0) + \
                solute_hard_todummy_intralj.components().lj(0) + solute_hard_fromdummy_intralj.components().lj(0) + \
                solute_todummy_fromdummy_intralj.components().lj(0) + \
                solventff.components().lj() + \
                solvent_intralj.components().lj() + \
                solute_hard_solventff.components().lj() + \
                solute_todummy_solventff.components().lj(0) + \
                solute_fromdummy_solventff.components().lj(0)

    e_total = system.totalComponent()

    lam = Symbol("lambda")

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


def addAnalyticalLRC(system, cutoff, bulk_density):

    # 1) Collect all solvent particles
    # 2) Average all solvent sigma/epsilon parameters
    solvent = system[MGName("solvent")]
    solvent_mols = solvent.molecules()
    solvent_molnums = solvent_mols.molNums()
    #
    # What if solvent contains more than one type of molecule?
    #
    #for molnum in solvent_molnums:
    avg_sigma = 0.0 * angstrom
    avg_epsilon = 0.0 * kcal_per_mol
    LJsites = 0
    # Is this always the right molecule?
    mol = solvent_mols.first()[0].molecule()
    if (mol.nAtoms() == 1):
        print ("This does not seem to be a solvent molecule. Picking up another one...")
        mol = solvent_mols.last().molecule()
        if (mol.nAtoms() == 1):
            print ("This also does not seem to be a solvent molecule ! Abort!")
            sys.exit(-1)
    #molecule(molnum).molecule()
    atoms = mol.atoms()
    natoms = mol.nAtoms()
    solv_mol_mass = 0 * g_per_mol
    for x in range(0,natoms):
        atom = atoms[x]
        at_mass = atom.property("mass")
        solv_mol_mass += at_mass
    LJparams = mol.property("LJ").toVector()
    for LJparam in LJparams:
        sigma = LJparam.sigma()
        epsilon = LJparam.epsilon()
        print (sigma, epsilon)
        # Are we supposed to skip null LJ sites?
        if epsilon.value() > 0:
            avg_epsilon += epsilon
            avg_sigma += sigma
            LJsites += 1
    avg_sigma /= LJsites
    avg_epsilon /= LJsites
    print (avg_sigma, avg_epsilon, LJsites)
    solv_sigma = avg_sigma
    solv_epsilon = avg_epsilon
    # Now for each LJsite in every molecule in the system, compute LRC term
    molecules_group = system[MGName("molecules")]
    # IS THE FORMULA CORRECT FOR SOFT-CORE POTENTIALS ? IN PRINCIPLE NO ERROR
    # AT THE END STATES
    molecules = molecules_group.molecules()
    molnums = molecules.molNums()
    E_lrc_full = 0.0 * kcal_per_mol
    for molnum in molnums:
        mol = molecules.molecule(molnum)[0].molecule()
        LJparams = mol.property("LJ").toVector()
        for LJparam in LJparams:
            sigma = LJparam.sigma()
            epsilon = LJparam.epsilon()
            #import pdb; pdb.set_trace()
            eps_ij = math.sqrt(epsilon.value()*solv_epsilon.value()) \
                     * kcal_per_mol
            if combining_rules.val == 'arithmetic':
                sig_ij = (0.5*(sigma+solv_sigma)).value()
            else:
                sig_ij = math.sqrt(sigma.value()*solv_sigma.value())
            sig_ij6 = sig_ij**6
            sig_ij12 = sig_ij6**2
            sig_ij6 = sig_ij6 #* angstrom
            sig_ij12 = sig_ij12 #* angstrom
            # units !!
            # density must be converted in molecule per cubic Angstrom
            rho = (bulk_density/solv_mol_mass).value()
            cutval = cutoff.value()
            #cutval = 3.15075
            E_lrc_pairwise = 8*pi*rho*\
                             ( (1/(9.*(cutval)**9))*(eps_ij*sig_ij12) -
                               (1/(3.*(cutval)**3))*(eps_ij*sig_ij6) )
            #print (E_lrc_pairwise)
            #import pdb; pdb.set_trace()
            #sys.exit(-1)
            E_lrc_full += E_lrc_pairwise
    return E_lrc_full

def zeroCharges(system):

    molecules = system[MGName("molecules")]
    molnums = molecules.molNums()
    #updated = []
    for molnum in molnums:
        mol = molecules.molecule(molnum)[0].molecule()
        editmol = mol.edit()
        for x in range(0,mol.nAtoms()):
            editatom = editmol.atom(AtomIdx(x))
            editatom.setProperty("charge", 0 * mod_electron)
            editmol = editatom.molecule()
        mol = editmol.commit()
        system.update(mol)
    return system

def updateSystemfromTraj(system, frame_xyz, cell_lengths, cell_angles):
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
        mol = system.molecule(molnum)[0].molecule()
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
    print("### Running LJ tail correction calculation on %s ###" % host)
    if verbose.val:
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

    # !!! NEED TO DISABLE CHANGE IN COULOMBIC CUTOFF !!
    #system = zeroCharges(system)

    #import pdb; pdb.set_trace()

    # THIS IS THE ONE WITH SHORT CUTOFF
    system_shortc = System()
    system_shortc.copy(system)
    #import pdb; pdb.set_trace()
    system_shortc = setupLJFF(system_shortc, space, \
                              cutoff=cutoff_dist.val)
    #import pdb; pdb.set_trace()

    # Determine longest cutoff that can be used. Take lowest space dimension,
    # and decrease by 5%
    dims = space.dimensions()
    mindim = dims.x()
    if mindim > dims.y():
        mindim = dims.y()
    if mindim > dims.z():
        mindim = dims.z()
    long_cutoff = (mindim/2.0 * 0.95) * angstrom
    print (long_cutoff)
    system_longc = System()
    system_longc.copy(system)
    system_longc = setupLJFF(system_longc, space, \
                             cutoff=long_cutoff)
    # NOW ADD ANALYTICAL CORRECTION TERM TO longc
    E_lrc_full = addAnalyticalLRC(system_longc, long_cutoff, bulk_rho.val)
    #import pdb; pdb.set_trace()
    # Now loop over snapshots in dcd and accumulate energies
    start_frame = 1
    end_frame = 1000000000
    step_frame = stepframe.val

    #mdtraj_top = mdtraj.load_prmtop(topfile.val)
    mdtraj_trajfile = mdtraj.open(trajfile.val,'r')
    nframes = len(mdtraj_trajfile)
    if end_frame > (nframes - 1):
        end_frame = nframes - 1
    mdtraj_trajfile.seek(start_frame)
    current_frame = start_frame

    delta_nrgs = []

    while (current_frame <= end_frame):
        print ("Processing frame %s " % current_frame)
        print ("CURRENT POSITION %s " % mdtraj_trajfile.tell() )
        frames_xyz, cell_lengths, cell_angles = mdtraj_trajfile.read(n_frames=1)
        #print (system_shortc.energy())
        #print (system_longc.energy())
        system_shortc = updateSystemfromTraj(system_shortc, frames_xyz, cell_lengths, cell_angles)
        system_longcc = updateSystemfromTraj(system_longc, frames_xyz, cell_lengths, cell_angles)
        #print (system_shortc.energy())
        #print (system_longc.energy())
        delta_nrg = (system_longc.energy()+E_lrc_full - system_shortc.energy())
        delta_nrgs.append(delta_nrg)
        current_frame += step_frame
        mdtraj_trajfile.seek(current_frame)
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
    print ("DG_LJ = %8.5f +/- %8.5f kcal/mol (1 sigma) " % (deltaG.value(), dev))
