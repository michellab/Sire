#
# Python script to perform a MD simulation in Sire with OpenMM
#

import os,re, sys, shutil
import math

from Sire.IO import *
from Sire.Mol import *
from Sire.CAS import *
from Sire.System import *
from Sire.Move import *
from Sire.MM import *
from Sire.FF import *
from Sire.Units import *
from Sire.Vol import *
from Sire.Maths import *
from Sire.Base import *
from Sire.Qt import *
from Sire.ID import *
from Sire.Config import *
from Sire.Analysis import *

from Sire.Tools.DCDFile import *

from Sire.Tools import Parameter, resolveParameters

import Sire.Stream

gpu = Parameter("gpu", 0, """The device ID of the GPU on which to run the simulation.""")

rf_dielectric = Parameter("reaction field dielectric", 78.3,
                          """Dielectric constant to use if the reaction field cutoff method is used.""")

temperature = Parameter("temperature", 25*celsius, """Simulation temperature""")

pressure = Parameter("pressure", 1*atm, """Simulation pressure""")

topfile = Parameter("topfile", "SYSTEM.top",
                    """Name of the topology file containing the system to be simulated.""")

crdfile = Parameter("crdfile", "SYSTEM.crd",
                    """Name of the coordinate file containing the coordinates of the 
                       system to be simulated.""")

s3file = Parameter("s3file", "SYSTEM.s3",
                    """Name to use for the intermediate s3 file that will contain the 
                       simulation system after it has been loaded from the top/crd files.""")

restart_file = Parameter("restart file", "sim_restart.s3",
                         """Name of the restart file to use to save progress during the simulation.""")

dcd_root = Parameter("dcd root", "traj", """Root of the filename of the output DCD trajectory files.""")

nmoves = Parameter("nmoves", 1000, """Number of Molecular Dynamics moves to perform during the simulation.""")

random_seed = Parameter("random seed", None, """Random number seed. Set this if you
                         want to have reproducible simulations.""")

ncycles = Parameter("ncycles", 1, """The number of MD cycles. The total elapsed time will be nmoves*ncycles*timestep""")

ncycles_per_snap = Parameter("ncycles_per_snap", 1, """Number of cycles between saving snapshots""")

save_coords = Parameter("save coordinates", True, """Whether or not to save coordinates.""")

buffered_coords_freq = Parameter("buffered coordinates frequency", 0, 
                                 """The number of time steps between saving of coordinates during
                                 a cycle of MD. 0 disables buffering.""")

time_to_skip = Parameter("time to skip", 0*picosecond, """Time to skip?""")

minimize = Parameter("minimize", False, """Whether or not to perform minimization before the simulation.""")

minimize_tol = Parameter("minimize tolerance", 1e-8, """Tolerance used to know when minimization is complete.""")

minimize_max_iter = Parameter("minimize maximum iterations", 1000, """Maximum number of iterations for minimization.""")

equilibrate = Parameter("equilibrate", False, """Whether or not to perform equilibration before dynamics.""")

equil_iterations = Parameter("equilibration iterations", 2000, """Number of equilibration steps to perform.""")

equil_timestep = Parameter("equilibration timestep", 0.5*femtosecond, """Timestep to use during equilibration.""")

combining_rules = Parameter("combining rules", "arithmetic", """Combining rules to use for the non-bonded interactions.""")

timestep = Parameter("timestep", 2 * femtosecond, """Timestep for the dynamics simulation.""")

platform = Parameter("platform", "CPU", """Which OpenMM platform should be used to perform the dynamics.""")

precision = Parameter("precision", "mixed", """The floating point precision model to use during dynamics.""")

constraint = Parameter("constraint", "none", """The constraint model to use during dynamics.""")

cutoff_type = Parameter("cutoff type", "cutoffperiodic", """The cutoff method to use during the simulation.""")

cutoff_dist = Parameter("cutoff distance", 10*angstrom, """The cutoff distance to use for the non-bonded interactions.""")

integrator_type = Parameter("integrator", "leapfrogverlet", """The integrator to use for dynamics.""")

inverse_friction = Parameter("inverse friction", 0.1*picosecond, """Inverse friction time for the Langevin thermostat.""")

andersen = Parameter("thermostat", True, """Whether or not to use the Andersen thermostat (needed for NVT or NPT simulation).""")

barostat = Parameter("barostat", False, """Whether or not to use a barostat (needed for NPT simulation).""")

andersen_frequency = Parameter("andersen frequency", 10.0, """Collision frequency in units of (1/ps)""")

barostat_frequency = Parameter("barostat frequency", 25, """Number of steps before attempting box changes if using the barostat.""")

lj_dispersion = Parameter("lj dispersion", False, """Whether or not to calculate and include the LJ dispersion term.""")

cmm_removal = Parameter("center of mass frequency", 10, "Frequency of which the system center of mass motion is removed.""")

center_solute = Parameter("center solute", True, """Whether or not to centre the centre of geometry of the solute in the box.""")

use_restraints = Parameter("use restraints", False, """Whether or not to use harmonic restaints on the solute atoms.""")

k_restraint = Parameter("restraint force constant", 100.0, """Force constant to use for the harmonic restaints.""")

heavy_mass_restraint = Parameter("heavy mass restraint", 1.10, """Only restrain solute atoms whose mass is greater than this value.""")

unrestrained_residues = Parameter("unrestrained residues", ["WAT", "HOH"], """Names of residues that are never restrained.""")

freeze_residues = Parameter("freeze residues", True, """Whether or not to freeze certain residues.""")

frozen_residues = Parameter("frozen residues", ["LGR", "SIT", "NEG", "POS"], """List of residues to freeze if 'freeze residues' is True.""")

#use_distance_restraints = Parameter("use distance restraints",True, """Whether or not to use restraints distances between pairs of atoms.""") 
use_distance_restraints = Parameter("use distance restraints",False, """Whether or not to use restraints distances between pairs of atoms.""")

#distance_restraints_dict = Parameter("distance restraints dictionary",{ (17,12):(3.0,10.0, 0.2) }, """Dictionnary of pair of atoms whose distance is restrained, and restraint parameters. Syntax is {(atom0,atom1):(reql, kl, Dl)} where atom0, atom1 are atomic indices. reql the equilibrium distance. Kl the force constant of the restraint. D the flat bottom radius.""")

distance_restraints_dict = Parameter("distance restraints dictionary",{ }, """Dictionnary of pair of atoms whose distance is restrained, and restraint parameters. Syntax is {(atom0,atom1):(reql, kl, Dl)} where atom0, atom1 are atomic indices. reql the equilibrium distance. Kl the force constant of the restraint. D the flat bottom radius. WARNING: PBC distance checks not implemented, avoid restraining pair of atoms that may diffuse out of the box.""")

## Free energy specific keywords 

morphfile = Parameter("morphfile", "SYSTEM.morph",
                      """Name of the morph file containing the perturbation to apply to the system.""")

lambda_val = Parameter("lambda_val", 0.0, 
                       """Value of the lambda parameter at which to evaluate free energy gradients.""")

delta_lambda = Parameter("delta_lambda", 0.001,
                         """Value of the lambda interval used to evaluate free energy gradients by finite difference.""")

shift_delta = Parameter("shift delta", 2.0,
                        """Value of the Lennard-Jones softcore parameter.""")

coulomb_power = Parameter("coulomb power", 0,
                          """Value of the Coulombic softcore parameter.""")

energy_frequency = Parameter("energy frequency", 10,
                             """The number of time steps between evaluation of free energy gradients.""")


#####################################

def setupDCD(DCD_root, system):

    files = os.listdir(os.getcwd())

    #if (save_coords.val):
    #    buffered_coords_freq = 500
    #else:
    #    buffer_freq = 0
    
    dcds = []
    for f in files:
        if f.endswith(".dcd"):
            dcds.append(f)
    
    dcds.sort()
    
    index = len(dcds) + 1
    
    dcd_filename = dcd_root.val + "%0009d" % index + ".dcd"

    Trajectory = DCDFile(dcd_filename, system[MGName("all")], system.property("space"), timestep.val, interval=buffered_coords_freq.val*ncycles_per_snap.val)

    return Trajectory


def writeSystemData( system, moves, Trajectory, block):

    localtimer = QTime()
    localtimer.start()

    #if (save_coords.val):
    #    buffer_freq = 500
    #else:
    #    buffer_freq = 0

    if (block % ncycles_per_snap.val == 0):
        #PDB().write(system[MGName("all")], "output%0009d.pdb" % block)

        if buffered_coords_freq.val > 0:
            dimensions = {}
            sysprops = system.propertyKeys()
            for prop in sysprops:
                if prop.startswith("buffered_space"):
                    dimensions[str(prop)] = system.property(prop) 
            Trajectory.writeBufferedModels( system[MGName("all")], dimensions )
        else:
            Trajectory.writeModel( system[MGName("all")], system.property("space") )

    FILE = open("moves.dat", "w")
    print("%s" % moves, file=FILE)

    #print(" Time to write coordinates %s ms " % localtimer.elapsed())


def centerSolute(system, space):

    # ! Assuming first molecule in the system is the solute ! 
    
    if space.isPeriodic():
        box_center = space.dimensions()/2
    else:
        box_center = Vector( 0.0, 0.0, 0.0 )

    solute = system.molecules().at(MolNum(1)) #first().molecule()

    solute_cog = CenterOfGeometry(solute).point()

    delta = box_center - solute_cog

    molNums = system.molNums()

    for molnum in molNums:
        mol = system.molecule(molnum).molecule()
        molcoords = mol.property("coordinates")
        molcoords.translate(delta)
        mol = mol.edit().setProperty("coordinates",molcoords).commit()
        system.update(mol)

    return system


def createSystem(molecules):
    #print("Applying flexibility and zmatrix templates...")
    print("Creating the system...")

    moleculeNumbers = molecules.molNums()
    moleculeList = []

    for moleculeNumber in moleculeNumbers:
        molecule = molecules.molecule(moleculeNumber).molecule()
        moleculeList.append(molecule)

    molecules = MoleculeGroup("molecules")
    ions = MoleculeGroup("ions")
    
    for molecule in moleculeList:
        natoms = molecule.nAtoms()
        if natoms == 1:
            ions.add(molecule)
        else:
            molecules.add(molecule)

    all = MoleculeGroup("all")
    all.add(molecules)
    all.add(ions)

    # Add these groups to the System
    system = System()

    system.add(all)
    system.add(molecules)
    system.add(ions)
       
    return system

def setupForcefields(system, space):

    print("Creating force fields... ")

    all = system[ MGName("all") ]
    molecules = system[ MGName("molecules")]
    ions = system[ MGName("ions") ]
    
    # - first solvent-solvent coulomb/LJ (CLJ) energy
    internonbondedff = InterCLJFF("molecules:molecules")
    if (cutoff_type.val != "nocutoff") :
        internonbondedff.setUseReactionField(True)
        internonbondedff.setReactionFieldDielectric(rf_dielectric.val)
    internonbondedff.add(molecules)

    inter_ions_nonbondedff = InterCLJFF("ions:ions")
    if (cutoff_type.val != "nocutoff") :
        inter_ions_nonbondedff.setUseReactionField(True)
        inter_ions_nonbondedff.setReactionFieldDielectric(rf_dielectric.val)

    inter_ions_nonbondedff.add(ions)
    
    inter_ions_molecules_nonbondedff = InterGroupCLJFF("ions:molecules")
    if (cutoff_type.val != "nocutoff") :
        inter_ions_molecules_nonbondedff.setUseReactionField(True)
        inter_ions_molecules_nonbondedff.setReactionFieldDielectric(rf_dielectric.val)

    inter_ions_molecules_nonbondedff.add(ions, MGIdx(0) )
    inter_ions_molecules_nonbondedff.add(molecules, MGIdx(1) )
    
    # Now solute bond, angle, dihedral energy
    intrabondedff = InternalFF("molecules-intrabonded")
    intrabondedff.add(molecules)

    # Now solute intramolecular CLJ energy
    intranonbondedff = IntraCLJFF("molecules-intranonbonded")

    if (cutoff_type.val != "nocutoff") :
        intranonbondedff.setUseReactionField(True)
        intranonbondedff.setReactionFieldDielectric(rf_dielectric.val)

    intranonbondedff.add(molecules)

    # solute restraint energy
    # 
    # We restrain atoms based ont he contents of the property "restrainedatoms"
    #
    restraintff = RestraintFF("restraint")

    if use_restraints.val:
        molnums = molecules.molecules().molNums()

        for molnum in molnums:
            mol = molecules.molecule(molnum).molecule()
            try:
                mol_restrained_atoms = propertyToAtomNumVectorList(mol.property("restrainedatoms"))
            except UserWarning as error:
                error_type = re.search(r"(Sire\w*::\w*)", str(error)).group(0)
                if error_type == "SireBase::missing_property":
                    continue
                else:
                    raise error

            for restrained_line in mol_restrained_atoms:
                atnum = restrained_line[0]
                restraint_atom = mol.select(atnum)
                restraint_coords = restrained_line[1]
                restraint_k = restrained_line[2] * kcal_per_mol / (angstrom*angstrom)
            
                restraint = DistanceRestraint.harmonic(restraint_atom, restraint_coords, restraint_k)

                restraintff.add(restraint)

    # Here is the list of all forcefields
    forcefields = [ internonbondedff, intrabondedff, intranonbondedff, inter_ions_nonbondedff, inter_ions_molecules_nonbondedff, restraintff ] 
    
    for forcefield in forcefields:
        system.add(forcefield)

    system.setProperty( "space", space )
    system.setProperty("switchingFunction", CHARMMSwitchingFunction(cutoff_dist.val) )
    system.setProperty( "combiningRules", VariantProperty(combining_rules.val) )
    #system.setProperty( "useReactionField", VariantProperty(True) )
    #system.setProperty( "reactionFieldDielectric", VariantProperty(rf_dielectric.val) )

    total_nrg = internonbondedff.components().total() +\
                intranonbondedff.components().total() + intrabondedff.components().total() +\
                inter_ions_nonbondedff.components().total() + inter_ions_molecules_nonbondedff.components().total() +\
                restraintff.components().total()

    e_total = system.totalComponent()

    system.setComponent( e_total, total_nrg )

    # Add a monitor that calculates the average total energy and average energy
    # deltas - we will collect both a mean average and an zwanzig average
    system.add( "total_energy", MonitorComponent(e_total, Average()) )

    return system

def setupMoves(system, random_seed, GPUS):

    print("Setting up moves...")


    molecules = system[ MGName("all") ]

    Integrator_OpenMM = OpenMMMDIntegrator( molecules )
  
    Integrator_OpenMM.setPlatform(platform.val)
    Integrator_OpenMM.setConstraintType(constraint.val)
    Integrator_OpenMM.setCutoffType(cutoff_type.val)
    Integrator_OpenMM.setIntegrator(integrator_type.val)
    Integrator_OpenMM.setFriction( inverse_friction.val )# Only meaningful for Langevin/Brownian integrators
    Integrator_OpenMM.setPrecision(precision.val)
    Integrator_OpenMM.setTimetoSkip(time_to_skip.val)
    Integrator_OpenMM.setMinimization(minimize.val)
    Integrator_OpenMM.setMinimizeTol(minimize_tol.val)
    Integrator_OpenMM.setMinimizeIterations(minimize_max_iter.val)

    if equilibrate.val:
        Integrator_OpenMM.setEquilib_iterations(equil_iterations.val)
    else:
        Integrator_OpenMM.setEquilib_iterations(0)

    Integrator_OpenMM.setEquilib_time_step(equil_timestep.val)
    Integrator_OpenMM.setDeviceIndex(str(GPUS))
    Integrator_OpenMM.setLJDispersion(lj_dispersion.val)

    if cutoff_type.val != "nocutoff":
        Integrator_OpenMM.setCutoff_distance(cutoff_dist.val)
    if cutoff_type.val == "cutoffperiodic":
        Integrator_OpenMM.setField_dielectric(rf_dielectric.val)

    Integrator_OpenMM.setCMMremoval_frequency(cmm_removal.val)

    #if (save_coords.val):
    #    buffer_freq = 500
    #else:
    #    buffer_freq = 0

    Integrator_OpenMM.setBufferFrequency(buffered_coords_freq.val)

    if use_restraints.val:
        Integrator_OpenMM.setRestraint(True)

    if andersen.val:
        Integrator_OpenMM.setTemperature(temperature.val)
        Integrator_OpenMM.setAndersen(andersen.val)
        Integrator_OpenMM.setAndersen_frequency(andersen_frequency.val)

    if barostat.val:
        Integrator_OpenMM.setPressure(pressure.val)
        Integrator_OpenMM.setMCBarostat(barostat.val)
        Integrator_OpenMM.setMCBarostat_frequency(barostat_frequency.val)

    #print Integrator_OpenMM.getDeviceIndex()
    Integrator_OpenMM.initialise()

    mdmove = MolecularDynamics(molecules, Integrator_OpenMM, timestep.val, {"velocity generator":MaxwellBoltzmann(temperature.val)})

    print("Created a MD move that uses OpenMM for all molecules on %s " % GPUS)

    moves = WeightedMoves()
    moves.add(mdmove, 1)
    
    if (not random_seed):
        random_seed = RanGenerator().randInt(100000,1000000)
    print("Generated random seed number %d " % random_seed)

    moves.setGenerator( RanGenerator(random_seed) )

    return moves

def atomNumListToProperty( list ):
    prop = Properties()

    i = 0

    for value in list:
        prop.setProperty(str(i), VariantProperty(value.value()))
        i += 1

    return prop

def atomNumVectorListToProperty( list ):
    prop = Properties()

    i = 0
                
    for value in list:
        prop.setProperty("AtomNum(%d)" % i, VariantProperty(value[0].value()))
        prop.setProperty("x(%d)" % i, VariantProperty(value[1].x()))
        prop.setProperty("y(%d)" % i, VariantProperty(value[1].y()))
        prop.setProperty("z(%d)" % i, VariantProperty(value[1].z()))
        prop.setProperty("k(%d)" % i, VariantProperty(value[2].val ) )
        i += 1
    
    prop.setProperty("nrestrainedatoms", VariantProperty(i) );

    return prop

def linkbondVectorListToProperty( list ):

    prop = Properties()

    i = 0

    for value in list:
        prop.setProperty("AtomNum0(%d)" % i, VariantProperty(value[0]))
        prop.setProperty("AtomNum1(%d)" % i, VariantProperty(value[1]))
        prop.setProperty("reql(%d)" % i, VariantProperty(value[2]))
        prop.setProperty("kl(%d)" % i, VariantProperty(value[3]))
        prop.setProperty("dl(%d)" % i, VariantProperty(value[4]))
        i += 1

    prop.setProperty("nbondlinks", VariantProperty(i) );

    return prop


def propertyToAtomNumList( prop ):
    list = []

    i = 0

    try:
        while True:
            list.append( AtomNum(prop[str(i)].toInt()) )
            i += 1
    except:
        pass

    return list

def propertyToAtomNumVectorList( prop ):
    list = []
            
    i = 0
        
    try:
        while True:
            num = AtomNum(prop["AtomNum(%d)" % i].toInt())
            x = prop["x(%d)" % i].toDouble()
            y = prop["y(%d)" % i].toDouble()
            z = prop["z(%d)" % i].toDouble()
            k = prop["k(%d)" % i].toDouble()

            list.append( (num, Vector(x,y,z), k ) )
 
            i += 1
    except:
        pass

    return list


def setupRestraints(system):

    molecules = system[ MGName("all") ].molecules()

    molnums = molecules.molNums()

    for molnum in molnums:
        mol = molecules.molecule(molnum).molecule()
        nats = mol.nAtoms()
        atoms = mol.atoms()
        
        restrainedAtoms = []

        #
        # This will apply a restraint to every atom that is 
        # A) NOT a hydrogen
        # B) NOT in an unrestrained residue.
        #
        for x in range(0,nats):
            at = atoms[x]
            atnumber = at.number()
            #print at, atnumber
            if at.residue().name().value() in unrestrained_residues.val:
                continue
            #print at, at.property("mass"), heavyMass
            if ( at.property("mass").value() < heavy_mass_restraint.val ):
                #print "LIGHT, skip"
                continue
            atcoords = at.property("coordinates")
            #print at
            restrainedAtoms.append( ( atnumber , atcoords, k_restraint) )
            #restrainedAtoms.append( atnumber )
        
        if len(restrainedAtoms) > 0:
            mol = mol.edit().setProperty("restrainedatoms", atomNumVectorListToProperty(restrainedAtoms)).commit()
            #print restrainedAtoms
            #print propertyToAtomNumVectorList( mol.property("restrainedatoms") )
            system.update(mol)

    return system

def setupDistanceRestraints(system):
    prop_list = []

    molecules = system[ MGName("all") ].molecules()

    dic_items = list( distance_restraints_dict.val.items() )

    for i in range(0,molecules.nMolecules()):
        mol = molecules.molecule(MolNum(i+1)).molecule()
        atoms_mol = mol.atoms()
        natoms_mol = mol.nAtoms()
        for j in range(0,natoms_mol):
            at = atoms_mol[j]
            atnumber = at.number()
            for k in range(len(dic_items)):
                if dic_items[k][0][0] == dic_items[k][0][1]:
                    print ("Error! It is not possible to place a distance restraint on the same atom")
                    sys.exit(-1)
                if atnumber.value() - 1 in dic_items[k][0]:
                    print (at)
                    # atom0index atom1index, reql, kl, dl 
                    prop_list.append((dic_items[k][0][0]+1, dic_items[k][0][1]+1,dic_items[k][1][0],dic_items[k][1][1], dic_items[k][1][2]))

    unique_prop_list = []

    [unique_prop_list.append(item) for item in prop_list if item not in unique_prop_list]

    print (unique_prop_list)

    #Mol number 0 will store all the information related to the bond-links in the system
    mol0 = molecules.molecule(MolNum(1)).molecule()

    mol0 = mol0.edit().setProperty("linkbonds", linkbondVectorListToProperty( unique_prop_list )).commit()

    system.update(mol0)

    return system



def freezeResidues(system):
   
    molecules = system[ MGName("all") ].molecules()

    molnums = molecules.molNums()

    for molnum in molnums:
        mol = molecules.molecule(molnum).molecule()
        nats = mol.nAtoms()
        atoms = mol.atoms()

        for x in range(0,nats):
            at = atoms[x]
            atnumber = at.number()
            if at.residue().name().value() in frozen_residues.val:
                print("Freezing %s %s %s " % (at, atnumber, at.residue().name().value() ))
                mol = at.edit().setProperty("mass", 0 * g_per_mol).molecule()
                    
        system.update(mol)

    return system

def getDummies(molecule):
    print ("Selecting dummy groups")
    natoms = molecule.nAtoms()
    atoms = molecule.atoms()

    from_dummies = None
    to_dummies = None

    for x in range(0,natoms):
        atom = atoms[x]
        if atom.property("initial_ambertype") == "du":
            if from_dummies is None:
                from_dummies = molecule.selectAll( atom.index() )
            else:
                from_dummies += molecule.selectAll( atom.index() )
        elif atom.property("final_ambertype") == "du":
            if to_dummies is None:
                to_dummies = molecule.selectAll( atom.index() )
            else:
                to_dummies += molecule.selectAll( atom.index() )
    
    return to_dummies, from_dummies

def createSystemFreeEnergy(molecules):
    print ("Create the System...")

    moleculeNumbers = molecules.molNums()
    moleculeList = []

    for moleculeNumber in moleculeNumbers:
        molecule = molecules.molecule(moleculeNumber).molecule()
        moleculeList.append(molecule)

    #
    # The code below assumes that the solute to be perturbed is 
    # the first molecule in the top file.
    # The residue name of the first residue in this molecule is 
    # used to name the solute. This is used later to match 
    # templates in the flex/pert files.

    solute = moleculeList[0]
    lig_name = solute.residue( ResIdx(0) ).name().value()
        
    solute = solute.edit().rename(lig_name).commit()
                    
    perturbations_lib = PerturbationsLibrary(morphfile.val)
    solute = perturbations_lib.applyTemplate(solute)

    perturbations = solute.property("perturbations")

    lam = Symbol("lambda")

    initial = Perturbation.symbols().initial()
    final = Perturbation.symbols().final()

    solute = solute.edit().setProperty("perturbations",
                perturbations.recreate( (1-lam)*initial + lam*final ) ).commit()

    # We put atoms in three groups depending on what happens in the perturbation
    # non dummy to non dummy --> the hard group, uses a normal intermolecular FF
    # non dummy to dummy --> the todummy group, uses SoftFF with alpha = Lambda
    # dummy to non dummy --> the fromdummy group, uses SoftFF with alpha = 1 - Lambda
    # We start assuming all atoms are hard atoms. Then we call getDummies to find which atoms 
    # start/end as dummies and update the hard, todummy and fromdummy groups accordingly

    solute_grp_ref = MoleculeGroup("solute_ref", solute)
    solute_grp_ref_hard = MoleculeGroup("solute_ref_hard")
    solute_grp_ref_todummy = MoleculeGroup("solute_ref_todummy")
    solute_grp_ref_fromdummy = MoleculeGroup("solute_ref_fromdummy")

    solute_ref_hard = solute.selectAllAtoms()
    solute_ref_todummy = solute_ref_hard.invert()
    solute_ref_fromdummy = solute_ref_hard.invert()

    to_dummies, from_dummies = getDummies(solute)

    if to_dummies is not None:
        ndummies = to_dummies.count()
        dummies = to_dummies.atoms()

        for x in range(0,ndummies):
            dummy_index = dummies[x].index()
            solute_ref_hard = solute_ref_hard.subtract( solute.select( dummy_index ) )
            solute_ref_todummy = solute_ref_todummy.add( solute.select( dummy_index ) )

    if from_dummies is not None:
        ndummies = from_dummies.count()
        dummies = from_dummies.atoms()

        for x in range(0,ndummies):
            dummy_index = dummies[x].index()
            solute_ref_hard = solute_ref_hard.subtract( solute.select( dummy_index ) )
            solute_ref_fromdummy = solute_ref_fromdummy.add( solute.select( dummy_index ) )

    solute_grp_ref_hard.add(solute_ref_hard)
    solute_grp_ref_todummy.add(solute_ref_todummy)
    solute_grp_ref_fromdummy.add(solute_ref_fromdummy)

    solutes = MoleculeGroup("solutes")
    solutes.add(solute)

    molecules = MoleculeGroup("molecules")
    molecules.add(solute)

    solvent = MoleculeGroup("solvent")

    for molecule in moleculeList[1:]:
        molecules.add(molecule)
        solvent.add(molecule)

    all = MoleculeGroup("all")

    all.add(molecules)
    all.add(solvent)

    all.add(solutes)
    all.add(solute_grp_ref)
    all.add(solute_grp_ref_hard)
    all.add(solute_grp_ref_todummy)
    all.add(solute_grp_ref_fromdummy)

    # Add these groups to the System
    system = System()

    system.add(solutes)
    system.add(solute_grp_ref)
    system.add(solute_grp_ref_hard)
    system.add(solute_grp_ref_todummy)
    system.add(solute_grp_ref_fromdummy)

    system.add(molecules)

    system.add(solvent)

    system.add(all)

    return system

def setupForcefieldsFreeEnergy(system, space ):

    print ("Creating force fields... ")

    solutes = system[ MGName("solutes") ]

    solute = system[ MGName("solute_ref") ]
    solute_hard = system[ MGName("solute_ref_hard") ]
    solute_todummy = system[ MGName("solute_ref_todummy") ]
    solute_fromdummy = system[ MGName("solute_ref_fromdummy") ]

    solvent = system [ MGName("solvent") ]

    all = system[ MGName("all") ]

    # ''solvent'' is actually every molecule that isn't perturbed ! 
    solvent_intraff = InternalFF("solvent_intraff")
    solvent_intraff.add(solvent)

    # Solute bond, angle, dihedral energy
    solute_intraff = InternalFF("solute_intraff")
    solute_intraff.add(solute)

    # Solvent-solvent coulomb/LJ (CLJ) energy
    solventff = InterCLJFF("solvent:solvent")
    if (cutoff_type.val != "nocutoff") :
        solventff.setUseReactionField(True)
        solventff.setReactionFieldDielectric(rf_dielectric.val)
    solventff.add(solvent)

    #Solvent intramolecular CLJ energy
    solvent_intraclj = IntraCLJFF("solvent_intraclj")
    if (cutoff_type.val != "nocutoff") :
        solvent_intraclj.setUseReactionField(True)
        solvent_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    solvent_intraclj.add(solvent)

    # Solute intramolecular CLJ energy
    solute_hard_intraclj = IntraCLJFF("solute_hard_intraclj")
    if (cutoff_type.val != "nocutoff") :
        solute_hard_intraclj.setUseReactionField(True)
        solute_hard_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    solute_hard_intraclj.add(solute_hard)

    solute_todummy_intraclj = IntraSoftCLJFF("solute_todummy_intraclj")
    solute_todummy_intraclj.setShiftDelta(shift_delta.val)
    solute_todummy_intraclj.setCoulombPower(coulomb_power.val)
    if (cutoff_type.val != "nocutoff") :
        solute_todummy_intraclj.setUseReactionField(True)
        solute_todummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    solute_todummy_intraclj.add(solute_todummy)

    solute_fromdummy_intraclj = IntraSoftCLJFF("solute_fromdummy_intraclj")
    solute_fromdummy_intraclj.setShiftDelta(shift_delta.val)
    solute_fromdummy_intraclj.setCoulombPower(coulomb_power.val)
    if (cutoff_type.val != "nocutoff") :
        solute_fromdummy_intraclj.setUseReactionField(True)
        solute_fromdummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    solute_fromdummy_intraclj.add(solute_fromdummy)

    solute_hard_todummy_intraclj = IntraGroupSoftCLJFF("solute_hard:todummy_intraclj")
    solute_hard_todummy_intraclj.setShiftDelta(shift_delta.val)
    solute_hard_todummy_intraclj.setCoulombPower(coulomb_power.val)
    if (cutoff_type.val != "nocutoff") :
        solute_hard_todummy_intraclj.setUseReactionField(True)
        solute_hard_todummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    solute_hard_todummy_intraclj.add(solute_hard, MGIdx(0))
    solute_hard_todummy_intraclj.add(solute_todummy, MGIdx(1))

    solute_hard_fromdummy_intraclj = IntraGroupSoftCLJFF("solute_hard:fromdummy_intraclj")
    solute_hard_fromdummy_intraclj.setShiftDelta(shift_delta.val)
    solute_hard_fromdummy_intraclj.setCoulombPower(coulomb_power.val)
    if (cutoff_type.val != "nocutoff") :
        solute_hard_fromdummy_intraclj.setUseReactionField(True)
        solute_hard_fromdummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    solute_hard_fromdummy_intraclj.add(solute_hard, MGIdx(0))
    solute_hard_fromdummy_intraclj.add(solute_fromdummy, MGIdx(1))

    solute_todummy_fromdummy_intraclj = IntraGroupSoftCLJFF("solute_todummy:fromdummy_intraclj")
    solute_todummy_fromdummy_intraclj.setShiftDelta(shift_delta.val)
    solute_todummy_fromdummy_intraclj.setCoulombPower(coulomb_power.val)
    if (cutoff_type.val != "nocutoff") :
        solute_todummy_fromdummy_intraclj.setUseReactionField(True)
        solute_todummy_fromdummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    solute_todummy_fromdummy_intraclj.add(solute_todummy, MGIdx(0))
    solute_todummy_fromdummy_intraclj.add(solute_fromdummy, MGIdx(1))

    #Solute-solvent CLJ energy
    solute_hard_solventff = InterGroupCLJFF("solute_hard:solvent")
    if (cutoff_type.val != "nocutoff") :
        solute_hard_solventff.setUseReactionField(True)
        solute_hard_solventff.setReactionFieldDielectric(rf_dielectric.val)
    solute_hard_solventff.add(solute_hard, MGIdx(0))
    solute_hard_solventff.add(solvent, MGIdx(1))

    solute_todummy_solventff = InterGroupSoftCLJFF("solute_todummy:solvent")
    if (cutoff_type.val != "nocutoff") :
        solute_todummy_solventff.setUseReactionField(True)
        solute_todummy_solventff.setReactionFieldDielectric(rf_dielectric.val)
    solute_todummy_solventff.add(solute_todummy, MGIdx(0))
    solute_todummy_solventff.add(solvent, MGIdx(1))


    solute_fromdummy_solventff = InterGroupSoftCLJFF("solute_fromdummy:solvent")
    if (cutoff_type.val != "nocutoff") :
        solute_fromdummy_solventff.setUseReactionField(True)
        solute_fromdummy_solventff.setReactionFieldDielectric(rf_dielectric.val)
    solute_fromdummy_solventff.add(solute_fromdummy, MGIdx(0))
    solute_fromdummy_solventff.add(solvent, MGIdx(1))


    # TOTAL
    forcefields =  [ solute_intraff,
                  solute_hard_intraclj, solute_todummy_intraclj, solute_fromdummy_intraclj,
                  solute_hard_todummy_intraclj, solute_hard_fromdummy_intraclj, 
                  solute_todummy_fromdummy_intraclj, 
                  solvent_intraff,
                  solventff,solvent_intraclj,
                  solute_hard_solventff, solute_todummy_solventff, solute_fromdummy_solventff ]

    # BONDED
    #   forcefields =  [ solute_intraff, solvent_intraff ]

    # NON BONDED
    #    forcefields =  [solute_hard_intraclj, solute_todummy_intraclj, solute_fromdummy_intraclj,
    #                 solute_hard_todummy_intraclj, solute_hard_fromdummy_intraclj, 
    #                 solute_todummy_fromdummy_intraclj, 
    #                 solventff, solvent_intraclj,
    #                 solute_hard_solventff, solute_todummy_solventff, solute_fromdummy_solventff ]

    for forcefield in forcefields:
        system.add(forcefield)

    system.setProperty( "space", space )

    if (cutoff_type.val != "nocutoff") :
        system.setProperty( "switchingFunction", CHARMMSwitchingFunction(cutoff_dist.val) )
    else :
        system.setProperty( "switchingFunction", NoCutoff())

    system.setProperty( "combiningRules", VariantProperty(combining_rules.val) )
    system.setProperty( "coulombPower", VariantProperty(coulomb_power.val) )
    system.setProperty( "shiftDelta", VariantProperty(shift_delta.val) )
    #system.setProperty( "useReactionField", VariantProperty(True) )
    #system.setProperty( "reactionFieldDielectric", VariantProperty(rf_dielectric.val) )
    
    # TOTAL
    total_nrg = solute_intraff.components().total() + solute_hard_intraclj.components().total() +\
        solute_todummy_intraclj.components().total(0) + solute_fromdummy_intraclj.components().total(0) +\
        solute_hard_todummy_intraclj.components().total(0) + solute_hard_fromdummy_intraclj.components().total(0) +\
        solute_todummy_fromdummy_intraclj.components().total(0) +\
        solvent_intraff.components().total() + solventff.components().total() +\
        solvent_intraclj.components().total() +\
        solute_hard_solventff.components().total() +\
        solute_todummy_solventff.components().total(0) +\
        solute_fromdummy_solventff.components().total(0)

    # BONDED
    #   total_nrg = solute_intraff.components().total() + solvent_intraff.components().total()


    # NON BONDED
    #   total_nrg = solute_hard_intraclj.components().total() +\
    #     solute_todummy_intraclj.components().total(0) + solute_fromdummy_intraclj.components().total(0) +\
    #     solute_hard_todummy_intraclj.components().total(0) + solute_hard_fromdummy_intraclj.components().total(0) +\
    #     solute_todummy_fromdummy_intraclj.components().total(0) +\
    #     solventff.components().total() +\
    #     solute_hard_solventff.components().total() +\
    #     solute_todummy_solventff.components().total(0) +\
    #     solute_fromdummy_solventff.components().total(0)
    

    e_total = system.totalComponent()

    lam = Symbol("lambda")

    system.setComponent( e_total, total_nrg )

    system.setConstant(lam, 0.0)

    system.add( PerturbationConstraint(solutes) )

    # NON BONDED Alpha constraints for the soft force fields

    system.add( PropertyConstraint( "alpha0", FFName("solute_todummy_intraclj"), lam ) )
    system.add( PropertyConstraint( "alpha0", FFName("solute_fromdummy_intraclj"), 1 - lam ) )
    system.add( PropertyConstraint( "alpha0", FFName("solute_hard:todummy_intraclj"), lam ) )
    system.add( PropertyConstraint( "alpha0", FFName("solute_hard:fromdummy_intraclj"), 1 - lam ) )
    system.add( PropertyConstraint( "alpha0", FFName("solute_todummy:fromdummy_intraclj"), Max( lam, 1 - lam )  ) )
    system.add( PropertyConstraint( "alpha0", FFName("solute_todummy:solvent"), lam ) )
    system.add( PropertyConstraint( "alpha0", FFName("solute_fromdummy:solvent"), 1 - lam ) )

    system.setComponent( lam, lambda_val.val )

    # printEnergies( system.componentValues() )

    return system

def setupMovesFreeEnergy(system,random_seed,GPUS,lam_val):

    print ("Setting up moves...")

    molecules = system[ MGName("molecules") ]
    solute = system[ MGName("solute_ref") ]
    solute_hard = system[ MGName("solute_ref_hard") ]
    solute_todummy = system[ MGName("solute_ref_todummy") ]
    solute_fromdummy = system[ MGName("solute_ref_fromdummy") ]

    Integrator_OpenMM = OpenMMFrEnergyST(molecules,solute,solute_hard,solute_todummy,solute_fromdummy)

    Integrator_OpenMM.setIntegrator(integrator_type.val)
    Integrator_OpenMM.setFriction(inverse_friction.val )# Only meaningful for Langevin/Brownian integrators
    Integrator_OpenMM.setPlatform(platform.val)
    Integrator_OpenMM.setConstraintType(constraint.val)
    Integrator_OpenMM.setCutoffType(cutoff_type.val)
    Integrator_OpenMM.setField_dielectric(rf_dielectric.val)
    Integrator_OpenMM.setAlchemical_value(lambda_val.val)
    Integrator_OpenMM.setDeviceIndex(str(GPUS))
    Integrator_OpenMM.setCoulomb_power(coulomb_power.val)
    Integrator_OpenMM.setShift_delta(shift_delta.val)
    Integrator_OpenMM.setDeltatAlchemical(delta_lambda.val)
    Integrator_OpenMM.setPrecision(precision.val)
    Integrator_OpenMM.setTimetoSkip(time_to_skip.val)
    Integrator_OpenMM.setMinimization(minimize.val)
    Integrator_OpenMM.setMinimizeTol(minimize_tol.val)
    Integrator_OpenMM.setMinimizeIterations(minimize_max_iter.val)
 
    
    if equilibrate.val:
        Integrator_OpenMM.setEquilib_iterations(equil_iterations.val)
    else:
        Integrator_OpenMM.setEquilib_iterations(0)
  
    Integrator_OpenMM.setEquilib_time_step(equil_timestep.val)

    Integrator_OpenMM.setBufferFrequency(buffered_coords_freq.val)

    if cutoff_type != "nocutoff":
        Integrator_OpenMM.setCutoff_distance(cutoff_dist.val)

    Integrator_OpenMM.setCMMremoval_frequency(cmm_removal.val)
    
    Integrator_OpenMM.setEnergyFrequency(energy_frequency.val)


    if use_restraints.val:
        Integrator_OpenMM.setRestraint(True)

    if andersen.val:
        Integrator_OpenMM.setTemperature(temperature.val)
        Integrator_OpenMM.setAndersen(andersen.val)
        Integrator_OpenMM.setAndersen_frequency(andersen_frequency.val)

    if barostat.val:
        Integrator_OpenMM.setPressure(pressure.val)
        Integrator_OpenMM.setMCBarostat(barostat.val)
        Integrator_OpenMM.setMCBarostat_frequency(barostat_frequency.val)


    Integrator_OpenMM.initialise()

    mdmove = MolecularDynamics(molecules, Integrator_OpenMM, timestep.val , {"velocity generator":MaxwellBoltzmann(temperature.val)})

    #mdmove = MolecularDynamics(molecules, Integrator_OpenMM, time_step)

    print("Created a MD move that uses OpenMM for all molecules on %s " % GPUS)

    moves = WeightedMoves()
    moves.add(mdmove, 1)
    
    if (not random_seed):
        random_seed = RanGenerator().randInt(100000,1000000)

    print("Generated random seed number %d " % random_seed)

    moves.setGenerator( RanGenerator(random_seed) )
 
    return moves

def clearBuffers( system ):

    print ("Clearing buffers...")

    mols = system[MGName("all")].molecules()
    molnums = mols.molNums()

    changedmols = MoleculeGroup("changedmols")

    for molnum in molnums:
        mol = mols.molecule(molnum).molecule()
        molprops = mol.propertyKeys()
        editmol = mol.edit()
        for molprop in molprops:
            if molprop.startswith("buffered_"):
                #print "Removing property %s " % molprop
                editmol.removeProperty( PropertyName(molprop) )
        mol = editmol.commit()
        changedmols.add(mol)
        #system.update(mol)

    system.update(changedmols)

    return system



######## MAIN SCRIPTS  #############

@resolveParameters
def run():

    try:
        host = os.environ['HOSTNAME']
    except KeyError:
        host = "unknown"

    print(" ### Running Molecular Dynamics on %s ### " % host)

    timer = QTime()
    timer.start()

    # Setup the system from scratch if no restart file is available

    if not os.path.exists(restart_file.val):

        print("New run. Loading input and creating restart")
        
        amber = Amber()
        
        if os.path.exists(s3file.val):
            (molecules, space) = Sire.Stream.load(s3file.val)
        else:
            (molecules, space) = amber.readCrdTop(crdfile.val, topfile.val)
            Sire.Stream.save( (molecules,space), s3file.val )

        system = createSystem(molecules)

        if (center_solute.val):
            system = centerSolute(system, space)

        if use_restraints.val:
            system = setupRestraints(system)
        
        # Note that this just set the mass to zero which freezes residues in OpenMM but Sire doesn't known that
        if freeze_residues.val:
            system = freezeResidues(system)

        system = setupForcefields(system, space)

        if random_seed.val:
            ranseed = random_seed.val
        else:
            ranseed = RanGenerator().randInt(100000,1000000)

        print("Setting up the simulation with random seed %s" % ranseed)

        moves = setupMoves(system, ranseed, gpu.val)

        print("Saving restart")
        Sire.Stream.save( [system, moves], restart_file.val )
    else:
        system, moves = Sire.Stream.load( restart_file.val )
        move0 =  moves.moves()[0]
        integrator = move0.integrator()
        integrator.setDeviceIndex(str(gpu.val))
        move0.setIntegrator(integrator)
        moves = WeightedMoves()
        moves.add(move0)
        print("Index GPU = %s " % moves.moves()[0].integrator().getDeviceIndex())
        print("Loaded a restart file on wich we have performed %d moves." % moves.nMoves())

    cycle_start = int(moves.nMoves() / nmoves.val)  + 1
    cycle_end = cycle_start + ncycles.val 

    if (save_coords.val):
        trajectory = setupDCD(dcd_root.val, system)

    s1 = timer.elapsed()/1000.

    # Jm debug 14/10/14
    PDB().write(system[MGName("all")], "frame-0.pdb" )
    
    print("Running MD simulation ")

    for i in range(cycle_start,cycle_end):
        print("\nCycle = ",i,"\n")

        #print("Energy before = %s kJ mol-1" % (system.energy().to(kJ_per_mol)))
        # import ipdb; ipdb.set_trace()
        system = moves.move(system, nmoves.val, True)
        #print("Energy after = %s kJ mol-1" % (system.energy().to(kJ_per_mol)))

        if (save_coords.val):
            writeSystemData(system, moves, trajectory, i)

    s2 = timer.elapsed()/1000.
    print("Simulation took %d s " % ( s2 - s1))

    print("Saving restart")
    Sire.Stream.save( [system, moves], restart_file.val )


@resolveParameters
def runFreeNrg():

    #if (save_coords.val):
    #    buffer_freq = 500
    #else:
    #    buffer_freq = 0

    try:
        host = os.environ['HOSTNAME']
    except KeyError:
        host = "unknown"

    print(" ### Running Single Topology Molecular Dynamics Free Energy on %s ### " % host)
    print ("### Simulation Parameters ### ") 
    Parameter.printAll()
    print ("### ###")

    timer = QTime()
    timer.start()

    # Setup the system from scratch if no restart file is available

    if not os.path.exists(restart_file.val):

        print("New run. Loading input and creating restart")
                
        print("lambda is %s" % lambda_val.val)

        amber = Amber()
        
        if os.path.exists(s3file.val):
            (molecules, space) = Sire.Stream.load(s3file.val)
        else:
            (molecules, space) = amber.readCrdTop(crdfile.val, topfile.val)
            Sire.Stream.save( (molecules,space), s3file.val )

        system = createSystemFreeEnergy(molecules)

        if (center_solute.val):
            system = centerSolute(system, space)

        if use_restraints.val:
            system = setupRestraints(system)
        
        if use_distance_restraints.val:
            system = setupDistanceRestraints(system)

        # Note that this just set the mass to zero which freezes residues in OpenMM but Sire doesn't known that
        if freeze_residues.val:
            system = freezeResidues(system)

        system = setupForcefieldsFreeEnergy(system, space)

        if random_seed.val:
            ranseed = random_seed.val
        else:
            ranseed = RanGenerator().randInt(100000,1000000)

        print("Setting up the simulation with random seed %s" % ranseed)

        moves = setupMovesFreeEnergy(system, ranseed, gpu.val, lambda_val.val)

        print("Saving restart")
        Sire.Stream.save( [system, moves], restart_file.val )
    else:
        system, moves = Sire.Stream.load( restart_file.val )
        move0 =  moves.moves()[0]
        integrator = move0.integrator()
        integrator.setDeviceIndex(str(gpu.val))
        move0.setIntegrator(integrator)
        moves = WeightedMoves()
        moves.add(move0)
        print("Index GPU = %s " % moves.moves()[0].integrator().getDeviceIndex())
        print("Loaded a restart file on wich we have performed %d moves." % moves.nMoves())

    cycle_start = int(moves.nMoves() / nmoves.val)  + 1
    cycle_end = cycle_start + ncycles.val 
   
    lam_str = "%7.5f" % lambda_val.val
    outgradients = open("gradients.dat","a", 1)
    outgradients.write("# lambba_val.val %s\n" % lam_str)

    if (save_coords.val):
        trajectory = setupDCD(dcd_root.val, system)

    s1 = timer.elapsed()/1000.
    
    print("Running MD simulation ")

    grads = {}
    grads[lambda_val.val] = AverageAndStddev()
    for i in range(cycle_start,cycle_end):
        print("\nCycle = ",i,"\n")

        #print("Energy before = %s kJ mol-1" % (system.energy().to(kJ_per_mol)))
        # import ipdb; ipdb.set_trace()
        system = moves.move(system, nmoves.val, True)
        #print("Energy after = %s kJ mol-1" % (system.energy().to(kJ_per_mol)))

        if (save_coords.val):
            writeSystemData(system, moves, trajectory, i)

        mdmoves = moves.moves()[0]
        integrator = mdmoves.integrator()
        gradients = integrator.getGradients()
        outgradients.write("%5d %20.10f\n" % (i, gradients[i-1]))
        grads[lambda_val.val].accumulate( gradients[i-1] )

    s2 = timer.elapsed()/1000.
    print("Simulation took %d s " % ( s2 - s1))

    if os.path.exists("gradients.s3"):
        siregrads = Sire.Stream.load("gradients.s3")
    else:
        siregrads = Gradients()
    siregrads = siregrads + Gradients(grads) 
   
    Sire.Stream.save(siregrads, "gradients.s3")

    if buffered_coords_freq.val > 0:
        system = clearBuffers(system)
        # Necessary to write correct restart
        system.mustNowRecalculateFromScratch()

    print("Backing up previous restart")
    cmd = "cp %s %s.previous" % (restart_file.val, restart_file.val)
    os.system(cmd)
    print ("Saving new restart")
    Sire.Stream.save( [system, moves], restart_file.val )


