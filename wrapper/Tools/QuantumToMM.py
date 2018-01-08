#!/bin/env python
# -*- coding: utf-8 -*-

from Sire.IO import *
from Sire.System import *
from Sire.Mol import *
from Sire.MM import *
from Sire.FF import *
from Sire.CAS import *
from Sire.Maths import *
from Sire.Analysis import *
from Sire.System import *
from Sire.Base import *
from Sire.Units import *
from Sire.Squire import *

import Sire.Config
import Sire.Stream

from Sire.Tools.AmberLoader import *
from Sire.Tools import Parameter, resolveParameters

import os
import shutil
import copy

########################################################
# ALL OF THE GLOBAL USER-AVAILABLE QuanToMM PARAMETERS #
########################################################

cutoff_method = Parameter("cutoff method", "shift electrostatics",
                          """Method used to apply the non-bonded electrostatic cutoff.""")

rf_dielectric = Parameter("reaction field dielectric", 78.3,
                          """Dielectric constant to use if the reaction field cutoff method is used.""")

coul_cutoff = Parameter("coulomb cutoff", 15*angstrom,
                        """Coulomb cutoff length""")

lj_cutoff = Parameter("LJ cutoff", 15*angstrom,
                      """Lennard Jones cutoff length""")

grid_spacing = Parameter("grid spacing", 1.0*angstrom,
                         """Grid spacing used for the grid-based forcefields""")
grid_buffer = Parameter("grid buffer", 2*angstrom,
                        """Buffer around the grid used to prevent recalculation
                           in the grid-based forcefields.""")

disable_grid = Parameter("disable grid", False, """Whether or not to disable use of the grid""")

temperature = Parameter("temperature", 25*celsius, """Simulation temperature""")
pressure = Parameter("pressure", 1*atm, """Simulation pressure. Ignored if a reflection sphere
                                           is used.""")

random_seed = Parameter("random seed", None, """Random number seed. Set this if you
                         want to have reproducible simulations.""")

lambda_values = Parameter("lambda values", [ 0.0, 0.142, 0.285, 0.429, 0.571, 0.714, 0.857, 1.0 ],
                          """The values of lambda to use in the RETI free energy simulation.""")

ligand_name = Parameter("ligand name", "LIG",
                        """The name of the ligand. This should be the name of one of the residues
                           in the ligand, so that this program can find the correct molecule.""")

reflection_radius = Parameter("reflection radius", None,
                              """The radius of the reflection sphere""")

topfile = Parameter("topfile", "system.top",
                    """Name of the topology file containing the system.""")

crdfile = Parameter("crdfile", "system.crd",
                    """Name of the coordinate file containing the coordinates of the 
                       system.""")

s3file = Parameter("s3file", "system.s3",
                   """Name to use for the intermediate s3 file that will contain the 
                      system after it has been loaded from the top/crd files.""")

outdir = Parameter("output directory", "output",
                   """Name of the directory in which to place all of the output files.""")

restart_file = Parameter("restart file", "quantomm_restart.s3",
                         """Name of the restart file to use to save progress during the simulation.""")

sysmoves_file = Parameter("sysmoves file", "quantomm_sysmoves.s3",
                          """Name of the file to save the initial QM/MM pre-simulation system.""")

nmoves = Parameter("nmoves", 200, """Number of RETI moves to perform during the simulation.""")

save_pdb = Parameter("save pdb", True,
                     """Whether or not to write a PDB of the system after each iteration.""")

save_all_pdbs = Parameter("save all pdbs", False,
                          """Whether or not to write all of the PDBs. If not, only PDBs at the two 
                             end points of the simulation will be written.""")

pdb_frequency = Parameter("pdb frequency", 10,
                          """The frequency (number of iterations between) saving PDBs""")

restart_frequency = Parameter("restart frequency", 1,
                              """The frequency (number of iterations between) saving the restart file for the simulation.""")

nslow = Parameter("nslow", 50,
                  """The number of 'slow' moves to perform per RETI iteration.""")

nfast = Parameter("nfast", 1000,
                  """The number of 'fast' moves to perform per slow move.""")

scale_charges = Parameter("scale charges", 1.0,
                          """The amount by which to scale MM charges in the QM/MM calculation.""")

intermolecular_only = Parameter("intermolecular only", False,
                                """Only calculate QM/MM intermolecular energies. This intramolecular energy of the QM atoms 
                                   is calculated using the MM forcefield.""")

amberhome = Parameter("amberhome", None,
                      """Use this to specify the installation directory of Amber / AmberTools, if you are using
                         the 'sqm' program. If this is not set, then you must have set this location using
                         the AMBERHOME environmental variable.""")

qm_program = Parameter("qm program", "sqm",
                       """The name of the program to use to calculate QM energies.""")

qm_executable = Parameter("qm executable", None,
                          """Exact path to the QM executable used to calculate the QM energies. If 
                             this is not set, then the QM executable will be searched from the path.""")

qm_method = Parameter("qm method", "AM1/d",
                      """The string passed to the QM program to specify the QM method to use to model the ligand.""")

basis_set = Parameter("basis set", "VDZ",
                      """The string passed to the QM program (molpro, as ab-initio only) to specify the basis set used to model the ligand.""")

qm_zero_energy = Parameter("qm zero energy", None,
                           """The value of 'zero' for the QM energy. This is used to shift the QM energy so that it
                              has the same comparable value as the MM energy. This is normally determined automatically,
                              but this option allows you to set the value manually.""")


def createQMMMMoves(system):
    # pull out all of the molecule groups for the mobile parts of the system
    print("Creating the Monte Carlo moves to sample the QM/MM system...")

    # create the global set of moves that will be applied to
    # the system
    moves = WeightedMoves()

    # create zmatrix moves to move the protein sidechains
    try:
        mobile_sidechains = system[MGName("mobile_sidechains")]

        if mobile_sidechains.nViews() > 0:
            sc_moves = ZMatMove(mobile_sidechains)
            moves.add( sc_moves, mobile_sidechains.nViews() )
    except:
        pass

    try:
        mobile_backbones = system[MGName("mobile_backbones")]

        if mobile_backbones.nViews() > 0:
            bb_moves = RigidBodyMC(mobile_backbones)
            bb_moves.setCenterOfRotation( GetCOGPoint( AtomName("CA", CaseInsensitive),
                                                       AtomName("N", CaseInsensitive) ) )

            bb_moves.setMaximumTranslation(0.030*angstrom)
            bb_moves.setMaximumRotation(1.0*degrees)
            moves.add( bb_moves, mobile_backbones.nViews() )
    except:
        pass

    use_reflection_sphere = False

    try:
        mobile_ligand = system[MGName("ligand")]

        if mobile_ligand.nViews() > 0:
            scale_moves = 10

            # get the amount to translate and rotate from the ligand's flexibility object
            flex = mobile_ligand.moleculeAt(0)[0].molecule().property("flexibility")

            # only move the solute if it is not the only molecule in the system
            if system.nMolecules() > 1 and (flex.translation().value() != 0 or flex.rotation().value() != 0):
                rb_moves = RigidBodyMC(mobile_ligand)
                rb_moves.setMaximumTranslation(flex.translation())
                rb_moves.setMaximumRotation(flex.rotation())

                if system.containsProperty("reflection sphere radius"):
                    reflection_radius = float(str(system.property("reflection sphere radius"))) * angstroms
                    reflection_center = system.property("reflection center").toVector()[0]
                    rb_moves.setReflectionSphere(reflection_center, reflection_radius)
                    use_reflection_sphere = True
                    print("Using the reflection sphere to constrain solvent moves...")

                scale_moves = scale_moves / 2
                moves.add( rb_moves, scale_moves * mobile_ligand.nViews() )

            intra_moves = InternalMove(mobile_ligand)
            moves.add( intra_moves, scale_moves * mobile_ligand.nViews() )
    except:
        pass

    try:
        mobile_solutes = system[MGName("mobile_solutes")]

        if mobile_solutes.nViews() > 0:
            rb_moves = RigidBodyMC(mobile_solutes)

            if system.containsProperty("average solute translation delta"):
                translation_delta = float(str(system.property("average solute translation delta")))
            else:
                translation_delta = 0

            if system.containsProperty("average solute rotation delta"):
                rotation_delta = float(str(system.property("average solute rotation delta")))
            else:
                rotation_delta = 0

            if translation_delta > 0 and rotation_delta > 0:
                rb_moves.setMaximumTranslation(translation_delta * angstroms)
                rb_moves.setMaximumRotation(rotation_delta * degrees)

                if system.containsProperty("reflection sphere radius"):
                    reflection_radius = float(str(system.property("reflection sphere radius"))) * angstroms
                    reflection_center = system.property("reflection center").toVector()[0]
                    rb_moves.setReflectionSphere(reflection_center, reflection_radius)

                moves.add(rb_moves, 4 * mobile_solutes.nViews())

            intra_moves = InternalMove(solute_group)
            moves.add(intra_moves, 4 * mobile_solutes.nViews())
    except:
        pass

    try:
        mobile_solvent = system[MGName("mobile_solvents")]

        max_water_translation = 0.15 * angstroms
        max_water_rotation = 15 * degrees

        if mobile_solvent.nViews() > 0:
            rb_moves = RigidBodyMC(mobile_solvent)
            rb_moves.setMaximumTranslation(max_water_translation)
            rb_moves.setMaximumRotation(max_water_rotation)

            if system.containsProperty("reflection sphere radius"):
                reflection_radius = float(str(system.property("reflection sphere radius"))) * angstroms
                reflection_center = system.property("reflection center").toVector()[0]
                rb_moves.setReflectionSphere(reflection_center, reflection_radius)

            moves.add(rb_moves, 4 * mobile_solvent.nViews())
    except:
        pass

    moves.setTemperature(temperature.val)

    try:
        if pressure.val and system.nMolecules() > 1 and system.property("space").isPeriodic():
            if not use_reflection_sphere:
                print("Running a constant pressure calculation")
                all = system[MGName("all")]
                volume_move = VolumeMove(all)
                volume_move.setTemperature(temperature.val)
                volume_move.setPressure(pressure.val)
                volume_move.setMaximumVolumeChange( 0.1 * all.nMolecules() * angstrom3 )
                moves.add(volume_move, 1)
    except:
        pass

    # Now create a multiple-timestep Monte Carlo move that
    # uses the above weighted moves for the fast energy
    mtsmc = MTSMC( moves, nfast.val )
    mtsmc.setFastEnergyComponent( Symbol("E_{fast}") )
    mtsmc.setSlowEnergyComponent( Symbol("E_{slow}") )

    seed = random_seed.val

    if seed is None:
        seed = RanGenerator().randInt(100000,1000000)
        print("Using generated random number seed %d" % seed)
    else:
        print("Using supplied random number seed %d" % seed)

    mtsmc.setGenerator( RanGenerator(seed) )

    return SameMoves(mtsmc)    


def printMoveInfo(moves):
    """Use this function to print out the specific details of the individual moves"""
    mtsmc = moves[0]
    fast_moves = mtsmc.fastMoves()

    for i in range(0, fast_moves.count()):
        move = fast_moves[i]
        print("%d : %s" % (i, move))


def getMinimumDistance(mol0, mol1):
    space = Cartesian()
    return space.minimumDistance(CoordGroup(mol0.molecule().property("coordinates").array()), \
                                 CoordGroup(mol1.molecule().property("coordinates").array()))


def setCLJProperties(forcefield, space):
    if cutoff_method.val.find("shift electrostatics") != -1:
        forcefield.setShiftElectrostatics(True)

    elif cutoff_method.val.find("reaction field") != -1:
        forcefield.setUseReactionField(True)
        forcefield.setReactionFieldDielectric(rf_dielectric.val)

    else:
        print("Cannot interpret the cutoff method from \"%s\"" % cutoff_method.val, file=sys.stderr)

    forcefield.setSpace(space)
    forcefield.setSwitchingFunction( HarmonicSwitchingFunction(coul_cutoff.val,coul_cutoff.val,
                                                               lj_cutoff.val,lj_cutoff.val) )

    return forcefield


def setFakeGridProperties(forcefield, space):
    forcefield.setSwitchingFunction( HarmonicSwitchingFunction(coul_cutoff.val,coul_cutoff.val,
                                                               lj_cutoff.val,lj_cutoff.val) )
    forcefield.setSpace(space)

    return forcefield


def setGridProperties(forcefield, extra_buffer=0*angstrom):
    forcefield.setGridSpacing(grid_spacing.val)
    forcefield.setBuffer(grid_buffer.val + extra_buffer)
    forcefield.setLJCutoff(lj_cutoff.val)
    forcefield.setCoulombCutoff(coul_cutoff.val)

    return forcefield


def setQMProperties(forcefield, space):
    forcefield.setSpace(space)
    forcefield.setSwitchingFunction( HarmonicSwitchingFunction(coul_cutoff.val,coul_cutoff.val,
                                                               lj_cutoff.val,lj_cutoff.val) )

    # calculate the total charge on the QM atoms
    total_charge = 0.0
    for molnum in forcefield.molecules().molNums():
        total_charge += forcefield[molnum].evaluate().charge().value()

    # round the charge to the nearest integer
    print("Charge on QM atoms is %s" % total_charge)
    
    total_charge = round(total_charge)
    print("Integer charge on QM atoms is %s" % total_charge)

    # Define the QM program used to calculate the QM part of the QM/MM energy
    if qm_program.val == "sqm":
        qm_prog = SQM()
        qm_prog.setMethod(qm_method.val)
        qm_prog.setTotalCharge(total_charge)

        if amberhome.val is None:
            ahome = os.getenv("AMBERHOME")

            if ahome is None:
                print("ERROR: YOU MUST SET THE AMBERHOME ENVIRONMENTAL VARIABLE!!!")
                print("This must point to the installation directory of Amber / AmberTools")
                sys.exit(0)
        else:
            ahome = amberhome.val

        if qm_executable.val is None:
            qm_prog.setExecutable( findExe("%s/bin/sqm" % ahome).absoluteFilePath() )
        else:
            qm_prog.setExecutable( findExe(qm_executable.val).absoluteFilePath() )

        qm_prog.setEnvironment("AMBERHOME", ahome)

        forcefield.setQuantumProgram(qm_prog)

    elif qm_program.val == "molpro":
        qm_prog = Molpro()
        qm_prog.setBasisSet(basis_set.val)
        qm_prog.setMethod(qm_method.val)
        qm_prog.setTotalCharge(total_charge)

        if qm_executable.val is None:
            qm_prog.setExecutable( findExe("molpro").absoluteFilePath() )
        else:
            qm_prog.setExecutable( findExe(qm_executable.val).absoluteFilePath() )

        forcefield.setQuantumProgram(qm_prog)

    elif qm_program.val is None:
        print("WARNING: Using a null (non-existant) QM program!")
        qm_prog = NullQM()
        forcefield.setQuantumProgram(qm_prog)

    else:
        print("WARNING: Could not recognise the QM program from %s. Will use the NULL program." \
                         % qm_program.val)

        forcefield.setQuantumProgram(NullQM())

    forcefield.setIntermolecularOnly( intermolecular_only.val )

    forcefield.setChargeScalingFactor( scale_charges.val )

    print("Using QM program %s, MM charges scaled by %s" % (qm_prog, scale_charges.val))

    return forcefield


def printEnergies(nrgs):
    keys = nrgs.keys()
    keys.sort()

    for key in keys:
        print("%s : %s kcal mol-1" % (key, nrgs[key]))


def loadQMMMSystem():
    """This function is called to set up the system. It sets everything
       up, then returns a System object that holds the configured system"""

    print("Loading the system...")

    t = QTime()

    if os.path.exists(s3file.val):
        print("Loading existing s3 file %s..." % s3file.val)
        loadsys = Sire.Stream.load(s3file.val)

    else:
        print("Loading from Amber files %s / %s..." % (topfile.val, crdfile.val))
        # Add the name of the ligand to the list of solute molecules
        sys_scheme = NamingScheme()
        sys_scheme.addSoluteResidueName(ligand_name.val)

        # Load up the system. This will automatically find the protein, solute, water, solvent
        # and ion molecules and assign them to different groups
        loadsys = createSystem(topfile.val, crdfile.val, sys_scheme)
        ligand_mol = findMolecule(loadsys, ligand_name.val)

        if ligand_mol is None:
            print("Cannot find the ligand (%s) in the set of loaded molecules!" % ligand_name.val)
            sys.exit(-1)

        # Center the system with the ligand at (0,0,0)
        loadsys = centerSystem(loadsys, ligand_mol)
        ligand_mol = loadsys[ligand_mol.number()][0].molecule()

        if reflection_radius.val is None:
            loadsys = addFlexibility(loadsys, naming_scheme=sys_scheme )
        else:
            loadsys = addFlexibility(loadsys, Vector(0), reflection_radius.val, naming_scheme=sys_scheme)

        Sire.Stream.save(loadsys, s3file.val)

    ligand_mol = findMolecule(loadsys, ligand_name.val)

    if ligand_mol is None:
        print("Cannot find the ligand (%s) in the set of loaded molecules!" % ligand_name.val)
        sys.exit(-1)

    # Now build the QM/MM system
    system = System("QMMM system")

    if loadsys.containsProperty("reflection center"):
        reflect_center = loadsys.property("reflection center").toVector()[0]
        reflect_radius = float(str(loadsys.property("reflection sphere radius")))

        system.setProperty("reflection center", AtomCoords(CoordGroup(1,reflect_center)))
        system.setProperty("reflection sphere radius", VariantProperty(reflect_radius))
        space = Cartesian()
    else:
        space = loadsys.property("space")

    if loadsys.containsProperty("average solute translation delta"):
        system.setProperty("average solute translation delta", \
                           loadsys.property("average solute translation delta"))

    if loadsys.containsProperty("average solute rotation delta"):
        system.setProperty("average solute rotation delta", \
                           loadsys.property("average solute rotation delta"))

    # create a molecule group to hold all molecules
    all_group = MoleculeGroup("all")

    # create a molecule group for the ligand
    ligand_group = MoleculeGroup("ligand")
    ligand_group.add(ligand_mol)
    all_group.add(ligand_mol)

    groups = []
    groups.append(ligand_group)

    # pull out the groups that we want from the two systems

    # create a group to hold all of the fixed molecules in the bound leg
    fixed_group = MoleculeGroup("fixed_molecules")
    if MGName("fixed_molecules") in loadsys.mgNames():
        fixed_group.add( loadsys[ MGName("fixed_molecules") ] )

    if save_pdb.val:
        # write a PDB of the fixed atoms in the bound and free legs
        if not os.path.exists(outdir.val):
            os.makedirs(outdir.val)

        PDB().write(fixed_group, "%s/fixed.pdb" % outdir.val)

    # create a group to hold all of the mobile solute molecules
    mobile_solutes_group = MoleculeGroup("mobile_solutes")
    if MGName("mobile_solutes") in loadsys.mgNames():
        mobile_solutes_group.add( loadsys[MGName("mobile_solutes")] )
        mobile_solutes_group.remove(ligand_mol)
        if mobile_solutes_group.nMolecules() > 0:
            all_group.add(mobile_solutes_group)
    
    groups.append(mobile_solutes_group)

    # create a group to hold all of the mobile solvent molecules
    mobile_solvents_group = MoleculeGroup("mobile_solvents")
    if MGName("mobile_solvents") in loadsys.mgNames():
        mols = loadsys[MGName("mobile_solvents")]
        for molnum in mols.molNums():
            solvent_mol = mols[molnum][0].molecule()
            mobile_solvents_group.add(solvent_mol)

        all_group.add(mobile_solvents_group)

        print("The number of mobile solvent molecules is %d." % mobile_solvents_group.nMolecules())

    groups.append(mobile_solvents_group)

    # create the groups to hold all of the protein molecules. We will use "extract" to 
    # pull out only those protein atoms that are in the mobile region
    protein_intra_group = MoleculeGroup("protein_intra_group")
    mobile_proteins_group = MoleculeGroup("proteins")
    mobile_protein_sidechains_group = MoleculeGroup("mobile_sidechains")
    mobile_protein_backbones_group = MoleculeGroup("mobile_backbones")

    if MGName("protein_sidechains") in loadsys.mgNames() or \
       MGName("protein_backbones") in loadsys.mgNames():

        all_proteins = Molecules()

        try:
            protein_sidechains = loadsys[MGName("protein_sidechains")]
            all_proteins.add(protein_sidechains.molecules())
        except:
            protein_sidechains = MoleculeGroup()

        try:
            protein_backbones = loadsys[MGName("protein_backbones")]
            all_proteins.add(protein_backbones.molecules())
        except:
            protein_backbones = MoleculeGroup()

        try:
            boundary_molecules = loadsys[MGName("boundary_molecules")]
            all_proteins.add(boundary_molecules.molecules())
        except:
            boundary_molecules = MoleculeGroup()

        for molnum in all_proteins.molNums():
            protein_mol = Molecule.join(all_proteins[molnum])
            
            if protein_mol.selectedAll():
                protein_intra_group.add(protein_mol)
                all_group.add(protein_mol)

                mobile_protein = []

                if protein_sidechains.contains(molnum):
                    sidechains = protein_sidechains[molnum]
                    for sidechain in sidechains:
                        mobile_protein_sidechains_group.add( sidechain )

                    mobile_protein += sidechains

                if protein_backbones.contains(molnum):
                    backbones = protein_backbones[molnum]
                    for backbone in backbones:
                        mobile_protein_backbones_group.add( backbone )

                    mobile_protein += backbones

                if len(mobile_protein) > 0:
                    mobile_proteins_group.add( Molecule.join(mobile_protein) )

            else:
                # only some of the atoms have been selected. We will extract
                # the mobile atoms and will then update all of the other selections
                print("Extracting the mobile atoms of protein %s" % protein_mol.molecule())
                new_protein_mol = protein_mol.extract()
                print("Extracted %d mobile atoms from %d total atoms..." % \
                                        (new_protein_mol.nAtoms(), protein_mol.molecule().nAtoms()))

                protein_intra_group.add(new_protein_mol)
                all_group.add( new_protein_mol )

                mobile_protein_view = new_protein_mol.selection()
                mobile_protein_view = mobile_protein_view.selectNone()

                if protein_sidechains.contains(molnum):
                    sidechains = protein_sidechains[molnum]

                    for sidechain in sidechains:
                        view = new_protein_mol.selection()
                        view = view.selectNone()

                        for atomid in sidechain.selection().selectedAtoms():
                            atom = protein_mol.atom(atomid)
                            resatomid = ResAtomID( atom.residue().number(), atom.name() )
                            view = view.select( resatomid )
                            mobile_protein_view = mobile_protein_view.select( resatomid )

                        if view.nSelected() > 0:
                            mobile_protein_sidechains_group.add( PartialMolecule(new_protein_mol, view) )

                if protein_backbones.contains(molnum):
                    backbones = protein_backbones[molnum]

                    for backbone in backbones:
                        view = new_protein_mol.selection()
                        view = view.selectNone()

                        for atomid in backbone.selection().selectedAtoms():
                            atom = protein_mol.atom(atomid)
                            resatomid = ResAtomID( atom.residue().number(), atom.name() )
                            view = view.select( resatomid )
                            mobile_protein_view = mobile_protein_view.select( resatomid )

                        if view.nSelected() > 0:
                            mobile_protein_backbones_group.add( PartialMolecule(new_protein_mol, view) )

                print("Number of moved protein sidechain residues = %s" % mobile_protein_sidechains_group.nViews())
                print("Number of moved protein backbone residues = %s" % mobile_protein_backbones_group.nViews())

                if mobile_protein_view.nSelected() > 0:
                    mobile_proteins_group.add( PartialMolecule(new_protein_mol, mobile_protein_view) )

    groups.append(mobile_protein_backbones_group)
    groups.append(mobile_protein_sidechains_group)
    groups.append(all_group)

    # finished added in all of the proteins
    for group in groups:
        if group.nMolecules() > 0:
            print("Adding group %s" % group.name())
            system.add(group)

    # now add in the forcefields for the system...
    print("Creating the forcefields for the QM/MM system...")

    # first, group together the molecules grouped above into convenient
    # groups for the forcefields

    # group holding just the ligand
    ligand_mols = ligand_group.molecules()

    # group holding all of the mobile atoms
    mobile_mols = mobile_solvents_group.molecules()
    mobile_mols.add( mobile_solutes_group.molecules() )
    mobile_mols.add( protein_intra_group.molecules() )

    # group holding all of the mobile atoms in the bound leg, excluding the 
    # buffer atoms that are fixed, but bonded to mobile atoms
    mobile_buffered_mols = mobile_solvents_group.molecules()
    mobile_buffered_mols.add( mobile_solutes_group.molecules() )
    mobile_buffered_mols.add( mobile_proteins_group.molecules() )

    # group holding all of the protein molecules that need intramolecular terms calculated
    protein_intra_mols = protein_intra_group.molecules()

    # group holding all of the solute molecules that nede intramolecular terms calculated
    solute_intra_mols = mobile_solutes_group.molecules()

    forcefields = []

    ###
    ### INTRA-ENERGY OF THE LIGAND AND CLUSTER
    ###
    
    # intramolecular energy of the ligand
    ligand_intraclj = IntraCLJFF("ligand:intraclj")
    ligand_intraclj = setCLJProperties(ligand_intraclj, space)
    ligand_intraclj.add(ligand_mols)

    ligand_intraff = InternalFF("ligand:intra")
    ligand_intraff.add(ligand_mols)

    forcefields.append(ligand_intraclj)
    forcefields.append(ligand_intraff)

    ligand_mm_nrg = ligand_intraclj.components().total() + ligand_intraff.components().total()

    ###
    ### FORCEFIELDS INVOLVING THE LIGAND/CLUSTER AND OTHER ATOMS
    ###

    # forcefield holding the energy between the ligand and the mobile atoms in the
    # bound leg
    ligand_mobile = InterGroupCLJFF("system:ligand-mobile")
    ligand_mobile = setCLJProperties(ligand_mobile, space)

    ligand_mobile.add(ligand_mols, MGIdx(0))
    ligand_mobile.add(mobile_mols, MGIdx(1))

    qm_ligand = QMMMFF("system:ligand-QM")    
    qm_ligand.add(ligand_mols, MGIdx(0))
    qm_ligand = setQMProperties(qm_ligand, space)

    zero_energy = 0

    if not intermolecular_only.val:
        if qm_zero_energy.val is None:
            # calculate the delta value for the system - this is the difference between
            # the MM and QM intramolecular energy of the ligand
            t.start()
            print("\nComparing the MM and QM energies of the ligand...")
            mm_intra = ligand_intraclj.energy().value() + ligand_intraff.energy().value()
            print("MM energy = %s kcal mol-1 (took %s ms)" % (mm_intra, t.elapsed()))

            t.start()
            zero_sys = System()
            zero_sys.add(qm_ligand)
            qm_intra = zero_sys.energy().value()
            print("QM energy = %s kcal mol-1 (took %s ms)" % (qm_intra, t.elapsed()))

            print("\nSetting the QM zero energy to %s kcal mol-1" % (qm_intra - mm_intra))
            qm_ligand.setZeroEnergy( (qm_intra-mm_intra) * kcal_per_mol )
            zero_energy = qm_intra - mm_intra
        else:
            print("\nManually setting the QM zero energy to %s" % qm_zero_energy.val)
            qm_ligand.setZeroEnergy( qm_zero_energy.val )
            zero_energy = qm_zero_energy.val

    qm_ligand.add(mobile_mols, MGIdx(1))

    ligand_mm_nrg += ligand_mobile.components().total()
    ligand_qm_nrg = qm_ligand.components().total() + ligand_mobile.components().lj()

    if intermolecular_only.val:
        # the QM model still uses the MM intramolecular energy of the ligand
        ligand_qm_nrg += ligand_intraclj.components().total() + ligand_intraff.components().total()

    forcefields.append(ligand_mobile)
    forcefields.append(qm_ligand)

    if fixed_group.nMolecules() > 0:
        # there are fixed molecules

        # Whether or not to disable the grid and calculate all energies atomisticly
        if disable_grid:
            # we need to renumber all of the fixed molecules so that they don't clash
            # with the mobile molecules
            print("Renumbering fixed molecules...")
            fixed_group = renumberMolecules(fixed_group)

        # forcefield holding the energy between the ligand and the fixed atoms in the bound leg
        if disable_grid:
            ligand_fixed = InterGroupCLJFF("system:ligand-fixed")
            ligand_fixed = setCLJProperties(ligand_fixed, space)
            ligand_fixed = setFakeGridProperties(ligand_fixed, space)

            ligand_fixed.add(ligand_mols, MGIdx(0))
            ligand_fixed.add(fixed_group, MGIdx(1))

            qm_ligand.add(fixed_group, MGIdx(1))

            ligand_mm_nrg += ligand_fixed.components().total()
            ligand_qm_nrg += ligand_fixed.components().lj()

            forcefields.append(ligand_fixed)

        else:
            ligand_fixed = GridFF2("system:ligand-fixed")
            ligand_fixed = setCLJProperties(ligand_fixed, space)
            ligand_fixed = setGridProperties(ligand_fixed)

            ligand_fixed.add(ligand_mols, MGIdx(0))
            ligand_fixed.addFixedAtoms( fixed_group )

            qm_ligand.addFixedAtoms( fixed_group )

            ligand_mm_nrg += ligand_fixed.components().total()
            ligand_qm_nrg += ligand_fixed.components().lj()

            forcefields.append(ligand_fixed)

    ###
    ### FORCEFIELDS NOT INVOLVING THE LIGAND
    ###

    # forcefield holding the intermolecular energy between all molecules
    mobile_mobile = InterCLJFF("mobile-mobile")
    mobile_mobile = setCLJProperties(mobile_mobile, space)

    mobile_mobile.add(mobile_mols)

    other_nrg = mobile_mobile.components().total()
    forcefields.append(mobile_mobile)

    # forcefield holding the energy between the mobile atoms and  
    # the fixed atoms
    if disable_grid.val:
        mobile_fixed = InterGroupCLJFF("mobile-fixed")
        mobile_fixed = setCLJProperties(mobile_fixed)
        mobile_fixed = setFakeGridProperties(mobile_fixed, space)
        mobile_fixed.add(mobile_buffered_mols, MGIdx(0))
        mobile_fixed.add(fixed_group, MGIdx(1))
        other_nrg += mobile_fixed.components().total()
        forcefields.append(mobile_fixed)
    else:
        mobile_fixed = GridFF2("mobile-fixed")
        mobile_fixed = setCLJProperties(mobile_fixed, space)
        mobile_fixed = setGridProperties(mobile_fixed)

        # we use mobile_buffered_group as this group misses out atoms that are bonded
        # to fixed atoms (thus preventing large energies caused by incorrect non-bonded calculations)
        mobile_fixed.add(mobile_buffered_mols, MGIdx(0))
        mobile_fixed.addFixedAtoms(fixed_group)
        other_nrg += mobile_fixed.components().total()
        forcefields.append(mobile_fixed)

    # intramolecular energy of the protein
    if protein_intra_mols.nMolecules() > 0:
        protein_intraclj = IntraCLJFF("protein_intraclj")
        protein_intraclj = setCLJProperties(protein_intraclj, space)

        protein_intraff = InternalFF("protein_intra")

        for molnum in protein_intra_mols.molNums():
            protein_mol = Molecule.join(protein_intra_mols[molnum])
            protein_intraclj.add(protein_mol)
            protein_intraff.add(protein_mol)

        other_nrg += protein_intraclj.components().total()
        other_nrg += protein_intraff.components().total()
        forcefields.append(protein_intraclj)
        forcefields.append(protein_intraff)

    # intramolecular energy of any other solutes
    if solute_intra_mols.nMolecules() > 0:
        solute_intraclj = IntraCLJFF("solute_intraclj")
        solute_intraclj = setCLJProperties(solute_intraclj, space)

        solute_intraff = InternalFF("solute_intra")

        for molnum in solute_intra_mols.molNums():
            solute_mol = Molecule.join(solute_intra_mols[molnum])
            solute_intraclj.add(solute_mol)
            solute_intraff.add(solute_mol)

        other_nrg += solute_intraclj.components().total()
        other_nrg += solute_intraff.components().total()
        forcefields.append(solute_intraclj)
        forcefields.append(solute_intraff)

    ###
    ### NOW ADD THE FORCEFIELDS TO THE SYSTEM
    ###
    ###
    ### SETTING THE FORCEFIELD EXPRESSIONS
    ###

    lam = Symbol("lambda")

    e_slow = ((1-lam) * ligand_qm_nrg) + (lam * ligand_mm_nrg) + other_nrg
    e_fast = ligand_mm_nrg + other_nrg

    de_by_dlam = ligand_mm_nrg - ligand_qm_nrg

    for forcefield in forcefields:
        system.add(forcefield)

    system.setConstant(lam, 0.0)

    system.setComponent(Symbol("E_{fast}"), e_fast)
    system.setComponent(Symbol("E_{slow}"), e_slow)
    system.setComponent(Symbol("dE/dlam"), de_by_dlam)
    system.setComponent( system.totalComponent(), e_slow )
 
    system.setProperty("space", space)
    
    if space.isPeriodic():
        # ensure that all molecules are wrapped into the space with the ligand at the center
        print("Adding in a space wrapper constraint %s, %s" % (space, ligand_mol.evaluate().center()))
        system.add( SpaceWrapper( ligand_mol.evaluate().center(), all_group ) )
        system.applyConstraints()

    print("\nHere are the values of all of the initial energy components...")
    t.start()
    printEnergies(system.energies())
    print("(these took %d ms to evaluate)\n" % t.elapsed())

    # Create a monitor to monitor the free energy average
    system.add( "dG/dlam", MonitorComponent(Symbol("dE/dlam"), AverageAndStddev()) )

    if intermolecular_only.val:
        print("\n\n## This simulation uses QM to model *only* the intermolecular energy between")
        print("## the QM and MM atoms. The intramolecular energy of the QM atoms is still")
        print("## modelled using MM.\n")
    else:
        print("\n\n## This simulation uses QM to model both the intermolecular and intramolecular")
        print("## energies of the QM atoms. Because the this, we have to adjust the 'zero' point")
        print("## of the QM potential. You need to add the value %s kcal mol-1 back onto the" % zero_energy)
        print("## QM->MM free energy calculated using this program.\n")

    return system


def makeRETI(system, moves):
    """This function replicates 'system' over each of the supplied lambda values
       and uses 'moves' to sample each of the replicated systems. This uses RETI
       to perform replica exchange moves across lambda"""

    lam = Symbol("lambda")

    replicas = Replicas( len(lambda_values.val) )

    replicas.setSubSystem(system)
    replicas.setSubMoves(moves)
    replicas.setNSubMoves(nslow.val)
    replicas.setLambdaComponent(lam)
    replicas.setRecordAllStatistics(True)        

    seed = random_seed.val
    
    if seed is None:
        seed = RanGenerator().randInt(100000,1000000)
        print("RETI system using generated random number seed %d" % seed)
    else:
        print("RETI system using supplied random number seed %d" % seed)
    
    replicas.setGenerator( RanGenerator(seed+5) )

    for i in range(0, len(lambda_values.val)):
        # set the initial lambda value for this replica
        replicas.setLambdaValue(i, lambda_values.val[i])

    for i in range(0, len(lambda_values.val)):
        print(lambda_values.val[i])
        print(replicas[i].subSystem().constants())

    # Now add monitors for each replica that will copy back
    nrgmons = [ "dG/dlam" ]

    for nrgmon in nrgmons:
        replicas.add( nrgmon, MonitorMonitor(MonitorName(nrgmon), True) )

    # now create the replica exchange moves for the replicas
    replica_moves = RepExMove()
    replica_moves.setDisableSwaps(False)
    replica_moves.setGenerator( RanGenerator(seed+7) )

    print("\nReturning the QM/MM RETI replicas and moves...")
    return (replicas, replica_moves)


def loadQMMM():
    """This loads the QM/MM system"""
    system = loadQMMMSystem()
    moves = createQMMMMoves(system)

    print("\nUsing moves:")
    printMoveInfo(moves)

    return (system, moves)


def analyseQMMM(replicas, iteration, ti_freenrgs):
    """This function is used to perform all analysis of iteration 'it' of the passed QMMM system"""

    print("Analysing iteration %d..." % iteration)

    if not os.path.exists(outdir.val):
        os.makedirs(outdir.val)

    # read the value of delta_lambda from the first system
    system = replicas[0].subSystem()

    logfile = "%s/results_%0004d.log" % (outdir.val, iteration)

    FILE = open(logfile, "w")

    print("===========================", file=FILE)
    print(" Results for iteration %d" % iteration, file=FILE)
    print("===========================", file=FILE)

    print("temperature == %f K\n" % replicas[0].subMoves().temperature().to(kelvin), file=FILE) 

    nreplicas = replicas.nReplicas()

    # extract all of the monitors from the replicas
    lambda_values = []

    dg = {}

    write_pdbs = (save_pdb.val) and (iteration % pdb_frequency.val == 0)

    if write_pdbs:
        print("Saving PDBs of the system at iteration %d..." % iteration)

    for i in range(0, nreplicas):
        replica = replicas[i]
        monitors = replica.monitors()
        lamval = replica.lambdaValue()
        lambda_values.append(lamval)

        if write_pdbs:
            if save_all_pdbs.val or (i == 0) or (i == nreplicas-1):
                # Save a PDB of the final configuration for the bound and free legs for each lambda value
                system = replica.subSystem()
                all = system[MGName("all")]

                PDB().write(all, "%s/coords_%000006d_%.5f.pdb" % (outdir.val, iteration, lamval))

        dg[lamval] = monitors[MonitorName("dG/dlam")][-1].accumulator()

    ti_freenrgs.set( iteration, dg )

    print("\n=============", file=FILE)
    print("Correction free energy for iteration %d equals %s" % (iteration, \
                        ti_freenrgs[iteration].integrate().values()[-1].y()), file=FILE)
    print("==============", file=FILE)


@resolveParameters
def run():
    """This function is used to actually run the simulation according to the 
       parameters passed in in 'params'"""

    t = QTime()
    total_t = QTime()
    total_t.start()


    if os.path.exists(restart_file.val):
        t.start()
        (qm_system,qm_moves) = Sire.Stream.load(restart_file.val)
        print("Loading the restart file took %d ms" % t.elapsed())
    else:
        # Load the system from the files
        t.start()
        if os.path.exists(sysmoves_file.val):
            (qm_system,qm_moves) = Sire.Stream.load(sysmoves_file.val)
        else:
            (qm_system,qm_moves) = loadQMMM()            
            Sire.Stream.save( (qm_system,qm_moves), sysmoves_file.val)

        # now replicate the system across all lambda values so that we can
        # run a simulation
        (qm_system,qm_moves) = makeRETI(qm_system,qm_moves)

        Sire.Stream.save( (qm_system,qm_moves), restart_file.val )

        print("Initialising the system took %d ms" % t.elapsed())

    # see how many blocks of moves we still need to perform...
    nattempted = qm_moves.nMoves()

    print("Number of iterations to perform: %d. Number of iterations completed: %d." % (nmoves.val, nattempted))

    # See if we have any existing free energy statistics files...
    t.start()
    freenrgs_file = "%s/freenrgs.s3" % outdir.val

    if not os.path.exists(freenrgs_file):
        ti_freenrgs = TI()
    else:
        ti_freenrgs = Sire.Stream.load(freenrgs_file)

    print("Initialising / loading the free energy files took %d ms" % t.elapsed())

    for i in range(nattempted+1, nmoves.val+1):
        t.start()
        print("Performing iteration %d..." % i)
        sim = SupraSim.run( qm_system, qm_moves, 1, True )
        sim.wait()

        qm_system = sim.system()
        qm_moves = sim.moves()

        print("...iteration complete (took %d ms)" % t.elapsed())

        # write a restart file every N moves in case of crash or run out of time
        if i % restart_frequency.val == 0 or i == nmoves.val:
            t.start()
            print("Saving the restart file from iteration %d..." % i)
            # save the old file to a backup
            try:
                shutil.copy(restart_file.val, "%s.bak" % restart_file.val)
            except:
                pass

            Sire.Stream.save( (qm_system, qm_moves), restart_file.val )
            print("...save complete (took %d ms)" % t.elapsed())

        t.start()
        print("Analysing iteration %d..." % i)
        analyseQMMM(qm_system, i, ti_freenrgs)
        qm_system.clearAllStatistics()
        print("...analysis complete (took %d ms)" % t.elapsed())

        if i % restart_frequency.val == 0 or i == nmoves.val:
            t.start()
            print("Saving the free energy analysis files from iteration %d..." % i)
            # save the old file to a backup
            try:
                shutil.copy(freenrgs_file, "%s.bak" % freenrgs_file)
            except:
                pass

            Sire.Stream.save( ti_freenrgs, freenrgs_file )

            print("...save complete (took %d ms)" % t.elapsed())

    print("All iterations complete. Total runtime was %d ms" % total_t.elapsed())
