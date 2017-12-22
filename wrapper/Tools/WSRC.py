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

import Sire.Config
import Sire.Stream

from Sire.Tools.AmberLoader import *
from Sire.Tools import Parameter, resolveParameters

import os
import shutil
import copy

wsrc_tools_dir = "%s/Tools/WSRC" % Sire.Config.share_directory

####################################################
# ALL OF THE GLOBAL USER-AVAILABLE WSRC PARAMETERS #
####################################################

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

fast_sim = Parameter("fast simulation", False,
                     """Switch on options that simplify the calculation so that it is significantly
                        sped up. The calculation has much lower accuracy, but could be useful
                        as a way to quickly get residue values, or to equilibrate a longer calculation""")

use_oldff = Parameter("use old forcefields", False, """For debugging, use the old forcefields rather than the 
                                                       new forcefields""")

temperature = Parameter("temperature", 25*celsius, """Simulation temperature""")
random_seed = Parameter("random seed", None, """Random number seed. Set this if you
                         want to have reproducible simulations.""")

identity_atoms = Parameter("identity atoms", None,
                           """The list of atom names in the ligand on which to place
                              identity points. If this is not set, then the identity atoms 
                              will be generated automatically.""")

use_fixed_points = Parameter("fixed points", False,
                             """Whether or not to use fixed identity points based on looking at
                                the overlap with the atoms""")

use_water_points = Parameter("water points", False,
                             """Whether or not to move the identity points to the oxygens of
                                the swap water molecules, and so keep them fixed in space during
                                the simulation""")

use_fixed_ligand = Parameter("fixed ligand", False,
                             """Whether or not to completely fix the ligand during the simulation.""")

use_rot_trans_ligand = Parameter("ligand rot-trans", True,
                                 """Whether or not the ligand is free to translate and rotate.""")

use_reflect_volume = Parameter("reflect volume", False,
                               """Use the reflection volume instead of the identity constraint to hold
                                  the swap water cluster in place.""")

reflect_volume_radius = Parameter("reflect volume radius", 1.75*angstrom,
                                  """The radius of the reflection volume, used only if the reflection volume
                                     is used to hold the swap water cluster in place.""")

reflect_volume_buffer = Parameter("reflect volume buffer", 0*angstrom,
                                  """The buffer beyond the reflection volume radius that is used when selecting
                                     the water molecules that will be swapped. Swapped water molecules will be those
                                     that are within 'reflect volume radius + reflect volume buffer' of any of 
                                     the heavy atoms of the swapped ligand.""")

n_equil_swap = Parameter("swap water nequilmoves", 5000,
                         """The number of moves to equilibrate the swap water cluster before applying
                            the identity or reflection volume constraint.""")

alpha_scale = Parameter("alpha_scale", 1.0,
                        """Amount by which to scale the alpha parameter. The lower the value,
                           the less softening with lambda, while the higher the value, the
                           more softening""")

delta_lambda = Parameter("delta_lambda", 0.001,
                         """Value of delta lambda used in the finite difference thermodynamic
                            integration algorithm used to calculate the free energy""")

water_monitor_distance = Parameter("water monitor distance", 5.0*angstrom,
                                   """The distance up to which the free energy of water molecules
                                      interacting with the ligand should be recorded.""")

waterbox_only = Parameter("waterbox only", False,
                          """Whether or not to select water molecules only from the water box.""")

nrgmon_frequency = Parameter("energy monitor frequency", 1000, 
                             """The number of steps between each evaluation of the energy monitors.""")

save_all_nrgmons = Parameter("save energy monitors", False,
                             """When debugging, you may want to switch on the saving of energy
                                monitors. Normally you shouldn't need to save these.""")

lambda_values = Parameter("lambda values", [ 0.005, 0.071, 0.137, 0.203, 0.269, 0.335, 0.401, 0.467, 0.533, 0.599, 0.665, 0.731, 0.797, 0.863, 0.929, 0.995 ],
                          """The values of lambda to use in the RETI free energy simulation. Note that it is not a good idea
                             to use lambda values 0 or 1 as this will introduce discontinuities into the PMF.""")

nsubmoves = Parameter("nsubmoves", 50000,
                      """The number of moves to perform between each RETI move.""")

ligand_name = Parameter("ligand name", None,
                        """The name of the ligand. This should be the name of one of the residues
                           in the ligand, so that this program can find the correct molecule. If it is not set, then
                           the first non-protein, non-solvent molecule is used.""")

reflection_radius = Parameter("reflection radius", 15*angstrom,
                              """The radius of the reflection sphere""")

ligand_reflection_radius = Parameter("ligand reflection radius", 1*angstrom,
                                     """The reflection radius of the ligand. This is used to constrain the ligand
                                        to remain in the active site. This is needed to define the accessible volume
                                        of the bound state.""")

protein_topfile = Parameter("protein topfile", "proteinbox.top",
                            """Name of the topology file containing the solvated protein-ligand complex.""")

protein_crdfile = Parameter("protein crdfile", "proteinbox.crd",
                            """Name of the coordinate file containing the coordinates of the 
                               solvated protein-ligand complex.""")

protein_s3file = Parameter("protein s3file", "proteinbox.s3",
                           """Name to use for the intermediate s3 file that will contain the 
                              solvated protein-ligand complex after it has been loaded from the top/crd files.""")

water_topfile = Parameter("water topfile", "%s/waterbox.top" % wsrc_tools_dir,
                          """Name of the topology file containing the water box.""")

water_crdfile = Parameter("water crdfile", "%s/waterbox.crd" % wsrc_tools_dir,
                          """Name of the coordinate file containing the coordinates of the water box.""")

water_s3file = Parameter("water s3file", "waterbox.s3",
                         """Name to use for the intermediate s3 file that will contain the 
                            water box after it has been loaded from the top/crd files.""")

outdir = Parameter("output directory", "output",
                   """Name of the directory in which to place all of the output files.""")

restart_file = Parameter("restart file", "wsrc_restart.s3",
                         """Name of the restart file to use to save progress during the simulation.""")

sysmoves_file = Parameter("sysmoves file", "wsrc_sysmoves.s3",
                          """Name of the file to save the initial WSRC pre-simulation system.""")

nequilmoves = Parameter("nequilmoves", 50000,
                        """Number of equilibration moves to perform before setting up the free energy simulation.""")

nmoves = Parameter("nmoves", 1000, """Number of RETI moves to perform during the simulation.""")

move_backbone = Parameter("move backbone", True, 
                          """Whether or not to move the protein backbone.""")

coul_power = Parameter("coulomb power", 0,
                       """The soft-core coulomb power parameter""")

shift_delta = Parameter("shift delta", 1.2,
                        """The soft-core shift delta parameter""")

soften_water = Parameter("soften water", 1.1, 
                         """The amount by which to scale the water-water electrostatic interactions in 
                            the swap-water cluster between lambda=0 and lambda=1. This helps keep the cluster
                            together as it is swapped between the two boxes.""")

lj_buffer = Parameter("LJ buffer", 0.005,
                      """To prevent end-point effects, the scale factor for the LJ interactions cannot fully
                         move between 0 and 1, as configurations sampled at 0 or 1 will be invalid for other states.
                         To overcome this problem, the LJ lambda scale is moved from "LJ buffer" to "1 - LJ buffer",
                         e.g. for the default buffer of 0.001, the LJ lambda scale is from 0.001 to 0.999.""")

uncharge_ligand = Parameter("uncharge ligand", False,
                            """Whether or not to uncharge the ligand (and swap water cluster) before 
                               swapping them. They are then recharged at the end of the swap.""")

uncharge_ligand_max = Parameter("uncharge ligand max", 0.5,
                                """The maximum amount to uncharge the ligand (and swap water cluster) before
                                   swapping them. A value of 1.0 will not uncharge them at all, while a value of
                                   0.0 will uncharge them completely. The default value is 0.5, which will 50%
                                   uncharge the ligand and swap water cluster before swapping.""")

uncharge_lambda_values = Parameter("uncharge lambda values", [0.0, 0.1, 0.25, 0.45, 0.55, 0.75, 0.9, 1.0],
                                   """Lambda values to use when uncharging (and then recharging) the ligand. These will be 
                                      added onto the swapped lambda values to give a new range squashed between 0 and 1.""")

save_pdb = Parameter("save pdb", True,
                     """Whether or not to write a PDB of the system after each iteration.""")

save_all_pdbs = Parameter("save all pdbs", False,
                          """Whether or not to write all of the PDBs. If not, only PDBs at the two 
                             end points of the simulation will be written.""")

pdb_frequency = Parameter("pdb frequency", 100,
                          """The frequency (number of iterations between) saving PDBs""")

binwidth = Parameter("free energy bin width", 1 * kcal_per_mol,
                     """The size of the bins used in the histogram of energies collected
                        as part of creating the free energy average, in multiples of delta lambda""")

restart_frequency = Parameter("restart frequency", 10,
                              """The frequency (number of iterations between) saving the restart file for the simulation.""")

####################################################

def getOverlapWaters(ligand, waters, radius=2*angstrom):

    overlaps = []

    space = Cartesian()

    # get the coordinates of all heavy atoms
    coords = []

    for atom in ligand.atoms():
        try:
            if atom.property("element").nProtons() >= 6:
                coords.append( atom.property("coordinates") )
        except:
            if atom.property("mass").value() >= 12:
                coords.append( atom.property("coordinates") )

    coords = CoordGroup(coords)

    for molnum in waters.molNums():
        water = waters[molnum][0].molecule()

        oxygen = None

        for atom in water.atoms():
            if atom.property("element").nProtons() == 8:
                oxygen = atom
                break

        if oxygen is None:
            oxygen = water.atoms()[0]

        mindist = space.minimumDistance( CoordGroup([oxygen.property("coordinates")]), coords)

        if mindist < radius.value():
            overlaps.append(oxygen.molecule())

    return overlaps


def getIdentityPoints(ligand):
    
    atoms = ligand.atoms()

    have_point = {}

    for atom in atoms:
        # skip small atoms
        try:
            if atom.property("element").nProtons() >= 6:
                have_point[str(atom.name().value())] = True

            else:
                have_point[str(atom.name().value())] = False
        except:
            try:
                if atom.property("mass").value() >= 12:
                    have_point[str(atom.name().value())] = True

                else:
                    have_point[str(atom.name().value())] = False
            except:
                print("Atom %s has neither a mass or element. Cannot add an identity point." % str(atom))
                have_point[str(atom.name().value())] = False

    if ligand.hasProperty("connectivity"):
        connectivity = ligand.property("connectivity")
    else:
        connectivity = Connectivity(ligand)

    have_point_keys = list(have_point.keys())
    have_point_keys.sort()

    for key in list(have_point_keys):
        if have_point[key]:
            # if this is bonded to 3+ atoms that also have 
            # identity points, then get rid of this point
            atom = ligand.atom( AtomName(key) )

            bonded = connectivity.connectionsTo(atom.name())

            if len(bonded) >=3:
                n = 0

                for b in bonded:
                    if have_point[ str(ligand.atom(b).name().value()) ]:
                        n += 1

                if n >= 3:
                    print("Skipping %s as it is bonded to 3 or more points..." % atom.name())
                    have_point[key] = False

    identity_points = []
    k2 = []

    # skip every 8th point
    iskip = 0

    have_point_keys = list(have_point.keys())
    have_point_keys.sort()

    for key in list(have_point_keys):
        if have_point[key]:
            iskip += 1

            if iskip == 8:
                iskip = 0
            else:
                k2.append(key)
                identity_points.append( ligand.atom( AtomName(key) ) )

    k2.sort()

    print("Using %d identity points: %s" % (len(k2), str(k2)))

    return identity_points


def getMinimumDistance(mol0, mol1):
    space = Cartesian()
    return space.minimumDistance(CoordGroup(mol0.molecule().property("coordinates").array()), \
                                 CoordGroup(mol1.molecule().property("coordinates").array()))


def setCLJProperties(forcefield):
    if cutoff_method.val.find("shift electrostatics") != -1:
        forcefield.setShiftElectrostatics(True)

    elif cutoff_method.val.find("reaction field") != -1:
        forcefield.setUseReactionField(True)
        forcefield.setReactionFieldDielectric(rf_dielectric.val)

    else:
        print("Cannot interpret the cutoff method from \"%s\"" % cutoff_method.val, file=sys.stderr)

    forcefield.setSpace(Cartesian())

    if fast_sim.val:
        forcefield.setSwitchingFunction( HarmonicSwitchingFunction(7.5 * angstrom, 7.5 * angstrom,
                                                                   7.5 * angstrom, 7.5 * angstrom) )
    else:
        forcefield.setSwitchingFunction( HarmonicSwitchingFunction(coul_cutoff.val,coul_cutoff.val,
                                                                   lj_cutoff.val,lj_cutoff.val) )

    return forcefield


def setFakeGridProperties(forcefield):
    if fast_sim.val:
        forcefield.setSwitchingFunction( HarmonicSwitchingFunction(7.5*angstrom, 7.5*angstrom,
                                                                   7.5*angstrom, 7.5*angstrom) )
    else:
        forcefield.setSwitchingFunction( HarmonicSwitchingFunction(coul_cutoff.val,coul_cutoff.val,
                                                                   lj_cutoff.val,lj_cutoff.val) )
    forcefield.setSpace(Cartesian())

    return forcefield


def setGridProperties(forcefield, extra_buffer=0*angstrom):
    forcefield.setGridSpacing(grid_spacing.val)
    forcefield.setBuffer(grid_buffer.val + extra_buffer)

    if fast_sim.val:
        forcefield.setLJCutoff(7.5*angstrom)
        forcefield.setCoulombCutoff(7.5*angstrom)
    else:
        forcefield.setLJCutoff(lj_cutoff.val)
        forcefield.setCoulombCutoff(coul_cutoff.val)

    return forcefield


def setSoftCoreProperties(forcefield):
    forcefield.setCoulombPower(coul_power.val)
    forcefield.setShiftDelta(shift_delta.val)

    return forcefield


def setCLJFuncProperties(cljfunc):
    cljfunc.setSpace(Cartesian())

    if fast_sim.val:
        cljfunc.setCoulombCutoff(7.5*angstrom)
        cljfunc.setLJCutoff(7.5*angstrom)
    else:
        cljfunc.setCoulombCutoff(coul_cutoff.val)
        cljfunc.setLJCutoff(lj_cutoff.val)

    cljfunc.setArithmeticCombiningRules( True )

    return cljfunc

def getInterCLJFunction():
    if cutoff_method.val.find("shift electrostatics") != -1:
        cljfunc = CLJShiftFunction()

    elif cutoff_method.val.find("reaction field") != -1:
        cljfunc = CLJRFFunction()
        cljfunc.setDielectric(rf_dielectric.val)

    else:
        print("Cannot interpret the cutoff method from \"%s\"" % cutoff_method.val, file=sys.stderr)

    return setCLJFuncProperties(cljfunc)

def getSoftInterCLJFunction():
    if cutoff_method.val.find("shift electrostatics") != -1:
        cljfunc = CLJSoftShiftFunction()

    elif cutoff_method.val.find("reaction field") != -1:
        cljfunc = CLJSoftRFFunction()
        cljfunc.setDielectric(rf_dielectric.val)

    else:
        print("Cannot interpret the cutoff method from \"%s\"" % cutoff_method.val, file=sys.stderr)

    cljfunc.setAlpha(0.0)
    cljfunc.setShiftDelta(shift_delta.val)
    cljfunc.setCoulombPower(coul_power.val)

    return setCLJFuncProperties(cljfunc)

def getIntraCLJFunction():
    if cutoff_method.val.find("shift electrostatics") != -1:
        cljfunc = CLJIntraShiftFunction()

    elif cutoff_method.val.find("reaction field") != -1:
        cljfunc = CLJIntraRFFunction()
        cljfunc.setDielectric(rf_dielectric.val)

    else:
        print("Cannot interpret the cutoff method from \"%s\"" % cutoff_method.val, file=sys.stderr)

    return setCLJFuncProperties(cljfunc)

def getSoftIntraCLJFunction():
    if cutoff_method.val.find("shift electrostatics") != -1:
        cljfunc = CLJSoftIntraShiftFunction()

    elif cutoff_method.val.find("reaction field") != -1:
        cljfunc = CLJSoftIntraRFFunction()
        cljfunc.setDielectric(rf_dielectric.val)

    else:
        print("Cannot interpret the cutoff method from \"%s\"" % cutoff_method.val, file=sys.stderr)

    cljfunc.setAlpha(0.0)
    cljfunc.setShiftDelta(shift_delta.val)
    cljfunc.setCoulombPower(coul_power.val)

    return setCLJFuncProperties(cljfunc)

def setNewGridProperties(forcefield, extra_buffer=0*angstrom):
    if disable_grid.val:
        forcefield.disableGrid()
    else:
        forcefield.enableGrid()
        forcefield.setGridSpacing(grid_spacing.val)
        forcefield.setGridBuffer(grid_buffer.val + extra_buffer)

    return forcefield


def getAtomNearCOG( molecule ):
    """Find the atom that is closest to the center of geometry of the passed molecule"""

    mol_centre = molecule.evaluate().center()
    mindist = 99999.0

    for x in range(0, molecule.nAtoms()):
        atom = molecule.atoms()[x]
        at_coords = atom.property('coordinates')
        dist = Vector().distance2(at_coords, mol_centre)
        if dist < mindist:
            mindist = dist
            nearest_atom = atom

    return nearest_atom


def createWSRCMoves(system):
    # pull out all of the molecule groups for the mobile parts of the system
    mobile_solvent = system[MGName("mobile_solvent")]
    mobile_sidechains = system[MGName("mobile_sidechains")]
    mobile_backbones = system[MGName("mobile_backbones")]
    mobile_solutes = system[MGName("mobile_solutes")]
    mobile_ligand = system[MGName("mobile_ligand")]
    mobile_swap = system[MGName("mobile_swap_water")]

    print("Creating the Monte Carlo moves to sample the WSRC system...")

    # create the global set of moves that will be applied to
    # the system
    moves = WeightedMoves()

    if fast_sim.val:
        # we will only move the water
        max_water_translation = 0.15 * angstroms
        max_water_rotation = 15 * degrees

        if mobile_swap.nViews() > 0:
            rb_moves = RigidBodyMC(mobile_swap)
            rb_moves.setMaximumTranslation(max_water_translation)
            rb_moves.setMaximumRotation(max_water_rotation)

            if use_reflect_volume.val:
                rb_moves.setReflectionVolume( mobile_ligand[MolIdx(0)], reflect_volume_radius.val )

            moves.add(rb_moves, 4 * mobile_swap.nViews())


        if mobile_solvent.nViews() > 0:
            rb_moves = RigidBodyMC(mobile_solvent)
            rb_moves.setMaximumTranslation(max_water_translation)
            rb_moves.setMaximumRotation(max_water_rotation)

            if system.containsProperty("reflection sphere radius"):
                reflection_radius = float(str(system.property("reflection sphere radius"))) * angstroms
                reflection_center = system.property("reflection center").toVector()[0]
                rb_moves.setReflectionSphere(reflection_center, reflection_radius)

            moves.add(rb_moves, 4 * mobile_solvent.nViews())

        moves.setTemperature(temperature.val)

        seed = random_seed.val

        if seed is None:
            seed = RanGenerator().randInt(100000,1000000)
            print("Using generated random number seed %d" % seed)
        else:
            print("Using supplied random number seed %d" % seed)

        moves.setGenerator( RanGenerator(seed) )

        return moves
    
    # we are performing a normal simulation, moving everything
  
    # create zmatrix moves to move the protein sidechains
    if mobile_sidechains.nViews() > 0:
        sc_moves = ZMatMove(mobile_sidechains)
        moves.add( sc_moves, mobile_sidechains.nViews() )

    if mobile_backbones.nViews() > 0 and move_backbone.val:
        bb_moves = RigidBodyMC(mobile_backbones)
        bb_moves.setCenterOfRotation( GetCOGPoint( AtomName("CA", CaseInsensitive),
                                                   AtomName("N", CaseInsensitive) ) )

        bb_moves.setMaximumTranslation(0.030*angstrom)
        bb_moves.setMaximumRotation(1.0*degrees)
        moves.add( bb_moves, mobile_backbones.nViews() )

    if not use_fixed_ligand.val:
        if mobile_ligand.nViews() > 0:
            scale_moves = 10

            # get the amount to translate and rotate from the ligand's flexibility object
            flex = mobile_ligand.moleculeAt(0)[0].molecule().property("flexibility")

            if use_rot_trans_ligand.val:
                if (flex.translation().value() != 0 or flex.rotation().value() != 0):
                    rb_moves = RigidBodyMC(mobile_ligand)
                    rb_moves.setMaximumTranslation(flex.translation())
                    rb_moves.setMaximumRotation(flex.rotation())
                    rb_moves.setCenterOfRotation(GetCOMPoint())

                    # the ligand is not allowed to move away from its original position,
                    # as we don't want to sample "unbound" states
                    if not ligand_reflection_radius.val is None:
                        rb_moves.setReflectionSphere(mobile_ligand.moleculeAt(0)[0].molecule().evaluate().centerOfMass(), 
                                                     ligand_reflection_radius.val)

                    scale_moves = scale_moves / 2
                    moves.add( rb_moves, scale_moves * mobile_ligand.nViews() )

            intra_moves = InternalMove(mobile_ligand)
            intra_moves.setCenterOfMolecule(GetCOMPoint())
            moves.add( intra_moves, scale_moves * mobile_ligand.nViews() )

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

        intra_moves = InternalMove(mobile_solutes)
        moves.add(intra_moves, 4 * mobile_solutes.nViews())

    max_water_translation = 0.15 * angstroms
    max_water_rotation = 15 * degrees

    if mobile_swap.nViews() > 0:
        rb_moves = RigidBodyMC(mobile_swap)
        rb_moves.setMaximumTranslation(max_water_translation)
        rb_moves.setMaximumRotation(max_water_rotation)

        if use_reflect_volume.val:
            rb_moves.setReflectionVolume( mobile_ligand[MolIdx(0)], reflect_volume_radius.val )

        moves.add(rb_moves, 4 * mobile_swap.nViews())

    if mobile_solvent.nViews() > 0:
        rb_moves = RigidBodyMC(mobile_solvent)
        rb_moves.setMaximumTranslation(max_water_translation)
        rb_moves.setMaximumRotation(max_water_rotation)

        if system.containsProperty("reflection sphere radius"):
            reflection_radius = float(str(system.property("reflection sphere radius"))) * angstroms
            reflection_center = system.property("reflection center").toVector()[0]
            rb_moves.setReflectionSphere(reflection_center, reflection_radius)

        moves.add(rb_moves, 4 * mobile_solvent.nViews())

    moves.setTemperature(temperature.val)

    seed = random_seed.val

    if seed is None:
        seed = RanGenerator().randInt(100000,1000000)
        print("Using generated random number seed %d" % seed)
    else:
        print("Using supplied random number seed %d" % seed)
    
    moves.setGenerator( RanGenerator(seed) )

    return moves    


def renumberMolecules(molgroup):
    newgroup = MoleculeGroup(molgroup.name().value())
    for molnum in molgroup.molNums():
        mol = molgroup[molnum][0]
        newmol = mol.molecule().edit().renumber().commit()
        newgroup.add( ViewsOfMol(newmol,mol.selections()) )

    return newgroup

def pruneLambda(lamvals):
    """Half the passed number of lambda values"""

    n = len(lamvals)

    if n % 2 == 0:
        mid = int(n/2)
        return lamvals[0:mid:2] + lamvals[mid+1::2]
    else:
        return lamvals[0::2]

def getLambdaValues():
    """Return the lambda values to use for the simulation. Lambda scale from 0 to 1
       and will include the discharging and charging steps. The values are set such
       that, if no charging is used, then return lambda_values.val (the lambda values
       set by the user). If discharging / charging is used, then discharging is from
       lambda = 0-0.25, swapping from 0.25-0.75 and recharging from 0.75-1.0"""

    if uncharge_ligand.val:
        lamvals = []

        charge_lams = copy.deepcopy( uncharge_lambda_values.val )
        charge_lams.sort()

        swap_lams = copy.deepcopy( lambda_values.val )
        swap_lams.sort()

        for lam in charge_lams:
            if lam >= 0.0 and lam <= 1.0:
                lamvals.append( 0.25 * lam )

        for lam in swap_lams:
            if lam >= 0.0 and lam <= 1.0:
                lamvals.append( 0.25 + (0.5*lam) )

        charge_lams.reverse()

        for lam in charge_lams:
            if lam >= 0.0 and lam <= 1.0:
                lamvals.append( 0.75 + (0.25*(1.0-lam)) )

        if fast_sim.val:
            lamvals = pruneLambda(lamvals)

        return lamvals

    else:
        swap_lams = copy.deepcopy( lambda_values.val )
        swap_lams.sort()

        lamvals = []

        for lam in swap_lams:
            if lam >= 0.0 and lam <= 1.0:
                lamvals.append(lam)

        if fast_sim.val:
            lamvals = pruneLambda(lamvals)

        return lamvals


def printEnergies(nrgs, FILE):
    """This function prints all of the energies in 'nrgs' to the file 'FILE'"""
    keys = list(nrgs.keys())
    keys.sort()

    for key in keys:
        FILE.write("%s  ==  %s kcal mol-1\n" % (key, nrgs[key]))


def mergeSystems(protein_system, water_system, ligand_mol):
    
    print("Merging the protein box and water box to create the WSRC system...")

    system = System("WSRC system")

    if protein_system.containsProperty("reflection center"):
        prot_reflection_center = protein_system.property("reflection center").toVector()[0]
        prot_reflection_radius = float(str(protein_system.property("reflection sphere radius")))

        wat_reflection_center = water_system.property("reflection center").toVector()[0]
        wat_reflection_radius = float(str(water_system.property("reflection sphere radius")))

        if prot_reflection_center != wat_reflection_center or \
           prot_reflection_radius != wat_reflection_radius:
            print("Disagreement of the reflection sphere in the protein and water boxes!")
            print("Protein: %s and %s    Water: %s and %s" % \
                    (prot_reflection_center,prot_reflection_radius,
                     wat_reflection_center,wat_reflection_radius))

            sys.exit(-1)

        system.setProperty("reflection center", AtomCoords(CoordGroup(1,prot_reflection_center)))
        system.setProperty("reflection sphere radius", VariantProperty(prot_reflection_radius))

    if protein_system.containsProperty("average solute translation delta"):
        system.setProperty("average solute translation delta", \
                       protein_system.property("average solute translation delta"))

    if protein_system.containsProperty("average solute rotation delta"):
        system.setProperty("average solute rotation delta", \
                       protein_system.property("average solute rotation delta"))

    # create a molecule group for the ligand
    ligand_group = MoleculeGroup("ligand")
    ligand_group.add(ligand_mol)

    bound_leg = MoleculeGroup("bound_leg")
    free_leg = MoleculeGroup("free_leg")

    bound_leg.add(ligand_mol)
    free_leg.add(ligand_mol)

    # pull out the groups that we want from the two systems

    # create a group to hold all of the mobile water molecules in the free leg
    mobile_free_water_group = MoleculeGroup("mobile_free")
    water_mol = None
    if MGName("mobile_solvents") in water_system.mgNames():
        mols = water_system[MGName("mobile_solvents")].molecules()
        for molnum in mols.molNums():
            water_mol = mols[molnum][0].molecule().edit().renumber().commit()
            for j in range(0,water_mol.nResidues()):
                water_mol = water_mol.residue( ResIdx(j) ).edit() \
                                           .setProperty( PDB.parameters().pdbResidueName(), "FWT" ) \
                                           .commit().molecule()

            mobile_free_water_group.add(water_mol)

    # create a group to hold all of the fixed water molecules in the free leg
    fixed_free_water_group = MoleculeGroup("fixed_free")
    if MGName("fixed_molecules") in water_system.mgNames():
        mols = water_system[MGName("fixed_molecules")].molecules()
        for molnum in mols.molNums():
            fixed_free_water_group.add( mols[molnum][0].molecule().edit().renumber().commit() )

    # create a group to hold all of the fixed molecules in the bound leg
    fixed_bound_group = MoleculeGroup("fixed_bound")
    if MGName("fixed_molecules") in protein_system.mgNames():
        fixed_bound_group.add( protein_system[ MGName("fixed_molecules") ] )

    if save_pdb.val:
        # write a PDB of the fixed atoms in the bound and free legs
        if not os.path.exists(outdir.val):
            os.makedirs(outdir.val)

        PDB().write(fixed_bound_group, "%s/bound_fixed.pdb" % outdir.val)
        PDB().write(fixed_free_water_group, "%s/free_fixed.pdb" % outdir.val)

    # create a group to hold all of the mobile solute molecules in the bound leg
    mobile_bound_solutes_group = MoleculeGroup("mobile_bound_solutes")
    if MGName("mobile_solutes") in protein_system.mgNames():
        mobile_bound_solutes_group.add( protein_system[MGName("mobile_solutes")] )
        mobile_bound_solutes_group.remove(ligand_mol)
        if mobile_bound_solutes_group.nMolecules() > 0:
            bound_leg.add(mobile_bound_solutes_group)
    
    # create a group to hold all of the mobile solvent molecules in the bound leg
    mobile_bound_solvents_group = MoleculeGroup("mobile_bound_solvents")
    mobile_bound_water_group = MoleculeGroup("mobile_bound_water")
    if MGName("mobile_solvents") in protein_system.mgNames():
        mols = protein_system[MGName("mobile_solvents")]
        for molnum in mols.molNums():
            solvent_mol = mols[molnum][0].molecule()

            try:
                # this is a water molecule if we can swap the coordinates with the 
                # water molecule from teh water box
                water_mol.edit().setProperty("coordinates", \
                                     solvent_mol.property("coordinates"))
                
                for j in range(0,solvent_mol.nResidues()):
                    solvent_mol = solvent_mol.residue( ResIdx(j) ).edit() \
                                             .setProperty( PDB.parameters().pdbResidueName(), "BWT" ) \
                                             .commit().molecule()
                
                mobile_bound_solvents_group.add(solvent_mol)
                mobile_bound_water_group.add(solvent_mol)
            except:
                # the test molecule is not compatible, so it is not 
                # compatible with the water in the water box
                mobile_bound_solvents_group.add(solvent_mol)

        print("The number of bound leg mobile solvent molecules is %d." % mobile_bound_solvents_group.nMolecules())
        print("The number of these which are compatible water molecules is %d." % mobile_bound_water_group.nMolecules())

    # create the groups to hold all of the protein molecules. We will use "extract" to 
    # pull out only those protein atoms that are in the mobile region
    bound_protein_intra_group = MoleculeGroup("bound_protein_intra_group")
    mobile_bound_proteins_group = MoleculeGroup("mobile_bound_proteins")
    mobile_bound_protein_sidechains_group = MoleculeGroup("mobile_bound_protein_sidechains")
    mobile_bound_protein_backbones_group = MoleculeGroup("mobile_bound_protein_backbones")

    if MGName("protein_sidechains") in protein_system.mgNames() or \
       MGName("protein_backbones") in protein_system.mgNames():

        all_proteins = Molecules()

        try:
            protein_sidechains = protein_system[MGName("protein_sidechains")]
            all_proteins.add(protein_sidechains.molecules())
        except:
            protein_sidechains = MoleculeGroup()

        try:
            protein_backbones = protein_system[MGName("protein_backbones")]
            all_proteins.add(protein_backbones.molecules())
        except:
            protein_backbones = MoleculeGroup()

        try:
            boundary_molecules = protein_system[MGName("boundary_molecules")]
            all_proteins.add(boundary_molecules.molecules())
        except:
            boundary_molecules = MoleculeGroup()

        for molnum in all_proteins.molNums():
            protein_mol = Molecule.join(all_proteins[molnum])
            
            if protein_mol.selectedAll():
                bound_protein_intra_group.add(protein_mol)
                bound_leg.add(protein_mol)

                mobile_protein = None                

                try:
                    mobile_protein = protein_sidechains[molnum]
                    mobile_bound_protein_sidechains_group.add( mobile_protein )
                except:
                    pass

                try:
                    if mobile_protein is None:
                        mobile_protein = protein_backbones[molnum]
                        mobile_bound_protein_backbones_group.add( mobile_protein )
                    else:
                        mobile_protein.add( protein_backbones[molnum].selection() )
                        mobile_bound_protein_backbones_group.add( protein_backbones[molnum] )
                except:
                    pass

                if not (mobile_protein is None):
                    mobile_bound_proteins_group.add( Molecule.join(mobile_protein) )

            else:
                # only some of the atoms have been selected. We will extract
                # the mobile atoms and will then update all of the other selections
                print("Extracting the mobile atoms of protein %s" % protein_mol)
                new_protein_mol = protein_mol.extract()
                print("Extracted %d mobile atoms from %d total atoms..." % \
                                        (new_protein_mol.nAtoms(), protein_mol.molecule().nAtoms()))

                bound_protein_intra_group.add(new_protein_mol)
                bound_leg.add( new_protein_mol )

                mobile_protein_view = new_protein_mol.selection()
                mobile_protein_view = mobile_protein_view.selectNone()

                try:
                    sidechains = protein_sidechains[molnum]

                    for i in range(0,sidechains.nViews()):
                        view = new_protein_mol.selection()
                        view = view.selectNone()

                        for atomid in sidechains.viewAt(i).selectedAtoms():
                            atom = protein_mol.atom(atomid)
                            resatomid = ResAtomID( atom.residue().number(), atom.name() )
                            view = view.select( resatomid )
                            mobile_protein_view = mobile_protein_view.select( resatomid )

                        if view.nSelected() > 0:
                            mobile_bound_protein_sidechains_group.add( PartialMolecule(new_protein_mol, view) )
                except:
                    pass

                try:
                    backbones = protein_backbones[molnum]

                    for i in range(0,backbones.nViews()):
                        view = new_protein_mol.selection()
                        view = view.selectNone()

                        for atomid in backbones.viewAt(i).selectedAtoms():
                            atom = protein_mol.atom(atomid)
                            resatomid = ResAtomID( atom.residue().number(), atom.name() )
                            view = view.select( resatomid )
                            mobile_protein_view = mobile_protein_view.select( resatomid )

                        if view.nSelected() > 0:
                            mobile_bound_protein_backbones_group.add( PartialMolecule(new_protein_mol, view) )
                except:
                    pass

                if mobile_protein_view.nSelected() > 0:
                    mobile_bound_proteins_group.add( PartialMolecule(new_protein_mol, mobile_protein_view) )

    # finished adding in all of the protein groups

    use_identity_constraint = not use_reflect_volume.val

    if use_identity_constraint:
        # get the identity points for the ligand
        print("\nObtaining the identity points...")

        if identity_atoms.val is None:
            print("Auto-identifying the identity atoms...")
            identity_points = getIdentityPoints(ligand_mol)
        else:
            identity_points = []
            for identity_atom in identity_atoms.val:
                identity_points.append( ligand_mol.atom( AtomName(identity_atom) ) )

        print("Using identity points:")
        print(identity_points)

        if use_fixed_points.val:
            print("\nUsing fixed identity points...")
            fixed_points = []
            for point in identity_points:
                fixed_points.append( point.property("coordinates") )
        
            identity_points = fixed_points
            print(identity_points) 

        print("\nIdentifying the swap-water cluster...")
        swap_water_group = MoleculeGroup("swap water")
        mobile_free_water_group = IdentityConstraint.constrain( mobile_free_water_group, identity_points )

        # Rename the residues of the swap solvent so that they are easy
        # to find in the output PDBs
        for i in range(0,len(identity_points)):
            swap_water_mol = mobile_free_water_group.moleculeAt(i)[0].molecule()
            mobile_free_water_group.remove(swap_water_mol)

            for j in range(0,swap_water_mol.nResidues()):
                swap_water_mol = swap_water_mol.residue( ResIdx(j) ).edit() \
                                               .setProperty( PDB.parameters().pdbResidueName(), "SWP" ) \
                                               .commit().molecule()

            swap_water_group.add(swap_water_mol)

        print("found %d molecules that are now part of the swap water cluster" % swap_water_group.nMolecules())
        PDB().write(swap_water_group.molecules(), "swapcluster00.pdb")

        # now equilibrate the swap water cluster, if requested
        if n_equil_swap.val:
            move = RigidBodyMC(swap_water_group)
            move.setMaximumTranslation(0.15*angstrom)
            move.setMaximumRotation(15*degrees)

            # use the same random number seed so that the swap water cluster is reproducible
            move.setGenerator( RanGenerator(4039251) )

            equil_system = System()

            ff = InterFF("interff")
            ff.setCLJFunction( getInterCLJFunction() )
            ff = setNewGridProperties(ff)

            ff.add(swap_water_group)

            fixed_mols = bound_leg.molecules()
            fixed_mols.remove(ligand_mol.number())

            ff.addFixedAtoms(fixed_mols)

            equil_system.add(swap_water_group)
            equil_system.add(ff)

            n = n_equil_swap.val
            print("Equilibrating the swap water cluster (moves = %d)" % n)

            if n > 10:
                for i in range(1,11):
                    move.move(equil_system, int(n / 10), False)
            else:
                move.move(equil_system, n, False)

            swap_water_group = equil_system[ swap_water_group.name() ]

            PDB().write(swap_water_group, "swapcluster01.pdb")

            print("Complete. Equilibrated water molecules in file 'swapcluster01.pdb'")

        if use_water_points.val:
            # use identity points that are fixed in space based on the current position
            # of the center of each molecule of the swap water cluster. This prevents the swap cluster
            # from changing shape too much during the calculation
            print("\nUsing identity points based on fixed positions in space from the swap water cluster...")
            
            fixed_points = []

            for i in range(0,swap_water_group.nMolecules()):
                sw = swap_water_group.molecule(MolIdx(i))[0].molecule()
                fixed_points.append( VectorPoint(sw.evaluate().centerOfMass()) )

            print("Using fixed identity points %s" % fixed_points)
            identity_points = fixed_points

    else:
        # we will be using the reflection volume to get the swap water cluster
        swap_water_group = MoleculeGroup("swap water")
        move = RigidBodyMC(swap_water_group)
        move.setMaximumTranslation(0.15*angstrom)
        move.setMaximumRotation(15*degrees)
        move.setReflectionVolume(ligand_mol, reflect_volume_radius.val)

        # find the swap water cluster by looking at overlapping waters
        swap_waters = move.extract(mobile_free_water_group.molecules(), reflect_volume_buffer.val)

        print("Swap water cluster based on the %d water molecules overlapping with the ligand." % swap_waters.nMolecules())

        for molnum in swap_waters.molNums():
            swap_water = swap_waters[molnum][0].molecule()

            for j in range(0,swap_water.nResidues()):
                swap_water = swap_water.residue( ResIdx(j) ).edit() \
                                       .setProperty( PDB.parameters().pdbResidueName(), "SWP" ) \
                                       .commit().molecule()

            swap_water_group.add( swap_water )
            mobile_free_water_group.remove(swap_water.number())

        # now equilibrate the swap cluster, if requested
        if n_equil_reflect.val:
             equil_system = System()
             ff = InterCLJFF("ff")
             ff.add(swap_water_group)
             liggroup = MoleculeGroup("liggroup")
             liggroup.add(ligand_mol)
             equil_system.add(ff)
             equil_system.add(swap_water_group)
             equil_system.add(liggroup)

             PDB().write(equil_system.molecules(), "swapcluster00.pdb")

             n = n_equil_reflect.val
             if n > 10:
                 for i in range(1,11):
                     move.move(equil_system, int(n / 10), False)
                     PDB().write(equil_system.molecules(), "swapcluster%02d.pdb" % i)
             else:
                 move.move(equil_system, n, False)
                 PDB().write(equil_system.molecules(), "swapcluster01.pdb")

             swap_water_group = equil_system[ swap_water_group.name() ]

    bound_leg.add(swap_water_group)
    bound_leg.add(mobile_bound_solvents_group)

    free_leg.add(swap_water_group)
    free_leg.add(mobile_free_water_group)

    system.add(bound_leg)
    system.add(free_leg)

    # now add in the forcefields for the system...
    print("Creating the forcefields for the WSRC system...")

    # first, group together the molecules grouped above into convenient
    # groups for the forcefields

    # group holding just the ligand
    ligand_mols = ligand_group.molecules()

    # group holding just the swap water cluster
    swap_water_mols = swap_water_group.molecules()

    # group holding all of the mobile atoms in the bound leg
    mobile_bound_mols = mobile_bound_solvents_group.molecules()
    mobile_bound_mols.add( mobile_bound_solutes_group.molecules() )
    mobile_bound_mols.add( bound_protein_intra_group.molecules() )

    # group holding all of the mobile atoms in the bound leg, excluding the 
    # buffer atoms that are fixed, but bonded to mobile atoms
    mobile_buffered_bound_mols = mobile_bound_solvents_group.molecules()
    mobile_buffered_bound_mols.add( mobile_bound_solutes_group.molecules() )
    mobile_buffered_bound_mols.add( mobile_bound_proteins_group.molecules() )

    # group holding all of the mobile water molecules in the free leg
    mobile_free_mols = mobile_free_water_group.molecules()

    # group holding all of the fixed water molecules in the free leg
    fixed_free_group = fixed_free_water_group

    # group holding all of the protein molecules that need intramolecular terms calculated
    bound_protein_intra_mols = bound_protein_intra_group.molecules()

    # group holding all of the solute molecules that nede intramolecular terms calculated
    bound_solute_intra_mols = mobile_bound_solutes_group.molecules()

    ###
    ### INTRA-ENERGY OF THE LIGAND AND CLUSTER
    ###
    
    # intramolecular energy of the ligand
    if use_oldff.val:
        ligand_intraclj = IntraCLJFF("ligand:intraclj")
        ligand_intraclj = setCLJProperties(ligand_intraclj)
        ligand_intraclj.add(ligand_mols)

        ligand_intraff = InternalFF("ligand:intra")
        ligand_intraff.setUse14Calculation(False)
        ligand_intraff.add(ligand_mols)
    else:
        print("Using the NEW PARALLEL FORCEFIELDS :-)")
        # intramolecular energy of the ligands
        ligand_intraclj = IntraFF("ligand:intraclj")
        ligand_intraclj.setCLJFunction( getIntraCLJFunction() )
        ligand_intraclj.add(ligand_mols)

        ligand_intraff = InternalFF("ligand:intra")
        ligand_intraff.setUse14Calculation(True)
        ligand_intraff.add(ligand_mols)

    # intramolecular energy of the swap water cluster
    if use_oldff.val:
        swap_interclj = InterSoftCLJFF("swap:interclj")
        swap_interclj = setCLJProperties(swap_interclj)
        swap_interclj = setSoftCoreProperties(swap_interclj)
        swap_interclj.add(swap_water_mols)
    else:
        swap_interclj = InterFF("swap:interclj")
        swap_interclj.setCLJFunction( getSoftInterCLJFunction() )
        swap_interclj.setCLJFunction( "f", getSoftInterCLJFunction() )
        swap_interclj.setCLJFunction( "b", getSoftInterCLJFunction() )
        swap_interclj.setCLJFunction( "next", getSoftInterCLJFunction() )
        swap_interclj.setCLJFunction( "prev", getSoftInterCLJFunction() )
        swap_interclj.add(swap_water_mols)

    ###
    ### FORCEFIELDS INVOLVING THE LIGAND/CLUSTER BOUND LEG
    ###

    # forcefield holding the energy between the ligand and the mobile atoms in the
    # bound leg
    if use_oldff.val:
        bound_ligand_mobile = InterGroupSoftCLJFF("bound:ligand-mobile")
        bound_ligand_mobile = setCLJProperties(bound_ligand_mobile)
        bound_ligand_mobile = setSoftCoreProperties(bound_ligand_mobile)

        bound_ligand_mobile.add(ligand_mols, MGIdx(0))
        bound_ligand_mobile.add(mobile_bound_mols, MGIdx(1))

        # Whether or not to disable the grid and calculate all energies atomisticly
        if disable_grid:
            # we need to renumber all of the fixed molecules so that they don't clash
            # with the mobile molecules
            print("Renumbering fixed molecules...")
            fixed_bound_group = renumberMolecules(fixed_bound_group)
            fixed_free_group = renumberMolecules(fixed_free_group)

        # forcefield holding the energy between the ligand and the fixed atoms in the bound leg
        if disable_grid:
            bound_ligand_fixed = InterGroupCLJFF("bound:ligand-fixed")
            bound_ligand_fixed = setCLJProperties(bound_ligand_fixed)
            bound_ligand_fixed = setFakeGridProperties(bound_ligand_fixed)

            bound_ligand_fixed.add(ligand_mols, MGIdx(0))
            bound_ligand_fixed.add(fixed_bound_group, MGIdx(1))
        else:
            bound_ligand_fixed = GridFF2("bound:ligand-fixed")
            bound_ligand_fixed = setCLJProperties(bound_ligand_fixed)
            bound_ligand_fixed = setGridProperties(bound_ligand_fixed)

            bound_ligand_fixed.add(ligand_mols, MGIdx(0))
            bound_ligand_fixed.addFixedAtoms( fixed_bound_group )
    else:
        bound_ligand_mobile = InterGroupFF("bound:ligand-mobile")
        bound_ligand_mobile.setCLJFunction( getSoftInterCLJFunction() )
        bound_ligand_mobile.setCLJFunction( "f", getSoftInterCLJFunction() )
        bound_ligand_mobile.setCLJFunction( "b", getSoftInterCLJFunction() )
        bound_ligand_mobile.setCLJFunction( "next", getSoftInterCLJFunction() )
        bound_ligand_mobile.setCLJFunction( "prev", getSoftInterCLJFunction() )

        bound_ligand_mobile.add(ligand_mols, MGIdx(0))
        bound_ligand_mobile.add(mobile_bound_mols, MGIdx(1))

        bound_ligand_fixed = InterGroupFF("bound:ligand0-fixed")
        bound_ligand_fixed.setCLJFunction( getInterCLJFunction() )
        bound_ligand_fixed = setNewGridProperties(bound_ligand_fixed)

        bound_ligand_fixed.add(ligand_mols, MGIdx(0))
        bound_ligand_fixed.addFixedAtoms(fixed_bound_group.molecules())

    # forcefield holding the energy between the swap water cluster and the mobile
    # atoms in the bound leg
    if use_oldff.val:
        bound_swap_mobile = InterGroupSoftCLJFF("bound:swap-mobile")
        bound_swap_mobile = setCLJProperties(bound_swap_mobile)
        bound_swap_mobile = setSoftCoreProperties(bound_swap_mobile)

        bound_swap_mobile.add(swap_water_mols, MGIdx(0))
        bound_swap_mobile.add(mobile_bound_mols, MGIdx(1))

        # forcefield holding the energy between the swap water cluster and the
        # fixed atoms in the bound leg
        if disable_grid:
            bound_swap_fixed = InterGroupCLJFF("bound:swap-fixed")
            bound_swap_fixed = setCLJProperties(bound_swap_fixed)
            bound_swap_fixed = setFakeGridProperties(bound_swap_fixed)

            bound_swap_fixed.add(swap_water_mols, MGIdx(0))
            bound_swap_fixed.add( fixed_bound_group, MGIdx(1) )
        else:
            bound_swap_fixed = GridFF2("bound:swap-fixed")
            bound_swap_fixed = setCLJProperties(bound_swap_fixed)
            # The swap water cluster is more mobile, so needs a bigger grid buffer
            bound_swap_fixed = setGridProperties(bound_swap_fixed, 6*angstrom)

            bound_swap_fixed.add(swap_water_mols, MGIdx(0))
            bound_swap_fixed.addFixedAtoms(fixed_bound_group)
    else:
        bound_swap_mobile = InterGroupFF("bound:swap-mobile")
        bound_swap_mobile.setCLJFunction( getSoftInterCLJFunction() )
        bound_swap_mobile.setCLJFunction( "f", getSoftInterCLJFunction() )
        bound_swap_mobile.setCLJFunction( "b", getSoftInterCLJFunction() )
        bound_swap_mobile.setCLJFunction( "next", getSoftInterCLJFunction() )
        bound_swap_mobile.setCLJFunction( "prev", getSoftInterCLJFunction() )

        bound_swap_mobile.add(swap_water_mols, MGIdx(0))
        bound_swap_mobile.add(mobile_bound_mols, MGIdx(1))
 
        bound_swap_fixed = InterGroupFF("bound:swap-fixed")
        bound_swap_fixed.setCLJFunction( getInterCLJFunction() )
        bound_swap_fixed = setNewGridProperties(bound_swap_fixed)

        bound_swap_fixed.add(swap_water_mols, MGIdx(0))
        bound_swap_fixed.addFixedAtoms(fixed_bound_group.molecules())

    ###
    ### FORCEFIELDS INVOLVING THE LIGAND/CLUSTER FREE LEG
    ###

    # forcefield holding the energy between the ligand and the mobile atoms
    # in the free leg
    if use_oldff.val:
        free_ligand_mobile = InterGroupSoftCLJFF("free:ligand-mobile")
        free_ligand_mobile = setCLJProperties(free_ligand_mobile)
        free_ligand_mobile = setSoftCoreProperties(free_ligand_mobile)

        free_ligand_mobile.add(ligand_mols, MGIdx(0))
        free_ligand_mobile.add(mobile_free_mols, MGIdx(1))

        # forcefield holding the energy between the ligand and the fixed atoms
        # in the free leg
        if disable_grid.val:
            free_ligand_fixed = InterGroupCLJFF("free:ligand_fixed")
            free_ligand_fixed = setCLJProperties(free_ligand_fixed)
            free_ligand_fixed = setFakeGridProperties(free_ligand_fixed)

            free_ligand_fixed.add(ligand_mols, MGIdx(0))
            free_ligand_fixed.add(fixed_free_group, MGIdx(1))
        else:
            free_ligand_fixed = GridFF2("free:ligand-fixed")
            free_ligand_fixed = setCLJProperties(free_ligand_fixed)
            free_ligand_fixed = setGridProperties(free_ligand_fixed)

            free_ligand_fixed.add(ligand_mols, MGIdx(0))
            free_ligand_fixed.addFixedAtoms(fixed_free_group)
    else:
        free_ligand_mobile = InterGroupFF("free:ligand-mobile")
        free_ligand_mobile.setCLJFunction( getSoftInterCLJFunction() )
        free_ligand_mobile.setCLJFunction( "f", getSoftInterCLJFunction() )
        free_ligand_mobile.setCLJFunction( "b", getSoftInterCLJFunction() )
        free_ligand_mobile.setCLJFunction( "next", getSoftInterCLJFunction() )
        free_ligand_mobile.setCLJFunction( "prev", getSoftInterCLJFunction() )

        free_ligand_mobile.add(ligand_mols, MGIdx(0))
        free_ligand_mobile.add(mobile_free_mols, MGIdx(1))

        free_ligand_fixed = InterGroupFF("free:ligand-fixed")
        free_ligand_fixed.setCLJFunction( getInterCLJFunction() )
        free_ligand_fixed = setNewGridProperties(free_ligand_fixed)

        free_ligand_fixed.add(ligand_mols, MGIdx(0))
        free_ligand_fixed.addFixedAtoms(fixed_free_group.molecules())

    # forcefield holding the energy between the swap water cluster and the
    # mobile atoms of the free leg
    if use_oldff.val:
        free_swap_mobile = InterGroupSoftCLJFF("free:swap-mobile")
        free_swap_mobile = setCLJProperties(free_swap_mobile)
        free_swap_mobile = setSoftCoreProperties(free_swap_mobile)

        free_swap_mobile.add(swap_water_mols, MGIdx(0))
        free_swap_mobile.add(mobile_free_mols, MGIdx(1))

        # forcefield holding the energy between the swap water cluster and the 
        # fixed atoms in the free leg
        if disable_grid.val:
            free_swap_fixed = InterGroupCLJFF("free:swap-fixed")
            free_swap_fixed = setCLJProperties(free_swap_fixed)
            free_swap_fixed = setFakeGridProperties(free_swap_fixed)

            free_swap_fixed.add(swap_water_mols, MGIdx(0))
            free_swap_fixed.add(fixed_free_group, MGIdx(1))
        else:
            free_swap_fixed = GridFF2("free:swap-fixed")
            free_swap_fixed = setCLJProperties(free_swap_fixed)
            # The swap water cluster is more mobile so needs a bigger grid buffer
            free_swap_fixed = setGridProperties(free_swap_fixed, 6*angstrom)

            free_swap_fixed.add(swap_water_mols, MGIdx(0))
            free_swap_fixed.addFixedAtoms(fixed_free_group)
    else:
        free_swap_mobile = InterGroupFF("free:swap-mobile")
        free_swap_mobile.setCLJFunction( getSoftInterCLJFunction() )
        free_swap_mobile.setCLJFunction( "f", getSoftInterCLJFunction() )
        free_swap_mobile.setCLJFunction( "b", getSoftInterCLJFunction() )
        free_swap_mobile.setCLJFunction( "next", getSoftInterCLJFunction() )
        free_swap_mobile.setCLJFunction( "prev", getSoftInterCLJFunction() )

        free_swap_mobile.add(swap_water_mols, MGIdx(0))
        free_swap_mobile.add(mobile_free_mols, MGIdx(1))

        free_swap_fixed = InterGroupFF("free:swap-fixed")
        free_swap_fixed.setCLJFunction( getInterCLJFunction() )
        free_swap_fixed = setNewGridProperties(free_swap_fixed)

        free_swap_fixed.add(swap_water_mols, MGIdx(0))
        free_swap_fixed.addFixedAtoms(fixed_free_group.molecules())

    ###
    ### FORCEFIELDS LOCAL ONLY TO THE BOUND LEG
    ###
    bound_forcefields = []

    if use_oldff.val:
        # forcefield holding the energy between the bound leg mobile atoms and
        # the bound leg fixed atoms
        if disable_grid.val:
            bound_mobile_fixed = InterGroupCLJFF("bound:mobile-fixed")
            bound_mobile_fixed = setCLJProperties(bound_mobile_fixed)
            bound_mobile_fixed = setFakeGridProperties(bound_mobile_fixed)
            bound_mobile_fixed.add(mobile_buffered_bound_mols, MGIdx(0))
            bound_mobile_fixed.add(fixed_bound_group, MGIdx(1))
            bound_forcefields.append(bound_mobile_fixed)
        else:
            bound_mobile_fixed = GridFF2("bound:mobile-fixed")
            bound_mobile_fixed = setCLJProperties(bound_mobile_fixed)
            bound_mobile_fixed = setGridProperties(bound_mobile_fixed)

            # we use mobile_buffered_bound_group as this group misses out atoms that are bonded
            # to fixed atoms (thus preventing large energies caused by incorrect non-bonded calculations)
            bound_mobile_fixed.add(mobile_buffered_bound_mols, MGIdx(0))
            bound_mobile_fixed.addFixedAtoms(fixed_bound_group)

            bound_forcefields.append(bound_mobile_fixed)

        # forcefield holding the intermolecular energy between all bound molecules
        bound_mobile_mobile = InterCLJFF("bound:mobile-mobile")
        bound_mobile_mobile = setCLJProperties(bound_mobile_mobile)

        bound_mobile_mobile.add(mobile_bound_mols)

        bound_forcefields.append(bound_mobile_mobile)
    else:
        # forcefield holding the energy between
        # the bound molecules and bound fixed atoms
        bound_mobile_fixed = InterGroupFF("bound:mobile-fixed")
        bound_mobile_fixed.setCLJFunction( getInterCLJFunction() )
        bound_mobile_fixed = setNewGridProperties(bound_mobile_fixed)

        # we use mobile_buffered_bound_group as this group misses out atoms that are bonded
        # to fixed atoms (thus preventing large energies caused by incorrect non-bonded calculations)
        bound_mobile_fixed.add(mobile_buffered_bound_mols, MGIdx(0))
        bound_mobile_fixed.addFixedAtoms(fixed_bound_group.molecules())
        bound_forcefields.append(bound_mobile_fixed)

        # forcefield holding the energy between all bound mobile molecules
        bound_mobile_mobile = InterFF("bound:mobile-mobile")
        bound_mobile_mobile.setCLJFunction( getInterCLJFunction() )
        bound_mobile_mobile.add(mobile_bound_mols)
        bound_forcefields.append(bound_mobile_mobile)

    # intramolecular energy of the protein
    if bound_protein_intra_mols.nMolecules() > 0:
        if use_oldff.val:
            protein_intraclj = IntraCLJFF("bound:protein_intraclj")
            protein_intraclj = setCLJProperties(protein_intraclj)

            protein_intraff = InternalFF("bound:protein_intra")
            protein_intraff.setUse14Calculation(False)

            for molnum in bound_protein_intra_mols.molNums():
                protein_mol = Molecule.join(bound_protein_intra_mols[molnum])
                protein_intraclj.add(protein_mol)
                protein_intraff.add(protein_mol)

            bound_forcefields.append(protein_intraclj)
            bound_forcefields.append(protein_intraff)
        else:
            protein_intraclj = IntraFF("bound:protein_intraclj")
            protein_intraclj.setCLJFunction( getIntraCLJFunction() )
            protein_intraff = InternalFF("bound:protein_intra")
            protein_intraff.setUse14Calculation(True)

            for molnum in bound_protein_intra_mols.molNums():
                protein_mol = Molecule.join(bound_protein_intra_mols[molnum])
                protein_intraclj.add(protein_mol)
                protein_intraff.add(protein_mol)

            bound_forcefields.append(protein_intraclj)
            bound_forcefields.append(protein_intraff)

    # intramolecular energy of any other solutes
    if bound_solute_intra_mols.nMolecules() > 0:
        if use_oldff.val:
            solute_intraclj = IntraCLJFF("bound:solute_intraclj")
            solute_intraclj = setCLJProperties(solute_intraclj)

            solute_intraff = InternalFF("bound:solute_intra")
            solute_intraff.setUse14Calculation(False)

            for molnum in bound_solute_intra_mols.molNums():
                solute_mol = Molecule.join(bound_solute_intra_mols[molnum])
                solute_intraclj.add(solute_mol)
                solute_intraff.add(solute_mol)

            bound_forcefields.append(solute_intraclj)
            bound_forcefields.append(solute_intraff)
        else:
            solute_intraclj = IntraFF("bound:solute_intraclj")
            solute_intraclj.setCLJFunction( getIntraCLJFunction() )
            solute_intraff = InternalFF("bound:solute_intra")
            solute_intraff.setUse14Calculation(True)

            for molnum in bound_solute_intra_mols.molNums():
                solute_mol = Molecule.join(bound_solute_intra_mols[molnum])
                solute_intraclj.add(solute_mol)
                solute_intraff.add(solute_mol)

            bound_forcefields.append(solute_intraclj)
            bound_forcefields.append(solute_intraff)

    ###
    ### FORCEFIELDS LOCAL ONLY TO THE FREE LEG
    ###
    free_forcefields = []

    if use_oldff.val:
        # forcefield holding the energy between the mobile free molecules and the
        # fixed free molecules
        if disable_grid:
            free_mobile_fixed = InterGroupCLJFF("free:mobile-fixed")
            free_mobile_fixed = setCLJProperties(free_mobile_fixed)
            free_mobile_fixed = setFakeGridProperties(free_mobile_fixed)

            free_mobile_fixed.add(mobile_free_mols, MGIdx(0))
            free_mobile_fixed.add(fixed_free_group, MGIdx(1))
            free_forcefields.append(free_mobile_fixed)
        else:
            free_mobile_fixed = GridFF2("free:mobile-fixed")
            free_mobile_fixed = setCLJProperties(free_mobile_fixed)
            free_mobile_fixed = setGridProperties(free_mobile_fixed)

            free_mobile_fixed.add(mobile_free_mols, MGIdx(0))
            free_mobile_fixed.addFixedAtoms(fixed_free_group)

            free_forcefields.append(free_mobile_fixed)

        # forcefield holding the intermolecular energy between the mobile free molecules
        free_mobile_mobile = InterCLJFF("free:mobile-mobile")
        free_mobile_mobile = setCLJProperties(free_mobile_mobile)

        free_mobile_mobile.add(mobile_free_mols)

        free_forcefields.append(free_mobile_mobile)

    else:
        # forcefield holding the energy between the mobile free molecules, and their
        # interaction with the fixed free molecules
        free_mobile = InterFF("free:mobile")
        free_mobile.setCLJFunction( getInterCLJFunction() )
        free_mobile = setNewGridProperties( free_mobile )

        free_mobile.add(mobile_free_mols)
        free_mobile.addFixedAtoms(fixed_free_group.molecules())

        free_forcefields.append(free_mobile)


    ###
    ### NOW ADD THE FORCEFIELDS TO THE SYSTEM
    ###
    ###
    ### SETTING THE FORCEFIELD EXPRESSIONS
    ###

    ligand_int_nrg_sym = Symbol("E_{ligand:internal}")

    ligand_int_nrg_f_sym = Symbol("E_{ligand:internal_{f}}")
    ligand_int_nrg_b_sym = Symbol("E_{ligand:internal_{b}}")
    ligand_int_nrg_next_sym = Symbol("E_{ligand:internal_{next}}")
    ligand_int_nrg_prev_sym = Symbol("E_{ligand:internal_{prev}}")

    ligand_bound_coul_nrg_sym = Symbol("E_{ligand:bound_coul}")
    ligand_bound_lj_nrg_sym = Symbol("E_{ligand:bound_lj}")
    ligand_bound_coul_nrg_f_sym = Symbol("E_{ligand:bound_coul{f}}")
    ligand_bound_lj_nrg_f_sym = Symbol("E_{ligand:bound_lj{f}}")
    ligand_bound_coul_nrg_b_sym = Symbol("E_{ligand:bound_coul{b}}")
    ligand_bound_lj_nrg_b_sym = Symbol("E_{ligand:bound_lj{b}}")
    ligand_bound_coul_nrg_next_sym = Symbol("E_{ligand:bound_coul{next}}")
    ligand_bound_lj_nrg_next_sym = Symbol("E_{ligand:bound_lj{next}}")
    ligand_bound_coul_nrg_prev_sym = Symbol("E_{ligand:bound_coul{prev}}")
    ligand_bound_lj_nrg_prev_sym = Symbol("E_{ligand:bound_lj{prev}}")

    ligand_free_coul_nrg_sym = Symbol("E_{ligand:free_coul}")
    ligand_free_lj_nrg_sym = Symbol("E_{ligand:free_lj}")
    ligand_free_coul_nrg_f_sym = Symbol("E_{ligand:free_coul{f}}")
    ligand_free_lj_nrg_f_sym = Symbol("E_{ligand:free_lj{f}}")
    ligand_free_coul_nrg_b_sym = Symbol("E_{ligand:free_coul{b}}")
    ligand_free_lj_nrg_b_sym = Symbol("E_{ligand:free_lj{b}}")
    ligand_free_coul_nrg_next_sym = Symbol("E_{ligand:free_coul{next}}")
    ligand_free_lj_nrg_next_sym = Symbol("E_{ligand:free_lj{next}}")
    ligand_free_coul_nrg_prev_sym = Symbol("E_{ligand:free_coul{prev}}")
    ligand_free_lj_nrg_prev_sym = Symbol("E_{ligand:free_lj{prev}}")

    ligand_int_nrg = ligand_intraclj.components().total() + \
                     ligand_intraff.components().total()

    ligand_int_nrg_f = ligand_intraclj.components().total() + \
                       ligand_intraff.components().total()

    ligand_int_nrg_b = ligand_intraclj.components().total() + \
                       ligand_intraff.components().total()

    ligand_int_nrg_next = ligand_intraclj.components().total() + \
                          ligand_intraff.components().total()

    ligand_int_nrg_prev = ligand_intraclj.components().total() + \
                          ligand_intraff.components().total()

    bound_ligand_fixed_coul_nrg = bound_ligand_fixed.components().coulomb()
    bound_ligand_fixed_lj_nrg = bound_ligand_fixed.components().lj()

    free_ligand_fixed_coul_nrg = free_ligand_fixed.components().coulomb()
    free_ligand_fixed_lj_nrg = free_ligand_fixed.components().lj()

    bound_swap_fixed_coul_nrg = bound_swap_fixed.components().coulomb()
    bound_swap_fixed_lj_nrg = bound_swap_fixed.components().lj()

    free_swap_fixed_coul_nrg = free_swap_fixed.components().coulomb()
    free_swap_fixed_lj_nrg = free_swap_fixed.components().lj()

    if use_oldff.val:
        ligand_bound_coul_nrg = bound_ligand_mobile.components().coulomb(0) + \
                                bound_ligand_fixed_coul_nrg

        ligand_bound_lj_nrg = bound_ligand_mobile.components().lj(0) + \
                              bound_ligand_fixed_lj_nrg

        ligand_bound_coul_nrg_f = bound_ligand_mobile.components().coulomb(1) + \
                                  bound_ligand_fixed_coul_nrg

        ligand_bound_lj_nrg_f = bound_ligand_mobile.components().lj(1) + \
                                bound_ligand_fixed_lj_nrg

        ligand_bound_coul_nrg_b = bound_ligand_mobile.components().coulomb(2) + \
                                  bound_ligand_fixed_coul_nrg

        ligand_bound_lj_nrg_b = bound_ligand_mobile.components().lj(2) + \
                                bound_ligand_fixed_lj_nrg

        ligand_bound_coul_nrg_next = bound_ligand_mobile.components().coulomb(3) + \
                                     bound_ligand_fixed_coul_nrg

        ligand_bound_lj_nrg_next = bound_ligand_mobile.components().lj(3) + \
                                   bound_ligand_fixed_lj_nrg

        ligand_bound_coul_nrg_prev = bound_ligand_mobile.components().coulomb(4) + \
                                     bound_ligand_fixed_coul_nrg

        ligand_bound_lj_nrg_prev = bound_ligand_mobile.components().lj(4) + \
                                   bound_ligand_fixed_lj_nrg

        ligand_free_coul_nrg = free_ligand_mobile.components().coulomb(0) + \
                               free_ligand_fixed_coul_nrg
 
        ligand_free_lj_nrg = free_ligand_mobile.components().lj(0) + \
                             free_ligand_fixed_lj_nrg

        ligand_free_coul_nrg_f = free_ligand_mobile.components().coulomb(1) + \
                                 free_ligand_fixed_coul_nrg

        ligand_free_lj_nrg_f = free_ligand_mobile.components().lj(1) + \
                               free_ligand_fixed_lj_nrg

        ligand_free_coul_nrg_b = free_ligand_mobile.components().coulomb(2) + \
                                 free_ligand_fixed_coul_nrg
 
        ligand_free_lj_nrg_b = free_ligand_mobile.components().lj(2) + \
                               free_ligand_fixed_lj_nrg

        ligand_free_coul_nrg_next = free_ligand_mobile.components().coulomb(3) + \
                                    free_ligand_fixed_coul_nrg

        ligand_free_lj_nrg_next = free_ligand_mobile.components().lj(3) + \
                                  free_ligand_fixed_lj_nrg

        ligand_free_coul_nrg_prev = free_ligand_mobile.components().coulomb(4) + \
                                    free_ligand_fixed_coul_nrg

        ligand_free_lj_nrg_prev = free_ligand_mobile.components().lj(4) + \
                                  free_ligand_fixed_lj_nrg
    else:
        ligand_bound_coul_nrg = bound_ligand_mobile.components().coulomb() + \
                                bound_ligand_fixed_coul_nrg

        ligand_bound_lj_nrg = bound_ligand_mobile.components().lj() + \
                              bound_ligand_fixed_lj_nrg

        ligand_bound_coul_nrg_f = bound_ligand_mobile.components().coulomb("f") + \
                                  bound_ligand_fixed_coul_nrg

        ligand_bound_lj_nrg_f = bound_ligand_mobile.components().lj("f") + \
                                bound_ligand_fixed_lj_nrg

        ligand_bound_coul_nrg_b = bound_ligand_mobile.components().coulomb("b") + \
                                  bound_ligand_fixed_coul_nrg

        ligand_bound_lj_nrg_b = bound_ligand_mobile.components().lj("b") + \
                                bound_ligand_fixed_lj_nrg

        ligand_bound_coul_nrg_next = bound_ligand_mobile.components().coulomb("next") + \
                                     bound_ligand_fixed_coul_nrg

        ligand_bound_lj_nrg_next = bound_ligand_mobile.components().lj("next") + \
                                   bound_ligand_fixed_lj_nrg

        ligand_bound_coul_nrg_prev = bound_ligand_mobile.components().coulomb("prev") + \
                                     bound_ligand_fixed_coul_nrg

        ligand_bound_lj_nrg_prev = bound_ligand_mobile.components().lj("prev") + \
                                   bound_ligand_fixed_lj_nrg

        ligand_free_coul_nrg = free_ligand_mobile.components().coulomb() + \
                               free_ligand_fixed_coul_nrg
 
        ligand_free_lj_nrg = free_ligand_mobile.components().lj() + \
                             free_ligand_fixed_lj_nrg

        ligand_free_coul_nrg_f = free_ligand_mobile.components().coulomb("f") + \
                                 free_ligand_fixed_coul_nrg

        ligand_free_lj_nrg_f = free_ligand_mobile.components().lj("f") + \
                               free_ligand_fixed_lj_nrg

        ligand_free_coul_nrg_b = free_ligand_mobile.components().coulomb("b") + \
                                 free_ligand_fixed_coul_nrg
 
        ligand_free_lj_nrg_b = free_ligand_mobile.components().lj("b") + \
                               free_ligand_fixed_lj_nrg

        ligand_free_coul_nrg_next = free_ligand_mobile.components().coulomb("next") + \
                                    free_ligand_fixed_coul_nrg

        ligand_free_lj_nrg_next = free_ligand_mobile.components().lj("next") + \
                                  free_ligand_fixed_lj_nrg

        ligand_free_coul_nrg_prev = free_ligand_mobile.components().coulomb("prev") + \
                                    free_ligand_fixed_coul_nrg

        ligand_free_lj_nrg_prev = free_ligand_mobile.components().lj("prev") + \
                                  free_ligand_fixed_lj_nrg

    lam = Symbol("lambda")
    lam_f = Symbol("lambda_{f}")
    lam_b = Symbol("lambda_{b}")
    lam_next = Symbol("lambda_{next}")
    lam_prev = Symbol("lambda_{prev}")

    lam_coul_on = Symbol("lambda_coul_on")
    lam_coul_on_f = Symbol("lambda_coul_on_f")
    lam_coul_on_b = Symbol("lambda_coul_on_b")
    lam_coul_on_next = Symbol("lambda_coul_on_next")
    lam_coul_on_prev = Symbol("lambda_coul_on_prev")

    lam_coul_off = Symbol("lambda_coul_off")
    lam_coul_off_f = Symbol("lambda_coul_off_f")
    lam_coul_off_b = Symbol("lambda_coul_off_b")
    lam_coul_off_next = Symbol("lambda_coul_off_next")
    lam_coul_off_prev = Symbol("lambda_coul_off_prev")

    lam_lj_on = Symbol("lambda_lj_on")
    lam_lj_on_f = Symbol("lambda_lj_on_f")
    lam_lj_on_b = Symbol("lambda_lj_on_b")
    lam_lj_on_next = Symbol("lambda_lj_on_next")
    lam_lj_on_prev = Symbol("lambda_lj_on_prev")

    lam_lj_off = Symbol("lambda_lj_off")
    lam_lj_off_f = Symbol("lambda_lj_off_f")
    lam_lj_off_b = Symbol("lambda_lj_off_b")
    lam_lj_off_next = Symbol("lambda_lj_off_next")
    lam_lj_off_prev = Symbol("lambda_lj_off_prev")

    lam_coul_swap = Symbol("lambda_coul_swap")
    lam_coul_swap_f = Symbol("lambda_coul_swap_f")
    lam_coul_swap_b = Symbol("lambda_coul_swap_b")
    lam_coul_swap_next = Symbol("lambda_coul_swap_next")
    lam_coul_swap_prev = Symbol("lambda_coul_swap_prev")

    lam_lj_swap = Symbol("lambda_lj_swap")
    lam_lj_swap_f = Symbol("lambda_lj_swap_f")
    lam_lj_swap_b = Symbol("lambda_lj_swap_b")
    lam_lj_swap_next = Symbol("lambda_lj_swap_next")
    lam_lj_swap_prev = Symbol("lambda_lj_swap_prev")

    S_sym = Symbol("S")
    S_scl = S_sym - 4*(S_sym-1)*(lam-0.5)**2
    S_scl_f = S_sym - 4*(S_sym-1)*(lam_f-0.5)**2
    S_scl_b = S_sym - 4*(S_sym-1)*(lam_b-0.5)**2
    S_scl_next = S_sym - 4*(S_sym-1)*(lam_next-0.5)**2
    S_scl_prev = S_sym - 4*(S_sym-1)*(lam_prev-0.5)**2

    swap_int_nrg_sym = Symbol("E_{swap:internal}")
    swap_int_nrg_f_sym = Symbol("E_{swap:internal_{f}}")
    swap_int_nrg_b_sym = Symbol("E_{swap:internal_{b}}")
    swap_int_nrg_next_sym = Symbol("E_{swap:internal_{next}}")
    swap_int_nrg_prev_sym = Symbol("E_{swap:internal_{prev}}")

    if use_oldff.val:
        swap_int_nrg = (lam_coul_swap * S_scl * swap_interclj.components().coulomb(0)) + \
                       (lam_lj_swap * swap_interclj.components().lj(0))

        swap_int_nrg_f = (lam_coul_swap_f * S_scl_f * swap_interclj.components().coulomb(1)) + \
                         (lam_lj_swap_f * swap_interclj.components().lj(1))

        swap_int_nrg_b = (lam_coul_swap_b * S_scl_b * swap_interclj.components().coulomb(2)) + \
                         (lam_lj_swap_b * swap_interclj.components().lj(2))

        swap_int_nrg_next = (lam_coul_swap_next * S_scl_next * swap_interclj.components().coulomb(3)) + \
                            (lam_lj_swap_next * swap_interclj.components().lj(3))

        swap_int_nrg_prev = (lam_coul_swap_prev * S_scl_prev * swap_interclj.components().coulomb(4)) + \
                            (lam_lj_swap_prev * swap_interclj.components().lj(4))
    else:
        swap_int_nrg = (lam_coul_swap * S_scl * swap_interclj.components().coulomb()) + \
                       (lam_lj_swap * swap_interclj.components().lj())

        swap_int_nrg_f = (lam_coul_swap_f * S_scl_f * swap_interclj.components().coulomb("f")) + \
                         (lam_lj_swap_f * swap_interclj.components().lj("f"))

        swap_int_nrg_b = (lam_coul_swap_b * S_scl_b * swap_interclj.components().coulomb("b")) + \
                         (lam_lj_swap_b * swap_interclj.components().lj("b"))

        swap_int_nrg_next = (lam_coul_swap_next * S_scl_next * swap_interclj.components().coulomb("next")) + \
                            (lam_lj_swap_next * swap_interclj.components().lj("next"))

        swap_int_nrg_prev = (lam_coul_swap_prev * S_scl_prev * swap_interclj.components().coulomb("prev")) + \
                            (lam_lj_swap_prev * swap_interclj.components().lj("prev"))

    swap_bound_coul_nrg_sym = Symbol("E_{swap:bound_coul}")
    swap_bound_lj_nrg_sym = Symbol("E_{swap:bound_lj}")
    swap_bound_coul_nrg_f_sym = Symbol("E_{swap:bound_coul{f}}")
    swap_bound_lj_nrg_f_sym = Symbol("E_{swap:bound_lj{f}}")
    swap_bound_coul_nrg_b_sym = Symbol("E_{swap:bound_coul{b}}")
    swap_bound_lj_nrg_b_sym = Symbol("E_{swap:bound_lj{b}}")
    swap_bound_coul_nrg_next_sym = Symbol("E_{swap:bound_coul{next}}")
    swap_bound_lj_nrg_next_sym = Symbol("E_{swap:bound_lj{next}}")
    swap_bound_coul_nrg_prev_sym = Symbol("E_{swap:bound_coul{prev}}")
    swap_bound_lj_nrg_prev_sym = Symbol("E_{swap:bound_lj{prev}}")

    swap_free_coul_nrg_sym = Symbol("E_{swap:free_coul}")
    swap_free_lj_nrg_sym = Symbol("E_{swap:free_lj}")
    swap_free_coul_nrg_f_sym = Symbol("E_{swap:free_coul{f}}")
    swap_free_lj_nrg_f_sym = Symbol("E_{swap:free_lj{f}}")
    swap_free_coul_nrg_b_sym = Symbol("E_{swap:free_coul{b}}")
    swap_free_lj_nrg_b_sym = Symbol("E_{swap:free_lj{b}}")
    swap_free_coul_nrg_next_sym = Symbol("E_{swap:free_coul{next}}")
    swap_free_lj_nrg_next_sym = Symbol("E_{swap:free_lj{next}}")
    swap_free_coul_nrg_prev_sym = Symbol("E_{swap:free_coul{prev}}")
    swap_free_lj_nrg_prev_sym = Symbol("E_{swap:free_lj{prev}}")

    if use_oldff.val:
        swap_bound_coul_nrg = bound_swap_mobile.components().coulomb(0) + \
                              bound_swap_fixed_coul_nrg

        swap_bound_lj_nrg = bound_swap_mobile.components().lj(0) + \
                            bound_swap_fixed_lj_nrg

        swap_bound_coul_nrg_f = bound_swap_mobile.components().coulomb(1) + \
                                bound_swap_fixed_coul_nrg

        swap_bound_lj_nrg_f = bound_swap_mobile.components().lj(1) + \
                              bound_swap_fixed_lj_nrg

        swap_bound_coul_nrg_b = bound_swap_mobile.components().coulomb(2) + \
                                bound_swap_fixed_coul_nrg

        swap_bound_lj_nrg_b = bound_swap_mobile.components().lj(2) + \
                              bound_swap_fixed_lj_nrg

        swap_bound_coul_nrg_next = bound_swap_mobile.components().coulomb(3) + \
                                   bound_swap_fixed_coul_nrg

        swap_bound_lj_nrg_next = bound_swap_mobile.components().lj(3) + \
                                 bound_swap_fixed_lj_nrg

        swap_bound_coul_nrg_prev = bound_swap_mobile.components().coulomb(4) + \
                                   bound_swap_fixed_coul_nrg

        swap_bound_lj_nrg_prev = bound_swap_mobile.components().lj(4) + \
                                 bound_swap_fixed_lj_nrg

        swap_free_coul_nrg = free_swap_mobile.components().coulomb(0) + \
                             free_swap_fixed_coul_nrg

        swap_free_lj_nrg = free_swap_mobile.components().lj(0) + \
                           free_swap_fixed_lj_nrg

        swap_free_coul_nrg_f = free_swap_mobile.components().coulomb(1) + \
                               free_swap_fixed_coul_nrg

        swap_free_lj_nrg_f = free_swap_mobile.components().lj(1) + \
                             free_swap_fixed_lj_nrg

        swap_free_coul_nrg_b = free_swap_mobile.components().coulomb(2) + \
                               free_swap_fixed_coul_nrg

        swap_free_lj_nrg_b = free_swap_mobile.components().lj(2) + \
                             free_swap_fixed_lj_nrg

        swap_free_coul_nrg_next = free_swap_mobile.components().coulomb(3) + \
                                  free_swap_fixed_coul_nrg
 
        swap_free_lj_nrg_next = free_swap_mobile.components().lj(3) + \
                                free_swap_fixed_lj_nrg

        swap_free_coul_nrg_prev = free_swap_mobile.components().coulomb(4) + \
                                  free_swap_fixed_coul_nrg

        swap_free_lj_nrg_prev = free_swap_mobile.components().lj(4) + \
                                free_swap_fixed_lj_nrg
    else:
        swap_bound_coul_nrg = bound_swap_mobile.components().coulomb() + \
                              bound_swap_fixed_coul_nrg

        swap_bound_lj_nrg = bound_swap_mobile.components().lj() + \
                            bound_swap_fixed_lj_nrg

        swap_bound_coul_nrg_f = bound_swap_mobile.components().coulomb("f") + \
                                bound_swap_fixed_coul_nrg

        swap_bound_lj_nrg_f = bound_swap_mobile.components().lj("f") + \
                              bound_swap_fixed_lj_nrg

        swap_bound_coul_nrg_b = bound_swap_mobile.components().coulomb("b") + \
                                bound_swap_fixed_coul_nrg

        swap_bound_lj_nrg_b = bound_swap_mobile.components().lj("b") + \
                              bound_swap_fixed_lj_nrg

        swap_bound_coul_nrg_next = bound_swap_mobile.components().coulomb("next") + \
                                   bound_swap_fixed_coul_nrg

        swap_bound_lj_nrg_next = bound_swap_mobile.components().lj("next") + \
                                 bound_swap_fixed_lj_nrg

        swap_bound_coul_nrg_prev = bound_swap_mobile.components().coulomb("prev") + \
                                   bound_swap_fixed_coul_nrg

        swap_bound_lj_nrg_prev = bound_swap_mobile.components().lj("prev") + \
                                 bound_swap_fixed_lj_nrg

        swap_free_coul_nrg = free_swap_mobile.components().coulomb() + \
                             free_swap_fixed_coul_nrg

        swap_free_lj_nrg = free_swap_mobile.components().lj() + \
                           free_swap_fixed_lj_nrg

        swap_free_coul_nrg_f = free_swap_mobile.components().coulomb("f") + \
                               free_swap_fixed_coul_nrg

        swap_free_lj_nrg_f = free_swap_mobile.components().lj("f") + \
                             free_swap_fixed_lj_nrg

        swap_free_coul_nrg_b = free_swap_mobile.components().coulomb("b") + \
                               free_swap_fixed_coul_nrg

        swap_free_lj_nrg_b = free_swap_mobile.components().lj("b") + \
                             free_swap_fixed_lj_nrg

        swap_free_coul_nrg_next = free_swap_mobile.components().coulomb("next") + \
                                  free_swap_fixed_coul_nrg
 
        swap_free_lj_nrg_next = free_swap_mobile.components().lj("next") + \
                                free_swap_fixed_lj_nrg

        swap_free_coul_nrg_prev = free_swap_mobile.components().coulomb("prev") + \
                                  free_swap_fixed_coul_nrg

        swap_free_lj_nrg_prev = free_swap_mobile.components().lj("prev") + \
                                free_swap_fixed_lj_nrg

    system.add(ligand_intraclj)
    system.add(ligand_intraff)
    system.add(swap_interclj)
    system.add(bound_ligand_mobile)
    system.add(bound_swap_mobile)
    system.add(free_ligand_mobile)
    system.add(free_swap_mobile)
    system.add(bound_ligand_fixed)
    system.add(bound_swap_fixed)
    system.add(free_ligand_fixed)
    system.add(free_swap_fixed)

    system.setConstant(lam, 0.0)
    system.setConstant(lam_f, 0.0)
    system.setConstant(lam_b, 0.0)
    system.setConstant(lam_next, 0.0)
    system.setConstant(lam_prev, 0.0)

    system.setConstant(lam_coul_on, 1.0)
    system.setConstant(lam_coul_on_f, 1.0)
    system.setConstant(lam_coul_on_b, 1.0)
    system.setConstant(lam_coul_on_next, 1.0)
    system.setConstant(lam_coul_on_prev, 1.0)

    system.setConstant(lam_coul_off, 0.0)
    system.setConstant(lam_coul_off_f, 0.0)
    system.setConstant(lam_coul_off_b, 0.0)
    system.setConstant(lam_coul_off_next, 0.0)
    system.setConstant(lam_coul_off_prev, 0.0)

    system.setConstant(lam_lj_on, 1.0)
    system.setConstant(lam_lj_on_f, 1.0)
    system.setConstant(lam_lj_on_b, 1.0)
    system.setConstant(lam_lj_on_next, 1.0)
    system.setConstant(lam_lj_on_prev, 1.0)

    system.setConstant(lam_lj_off, 0.0)
    system.setConstant(lam_lj_off_f, 0.0)
    system.setConstant(lam_lj_off_b, 0.0)
    system.setConstant(lam_lj_off_next, 0.0)
    system.setConstant(lam_lj_off_prev, 0.0)

    system.setConstant(lam_coul_swap, 1.0)
    system.setConstant(lam_coul_swap_f, 1.0)
    system.setConstant(lam_coul_swap_b, 1.0)
    system.setConstant(lam_coul_swap_next, 1.0)
    system.setConstant(lam_coul_swap_prev, 1.0)

    system.setConstant(lam_lj_swap, 1.0)
    system.setConstant(lam_lj_swap_f, 1.0)
    system.setConstant(lam_lj_swap_b, 1.0)
    system.setConstant(lam_lj_swap_next, 1.0)
    system.setConstant(lam_lj_swap_prev, 1.0)

    if uncharge_ligand.val:
        system.setComponent(S_sym, 1.0)
    else:
        system.setComponent(S_sym, soften_water.val)

    system.setComponent(ligand_int_nrg_sym, ligand_int_nrg)
    system.setComponent(ligand_int_nrg_f_sym, ligand_int_nrg_f)
    system.setComponent(ligand_int_nrg_b_sym, ligand_int_nrg_b)
    system.setComponent(ligand_int_nrg_next_sym, ligand_int_nrg_next)
    system.setComponent(ligand_int_nrg_prev_sym, ligand_int_nrg_prev)

    system.setComponent(ligand_bound_coul_nrg_sym, ligand_bound_coul_nrg)
    system.setComponent(ligand_bound_coul_nrg_f_sym, ligand_bound_coul_nrg_f)
    system.setComponent(ligand_bound_coul_nrg_b_sym, ligand_bound_coul_nrg_b)
    system.setComponent(ligand_bound_coul_nrg_next_sym, ligand_bound_coul_nrg_next)
    system.setComponent(ligand_bound_coul_nrg_prev_sym, ligand_bound_coul_nrg_prev)

    system.setComponent(ligand_bound_lj_nrg_sym, ligand_bound_lj_nrg)
    system.setComponent(ligand_bound_lj_nrg_f_sym, ligand_bound_lj_nrg_f)
    system.setComponent(ligand_bound_lj_nrg_b_sym, ligand_bound_lj_nrg_b)
    system.setComponent(ligand_bound_lj_nrg_next_sym, ligand_bound_lj_nrg_next)
    system.setComponent(ligand_bound_lj_nrg_prev_sym, ligand_bound_lj_nrg_prev)

    system.setComponent(ligand_free_coul_nrg_sym, ligand_free_coul_nrg)
    system.setComponent(ligand_free_coul_nrg_f_sym, ligand_free_coul_nrg_f)
    system.setComponent(ligand_free_coul_nrg_b_sym, ligand_free_coul_nrg_b)
    system.setComponent(ligand_free_coul_nrg_next_sym, ligand_free_coul_nrg_next)
    system.setComponent(ligand_free_coul_nrg_prev_sym, ligand_free_coul_nrg_prev)

    system.setComponent(ligand_free_lj_nrg_sym, ligand_free_lj_nrg)
    system.setComponent(ligand_free_lj_nrg_f_sym, ligand_free_lj_nrg_f)
    system.setComponent(ligand_free_lj_nrg_b_sym, ligand_free_lj_nrg_b)
    system.setComponent(ligand_free_lj_nrg_next_sym, ligand_free_lj_nrg_next)
    system.setComponent(ligand_free_lj_nrg_prev_sym, ligand_free_lj_nrg_prev)

    system.setComponent(swap_int_nrg_sym, swap_int_nrg)
    system.setComponent(swap_int_nrg_f_sym, swap_int_nrg_f)
    system.setComponent(swap_int_nrg_b_sym, swap_int_nrg_b)
    system.setComponent(swap_int_nrg_next_sym, swap_int_nrg_next)
    system.setComponent(swap_int_nrg_prev_sym, swap_int_nrg_prev)

    system.setComponent(swap_bound_coul_nrg_sym, swap_bound_coul_nrg)
    system.setComponent(swap_bound_coul_nrg_f_sym, swap_bound_coul_nrg_f)
    system.setComponent(swap_bound_coul_nrg_b_sym, swap_bound_coul_nrg_b)
    system.setComponent(swap_bound_coul_nrg_next_sym, swap_bound_coul_nrg_next)
    system.setComponent(swap_bound_coul_nrg_prev_sym, swap_bound_coul_nrg_prev)

    system.setComponent(swap_bound_lj_nrg_sym, swap_bound_lj_nrg)
    system.setComponent(swap_bound_lj_nrg_f_sym, swap_bound_lj_nrg_f)
    system.setComponent(swap_bound_lj_nrg_b_sym, swap_bound_lj_nrg_b)
    system.setComponent(swap_bound_lj_nrg_next_sym, swap_bound_lj_nrg_next)
    system.setComponent(swap_bound_lj_nrg_prev_sym, swap_bound_lj_nrg_prev)

    system.setComponent(swap_free_coul_nrg_sym, swap_free_coul_nrg)
    system.setComponent(swap_free_coul_nrg_f_sym, swap_free_coul_nrg_f)
    system.setComponent(swap_free_coul_nrg_b_sym, swap_free_coul_nrg_b)
    system.setComponent(swap_free_coul_nrg_next_sym, swap_free_coul_nrg_next)
    system.setComponent(swap_free_coul_nrg_prev_sym, swap_free_coul_nrg_prev)

    system.setComponent(swap_free_lj_nrg_sym, swap_free_lj_nrg)
    system.setComponent(swap_free_lj_nrg_f_sym, swap_free_lj_nrg_f)
    system.setComponent(swap_free_lj_nrg_b_sym, swap_free_lj_nrg_b)
    system.setComponent(swap_free_lj_nrg_next_sym, swap_free_lj_nrg_next)
    system.setComponent(swap_free_lj_nrg_prev_sym, swap_free_lj_nrg_prev)

    bound_bound_nrg_sym = Symbol("E_{bound-bound}")
    bound_bound_nrg = None
    
    for bound_forcefield in bound_forcefields:
        if bound_bound_nrg is None:
            bound_bound_nrg = bound_forcefield.components().total()
        else:
            bound_bound_nrg = bound_bound_nrg + bound_forcefield.components().total()

        system.add(bound_forcefield)

    system.setComponent(bound_bound_nrg_sym, bound_bound_nrg)

    free_free_nrg_sym = Symbol("E_{free-free}")
    free_free_nrg = None

    for free_forcefield in free_forcefields:
        if free_free_nrg is None:
            free_free_nrg = free_forcefield.components().total()
        else:
            free_free_nrg = free_free_nrg + free_forcefield.components().total()

        system.add(free_forcefield)

    system.setComponent(free_free_nrg_sym, free_free_nrg)

    bound_nrg_sym = Symbol("E_{bound}")
    bound_nrg = (lam_coul_on * ligand_bound_coul_nrg_sym) + (lam_coul_off * swap_bound_coul_nrg_sym) + \
                (lam_lj_on * ligand_bound_lj_nrg_sym) + (lam_lj_off * swap_bound_lj_nrg_sym)

    bound_nrg_f_sym = Symbol("E_{bound_{f}}")
    bound_nrg_f = (lam_coul_on_f * ligand_bound_coul_nrg_f_sym) + (lam_coul_off_f * swap_bound_coul_nrg_f_sym) + \
                  (lam_lj_on_f * ligand_bound_lj_nrg_f_sym) + (lam_lj_off_f * swap_bound_lj_nrg_f_sym)

    bound_nrg_b_sym = Symbol("E_{bound_{b}}")
    bound_nrg_b = (lam_coul_on_b * ligand_bound_coul_nrg_b_sym) + (lam_coul_off_b * swap_bound_coul_nrg_b_sym) + \
                  (lam_lj_on_b * ligand_bound_lj_nrg_b_sym) + (lam_lj_off_b * swap_bound_lj_nrg_b_sym)

    bound_nrg_next_sym = Symbol("E_{bound_{next}}")
    bound_nrg_next = (lam_coul_on_next * ligand_bound_coul_nrg_next_sym) + (lam_coul_off_next * swap_bound_coul_nrg_next_sym) + \
                     (lam_lj_on_next * ligand_bound_lj_nrg_next_sym) + (lam_lj_off_next * swap_bound_lj_nrg_next_sym)

    bound_nrg_prev_sym = Symbol("E_{bound_{prev}}")
    bound_nrg_prev = (lam_coul_on_prev * ligand_bound_coul_nrg_prev_sym) + (lam_coul_off_prev * swap_bound_coul_nrg_prev_sym) + \
                     (lam_lj_on_prev * ligand_bound_lj_nrg_prev_sym) + (lam_lj_off_prev * swap_bound_lj_nrg_prev_sym)
    
    free_nrg_sym = Symbol("E_{free}")
    free_nrg = (lam_coul_off * ligand_free_coul_nrg_sym) + (lam_coul_on * swap_free_coul_nrg_sym) + \
               (lam_lj_off * ligand_free_lj_nrg_sym) + (lam_lj_on * swap_free_lj_nrg_sym)

    free_nrg_f_sym = Symbol("E_{free_{f}}")
    free_nrg_f = (lam_coul_off_f * ligand_free_coul_nrg_f_sym) + (lam_coul_on_f * swap_free_coul_nrg_f_sym) + \
                 (lam_lj_off_f * ligand_free_lj_nrg_f_sym) + (lam_lj_on_f * swap_free_lj_nrg_f_sym)

    free_nrg_b_sym = Symbol("E_{free_{b}}")
    free_nrg_b = (lam_coul_off_b * ligand_free_coul_nrg_b_sym) + (lam_coul_on_b * swap_free_coul_nrg_b_sym) + \
                 (lam_lj_off_b * ligand_free_lj_nrg_b_sym) + (lam_lj_on_b * swap_free_lj_nrg_b_sym)

    free_nrg_next_sym = Symbol("E_{free_{next}}")
    free_nrg_next = (lam_coul_off_next * ligand_free_coul_nrg_next_sym) + (lam_coul_on_next * swap_free_coul_nrg_next_sym) + \
                    (lam_lj_off_next * ligand_free_lj_nrg_next_sym) + (lam_lj_on_next * swap_free_lj_nrg_next_sym)

    free_nrg_prev_sym = Symbol("E_{free_{prev}}")
    free_nrg_prev = (lam_coul_off_prev * ligand_free_coul_nrg_prev_sym) + (lam_coul_on_prev * swap_free_coul_nrg_prev_sym) + \
                    (lam_lj_off_prev * ligand_free_lj_nrg_prev_sym) + (lam_lj_on_prev * swap_free_lj_nrg_prev_sym)

    box_nrg_sym = Symbol("E_{box}")
    box_nrg = bound_bound_nrg_sym + free_free_nrg_sym + ligand_int_nrg_sym + swap_int_nrg_sym

    box_nrg_f_sym = Symbol("E_{box_{f}}")
    box_nrg_f = bound_bound_nrg_sym + free_free_nrg_sym + ligand_int_nrg_f_sym + swap_int_nrg_f_sym

    box_nrg_b_sym = Symbol("E_{box_{b}}")
    box_nrg_b = bound_bound_nrg_sym + free_free_nrg_sym + ligand_int_nrg_b_sym + swap_int_nrg_b_sym

    box_nrg_next_sym = Symbol("E_{box_{next}}")
    box_nrg_next = bound_bound_nrg_sym + free_free_nrg_sym + ligand_int_nrg_next_sym + swap_int_nrg_next_sym

    box_nrg_prev_sym = Symbol("E_{box_{prev}}")
    box_nrg_prev = bound_bound_nrg_sym + free_free_nrg_sym + ligand_int_nrg_prev_sym + swap_int_nrg_prev_sym

    total_nrg_sym = system.totalComponent()
    total_nrg = bound_nrg_sym + free_nrg_sym + box_nrg_sym

    total_nrg_f_sym = Symbol("E_{total_{f}}")
    total_nrg_f = bound_nrg_f_sym + free_nrg_f_sym + box_nrg_f_sym

    total_nrg_b_sym = Symbol("E_{total_{b}}")
    total_nrg_b = bound_nrg_b_sym + free_nrg_b_sym + box_nrg_b_sym

    total_nrg_next_sym = Symbol("E_{total_{next}}")
    total_nrg_next = bound_nrg_next_sym + free_nrg_next_sym + box_nrg_next_sym

    total_nrg_prev_sym = Symbol("E_{total_{prev}}")
    total_nrg_prev = bound_nrg_prev_sym + free_nrg_prev_sym + box_nrg_prev_sym

    system.setComponent(bound_nrg_sym, bound_nrg)
    system.setComponent(bound_nrg_f_sym, bound_nrg_f)
    system.setComponent(bound_nrg_b_sym, bound_nrg_b)
    system.setComponent(bound_nrg_next_sym, bound_nrg_next)
    system.setComponent(bound_nrg_prev_sym, bound_nrg_prev)

    system.setComponent(free_nrg_sym, free_nrg)
    system.setComponent(free_nrg_f_sym, free_nrg_f)
    system.setComponent(free_nrg_b_sym, free_nrg_b)
    system.setComponent(free_nrg_next_sym, free_nrg_next)
    system.setComponent(free_nrg_prev_sym, free_nrg_prev)

    system.setComponent(box_nrg_sym, box_nrg)
    system.setComponent(box_nrg_f_sym, box_nrg_f)
    system.setComponent(box_nrg_b_sym, box_nrg_b)
    system.setComponent(box_nrg_next_sym, box_nrg_next)
    system.setComponent(box_nrg_prev_sym, box_nrg_prev)

    system.setComponent(total_nrg_sym, total_nrg)
    system.setComponent(total_nrg_f_sym, total_nrg_f)
    system.setComponent(total_nrg_b_sym, total_nrg_b)
    system.setComponent(total_nrg_next_sym, total_nrg_next)
    system.setComponent(total_nrg_prev_sym, total_nrg_prev)

    system.setComponent( Symbol("delta_nrg^{F}"), (total_nrg_f_sym - total_nrg_sym) )
    system.setComponent( Symbol("delta_nrg^{B}"), (total_nrg_b_sym - total_nrg_sym) )
    system.setComponent( Symbol("delta_nrg^{next}"), (total_nrg_next_sym - total_nrg_sym) )
    system.setComponent( Symbol("delta_nrg^{prev}"), (total_nrg_prev_sym - total_nrg_sym) )

    system.setComponent( Symbol("delta_bound_nrg^{F}"), (bound_nrg_f_sym - bound_nrg_sym) )
    system.setComponent( Symbol("delta_bound_nrg^{B}"), (bound_nrg_b_sym - bound_nrg_sym) )
    system.setComponent( Symbol("delta_free_nrg^{F}"), (free_nrg_f_sym - free_nrg_sym) )
    system.setComponent( Symbol("delta_free_nrg^{B}"), (free_nrg_b_sym - free_nrg_sym) )

    # Now add constraints. These are used to keep the identity of the 
    # swap water, to keep all lambda values between 0 and 1, and to
    # map the alpha values of the softcore forcefields to lambda
    print("\nCreating WSRC system constraints...\n")

    # Add the constraint that lambda_f is lambda + delta_lambda and
    # lambda_b is lambda - delta_lambda (kept to between 0 and 1)
    dlam = delta_lambda.val

    if dlam > 1 or dlam < 0.0000001:
        print("WARNING: Weird value of delta_lambda (%s). Setting it to 0.001" % dlam)
        dlam = 0.001

    # Constrain lam_f and lam_b to lie with delta_lambda of lambda
    dlam_sym = Symbol("delta_lambda")
    system.setConstant( dlam_sym, dlam )
    system.add( ComponentConstraint( lam_f, Min( lam + dlam_sym, 1 ) ) )
    system.add( ComponentConstraint( lam_b, Max( lam - dlam_sym, 0 ) ) )

    # Constrain lam_next and lam_prev to be equal to the next and previous
    # windows lambda values
    lamvals = getLambdaValues()

    print("Using lambda values: %s" % lamvals)

    if lamvals[-1] != 1:
        lamvals.append(1)

    if lamvals[0] != 0:
        lamvals.insert(0,0)

    system.add( WindowedComponent( lam_next, lam, lamvals, 1 ) )
    system.add( WindowedComponent( lam_prev, lam, lamvals, -1 ) )
    system.setConstant( lam, lambda_values.val[0] )

    # work out the maximum and minimum permissable values of lam_lj
    lam_lj_min = lj_buffer.val
    lam_lj_max = 1.0 - lj_buffer.val

    # now constrain lam_coul_on, lam_coul_off, lam_lj_on and lam_lj_off to follow lambda
    if uncharge_ligand.val:
        v = uncharge_ligand_max.val

        vfunc = Max
        if v > 0.75:
            vfunc = Min

        system.add( ComponentConstraint( lam_coul_on, vfunc( 1 + (4*lam*(v-1.0)), (v/0.75)*(1.0-lam) ) ) )
        system.add( ComponentConstraint( lam_coul_off, vfunc( 1 + (4*(1-lam)*(v-1.0)), (v/0.75)*(1.0-(1-lam)) ) ) )
        system.add( ComponentConstraint( lam_coul_on_f, vfunc( 1 + (4*lam_f*(v-1.0)), (v/0.75)*(1.0-lam_f) ) ) )
        system.add( ComponentConstraint( lam_coul_off_f, vfunc( 1 + (4*(1-lam_f)*(v-1.0)), (v/0.75)*(1.0-(1-lam_f)) ) ) )
        system.add( ComponentConstraint( lam_coul_on_b, vfunc( 1 + (4*lam_b*(v-1.0)), (v/0.75)*(1.0-lam_b) ) ) )
        system.add( ComponentConstraint( lam_coul_off_b, vfunc( 1 + (4*(1-lam_b)*(v-1.0)), (v/0.75)*(1.0-(1-lam_b)) ) ) )
        system.add( ComponentConstraint( lam_coul_on_next, vfunc( 1 + (4*lam_next*(v-1.0)), (v/0.75)*(1.0-lam_next) ) ) )
        system.add( ComponentConstraint( lam_coul_off_next, vfunc( 1 + (4*(1-lam_next)*(v-1.0)), (v/0.75)*(1.0-(1-lam_next)) ) ) )
        system.add( ComponentConstraint( lam_coul_on_prev, vfunc( 1 + (4*lam_prev*(v-1.0)), (v/0.75)*(1.0-lam_prev) ) ) )
        system.add( ComponentConstraint( lam_coul_off_prev, vfunc( 1 + (4*(1-lam_prev)*(v-1.0)), (v/0.75)*(1.0-(1-lam_prev)) ) ) )

        system.add( ComponentConstraint( lam_lj_on, Max( Min( 2 * ((1-lam)-0.25), lam_lj_max ), lam_lj_min ) ) ) # scale from 1 to 0 from lam=0.25 to 0.75
        system.add( ComponentConstraint( lam_lj_off, Max( Min( 2 * (lam-0.25), lam_lj_max ), lam_lj_min ) ) )    # scale from 0 to 1 from lam=0.25 to 0.75
        system.add( ComponentConstraint( lam_lj_on_f, Max( Min( 2 * ((1-lam_f)-0.25), lam_lj_max ), lam_lj_min ) ) ) # scale from 1 to 0 from lam=0.25 to 0.75
        system.add( ComponentConstraint( lam_lj_off_f, Max( Min( 2 * (lam_f-0.25), lam_lj_max ), lam_lj_min ) ) )    # scale from 0 to 1 from lam=0.25 to 0.75
        system.add( ComponentConstraint( lam_lj_on_b, Max( Min( 2 * ((1-lam_b)-0.25), lam_lj_max ), lam_lj_min ) ) ) # scale from 1 to 0 from lam=0.25 to 0.75
        system.add( ComponentConstraint( lam_lj_off_b, Max( Min( 2 * (lam_b-0.25), lam_lj_max ), lam_lj_min ) ) )    # scale from 0 to 1 from lam=0.25 to 0.75
        system.add( ComponentConstraint( lam_lj_on_next, Max( Min( 2 * ((1-lam_next)-0.25), lam_lj_max ), lam_lj_min ) ) ) # scale from 1 to 0 from lam=0.25 to 0.75
        system.add( ComponentConstraint( lam_lj_off_next, Max( Min( 2 * (lam_next-0.25), lam_lj_max ), lam_lj_min ) ) )    # scale from 0 to 1 from lam=0.25 to 0.75
        system.add( ComponentConstraint( lam_lj_on_prev, Max( Min( 2 * ((1-lam_prev)-0.25), lam_lj_max ), lam_lj_min ) ) ) # scale from 1 to 0 from lam=0.25 to 0.75
        system.add( ComponentConstraint( lam_lj_off_prev, Max( Min( 2 * (lam_prev-0.25), lam_lj_max ), lam_lj_min ) ) )    # scale from 0 to 1 from lam=0.25 to 0.75

        system.add( ComponentConstraint( lam_coul_swap, Max( v, Max(4.0*(v-1.0)*lam + 1, 4.0*(v-1.0)*(1-lam) + 1 ) ) ) )
        system.add( ComponentConstraint( lam_coul_swap_f, Max( v, Max(4.0*(v-1.0)*lam_f + 1, 4.0*(v-1.0)*(1-lam_f) + 1 ) ) ) )
        system.add( ComponentConstraint( lam_coul_swap_b, Max( v, Max(4.0*(v-1.0)*lam_b + 1, 4.0*(v-1.0)*(1-lam_b) + 1 ) ) ) )
        system.add( ComponentConstraint( lam_coul_swap_next, Max( v, Max(4.0*(v-1.0)*lam_next + 1, 4.0*(v-1.0)*(1-lam_next) + 1 ) ) ) )
        system.add( ComponentConstraint( lam_coul_swap_prev, Max( v, Max(4.0*(v-1.0)*lam_prev + 1, 4.0*(v-1.0)*(1-lam_prev) + 1 ) ) ) )
    else:
        system.add( ComponentConstraint( lam_coul_on, 1-lam ) )
        system.add( ComponentConstraint( lam_coul_off, lam ) )
        system.add( ComponentConstraint( lam_lj_on, Max( Min( 1-lam, lam_lj_max ), lam_lj_min ) ) )
        system.add( ComponentConstraint( lam_lj_off, Max( Min( lam, lam_lj_max ), lam_lj_min ) ) )

        system.add( ComponentConstraint( lam_coul_on_f, 1-lam_f ) )
        system.add( ComponentConstraint( lam_coul_off_f, lam_f ) )
        system.add( ComponentConstraint( lam_lj_on_f, Max( Min( 1-lam_f, lam_lj_max ), lam_lj_min ) ) )
        system.add( ComponentConstraint( lam_lj_off_f, Max( Min( lam_f, lam_lj_max ), lam_lj_min ) ) )

        system.add( ComponentConstraint( lam_coul_on_b, 1-lam_b ) )
        system.add( ComponentConstraint( lam_coul_off_b, lam_b ) )
        system.add( ComponentConstraint( lam_lj_on_b, Max( Min( 1-lam_b, lam_lj_max ), lam_lj_min ) ) )
        system.add( ComponentConstraint( lam_lj_off_b, Max( Min( lam_b, lam_lj_max ), lam_lj_min ) ) )

        system.add( ComponentConstraint( lam_coul_on_next, 1-lam_next ) )
        system.add( ComponentConstraint( lam_coul_off_next, lam_next ) )
        system.add( ComponentConstraint( lam_lj_on_next, Max( Min( 1-lam_next, lam_lj_max ), lam_lj_min ) ) )
        system.add( ComponentConstraint( lam_lj_off_next, Max( Min( lam_next, lam_lj_max ), lam_lj_min ) ) )

        system.add( ComponentConstraint( lam_coul_on_prev, 1-lam_prev ) )
        system.add( ComponentConstraint( lam_coul_off_prev, lam_prev ) )
        system.add( ComponentConstraint( lam_lj_on_prev, Max( Min( 1-lam_prev, lam_lj_max ), lam_lj_min ) ) )
        system.add( ComponentConstraint( lam_lj_off_prev, Max( Min( lam_prev, lam_lj_max ), lam_lj_min ) ) )

    # now add alpha variables that can be used by the EnergyMonitors
    alpha_on = Symbol("alpha_on")
    alpha_off = Symbol("alpha_off")

    system.setConstant(alpha_on, 0)
    system.setConstant(alpha_off, 1)

    system.setConstant( Symbol("alpha_scale"), alpha_scale.val )
    system.add( ComponentConstraint( alpha_on, alpha_scale.val * lam ) )
    system.add( ComponentConstraint( alpha_off, alpha_scale.val * (1-lam) ) )

    # Now constrain alpha to follow lambda
    # First, the reference state (alpha0)
    if use_oldff.val:
        system.add( PropertyConstraint( "alpha0", FFName("free:swap-mobile"), alpha_scale.val * lam ) )
        system.add( PropertyConstraint( "alpha0", FFName("bound:swap-mobile"), alpha_scale.val * (1 - lam) ) )

        system.add( PropertyConstraint( "alpha0", FFName("bound:ligand-mobile"), alpha_scale.val * lam ) )
        system.add( PropertyConstraint( "alpha0", FFName("free:ligand-mobile"), alpha_scale.val * (1 - lam) ) )

        # Now the forwards perturbed state (alpha1)
        system.add( PropertyConstraint( "alpha1", FFName("free:swap-mobile"),  alpha_scale.val * lam_f ) )
        system.add( PropertyConstraint( "alpha1", FFName("bound:swap-mobile"),  alpha_scale.val * (1 - lam_f) ) )

        system.add( PropertyConstraint( "alpha1", FFName("bound:ligand-mobile"),  alpha_scale.val * lam_f ) )
        system.add( PropertyConstraint( "alpha1", FFName("free:ligand-mobile"),  alpha_scale.val * (1 - lam_f) ) )

        # Now the backwards perturbed state (alpha2)
        system.add( PropertyConstraint( "alpha2", FFName("free:swap-mobile"),  alpha_scale.val * lam_b ) )
        system.add( PropertyConstraint( "alpha2", FFName("bound:swap-mobile"),  alpha_scale.val * (1 - lam_b) ) )

        system.add( PropertyConstraint( "alpha2", FFName("bound:ligand-mobile"),  alpha_scale.val * lam_b ) )
        system.add( PropertyConstraint( "alpha2", FFName("free:ligand-mobile"),  alpha_scale.val * (1 - lam_b) ) )

        # Now the next window perturbed state (alpha3)
        system.add( PropertyConstraint( "alpha3", FFName("free:swap-mobile"),  alpha_scale.val * lam_next ) )
        system.add( PropertyConstraint( "alpha3", FFName("bound:swap-mobile"),  alpha_scale.val * (1 - lam_next) ) )

        system.add( PropertyConstraint( "alpha3", FFName("bound:ligand-mobile"),  alpha_scale.val * lam_next ) )
        system.add( PropertyConstraint( "alpha3", FFName("free:ligand-mobile"),  alpha_scale.val * (1 - lam_next) ) )

        # Now the previous window perturbed state (alpha4)
        system.add( PropertyConstraint( "alpha4", FFName("free:swap-mobile"),  alpha_scale.val * lam_prev ) )
        system.add( PropertyConstraint( "alpha4", FFName("bound:swap-mobile"),  alpha_scale.val * (1 - lam_prev) ) )

        system.add( PropertyConstraint( "alpha4", FFName("bound:ligand-mobile"),  alpha_scale.val * lam_prev ) )
        system.add( PropertyConstraint( "alpha4", FFName("free:ligand-mobile"),  alpha_scale.val * (1 - lam_prev) ) )
    else:
        system.add( PropertyConstraint( "alpha", FFName("free:swap-mobile"), alpha_scale.val * lam ) )
        system.add( PropertyConstraint( "alpha", FFName("bound:swap-mobile"), alpha_scale.val * (1 - lam) ) )

        system.add( PropertyConstraint( "alpha", FFName("bound:ligand-mobile"), alpha_scale.val * lam ) )
        system.add( PropertyConstraint( "alpha", FFName("free:ligand-mobile"), alpha_scale.val * (1 - lam) ) )

        # Now the forwards perturbed state (alpha1)
        system.add( PropertyConstraint( "alpha[f]", FFName("free:swap-mobile"),  alpha_scale.val * lam_f ) )
        system.add( PropertyConstraint( "alpha[f]", FFName("bound:swap-mobile"),  alpha_scale.val * (1 - lam_f) ) )

        system.add( PropertyConstraint( "alpha[f]", FFName("bound:ligand-mobile"),  alpha_scale.val * lam_f ) )
        system.add( PropertyConstraint( "alpha[f]", FFName("free:ligand-mobile"),  alpha_scale.val * (1 - lam_f) ) )

        # Now the backwards perturbed state (alpha2)
        system.add( PropertyConstraint( "alpha[b]", FFName("free:swap-mobile"),  alpha_scale.val * lam_b ) )
        system.add( PropertyConstraint( "alpha[b]", FFName("bound:swap-mobile"),  alpha_scale.val * (1 - lam_b) ) )

        system.add( PropertyConstraint( "alpha[b]", FFName("bound:ligand-mobile"),  alpha_scale.val * lam_b ) )
        system.add( PropertyConstraint( "alpha[b]", FFName("free:ligand-mobile"),  alpha_scale.val * (1 - lam_b) ) )

        # Now the next window perturbed state (alpha3)
        system.add( PropertyConstraint( "alpha[next]", FFName("free:swap-mobile"),  alpha_scale.val * lam_next ) )
        system.add( PropertyConstraint( "alpha[next]", FFName("bound:swap-mobile"),  alpha_scale.val * (1 - lam_next) ) )

        system.add( PropertyConstraint( "alpha[next]", FFName("bound:ligand-mobile"),  alpha_scale.val * lam_next ) )
        system.add( PropertyConstraint( "alpha[next]", FFName("free:ligand-mobile"),  alpha_scale.val * (1 - lam_next) ) )

        # Now the previous window perturbed state (alpha4)
        system.add( PropertyConstraint( "alpha[prev]", FFName("free:swap-mobile"),  alpha_scale.val * lam_prev ) )
        system.add( PropertyConstraint( "alpha[prev]", FFName("bound:swap-mobile"),  alpha_scale.val * (1 - lam_prev) ) )

        system.add( PropertyConstraint( "alpha[prev]", FFName("bound:ligand-mobile"),  alpha_scale.val * lam_prev ) )
        system.add( PropertyConstraint( "alpha[prev]", FFName("free:ligand-mobile"),  alpha_scale.val * (1 - lam_prev) ) )

    # Now soften the swap-water-swap-water interactions around lambda = 0.5 (used if decharging the ligand)
    if uncharge_ligand.val:
        s_scl = soften_water.val
    else:
        s_scl = 0

    if use_oldff.val:
        system.add( PropertyConstraint( "alpha0", FFName("swap:interclj"), s_scl * (1 - 2*Abs(lam - 0.5)) ) )
        system.add( PropertyConstraint( "alpha1", FFName("swap:interclj"), s_scl * (1 - 2*Abs(lam_f - 0.5)) ) )
        system.add( PropertyConstraint( "alpha2", FFName("swap:interclj"), s_scl * (1 - 2*Abs(lam_b - 0.5)) ) )
        system.add( PropertyConstraint( "alpha3", FFName("swap:interclj"), s_scl * (1 - 2*Abs(lam_next - 0.5)) ) )
        system.add( PropertyConstraint( "alpha4", FFName("swap:interclj"), s_scl * (1 - 2*Abs(lam_prev - 0.5)) ) )
    else:
        system.add( PropertyConstraint( "alpha", FFName("swap:interclj"), s_scl * (1 - 2*Abs(lam - 0.5)) ) )
        system.add( PropertyConstraint( "alpha[f]", FFName("swap:interclj"), s_scl * (1 - 2*Abs(lam_f - 0.5)) ) )
        system.add( PropertyConstraint( "alpha[b]", FFName("swap:interclj"), s_scl * (1 - 2*Abs(lam_b - 0.5)) ) )
        system.add( PropertyConstraint( "alpha[next]", FFName("swap:interclj"), s_scl * (1 - 2*Abs(lam_next - 0.5)) ) )
        system.add( PropertyConstraint( "alpha[prev]", FFName("swap:interclj"), s_scl * (1 - 2*Abs(lam_prev - 0.5)) ) )

    # Now lets create all of the groups for moves based on the above

    # All solvent molecules in the bound and free legs are moved together
    mobile_solvent = MoleculeGroup("mobile_solvent")
    mobile_solvent.add( mobile_bound_solvents_group.molecules() )
    mobile_solvent.add( mobile_free_water_group.molecules() )

    system.add( mobile_solvent )

    # All protein sidechains are moved together
    mobile_sidechains = MoleculeGroup("mobile_sidechains")
    mobile_sidechains.add(mobile_bound_protein_sidechains_group.molecules())
    
    system.add( mobile_sidechains )

    # All protein backbones are moved together
    mobile_backbones = MoleculeGroup("mobile_backbones")
    mobile_backbones.add(mobile_bound_protein_backbones_group.molecules())

    system.add( mobile_backbones )

    # All solute molecules are moved together
    mobile_solutes = MoleculeGroup("mobile_solutes")
    mobile_solutes.add(mobile_bound_solutes_group.molecules())

    system.add( mobile_solutes )

    # The ligand is moved in its own group
    mobile_ligand = MoleculeGroup("mobile_ligand")
    mobile_ligand.add(ligand_mol)

    system.add( mobile_ligand )

    # The swap water cluster is moved in its own group
    mobile_swap = MoleculeGroup("mobile_swap_water")
    mobile_swap.add(swap_water_group.molecules())

    system.add( mobile_swap )

    if use_identity_constraint:
        print("Adding the identity constraint...")

        # Now add the constraint that keeps the identities of the
        # swap molecules. The swap molecules are chosen from all available mobile
        # water molecules. We need to build a group of all mobile water molecules that
        # are waters (as opposed to ions, as other molecules may be in mobile_solvent)
        mobile_water = MoleculeGroup("mobile_water")

        # The mobile water *must* contain the swap waters, so that they can be swapped
        mobile_water.add(swap_water_group)
        mobile_water.add(mobile_free_water_group)

        if waterbox_only.val:
            print("Choosing water molecules only from the free water box.")
        else:
            print("Choosing swap waters from both the protein box and water box.")
            mobile_water.add(mobile_bound_water_group)

        print("The number of candidates for the swap water equals: %d" % mobile_water.nMolecules())

        system.add(mobile_water)
        system.add( IdentityConstraint(identity_points, mobile_water, { "space" : Cartesian() } ) )

        # Apply all of the constraints
        system.applyConstraints()
    else:
        print("Using the reflection volume to hold the swap water in place.")
        system.applyConstraints()

    ###
    ### ADD THE SYSTEM MONITORS
    ###

    # Now we need to add the monitors...
    print("\nAdding system monitors...")

    system.add( "delta_g^{F}", MonitorComponent( Symbol("delta_nrg^{F}"),
                                                 FreeEnergyAverage(temperature.val,
                                                                   dlam * binwidth.val) ) )

    system.add( "delta_g^{B}", MonitorComponent( Symbol("delta_nrg^{B}"),
                                                 FreeEnergyAverage(temperature.val,
                                                                   dlam * binwidth.val, False) ) )

    system.add( "delta_g^{next}", MonitorComponent( Symbol("delta_nrg^{next}"),
                                                    BennettsFreeEnergyAverage(0 * kcal_per_mol,
                                                                              temperature.val,
                                                                              0.1 * binwidth.val) ) )

    system.add( "delta_g^{prev}", MonitorComponent( Symbol("delta_nrg^{prev}"),
                                                    BennettsFreeEnergyAverage(0 * kcal_per_mol,
                                                                              temperature.val,
                                                                              0.1 * binwidth.val, False) ) )
    
    system.add( "delta_bound_g^{F}", MonitorComponent( Symbol("delta_bound_nrg^{F}"),
                                                       FreeEnergyAverage(temperature.val,
                                                                         dlam * binwidth.val) ) )
    system.add( "delta_bound_g^{B}", MonitorComponent( Symbol("delta_bound_nrg^{B}"),
                                                       FreeEnergyAverage(temperature.val,
                                                                         dlam * binwidth.val, False) ) )
    system.add( "delta_free_g^{F}", MonitorComponent( Symbol("delta_free_nrg^{F}"),
                                                      FreeEnergyAverage(temperature.val,
                                                                        dlam * binwidth.val) ) )
    system.add( "delta_free_g^{B}", MonitorComponent( Symbol("delta_free_nrg^{B}"),
                                                      FreeEnergyAverage(temperature.val,
                                                                        dlam * binwidth.val, False) ) )

    # we will monitor the average energy between the swap cluster/ligand and each
    # residue with mobile sidechain, and each mobile solute
    monitor_prosol = None

    if mobile_solutes.isEmpty():
        monitor_prosol = mobile_sidechains
    elif mobile_sidechains.isEmpty():
        monitor_prosol = mobile_solutes
    else:
        monitor_prosol = MoleculeGroup("monitored_protein_solute")
        monitor_prosol.add(mobile_sidechains)
        monitor_prosol.add(mobile_solutes)
        system.add(monitor_prosol)

    residue_nrgmon = FreeEnergyMonitor(monitor_prosol, ligand_group, mobile_swap)

    nrgmons = {}
    nrgmons["residue_nrgmon"] = residue_nrgmon

    # because the water molecules can diffuse, we find all waters within
    # a certain distance of the ligand, and then identify them using identity
    # points (placed at the center of the initial positions of the waters),
    # and then monitor those...
    boundwater_points = []
    freewater_points = []

    if water_monitor_distance.val:
        dist = water_monitor_distance.val.to(angstrom)

        for molnum in mobile_bound_water_group.molNums():
            water_mol = mobile_bound_water_group[molnum].molecule()
            if getMinimumDistance(ligand_mol,water_mol) < dist:
                # we should monitor this water
                boundwater_points.append( VectorPoint(water_mol.evaluate().center()) )
    
        for molnum in mobile_free_water_group.molNums():
            #this is a mobile water, so a candidate for monitoring
            water_mol = mobile_free_water_group[molnum].molecule()
            if getMinimumDistance(ligand_mol,water_mol) < dist:
                # we should monitor this water
                freewater_points.append( VectorPoint(water_mol.evaluate().center()) )

        system.add(mobile_bound_water_group)
        system.add(mobile_free_water_group)

        boundwater_assigner = IDAssigner(boundwater_points, mobile_bound_water_group,
                                         {"space" : Cartesian()})

        boundwater_assigner.update(system)

        freewater_assigner = IDAssigner(freewater_points, mobile_free_water_group,
                                        {"space" : Cartesian()})

        freewater_assigner.update(system)

        boundwater_nrgmon = FreeEnergyMonitor(boundwater_assigner, ligand_group, mobile_swap)
        freewater_nrgmon = FreeEnergyMonitor(freewater_assigner, mobile_swap, ligand_group)

        nrgmons["boundwater_nrgmon"] = boundwater_nrgmon
        nrgmons["freewater_nrgmon"] = freewater_nrgmon

    for key in list(nrgmons.keys()):
        nrgmons[key].setCoulombPower(coul_power.val)
        nrgmons[key].setShiftDelta(shift_delta.val)
        nrgmons[key].setTemperature(temperature.val)

        system.add(key, nrgmons[key], nrgmon_frequency.val)

    # now calculate the total energy of the system - this initialises grids etc.
    # ensuring that, when we make the replicas, the maximum amount of sharing between
    # replicas occurs
    print("\nEnergies of this system at lambda == 0...")
    system.setConstant(lam, 0.0)
    printEnergies(system.energies(), sys.stdout)

    print("\nEnergies of this system at lambda == 1...")
    system.setConstant(lam, 1.0)
    printEnergies(system.energies(), sys.stdout)

    system.setConstant(lam, 0.0)

    return system


def makeRETI(system, moves):
    """This function replicates 'system' over each of the supplied lambda values
       and uses 'moves' to sample each of the replicated systems. This uses RETI
       to perform replica exchange moves across lambda"""

    lam = Symbol("lambda")

    lamvals = getLambdaValues()

    replicas = Replicas( len(lamvals) )

    replicas.setSubSystem(system)
    replicas.setSubMoves(moves)
    replicas.setNSubMoves(nsubmoves.val)
    replicas.setLambdaComponent(lam)
    replicas.setRecordAllStatistics(True)        

    seed = random_seed.val
    
    if seed is None:
        seed = RanGenerator().randInt(100000,1000000)
        print("RETI system using generated random number seed %d" % seed)
    else:
        print("RETI system using supplied random number seed %d" % seed)
    
    replicas.setGenerator( RanGenerator(seed+5) )

    for i in range(0, len(lamvals)):
        # set the initial lambda value for this replica
        replicas.setLambdaValue(i, lamvals[i])

    for i in range(0, len(lamvals)):
        print(lamvals[i])
        print(replicas[i].subSystem().constants())

    # Now add monitors for each replica that will copy back
    nrgmons = [ "delta_g^{F}", "delta_g^{B}", "delta_g^{next}", "delta_g^{prev}",
                "delta_bound_g^{F}", "delta_bound_g^{B}",
                "delta_free_g^{F}", "delta_free_g^{B}",
                "residue_nrgmon", "boundwater_nrgmon", "freewater_nrgmon" ]

    for nrgmon in nrgmons:
        replicas.add( nrgmon, MonitorMonitor(MonitorName(nrgmon), True) )

    # now create the replica exchange moves for the replicas
    replica_moves = RepExMove2()
    replica_moves.setDisableSwaps(False)
    replica_moves.setGenerator( RanGenerator(seed+7) )

    print("\nReturning the WSRC RETI replicas and moves...")
    return (replicas, replica_moves)


def getName(view):
   try:
       residue = view.residue()
       return "%s:%s" % (residue.name().value(), residue.number().value())
   except:
       return "%s:%s" % (view.name().value(), view.number().value())


def loadWSRC():
    # Load the WSRC system and moves using the passed parameters
    # This returns (wsrc_system, wsrc_moves), ready for simulation

    print("Loading the protein box system...")

    if os.path.exists(protein_s3file.val):
        print("Loading existing s3 file %s..." % protein_s3file.val)
        proteinsys = Sire.Stream.load(protein_s3file.val)

    else:
        print("Loading from Amber files %s / %s..." % (protein_topfile.val, protein_crdfile.val))
        # Add the name of the ligand to the list of solute molecules
        proteinsys_scheme = NamingScheme()
        proteinsys_scheme.addSoluteResidueName( ligand_name.val )

        # Load up the system. This will automatically find the protein, solute, water, solvent
        # and ion molecules and assign them to different groups
        proteinsys = createSystem(protein_topfile.val, protein_crdfile.val, proteinsys_scheme)
        ligand_mol = findMolecule(proteinsys, ligand_name.val)

        if ligand_mol is None:
            print("Cannot find the ligand (%s) in the set of loaded molecules!" % ligand_name.val)
            sys.exit(-1)

        # Center the system with the ligand at (0,0,0)
        proteinsys = centerSystem(proteinsys, ligand_mol)
        ligand_mol = proteinsys[ligand_mol.number()][0].molecule()

        proteinsys = addFlexibility(proteinsys, Vector(0,0,0), reflection_radius.val, proteinsys_scheme )
        Sire.Stream.save(proteinsys, protein_s3file.val)

    ligand_mol = findMolecule(proteinsys, ligand_name.val)

    if ligand_mol is None:
        print("Cannot find the ligand (%s) in the set of loaded molecules!" % ligand_name.val)
        sys.exit(-1)

    print("Loading the water box system...")

    if os.path.exists(water_s3file.val):
        print("Loading from existing s3 file %s..." % water_s3file.val)
        watersys = Sire.Stream.load(water_s3file.val)
    else:
        print("Loading from Amber files %s / %s..." % (water_topfile.val, water_crdfile.val))
        watersys = createSystem(water_topfile.val, water_crdfile.val)
        watersys = addFlexibility(watersys, Vector(0,0,0), reflection_radius.val)
        Sire.Stream.save(watersys, water_s3file.val)

    print("Building the WSRC forcefields")

    wsrc_system = mergeSystems(proteinsys, watersys, ligand_mol)
    wsrc_moves = createWSRCMoves(wsrc_system)

    return (wsrc_system, wsrc_moves)


def printComponents(comps, FILE):
    """This function prints out all of the free energy components in the passed object"""
    print("RESIDUE    TOTAL    COULOMB    LJ", file=FILE)
    for i in range(0, comps.nComponents()):
        print("%s  %s  %s  %s" % (comps.viewAt(i).residue(), \
                                  -comps.integrate(i).values()[-1].y(), \
                                  -comps.integrateCoulomb(i).values()[-1].y(), \
                                  -comps.integrateLJ(i).values()[-1].y()), file=FILE)


def printFreeEnergy(total, bound, free, FILE):
    """This function prints out the total, bound and free free energies"""
    print("TOTAL   BOUND    FREE", file=FILE)
    print("%s   %s   %s" % (-total.integrate().values()[-1].y(), \
                            -bound.integrate().values()[-1].y(), \
                            -free.integrate().values()[-1].y()), file=FILE)


def analyseWSRC(replicas, iteration, bennetts_freenrgs, fep_freenrgs, ti_freenrgs, bound_freenrgs, free_freenrgs,
                res_freenrgs, bound_water_freenrgs, free_water_freenrgs):
    """This function is used to perform all analysis of iteration 'it' of the passed WSRC system"""

    print("Analysing iteration %d..." % iteration)

    if not os.path.exists(outdir.val):
        os.makedirs(outdir.val)

    # read the value of delta_lambda from the first system
    system = replicas[0].subSystem()
    delta_lambda = system.constant(Symbol("delta_lambda"))

    logfile = "%s/results_%0004d.log" % (outdir.val, iteration)

    FILE = open(logfile, "w")

    print("===========================", file=FILE)
    print(" Results for iteration %d" % iteration, file=FILE)
    print("===========================", file=FILE)

    print("\ndelta_lambda == %f" % delta_lambda, file=FILE)
    print("temperature == %f K\n" % replicas[0].subMoves().temperature().to(kelvin), file=FILE) 

    nreplicas = replicas.nReplicas()

    # extract all of the monitors from the replicas
    lambda_values = []

    dg_f = {}
    dg_b = {}

    dg_next = {}
    dg_prev = {}

    dg_bound_f = {}
    dg_bound_b = {}
    
    dg_free_f = {}
    dg_free_b = {}

    dg_residue = {}
    dg_boundwater = {}
    dg_freewater = {}

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
                bound_leg = system[MGName("bound_leg")]
                free_leg = system[MGName("free_leg")]

                PDB().write(bound_leg, "%s/bound_mobile_%000006d_%.5f.pdb" % (outdir.val, iteration, lamval))
                PDB().write(free_leg, "%s/free_mobile_%000006d_%.5f.pdb" % (outdir.val, iteration, lamval))

        dg_f[lamval] = monitors[MonitorName("delta_g^{F}")][-1].accumulator()
        dg_b[lamval] = monitors[MonitorName("delta_g^{B}")][-1].accumulator()
        dg_next[lamval] = monitors[MonitorName("delta_g^{next}")][-1].accumulator()
        dg_prev[lamval] = monitors[MonitorName("delta_g^{prev}")][-1].accumulator()
        dg_bound_f[lamval] = monitors[MonitorName("delta_bound_g^{F}")][-1].accumulator()
        dg_bound_b[lamval] = monitors[MonitorName("delta_bound_g^{B}")][-1].accumulator()
        dg_free_f[lamval] = monitors[MonitorName("delta_free_g^{F}")][-1].accumulator()
        dg_free_b[lamval] = monitors[MonitorName("delta_free_g^{B}")][-1].accumulator()

        dg_residue[lamval] = monitors[MonitorName("residue_nrgmon")][-1]
        dg_boundwater[lamval] = monitors[MonitorName("boundwater_nrgmon")][-1]
        dg_freewater[lamval] = monitors[MonitorName("freewater_nrgmon")][-1]

    windows = copy.deepcopy(lambda_values)
    windows.sort()

    if windows[-1] != 1:
        windows.append(1)

    if windows[0] != 0:
        windows.insert(0,0)

    bennetts_freenrgs.set( iteration, windows, dg_next, dg_prev )
    fep_freenrgs.set( iteration, windows, dg_next, dg_prev )
    ti_freenrgs.set( iteration, dg_f, dg_b, delta_lambda )

    bound_freenrgs.set( iteration, dg_bound_f, dg_bound_b, delta_lambda )
    free_freenrgs.set( iteration, dg_free_f, dg_free_b, delta_lambda )

    print("\nTOTAL BINDING FREE ENERGY\n", file=FILE)
    printFreeEnergy(ti_freenrgs[iteration], bound_freenrgs[iteration], free_freenrgs[iteration], FILE)

    res_freenrgs.set( iteration, dg_residue )
    bound_water_freenrgs.set( iteration, dg_boundwater )
    free_water_freenrgs.set( iteration, dg_freewater )

    print("\nRESIDUE FREE ENERGY COMPONENTS\n", file=FILE)
    printComponents(res_freenrgs[iteration], FILE)

    print("\nPROTEIN BOX WATER FREE ENERGY COMPONENTS\n", file=FILE)
    printComponents(bound_water_freenrgs[iteration], FILE)

    print("\nWATER BOX WATER FREE ENERGY COMPONENTS\n", file=FILE)
    printComponents(free_water_freenrgs[iteration], FILE)

    print("\n=============", file=FILE)
    print("Binding free energy for iteration %d equals %s" % (iteration, \
                        -ti_freenrgs[iteration].integrate().values()[-1].y()), file=FILE)
    print("==============", file=FILE)


@resolveParameters
def run():
    """This is a very high level function that does everything to run a WSRC simulation"""

    t = QTime()
    total_t = QTime()
    total_t.start()

    if os.path.exists(restart_file.val):
        t.start()
        (wsrc_system, wsrc_moves) = Sire.Stream.load(restart_file.val)
        print("Loading the restart file took %d ms" % t.elapsed())
    else:
        # Load the WSRC protein and water boxes from the topology and coordinate
        # files and merge together into the WSRC system and moves object
        t.start()
        if os.path.exists(sysmoves_file.val):
            (wsrc_system, wsrc_moves) = Sire.Stream.load(sysmoves_file.val)
        else:
            (wsrc_system, wsrc_moves) = loadWSRC()

            Sire.Stream.save( (wsrc_system, wsrc_moves), "pre_equil.s3")

            # Should add in some equilibration here...
            if nequilmoves.val:
                print("Equilibrating the system (number of moves: %d)..." % nequilmoves.val)
                wsrc_system = wsrc_moves.move(wsrc_system, nequilmoves.val, False)
                print("...equilibration complete")
            
            Sire.Stream.save( (wsrc_system, wsrc_moves), sysmoves_file.val)

        # Now replicate the WSRC system across all lambda values so that we
        # can run a RETI simulation
        (wsrc_system, wsrc_moves) = makeRETI(wsrc_system, wsrc_moves)

        Sire.Stream.save( (wsrc_system, wsrc_moves), restart_file.val )

        print("Initialising the simulation took %d ms" % t.elapsed())

    # see how many blocks of moves we still need to perform...
    nattempted = wsrc_moves.nMoves()

    print("Number of iterations to perform: %d. Number of iterations completed: %d." % (nmoves.val, nattempted))

    # See if we have any existing free energy statistics files...
    t.start()
    freenrgs_file = "%s/freenrgs.s3" % outdir.val

    if not os.path.exists(freenrgs_file):
        bennetts_freenrgs = Bennetts()
        fep_freenrgs = FEP()
        ti_freenrgs = TI()
    else:
        [bennetts_freenrgs, fep_freenrgs, ti_freenrgs] = Sire.Stream.load(freenrgs_file)

    freenrg_parts_file = "%s/freenrg_parts.s3" % outdir.val

    if not os.path.exists(freenrg_parts_file):
        bound_freenrgs = TI()
        free_freenrgs = TI()
    else:
        [bound_freenrgs, free_freenrgs] = Sire.Stream.load(freenrg_parts_file)

    freenrg_components_file = "%s/freenrg_components.s3" % outdir.val

    if not os.path.exists(freenrg_components_file):
        res_freenrgs = TIComponents()
        bound_water_freenrgs = TIComponents()
        free_water_freenrgs = TIComponents()
    else:
        [res_freenrgs, bound_water_freenrgs, free_water_freenrgs] = Sire.Stream.load(freenrg_components_file)

    print("Initialising / loading the free energy files took %d ms" % t.elapsed())

    for i in range(nattempted+1, nmoves.val+1):
        t.start()
        print("Performing iteration %d..." % i)
        wsrc_moves.move(wsrc_system, 1, True)

        print("...iteration complete (took %d ms)" % t.elapsed())

        t.start()
        print("Analysing iteration %d..." % i)
        analyseWSRC(wsrc_system, i, bennetts_freenrgs, fep_freenrgs, ti_freenrgs, bound_freenrgs, free_freenrgs,
                    res_freenrgs, bound_water_freenrgs, free_water_freenrgs)
        wsrc_system.clearAllStatistics()
        print("...analysis complete (took %d ms)" % t.elapsed())

        if i % restart_frequency.val == 0 or i == nmoves.val:
            t.start()
            print("Saving the free energy analysis files from iteration %d..." % i)
            # save the old file to a backup
            try:
                shutil.copy(freenrgs_file, "%s.bak" % freenrgs_file)
            except:
                pass

            try:
                shutil.copy(freenrg_components_file, "%s.bak" % freenrg_components_file)
            except:
                pass

            try:
                shutil.copy(freenrg_parts_file, "%s.bak" % freenrg_parts_file)
            except:
                pass

            Sire.Stream.save( [bennetts_freenrgs, fep_freenrgs, ti_freenrgs], freenrgs_file )
            Sire.Stream.save( [bound_freenrgs, free_freenrgs], freenrg_parts_file )
            Sire.Stream.save( [res_freenrgs, bound_water_freenrgs, free_water_freenrgs], freenrg_components_file )
            print("...save complete (took %d ms)" % t.elapsed())

        # write a restart file every N moves in case of crash or run out of time
        if i % restart_frequency.val == 0 or i == nmoves.val:
            t.start()
            print("Saving the restart file from iteration %d..." % i)
            # save the old file to a backup
            try:
                shutil.copy(restart_file.val, "%s.bak" % restart_file.val)
            except:
                pass

            Sire.Stream.save( (wsrc_system, wsrc_moves), restart_file.val )
            print("...save complete (took %d ms)" % t.elapsed())

    print("All iterations complete. Total runtime was %d ms" % total_t.elapsed())
