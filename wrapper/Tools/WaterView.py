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
from Sire.Cluster import *
from Sire.Units import *

import Sire.Config
import Sire.Stream

from Sire.Tools.AmberLoader import *
from Sire.Tools import Parameter, resolveParameters

import os
import shutil
import copy

# We will use the waterbox held in the WSRC tools directory
wsrc_tools_dir = "%s/Tools/WSRC" % Sire.Config.share_directory

####################################################
# ALL OF THE GLOBAL USER-AVAILABLE LSRC PARAMETERS #
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

use_fast_ff = Parameter("fast forcefield", False, """Whether or not to use the experimental fast forcefield""")

temperature = Parameter("temperature", 25*celsius, """Simulation temperature""")
random_seed = Parameter("random seed", None, """Random number seed. Set this if you
                         want to have reproducible simulations.""")

overlap_dist = Parameter("overlap distance", 1*angstrom,
                         """Maximum distance between the ligand and an overlapping water molecule""")

reflection_radius = Parameter("reflection radius", 20*angstrom,
                              """The radius of the reflection sphere""")

ligand_name = Parameter("ligand name", None,
                         """The name of ligand. This should be the name of one of the residues
                            in the ligand, so that this program can find the correct molecule. If it is not set, then 
                            the first non-protein, non solvent molecule is used.""")

topfile = Parameter("topfile", "complex.top",
                     """Name of the topology file containing the solvated protein-ligand complex.""")

crdfile = Parameter("crdfile", "complex.crd",
                     """Name of the coordinate file containing the coordinates of the 
                        solvated protein-ligand complex.""")

s3file = Parameter("s3file", "complex.s3",
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

restart_file = Parameter("restart file", "restart.s3",
                          """Name of the restart file to use to save progress during calculation.""")

nmoves = Parameter("nmoves", 100, """Number of blocks of moves to perform during the simulation.""")

nsubmoves = Parameter("nsubmoves", 50000,
                      """The number of moves to perform in each block.""")

nequilmoves = Parameter("nequilmoves", 50000,
                        """Number of equilibration moves to perform before monitoring the water distribution.""")

save_pdb = Parameter("save pdb", True,
                     """Whether or not to write a PDB of the system after each block.""")

restart_frequency = Parameter("restart frequency", 10,
                              """The frequency (number of blocks between) saving the restart file for the simulation.""")


####################################################


def setCLJProperties(forcefield):
    if cutoff_method.val.find("shift electrostatics") != -1:
        forcefield.setShiftElectrostatics(True)

    elif cutoff_method.val.find("reaction field") != -1:
        forcefield.setUseReactionField(True)
        forcefield.setReactionFieldDielectric(rf_dielectric.val)

    else:
        print("Cannot interpret the cutoff method from \"%s\"" % cutoff_method.val, file=sys.stderr)

    forcefield.setSpace(Cartesian())
    forcefield.setSwitchingFunction( HarmonicSwitchingFunction(coul_cutoff.val,coul_cutoff.val,
                                                               lj_cutoff.val,lj_cutoff.val) )

    return forcefield


def setFakeGridProperties(forcefield):
    forcefield.setSwitchingFunction( HarmonicSwitchingFunction(coul_cutoff.val,coul_cutoff.val,
                                                               lj_cutoff.val,lj_cutoff.val) )
    forcefield.setSpace(Cartesian())

    return forcefield


def setGridProperties(forcefield, extra_buffer=0*angstrom):
    forcefield.setGridSpacing(grid_spacing.val)
    forcefield.setBuffer(grid_buffer.val + extra_buffer)
    forcefield.setLJCutoff(lj_cutoff.val)
    forcefield.setCoulombCutoff(coul_cutoff.val)

    return forcefield


def loadWater():
    """Load the the water box to substitute for the ligand"""

    if os.path.exists(water_s3file.val):
        print("Restoring from Sire Streamed Save file %s..." % water_s3file.val)
        watersys = Sire.Stream.load(water_s3file.val)
    else:
        print("Loading from Amber files %s / %s..." % (water_topfile.val, water_crdfile.val))
        watersys = createSystem(water_topfile.val, water_crdfile.val)
        watersys = addFlexibility(watersys, Vector(0,0,0), reflection_radius.val)
        Sire.Stream.save(watersys, water_s3file.val)

    return watersys


def loadSystem(topfile, crdfile, s3file, ligand_name):
    """Load the solvated protein-ligand system from topfile/crdfile, saving into the s3 file 's3file'
       and locating the ligand molecule called 'ligand_name'"""

    if os.path.exists(s3file):
        print("Restoring from Sire Streamed Save file %s..." % s3file)
        system = Sire.Stream.load(s3file)
    else:
        print("Loading from Amber files %s / %s..." % (topfile,crdfile))
        # Add the name of the ligand to the list of solute molecules
        scheme = NamingScheme()

        if ligand_name:
            scheme.addSoluteResidueName(ligand_name)

        # Load up the system. This will automatically find the protein, solute, water, solvent
        # and ion molecules and assign them to different groups
        system = createSystem(topfile, crdfile, scheme)
        ligand_mol = findMolecule(system, ligand_name)

        if ligand_mol is None:
            print("Cannot find the ligand (%s) in the set of loaded molecules!" % ligand_name)
            sys.exit(-1)

        # Center the system with the ligand at (0,0,0)
        system = centerSystem(system, ligand_mol)
        ligand_mol = system[ligand_mol.number()][0].molecule()

        system = addFlexibility(system, Vector(0,0,0), reflection_radius.val, scheme )
        Sire.Stream.save(system, s3file)

    ligand_mol = findMolecule(system, ligand_name)

    if ligand_mol is None:
        print("Cannot find the ligand (%s) in the set of loaded molecules!" % ligand_name)
        sys.exit(-1)

    return (system, ligand_mol)


def getMinimumDistance(mol0, mol1):
    space = Cartesian()
    return space.minimumDistance(CoordGroup(mol0.molecule().property("coordinates").array()), \
                                 CoordGroup(mol1.molecule().property("coordinates").array()))


def findOverlappingWaters(ligand_mol, watersys):
    """Find all of the water molecules in watersys that overlap with ligand_mol"""
    waters = Molecules()

    if MGName("mobile_solvents") in watersys.mgNames():
        mols = watersys[MGName("mobile_solvents")].molecules()
        for molnum in mols.molNums():
            # only add this water if it overlaps with the ligand
            water_mol = mols[molnum].molecule().edit().renumber().commit()

            if getMinimumDistance(ligand_mol,water_mol) < overlap_dist.val.to(angstrom):
                for j in range(0,water_mol.nResidues()):
                    water_mol = water_mol.residue( ResIdx(j) ).edit() \
                                                   .setProperty( PDB.parameters().pdbResidueName(), "SWP" ) \
                                                   .commit().molecule()

                waters.add(water_mol)

    return waters

def makeSim(system, ligand_mol, watersys):
    """Create simulation systems with and without the ligand and return those systems together
       with the moves"""

    stage1 = System("with_ligand")
    stage2 = System("without_ligand")

    if system.containsProperty("reflection center"):
        reflection_center = system.property("reflection center").toVector()[0]
        reflection_radius = float(str(system.property("reflection sphere radius")))

        stage1.setProperty("reflection center", AtomCoords(CoordGroup(1,reflection_center)))
        stage1.setProperty("reflection sphere radius", VariantProperty(reflection_radius))

        stage2.setProperty("reflection center", AtomCoords(CoordGroup(1,reflection_center)))
        stage2.setProperty("reflection sphere radius", VariantProperty(reflection_radius))

    # create a molecule group for fixed atoms (everything except the mobile water)
    fixed_group = MoleculeGroup("fixed")

    if MGName("fixed_molecules") in system.mgNames():
        fixed_group.add( system[ MGName("fixed_molecules") ] )

    if MGName("mobile_solutes") in system.mgNames():
        fixed_group.add( system[MGName("mobile_solutes")] )

    if MGName("protein_sidechains") in system.mgNames() or \
       MGName("protein_backbones") in system.mgNames():

        all_proteins = Molecules()

        try:
            protein_sidechains = system[MGName("protein_sidechains")]
            all_proteins.add(protein_sidechains.molecules())
        except:
            pass

        try:
            protein_backbones = system[MGName("protein_backbones")]
            all_proteins.add(protein_backbones.molecules())
        except:
            pass

        try:
            boundary_molecules = system[MGName("boundary_molecules")]
            all_proteins.add(boundary_molecules.molecules())
        except:
            pass

        for molnum in all_proteins.molNums():
            protein_mol = Molecule.join(all_proteins[molnum])
            fixed_group.add(protein_mol)

    stage1_fixed_group = MoleculeGroup(fixed_group)
    stage2_fixed_group = MoleculeGroup(fixed_group)

    stage1_fixed_group.add(ligand_mol)
    stage2_fixed_group.remove(ligand_mol)

    mobile_group = MoleculeGroup("mobile_group")
    if MGName("mobile_solvents") in system.mgNames():
        mobile_group.add( system[MGName("mobile_solvents")] )

    stage1_mobile_group = MoleculeGroup(mobile_group)
    stage2_mobile_group = MoleculeGroup(mobile_group)

    # now find water molecules from the water system that can be substituted for the ligand
    watermols = findOverlappingWaters(ligand_mol, watersys)

    stage2_mobile_group.add(watermols)

    print("The number of stage 1 fixed non-solvent molecules is %d." % stage1_fixed_group.nMolecules())
    print("The number of stage 1 mobile solvent molecules is %d." % stage1_mobile_group.nMolecules())

    print("The number of stage 2 fixed non-solvent molecules is %d." % stage2_fixed_group.nMolecules())
    print("The number of stage 2 mobile solvent molecules is %d." % stage2_mobile_group.nMolecules())

    # write a PDB of all of the fixed molecules
    PDB().write(stage1_mobile_group, "stage1_mobile_atoms.pdb")
    PDB().write(stage2_mobile_group, "stage2_mobile_atoms.pdb")
    PDB().write(stage1_fixed_group, "stage1_fixed_atoms.pdb")
    PDB().write(stage2_fixed_group, "stage2_fixed_atoms.pdb")

    # create the forcefields

    if use_fast_ff.val:
        stage1_ff = InterFF("ff")
        stage2_ff = InterFF("ff")
        stage1_ff.setCLJFunction( CLJShiftFunction(Cartesian(), coul_cutoff.val, lj_cutoff.val) )
        stage2_ff.setCLJFunction( CLJShiftFunction(Cartesian(), coul_cutoff.val, lj_cutoff.val) )
        
        if disable_grid.val:
            stage1_ff.disableGrid()
            stage2_ff.disableGrid()
        else:
            stage1_ff.enableGrid()
            stage1_ff.setGridSpacing(grid_spacing.val)
            stage1_ff.setGridBuffer(grid_buffer.val)
            stage2_ff.enableGrid()
            stage2_ff.setGridSpacing(grid_spacing.val)
            stage2_ff.setGridBuffer(grid_buffer.val)

        stage1_ff.add(stage1_mobile_group)
        stage1_ff.setFixedAtoms(stage1_fixed_group.molecules())
        stage2_ff.add(stage2_mobile_group)
        stage2_ff.setFixedAtoms(stage2_fixed_group.molecules())

        stage1.add(stage1_ff)
        stage1.setComponent(stage1.totalComponent(), stage1_ff.components().total())
        stage2.add(stage1_ff)
        stage2.setComponent(stage2.totalComponent(), stage2_ff.components().total())

    else:
        # forcefield holding the energy between the mobile atoms and  
        # the fixed atoms
        if disable_grid.val:
            stage1_mobile_fixed = InterGroupCLJFF("mobile-fixed")
            stage1_mobile_fixed = setCLJProperties(stage1_mobile_fixed)
            stage1_mobile_fixed = setFakeGridProperties(stage1_mobile_fixed)
            stage1_mobile_fixed.add(stage1_mobile_group, MGIdx(0))
            stage1_mobile_fixed.add(stage1_fixed_group, MGIdx(1))

            stage2_mobile_fixed = InterGroupCLJFF("mobile-fixed")
            stage2_mobile_fixed = setCLJProperties(stage2_mobile_fixed)
            stage2_mobile_fixed = setFakeGridProperties(stage2_mobile_fixed)
            stage2_mobile_fixed.add(stage2_mobile_group, MGIdx(0))
            stage2_mobile_fixed.add(stage2_fixed_group, MGIdx(1))
        else:
            stage1_mobile_fixed = GridFF2("mobile-fixed")
            stage1_mobile_fixed = setCLJProperties(stage1_mobile_fixed)
            stage1_mobile_fixed = setGridProperties(stage1_mobile_fixed)

            stage1_mobile_fixed.add(stage1_mobile_group, MGIdx(0))
            stage1_mobile_fixed.addFixedAtoms(stage1_fixed_group)

            stage2_mobile_fixed = GridFF2("mobile-fixed")
            stage2_mobile_fixed = setCLJProperties(stage2_mobile_fixed)
            stage2_mobile_fixed = setGridProperties(stage2_mobile_fixed)

            stage2_mobile_fixed.add(stage2_mobile_group, MGIdx(0))
            stage2_mobile_fixed.addFixedAtoms(stage2_fixed_group)

        # forcefield holding the energy between fixed atoms
        stage1_mobile_mobile = InterCLJFF("mobile-mobile")
        stage1_mobile_mobile = setCLJProperties(stage1_mobile_mobile)
        stage1_mobile_mobile.add(stage1_mobile_group)

        stage2_mobile_mobile = InterCLJFF("mobile-mobile")
        stage2_mobile_mobile = setCLJProperties(stage2_mobile_mobile)
        stage2_mobile_mobile.add(stage2_mobile_group)

        stage1.add(stage1_mobile_group)
        stage1.add(stage1_mobile_fixed)
        stage1.add(stage1_mobile_mobile)

        stage2.add(stage2_mobile_group)
        stage2.add(stage2_mobile_fixed)
        stage2.add(stage2_mobile_mobile)

        stage1.setComponent(stage1.totalComponent(), stage1_mobile_mobile.components().total() + stage1_mobile_fixed.components().total())
        stage2.setComponent(stage2.totalComponent(), stage2_mobile_mobile.components().total() + stage2_mobile_fixed.components().total())

    stage1.add(stage1_mobile_group)
    stage2.add(stage2_mobile_group)

    stage1.add("volume_map", VolMapMonitor(stage1_mobile_group), 1000)
    stage2.add("volume_map", VolMapMonitor(stage2_mobile_group), 1000)

    max_water_translation = 0.15 * angstroms
    max_water_rotation = 15 * degrees

    stage1_moves = WeightedMoves()

    if stage1_mobile_group.nViews() > 0:
        rb_moves = RigidBodyMC(stage1_mobile_group)
        rb_moves.setMaximumTranslation(max_water_translation)
        rb_moves.setMaximumRotation(max_water_rotation)

        if stage1.containsProperty("reflection sphere radius"):
            reflection_radius = float(str(stage1.property("reflection sphere radius"))) * angstroms
            reflection_center = stage1.property("reflection center").toVector()[0]
            rb_moves.setReflectionSphere(reflection_center, reflection_radius)

        stage1_moves.add(rb_moves, stage1_mobile_group.nViews())

    stage2_moves = WeightedMoves()

    if stage2_mobile_group.nViews() > 0:
        rb_moves = RigidBodyMC(stage2_mobile_group)
        rb_moves.setMaximumTranslation(max_water_translation)
        rb_moves.setMaximumRotation(max_water_rotation)

        if stage2.containsProperty("reflection sphere radius"):
            reflection_radius = float(str(stage2.property("reflection sphere radius"))) * angstroms
            reflection_center = stage2.property("reflection center").toVector()[0]
            rb_moves.setReflectionSphere(reflection_center, reflection_radius)

        stage2_moves.add(rb_moves, stage2_mobile_group.nViews())

    stage1_moves.setTemperature(temperature.val)
    stage2_moves.setTemperature(temperature.val)

    seed = random_seed.val

    if seed is None:
        seed = RanGenerator().randInt(100000,1000000)
        print("Using generated random number seed %d" % seed)
    else:
        print("Using supplied random number seed %d" % seed)
    
    stage1_moves.setGenerator( RanGenerator(seed) )
    stage2_moves.setGenerator( RanGenerator(seed+7) )

    return ( (stage1,stage1_moves), (stage2,stage2_moves) )


def tryBackup(filename):
    # save the old file to a backup
    try:
        shutil.copy(filename, "%s.bak" % filename)
    except:
        pass


@resolveParameters
def run():

    print("Getting hold of (up to) 2 processor cores to run the calculation...")
    nodes = Nodes()
    this_thread = nodes.borrowThisThread()
    nodes.addNode()
    print(nodes)

    if os.path.exists(restart_file.val):
        (stage1, stage1_moves, stage2, stage2_moves) = Sire.Stream.load(restart_file.val)
        nblocks = int( stage1.property("nblocks").value() )
    else:
        (system, ligand_mol) = loadSystem(topfile.val, crdfile.val, s3file.val, ligand_name.val)

        watersys = loadWater()

        ( (stage1,stage1_moves), (stage2,stage2_moves) ) = makeSim(system, ligand_mol, watersys)

        if nequilmoves.val:
            print("Equilibrating for %d moves...." % nequilmoves.val)
            node1 = nodes.getNode()
            sim1 = Simulation.run(node1, stage1, stage1_moves, nequilmoves.val, False)
            node2 = nodes.getNode()
            sim2 = Simulation.run(node2, stage2, stage2_moves, nequilmoves.val, False)

            sim1.wait()
            stage1 = sim1.system()
            stage1_moves = sim1.moves()

            sim2.wait()
            stage2 = sim2.system()
            stage2_moves = sim2.moves()

            node1.release()
            node2.release() 

            sim1 = None
            sim2 = None

            stage1 = stage1_moves.move(stage1, nequilmoves.val, False)
            stage2 = stage2_moves.move(stage2, nequilmoves.val, False)

        nblocks = 0
        stage1.setProperty("nblocks", NumberProperty(nblocks))

        Sire.Stream.save( (stage1,stage1_moves,stage2,stage2_moves), restart_file.val )

    print("Number of blocks completed = %s" % nblocks)
    print("Number of blocks to perform = %s" % nmoves.val)

    for i in range(nblocks, nmoves.val):
        print("Performing block %s..." % (i+1))
        node1 = nodes.getNode()
        sim1 = Simulation.run(node1, stage1, stage1_moves, nsubmoves.val, True)
        node2 = nodes.getNode()
        sim2 = Simulation.run(node2, stage2, stage2_moves, nsubmoves.val, True)

        sim1.wait()
        stage1 = sim1.system()
        stage1_moves = sim1.moves()

        sim2.wait()
        stage2 = sim2.system()
        stage2_moves = sim2.moves()

        sim1 = None
        sim2 = None

        node1.release()
        node2.release()

        stage1.setProperty("nblocks", NumberProperty(i+1))

        if save_pdb.val:
            PDB().write(stage1.molecules(), "stage1_mobile_%0004d.pdb" % (i+1))
            PDB().write(stage2.molecules(), "stage2_mobile_%0004d.pdb" % (i+1))

        if i % restart_frequency.val == 0 or i == (nmoves.val - 1):
            tryBackup(restart_file.val)
            tryBackup("stage1_vol.dx")
            tryBackup("stage2_vol.dx")
                
            Sire.Stream.save( (stage1,stage1_moves,stage2,stage2_moves), restart_file.val )
            DX().write( stage1[MonitorName("volume_map")].volumeMap(), "stage1_vol.dx" )
            DX().write( stage2[MonitorName("volume_map")].volumeMap(), "stage2_vol.dx" )

    print("Simulation complete")
