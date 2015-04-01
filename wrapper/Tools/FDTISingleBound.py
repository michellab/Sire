#!/bin/env python
# -*- coding: utf-8 -*-

import os,re, sys, shutil

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

from Sire.Tools.AmberLoader import *
from Sire.Tools import Parameter, resolveParameters

import Sire.Stream

##### This is how we can have the script specify all of the 
##### user-controllable parameters

use_sphere = Parameter("use sphere", False,
                       """Whether or not to use sphereical boundary conditions""")

sphere_radius = Parameter("spherical boundary radius", 10*angstrom,
                          """The radius for the spherical boundary conditions.""")

sphere_center = None # this parameter will be calculated and set in the script

use_softcore = Parameter("use softcore", True,
                       """Whether or not to use a soft-core potential for the perturbed solute.""")

use_grid = Parameter("use grid", False,
                     """Whether or not to use a grid for the interactions with atoms 
                        that are beyond the spherical boundary""")

grid_spacing = Parameter("grid spacing", 0.5*angstrom,
                         """The spacing between grid points if a grid is used""")

grid_buffer = Parameter("grid buffer", 3*angstrom,
                        """The grid is generated to enclose all of the molecules in group 0,
                           plus a buffer specified by this parameter. The larger this buffer,
                           the larger the grid, but also the lower the chance that the grid
                           will need to be recalculated as the molecules in group 0 move.""")

cutoff_scheme = Parameter("cutoff scheme", "group",
                          """The method used to apply the non-bonded cutoff. Choices are;
                             (1) shift_electrostatics : This should become the default, and uses an atomistic cutoff
                                                        with a force-shifted cutoff.
                             (2) reaction_field : This uses the atomistic reaction field cutoff. You can
                                                  set the reaction field dielectric using the "dielectric"
                                                  parameter.
                             (3) group : This is the default, and uses a group-based cutoff with a feather. Note that this is 
                                         incompatible with a grid, so an error will be raised if you try
                                         to use a group-based cutoff with a grid.""")

rf_dielectric = Parameter("dielectric", 78.3,
                          """The dielectric constant to use with the reaction field cutoff method.""")

out_dir = Parameter("output directory", "output",
                    """The directory in which to place all output files.""")

top_file = Parameter("topology file", "../../SYSTEM.top",
                     """The name of the topology file that contains the solvated solute.""")

crd_file = Parameter("coordinate file", "../../SYSTEM.crd",
                     """The name of the coordinate file that contains the solvated solute.""")

ligand_flex_file = Parameter("ligand flex file", "../../MORPH.flex",
                             """Flexibility file describing how the morph is perturbed.""")

ligand_pert_file = Parameter("ligand perturbation file", "../../MORPH.pert",
                             """Perturbation file describing how the morph is perturbed.""")

protein_flex_file = Parameter("protein flex file", "../../PROTEIN.flex",
                              """Flexibility file describing which residues of the protein should be moved.""")

lig_name = Parameter("ligand name", "MORPH",
                     """Optional, the name of the ligand used in the flexibility file.
                        If the ligand has a single residue, the program will use the residue name
                        by default to look up the flexibility template""")

restart_file = Parameter("restart file", "sim_restart.s3",
                         """The name of the restart file.""")

random_seed = Parameter("random seed", 0, """The random number seed""")

nmoves = Parameter("number of moves", 50, """The number of moves per block""")

nmoves_per_energy = Parameter("number of energy snapshots", 1,
                              """The number of times during the simulation that you want the 
                                 energy to be recorded.""")

nmoves_per_pdb = Parameter("number of structure snapshots", 1,
                           """The number of times during the simulation that you want the 
                                    structure to be recorded (as a PDB).""")

nmoves_per_pdb_intermediates = Parameter("number of intermediate structure snapshots", None,
                                         """The number of times during an intermediate simulation to save 
                                            the structure (as a PDB).""")

temperature = Parameter("temperature", 25 * celsius, """The temperature of the simulation""")

pressure = Parameter("pressure", 1 * atm,
                     """The pressure of the simulation. Note that this is ignored if you
                        are using spherical boundary conditions.""")

coul_cutoff = Parameter("coulomb cutoff", 10*angstrom,
                        """The cutoff radius for non-bonded coulomb interactions""")

coul_feather = Parameter("coulomb feather", 0.5*angstrom,
                         """The feather radius for the non-bonded coulomb interactions
                            (only needed if a group-based cutoff is used)""")

lj_cutoff = Parameter("lj cutoff", 10*angstrom,
                      """The cutoff radius for non-bonded LJ interactions""")

lj_feather = Parameter("lj feather", 0.5*angstrom,
                       """The feather radius for the non-bonded LJ interactions
                          (only needed if a group-based cutoff is used)""")

coulomb_power = Parameter("coulomb power", 0,
                          """The soft-core coulomb power parameter""")

shift_delta = Parameter("shift delta", 2.0,
                        """The soft-core shift delta parameter""")

combining_rules = Parameter("combining rules", "arithmetic",
                            """The combinining rules for LJ interactions""")

pref_constant = Parameter("preferential constant", 200 * angstrom2,
                          """The preferential sampling constant""")

max_solvent_translation = Parameter("maximum solvent translation", 0.15*angstrom,
                                    """Maximum amount to translate the solvent""")

max_solvent_rotation = Parameter("maximum solvent rotation", 15*degrees,
                                 """Maximum amount to rotate the solvent""")
  
protein_mc_weight = Parameter("protein move weight", 1000,
                              """Factor used to multiply the weight of the protein moves.""")
                           
solvent_mc_weight_factor = Parameter("solvent move weight", 5,
                                     """Factor used to multiply the weight of the solvent moves.""")

solute_mc_weight = Parameter("solute move weight", 100,
                             """Factor used to multiply the weight of the solute moves.""")

volume_mc_weight = Parameter("volume move weight", 1,
                             """Factor used to multiply the weight of the volume moves.""")

delta_lambda = Parameter("delta lambda", 0.001,
                         """Delta lambda for finite difference gradients.""")

compress = Parameter("compression method", "bzip2 -f",
                     """Command used to compress output files.""")

lam_val = Parameter("lambda", 0.0, """Value of lambda for the simulation""")

print_nrgs = Parameter("print energies", None, 
                       """Whether or not to print all energy components after loading 
                          the restart file or starting the simulation. Useful for debugging.""")


####### FUNCTIONS  ###############

def readProteinFlexibility(protein_flex_file):

    FILE = open(protein_flex_file,'r')
    buffer = FILE.readlines()

    sc_flex = []
    bb_flex = []

    scmode = False
    bbmode = False

    for line in buffer:
        if line.startswith("#"):
            continue
        if line.startswith("flexible sidechain"):
            scmode = True
            continue
        if line.startswith("flexible backbone"):
            scmode = False
            bbmode = True
            continue

        elems = line.split()

        for elem in elems:
            res_number = ResNum( int(elem) )
            if scmode:
                sc_flex.append(res_number)
            elif bbmode:
                bb_flex.append(res_number)

    return sc_flex, bb_flex


def createBBMoveGroup(protein, bbgroup, flex_list):

    hn_atoms = AtomName("N", CaseInsensitive) * \
               AtomName("HN", CaseInsensitive) * AtomName("HN1", CaseInsensitive) * \
               AtomName("HN2", CaseInsensitive) * AtomName("HN3", CaseInsensitive)

    for i in range(0, protein.nResidues()):
        residue = protein.select(ResIdx(i))
        # Create only if mentioned as flexible
        if (not residue.number() in flex_list):
            continue
        atoms = protein.select(ResIdx(i)).selection()

        if i < (protein.nResidues()-1):
           try:
                atoms.deselect( hn_atoms + ResIdx(i) )
           except:
               pass

        if i > 0:
           try:
               atoms.select( hn_atoms + ResIdx(i-1) )
           except:
               pass

        bbgroup.add( PartialMolecule(protein, atoms) )


def adjustPerturbedDOFs( molecule ):
    
    perturbations = molecule.property("perturbations").perturbations()
    
    r0 = Symbol("r0")
    theta0 = Symbol("theta0")

    for pert in perturbations:
        if ( pert.what() == 'SireMM::TwoAtomPerturbation'):
            ri = pert.initialForms()[r0].toString().toDouble()
            rf = pert.finalForms()[r0].toString().toDouble()
            if (abs(ri-rf) > 0.001):
                #rint ri,rf
                r = (1-lam_val.val) * ri + lam_val.val * rf
                r = r * angstrom
                bond = BondID(pert.atom0(), pert.atom1() )
                mover = molecule.move()
                try:
                    mover.set(bond, r)
                except UserWarning:
                    # extract the type of the errror
                    _, error, _ = sys.exc_info()
                    error_type = re.search(r"(Sire\w*::\w*)", str(error)).group(0)
                    if error_type == "SireMol::ring_error":
                        continue
                molecule = mover.commit()
        elif ( pert.what() == 'SireMM::ThreeAtomPerturbation'):
            thetai = pert.initialForms()[theta0].toString().toDouble()
            thetaf = pert.finalForms()[theta0].toString().toDouble()
            if (abs(thetai-thetaf) > 0.001):
                #rint thetai,thetaf
                theta = (1-lam_val.val) * thetai + lam_val.val * thetaf
                theta = theta * radians
                angle = AngleID(pert.atom0(), pert.atom1(), pert.atom2() )
                mover = molecule.move()
                try:                
                    mover.set(angle, theta)
                except UserWarning:
                    # extract the type of the errror
                    _, err, _ = sys.exc_info()
                    error_type = re.search(r"(Sire\w*::\w*)", str(error)).group(0)
                    if error_type == "SireMol::ring_error":
                        continue
                molecule = mover.commit()
    return molecule


def getDummies(molecule):
    print "Selecting dummy groups"
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


def createSystem(molecules, space, naming_scheme=NamingScheme()):
    # First, sanity check that the cutoff is not greater than half the box length for
    # periodic spaces...
    radius = sphere_radius.val.to(angstrom)

    cutoff = coul_cutoff.val.to(angstrom)

    if lj_cutoff.val.to(angstrom) > cutoff:
        cutoff = lj_cutoff.val.to(angstrom)

    if space.isPeriodic():
        eps_cutoff = cutoff - 1e-6

        ok_x = (space.getMinimumImage(Vector(eps_cutoff,0,0), Vector(0)).length() <= cutoff)
        ok_y = (space.getMinimumImage(Vector(0,eps_cutoff,0), Vector(0)).length() <= cutoff)
        ok_z = (space.getMinimumImage(Vector(0,0,eps_cutoff), Vector(0)).length() <= cutoff)

        if not (ok_x and ok_y and ok_z):
            print >>sys.stderr,"The cutoff (%f A) is too large for periodic box %s" % \
                        (cutoff, space)
            raise RuntimeError()

        if use_grid.val:
            eps_radius = cutoff + radius - 1e-6

            ok_x = (space.getMinimumImage(Vector(eps_radius,0,0), Vector(0)).length() > radius)
            ok_y = (space.getMinimumImage(Vector(0,eps_radius,0), Vector(0)).length() > radius)
            ok_z = (space.getMinimumImage(Vector(0,0,eps_radius), Vector(0)).length() > radius)

            if not (ok_x and ok_y and ok_z):
                print >>sys.stderr,"The sphere radius (%f A) plus non-bonded cutoff (%f A) is too large for periodic box %s" \
                                 % (radius, cutoff, space)
                print >>sys.stderr, \
                         "Two times the sphere radius plus the cutoff distance cannot exceed the dimension of the box."

                raise RuntimeError()

    # add the ligand name to the naming scheme
    naming_scheme.addSoluteResidueName(lig_name.val)

    print "Applying flexibility and zmatrix templates..."

    moleculeNumbers = molecules.molNums()
    moleculeNumbers.sort()
    moleculeList = []

    for moleculeNumber in moleculeNumbers:
        molecule = molecules.molecule(moleculeNumber).molecule()
        moleculeList.append(molecule)

    # Create molecule groups and system

    system = System()

    all = MoleculeGroup("all")

    solutes = MoleculeGroup("solutes")
    perturbed_solutes = MoleculeGroup("perturbed_solutes")    

    protein = MoleculeGroup("protein")
    bbgroup = MoleculeGroup("bbresidues")
    residues = MoleculeGroup("residues")
    protein_and_buffer = MoleculeGroup("protein_and_buffer")

    solvent = MoleculeGroup("solvent")
    mobile_solvent = MoleculeGroup("solvent")
    water = MoleculeGroup( "water")
    ion = MoleculeGroup( "ion")

    mobilewater = MoleculeGroup("mobilewater")
    fixwater = MoleculeGroup("fixwater")

    mobileion = MoleculeGroup("mobileion")
    fixion = MoleculeGroup("fixion")

    fixed_atoms = MoleculeGroup("fixed_atoms")

    # We will need those to set properties
    zmat_maker = ZmatrixMaker()
    zmat_maker.loadTemplates(os.path.join(parameter_directory, "amber.zmatrices"))
    
    flexibility_lib = FlexibilityLibrary(ligand_flex_file.val)
    perturbations_lib = PerturbationsLibrary(ligand_pert_file.val)

    # Read the list of flexible sidechains and backbones
    sc_flex , bb_flex = readProteinFlexibility(protein_flex_file.val)

    # First, we need to find the solute molecule
    for molecule in moleculeList:
        if naming_scheme.isSolute(molecule):
            #
            # SOLUTE SETUP
            #
            if not solutes.isEmpty():
                print "This simulation script only supports one solute molecule."
                sys.exit(-1)

            # Solute...
            solute = molecule
            # If ligname has not been defined, and there is a single residue,
            # use the residue name
            ligand_name = lig_name.val

            if ligand_name is None:
                if ( solute.nResidues() == 1 ):
                    ligand_name = solute.residue( ResIdx(0) ).name().value()
                else:
                    ligand_name = "ligand" # Likely not good...

            #print lig_name
            solute = solute.edit().rename(ligand_name).commit()
            # print solute
            # This will add the property "flexibility" to the solute

            flexibility = flexibility_lib.getFlexibility(solute)
            solute = solute.edit().setProperty("flexibility", flexibility).commit()

            solute = perturbations_lib.applyTemplate(solute)

            perturbations = solute.property("perturbations")

            # print lam_val
            lam = Symbol("lambda")
            lam_fwd = Symbol("lambda_{fwd}")
            lam_bwd = Symbol("lambda_{bwd}")
    
            initial = Perturbation.symbols().initial()
            final = Perturbation.symbols().final()

            solute = solute.edit().setProperty("perturbations",
                       perturbations.recreate( (1-lam)*initial + lam*final ) ).commit()

            # Set the geometry of perturbed bonds/angles to match the corresponding equilibrium value 
            solute = adjustPerturbedDOFs( solute )

            solute_fwd = solute.edit().renumber().setProperty("perturbations",
                            perturbations.substitute( lam, lam_fwd ) ).commit()
            solute_bwd = solute.edit().renumber().setProperty("perturbations",
                            perturbations.substitute( lam, lam_bwd ) ).commit()

            # print solute
            # print solute_fwd
            # print solute_bwd

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

            solute_grp_fwd = MoleculeGroup("solute_fwd", solute_fwd)
            solute_grp_fwd_hard = MoleculeGroup("solute_fwd_hard")
            solute_grp_fwd_todummy = MoleculeGroup("solute_fwd_todummy")
            solute_grp_fwd_fromdummy = MoleculeGroup("solute_fwd_fromdummy")

            solute_grp_bwd = MoleculeGroup("solute_bwd", solute_bwd)
            solute_grp_bwd_hard = MoleculeGroup("solute_bwd_hard")
            solute_grp_bwd_todummy = MoleculeGroup("solute_bwd_todummy")
            solute_grp_bwd_fromdummy = MoleculeGroup("solute_bwd_fromdummy")
    
            solute_ref_hard = solute.selectAllAtoms()
            solute_ref_todummy = solute_ref_hard.invert()
            solute_ref_fromdummy = solute_ref_hard.invert()
            solute_fwd_hard = solute_fwd.selectAllAtoms()
            solute_fwd_todummy = solute_fwd_hard.invert()
            solute_fwd_fromdummy = solute_fwd_hard.invert()
            solute_bwd_hard = solute_bwd.selectAllAtoms()
            solute_bwd_todummy = solute_bwd_hard.invert()
            solute_bwd_fromdummy = solute_bwd_hard.invert()

            to_dummies, from_dummies = getDummies(solute)

            # print to_dummies
            # print from_dummies

            if to_dummies is not None:
                ndummies = to_dummies.count()
                dummies = to_dummies.atoms()

                for x in range(0,ndummies):
                    dummy_index = dummies[x].index()
                    solute_ref_hard = solute_ref_hard.subtract( solute.select( dummy_index ) )
                    solute_fwd_hard = solute_fwd_hard.subtract( solute_fwd.select( dummy_index ) )
                    solute_bwd_hard = solute_bwd_hard.subtract( solute_bwd.select( dummy_index ) )
                    solute_ref_todummy = solute_ref_todummy.add( solute.select( dummy_index ) )
                    solute_fwd_todummy = solute_fwd_todummy.add( solute_fwd.select( dummy_index ) )
                    solute_bwd_todummy = solute_bwd_todummy.add( solute_bwd.select( dummy_index ) )

            if from_dummies is not None:
                ndummies = from_dummies.count()
                dummies = from_dummies.atoms()

                for x in range(0,ndummies):
                    dummy_index = dummies[x].index()
                    solute_ref_hard = solute_ref_hard.subtract( solute.select( dummy_index ) )
                    solute_fwd_hard = solute_fwd_hard.subtract( solute_fwd.select( dummy_index ) )
                    solute_bwd_hard = solute_bwd_hard.subtract( solute_bwd.select( dummy_index ) )
                    solute_ref_fromdummy = solute_ref_fromdummy.add( solute.select( dummy_index ) )
                    solute_fwd_fromdummy = solute_fwd_fromdummy.add( solute_fwd.select( dummy_index ) )
                    solute_bwd_fromdummy = solute_bwd_fromdummy.add( solute_bwd.select( dummy_index ) )

            solute_grp_ref_hard.add(solute_ref_hard)
            solute_grp_fwd_hard.add(solute_fwd_hard)
            solute_grp_bwd_hard.add(solute_bwd_hard)
            
            solute_grp_ref_todummy.add(solute_ref_todummy)
            solute_grp_fwd_todummy.add(solute_fwd_todummy)
            solute_grp_bwd_todummy.add(solute_bwd_todummy)

            solute_grp_ref_fromdummy.add(solute_ref_fromdummy)
            solute_grp_fwd_fromdummy.add(solute_fwd_fromdummy)
            solute_grp_bwd_fromdummy.add(solute_bwd_fromdummy)
 
            solutes.add(solute)
            solutes.add(solute_fwd)
            solutes.add(solute_bwd)

            perturbed_solutes.add(solute_fwd)
            perturbed_solutes.add(solute_bwd)

    ligand = solutes[MolIdx(0)].molecule()
    
    global sphere_center
    sphere_center = ligand.evaluate().center()

    # make sure that this center is used by the ligand from now on
    ligand = ligand.edit().setProperty("center", wrap(sphere_center)).commit()
    solutes.update(ligand)
    molecules.update(ligand)

    num_images = 0

    for molecule in moleculeList:
        # get the center of the solvent
        mol_center = molecule.evaluate().center()

        if not naming_scheme.isSolute(molecule):
            # first, if this is a periodic space, then wrap the molecule
            # into the same periodic box as the ligand
            if space.isPeriodic():
                # wrap the molecule into the same space as the solute
                wrapped_mol_center = space.getMinimumImage(mol_center, sphere_center)

                if wrapped_mol_center != mol_center:
                    molecule = molecule.move().translate( wrapped_mol_center - mol_center ).commit()
                    mol_center = molecule.evaluate().center()

                # next, if we are using a grid, then we need to manually add the images of this
                # molecule from the neighbouring periodic boxes

                if use_grid.val:
                    image_cutoff = cutoff + radius

                    if Vector.distance(mol_center, sphere_center) > radius:
                        # this will be one of the fixed molecules, so may will need to be manually wrapped
                        for i in (-1,0,1):
                            for j in (-1,0,1):
                                for k in (-1,0,1):
                                    delta = Vector(i * radius, j * radius, k * radius)

                                    if delta.length() != 0:
                                        # get the image of the solvent molecule in this box
                                        image_solv_center = space.getMinimumImage(solv_center, sphere_center+delta)

                                        delta = image_solv_center - solv_center

                                        if delta.length() > 0:
                                            #there is a periodic image available here - is it within non-bonded cutoff?
                                            if (image_solv_center - sphere_center).length() <= image_cutoff:
                                               # it is within cutoff, so a copy of this molecule should be added
                                               image = molecule.edit().renumber().move().translate(delta).commit()
                                               fixed_atoms.add(image)
                                               num_images += 1

        # Ok, we have wrapped the molecule into the same box as the solute and have added any necessary images.
        # Now lets process each molecule according to its type...
        if naming_scheme.isProtein(molecule):
            # This will add the property "z-matrix" to the protein

            molecule = zmat_maker.applyTemplates( molecule )

            if use_grid.val:
                # loop through all residues and find all those with at least one
                # atom in the reflection sphere
                space = Cartesian()
                sphere_resnums = []

                # the extra atoms moved as part of a backbone move
                hn_atoms = AtomName("N", CaseInsensitive) * AtomName("H", CaseInsensitive) * \
                           AtomName("HN", CaseInsensitive) * AtomName("HN1", CaseInsensitive) * \
                           AtomName("HN2", CaseInsensitive) * AtomName("HN3", CaseInsensitive)

                for i in range(0, molecule.nResidues()):
                    res = molecule.residue( ResIdx(i) )
                    distance = space.minimumDistance(CoordGroup(1,sphere_center), getCoordGroup(res.atoms()))
          
                    if distance < radius:
                        sphere_resnums.append( res.number() )
                        protein_and_buffer.add(res)
                        protein.add(res)

                        if ( res.number() in sc_flex ):
                            # add the residue to the mobile sidechains group
                            residues.add(res)

                        if ( res.number() in bb_flex ):
                            # now add the atoms needed from the residue to the mobile backbones group
                            atoms = molecule.select(ResIdx(i)).selection()
    
                            if i < (molecule.nResidues()-1):
                                try:
                                    atoms.deselect( hn_atoms + ResIdx(i) )
                                except:
                                    pass

                            if i > 0:
                                try:
                                    atoms.select( hn_atoms + ResIdx(i+1) )
                                except:
                                    pass

                            bbgroup.add( PartialMolecule(molecule, atoms) )

                # now loop over all of the residues and work out which ones are fixed, and which ones
                # are bonded to sphere residues
                connectivity = molecule.property("connectivity")

                for i in range(0, molecule.nResidues()):
                    res = molecule.residue( ResIdx(i) )

                    if not res.number() in sphere_resnums:
                        # is this residue bonded to any of the sphere residues? If so, then it is a boundary residue
                        is_boundary = False

                        for bonded_res in connectivity.connectionsTo( res.number() ):
                            bonded_resnum = molecule.residue(bonded_res).number()

                            if bonded_resnum in sphere_resnums:
                                is_boundary = True
                                break

                        if is_boundary:
                            protein_and_buffer.add(res)
                        else:
                            fixed_atoms.add(res)

            else:
                # Update the MoleculeGroup used to perform backbone moves
                # on the protein
                createBBMoveGroup(molecule, bbgroup, bb_flex)

                for i in range(0, molecule.nResidues()):
                    res = molecule.residue( ResIdx(i) )
                    # Skip residues that have not been flagged as flexible
                    if ( not res.number() in sc_flex ):
                        continue
                    residues.add(res)

                protein.add(molecule)

        elif naming_scheme.isWater(molecule):
            # Separate water from ions
            water.add(molecule)
            solvent.add(molecule)

            if Vector.distance(mol_center, sphere_center) < radius:
                mobilewater.add(molecule)
                mobile_solvent.add(molecule)
            else:
                fixwater.add(molecule)
                fixed_atoms.add(molecule)

        elif naming_scheme.isIon(molecule):
            ion.add(molecule)
            solvent.add(molecule)

            if Vector.distance(mol_center, sphere_center) < radius:
                mobileions.add(molecule)
                mobile_solvent.add(molecule)
            else:
                fixion.add(molecule)
                fixed_atoms.add(molecule)
        
        elif naming_scheme.isSolute(molecule):
            # already processed the solute
            pass

        else:
            print "Oops, script does not know in what kind of molecule this residue is %s " % (molecule.residues()[0].name().value())
            raise RunTimeError()


    if use_grid.val:
        # now the fun part - we need to split the protein into two parts so that the fixed
        # atoms can be removed from the main code. This significantly improves the efficiency
        # of the code

        new_protein = MoleculeGroup("protein")
        new_bbgroup = MoleculeGroup("bbresidues")
        new_residues = MoleculeGroup("residues")
        new_protein_and_buffer = MoleculeGroup("protein_and_buffer")

        for molnum in protein_and_buffer.molNums():
            protein_mol = protein_and_buffer[molnum].join()

            if protein_mol.selectedAll():
                new_protein_and_buffer.add( protein_mol )

                # the protein hasn't changed, so just copy across the existing protein into the new groups
                try:
                    new_protein.add( protein[molnum] )
                except:
                    pass

                try:
                    new_residues.add( residues[molnum] )
                except:
                    pass

                try:
                    new_bbgroup.add( bbgroup[molnum] )
                except:
                    pass

            else:
                # some of the protein is fixed, so we need to extract only the part
                # that is in the reflection sphere. Also renumber the molecule, to prevent
                # clashes with the original
                new_protein_mol = protein_mol.extract().molecule().edit().renumber().commit()

                new_protein_and_buffer.add(new_protein_mol)

                # copy the selection from the old groups using the new_protein extracted molecule
                try:
                    protein_views = protein[molnum]

                    for i in range(0,protein_views.nViews()):
                        view = new_protein_mol.selection()
                        view = view.selectNone()

                        for atomid in protein_views.viewAt(i).selectedAtoms():
                            atom = protein_mol.atom(atomid)
                            resatomid = ResAtomID( atom.residue().number(), atom.name() )
                            view = view.select( resatomid )

                        if view.nSelected() > 0:
                            new_protein.add( PartialMolecule(new_protein_mol, view) )
                except:
                    print "ERROR IN PROTEIN"
                    _, error, _ = sys.exc_info()
                    print error
                    pass


                try:
                    residue_views = residues[molnum]

                    for i in range(0,residue_views.nViews()):
                        view = new_protein_mol.selection()
                        view = view.selectNone()

                        for atomid in residue_views.viewAt(i).selectedAtoms():
                            atom = protein_mol.atom(atomid)
                            resatomid = ResAtomID( atom.residue().number(), atom.name() )
                            view = view.select( resatomid )

                        if view.nSelected() > 0:
                            new_residues.add( PartialMolecule(new_protein_mol, view) )
                except:
                    print "ERROR IN RESIDUES"
                    _, error, _ = sys.exc_info()
                    print error
                    pass

                try:
                    bb_views = bbgroup[molnum]

                    for i in range(0,bb_views.nViews()):
                        view = new_protein_mol.selection()
                        view = view.selectNone()

                        for atomid in bb_views.viewAt(i).selectedAtoms():
                            atom = protein_mol.atom(atomid)
                            resatomid = ResAtomID( atom.residue().number(), atom.name() )
                            view = view.select( resatomid )

                        if view.nSelected() > 0:
                            new_bbgroup.add( PartialMolecule(new_protein_mol, view) )
                except:
                    print "ERROR IN BBGROUP"
                    _, error, _ = sys.exc_info()
                    print error
                    pass

        protein = new_protein
        bbgroup = new_bbgroup
        residues = new_residues
        protein_and_buffer = new_protein_and_buffer
        print "Number of images == %d" % num_images

    all.add(solutes)
    
    all.add(perturbed_solutes)

    all.add(solute_grp_ref)
    all.add(solute_grp_fwd)
    all.add(solute_grp_bwd)

    all.add(solute_grp_ref_hard)
    all.add(solute_grp_ref_todummy)
    all.add(solute_grp_ref_fromdummy)

    all.add(solute_grp_fwd_hard)
    all.add(solute_grp_fwd_todummy)
    all.add(solute_grp_fwd_fromdummy)

    all.add(solute_grp_bwd_hard)
    all.add(solute_grp_bwd_todummy)
    all.add(solute_grp_bwd_fromdummy)

    traj = MoleculeGroup("traj")

    traj.add(solute)

    if use_grid.val:
        # don't save the fixed molecules as these will take up space
        # and slow down the simulation
        all.add(mobilewater)
        all.add(mobileion)
        all.add(protein_and_buffer)

        traj.add(mobilewater)
        traj.add(mobileion)
        traj.add(protein_and_buffer)

        system.add(mobile_solvent)
        system.add(mobilewater)
        system.add(mobileion)
        system.add(protein_and_buffer)
        system.add(protein)
    else:
        all.add(water)
        all.add(ion) 
        all.add(protein)

        traj.add(water)
        traj.add(ion)
        traj.add(protein)

        system.add(solvent)
        system.add(water)
        system.add(ion)

        system.add(mobilewater)
        system.add(fixwater)

        system.add(mobileion)
        system.add(fixion)
    
        system.add(protein)

    # NOT SURE NEED TO ADD REDUNDANT GRPS TO ALL...

    # Add these groups to the System
    system.add(solutes)

    system.add(perturbed_solutes)

    system.add(solute_grp_ref)
    system.add(solute_grp_fwd)
    system.add(solute_grp_bwd)

    system.add(solute_grp_ref_hard)
    system.add(solute_grp_ref_todummy)
    system.add(solute_grp_ref_fromdummy)

    system.add(solute_grp_fwd_hard)
    system.add(solute_grp_fwd_todummy)
    system.add(solute_grp_fwd_fromdummy)

    system.add(solute_grp_bwd_hard)
    system.add(solute_grp_bwd_todummy)
    system.add(solute_grp_bwd_fromdummy)

    system.add(residues)
    system.add(bbgroup)

    system.add(all)

    system.add(traj)

    if use_grid.val:
        # add the fixed atoms as a system property, so that they are not included
        # in the main simulated system
        system.setProperty("fixed_atoms", fixed_atoms)

    return system


def setupForcefields(system, space):

    print "Creating force fields... "

    solutes = system[ MGName("solutes") ]

    solute = system[ MGName("solute_ref") ]
    solute_hard = system[ MGName("solute_ref_hard") ]
    solute_todummy = system[ MGName("solute_ref_todummy") ]
    solute_fromdummy = system[ MGName("solute_ref_fromdummy") ]

    solute_fwd = system[ MGName("solute_fwd") ]
    solute_fwd_hard = system[ MGName("solute_fwd_hard") ]
    solute_fwd_todummy = system[ MGName("solute_fwd_todummy") ]
    solute_fwd_fromdummy = system[ MGName("solute_fwd_fromdummy") ]

    solute_bwd = system[ MGName("solute_bwd") ]
    solute_bwd_hard = system[ MGName("solute_bwd_hard") ]
    solute_bwd_todummy = system[ MGName("solute_bwd_todummy") ]
    solute_bwd_fromdummy = system[ MGName("solute_bwd_fromdummy") ]

    solvent = system[ MGName("solvent") ]

    protein = system[ MGName("protein") ]

    # As we build forcefields, keep a list of each type so we can 
    # easily set parameters etc.
    cljffs = []
    gridffs = []
    intraffs = []

    if use_grid.val:
        protein = system[ MGName("protein_and_buffer") ]
        fixed_atoms = system.property("fixed_atoms")
            
        # Start by creating a template GridFF forcefield, which can be duplicated
        # for each grid. This ensures that only a single copy of the fixed atoms
        # will be saved in the system, saving space and improving efficiency
        gridff = GridFF("template")
        gridff.addFixedAtoms(fixed_atoms)
        gridff.setGridSpacing( grid_spacing.val )
        gridff.setBuffer( grid_buffer.val )
        gridff.setCoulombCutoff( coul_cutoff.val )
        gridff.setLJCutoff( lj_cutoff.val )
        gridff.setProperty("combiningRules", VariantProperty(combining_rules.val) )
        gridff.setProperty("space", Cartesian())
           
        if cutoff_scheme.val == "shift_electrostatics":
            gridff.setShiftElectrostatics(True)

        elif cutoff_scheme.val == "reaction_field":
            gridff.setUseReactionField(True)
            gridff.setReactionFieldDielectric(rf_dielectric.val)

        elif cutoff_scheme.val == "group":
            print >>sys.stderr,"You cannot use a group-based cutoff with a grid!"
            print >>sys.stderr,"Please choose either the shift_electrostatics or reaction_field cutoff schemes."
            raise RuntimeError()

        else:
            print "WARNING. Unrecognised cutoff scheme. Using \"shift_electrostatics\"."
            gridff.setShiftElectrostatics(True)

    residues = system[ MGName("residues") ]

    # - first solvent-solvent coulomb/LJ (CLJ) energy
    solventff = InterCLJFF("solvent:solvent")
    solventff.add(solvent)
    cljffs.append(solventff)

    # The protein bond, angle, dihedral energy
    protein_intraff = InternalFF("protein_intraff")
    protein_intraff.add(protein)
    intraffs.append(protein_intraff)

    # The protein intramolecular CLJ energy
    protein_intraclj = IntraCLJFF("protein_intraclj")
    protein_intraclj.add(protein)
    cljffs.append(protein_intraclj)

    # The protein-solvent energy 
    protein_solventff = InterGroupCLJFF("protein:solvent")
    protein_solventff.add(protein, MGIdx(0))
    protein_solventff.add(solvent, MGIdx(1))
    cljffs.append(protein_solventff)

    # now the energy between the protein+solvent and the fixed atoms
    if use_grid.val:
        protein_solvent_fixedff = gridff.clone()
        protein_solvent_fixedff.setName("protein_solvent:fixed")
        protein_solvent_fixedff.add(protein, MGIdx(0))
        protein_solvent_fixedff.add(solvent, MGIdx(0))
        gridffs.append(protein_solvent_fixedff)

    # Now solute bond, angle, dihedral energy
    solute_intraff = InternalFF("solute_intraff")
    solute_intraff.add(solute)
    intraffs.append(solute_intraff)    

    solute_fwd_intraff = InternalFF("solute_fwd_intraff")
    solute_fwd_intraff.add(solute_fwd)
    intraffs.append(solute_fwd_intraff)

    solute_bwd_intraff = InternalFF("solute_bwd_intraff")
    solute_bwd_intraff.add(solute_bwd)
    intraffs.append(solute_bwd_intraff)

    # Now solute intramolecular CLJ energy
    solute_hard_intraclj = IntraCLJFF("solute_hard_intraclj")
    solute_hard_intraclj.add(solute_hard)
    cljffs.append(solute_hard_intraclj)

    solute_todummy_intraclj = IntraSoftCLJFF("solute_todummy_intraclj")
    solute_todummy_intraclj.add(solute_todummy)
    cljffs.append(solute_todummy_intraclj)

    solute_fromdummy_intraclj = IntraSoftCLJFF("solute_fromdummy_intraclj")
    solute_fromdummy_intraclj.add(solute_fromdummy)
    cljffs.append(solute_fromdummy_intraclj)

    solute_hard_todummy_intraclj = IntraGroupSoftCLJFF("solute_hard:todummy_intraclj")
    solute_hard_todummy_intraclj.add(solute_hard, MGIdx(0))
    solute_hard_todummy_intraclj.add(solute_todummy, MGIdx(1))
    cljffs.append(solute_hard_todummy_intraclj)

    solute_hard_fromdummy_intraclj = IntraGroupSoftCLJFF("solute_hard:fromdummy_intraclj")
    solute_hard_fromdummy_intraclj.add(solute_hard, MGIdx(0))
    solute_hard_fromdummy_intraclj.add(solute_fromdummy, MGIdx(1))
    cljffs.append(solute_hard_fromdummy_intraclj)

    solute_todummy_fromdummy_intraclj = IntraGroupSoftCLJFF("solute_todummy:fromdummy_intraclj")
    solute_todummy_fromdummy_intraclj.add(solute_todummy, MGIdx(0))
    solute_todummy_fromdummy_intraclj.add(solute_fromdummy, MGIdx(1))
    cljffs.append(solute_todummy_fromdummy_intraclj)

    # The forwards intramolecular CLJ energy

    solute_fwd_hard_intraclj = IntraCLJFF("solute_fwd_hard_intraclj")
    solute_fwd_hard_intraclj.add(solute_fwd_hard)
    cljffs.append(solute_fwd_hard_intraclj)

    solute_fwd_todummy_intraclj = IntraSoftCLJFF("solute_fwd_todummy_intraclj")
    solute_fwd_todummy_intraclj.add(solute_fwd_todummy)
    cljffs.append(solute_fwd_todummy_intraclj)

    solute_fwd_fromdummy_intraclj = IntraSoftCLJFF("solute_fwd_fromdummy_intraclj")
    solute_fwd_fromdummy_intraclj.add(solute_fwd_fromdummy)
    cljffs.append(solute_fwd_fromdummy_intraclj)

    solute_fwd_hard_todummy_intraclj = IntraGroupSoftCLJFF("solute_fwd_hard:todummy_intraclj")
    solute_fwd_hard_todummy_intraclj.add(solute_fwd_hard, MGIdx(0))
    solute_fwd_hard_todummy_intraclj.add(solute_fwd_todummy, MGIdx(1))
    cljffs.append(solute_fwd_hard_todummy_intraclj)

    solute_fwd_hard_fromdummy_intraclj = IntraGroupSoftCLJFF("solute_fwd_hard:fromdummy_intraclj")
    solute_fwd_hard_fromdummy_intraclj.add(solute_fwd_hard, MGIdx(0))
    solute_fwd_hard_fromdummy_intraclj.add(solute_fwd_fromdummy, MGIdx(1))
    cljffs.append(solute_fwd_hard_fromdummy_intraclj)

    solute_fwd_todummy_fromdummy_intraclj = IntraGroupSoftCLJFF("solute_fwd_todummy:fromdummy_intraclj")
    solute_fwd_todummy_fromdummy_intraclj.add(solute_fwd_todummy, MGIdx(0))
    solute_fwd_todummy_fromdummy_intraclj.add(solute_fwd_fromdummy, MGIdx(1))
    cljffs.append(solute_fwd_todummy_fromdummy_intraclj)

    # The backwards intramolecular CLJ energy

    solute_bwd_hard_intraclj = IntraCLJFF("solute_bwd_hard_intraclj")
    solute_bwd_hard_intraclj.add(solute_bwd_hard)
    cljffs.append(solute_bwd_hard_intraclj)

    solute_bwd_todummy_intraclj = IntraSoftCLJFF("solute_bwd_todummy_intraclj")
    solute_bwd_todummy_intraclj.add(solute_bwd_todummy)
    cljffs.append(solute_bwd_todummy_intraclj)

    solute_bwd_fromdummy_intraclj = IntraSoftCLJFF("solute_bwd_fromdummy_intraclj")
    solute_bwd_fromdummy_intraclj.add(solute_bwd_fromdummy)
    cljffs.append(solute_bwd_fromdummy_intraclj)

    solute_bwd_hard_todummy_intraclj = IntraGroupSoftCLJFF("solute_bwd_hard:todummy_intraclj")
    solute_bwd_hard_todummy_intraclj.add(solute_bwd_hard, MGIdx(0))
    solute_bwd_hard_todummy_intraclj.add(solute_bwd_todummy, MGIdx(1))
    cljffs.append(solute_bwd_hard_todummy_intraclj)

    solute_bwd_hard_fromdummy_intraclj = IntraGroupSoftCLJFF("solute_bwd_hard:fromdummy_intraclj")
    solute_bwd_hard_fromdummy_intraclj.add(solute_bwd_hard, MGIdx(0))
    solute_bwd_hard_fromdummy_intraclj.add(solute_bwd_fromdummy, MGIdx(1))
    cljffs.append(solute_bwd_hard_fromdummy_intraclj)

    solute_bwd_todummy_fromdummy_intraclj = IntraGroupSoftCLJFF("solute_bwd_todummy:fromdummy_intraclj")
    solute_bwd_todummy_fromdummy_intraclj.add(solute_bwd_todummy, MGIdx(0))
    solute_bwd_todummy_fromdummy_intraclj.add(solute_bwd_fromdummy, MGIdx(1))
    cljffs.append(solute_bwd_todummy_fromdummy_intraclj)

    # Now the solute-solvent CLJ energy
    solute_hard_solventff = InterGroupCLJFF("solute_hard:solvent")
    solute_hard_solventff.add(solute_hard, MGIdx(0))
    solute_hard_solventff.add(solvent, MGIdx(1))
    cljffs.append(solute_hard_solventff)

    solute_todummy_solventff = InterGroupSoftCLJFF("solute_todummy:solvent")
    solute_todummy_solventff.add(solute_todummy, MGIdx(0))
    solute_todummy_solventff.add(solvent, MGIdx(1))
    cljffs.append(solute_todummy_solventff)

    solute_fromdummy_solventff = InterGroupSoftCLJFF("solute_fromdummy:solvent")
    solute_fromdummy_solventff.add(solute_fromdummy, MGIdx(0))
    solute_fromdummy_solventff.add(solvent, MGIdx(1))
    cljffs.append(solute_fromdummy_solventff)

    if use_grid.val:
        solute_fixedff = gridff.clone()
        solute_fixedff.setName("solute:fixed")
        solute_fixedff.add(solute_hard, MGIdx(0))
        solute_fixedff.add(solute_todummy, MGIdx(0))
        solute_fixedff.add(solute_fromdummy, MGIdx(0))
        gridffs.append(solute_fixedff)

    # Now the forwards solute-solvent CLJ energy
    solute_fwd_hard_solventff = InterGroupCLJFF("solute_fwd_hard:solvent")
    solute_fwd_hard_solventff.add(solute_fwd_hard, MGIdx(0))
    solute_fwd_hard_solventff.add(solvent, MGIdx(1))
    cljffs.append(solute_fwd_hard_solventff)

    solute_fwd_todummy_solventff = InterGroupSoftCLJFF("solute_fwd_todummy:solvent")
    solute_fwd_todummy_solventff.add(solute_fwd_todummy, MGIdx(0))
    solute_fwd_todummy_solventff.add(solvent, MGIdx(1))
    cljffs.append(solute_fwd_todummy_solventff)

    solute_fwd_fromdummy_solventff = InterGroupSoftCLJFF("solute_fwd_fromdummy:solvent")
    solute_fwd_fromdummy_solventff.add(solute_fwd_fromdummy, MGIdx(0))
    solute_fwd_fromdummy_solventff.add(solvent, MGIdx(1))
    cljffs.append(solute_fwd_fromdummy_solventff)

    if use_grid.val:
        solute_fwd_fixedff = gridff.clone()
        solute_fwd_fixedff.setName("solute_fwd:fixed")
        solute_fwd_fixedff.add(solute_fwd_hard, MGIdx(0))
        solute_fwd_fixedff.add(solute_fwd_todummy, MGIdx(0))
        solute_fwd_fixedff.add(solute_fwd_fromdummy, MGIdx(0))
        gridffs.append(solute_fwd_fixedff)

    # Now the backwards solute-solvent CLJ energy
    solute_bwd_hard_solventff = InterGroupCLJFF("solute_bwd_hard:solvent")
    solute_bwd_hard_solventff.add(solute_bwd_hard, MGIdx(0))
    solute_bwd_hard_solventff.add(solvent, MGIdx(1))
    cljffs.append(solute_bwd_hard_solventff)

    solute_bwd_todummy_solventff = InterGroupSoftCLJFF("solute_bwd_todummy:solvent")
    solute_bwd_todummy_solventff.add(solute_bwd_todummy, MGIdx(0))
    solute_bwd_todummy_solventff.add(solvent, MGIdx(1))
    cljffs.append(solute_bwd_todummy_solventff)

    solute_bwd_fromdummy_solventff = InterGroupSoftCLJFF("solute_bwd_fromdummy:solvent")
    solute_bwd_fromdummy_solventff.add(solute_bwd_fromdummy, MGIdx(0))
    solute_bwd_fromdummy_solventff.add(solvent, MGIdx(1))
    cljffs.append(solute_bwd_fromdummy_solventff)

    if use_grid.val:
        solute_bwd_fixedff = gridff.clone()
        solute_bwd_fixedff.setName("solute_bwd:fixed")
        solute_bwd_fixedff.add(solute_bwd_hard, MGIdx(0))
        solute_bwd_fixedff.add(solute_bwd_todummy, MGIdx(0))
        solute_bwd_fixedff.add(solute_bwd_fromdummy, MGIdx(0))
        gridffs.append(solute_bwd_fixedff)

    # Now the solute-protein ( soft-core ) CLJ energy
    solute_hard_proteinff = InterGroupCLJFF("solute_hard:protein")
    solute_hard_proteinff.add(solute_hard, MGIdx(0))
    solute_hard_proteinff.add(protein, MGIdx(1))
    cljffs.append(solute_hard_proteinff)

    solute_todummy_proteinff = InterGroupSoftCLJFF("solute_todummy:protein")
    solute_todummy_proteinff.add(solute_todummy, MGIdx(0))
    solute_todummy_proteinff.add(protein, MGIdx(1))
    cljffs.append(solute_todummy_proteinff)

    solute_fromdummy_proteinff = InterGroupSoftCLJFF("solute_fromdummy:protein")
    solute_fromdummy_proteinff.add(solute_fromdummy, MGIdx(0))
    solute_fromdummy_proteinff.add(protein, MGIdx(1))
    cljffs.append(solute_fromdummy_proteinff)

    # Now the forwards solute-protein CLJ energy
    solute_fwd_hard_proteinff = InterGroupCLJFF("solute_fwd_hard:protein")
    solute_fwd_hard_proteinff.add(solute_fwd_hard, MGIdx(0))
    solute_fwd_hard_proteinff.add(protein, MGIdx(1))
    cljffs.append(solute_fwd_hard_proteinff)

    solute_fwd_todummy_proteinff = InterGroupSoftCLJFF("solute_fwd_todummy:protein")
    solute_fwd_todummy_proteinff.add(solute_fwd_todummy, MGIdx(0))
    solute_fwd_todummy_proteinff.add(protein, MGIdx(1))
    cljffs.append(solute_fwd_todummy_proteinff)

    solute_fwd_fromdummy_proteinff = InterGroupSoftCLJFF("solute_fwd_fromdummy:protein")
    solute_fwd_fromdummy_proteinff.add(solute_fwd_fromdummy, MGIdx(0))
    solute_fwd_fromdummy_proteinff.add(protein, MGIdx(1))
    cljffs.append(solute_fwd_fromdummy_proteinff)

    # Now the backwards solute-protein CLJ energy
    solute_bwd_hard_proteinff = InterGroupCLJFF("solute_bwd_hard:protein")
    solute_bwd_hard_proteinff.add(solute_bwd_hard, MGIdx(0))
    solute_bwd_hard_proteinff.add(protein, MGIdx(1))
    cljffs.append(solute_bwd_hard_proteinff)

    solute_bwd_todummy_proteinff = InterGroupSoftCLJFF("solute_bwd_todummy:protein")
    solute_bwd_todummy_proteinff.add(solute_bwd_todummy, MGIdx(0))
    solute_bwd_todummy_proteinff.add(protein, MGIdx(1))
    cljffs.append(solute_bwd_todummy_proteinff)

    solute_bwd_fromdummy_proteinff = InterGroupSoftCLJFF("solute_bwd_fromdummy:protein")
    solute_bwd_fromdummy_proteinff.add(solute_bwd_fromdummy, MGIdx(0))
    solute_bwd_fromdummy_proteinff.add(protein, MGIdx(1))
    cljffs.append(solute_bwd_fromdummy_proteinff)

    for cljff in cljffs:
        try:
            if use_grid.val:
                cljff.setProperty("space", Cartesian())
            else:
                cljff.setProperty("space", space)

            cljff.setProperty("switchingFunction", HarmonicSwitchingFunction(coul_cutoff.val, coul_feather.val,
                                                                             lj_cutoff.val, lj_feather.val) )
            cljff.setProperty("combiningRules", VariantProperty(combining_rules.val) )
            cljff.setProperty("coulombPower", VariantProperty(coulomb_power.val) )
            cljff.setProperty("shiftDelta", VariantProperty(shift_delta.val) )
        except:
            pass

        try:
            if cutoff_scheme.val == "shift_electrostatics":
                cljff.setShiftElectrostatics(True)

            elif cutoff_scheme.val == "reaction_field":
                cljff.setUseReactionField(True)
                cljff.setReactionFieldDielectric(rf_dielectric.val)

            elif cutoff_scheme.val == "group":
                cljff.setUseGroupCutoff(True)

            else:
                print "WARNING. Unrecognised cutoff scheme. Using \"shift_electrostatics\"."
                gridff.setShiftElectrostatics(True)
        except:
            pass

        system.add(cljff)

    for gridff in gridffs:
        system.add(gridff)

    for intraff in intraffs:
        system.add(intraff)
    
    total_nrg = solute_intraff.components().total() + solute_hard_intraclj.components().total() +\
        solute_todummy_intraclj.components().total(0) + solute_fromdummy_intraclj.components().total(0) +\
        solute_hard_todummy_intraclj.components().total(0) + solute_hard_fromdummy_intraclj.components().total(0) +\
        solute_todummy_fromdummy_intraclj.components().total(0) +\
        solventff.components().total() +\
        solute_hard_solventff.components().total() +\
        solute_todummy_solventff.components().total(0) +\
        solute_fromdummy_solventff.components().total(0) +\
        protein_intraff.components().total() + protein_intraclj.components().total() + protein_solventff.components().total() +\
        solute_hard_proteinff.components().total() +\
        solute_todummy_proteinff.components().total(0) +\
        solute_fromdummy_proteinff.components().total(0) 

    fwd_nrg = solute_fwd_intraff.components().total() + solute_fwd_hard_intraclj.components().total() +\
        solute_fwd_todummy_intraclj.components().total(0) + solute_fwd_fromdummy_intraclj.components().total(0) +\
        solute_fwd_hard_todummy_intraclj.components().total(0) + solute_fwd_hard_fromdummy_intraclj.components().total(0) +\
        solute_fwd_todummy_fromdummy_intraclj.components().total(0) +\
        solventff.components().total() +\
        solute_fwd_hard_solventff.components().total() +\
        solute_fwd_todummy_solventff.components().total(0) +\
        solute_fwd_fromdummy_solventff.components().total(0) +\
        protein_intraff.components().total() + protein_intraclj.components().total() + protein_solventff.components().total() +\
        solute_fwd_hard_proteinff.components().total() +\
        solute_fwd_todummy_proteinff.components().total(0) +\
        solute_fwd_fromdummy_proteinff.components().total(0) 
        
    bwd_nrg = solute_bwd_intraff.components().total() + solute_bwd_hard_intraclj.components().total() +\
        solute_bwd_todummy_intraclj.components().total(0) + solute_bwd_fromdummy_intraclj.components().total(0) +\
        solute_bwd_hard_todummy_intraclj.components().total(0) + solute_bwd_hard_fromdummy_intraclj.components().total(0) +\
        solute_bwd_todummy_fromdummy_intraclj.components().total(0) +\
        solventff.components().total() +\
        solute_bwd_hard_solventff.components().total() +\
        solute_bwd_todummy_solventff.components().total(0) +\
        solute_bwd_fromdummy_solventff.components().total(0) +\
        protein_intraff.components().total() + protein_intraclj.components().total() + protein_solventff.components().total() +\
        solute_bwd_hard_proteinff.components().total() +\
        solute_bwd_todummy_proteinff.components().total(0) +\
        solute_bwd_fromdummy_proteinff.components().total(0) 

    if use_grid.val:
        total_nrg += protein_solvent_fixedff.components().total() + solute_fixedff.components().total()
        fwd_nrg += protein_solvent_fixedff.components().total() + solute_fwd_fixedff.components().total()
        bwd_nrg += protein_solvent_fixedff.components().total() + solute_bwd_fixedff.components().total()

    e_total = system.totalComponent()
    e_fwd = Symbol("E_{fwd}")
    e_bwd = Symbol("E_{bwd}")

    lam = Symbol("lambda")
    lam_fwd = Symbol("lambda_{fwd}")
    lam_bwd = Symbol("lambda_{bwd}")

    system.setComponent( e_total, total_nrg )
    system.setComponent( e_fwd, fwd_nrg )
    system.setComponent( e_bwd, bwd_nrg )

    system.setConstant(lam, 0.0)
    system.setConstant(lam_fwd, 0.0)
    system.setConstant(lam_bwd, 0.0)

    de_fwd = Symbol("dE_{fwd}")
    de_bwd = Symbol("dE_{bwd}")

    system.setComponent( de_fwd, fwd_nrg - total_nrg )
    system.setComponent( de_bwd, total_nrg - bwd_nrg )

    # Add a monitor that calculates the average total energy and average energy
    # deltas - we will collect both a mean average and an zwanzig average
    system.add( "total_energy", MonitorComponent(e_total, Average()) )
    system.add( "dg_fwd", MonitorComponent(de_fwd, FreeEnergyAverage(temperature.val)) )
    system.add( "dg_bwd", MonitorComponent(de_bwd, FreeEnergyAverage(temperature.val)) )

    system.add( PerturbationConstraint(solutes) )

    system.add( ComponentConstraint( lam_fwd, Min( lam + delta_lambda.val, 1 ) ) )
    system.add( ComponentConstraint( lam_bwd, Max( lam - delta_lambda.val, 0 ) ) )

    # Add a monitor that records the value of all energy components
    if nmoves_per_energy.val:
        if nmoves_per_energy.val > 0:
            system.add( "energies", MonitorComponents(RecordValues()), nmoves.val / nmoves_per_energy.val )
    
    # Add a monitor that records the coordinates of the system
    if (lam_val.val < 0.001 or lam_val.val > 0.999):
        if nmoves_per_pdb.val:
            if nmoves_per_pdb.val > 0:
                system.add( "trajectory", TrajectoryMonitor(MGName("traj")), nmoves.val / nmoves_per_pdb.val )
    elif not (nmoves_per_pdb_intermediates.val is None):
        if nmoves_per_pdb_intermediates.val > 0:
            system.add( "trajectory", TrajectoryMonitor(MGName("traj")), nmoves.val / nmoves_per_pdb_intermediates.val )

    # Alpha constraints for the soft force fields

    system.add( PropertyConstraint( "alpha0", FFName("solute_todummy_intraclj"), lam ) )
    system.add( PropertyConstraint( "alpha0", FFName("solute_fromdummy_intraclj"), 1 - lam ) )
    system.add( PropertyConstraint( "alpha0", FFName("solute_hard:todummy_intraclj"), lam ) )
    system.add( PropertyConstraint( "alpha0", FFName("solute_hard:fromdummy_intraclj"), 1 - lam ) )
    system.add( PropertyConstraint( "alpha0", FFName("solute_todummy:fromdummy_intraclj"), Min( lam, 1 - lam )  ) ) 
    system.add( PropertyConstraint( "alpha0", FFName("solute_todummy:solvent"), lam ) )
    system.add( PropertyConstraint( "alpha0", FFName("solute_fromdummy:solvent"), 1 - lam ) )
    system.add( PropertyConstraint( "alpha0", FFName("solute_todummy:protein"), lam ) )
    system.add( PropertyConstraint( "alpha0", FFName("solute_fromdummy:protein"), 1 - lam ) )

    system.add( PropertyConstraint( "alpha0", FFName("solute_fwd_todummy_intraclj"), lam_fwd ) )
    system.add( PropertyConstraint( "alpha0", FFName("solute_fwd_fromdummy_intraclj"), 1 - lam_fwd ) )
    system.add( PropertyConstraint( "alpha0", FFName("solute_fwd_hard:todummy_intraclj"), lam_fwd ) )
    system.add( PropertyConstraint( "alpha0", FFName("solute_fwd_hard:fromdummy_intraclj"), 1 - lam_fwd ) )
    system.add( PropertyConstraint( "alpha0", FFName("solute_fwd_todummy:fromdummy_intraclj"), Min( lam_fwd, 1 - lam_fwd ) ) ) 
    system.add( PropertyConstraint( "alpha0", FFName("solute_fwd_todummy:solvent"), lam_fwd ) )
    system.add( PropertyConstraint( "alpha0", FFName("solute_fwd_fromdummy:solvent"), 1 - lam_fwd ) )
    system.add( PropertyConstraint( "alpha0", FFName("solute_fwd_todummy:protein"), lam_fwd ) )
    system.add( PropertyConstraint( "alpha0", FFName("solute_fwd_fromdummy:protein"), 1 - lam_fwd ) )

    system.add( PropertyConstraint( "alpha0", FFName("solute_bwd_todummy_intraclj"), lam_bwd ) )
    system.add( PropertyConstraint( "alpha0", FFName("solute_bwd_fromdummy_intraclj"), 1 - lam_bwd ) )
    system.add( PropertyConstraint( "alpha0", FFName("solute_bwd_hard:todummy_intraclj"), lam_bwd ) )
    system.add( PropertyConstraint( "alpha0", FFName("solute_bwd_hard:fromdummy_intraclj"), 1 - lam_bwd ) )
    system.add( PropertyConstraint( "alpha0", FFName("solute_bwd_todummy:fromdummy_intraclj"), Min( lam_bwd, 1 - lam_bwd ) ) ) 
    system.add( PropertyConstraint( "alpha0", FFName("solute_bwd_todummy:solvent"), lam_bwd ) )
    system.add( PropertyConstraint( "alpha0", FFName("solute_bwd_fromdummy:solvent"), 1 - lam_bwd ) )
    system.add( PropertyConstraint( "alpha0", FFName("solute_bwd_todummy:protein"), lam_bwd ) )
    system.add( PropertyConstraint( "alpha0", FFName("solute_bwd_fromdummy:protein"), 1 - lam_bwd ) )

    system.setComponent( lam, lam_val.val )

    return system


def getAtomNearCOG( molecule ):

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


def setupMoves(system, random_seed):

    moves = WeightedMoves()

    solutes = system[ MGName("solutes") ]
    solute_ref = system[ MGName("solute_ref") ]

    mobilewater = system[ MGName("mobilewater") ]
    mobileion = system[ MGName("mobileion") ]

    protein = system[ MGName("protein") ]
    residues = system[ MGName("residues") ]
    bbresidues = system[ MGName("bbresidues") ]
    
    print "Setting up moves..."
    # Setup Moves
    solute_moves = RigidBodyMC( solutes ) 
    solute_moves.setMaximumTranslation(solutes[MolIdx(0)].molecule().property('flexibility').translation() )
    solute_moves.setMaximumRotation(solutes[MolIdx(0)].molecule().property('flexibility').rotation() )

    # Find solute atom nearest to the center of geometry
    # (note that in new Sire code, the solute will be rotated using its "center" property, so
    #  this nearestcog_atom center will not be used)
    nearestcog_atom = getAtomNearCOG( solutes[MolIdx(0)].molecule() )

    #print nearestcog_atom
    solute_moves.setCenterOfRotation( GetCOGPoint( nearestcog_atom.name() ) )
    solute_moves.setSynchronisedTranslation(True)
    solute_moves.setSynchronisedRotation(True)
    #solute_moves.setSharedRotationCenter(True)
    moves.add( solute_moves, solute_mc_weight.val / 2 )

    solute_intra_moves = InternalMoveSingle( solute_ref )

    perturbed_solutes = system[ MGName("perturbed_solutes") ]
    # Each molecule in perturbed_solutes will have its coordinates set to those 
    # of solute_ref after the move
    solute_intra_moves.setSynchronisedCoordinates(perturbed_solutes)
    moves.add( solute_intra_moves, solute_mc_weight.val / 2)
    
    # Solvent moves, split in water and ions
    if mobilewater.nMolecules() > 0:
        water_moves = RigidBodyMC( PrefSampler(solute_ref[MolIdx(0)], 
                                           mobilewater, pref_constant.val) )    
        water_moves.setMaximumTranslation(max_solvent_translation.val)
        water_moves.setMaximumRotation(max_solvent_rotation.val)
        water_moves.setReflectionSphere(sphere_center, sphere_radius.val)
        moves.add( water_moves, mobilewater.nMolecules()*solvent_mc_weight_factor.val )

    #ion_moves = RigidBodyMC( PrefSampler(solute_ref[MolIdx(0)], ion, pref_constant.val) )
    if mobileion.nMolecules() > 0:
        ion_moves = RigidBodyMC( PrefSampler(solute_ref[MolIdx(0)], mobileion, pref_constant.val) )
        ion_moves.setMaximumTranslation(max_solvent_translation.val)
        ion_moves.setMaximumRotation(max_solvent_rotation.val)
        ion_moves.setReflectionSphere(sphere_center, sphere_radius.val)
        moves.add( ion_moves, ion.nMolecules() )

    # Protein intra moves
    if residues.nMolecules() > 0:
        protein_intra_moves = ZMatMove( PrefSampler(solute_ref[MolIdx(0)], residues, pref_constant.val) ) 
        moves.add( protein_intra_moves, protein_mc_weight.val / 2)
    
    # Now add protein backbone moves
    if bbresidues.nMolecules() > 0:
        bbmoves = RigidBodyMC(bbresidues)
        bbmoves.setMaximumTranslation(0.025*angstrom)
        bbmoves.setMaximumRotation(1*degrees)
        bbmoves.setCenterOfRotation( GetCOGPoint( AtomName("CA", CaseInsensitive),
                                                  AtomName("N", CaseInsensitive) ) )
    
        moves.add( bbmoves, protein_mc_weight.val / 2 )

    moves.setTemperature(temperature.val)

    if (not random_seed):
	random_seed = RanGenerator().randInt(100000,1000000)
	print "Generated random seed number %d " % random_seed

    moves.setGenerator( RanGenerator(random_seed) )
    
    return moves


def writeComponents(components, filename):
    """This function writes the energy components to the file 'filename'"""

    symbols = components.monitoredComponents()

    if len(symbols) == 0:
        return

    newrun = False
    if not os.path.exists(filename):
        newrun = True

    FILE = open(filename, "a")

    nrgs = {}

    for symbol in symbols:
        nrgs[str(symbol)] = components.accumulator(symbol).values()

    symbols = nrgs.keys()
    symbols.sort()

    if newrun:
        print >>FILE,"#step   ",

        for symbol in symbols:
            print >>FILE,"%s   " % symbol,

        print >>FILE,"\n",

    for i in range(0, len(nrgs[symbols[0]])):
        print >>FILE,"%d   " % i,

        for symbol in symbols:
            print >>FILE,"%f   " % nrgs[symbol][i],

        print >>FILE,"\n",


def writeSystemData( system, moves, block):

    nmoves = moves.nMoves()
    monitors = system.monitors()
    outdir = out_dir.val

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if use_grid.val and block == 1:
        # write out the coordinates of the fixed atoms to the output directory
        print "Saving the fixed atoms as %s/fixed_atoms.pdb" % outdir
        PDB().write( system.property("fixed_atoms"), "%s/fixed_atoms.pdb" % outdir )

    try:
        pdb = monitors[MonitorName("trajectory")]
        pdb.writeToDisk("%s/output%0009d.pdb" % (outdir,block))
    except:
        pass

    try:
        energies = monitors[MonitorName("energies")]
        if os.path.exists("%s/energies.dat.bz2" % outdir):
            os.system("bunzip2 -f %s/energies.dat.bz2" % outdir)

        writeComponents( energies, "%s/energies.dat" % outdir )
    except:
        pass

    total_energy = monitors[MonitorName("total_energy")]
    
    dg_fwd = monitors[MonitorName("dg_fwd")]
    dg_bwd = monitors[MonitorName("dg_bwd")]
    
    dg_fwd = dg_fwd.accumulator().average() / delta_lambda.val
    dg_bwd = dg_bwd.accumulator().average() / delta_lambda.val

    system.clearStatistics()
    
    # Ugly
    lam = system.constantExpression(Symbol("lambda")).toString().toDouble()    
    #print dg_bwd, dg_fwd, lam

    if lam < 0.0001:
        dg_bwd = dg_fwd
    elif lam > 0.9999:
        dg_fwd = dg_bwd

    dg_avg = 0.5 * ( dg_fwd + dg_bwd )

    #print dg_avg

    FILE = open("%s/gradients.dat" % outdir , 'a')
    print >>FILE, "%9d %12.8f " % ( block, dg_avg)
   
    FILE = open("moves.dat", "w")
    print >>FILE, "%s" % moves


def printComponents(nrgs):
    keys = nrgs.keys()
    keys.sort()

    for key in keys:
        print "%s     %s" % (key, nrgs[key])

    print "\n",


@resolveParameters
def run():
    print " ### Running a \"bound leg\" single topology free energy calculation ### "

    timer = QTime()
    timer.start()

    # Setup the system from scratch if no restart file is available

    if not os.path.exists("%s/%s" % (out_dir.val,restart_file.val)):
        print "New run. Loading input and creating restart"
        print "Lambda is %5.3f" % lam_val.val       

        amber = Amber()

        print "Reading in coordinate and topology file..."
        molecules, space = amber.readCrdTop(crd_file.val, top_file.val)

        print "Creating the simulation system..."
        system = createSystem(molecules, space)

        print "Setting up the forcefields..."
        system = setupForcefields(system, space)

        print "Setting up the moves..."
        moves = setupMoves(system, random_seed.val)
        print "Saving restart"

        if not os.path.exists(out_dir.val):
            os.makedirs(out_dir.val)

        Sire.Stream.save( [system, moves], "%s/%s" % (out_dir.val,restart_file.val) )

    system, moves = Sire.Stream.load("%s/%s" % (out_dir.val,restart_file.val))
    print "Loaded a restart file on wich we have performed %d moves." % moves.nMoves()
    block_number = moves.nMoves() / nmoves.val  + 1
    s1 = timer.elapsed()/1000.
    print "Setup took %d s " % ( s1 )

    # Run a short simulation

    print "Performing simulation for block number %d " % block_number

    if print_nrgs.val:
        printComponents(system.energies())
    
    system = moves.move(system, nmoves.val, True)

    s2 = timer.elapsed()/1000.
    print "Simulation took %d s " % ( s2 - s1)

    # Update statistics and save restart
    writeSystemData(system, moves, block_number) 

    Sire.Stream.save( [system, moves], "%s/%s" % (out_dir.val,restart_file.val) )

    # Compress some output files
    outpdb = "%s/output%0009d.pdb" % (out_dir.val,block_number)
    if os.path.exists(outpdb):
        os.system( "%s %s/output%0009d*" % (compress.val, out_dir.val, block_number) )
    if os.path.exists("energies.dat"):
        os.system(" %s %s/energies.dat" % (out_dir.val,compress.val) )
