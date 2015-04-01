#!/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import re

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

import Sire.Stream

from Sire.Tools import Parameter, resolveParameters
from Sire.Tools.WaterChanger import convertTip3PtoTip4P

###################################
# Parameters used by this module  #
###################################

dobonds = Parameter("move bonds", True, """Whether or not to move the ligands bonds""")
doangles = Parameter("move angles", True, """Whether or not to move the ligands angles""")
dodihedrals = Parameter("move dihedrals", True, """Whether or not to move the ligands dihedrals""")

water_model = Parameter("water model", None,
                        """The water model to use. Note, by default the water model is read from
                           the protein and water crd/top files. If you want to force a change
                           in water model, then set it here, e.g. if you are loading a TIP3P box
                           but want to use TIP4P, then set this parameter to "tip4p".""")

BASE_DIHEDRALH_FLEX = Parameter("h dihedral flex", 30*degrees, "Base dihedral rotation for H")
BASE_DIHEDRAL_FLEX = Parameter("dihedral flex", 20*degrees, "Base dihedral rotation")
BASE_ANGLE_FLEX = Parameter("angle flex", 0.25*degrees, "Base angle rotation")
BASE_BOND_FLEX = Parameter("bond flex", 0.025*angstroms, "Base bond stretch amount")
BASE_TRANSLATION = Parameter("translation", 0.75*angstroms, "Base translation delta amount")
BASE_ROTATION = Parameter("rotation", 30*degrees, "Base rigid body rotation")
BASE_MAXVAR = Parameter("maxvar", 10, "Maximum number of degrees of freedom to move at once")
BASE_MAXVAR_B = Parameter("maxvar bonds", 2, "Maximum number of bonds to move at once")
BASE_MAXVAR_A = Parameter("maxvar angles", 4, "Maximum number of angles to move at once")
BASE_MAXVAR_D = Parameter("maxvar dihedrals", 4, "Maximum number of dihedrals to move at once")

###################################

def getResidueNames(molecule):
    nres = molecule.nResidues()

    resnams = []

    for i in range(0, nres):
        resnams.append( str( molecule.residue(ResIdx(i)).name().value()).upper() )

    return resnams


class NamingScheme:
    def __init__(self):
        self._protein_names = ["GLH", "ILE", "GLN", "GLY", "GLU",
                               "CYS", "HIS", "HID", "SER", "LYS",
                               "LYN", "PRO", "CYX", "HIE", "ASH",
                               "ASN", "HIP", "VAL", "THR", "ASP",
                               "TRP", "PHE", "ALA", "MET", "LEU",
                               "ARG", "TYR", "NME", "ACE"]

        self._water_names = [ "WAT", "T3P", "T4P", "HOH" ]

        self._ion_names = [ "NA+", "Na+", "CA+", "Ca+", "CAL", "CL-", "Cl-" ]

        self._solute_names = [ "LIG" ]

    def proteinsGroupName(self):
        return MGName("protein")

    def solutesGroupName(self):
        return MGName("solute")

    def solventsGroupName(self):
        return MGName("solvent")

    def watersGroupName(self):
        return MGName("water")

    def ionsGroupName(self):
        return MGName("ions")

    def allMoleculesGroupName(self):
        return MGName("all")

    def fixedMoleculesGroupName(self):
        return MGName("fixed_molecules")

    def boundaryMoleculesGroupName(self):
        return MGName("boundary_molecules")

    def mobileProteinSidechainsGroupName(self):
        return MGName("protein_sidechains")

    def mobileProteinBackbonesGroupName(self):
        return MGName("protein_backbones")

    def mobileSolutesGroupName(self):
        return MGName("mobile_solutes")

    def mobileSolventsGroupName(self):
        return MGName("mobile_solvents")

    def addProteinResidueName(self, name):
        self._protein_names.append( name.upper() )

    def addWaterResidueName(self, name):
        self._water_names.append( name.upper() )

    def addSoluteResidueName(self, name):
        self._solute_names.append( name.upper() )

    def addIonResidueName(self, name):
        self._ion_names.append( name.upper() )

    def proteinResidueNames(self):
        return self._protein_names

    def waterResidueNames(self):
        return self._water_names

    def soluteResidueNames(self):
        return self._solute_names

    def ionResidueNames(self):
        return self._ion_names

    def setProteinResidueNames(self, names):
        self._protein_names = []
        for name in names:
            self.addProteinResidueName(name)

    def setWaterResidueNames(self, names):
        self._water_names = []
        for name in names:
            self.addWaterResidueName(name)

    def setSoluteResidueNames(self, name):
        self._solute_names = []
        for name in names:
            self.addSoluteResidueName(name)

    def setIonResidueNames(self, name):
        self._ion_names = []
        for name in names:
            self.addIonResidueName(name)

    def _isType(self, molecule, names, max_residues = None):
        try:
            resnams = getResidueNames(molecule)
        except:
            resnams = molecule

        if max_residues:
            if len(resnams) > max_residues:
                return False

        for resnam in resnams:
            if resnam in names:
                return True

        try:
            if str(molecule.name().value()).upper() in names:
                return True
            else:
                return False
        except:
            return False

    def isProtein(self, molecule):
        return self._isType(molecule, self._protein_names)

    def isWater(self, molecule):
        return self._isType(molecule, self._water_names, 1)

    def isIon(self, molecule):
        return self._isType(molecule, self._ion_names, 1)

    def isSolute(self, molecule):
        return self._isType(molecule, self._solute_names)


def findMolecule(system, molname):
    molecules = system.molecules()

    molname = molname.upper()

    for molnum in molecules.molNums():
        molecule = molecules[molnum].molecule()

        if str(molecule.name().value()).upper() == molname:
            return molecule

        resnams = getResidueNames(molecule)

        for resnam in resnams:
            if resnam == molname:
                return molecule

    return None


def createSystem(top_file, crd_file, naming_scheme = NamingScheme()):
    
    system = System(top_file)

    # Load all of the molecules and their parameters from
    # the topology and coordinate files
    amber = Amber()
    print("Loading the molecules from the Amber files \"%s\" and \"%s\"..." % \
                  (crd_file, top_file))
    (molecules, space) = amber.readCrdTop(crd_file, top_file)

    # If requested, change the water model for all water molecules
    if water_model.val == "tip4p":
        molnums = molecules.molNums()
        new_molecules = Molecules()

        print("Forcing all water molecules to use the %s water model..." % water_model.val)
        print("Converting %d molecules..." % len(molnums))
        i = 0
        for molnum in molnums:
            molecule = molecules[molnum].molecule()

            if i % 100 == 0:
                print("%d" % i)                
                sys.stdout.flush()

            elif i % 10 == 0:
                print(".", end=' ')
                sys.stdout.flush()

            i += 1

            if molecule.nAtoms() == 3:
                # this could be a TIP3P water
                resname =str(molecule.residue().name().value()).lower()

                if resname == "wat" or resname == "t3p":
                    new_molecule = convertTip3PtoTip4P(molecule)
                    if new_molecule:
                        molecule = new_molecule

            new_molecules.add(molecule)

        print("%d" % i)

        molecules = new_molecules

    nmols = molecules.nMolecules()

    print("Number of molecules == %s" % nmols)
    print("System space == %s" % space)

    if nmols == 0:
        return system

    print("Assigning molecules to molecule groups...")
    solute_group = MoleculeGroup(naming_scheme.solutesGroupName().value())
    protein_group = MoleculeGroup(naming_scheme.proteinsGroupName().value())
    solvent_group = MoleculeGroup(naming_scheme.solventsGroupName().value())
    water_group = MoleculeGroup(naming_scheme.watersGroupName().value())
    ion_group = MoleculeGroup(naming_scheme.ionsGroupName().value())
    all_group = MoleculeGroup(naming_scheme.allMoleculesGroupName().value())

    # The all molecules group has all of the molecules
    all_group.add(molecules)

    system.add(all_group)

    # Run through each molecule and decide what type it is...
    molnums = molecules.molNums()
    molnums.sort()

    central_molecule = None

    solutes = []
    proteins = []
    solvents = []
    waters = []
    ions = []

    for molnum in molnums:
        molecule = molecules[molnum].molecule()

        resnams = getResidueNames(molecule)

        if naming_scheme.isSolute(resnams):
            solutes.append(molecule)

        elif naming_scheme.isProtein(resnams):
            proteins.append(molecule)

        elif naming_scheme.isWater(resnams):
            waters.append(molecule)

        elif naming_scheme.isIon(resnams):
            ions.append(molecule)

        elif molecule.nResidues() == 1:
            solvents.append(molecule)

        else:
            solutes.append(molecule)

    # Ok - we have now divided everything up into groups
    for solute in solutes:
        solute_group.add(solute)

    for protein in proteins:
        protein_group.add(protein)

    for water in waters:
        solvent_group.add(water)
        water_group.add(water)

    for solvent in solvents:
        solvent_group.add(solvent)
    
    for ion in ions:
        solvent_group.add(ion)
        ion_group.add(ion)

    if solute_group.nMolecules() > 0:
        system.add(solute_group)

    if protein_group.nMolecules() > 0:
        system.add(protein_group)

    if solvent_group.nMolecules() > 0:
        system.add(solvent_group)

    if water_group.nMolecules() > 0:
        system.add(water_group)

    if ion_group.nMolecules() > 0:
        system.add(ion_group)    

    print("Number of solute molecules == %s" % solute_group.nMolecules()) 
    print("Number of protein molecules == %s" % protein_group.nMolecules())
    print("Number of ions == %s" % ion_group.nMolecules())
    print("Number of water molecules == %s" % water_group.nMolecules())
    print("Number of solvent molecules == %s" % solvent_group.nMolecules())
    print("(solvent group is waters + ions + unidentified single-residue molecules)")

    system.setProperty("space", space)
    system.add( SpaceWrapper( Vector(0), all_group ) )
    system.applyConstraints()

    print("Returning the constructed system")

    return system


def centerSystem(system, molecule):
    print("Setting the origin of the system to the center of molecule %s (%s)..." % (molecule, molecule.number()))
    center = molecule.evaluate().center()
    print("This requires translating everything by %s..." % (-center))    
    
    for molnum in system.molNums():
        molecule = system[molnum].molecule()
        molecule = molecule.move().translate(-center).commit()
        system.update(molecule)

    return system


def guessTranslation( solute ):
    natoms = solute.nAtoms()
    return (BASE_TRANSLATION.val) / ( natoms / 5 + 1)


def guessRotation( solute ):
    natoms = solute.nAtoms()
    sphere_radius = solute.evaluate().boundingSphere().radius()
    return (BASE_ROTATION.val) / ( sphere_radius ** 2)


def generateFlexibility(solute):

    connectivity = solute.property('connectivity')
    all_bonds = connectivity.getBonds()
    all_angles = connectivity.getAngles()
    all_dihedrals = connectivity.getDihedrals()

    flexibility = Flexibility(solute)
    flexibility.setRotation( guessRotation(solute) )
    flexibility.setTranslation( guessTranslation(solute) )

    try:
        flexibility.setMaximumVar( BASE_MAXVAR.val )
    except:
        flexibility.setMaximumBondVar( BASE_MAXVAR_B.val )
        flexibility.setMaximumAngleVar( BASE_MAXVAR_A.val )
        flexibility.setMaximumDihedralVar( BASE_MAXVAR_D.val )

    # Redundant torsions are discarded according to the following algorithm
    # 1) Do not sample a torsion at0-at1-at2-at3 if a variable torsion has 
    # already been defined around at1-at2 or at2-at1.
    # 2) Do not sample a torsion if it would break a ring
    #
    if dodihedrals.val:
        var_dihedrals = []

        for dihedral in all_dihedrals:
            #print dihedral
            tomove = True
            # print dihedral
            at0 = dihedral.atom0()
            at1 = dihedral.atom1()
            at2 = dihedral.atom2()
            at3 = dihedral.atom3()
            # See if a one of the variable dihedral 
            # already rotates around the same torsion
            for vardih in var_dihedrals:
                if ( ( at1 == vardih.atom1() and at2 == vardih.atom2() ) or 
                     ( at2 == vardih.atom1() and at1 == vardih.atom2() ) ):
                    # Yes so will not move this torsion 
                    tomove = False
                    break

            # If still wondering...See if a rotation around this dihedral would break a ring 
            if tomove:
                try:
                    dihbond = BondID(at1, at2)
                    #print dihbond
                    solute.move().change(dihbond,1*degrees)
                except UserWarning as error:
                    # extract the type of the errror
                    error_type = re.search(r"(Sire\w*::\w*)", str(error)).group(0)
                    if error_type == "SireMol::ring_error":
                        # print "This dof would move a ring and is therefore skipped"
                        tomove = False
                    else:
                        # re-throw the exception
                        raise error 

            if tomove:
                # Find out how many atoms would move 
                #print dihedral
                gr0, gr1 = connectivity.split(at1, at2)
                ngr0 = gr0.nSelected()
                ngr1 = gr1.nSelected()

                if (ngr0 <= ngr1):
                    smallgroup = gr0
                else:
                    smallgroup = gr1

                smallgroup = smallgroup.subtract(at1)
                smallgroup = smallgroup.subtract(at2)
                factor = smallgroup.nSelected()

                flexibility.add(dihedral, BASE_DIHEDRAL_FLEX.val/factor)
                var_dihedrals.append(dihedral)

    # And the angles ....
    if doangles.val:
        moved_atoms = []

        for angle in all_angles:
            # print angle
            at0 = angle.atom0()
            at2 = angle.atom2()
            # Do not sample that dof if an existing dof would already move this atom
            if ( ( at0 in moved_atoms) and (at2 in moved_atoms) ):
                continue

            # Test if the angle breaks a ring, if so do not sample it
            try:
                solute.move().change(angle,1*degrees)
            except UserWarning as error:
                # extract the type of the errror
                error_type = re.search(r"(Sire\w*::\w*)", str(error)).group(0)
                if error_type == "SireMol::ring_error":
                    # print "This dof would move a ring and is therefore skipped"
                    continue
                else:
                    # re-throw the exception
                    raise error 

            gr0, gr1 = connectivity.split(at0, angle.atom1(), at2)
            ngr0 = gr0.nSelected()
            ngr1 = gr1.nSelected()

            if (ngr0 <= ngr1):
                smallgroup = gr0
            else:
                smallgroup = gr1

            factor = smallgroup.nSelected()
            flexibility.add(angle, BASE_ANGLE_FLEX.val/factor)

            if at0 not in moved_atoms:
                moved_atoms.append(at0)
            if at2 not in moved_atoms:
                moved_atoms.append(at2)    

    # And the bonds...
    if dobonds.val:
        for bond in all_bonds:
            try:
                solute.move().change(bond,1*angstrom)
            except UserWarning as error:
                # extract the type of the errror
                error_type = re.search(r"(Sire\w*::\w*)", str(error)).group(0)
                if error_type == "SireMol::ring_error":
                    # print "This dof would move a ring and is therefore skipped"
                    continue
                else:
                    # re-throw the exception
                    raise error 

            gr0, gr1 = connectivity.split(bond.atom0(), bond.atom1() )
            ngr0 = gr0.nSelected()
            ngr1 = gr1.nSelected()
            if (ngr0 <= ngr1):
                smallgroup = gr0
            else:
                smallgroup = gr1

            factor = smallgroup.nSelected()
            flexibility.add(bond, BASE_BOND_FLEX.val/factor)

    return flexibility


def getCoordGroup(atoms, coords_property="coordinates"):
    coords = []

    for i in range(0, atoms.count()):
        atom = atoms[i]
        coords.append(atom.property(coords_property))

    return CoordGroup(coords)


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


def addFlexibility(system, reflection_center=None, reflection_radius=None, \
                           naming_scheme=NamingScheme()):
    
    print("Adding flexibility to the system...")

    # create a group for all of the fixed molecules and residues
    fixed_group = MoleculeGroup( naming_scheme.fixedMoleculesGroupName().value() )

    # create a group for the fixed residues that are bonded to the mobile residues
    boundary_group = MoleculeGroup( naming_scheme.boundaryMoleculesGroupName().value() )

    if reflection_center is None or reflection_radius is None:
        print ("No reflection radius or reflection molecule specified, so moving all "
               "molecules and residues in the system.")
        reflection_radius = None
        reflection_center = None

    else:
        print(("Only moving molecules/residues that are within a distance %s A "
               "of the point %s.") % (reflection_radius.value(), reflection_center))

        system.setProperty("reflection center", AtomCoords(CoordGroup(1,reflection_center)))
        system.setProperty("reflection sphere radius", VariantProperty(reflection_radius.to(angstroms)))

    # fit the protein z-matrix templates to all of the protein molecules and add the mobile
    # residues to the mobile_sc_group and mobile_bb_group for mobile sidechains and backbones
    if naming_scheme.proteinsGroupName() in system.mgNames():
        protein_group = system[naming_scheme.proteinsGroupName()]

        # create a zmatrix maker that will be used to build the z-matrices for each protein molecule
        zmat_maker = ZmatrixMaker()
        zmat_maker.loadTemplates( os.path.join(parameter_directory, "amber.zmatrices") )

        # now create the molecule groups that hold the flexible side chains and flexible backbone groups
        mobile_sc_group = MoleculeGroup(naming_scheme.mobileProteinSidechainsGroupName().value())
        mobile_bb_group = MoleculeGroup(naming_scheme.mobileProteinBackbonesGroupName().value())

        # the extra atoms moved as part of a backbone move
        hn_atoms = AtomName("N", CaseInsensitive) * AtomName("H", CaseInsensitive) * \
                   AtomName("HN", CaseInsensitive) * AtomName("HN1", CaseInsensitive) * \
                   AtomName("HN2", CaseInsensitive) * AtomName("HN3", CaseInsensitive)

        # loop over each protein molecule
        for molnum in protein_group.molNums():
            protein_mol = protein_group[molnum].molecule()

            print("Applying residue templates for protein %s" % molnum)

            protein_mol = zmat_maker.applyTemplates(protein_mol)

            system.update(protein_mol)

            if reflection_radius:
                space = Cartesian()

                mobile_resnums = []

                # only move side chains within "sc_radius" and backbones within "bb_radius" of the ligand molecule
                print("Looking for which residues are within the reflection sphere...")
                for i in range(0, protein_mol.nResidues()):
                    res = protein_mol.residue( ResIdx(i) )
                    distance = space.minimumDistance(CoordGroup(1,reflection_center), getCoordGroup(res.atoms()))
          
                    if distance < reflection_radius.value():
                        # add the residue to the mobile sidechains group
                        mobile_sc_group.add(res)
                        mobile_resnums.append( res.number() )

                        # now add the atoms needed from the residue to the mobile backbones group
                        atoms = protein_mol.select(ResIdx(i)).selection()
    
                        if i < (protein_mol.nResidues()-1):
                            try:
                                atoms.deselect( hn_atoms + ResIdx(i) )
                            except:
                                pass

                        if i > 0:
                            try:
                                atoms.select( hn_atoms + ResIdx(i+1) )
                            except:
                                pass

                        mobile_bb_group.add( PartialMolecule(protein_mol, atoms) )

                # now loop over all of the residues and work out which ones are fixed, and which ones
                # are bonded to fixed residues
                connectivity = protein_mol.property("connectivity")

                for i in range(0, protein_mol.nResidues()):
                    res = protein_mol.residue( ResIdx(i) )

                    if not res.number() in mobile_resnums:
                        # is this residue bonded to any of the mobile residues? If so, then it is a boundary residue
                        is_boundary = False

                        for bonded_res in connectivity.connectionsTo( res.number() ):
                            bonded_resnum = protein_mol.residue(bonded_res).number()

                            if bonded_resnum in mobile_resnums:
                                is_boundary = True
                                break

                        if is_boundary:
                            boundary_group.add(res)
                        else:
                            fixed_group.add(res)

            else:
                # assume that the backbone and side chains of all residues are flexible
                for i in range(0, protein_mol.nResidues()):
                    res = protein_mol.residue( ResIdx(i) )
                    mobile_sc_group.add(res)

                    atoms = protein_mol.select(ResIdx(i)).selection()
    
                    if i < (protein_mol.nResidues()-1):
                        try:
                            atoms.deselect( hn_atoms + ResIdx(i) )
                        except:
                            pass

                    if i > 0:
                        try:
                            atoms.select( hn_atoms + ResIdx(i+1) )
                        except:
                            pass

                    mobile_bb_group.add( PartialMolecule(protein_mol, atoms) )

        if mobile_sc_group.nMolecules() > 0:
            system.add(mobile_sc_group)

        if mobile_bb_group.nMolecules() > 0:
            system.add(mobile_bb_group)

        print("The number of residues with flexible sidechains equals %s" % mobile_sc_group.nViews())
        print("The number of residues with flexible backbones equals %s" % mobile_bb_group.nViews())
        print("The number of boundary residues equals %s" % boundary_group.nViews())
        print("The number of fixed residues equals %s" % fixed_group.nViews())

    # add all of the mobile solute molecules to the mobile_solute_group and auto-generate
    # the z-matricies of all of the mobile solutes
    if naming_scheme.solutesGroupName() in system.mgNames():
        solute_group = system[naming_scheme.solutesGroupName()]
        mobile_solute_group = MoleculeGroup( naming_scheme.mobileSolutesGroupName().value() )

        # store the average solute translation and rotation deltas
        avg_trans_delta = 0
        avg_rot_delta = 0

        for molnum in solute_group.molNums():
            solute_mol = solute_group[molnum].molecule()

            move_solute = True

            # Only move the solute if it is within the sphere cutoff of the ligand (if a ligand and solvent 
            # radius have been specified...)
            if reflection_radius:
                move_solute = (Vector.distance(reflection_center, \
                                  solute_mol.evaluate().center()) < reflection_radius.value())

            if move_solute:
                print("\nAuto-detecting the flexible degrees of freedom for solute %s" % molnum)

                # auto-generate the flexibility - bonds, angles and dihedrals
                flexibility = generateFlexibility(solute_mol)

                solute_mol = solute_mol.edit().setProperty("flexibility", flexibility).commit()

                print("\nFlexibility of solute %s equals:" % molnum)
                flex = solute_mol.property("flexibility")
                print(flex)

                avg_trans_delta += flex.translation().to(angstrom)
                avg_rot_delta += flex.rotation().to(degrees)

                system.update(solute_mol)

                mobile_solute_group.add(solute_mol)
            else:
                print("Not moving solute %s as it is outside the spherical solvent cutoff of the ligand." % solute_mol)
                fixed_group.add(solute_mol)

        if mobile_solute_group.nMolecules() > 0:
            system.add(mobile_solute_group)
            system.setProperty("average solute translation delta", \
                                 VariantProperty(avg_trans_delta / mobile_solute_group.nMolecules()))
            system.setProperty("average solute rotation delta", \
                                 VariantProperty(avg_rot_delta / mobile_solute_group.nMolecules()))

        print("\nNumber of mobile solute molecules equals %s" % mobile_solute_group.nMolecules())

    # add all of the mobile solvent molecules to the mobile_solvent_group
    if naming_scheme.solventsGroupName() in system.mgNames():
        solvent_group = system[ naming_scheme.solventsGroupName() ]        

        mobile_solvent_group = MoleculeGroup( naming_scheme.mobileSolventsGroupName().value() )

        print("Adding flexibility to the solvent...")

        if reflection_radius:
            for molnum in solvent_group.molNums():
                solvent_mol = solvent_group[molnum]
                if Vector.distance(reflection_center, solvent_mol.evaluate().center()) < reflection_radius.value():
                    mobile_solvent_group.add(solvent_mol)
                else:
                    fixed_group.add(solvent_mol)

        else:
            mobile_solvent_group.add( solvent_group.molecules() )

        if mobile_solvent_group.nMolecules() > 0:
            system.add(mobile_solvent_group)

        print("\nNumber of mobile solvent molecules equals %s" % mobile_solvent_group.nMolecules())

    # All finished - just need to add in the fixed and boundary groups
    if fixed_group.nMolecules() > 0:
        system.add(fixed_group)

    if boundary_group.nMolecules() > 0:
        system.add(boundary_group)    

    print("\nNumber of fixed (or partially fixed) molecules equals %s" % fixed_group.nMolecules())

    return system

def printGroupInfo(system, group_name):
    try:
        group = system[MGName(group_name)]
        print("%s : nMolecules() == %d" % (str(group), group.nMolecules()))
    except:
        print("There is no group called \"%s\"" % group_name)

