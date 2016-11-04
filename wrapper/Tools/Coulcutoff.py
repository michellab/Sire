#
# Evaluates electrostatics corrections to free energy changes
#
import os,sys, random
import math
from Sire.Tools.OpenMMMD import *
from Sire.Tools import Parameter, resolveParameters
from Sire.Tools.LJcutoff import getFreeEnergy, resample
from Sire.Units import *

# Convert from kT to kcal_per_mol
kT_to_kcal = 0.5933 # !!! THIS ASSUMES 1.99E-3 * 298.15 K
# Python dependencies
#
try:
    import mdtraj
except ImportError:
    print ("CoulCutoff.py depends on a working install of the python module mdtraj. Please install mdtraj in your sire python.")
    sys.exit(-1)

try:
    import numpy as np
except ImportError:
    print ("CoulCutoff.py depends on a working install of the python module mdtraj. Please install mdtraj in your sire python.")
    sys.exit(-1)


bulk_eps = Parameter("bulk_eps", 78.4,
                     """The dielectric constant of the bulk solvent.""")

model_eps = Parameter("model_eps", 78.4,
                     """The dielectric constant of the modelled solvent.""")

model_rho = Parameter("model_rho", 1.0 * gram/(centimeter*centimeter*centimeter)\
                     ,"""The density of buk solvent.""")

trajfile = Parameter("trajfile", "traj000000001.dcd",
                    """File name of the trajectory to process.""")

stepframe = Parameter("step_frame",1,
    """The number of frames to step to between two succcessive evaluations.""")

PoissonPBCSolverBin = Parameter("PoissonPBCSolverBin","/home/julien/local/bin/pb_generalT","""Path to the PBC Poisson solver.r""")

PoissonNPSolverBin = Parameter("PoissonNPSolverBin","/home/julien/local/APBS-1.4.1-binary/bin/apbs","""Path to the NP Poisson solver.""")

#### Hardcoded parameters (may need revision)
solvent_residues = ["WAT","ZBK","ZBT","CYC"]
ion_residues = ["Cl-"]#,"Na+"]  IMPORTANT RE-INSERT THE NA+
DIME = 97


def setupIntraCoulFF(system, space, cut_type="nocutoff", cutoff= 999* angstrom, dielectric=1.0):

    print ("Creating force fields... ")

    solutes = system[MGName("solutes")]
    solute = system[MGName("solute_ref")]
    solute_hard = system[MGName("solute_ref_hard")]
    solute_todummy = system[MGName("solute_ref_todummy")]
    solute_fromdummy = system[MGName("solute_ref_fromdummy")]

    # Solute intramolecular LJ energy
    solute_hard_intracoul = IntraCLJFF("solute_hard_intracoul")
    solute_hard_intracoul.add(solute_hard)
    if (cut_type != "nocutoff"):
        solute_hard_intracoul.setUseReactionField(True)
        solute_hard_intracoul.setReactionFieldDielectric(dielectric)

    solute_todummy_intracoul = IntraSoftCLJFF("solute_todummy_intracoul")
    solute_todummy_intracoul.setShiftDelta(shift_delta.val)
    solute_todummy_intracoul.setCoulombPower(coulomb_power.val)
    solute_todummy_intracoul.add(solute_todummy)
    if (cut_type != "nocutoff"):
        solute_todummy_intracoul.setUseReactionField(True)
        solute_todummy_intracoul.setReactionFieldDielectric(dielectric)

    solute_fromdummy_intracoul = IntraSoftCLJFF("solute_fromdummy_intracoul")
    solute_fromdummy_intracoul.setShiftDelta(shift_delta.val)
    solute_fromdummy_intracoul.setCoulombPower(coulomb_power.val)
    solute_fromdummy_intracoul.add(solute_fromdummy)
    if (cut_type != "nocutoff"):
        solute_fromdummy_intracoul.setUseReactionField(True)
        solute_fromdummy_intracoul.setReactionFieldDielectric(dielectric)

    solute_hard_todummy_intracoul = IntraGroupSoftCLJFF("solute_hard:todummy_intracoul")
    solute_hard_todummy_intracoul.setShiftDelta(shift_delta.val)
    solute_hard_todummy_intracoul.setCoulombPower(coulomb_power.val)
    solute_hard_todummy_intracoul.add(solute_hard, MGIdx(0))
    solute_hard_todummy_intracoul.add(solute_todummy, MGIdx(1))
    if (cut_type != "nocutoff"):
        solute_hard_todummy_intracoul.setUseReactionField(True)
        solute_hard_todummy_intracoul.setReactionFieldDielectric(dielectric)

    solute_hard_fromdummy_intracoul = IntraGroupSoftCLJFF("solute_hard:fromdummy_intracoul")
    solute_hard_fromdummy_intracoul.setShiftDelta(shift_delta.val)
    solute_hard_fromdummy_intracoul.setCoulombPower(coulomb_power.val)
    solute_hard_fromdummy_intracoul.add(solute_hard, MGIdx(0))
    solute_hard_fromdummy_intracoul.add(solute_fromdummy, MGIdx(1))
    if (cut_type != "nocutoff"):
        solute_hard_fromdummy_intracoul.setUseReactionField(True)
        solute_hard_fromdummy_intracoul.setReactionFieldDielectric(dielectric)

    solute_todummy_fromdummy_intracoul = IntraGroupSoftCLJFF("solute_todummy:fromdummy_intracoul")
    solute_todummy_fromdummy_intracoul.setShiftDelta(shift_delta.val)
    solute_todummy_fromdummy_intracoul.setCoulombPower(coulomb_power.val)
    solute_todummy_fromdummy_intracoul.add(solute_todummy, MGIdx(0))
    solute_todummy_fromdummy_intracoul.add(solute_fromdummy, MGIdx(1))
    if (cut_type != "nocutoff"):
        solute_todummy_fromdummy_intracoul.setUseReactionField(True)
        solute_todummy_fromdummy_intracoul.setReactionFieldDielectric(dielectric)

    # TOTAL
    forcefields = [solute_hard_intracoul, solute_todummy_intracoul, solute_fromdummy_intracoul,
                   solute_hard_todummy_intracoul, solute_hard_fromdummy_intracoul,
                   solute_todummy_fromdummy_intracoul]

    for forcefield in forcefields:
        system.add(forcefield)

    system.setProperty("space", space)
    system.setProperty("switchingFunction", CHARMMSwitchingFunction(cutoff))
    system.setProperty("combiningRules", VariantProperty(combining_rules.val))
    system.setProperty("coulombPower", VariantProperty(coulomb_power.val))
    system.setProperty("shiftDelta", VariantProperty(shift_delta.val))


    total_nrg = solute_hard_intracoul.components().coulomb() + \
                solute_todummy_intracoul.components().coulomb(0) + solute_fromdummy_intracoul.components().coulomb(0) + \
                solute_hard_todummy_intracoul.components().coulomb(0) + solute_hard_fromdummy_intracoul.components().coulomb(0) + \
                solute_todummy_fromdummy_intracoul.components().coulomb(0)

    e_total = system.totalComponent()

    lam = Symbol("lambda")

    system.setComponent(e_total, total_nrg)

    system.setConstant(lam, 0.0)

    system.add(PerturbationConstraint(solutes))

    # NON BONDED Alpha constraints for the soft force fields
    system.add(PropertyConstraint("alpha0", FFName("solute_todummy_intracoul"), lam))
    system.add(PropertyConstraint("alpha0", FFName("solute_fromdummy_intracoul"), 1 - lam))
    system.add(PropertyConstraint("alpha0", FFName("solute_hard:todummy_intracoul"), lam))
    system.add(PropertyConstraint("alpha0", FFName("solute_hard:fromdummy_intracoul"), 1 - lam))
    system.add(PropertyConstraint("alpha0", FFName("solute_todummy:fromdummy_intracoul"), Max(lam, 1 - lam)))

    system.setComponent(lam, lambda_val.val)

    return system

def setupInterCoulFF(system, space, cut_type="nocutoff", cutoff= 999* angstrom, dielectric=1.0):

    solute = system[MGName("solute_ref")]
    other_solutes = system[MGName("molecules")]
    mols = other_solutes.molecules()
    molnums = mols.molNums()

    #mutated_wats = [3333,5400,2682,729,4695,2612,3903,4707,5235,5248,4717]
    for molnum in molnums:
        mol = other_solutes.molecule(molnum).molecule()
        if solute.contains(mol):
            other_solutes.remove(mol)
        #if mol.residues()[0].number().value() in mutated_wats:
        #    continue
        if mol.residues()[0].name().value() in solvent_residues:
            other_solutes.remove(mol)
        if mol.residues()[0].name().value() in ion_residues:
            other_solutes.remove(mol)

    print ("There are %s mols left " % other_solutes.nMolecules())

    inter_nonbondedff = InterGroupCLJFF("solute_red:other_solutes")
    # JM bug in old code. Should be.
    #if (cutoff_type.val != "nocutoff"):
    #    inter_nonbondedff.setUseReactionField(True)
    #    inter_nonbondedff.setReactionFieldDielectric(rf_dielectric.val)
    # This should be correct
    if (cut_type != "nocutoff"):
        inter_nonbondedff.setUseReactionField(True)
        inter_nonbondedff.setReactionFieldDielectric(dielectric)

    inter_nonbondedff.add(solute, MGIdx(0))
    inter_nonbondedff.add(other_solutes, MGIdx(1))

    system.add(inter_nonbondedff)

    system.setProperty("space", space)
    system.setProperty("switchingFunction", CHARMMSwitchingFunction(cutoff))
    system.setProperty("combiningRules", VariantProperty(combining_rules.val))

    total_nrg = inter_nonbondedff.components().coulomb()
    e_total = system.totalComponent()
    system.setComponent(e_total, total_nrg)

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
        mol = system.molecule(molnum).molecule()
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
        mol = system.molecule(molnum).molecule()
        newmol_coords = newmols_coords[molnum]
        mol = mol.edit().setProperty("coordinates", newmol_coords).commit()
        changedmols.add(mol)
    system.update(changedmols)

    space = PeriodicBox(Vector( traj_box_x, traj_box_y, traj_box_z ) )
    system.setProperty("space",space)

    return system

def SplitSoluteSolvent(system):
    molecules = system.molecules()
    mol_numbers = molecules.molNums()
    solutes = MoleculeGroup("solutes")
    solvent = MoleculeGroup("solvent")
    ions = MoleculeGroup("ions")

    #mutated_wats = [3333,5400,2682,729,4695,2612,3903,4707,5235,5248,4717]

    for molnum in mol_numbers:
        mol = molecules.molecule(molnum).molecule()
        res0 = mol.residues()[0]
        #if res0.number().value() in mutated_wats:
        #    print ("Mutating wat %s" % res0.number().value())
        #    mol = mol.atom(AtomName("O")).edit().setProperty("charge",+1 * mod_electron).molecule().commit()
        #    mol = mol.atom(AtomName("H1")).edit().setProperty("charge",+0 * mod_electron).molecule().commit()
        #    mol = mol.atom(AtomName("H2")).edit().setProperty("charge",+0 * mod_electron).molecule().commit()            
        #    solutes.add(mol)
        #    continue
        if res0.name().value() in solvent_residues:
            solvent.add(mol)
        elif res0.name().value() in ion_residues:
            ions.add(mol)
        else:
            solutes.add(mol)
    return solutes, solvent, ions

def centerAll(solutes, solvent, space):

    if space.isPeriodic():
        box_center = space.dimensions()/2
    else:
        box_center = Vector(0.0, 0.0, 0.0)

    solutes_mols = solutes.molecules()
    solutes_cog = CenterOfGeometry(solutes_mols).point()

    delta = box_center - solutes_cog

    molNums = solutes_mols.molNums()
    for molnum in molNums:
        mol = solutes.molecule(molnum).molecule()
        molcoords = mol.property("coordinates")
        molcoords.translate(delta)
        mol = mol.edit().setProperty("coordinates", molcoords).commit()
        solutes.update(mol)

    solvent_mols = solvent.molecules()
    solvmolNums = solvent_mols.molNums()
    for molnum in solvmolNums:
        mol = solvent.molecule(molnum).molecule()
        molcoords = mol.property("coordinates")
        molcoords.translate(delta)
        mol = mol.edit().setProperty("coordinates",molcoords).commit()
        solvent.update(mol)

    return solutes, solvent

def PoissonPBC2(binary, solutes, space, cutoff, dielectric, framenum, solute_ref,\
                zerorefcharges=False, neutatm=None):

    mol_solref = solute_ref.molecules().first().molecule()
    space_x = space.dimensions()[0]/10.0 # nm
    space_y = space.dimensions()[1]/10.0 #
    space_z = space.dimensions()[2]/10.0 #
    nions = 0
    sol_mols = solutes.molecules()
    molnums = sol_mols.molNums()
    for molnum in molnums:
        mol = sol_mols.molecule(molnum).molecule()
        nions += mol.nAtoms()
    if neutatm is not None:
        nions += len(neutatm)

    dielec = dielectric
    cut = cutoff/10.0

    iterations=20# CHANGE ME BACK TO 200 !

    infilepart1 = """GRID
%s %s %s
%8.5f %8.5f %8.5f
4
END
ITERATION
%s 1.5 0.1 -0.001
0 50 50 50 1
END
ELECTRO
%s %8.5f
""" % (DIME, DIME, DIME, space_x, space_y, space_z, iterations, nions, dielec)

    infilepart2 = ""
    # Dict of charges to use for DG calcs, key is index
    probecharges = {}
    nzeroq = 0
    nions = 0
    srad = 0.14
    #srad = 0.0
    for molnum in molnums:
        mol = sol_mols.molecule(molnum).molecule()
        atoms = mol.atoms()
        for atom in atoms:
            if (mol.name() == mol_solref.name()):
                probecharge = atom.property("charge").value()
                probecharges[nions] = probecharge
                if (zerorefcharges == True):
                    charge = 0.0
                    radius = 0.0
                    nzeroq += 1
                else:
                    charge = probecharge
                    sigma = atom.property("LJ").sigma().value()
                    if (charge == 0 and sigma ==0):
                        radius = 0.0
                    else:
                        radius = 0.5*(sigma*2**(1/6.))/10.0 + srad # nm
            else:
                charge = atom.property("charge").value()
                sigma = atom.property("LJ").sigma().value()
                if (charge == 0 and sigma == 0):
                    radius = 0.0
                else:
                    radius = 0.5*(sigma*2**(1/6.))/10.0 + srad# nm. Is this a decently good approximation?
            coords = atom.property("coordinates")/10.0 # to nm
            line = "%8.6f %8.5f %8.5f %8.5f %8.5f\n" % (charge, radius, coords[0], coords[1], coords[2])
            nions += 1
            #import pdb; pdb.set_trace()
            infilepart2 += line
    if neutatm is not None:
        for ion in neutatm:
            line = "%8.6f %8.5f %8.5f %8.5f %8.5f\n" % (ion[1], 0.0, ion[0][0], ion[0][1], ion[0][2])
            nions += 1
            infilepart2 += line
    infilepart2 += "END\n"

    infilepart3 = """BOUNDARY
3
%s
END
""" % (cut)
    infile = infilepart1 + infilepart2 + infilepart3
    #lines = infile.split('\n')
    #for line in lines:
    #    print (line)
    #import pdb; pdb.set_trace()

    if (zerorefcharges == True):
        poissondir = "poisson-pbc-%s-zeroref" % framenum
    else:
        poissondir = "poisson-pbc-%s" % framenum

    cmd = " mkdir -p %s" % poissondir
    os.system(cmd)

    os.chdir(poissondir)

    wstream = open('inFile','w')
    wstream.write(infile)
    wstream.close()

    #import pdb; pdb.set_trace()

    allzeroq = False
    if (nions == nzeroq):
        allzeroq = True

    if ( not allzeroq ):
        cmd = "%s > PB.out" % binary
        os.system(cmd)
        # Or is it phi_het - phi_hom ?
        phi_het = loadPots("atompot.dat")
        # Compute free energies
        #import pdb; pdb.set_trace()
        DG_HET = elecEnergy(probecharges, phi_het)
    else:
        DG_HET = 0.0
    DG_HOM = 0.0
    os.chdir("../")

    cmd = "rm -rf %s" % poissondir
    #os.system(cmd)

    DF_BA_PBC = (DG_HET-DG_HOM) * kT_to_kcal * kcal_per_mol

    return DF_BA_PBC

def loadPots(txtfile):
    stream = open(txtfile,'r')
    buffer = stream.readlines()
    stream.close()
    c = 0
    pots = {}
    for line in buffer:
        if line.startswith("#"):
            continue
        val = float(line)
        pots[c] = val
        c += 1
    return pots

def elecEnergy(charges, potentials):
    # Charges assumed to be in fractional electron units
    # Pots assumed to be in kT/e
    indices = list( charges.keys() )
    indices.sort()
    nrg = 0.0
    for idx in indices:
        q = charges[idx]
        p = potentials[idx]
        nrg += q*p
    nrg *= 0.5
    return nrg

def PoissonNP2(binary, solutes, dielectric,framenum, space,\
               solute_ref, zerorefcharges=False, neutatm=None):

    mol_solref = solute_ref.molecules().first().molecule()
    # Here write input file for apbs...
    DF_CB_NP = 0.0 * kcal_per_mol
    space_x = space.dimensions()[0]
    space_y = space.dimensions()[1]
    space_z = space.dimensions()[2]

    nions = 0
    sol_mols = solutes.molecules()
    molnums = sol_mols.molNums()


    xmlinpart1 = """# Only the x, y, z, charge, and radii fields are read.
<ion-example>
  <residue>
     <resName>ION</resName>
     <chainID>A</chainID>
     <resSeq>1</resSeq>
"""

    srad = 1.4
    #srad = 0.0
    xmlinpart2 = ""
    # Dict of charges to use for DG calcs, key is index
    probecharges = {}
    nzeroq = 0
    for molnum in molnums:
        mol = sol_mols.molecule(molnum).molecule()
        atoms = mol.atoms()
        for atom in atoms:
            name = atom.name().value()
            if (mol.name() == mol_solref.name()):
                probecharge = atom.property("charge").value()
                probecharges[nions] = probecharge
                if (zerorefcharges == True):
                    charge = 0.0
                    radius = 0.0
                    nzeroq += 1
                else:
                    charge = probecharge
                    sigma = atom.property("LJ").sigma().value()
                    if (charge == 0 and sigma == 0):
                        radius = 0
                    else:
                        radius = 0.5*(sigma*2**(1/6.))+srad # Ang
            else:
                charge = atom.property("charge").value()
                sigma = atom.property("LJ").sigma().value()
                if (charge == 0 and sigma == 0):
                    radius = 0
                else:
                    radius = 0.5*(sigma*2**(1/6.)) + srad# Ang. Is this a decently good approximation?
            coords = atom.property("coordinates") # APBS uses angstroms
            xmlinpart2 += "     <atom>\n"
            xmlinpart2 += "        <serial>ATOM</serial>\n"
            xmlinpart2 += "        <name>%s</name>\n" % name
            xmlinpart2 += "        <x>%s</x>\n" % coords[0]
            xmlinpart2 += "        <y>%s</y>\n" % coords[1]
            xmlinpart2 += "        <z>%s</z>\n" % coords[2]
            xmlinpart2 += "        <charge>%s</charge>\n" % charge
            xmlinpart2 += "        <radius>%s</radius>\n" % radius
            xmlinpart2 += "     </atom>\n"
            nions += 1

    if neutatm is not None:
        for ion in neutatm:
            xmlinpart2 += "     <atom>\n"
            xmlinpart2 += "        <serial>ATOM</serial>\n"
            xmlinpart2 += "        <name>ATM</name>\n"
            xmlinpart2 += "        <x>%s</x>\n" % ion[0][0]
            xmlinpart2 += "        <y>%s</y>\n" % ion[0][1]
            xmlinpart2 += "        <z>%s</z>\n" % ion[0][2]
            xmlinpart2 += "        <charge>%s</charge>\n" % ion[1]
            xmlinpart2 += "        <radius>0.000</radius>\n"
            xmlinpart2 += "     </atom>\n"
            nions += 1

    xmlinpart3 = """  </residue>
</ion-example>
"""

    xmlin = xmlinpart1 + xmlinpart2 + xmlinpart3

    if (zerorefcharges ==True):
        poissondir = "poisson-np-%s-zeroref" % framenum
    else:
        poissondir = "poisson-np-%s" % framenum

    cmd = "mkdir -p %s" % poissondir
    os.system(cmd)

    os.chdir(poissondir)

    wstream = open("apbssystem.xml","w")
    wstream.write(xmlin)
    wstream.close()

    apbsin="""#############################################################################
### ATOMIC POTENTIALS
### $Id$
###
### Please see APBS documentation (http://apbs.sourceforge.net/doc/) for
### input file sytax.
#############################################################################

# READ IN MOLECULES
read
    mol xml apbssystem.xml
end

# COMPUTE POTENTIAL FOR SOLVATED STATE
elec name solvated
    mg-manual
    dime %s %s %s
    grid %8.5f %8.5f %8.5f
    gcent %8.5f %8.5f %8.5f
    mol 1
    lpbe
    bcfl mdh
    pdie 1.0
    sdie %s
    chgm spl2
    srfm mol
    srad 0.0
    swin 0.3
    sdens 10.0
    temp 298.15
    calcenergy total
    calcforce no
    write atompot flat atompotential
end

# COMPUTE POTENTIAL FOR REFERENCE STATE
elec name reference
    mg-manual
    dime %s %s %s
    grid %8.5f %8.5f %8.5f
    gcent %8.5f %8.5f %8.5f
    mol 1
    lpbe
    bcfl mdh
    pdie 1.0
    sdie 1.0
    chgm spl2
    srfm mol
    srad 0.0
    swin 0.3
    sdens 10.0
    temp 298.15
    calcenergy total
    calcforce no
    write atompot flat refatompotential
end

quit
    """ % (DIME, DIME, DIME,
           space_x/DIME, space_y/DIME, space_z/DIME,
           space_x/2.0, space_y/2.0, space_z/2.0, dielectric,
           DIME, DIME, DIME,
           space_x/DIME, space_y/DIME, space_z/DIME,
           space_x/2.0, space_y/2.0, space_z/2.0)

    wstream = open("apbs.in","w")
    wstream.write(apbsin)
    wstream.close()

    allzeroq = False
    if nions == nzeroq:
        allzeroq = True

    if (not allzeroq):
        cmd = "%s apbs.in 1> apbs.out 2> apbs.err" % (binary)
        os.system(cmd)
        # Now load electrostatic potentials
        phi_het = loadPots("atompotential.txt")
        phi_hom = loadPots("refatompotential.txt")
        #import pdb; pdb.set_trace()
        # Compute free energies
        DG_HET = elecEnergy(probecharges, phi_het)
        DG_HOM = elecEnergy(probecharges, phi_hom)
    else:
        DG_HET = 0.0
        DG_HOM = 0.0

    os.chdir("../")
    cmd = "rm -rf %s" % poissondir
    #os.system(cmd)

    DF_CB_NP = (DG_HET - DG_HOM) * kT_to_kcal * kcal_per_mol
    #if (zerorefcharges == True):
    #    import pdb; pdb.set_trace()

    return DF_CB_NP

def calcQuadrupoleTrace(molecule):
    atoms = molecule.atoms()
    com = [0.0, 0.0, 0.0]
    totmass = 0.0
    for atom in atoms:
        coords = atom.property("coordinates")
        mass = atom.property("mass").value()
        totmass += mass
        for i in range(0,3):
            com[i] += coords[i]*mass
    for i in range(0,3):
        com[i] /= totmass

    q_trace = 0.0
    for atom in atoms:
        coords = atom.property("coordinates")
        charge = atom.property("charge")
        ri_x = coords[0] - com[0]
        ri_y = coords[1] - com[1]
        ri_z = coords[2] - com[2]
        ri2 = ri_x**2 + ri_y**2 + ri_z**2
        cont = charge.value() * ri2
        #print ("cont %s " % cont)
        q_trace += charge.value() * ri2
    q_trace /= 100.0
    #print ("q_trace is %s (in e nm2)" % q_trace)

    return q_trace


def NPCoulPots(solutes, solute_ref):
    mol_solref = solute_ref.molecules().first().molecule()
    sol_mols = solutes.molecules()
    molnums = sol_mols.molNums()
    sol_atoms = mol_solref.atoms()
    nsolatoms = mol_solref.nAtoms()
    pots = {}
    nrg = 0.0
    for i in range(0,nsolatoms):
        ri = sol_atoms[i].property("coordinates")
        qi = sol_atoms[i].property("charge").value()
        poti = 0.0
        for molnum in molnums:
            mol = sol_mols.molecule(molnum).molecule()
            if mol.name() == mol_solref.name():
                continue
            molatoms = mol.atoms()
            molnats = mol.nAtoms()
            for j in range(0,molnats):
                rj = molatoms[j].property("coordinates")
                qj = molatoms[j].property("charge").value()
                rij2 = (rj[0]-ri[0])**2+(rj[1]-ri[1])**2+(rj[2]-ri[2])**2
                rij = math.sqrt(rij2)
                poti += one_over_four_pi_eps0 * (qj/rij)
        pots[i] = poti
        nrg += poti*qi
    print (pots)
    nrg *= 0.5
    print (nrg)
    #sys.exit(-1)

def DirectSummation(solutes, space,
                    cutoff, dielectric,
                    framenum,  solute_ref, zerorefcharges=False,
                    neutatm=None):
    # But really, shouldn't I just compute the intermolecular ligand energy +
    # the intramolecular ligand energy ?
    Udir_cb = 0.0
    Udir_pbc = 0.0
    # Initialise reaction-field parameters
    krf = ( 1 / cutoff**3 ) * ( dielectric - 1.0 ) / ( 2 * dielectric + 1.0 )
    crf = (1 / cutoff ) * (3 * dielectric / ( 2 * dielectric + 1.0 ) )

    # Select atoms for summation
    mol_solref = solute_ref.molecules().first().molecule()
    coords = []
    charges = []
    sol_mols = solutes.molecules()
    molnums = sol_mols.molNums()
    for molnum in molnums:
        mol = sol_mols.molecule(molnum).molecule()
        atoms = mol.atoms()
        for atom in atoms:
            coord = atom.property("coordinates")
            if (mol.name() == mol_solref.name() and zerorefcharges == True):
                charge = 0.0
            else:
                charge = atom.property("charge").value()
            coords.append(coord)
            charges.append(charge)
    if neutatm is not None:
        for ion in neutatm:
            coord = Vector(ion[0])
            charge = ion[1]
            coords.append(coord)
            charges.append(charge)
    # Now have array of atoms to consider for calculation, do double loop
    natoms = len(coords)
    for i in range(0,natoms):
        ri = coords[i]
        qi = charges[i]
        for j in range(i+1,natoms):
            rj = coords[j]
            qj = charges[j]
            # Compute Coulombic energy
            rij2_np = (rj[0]-ri[0])**2+(rj[1]-ri[1])**2+(rj[2]-ri[2])**2
            rij_np = math.sqrt(rij2_np)
            coul_nrg = one_over_four_pi_eps0 * ( (qi*qj) /rij_np )
            # Compute Barker-Watts reaction field energy
            rij_pbc = space.calcDist(ri,rj)
            if rij_pbc > cutoff:
                rf_nrg = 0.0
            else:
                rf_nrg =  qi*qj*one_over_four_pi_eps0 * ( 1/rij_pbc +\
                                                          krf*rij_pbc**2 - crf )
            Udir_cb += coul_nrg
            Udir_pbc += rf_nrg

    Udir_cb = Udir_cb * kcal_per_mol
    Udir_pbc = Udir_pbc * kcal_per_mol
    return Udir_cb, Udir_pbc

def ExcludedInteractions(solutes, space,
                         cutoff, dielectric,
                         framenum,  solute_ref, zerorefcharges=False):
    Udir_cb = 0.0
    Udir_pbc = 0.0
    # Initialise reaction-field parameters
    krf = ( 1 / cutoff**3 ) * ( dielectric - 1.0 ) / ( 2 * dielectric + 1.0 )
    crf = (1 / cutoff ) * (3 * dielectric / ( 2 * dielectric + 1.0 ) )
    # Select atoms for summation
    mol_solref = solute_ref.molecules().first().molecule()
    coords = []
    charges = []
    sol_mols = solutes.molecules()
    molnums = sol_mols.molNums()
    mol_atoms = []
    for molnum in molnums:
        mol = sol_mols.molecule(molnum).molecule()
        if (mol.name() != mol_solref.name()):
            continue
        atoms = mol.atoms()
        CLJNB = mol.property("intrascale")
        for atom in atoms:
            mol_atoms.append(atom)
            #coord = atom.property("coordinates")
            #charge = atom.property("charge").value()
            #coords.append(coord)
            #charges.append(charge)
    # Now have array of atoms to consider for calculation, do double loop
    natoms = len(mol_atoms)
    for i in range(0,natoms):
        cgi = mol_atoms[i].cgAtomIdx()
        ri = mol_atoms[i].property("coordinates")
        qi = mol_atoms[i].property("charge").value()
        for j in range(i+1,natoms):
            cgj = mol_atoms[j].cgAtomIdx()
            scale_fac = CLJNB( cgi, cgj)
            if (scale_fac.coulomb() == 0):
                # 1,2 or 1,3
                rj = mol_atoms[j].property("coordinates")
                qj = mol_atoms[j].property("charge").value()
                # Compute Coulombic energy
                rij2_np = (rj[0]-ri[0])**2+(rj[1]-ri[1])**2+(rj[2]-ri[2])**2
                rij_np = math.sqrt(rij2_np)
                coul_nrg = one_over_four_pi_eps0 * ( (qi*qj) /rij_np )
                # Compute Barker-Watts reaction field energy
                rij_pbc = space.calcDist(ri,rj)
                if rij_pbc > cutoff:
                    rf_nrg = 0.0
                else:
                    rf_nrg =  qi*qj*one_over_four_pi_eps0 * ( 1/rij_pbc +\
                                                              krf*rij_pbc**2 - crf )
            Udir_cb += coul_nrg
            Udir_pbc += rf_nrg
            # Also deal with 1,4 ?
    Udir_cb = Udir_cb * kcal_per_mol
    Udir_pbc = Udir_pbc * kcal_per_mol
    return Udir_cb, Udir_pbc

def SummationCorrection(solutes, solvent, solute_ref, space, rho_solvent_model,
                        eps_solvent_model, BAcutoff):
    nrg = 0.0
    mol_solref = solute_ref.molecules().first().molecule()
    sol_mols = solutes.molecules()
    molnums = sol_mols.molNums()

    for molnum in molnums:
        mol = sol_mols.molecule(molnum).molecule()
        if mol.name() != mol_solref.name():
            continue
        #print ("PSUM using mol %s (natoms %s )" % (mol, mol.nAtoms()) )
        #import pdb; pdb.set_trace()
        atoms = mol.atoms()
        mol_charge = 0.0
        mol_com = [0.0, 0.0, 0.0]
        mol_mass = 0.0
        for atom in atoms:
            name = atom.name().value()
            charge = atom.property("charge").value()
            mass = atom.property("mass").value()
            coords = atom.property("coordinates")
            # What about PBC ?
            mol_com[0] += coords[0]*mass
            mol_com[1] += coords[1]*mass
            mol_com[2] += coords[2]*mass
            mol_charge += charge
            mol_mass += mass
        for x in range(0,3):
            mol_com[x] = mol_com[x]/mol_mass
        #print (mol_com)
        #print (mol_charge)
        # JM 01/16 This code only deals with a single solute == solute_ref
        break

    nsolv = 0
    solv_mols = solvent.molecules()
    solvmolnums = solv_mols.molNums()


    # Now that we know the mol_com, compute the quadrupole trace
    quadrupole_trace = calcQuadrupoleTrace(solv_mols.first().molecule())

    for molnum in solvmolnums:
        solvent = solv_mols.molecule(molnum).molecule()
        # Find com
        solv_com = [0.0, 0.0, 0.0]
        solv_atoms = solvent.atoms()
        solv_mass = 0.0
        for solv_atom in solv_atoms:
            coords = solv_atom.property("coordinates")
            mass = solv_atom.property("mass").value()
            solv_mass += mass
            for i in range(0,3):
                solv_com[i] += coords[i]*mass
        for i in range(0,3):
            solv_com[i] /= solv_mass
        #print (mol_com)
        #print (solv_com)
        d = space.calcDist(Vector(mol_com), Vector(solv_com))
        #print (d)
        if d < BAcutoff:
            nsolv += 1
    #print (nsolv)
    ONE_OVER_6PI_EPS0 = 290.98622868361923
    nrg = -ONE_OVER_6PI_EPS0 * mol_charge * quadrupole_trace *\
        ( ( (2*(eps_solvent_model-1) / (2*eps_solvent_model+1) )*\
              nsolv/(4*pi*(BAcutoff/10.0)**3/3.0) ) +\
              (3/(2*eps_solvent_model+1)))
    #sys.exit(-1)
    DF_PSUM = nrg * kJ_per_mol
    return DF_PSUM

def SummationCorrection2(solutes, solvent, solute_ref, space, rho_solvent_model,
                         eps_solvent_model, BAcutoff):
    nrg_tot = 0.0
    solv_mols = solvent.molecules()
    solvmolnums = solv_mols.molNums()
    # Now that we know the mol_com, compute the quadrupole trace
    quadrupole_trace = calcQuadrupoleTrace(solv_mols.first().molecule())

    solv_coms = []

    for molnum in solvmolnums:
        solvent = solv_mols.molecule(molnum).molecule()
        # Find com
        solv_com = [0.0, 0.0, 0.0]
        solv_atoms = solvent.atoms()
        solv_mass = 0.0
        for solv_atom in solv_atoms:
            coords = solv_atom.property("coordinates")
            mass = solv_atom.property("mass").value()
            solv_mass += mass
            for i in range(0,3):
                solv_com[i] += coords[i]*mass
        for i in range(0,3):
            solv_com[i] /= solv_mass
        solv_coms.append(solv_com)

    mol_solref = solute_ref.molecules().first().molecule()
    sol_mols = solutes.molecules()
    molnums = sol_mols.molNums()

    for molnum in molnums:
        mol = sol_mols.molecule(molnum).molecule()
        if mol.name() != mol_solref.name():
            continue
        print ("PSUM using mol %s (natoms %s )" % (mol, mol.nAtoms()) )
        #import pdb; pdb.set_trace()
        atoms = mol.atoms()
        #mol_charge = 0.0
        #mol_com = [0.0, 0.0, 0.0]
        #mol_mass = 0.0
        for atom in atoms:
            nsolv = 0
            name = atom.name().value()
            sol_charge = atom.property("charge").value()
            #mass = atom.property("mass").value()
            sol_coords = atom.property("coordinates")
            for solv_com in solv_coms:
                d = space.calcDist(sol_coords, Vector(solv_com))
                if d < BAcutoff:
                    nsolv += 1
            #print (nsolv)
            ONE_OVER_6PI_EPS0 = 290.98622868361923
            nrg = -ONE_OVER_6PI_EPS0 * sol_charge * quadrupole_trace *\
                  ( ( (2*(eps_solvent_model-1) / (2*eps_solvent_model+1) )*\
                      nsolv/(4*pi*(BAcutoff/10.0)**3/3.0) ) +\
                    (3/(2*eps_solvent_model+1)))
            #print (sol_charge,nsolv,nrg)
            nrg_tot += nrg
            #print("Average molecules: %d" %nsolv)
    # JM 04/16 This code only deals with a single solute == solute_ref
        break

    DF_PSUM = nrg_tot * kJ_per_mol
    return DF_PSUM

def genNeutAtmosphere(solutes, solute_ref):
    # Compute net charge of all solutes excluding solute_ref
    sol_mols = solutes.molecules()
    mol_solref = solute_ref.molecules().first().molecule()
    molnums = sol_mols.molNums()
    netcharge = 0.0
    min_coord = [99999,99999,99999]
    max_coord = [-99999,-99999,-99999]
    sol_coords = []
    for molnum in molnums:
        mol = sol_mols.molecule(molnum).molecule()
        atoms = mol.atoms()
        for atom in atoms:
            coords = atom.property("coordinates")
            #import pdb; pdb.set_trace()
            sigma = atom.property("LJ").sigma().value()
            radius = max(0.5*(sigma*2**(1/6.)),1.0)#mininum radius is 1 Ang
            if coords.x() > max_coord[0]:
                max_coord[0] = coords.x()
            if coords.y() > max_coord[1]:
                max_coord[1] = coords.y()
            if coords.z() > max_coord[2]:
                max_coord[2] = coords.z()
            if coords.x() < min_coord[0]:
                min_coord[0] = coords.x()
            if coords.y() < min_coord[1]:
                min_coord[1] = coords.y()
            if coords.z() < min_coord[2]:
                min_coord[2] = coords.z()
            sol_coords.append( [[coords.x(), coords.y(), coords.z()], radius])
            if (mol.name() == mol_solref.name() ):
                continue
            else:
                netcharge += atom.property("charge").value()
    print ("The environment netcharge is %s " % netcharge)
    print ("min_coord %s max_coord %s " % (min_coord,max_coord))
    min_coord[0] = min_coord[0] - 5.0 # Space for extra points
    min_coord[1] = min_coord[1] - 5.0 # Space for extra points
    min_coord[2] = min_coord[2] - 5.0 # Space for extra points
    max_coord[0] = max_coord[0] + 5.0 # Space for extra points
    max_coord[1] = max_coord[1] + 5.0 # Space for extra points
    max_coord[2] = max_coord[2] + 5.0 # Space for extra points
    netcharge = int(netcharge)
    if (netcharge == 0):
        # No couter ions atmosphere needed
        return None
    # Embedd solutes into a 3D grid with a 0.5 Angstrom spacing
    step = 1.0
    nx = int( (max_coord[0] - min_coord[0])/step) + 1
    ny = int( (max_coord[1] - min_coord[1])/step) + 1
    nz = int( (max_coord[2] - min_coord[2])/step) + 1
    print ("nx %s ny %s nz %s" % (nx,ny,nz))
    grid = []
    for x in range(0,nx):
        for y in range(0,ny):
            for z in range(0,nz):
                gx = min_coord[0] + x * step
                gy = min_coord[1] + y * step
                gz = min_coord[2] + z * step
                # Datastructure is coordinates, active flag. All points are active and not within ionic radius initially
                grid.append( [ [gx, gy, gz],1,0 ] )

    # Compute the concave hull of solutes?
    # Disable grid points that are inside the hull
    # Disable grid points that are within VDW radius of a solute atom
    # Disable grid points that are more than VDW + cutoff_plus OR cutoff_minus of 
    # all solute atoms
    if netcharge < 0:# Use ionic radius of Sodium
        ionrad = 1.4
    elif netcharge > 0:# Use ionic radius of Chloride
        ionrad = 2.3
    c = 0
    for atom_coord in sol_coords:
        c += 1
        print ("Doing atom %s ..." % c)
        #if c == 100:
        #    break
        # if distance less than radmin assume steric clash
        radmin = atom_coord[1]+0.5*ionrad
        # if distance greater than radmax assume no contact
        radmax = atom_coord[1]+1.5*ionrad
        radmin2 = radmin**2
        radmax2 = radmax**2
        for point in grid:
            if (point[1] == 0):
                continue
            d2 = (point[0][0]-atom_coord[0][0])**2 +\
                 (point[0][1]-atom_coord[0][1])**2 +\
                 (point[0][2]-atom_coord[0][2])**2
            if d2 < radmin2:
                point[1] = 0
                continue
            if d2 < radmax2:
                point[2] = 1
            #print (d2,atom_coord,point)
            #import pdb; pdb.set_trace()
    new_points = []
    for point in grid:
        if (point[1] == 1 and point[2] == 1):
            new_points.append(point)
    # Create a particle at each remaining grid point with a zero radius and a
    # partial charge of -net_charge/ngrid_points
    pcharge= -netcharge/len(new_points)
    for point in new_points:
        point[1] = pcharge
    # Return particles list
    return new_points



@resolveParameters
def runLambda():
    try:
        host = os.environ['HOSTNAME']
    except KeyError:
        host = "unknown"
    print("### Running electrostatics correction calculation on %s ###" % host)
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

    # What to do with this...
    system = createSystemFreeEnergy(molecules)
    lam = Symbol("lambda")
    solutes = system[MGName("solutes")]
    solute_ref = system[MGName("solute_ref")]
    system.setConstant(lam, lambda_val.val)
    system.add(PerturbationConstraint(solutes))
    system.setComponent(lam, lambda_val.val)
    #system = createSystemFreeEnergy(molecules)

    system_solute_rf = System()
    system_solute_rf.add(solutes)
    system_solute_rf.add(system[MGName("solute_ref")])
    system_solute_rf.add(system[MGName("solute_ref_hard")])
    system_solute_rf.add(system[MGName("solute_ref_todummy")])
    system_solute_rf.add(system[MGName("solute_ref_fromdummy")])

    system_solute_rf = setupIntraCoulFF(system_solute_rf, space, \
                                        cut_type=cutoff_type.val,
                                        cutoff=cutoff_dist.val,
                                        dielectric=model_eps.val)
    #import pdb; pdb.set_trace()

    system_solute_cb = System()
    system_solute_cb.add(solutes)
    system_solute_cb.add(system[MGName("solute_ref")])
    system_solute_cb.add(system[MGName("solute_ref_hard")])
    system_solute_cb.add(system[MGName("solute_ref_todummy")])
    system_solute_cb.add(system[MGName("solute_ref_fromdummy")])

    system_solute_cb = setupIntraCoulFF(system_solute_cb, Cartesian(), \
                                        cut_type="nocutoff")


    system_solute_host_rf = System()
    system_solute_host_rf.add(system[MGName("molecules")])
    system_solute_host_rf.add(system[MGName("solute_ref")])
    system_solute_host_rf.add(system[MGName("solute_ref")])

    system_solute_host_rf = setupInterCoulFF(system_solute_host_rf, space, \
                                             cut_type=cutoff_type.val,
                                             cutoff=cutoff_dist.val,
                                             dielectric=model_eps.val)

    system_solute_host_cb = System()
    system_solute_host_cb.add(system[MGName("molecules")])
    system_solute_host_cb.add(system[MGName("solute_ref")])

    system_solute_host_cb = setupInterCoulFF(system_solute_host_cb, Cartesian(),\
                                             cut_type="nocutoff")#,cutoff=22*angstrom)#cutoff=cutoff_dist.val)

    #import pdb; pdb.set_trace()

    delta_func_nrgs = []
    delta_dir_nrgs = []
    DG_pols = []
    DG_dirs = []
    DG_psums = []
    DG_excluded = []

    # Now loop over snapshots in dcd and accumulate energies
    start_frame = 1
    end_frame = 1000000000
    step_frame = stepframe.val

    mdtraj_trajfile = mdtraj.open(trajfile.val,'r')
    nframes = len(mdtraj_trajfile)
    if end_frame > (nframes - 1):
        end_frame = nframes - 1
    mdtraj_trajfile.seek(start_frame)
    current_frame = start_frame

    neutatmosphere = True

    print ("#FrameNum DG_CH[P+L,L]^CB,NBC DG_CH[P,L]^CB,NBC UDIR[P+L]^CB,NBC UDIR[P]^CB,NBC"\
           "DG_CH[P+L,L]^BA,PBC DG_CH[P,L]^BA,PBC UDIR[P+L]^BA,PBC UDIR[P]^BA,PBC, DG_PSUM")

    while (current_frame <= end_frame):
        print ("#Processing frame %s " % current_frame)
        frames_xyz, cell_lengths, cell_angles = mdtraj_trajfile.read(n_frames=1)
        system = updateSystemfromTraj(system, frames_xyz, cell_lengths, cell_angles)
        #import pdb; pdb.set_trace()
        # Now filter out solvent molecules
        solutes, solvent, ions = SplitSoluteSolvent(system)
        solutes, solvent = centerAll(solutes, solvent, system.property("space"))
        # Now generate neutralising atmosphere
        neutatm = None
        if (neutatmosphere):
            neutatm = genNeutAtmosphere(solutes, solute_ref)
            if neutatm is not None:
                stream = open("ions.xyz","w")
                stream.write("%s\n \n" % len(neutatm))
                for ion in neutatm:
                    line = "ION %12.5f %12.5f %12.5f\n" % (ion[0][0],ion[0][1],ion[0][2])
                    stream.write(line)
                stream.close()
        #import pdb; pdb.set_trace()
        PDB().write(solutes, "solutes-%s.pdb" % current_frame )
        #print (solutes.molecules().first().molecule().property("coordinates").toVector()[0])
        #NPCoulPots(solutes, solute_ref)
        # ???Should we center solutes to the center of the box???
        # Compute free energy corrections
        # ############################################
       # Use APBS to compute DF^{CB}_{infinity}
        # ??? Should we use experimental dielectric constant rather than
        # the one from the model ???
        print ("#Poisson NP calculation... ")
        DG_CB_NP_HG = 0.0 * kcal_per_mol
        # First calculation - using all charges
        # JM 04/16 MODIFY TO GET ATOMIC POTENTIALS
        DG_CB_NP_HG = PoissonNP2(PoissonNPSolverBin.val, solutes, bulk_eps.val,\
                                current_frame, system.property("space"),\
                                 solute_ref, zerorefcharges=False,
                                 neutatm=neutatm)
        # NOTE ! Will generate nan if all the partial charges of the solutes are zero
        if math.isnan(DG_CB_NP_HG.value()):
            DG_CB_NP_HG = 0.0 * kcal_per_mol
        # Second calculation - setting charges of solute_ref to 0
        #DG_CB_NP_H = 0.0 * kcal_per_mol
        DG_CB_NP_H = PoissonNP2(PoissonNPSolverBin.val, solutes, bulk_eps.val,\
                               current_frame, system.property("space"),\
                                solute_ref, zerorefcharges=True,
                                neutatm=neutatm)
        #import pdb; pdb.set_trace()
        if math.isnan(DG_CB_NP_H.value()):
            DG_CB_NP_H = 0.0 * kcal_per_mol
        
        # Use PH' solver to compute DF^{BA}_{PBC}
        # First calculation - using all charges
        print ("#Poisson PBC calculation... ")
        DG_BA_PBC_HG = 0.0 * kcal_per_mol
        #DG_BA_PBC_HG = PoissonPBC(PoissonPBCSolverBin.val ,solutes, \
        #                          system.property("space"),cutoff_dist.val.value(),\
        #                          model_eps.val,
        #                          current_frame,solute_ref, zerorefcharges=False)
        DG_BA_PBC_HG = PoissonPBC2(PoissonPBCSolverBin.val ,solutes, \
                                   system.property("space"),cutoff_dist.val.value(),\
                                   model_eps.val,
                                   current_frame,solute_ref, zerorefcharges=False,
                                   neutatm=neutatm)
        # Second calculation - setting charges of solute_ref to 0.
        #DG_BA_PBC_H = 0.0 * kcal_per_mol
        #DG_BA_PBC_H = PoissonPBC(PoissonPBCSolverBin.val ,solutes, \
        #                          system.property("space"),cutoff_dist.val.value(),\
        #                          model_eps.val,
        #                          current_frame,solute_ref, zerorefcharges=True)
        DG_BA_PBC_H = PoissonPBC2(PoissonPBCSolverBin.val ,solutes, \
                                  system.property("space"),cutoff_dist.val.value(),\
                                  model_eps.val,
                                  current_frame,solute_ref, zerorefcharges=True,
                                  neutatm=neutatm)
        #import pdb; pdb.set_trace()
        system_solute_host_rf.update(solutes)
        system_solute_host_cb.update(solutes)
        Udirinter_cb = system_solute_host_cb.energy()
        Udirinter_rf = system_solute_host_rf.energy()
        # Direct summation code
        print ("#Direct Summation")
        #Udir_cb_hg = 0.0 * kcal_per_mol
        #Udir_pbc_hg = 0.0 * kcal_per_mol
        #Udir_cb_h = 0.0 * kcal_per_mol
        #Udir_pbc_h = 0.0 * kcal_per_mol
        Udir_cb_hg, Udir_pbc_hg = DirectSummation(solutes, system.property("space"),
                                                  cutoff_dist.val.value(), model_eps.val,
                                                  current_frame, solute_ref, zerorefcharges=False,
                                                  neutatm=neutatm)
        Udir_cb_h, Udir_pbc_h = DirectSummation(solutes, system.property("space"),
                                                cutoff_dist.val.value(), model_eps.val,
                                                current_frame, solute_ref, zerorefcharges=True,
                                                neutatm=neutatm)
        #delta_dir_nrgs.append(delta_dir_nrg)
        #DG_PB = (DG_CB_NP_HG + DG_CB_NP_H + Udir_cb_hg - Udir_cb_h) - \
        #         (DG_BA_PBC_HG + DG_BA_PBC_H + Udir_pbc_hg - Udir_pbc_h)
        #DG_POL = (DG_CB_NP_HG + DG_CB_NP_H) - (DG_BA_PBC_HG + DG_BA_PBC_H)
        #DG_pols.append(DG_POL)
        #DG_DIR = Udirinter_cb - Udirinter_rf
        #DG_dirs.append(DG_DIR)
        #print ("DG_CB_NP_HG IS %s " % DG_CB_NP_HG)
        #print ("DG_CB_NP_H IS %s " % DG_CB_NP_H)
        #print ("Udir_CB_hg IS %s " % Udir_cb_hg)
        #print ("Udir_CB_h IS %s " % Udir_cb_h)
        #print ("DG_BA_PBC_HG IS %s " % DG_BA_PBC_HG)
        #print ("DG_BA_PBC_H IS %s " % DG_BA_PBC_H)
        #print ("Udir_pbc_hg IS %s " % Udir_pbc_hg)
        #print ("Udir_pbc_h  IS %s " % Udir_pbc_h)
        #print ("Udirinter_cb is %s " %  Udirinter_cb)
        #print ("Udirinter_rf is %s " %  Udirinter_rf)
        #print ("DG_POL IS %s " % DG_POL)
        #print ("DG_DIR IS %s " % DG_DIR)
        #print ("DG_PB IS %s " % DG_PB)
        # Now compute term for excluded atoms
        #system_solute_rf.update(solutes)
        #system_solute_cb.update(solutes)
        #Uexc_cb_g, Uexc_pbc_g = ExcludedInteractions(solutes, system.property("space"),
        #                                        cutoff_dist.val.value(), model_eps.val,
        #                                        current_frame, solute_ref)
        #print ("Uexc_cb_g IS %s " % Uexc_cb_g)
        #print ("Uexc_pbc_g IS %s " % Uexc_pbc_g)
        #DG_EXC = Uexc_cb_g - Uexc_pbc_g
        #print ("DG_EXC IS %s " % DG_EXC)
        #DG_excluded.append(DG_EXC)
        #import pdb; pdb.set_trace()
        # ############################################
        # Compute psum
        print ("#Psum correction... ")
        #DG_PSUM = SummationCorrection(solutes, solvent, solute_ref,\
        #                              space, model_rho.val.value(),\
        #                              model_eps.val, cutoff_dist.val.value())
        DG_PSUM = SummationCorrection2(solutes, solvent, solute_ref,\
                                      space, model_rho.val.value(),\
                                      model_eps.val, cutoff_dist.val.value())
        #DG_psums.append(DG_PSUM)
        #print ("DG_PSUM IS %s" % DG_PSUM)
        # ############################################
        # Compute DG_func
        # Free energy change for changing from a reaction field cutoff to coulombic nocutoff
        # Update system_solute_rf
        #system_solute_rf.update(solutes)
        #system_solute_cb.update(solutes)
        #delta_func_nrg = (system_solute_cb.energy() - system_solute_rf.energy())
        #delta_func_nrgs.append(delta_func_nrg)
        #import pdb; pdb.set_trace()
        print ("%8d %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f" % \
               (current_frame,DG_CB_NP_HG.value(),DG_CB_NP_H.value(), Udir_cb_hg.value(), Udir_cb_h.value(),\
                DG_BA_PBC_HG.value(),DG_BA_PBC_H.value(),Udir_pbc_hg.value(),Udir_pbc_h.value(),
                Udirinter_cb.value(),Udirinter_rf.value(),DG_PSUM.value()))
        current_frame += step_frame
        mdtraj_trajfile.seek(current_frame)#step_frame, whence=1)
    # Now compute average POL term and uncertainty
    # Note that POL and DIR are correlated, so should evaluate sum of these terms
    # for bootstrap
    #nvals = len(DG_pols)
    #DG_POL_avg = 0.0 * kcal_per_mol
    #dev_POL = 0.0
    #for x in range(0,nvals):
    #    DG_POL_avg += DG_pols[x]
    #    dev_POL += (DG_pols[x].value())**2
    #DG_POL_avg /= nvals
    #dev_POL = (dev_POL / nvals) - (DG_POL_avg.value())**2
    #dev_POL = math.sqrt(dev_POL)
    ## Now do the same for DG_DIR
    #nvals = len(DG_dirs)
    #DG_DIR_avg = 0.0 * kcal_per_mol
    #dev_DIR = 0.0
    #for x in range(0,nvals):
    #    DG_DIR_avg += DG_dirs[x]
    #    dev_DIR += (DG_dirs[x].value())**2
    #DG_DIR_avg /= nvals
    #dev_DIR = (dev_DIR / nvals) - (DG_DIR_avg.value())**2
    #dev_DIR = math.sqrt(dev_DIR)
    ## Now do the same for PSUM
    #nvals = len(DG_psums)
    #DG_PSUM_avg = 0.0 * kcal_per_mol
    #dev_PSUM = 0.0
    #for x in range(0,nvals):
    #    DG_PSUM_avg += DG_psums[x]
    #    dev_PSUM += (DG_psums[x].value())**2
    #DG_PSUM_avg /= nvals
    #dev_PSUM = (dev_PSUM / nvals) - (DG_PSUM_avg.value())**2
    #dev_PSUM = math.sqrt(dev_PSUM)
    # Now do the same for excluded
    #nvals = len(DG_excluded)
    #DG_EXC_avg = 0.0 * kcal_per_mol
    #dev_EXC = 0.0
    #for x in range(0,nvals):
    #   DG_EXC_avg += DG_excluded[x]
    #    dev_EXC += (DG_excluded[x].value())**2
    #DG_EXC_avg /= nvals
    #dev_EXC = (dev_EXC / nvals) - (DG_EXC_avg.value())**2
    #dev_EXC = math.sqrt(dev_EXC)
    # Now do the same for 'ALL' = POL + DIR + PSUM
    #nvals = len(DG_pols)
    #DG_ALL_avg = 0.0 * kcal_per_mol
    #dev_ALL = 0.0
    #for x in range(0,nvals):
    #    val = DG_pols[x] + DG_dirs[x] + DG_psums[x]
    #    DG_ALL_avg += val
    #    dev_ALL += (val.value())**2
    #DG_ALL_avg /= nvals
    #dev_ALL = (dev_ALL / nvals) - (DG_ALL_avg.value())**2
    #dev_ALL = math.sqrt(dev_ALL)
    ## Now we do the FUNC free energy correction
    ##print (delta_func_nrgs)
    #DG_FUNC = getFreeEnergy(delta_func_nrgs)
    ##print (deltaG)
    #nbootstrap = 100
    #deltaG_bootstrap = np.zeros(nbootstrap)
    #for x in range(0,nbootstrap):
    #    resampled_nrgs = resample(delta_func_nrgs)
    #    dG = getFreeEnergy(resampled_nrgs)
    #    deltaG_bootstrap[x] = dG.value()
    #dev_FUNC = deltaG_bootstrap.std()
    #print ("DG_POL = %8.5f +/- %8.5f kcal/mol (1 sigma) " % (DG_POL_avg.value(), dev_POL))
    #print ("DG_DIR = %8.5f +/- %8.5f kcal/mol (1 sigma) " % (DG_DIR_avg.value(), dev_DIR))
    #print ("DG_PSUM = %8.5f +/- %8.5f kcal/mol (1 sigma) " % (DG_PSUM_avg.value(), dev_PSUM))
    #print ("DG_COR = %8.5f +/- %8.5f kcal/mol ( 1 sigma, nsamples %s) " % (DG_ALL_avg.value(),dev_ALL,len(DG_pols)))
    #print ("DG_EXC = %8.5f +/- %8.5f kcal/mol (1 sigma) " % (DG_EXC_avg.value(), dev_EXC))
    #print ("DG_FUNC = %8.5f +/- %8.5f kcal/mol (1 sigma) " % (DG_FUNC.value(), dev_FUNC))
