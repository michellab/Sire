#
# Evaluates electrostatics corrections to free energy changes
#
import os,sys, random
import math
from Sire.Tools.OpenMMMD import *
from Sire.Tools import Parameter, resolveParameters
from Sire.Tools.LJcutoff import getFreeEnergy, resample

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
ion_residues = ["Cl-","Na+"]
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

    for molnum in molnums:
        mol = other_solutes.molecule(molnum).molecule()
        if solute.contains(mol):
            other_solutes.remove(mol)
        if mol.residues()[0].name().value() in solvent_residues:
            other_solutes.remove(mol)
        if mol.residues()[0].name().value() in ion_residues:
            other_solutes.remove(mol)

    print ("There are %s mols left " % other_solutes.nMolecules())

    inter_nonbondedff = InterGroupCLJFF("solute_red:other_solutes")
    if (cutoff_type.val != "nocutoff"):
        inter_nonbondedff.setUseReactionField(True)
        inter_nonbondedff.setReactionFieldDielectric(rf_dielectric.val)

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
    for molnum in mol_numbers:
        mol = molecules.molecule(molnum).molecule()
        res0 = mol.residues()[0]
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

def PoissonPBC(binary, solutes, space, cutoff, dielectric, framenum, solute_ref, zerorefcharges=False):

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
    dielec = dielectric
    cut = cutoff/10.0

    iterations=200# CHANGE ME BACK TO 200 !

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
    for molnum in molnums:
        mol = sol_mols.molecule(molnum).molecule()
        atoms = mol.atoms()
        for atom in atoms:
            if (mol.name() == mol_solref.name() and zerorefcharges == True):
                charge = 0.0
            else:
                charge = atom.property("charge").value()
            sigma = atom.property("LJ").sigma().value()
            radius = 0.5*(sigma*2**(1/6.))/10.0 # to nm. Is this a decently good approximation?
            coords = atom.property("coordinates")/10.0 # to nm 
            line = "%8.5f %8.5f %8.5f %8.5f %8.5f\n" % (charge, radius, coords[0], coords[1], coords[2])
            #import pdb; pdb.set_trace()
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

    cmd = "%s > PB.out" % binary
    os.system(cmd)

    cmd = "tail -1 PB.out > temp"
    os.system(cmd)
    rstream = open('temp','r')
    buffer = rstream.readlines()
    elems = buffer[0].split()
    nrg = float(elems[1])

    os.chdir("../")

    cmd = "rm -rf %s" % poissondir
    #os.system(cmd)

    DF_BA_PBC = nrg * kJ_per_mol

    return DF_BA_PBC

def PoissonNP(binary, solutes, dielectric,framenum, space,\
              solute_ref, zerorefcharges=False):

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

    xmlinpart2 = ""

    for molnum in molnums:
        mol = sol_mols.molecule(molnum).molecule()
        nions += mol.nAtoms()
        atoms = mol.atoms()
        for atom in atoms:
            name = atom.name().value()
            if (mol.name() == mol_solref.name() and zerorefcharges == True):
                charge = 0.0
            else:
                charge = atom.property("charge").value()
            #charge = atom.property("charge").value()
            sigma = atom.property("LJ").sigma().value()
            radius = 0.5*(sigma*2**(1/6.)) # Ang. Is this a decently good approximation?
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
### BORN ION SOLVATION ENERGY
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
    srad 1.4
    swin 0.3
    sdens 10.0
    temp 298.15
    calcenergy total
    calcforce no
    write pot gz potential
    # write pot dx potential
    # write charge dx charge
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
    srad 1.4
    swin 0.3
    sdens 10.0
    temp 298.15
    calcenergy total
    calcforce no
end

# COMBINE TO GIVE SOLVATION ENERGY
print elecEnergy solvated - reference end

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

    cmd = "%s apbs.in 1> apbs.out 2> apbs.err" % (binary)
    os.system(cmd)

    rstream = open("apbs.out","r")
    buffer = rstream.readlines()

    nrg = 0.0
    nrgfound = False
    for line in buffer:
        if line.find("Global net ELEC" ) > 0:
            elems = line.split()
            #print (line)
            #print (elems)
            nrg = float(elems[5])
            nrgfound = True
            break
    if not nrgfound:
        print ("ERROR. Could not find electrostatic energy in apbs.out file. ABORT")
        sys.exit(-1)

    os.chdir("../")
    cmd = "rm -rf %s" % poissondir
    #os.system(cmd)

    DF_CB_NP = nrg * kJ_per_mol

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
    print ("q_trace is %s (in e nm2)" % q_trace)

    return q_trace

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
        print ("PSUM using mol %s (natoms %s )" % (mol, mol.nAtoms()) )
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

    system_solute_host_cb = setupInterCoulFF(system_solute_host_cb, Cartesian(), \
                                             cut_type="nocutoff")

    #import pdb; pdb.set_trace()

    delta_func_nrgs = []
    delta_dir_nrgs = []
    DG_pols = []
    DG_psums = []

    while (current_frame <= end_frame):
        print ("Processing frame %s " % current_frame)
        print ("CURRENT POSITION %s " % mdtraj_trajfile.tell() )
        frames_xyz, cell_lengths, cell_angles = mdtraj_trajfile.read(n_frames=1)
        system = updateSystemfromTraj(system, frames_xyz, cell_lengths, cell_angles)
        #import pdb; pdb.set_trace()
        # Now filter out solvent molecules
        solutes, solvent, ions = SplitSoluteSolvent(system)
        solutes, solvent = centerAll(solutes, solvent, system.property("space"))
        PDB().write(solutes, "solutes-%s.pdb" % current_frame )
        #print (solutes.molecules().first().molecule().property("coordinates").toVector()[0])
        #import pdb; pdb.set_trace()
        # ???Should we center solutes to the center of the box???
        # Compute free energy corrections
        # ############################################
        # Use PH' solver to compute DF^{BA}_{PBC}
        # First calculation - using all charges
        print ("Poisson PBC calculation... ")
        #DG_BA_PBC = 0.0 * kcal_per_mol
        DG_BA_PBC_HG = PoissonPBC(PoissonPBCSolverBin.val ,solutes, \
                                  system.property("space"),cutoff_dist.val.value(),\
                                  model_eps.val,
                                  current_frame,solute_ref, zerorefcharges=False)
        # Second calculation - setting charges of solute_ref to 0.
        DG_BA_PBC_H = PoissonPBC(PoissonPBCSolverBin.val ,solutes, \
                                  system.property("space"),cutoff_dist.val.value(),\
                                  model_eps.val,
                                  current_frame,solute_ref, zerorefcharges=True)

        # Use APBS to compute DF^{CB}_{infinity}
        # ??? Should we use experimental dielectric constant rather than
        # the one from the model ???
        print ("Poisson NP calculation... ")
        #DG_CB_NP = 0.0 * kcal_per_mol
        # First calculation - using all charges
        DG_CB_NP_HG = PoissonNP(PoissonNPSolverBin.val, solutes, bulk_eps.val,\
                                current_frame, system.property("space"),\
                                solute_ref, zerorefcharges=False)
        # NOTE ! Will generate nan if all the partial charges of the solutes are zero
        if math.isnan(DG_CB_NP_HG.value()):
            DG_CB_NP_PL = 0.0 * kcal_per_mol
        # Second calculation - setting charges of solute_ref to 0
        DG_CB_NP_H = PoissonNP(PoissonNPSolverBin.val, solutes, bulk_eps.val,\
                                current_frame, system.property("space"),\
                                solute_ref, zerorefcharges=True)
        if math.isnan(DG_CB_NP_H.value()):
            DG_CB_NP_H = 0.0 * kcal_per_mol
        #DG_POL = DG_CB_NP - DG_BA_PBC
        system_solute_host_rf.update(solutes)
        system_solute_host_cb.update(solutes)
        Udir_cb = system_solute_host_cb.energy()
        Udir_rf = system_solute_host_rf.energy() 
        #delta_dir_nrgs.append(delta_dir_nrg)
        DG_POL = (DG_CB_NP_HG-DG_CB_NP_H+Udir_cb) - \
                 (DG_BA_PBC_HG-DG_BA_PBC_H+Udir_rf)
        DG_pols.append(DG_POL)
        print ("DG_CB_NP_HG IS %s " % DG_CB_NP_HG)
        print ("DG_CB_NP_H IS %s " % DG_CB_NP_H)
        print ("Udir_CB IS %s " % Udir_cb)
        print ("DG_BA_PBC_HG IS %s " % DG_BA_PBC_HG)
        print ("DG_BA_PBC_H IS %s " % DG_BA_PBC_H)
        print ("Udir_rf IS %s " % Udir_rf)
        print ("DG_POL IS %s " % DG_POL)
        #import pdb; pdb.set_trace()
        # ############################################
        # Compute psum
        print ("Psum correction... ")
        DG_PSUM = SummationCorrection(solutes, solvent, solute_ref,\
                                      space, model_rho.val.value(),\
                                      model_eps.val, cutoff_dist.val.value())
        DG_psums.append(DG_PSUM)
        print ("DG_PSUM IS %s" % DG_PSUM)
        # ############################################
        # Compute DG_func
        # Free energy change for changing from a reaction field cutoff to coulombic nocutoff
        # Update system_solute_rf
        system_solute_rf.update(solutes)
        system_solute_cb.update(solutes)
        delta_func_nrg = (system_solute_cb.energy() - system_solute_rf.energy())
        delta_func_nrgs.append(delta_func_nrg)
        #import pdb; pdb.set_trace()
        current_frame += step_frame
        mdtraj_trajfile.seek(current_frame)#step_frame, whence=1)

    # Now compute average POL term and uncertainty
    # Note that POL and DIR are correlated, so should evaluate sum of these terms 
    # for bootstrap
    nvals = len(DG_pols)
    DG_POL_avg = 0.0 * kcal_per_mol
    dev_POL = 0.0
    for x in range(0,nvals):
        DG_POL_avg += DG_pols[x]
        dev_POL += (DG_pols[x].value())**2
    DG_POL_avg /= nvals
    dev_POL = (dev_POL / nvals) - (DG_POL_avg.value())**2
    dev_POL = math.sqrt(dev_POL)
    # Now do the same for PSUM
    nvals = len(DG_psums)
    DG_PSUM_avg = 0.0 * kcal_per_mol
    dev_PSUM = 0.0
    for x in range(0,nvals):
        DG_PSUM_avg += DG_psums[x]
        dev_PSUM += (DG_psums[x].value())**2
    DG_PSUM_avg /= nvals
    dev_PSUM = (dev_PSUM / nvals) - (DG_PSUM_avg.value())**2
    dev_PSUM = math.sqrt(dev_PSUM)
    # Now we do the FUNC free energy correction
    #print (delta_func_nrgs)
    DG_FUNC = getFreeEnergy(delta_func_nrgs)
    #print (deltaG)
    nbootstrap = 100
    deltaG_bootstrap = np.zeros(nbootstrap)
    for x in range(0,nbootstrap):
        resampled_nrgs = resample(delta_func_nrgs)
        dG = getFreeEnergy(resampled_nrgs)
        deltaG_bootstrap[x] = dG.value()
    dev_FUNC = deltaG_bootstrap.std()
    print ("DG_POL = %8.5f +/- %8.5f kcal/mol (1 sigma) " % (DG_POL_avg.value(), dev_POL))
    print ("DG_PSUM = %8.5f +/- %8.5f kcal/mol (1 sigma) " % (DG_PSUM_avg.value(), dev_PSUM))
    print ("DG_FUNC = %8.5f +/- %8.5f kcal/mol (1 sigma) " % (DG_FUNC.value(), dev_FUNC))
