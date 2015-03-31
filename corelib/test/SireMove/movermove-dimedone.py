from Sire.IO import *
from Sire.Mol import *
from Sire.CAS import *
from Sire.Move import *
from Sire.FF import *
from Sire.MM import *
from Sire.Units import *
from Sire.Vol import *
from Sire.Maths import *
from Sire.Base import *
from Sire.System import *
from Sire.Error import *

import os,sys,re
import shutil

temperature = 25*celsius

# Will adjust acceptace ratio using nblocks*nmoves
nmoves_eq = 500
nblocks_eq = 25
# For production
nmoves = 1000
nblocks = 10

solute_file = "test/io/DIMfree.pdb"
solute_name = "DIM"
solute_params = "test/io/DIMfree.par"

protoms_dir = "%s/software/protoms" % os.getenv("HOME")

ff_parameters = [ "%s/parameter/amber99.ff" % protoms_dir,
                  "%s/parameter/gaff.ff" % protoms_dir ]

combining_rules = "arithmetic"

compress = "none"

def createProtoMS():
    """Create and return a ProtoMS parameteriser object"""
    protoms = ProtoMS("%s/protoms2" % protoms_dir)
    for params in ff_parameters:
         protoms.addParameterFile(params)
    return protoms

def histoAnalysis(nbins=1, values=None, output="BOND.dat", mode="bond",avgnrg=0.0):
    """Histogram analysis"""
    if mode == "angle":
        unit = degrees
    else:
        unit = angstrom
    values.sort()
    minval = convertTo(values[0], unit)
    maxval = convertTo(values[-1], unit)
    histo = Histogram(minval,maxval,nbins=nbins)

    for val in values:
        histo.accumulate( convertTo(val, unit) )

    count = 0.0
    for x in range(0,nbins-1):
        bin = histo[x]
        count += bin.value()

    outfile = output
    stream =  open(outfile,'w')
    stream.write("# distribution Average energy %8.5f kcal.mol-1\n" % avgnrg)
    sumbin = 0.0
    for x in range(0,nbins-1):
        bin = histo[x]
        sumbin += bin.value()/count
        stream.write("%8.5f %8.5f %8.5f\n" % (bin.middle(),bin.value()/count,sumbin))   

    
def MoverInspector(solute=None):
    """
    1) Inspect a collection of dofs in a solute inferred from its connectivity
    2) Decide how to sample dofs
    """
    all_bonds = solute.property('connectivity').getBonds()
    all_angles = solute.property('connectivity').getAngles()
    all_dihedrals = solute.property('connectivity').getDihedrals()

    angle_deltas = {}
    bond_deltas = {}
    # Redundant torsions are discarded according to the following algorithm
    # 1) Do not sample a torsion at0-at1-at2-at3 if a variable torsion has 
    # already been defined around at1-at2. 
    # 2) Do not sample a torsion if it would break a ring
    print(all_dihedrals,len(all_dihedrals))
    var_dihedrals = []
    for dihedral in all_dihedrals:
        tomove = True
        print(dihedral)
        at0 = dihedral.atom0()
        at1 = dihedral.atom1()
        at2 = dihedral.atom2()
        at3 = dihedral.atom3()
        # See if a one of the variable dihedral 
        # already rotates around the same torsion
        for vardih in var_dihedrals:
            if ( at1 == vardih.atom1() and at2 == vardih.atom2() ):
                # Yes so will not move this torsion 
                tomove = False
                break
        # See if a rotation around this dihedral would break a ring 
        try:
            dihbond = BondID(at1, at2)
            solute.move().change(dihbond,1*degrees)
        except UserWarning as error:
            # extract the type of the errror
            error_type = re.search(r"(Sire\w*::\w*)", str(error)).group(0)
            if error_type == "SireMol::ring_error":
                print("This dof would move a ring and is therefore skipped")
                tomove = False
            else:
                # re-throw the exception
                raise error 
        if tomove:
            var_dihedrals.append(dihedral)
            angle_deltas[dihedral] = 5.0*degrees
    print("TORSIONS TO SAMPLE")
    for var in var_dihedrals:
        print(var)
    # And the angles ....
    print(all_angles,len(all_angles))
    var_angles = []
    moved_atoms = []
    for angle in all_angles:
        print(angle)
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
                print("This dof would move a ring and is therefore skipped")
                continue
            else:
                # re-throw the exception
                raise error 
        var_angles.append(angle)
        angle_deltas[angle] = 0.5*degrees
        if at0 not in moved_atoms:
            moved_atoms.append(at0)
        if at2 not in moved_atoms:
            moved_atoms.append(at2)    
    print(moved_atoms)
    print(var_angles,len(var_angles))
    # And the bonds...
    var_bonds = []
    for bond in all_bonds:
        try:
            solute.move().change(bond,1*angstrom)
        except UserWarning as error:
            # extract the type of the errror
            error_type = re.search(r"(Sire\w*::\w*)", str(error)).group(0)
            if error_type == "SireMol::ring_error":
                print("This dof would move a ring and is therefore skipped")
                continue
            else:
                # re-throw the exception
                raise error 
        var_bonds.append(bond)
        bond_deltas[bond] = 0.01*angstrom

    return var_bonds, var_angles, var_dihedrals, bond_deltas, angle_deltas

# Load and parameterise
solute = PDB().readMolecule(solute_file)
solute = solute.edit().rename(solute_name).commit()

protoms = createProtoMS()
protoms.addParameterFile(solute_params)
solute = protoms.parameterise(solute, ProtoMS.SOLUTE)

# Putting the solute in the system
system = System()
solutemolecules = MoleculeGroup("solute", solute)

system.add(solutemolecules)

# Setting a forcefield 
solute_intraff = InternalFF("solute_intraff")
solute_intraff.add(solutemolecules)
solute_intraclj = IntraCLJFF("solute_intraclj")
solute_intraclj.add(solute)

forcefields = [ solute_intraff, solute_intraclj ]

for forcefield in forcefields:
    system.add(forcefield)

# Setup the Moves 

solute_mover_moves = MoverMove(solutemolecules)
solute_mover_moves.setTemperature(temperature)

var_bonds, var_angles, var_dihedrals, bond_deltas, angle_deltas = MoverInspector(solute=solute)

solute_mover_moves.setBonds(var_bonds)
solute_mover_moves.setBondDeltas(bond_deltas)

solute_mover_moves.setAngles(var_angles)
solute_mover_moves.setDihedrals(var_dihedrals)
solute_mover_moves.setAngleDeltas(angle_deltas)

moves = SameMoves(solute_mover_moves)

#sys.exit(-1)
# Create monitors 

for bond in var_bonds:
    idx0 = bond.atom0().value()
    idx1 = bond.atom1().value()
    bcomp = Symbol("b_%s_%s" % ( bond.atom0().value(),bond.atom1().value() ) )
    system.add( DistanceComponent( bcomp, solute.atom(bond.atom0()), 
                                   solute.atom(bond.atom1()) ) )
    system.add( "b_%s_%s" % (bond.atom0().value(),bond.atom1().value()), 
                MonitorComponent(bcomp, RecordValues()), 10)

for angle in var_angles:
    idx0 = angle.atom0().value()
    idx1 = angle.atom1().value()
    idx2 = angle.atom2().value()
    acomp = Symbol("a_%s_%s_%s" % ( angle.atom0().value(),
                                    angle.atom1().value(),
                                    angle.atom2().value() ) )
    system.add( AngleComponent( acomp, 
                                solute.atom(angle.atom0()), 
                                solute.atom(angle.atom1()),
                                solute.atom(angle.atom2()) ) )
    system.add( "a_%s_%s_%s" % (angle.atom0().value(), 
                                angle.atom1().value(),
                                angle.atom2().value() ), 
                MonitorComponent(acomp, RecordValues()), 10)

for dihedral in var_dihedrals:
    idx0 = dihedral.atom0().value()
    idx1 = dihedral.atom1().value()
    idx2 = dihedral.atom2().value()
    idx3 = dihedral.atom3().value()
    dcomp = Symbol("d_%s_%s_%s_%s" % ( dihedral.atom0().value(),
                                       dihedral.atom1().value(),
                                       dihedral.atom2().value(),
                                       dihedral.atom3().value() ) )
    system.add( DihedralComponent( dcomp, 
                                   solute.atom(dihedral.atom0()), 
                                   solute.atom(dihedral.atom1()),
                                   solute.atom(dihedral.atom2()),
                                   solute.atom(dihedral.atom3()) ))
    system.add( "d_%s_%s_%s_%s" % (dihedral.atom0().value(), 
                                   dihedral.atom1().value(),
                                   dihedral.atom2().value(),
                                   dihedral.atom3().value() ), 
                MonitorComponent(dcomp, RecordValues()), 10)

e_total = system.totalComponent()
system.add( "total_energy", MonitorComponent(e_total, Average()) )

system.add( "trajectory", TrajectoryMonitor(system[MGIdx(0)]), 100 )

# Equilibrate the acceptance rate to lie within 30-50%
# In a more complex simulation, we should also monitor the acceptance rate of the 
# translations/rotations of the solute. Statistics should not be recorded during
# equilibration as changing the delta values breaks detailed balance.
print("Equilibration")
for i in range(0,nblocks_eq):
    print("Running block %d " % i)
    system = moves.move(system, nmoves_eq, False)
    solute_mover = moves[0]
    acc_rate = solute_mover.acceptanceRatio()
    solute_mover.clearStatistics()
    print(acc_rate)
    if (acc_rate > 0.5):
        print("Increasing deltas...")
        stability = 0
        solute_mover.changeDeltas(prob=0.25, scale=1.25)
    elif (acc_rate < 0.3):
        print("Decreasing deltas...")
        stability = 0
        solute_mover.changeDeltas(prob=0.25, scale=0.75)
    else:
        print("Acceptance rate was in range ")
        stability += 1
        if (stability > 5):
            print("The acceptance rate is stable.")
            break
    moves = SameMoves(solute_mover)

print("Production")
for i in range(0,nblocks):
    print("Running block %d " % i)
    system = moves.move(system, nmoves, True)
#sys.exit(-1)
print("Analysis...")
# make an histogram...
nbins=31

for bond in var_bonds:
    bmonitor = "b_%s_%s" % (bond.atom0().value(),bond.atom1().value()) 
    histoAnalysis(nbins=nbins,
                  values=list(system.monitor( MonitorName(bmonitor) ).accumulator().values()),
                  avgnrg = system.monitor( MonitorName("total_energy") ).accumulator().average(),
                  mode="bond",
                  output="BOND-%s.dat" % bmonitor)


for angle in var_angles:
    amonitor = "a_%s_%s_%s" % (angle.atom0().value(), 
                               angle.atom1().value(),
                               angle.atom2().value() ) 
    histoAnalysis(nbins=nbins,
                  values=list(system.monitor( MonitorName(amonitor) ).accumulator().values()),
                  avgnrg = system.monitor( MonitorName("total_energy") ).accumulator().average(),
                  mode="angle",
                  output="ANGLE-%s.dat" % amonitor)

for dihedral in var_dihedrals:
    dmonitor = "d_%s_%s_%s_%s" % (dihedral.atom0().value(), 
                                  dihedral.atom1().value(),
                                  dihedral.atom2().value(),
                                  dihedral.atom3().value() ) 
    histoAnalysis(nbins=nbins,
                  values=list(system.monitor( MonitorName(dmonitor) ).accumulator().values()),
                  avgnrg = system.monitor( MonitorName("total_energy") ).accumulator().average(),
                  mode="angle",
                  output="DIHEDRAL-%s.dat" % dmonitor)


#**  Output the trajectory
system.monitor( MonitorName("trajectory") ).writeToDisk("outputXXXXXX.pdb")
