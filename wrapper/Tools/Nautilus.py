"""Nautilus is a python module that implements the grid cell theory (GCT) method.
GCT enables analysis of explicit solvent molecular dynamics simulations to compute and spatially resolve enthalpies and entropies of
water molecules

Details of the methodology are available in the following publications:

Evaluation of Host/Guest Binding Thermodynamics of Model Cavities with Grid Cell Theory
Michel, J. ; Henchman, R.H. ; Gerogiokas, G. ; Southey, M. W. Y. ; Mazanetz, M. P. ; Law, R. J. ;
J. Chem. Theory Comput., 10 (9), 4055-4068, 2014

Prediction of Small Molecule Hydration Thermodynamics with Grid Cell Theory
Gerogiokas, G. ; Calabro, G. ; Henchman, R.H. ; Southey, M. W. Y. ; Law, R. J. ; Michel, J.
J. Chem. Theory Comput., 10 (1), 35 - 48, 2014

Gerogiokas, G. ; Southey, M. W. Y. ; Mazanetz, M. P. ; Heifetz, A.; Bodkin, M.;Law, R. J. ; Michel, J.
Phys. Chem. Chem. Phys., 2015. DOI: 10.1039/C4CP05572A

"""
#
# Nautilus module
#
# (C) Georgios Gerogiokas and Julien Michel 2014
#

import os, re, sys, shutil
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

from Sire.Move import *

import itertools
import operator
import io  # For compatibility with old open style for binary files
import copy
import csv


#
# Python dependencies
#
try:
    numpy = Sire.try_import("numpy")
except:
    pass

try:
    mdtraj = Sire.try_import("mdtraj")
except:
    pass

import numpy as np
#try:
#    import mdtraj
#except ImportError:
#    raise "Nautilus depends on a working install of the python module mdtraj. Please install mdtraj in your sire python."
#    sys.exit(-1)
#
#try:
#    import numpy as np
#except ImportError:
#    raise "Nautilus depends on a working install of the python module numpy. Please install numpy in your sire python"
#    sys.exit(-1)
#
# #### VARIABLES #######
#

grid_center_x = Parameter("grid_center_x", 0.0, """The x coordinate of the center of the grid.""")
grid_center_y = Parameter("grid_center_y", 0.0, """The y coordinate of the center of the grid.""")
grid_center_z = Parameter("grid_center_z", 0.0, """The z coordinate of the center of the grid.""")
grid_plus_x = Parameter("grid_plus_x", 99999.0, """The positive extent along the x coordinate.""")
grid_min_x = Parameter("grid_min_x", 99999.0, """The negative extent along the x coordinate.""")
grid_plus_y = Parameter("grid_plus_y", 99999.0, """The positive extent along the y coordinate.""")
grid_min_y = Parameter("grid_min_y", 99999.0, """The negative extent along the y coordinate.""")
grid_plus_z = Parameter("grid_plus_z", 99999.0, """The positive extent along the z coordinate.""")
grid_min_z = Parameter("grid_min_z", 99999.0, """The negative extent along the z coordinate.""")

topfile = Parameter("topfile", "SYSTEM.top", """Name of the topology file containing the system to be simulated.""")
crdfile = Parameter("crdfile", "SYSTEM.crd", """Name of the coordinate file containing the coordinates of the system to be simulated.""")
trajfile = Parameter("trajfile", "traj0000000001.dcd",
                     """Name of the topology file containing the system to be simulated.""")
startframe = Parameter("start_frame", 0, """The first frame to analyse.""")
endframe = Parameter("end_frame", 100000000, """The last frame to analyse.""")
cutoff = Parameter("cutoff", 10 * angstrom, """The cutoff distance to use for the non-bonded interactions.""")
rfdielectric = Parameter("rfdielectric", 78.3, """Dielectric constant to use with the reaction field cutoff method.""")
watermodel = Parameter("watermodel", "TIP4PEW-SireOpenMM",
                       """The cell theory parameterisation to use for bulk water.""")
cell_dir = Parameter("cell_dir", "cell", """The output folder for generated cell files.""")
benchmark = Parameter("benchmark", False, """Whether or not to benchmark the Nautilus scripts.""")

# ######## Parameters specific for grid-from-cell

grid_dir = Parameter("grid_dir", "grid", """The output folder for generated grid files.""")
cell_interval = Parameter("cell_interval", (1000, 22000), """The grid will be averaged using all the cell files in the time interval below. Set to -1 and 1e10 if you want to use all data. Note the interval is in frame numbers""")
grid_step = Parameter("grid_step", 1.0, """The grid step defines the grid density in x/y/z in Angstroms""")
grid_count_cutoff = Parameter("grid_count_cutoff", 0,
                              """Grid points with less than the count cutoff will be discarded""")
temperature = Parameter("temperature", 298, """Temperature in Kelvin""")
frequencyupdate = Parameter("frequencyupdate", 1000,
                            """running average update for entire computed thermodynamics of the grid""")
Overwrite = Parameter("Overwrite", True, """Overwrite an existing grid folder """)

######### Parameters specific for regionproperties
regionfile = Parameter("regionfile", "all.region", """Region file which specifies grid points to be averaged """)
gridforces = Parameter("gridforces", "grid.forces", """Grid.forces file which specifies average parameters of each grid point""")

######## Parameters specific for cluster.py
sites_dir = Parameter("sites_dir", "sites", """The output folder for the clustered sites""")
lowt = Parameter("lowt", 1.5, """The density threshold to terminate clustering""")
neighcut = Parameter("neighcut", 1.5, """The maximum distance (in Angstroms) between grid points to consider them neighbors""")

####### Parameters specific for averaging and subtracting grids
gridf = Parameter("gridf", "gridf.dx", """Grid dx file f to be subtracted""")
gridl = Parameter("gridl", "gridl.dx", """Grid dx file l used to subtract""")
diffdx = Parameter("diffdx", "diff.dx", """Difference grid dx file """)
avgdx = Parameter("avgdx", "avg.dx", """Average grid dx file """)
avggridfiles = Parameter("avggridfiles", [], """List of grid dx files to be averaged""")


extraacceptors = Parameter("extra_acceptors",[],"""List of extra atom types to consider as hydrogen bond acceptors.""")
extrapolarH = Parameter("extra_polar_hydrogens",[],"""List of extra atom types to consider as polar hydrogens.""")

##### CONSTANTS #########
combining_rules = "arithmetic"
# Name of water residues in input
Water_residue_names = ["T4P", "T3P", "Wat", "HOH", "WAT"]
# Amber/GAFF atom types classified as donors or acceptors
polarH = ["HW", "H", "HO"]
acceptors = ["OW", "O", "OH", "Cl", "N", "N3", "N2", "Na"]
bulk = {"TIP4PEW-SireOpenMM": 4.88,
        "TIP4PEW-RH": 4.67,
	"TIP3P-SireOpenMM": 5.02}

# The cutoff to compute coordination numbers of water to other waters
CUTOFF_CN_WAT = 3.4
# The curoff to use for coordination number of water to solute acceptors
CUTOFF_CN_SOL = 3.4  # This is for Cl, but do we need different numbers for each solute acceptor?
# GRID BUFFER. When parsing input, keep track of position of donors/acceptors that are within grid dimension + GRIDBUF
# This will cause an incorrect calculation of hydrogen-bonds and coordination numbers for waters that are at the edge of the extended
# grid, since they will miss some neighbors. It is important that GRIDBUF is "big enough" then.
GRIDBUFF = 5.0
NONPOLARFLAG = 0
POLARHFLAG = 1
ACCEPTORFLAG = 2
# Atom with masses ( in g.mol-1) above this are considered "HEAVY"
LIGHT = 1.10
CUTOFF_CN_WAT2 = CUTOFF_CN_WAT ** 2
CUTOFF_CN_SOL2 = CUTOFF_CN_SOL ** 2

#### constants used in grid-from-cell
#########################################
T = temperature.val
SMALL = 1.  #0.00000001
MININF = -99999999.0
NAVOGADRO = 6.022141e23
MASS_WATER = 18.01528  # g/mol
KCAL_TO_J = 4184
KB = 0.001987206500956023  # * 1000. * 4.184

## DICTIONARY OF BULK PROPERTIES PER WATER MODEL ###############
# For consistency, all should be derived using same cutoff/code.#
#################################################################
# TIP3P-SireOpenMM: 10 ns NPT, 10 Angstrom cutoff, reaction field, 10k snapshots
# RF seems to lower density and make
# TIP4P-RH From J. Chem. Phys. 126, 064504 2007
# Forces are in N x 10e-10
# Torques are in N m-1 x 10e-20
# Energies are in kcal/mol
# Densities are in kg/m3

waterProps = {"TIP4P-RH": {"BULKORI": 2.94,
                           "BULKFORCES": [1.627, 1.518, 1.284],
                           "BULKTORQUES": [1.157, 0.994, 1.371],
                           "BULKNRG": -9.92,
                           "BULKDENSITY": "X.XX"
                          },
              "TIP4PEW-SireOpenMM":
                          {"BULKORI": 3.305,
                           "BULKFORCES": [1.587, 1.735, 1.334],
                           "BULKTORQUES": [1.061, 1.194, 1.453],
                           "BULKNRG": -11.025,
                           "BULKDENSITY": 995.4
                          },
               "TIP3P-SireOpenMM":
                          { "BULKORI": 3.659,
                            "BULKFORCES": [ 1.377, 1.643, 1.105 ],
                            "BULKTORQUES": [ 0.985, 0.935, 1.283 ],
                            "BULKNRG": -9.520,
                            "BULKDENSITY":  982.2
                          }
}
###############################################


###### Internal functions ##########
def _createSystem(molecules):

    moleculeNumbers = molecules.molNums()
    moleculeList = []

    for moleculeNumber in moleculeNumbers:
        molecule = molecules.molecule(moleculeNumber).molecule()
        moleculeList.append(molecule)

    all = MoleculeGroup("all")

    for molecule in moleculeList[0:]:
        all.add(molecule)

    gridwater = MoleculeGroup("gridwater")
    otherwater = MoleculeGroup("otherwater")
    solutes = MoleculeGroup("solutes")

    # Add these groups to the System
    system = System()

    system.add(all)
    system.add(gridwater)
    system.add(otherwater)
    system.add(solutes)

    return system

def _initGrid( grid_center_x, grid_center_y, grid_center_z, grid_plus_x, grid_min_x, grid_plus_y, grid_min_y, grid_plus_z, grid_min_z ):

    grid = {}
    grid['center'] = ( grid_center_x, grid_center_y, grid_center_z )
    grid['step'] = ( grid_plus_x, grid_min_x, grid_plus_y, grid_min_y, grid_plus_z, grid_min_z )
    grid['min'] = ( grid['center'][0] - grid['step'][1], grid['center'][1] - grid['step'][3], grid['center'][2] - grid['step'][5] )
    grid['max'] = ( grid['center'][0] + grid['step'][0], grid['center'][1] + grid['step'][2], grid['center'][2] + grid['step'][4] )

    vol = ( grid['step'][0] + grid['step'][1] ) * ( grid['step'][2] + grid['step'][3] ) *  ( grid['step'][4] + grid['step'][5] )
    grid['volume'] = vol

    return grid

def _assignDonorsAcceptors( system ):

    # Update list of default polarH, acceptors with user passed data (if any)
    for atype in extrapolarH.val:
        polarH.append(atype)
    for atype in extraacceptors.val:
        acceptors.append(atype)

    print ("List of polarH types: %s " % polarH)
    print ("list of acceptor types: %s " % acceptors)

    atomsDAtype = {}

    # map molecule number with atomic numbers
    mol_to_atoms = {}

    molnums = system.molecules().molNums()
    for molnum in molnums:
        mol = system.molecule(molnum).molecule()

        if mol.residues()[0].name().value() in Water_residue_names:
            isawatmol = True
        else:
            isawatmol = False

        mol_to_atoms[ molnum.value() ] = []

        atoms = mol.atoms()
        for x in range(0, mol.nAtoms()):
            atom = atoms[x]
            at_num = atom.number().value()
            at_type = atom.property("ambertype")
            at_charge = atom.property("charge").value()

            mol_to_atoms[ molnum.value() ].append( at_num )
            #print at_type
            #sys.exit(-1)
            connected_nums = []
            hbondflag = NONPOLARFLAG

            if at_type in polarH:
                hbondflag = POLARHFLAG
                # Also ...if this is NOT a water molecule, we
                # need to know to which heavy atoms the hydrogen is bonded
                # Annoyingly amber potentially define bonds between hydrogens so ignore H
                #if not isawatmol:
                if mol.nAtoms() > 1:
                    connectivity = mol.property("connectivity")
                    connected_atoms_indices = connectivity.connectionsTo(atom.number())
                    for index in connected_atoms_indices:
                        c_at = mol.select(index)
                        if c_at.property("mass").value() > LIGHT:
                            connected_nums.append( c_at.number().value() )
                else:
                   connected_nums.append( -1 )

            elif at_type in acceptors:
                hbondflag = ACCEPTORFLAG
                doconnect = True
                # Here needs to know to which polarH hydrogens the acceptor is bonded to
                try:
                    connectivity = mol.property("connectivity")
                except UserWarning:
                    doconnect = False
                if doconnect:
                    connected_atoms_indices = connectivity.connectionsTo(atom.number())
                    for index in connected_atoms_indices:
                        c_at = mol.select(index)
                        c_at_type = c_at.property("ambertype")
                        if c_at_type in polarH:
                            connected_nums.append( c_at.number().value() )

            atomsDAtype[at_num] = ( hbondflag, at_charge, connected_nums, molnum.value(), isawatmol)

    return atomsDAtype, mol_to_atoms

def _updateSystemfromDCDgrid(system, frame_xyz, cell_lengths, cell_angles, atomsDAtype, grid):

    dcd_coordinates = frame_xyz[0]

    dcd_box_x = cell_lengths[0][0].tolist()
    dcd_box_y = cell_lengths[0][1].tolist()
    dcd_box_z = cell_lengths[0][2].tolist()

    dcd_natoms = len(dcd_coordinates)

    newmols_coords = {}

    dcd_index = 0
    mol_index = 0

    molnums = system.molNums()
    molnums.sort()

    # Dictionnary of coordinates for polarH atoms
    polarH_coords = {}
    # Dictionnary of coordinates for acceptor atoms
    acceptor_coords = {}
    # alternate coords
    alternate_acceptor_coords = {}
    tb = time.time()
    for molnum in molnums:
        mol = system.molecule(molnum).molecule()
        molatoms = mol.atoms()
        molnatoms = mol.nAtoms()
        # Create an empty coord group using molecule so we get the correct layout
        newmol_coords = AtomCoords( mol.property("coordinates") )
        #print "Processing molecule %s natoms %s " % ( molnum.value(), molnatoms )
        for x in range(0,molnatoms):
            tmparray = dcd_coordinates[dcd_index]
            atom_coord = Vector( tmparray[0].tolist() , tmparray[1].tolist() , tmparray[2].tolist() )
            atom = molatoms[x]
            cgatomidx = atom.cgAtomIdx()
            newmol_coords.set( cgatomidx, atom_coord)
            dcd_index += 1
        newmols_coords[molnum] = newmol_coords
        mol_index += 1

    if dcd_natoms != dcd_index:
        print ("The number of atoms in the system is not equal to the number of atoms in the DCD file ! Aborting.")
        sys.exit(-1)


    grid_box_min = Vector(grid['min'][0], grid['min'][1], grid['min'][2] )
    grid_box_max = Vector(grid['max'][0], grid['max'][1], grid['max'][2] )

    grid_buffer_box_min = Vector(grid['min'][0]-GRIDBUFF, grid['min'][1]-GRIDBUFF, grid['min'][2]-GRIDBUFF )
    grid_buffer_box_max = Vector(grid['max'][0]+GRIDBUFF, grid['max'][1]+GRIDBUFF, grid['max'][2]+GRIDBUFF )

    new_gridwater = MoleculeGroup("gridwater")
    new_otherwater = MoleculeGroup("otherwater")
    new_solutes = MoleculeGroup("solutes")

    changedmols = MoleculeGroup("changedmols")

    mol_index = 0
    for molnum in molnums:
        mol = system.molecule(molnum).molecule()
        newmol_coords = newmols_coords[molnum]
        mol = mol.edit().setProperty("coordinates", newmol_coords).commit()
        changedmols.add(mol)

        if mol.residues()[0].name() == ResName('WAT'):
            # If water. check if oxygen atom within grid
            # Assumes oxygen first atom
            o_at = mol.select(  AtomName("O") ).property("coordinates")
            #
            if not ( o_at.x() < grid_box_min.x() or
                     o_at.y() < grid_box_min.y() or
                     o_at.z() < grid_box_min.z() or
                     o_at.x() > grid_box_max.x() or
                     o_at.y() > grid_box_max.y() or
                     o_at.z() > grid_box_max.z() ):
                new_gridwater.add(mol)
            else:
                new_otherwater.add(mol)
        else:
            new_solutes.add(mol)

        # Update dictionnaries of donors/acceptors/alternate_acceptors and gridwater group based on new cartesian coordinates
        molatoms = mol.atoms()
        molnatoms = mol.nAtoms()

        for x in range(0,molnatoms):
            atom = molatoms[x]
            at_num = atom.number().value()
            atom_coord = atom.property("coordinates")

            # Complicated because grid_buffer could extend beyond the dcd box
            # This code assumes the origin of the dcd box is 0,0,0
            within_x = False
            within_y = False
            within_z = False

            if atom_coord.x() < grid_buffer_box_max.x() and atom_coord.x() > grid_buffer_box_min.x():
               within_x = True
            # Also check wrapping cases
            if grid_buffer_box_max.x() > dcd_box_x:
               if atom_coord.x() < ( grid_buffer_box_max.x() - dcd_box_x ):
                   within_x = True
            if grid_buffer_box_min.x() < 0:
               if atom_coord.x() > ( dcd_box_x + grid_buffer_box_min.x() ):
                   within_x = True

            if atom_coord.y() < grid_buffer_box_max.y() and atom_coord.y() > grid_buffer_box_min.y():
               within_y = True
            # Also check wrapping cases
            if grid_buffer_box_max.y() > dcd_box_y:
               if atom_coord.y() < ( grid_buffer_box_max.y() - dcd_box_y ):
                   within_y = True
            if grid_buffer_box_min.y() < 0:
               if atom_coord.y() > ( dcd_box_y + grid_buffer_box_min.y() ):
                   within_y = True

            if atom_coord.z() < grid_buffer_box_max.z() and atom_coord.z() > grid_buffer_box_min.z():
               within_z = True
            # Also check wrapping cases
            if grid_buffer_box_max.z() > dcd_box_z:
               if atom_coord.z() < ( grid_buffer_box_max.z() - dcd_box_z ):
                   within_z = True
            if grid_buffer_box_min.z() < 0:
               if atom_coord.z() > ( dcd_box_z + grid_buffer_box_min.z() ):
                   within_z = True

            if not (within_x and within_y and within_z):
                continue

            # Save coordinates for later h-bond analysis
            at_num = atom.number().value()
            if atomsDAtype[at_num][0] == POLARHFLAG:
                polarH_coords[at_num] = atom_coord
            elif atomsDAtype[at_num][0] == ACCEPTORFLAG:
                acceptor_coords[at_num] = atom_coord
            # Ugly, to improve
            if atom.name().value() == "EPW":
                # The 'O' atom is 3 at_num before
                alternate_acceptor_coords[ at_num - 3 ] = (at_num, atom_coord )
        mol_index +=1
        #sys.exit(-1)
    # Update everything
    system.update(changedmols)
    system.removeAllMoleculeGroups()
    system.removeAllForceFields()
    new_all = MoleculeGroup("all", changedmols)
    system.add(new_all)
    system.add(new_gridwater)
    system.add(new_otherwater)
    system.add(new_solutes)

    #if  dcd_traj.trajectory.periodic:
    if 1:# HOW TO DETECT NON PERIODIC SYSTEMS ?
        space = PeriodicBox(Vector( dcd_box_x, dcd_box_y, dcd_box_z ) )
    else:
        space = Cartesian()

    system = _setupForcefields(system, space)
    return system, polarH_coords, acceptor_coords, alternate_acceptor_coords

def _setupForcefields(system, space):

    all = system[ MGName("all") ]
    gridwater = system[ MGName("gridwater") ]
    otherwater = system[MGName("otherwater")]
    solutes = system[MGName("solutes")]

    # Intermolecular energy for entire system
    gridwaterff = InterCLJFF("gridwaterff")
    gridwaterff.add(gridwater)

    gridwater_otherwaterff = InterGroupCLJFF("gridwater:otherwaterff")
    gridwater_otherwaterff.add(gridwater, MGIdx(0) )
    gridwater_otherwaterff.add(otherwater, MGIdx(1) )

    gridwater_solutesff = InterGroupCLJFF("gridwater:solutesff")
    gridwater_solutesff.add(gridwater, MGIdx(0) )
    gridwater_solutesff.add(solutes, MGIdx(1) )

    forcefields = [ gridwaterff, gridwater_otherwaterff, gridwater_solutesff ]

    for forcefield in forcefields:
        system.add(forcefield)

    system.setProperty( "space", space )
    system.setProperty( "switchingFunction",
                        CHARMMSwitchingFunction(cutoff.val) )
    system.setProperty( "useReactionField", VariantProperty(True) )
    system.setProperty( "reactionFieldDielectric", VariantProperty(rfdielectric.val) )
    system.setProperty( "combiningRules", VariantProperty(combining_rules) )

    gridwater_nrg =  gridwaterff.components().total() + gridwater_otherwaterff.components().total() +\
                     gridwater_solutesff.components().total()

    gridwater_solutes_nrg = gridwater_solutesff.components().total()

    e_gridwater = Symbol("E_{gridwater}")

    e_gridwater_solutes = Symbol("E_{gridwater:solutes}")

    system.setComponent( e_gridwater, gridwater_nrg )

    system.setComponent( e_gridwater_solutes, gridwater_solutes_nrg )

    return system

def _defineFSW(space, acceptor_coords, atomsDAtype, grid ):

    cutoff_cn2 = CUTOFF_CN_SOL2 # Angstroms

    acceptor_atoms = list(acceptor_coords.keys())

    nacceptors = len(acceptor_atoms)

    waters_firstshell = {}

    # Loop over acceptor_atoms and split in solute/waters
    water_acceptors = []
    solute_acceptors = []

    nwat = 0
    nsol = 0
    for i in range(0,nacceptors):
        iat = acceptor_atoms[i]
        iat_inwat = atomsDAtype[iat][4]
        if iat_inwat:
            water_acceptors.append(iat)
            nwat += 1
        else:
            solute_acceptors.append(iat)
            nsol += 1

    water_coords = np.zeros([nwat,3])
    for i in range(0,nwat):
        iat = water_acceptors[i]
        water_coords[i] = [ acceptor_coords[iat][0], acceptor_coords[iat][1], acceptor_coords[iat][2] ]

    solute_coords = np.zeros([nsol,3])
    for i in range(0,nsol):
        iat = solute_acceptors[i]
        solute_coords[i] = [ acceptor_coords[iat][0], acceptor_coords[iat][1], acceptor_coords[iat][2] ]

    dimensions = np.array( [space.dimensions().x(), space.dimensions().y(), space.dimensions().z()] )

    fsw = np.zeros(nwat)

    for i in range(0,nsol):
        point = solute_coords[i]

        d2 = _npdistance2(point, water_coords, dimensions)

        contact = np.where(d2 < cutoff_cn2, 1, 0)
        fsw += contact

    fsw_indices = fsw.nonzero()[0]
    for idx in fsw_indices:
        iat = water_acceptors[idx]
        iatmolnum = atomsDAtype[iat][3]
        waters_firstshell[iatmolnum] = True

    return waters_firstshell

def _assignHbonds( space, polarH_coords, acceptor_coords, alternate_acceptor_coords, atomsDAtype, grid):
    # the key-value mark a h-bond between a polarH and an acceptor
    h_bonds = {}

    polarH_atoms = list(polarH_coords.keys())
    polarH_atoms.sort()
    ND = len(polarH_atoms)

    acceptor_atoms = list(acceptor_coords.keys())
    acceptor_atoms.sort()
    NA = len(acceptor_atoms)

    for polarH_at in polarH_atoms:
        h_bonds[polarH_at] = []

    acoords = np.zeros([NA,3])
    ac = np.zeros(NA)

    dimensions = np.array( [space.dimensions().x(), space.dimensions().y(), space.dimensions().z()] )

    for i in range(0,NA):
        acceptor = acceptor_atoms[i]
        h_bonds[acceptor] = []
        try:
            A_coords = alternate_acceptor_coords[acceptor][1]
            A_charge = atomsDAtype[ alternate_acceptor_coords[acceptor][0] ][1]
        except KeyError:
            A_coords = acceptor_coords[acceptor]
            A_charge = atomsDAtype[acceptor][1]
        acoords[i] = [ A_coords[0], A_coords[1], A_coords[2] ]
        ac[i] = A_charge

    for polarH_at in polarH_atoms:
        H_coords = polarH_coords[polarH_at]
        point = np.array( [ H_coords[0], H_coords[1], H_coords[2] ] )

        connected_acceptors = atomsDAtype[polarH_at][2]

        hforce = np.empty(NA)
        d2 = _npdistance2(point, acoords, dimensions)

        hforce = ac / d2
        # Have to mask acceptor(s) bonded to polarH_at
        for acc in connected_acceptors:
            try:
                acceptor_index = acceptor_atoms.index(acc)
                hforce[acceptor_index] = +99
            except ValueError:
                # This could happen if acceptor was not in the buffered grid
                continue

        # Then get index from hforce.argmax()
        hmax = hforce.argmin()
        hbonded = acceptor_atoms[hmax]
        h_bonds[polarH_at].append(hbonded)
        h_bonds[hbonded].append(polarH_at)

    return h_bonds

def _coordinationNumbers( space, acceptor_coords, h_bonds, atomsDAtype, firstshellwaters, grid):

    grid_minx = grid['min'][0]
    grid_miny = grid['min'][1]
    grid_minz = grid['min'][2]
    grid_maxx = grid['max'][0]
    grid_maxy = grid['max'][1]
    grid_maxz = grid['max'][2]

    #print (grid_minx, grid_min_y, grid_minz, grid_maxx, grid_maxy, grid_maxz)

    water_cn = {}
    acceptor_atoms = list(acceptor_coords.keys())
    #JM TEST
    #acceptor_atoms.sort()
    #
    nacceptors = len(acceptor_atoms)

    acoords = np.zeros([nacceptors,3])
    for i in range(0,nacceptors):
        iat = acceptor_atoms[i]
        acoords[i] = [ acceptor_coords[iat][0], acceptor_coords[iat][1], acceptor_coords[iat][2] ]

    dimensions = np.array( [space.dimensions().x(), space.dimensions().y(), space.dimensions().z()] )

    for i in range(0,nacceptors-1):
        point = acoords[i]
        iat = acceptor_atoms[i]
        iat_inwat = atomsDAtype[iat][4]
        #print (point)
        #print (iat)
        #print (iat_inwat)

        iat_ingrid = True
        if ( acoords[i][0] < grid_minx or
             acoords[i][1] < grid_miny or
             acoords[i][2] < grid_minz or
             acoords[i][0] > grid_maxx or
             acoords[i][1] > grid_maxy or
             acoords[i][2] > grid_maxz ):
            # iat is outside of the grid...
            iat_ingrid = False

        d2 = _npdistance2( point, acoords[i+1:], dimensions)
        # Select acceptors within cutoff
        check = np.where(d2 < CUTOFF_CN_WAT2, 1, 0)
        coordinated = check.nonzero()[0]
        #print (coordinated)
        # Need to consider case where there are no coordinated particles
        # but particle is in grid then water_cn has to be initialised

        imolnum = atomsDAtype[iat][3]
        #print (imolnum)

        #if iat_ingrid and not water_cn.has_key(imolnum):
        if iat_ingrid and (not imolnum in water_cn):
            water_cn[imolnum] = {}
            water_cn[imolnum]['fsw'] = []
            water_cn[imolnum]['bulk'] = []
            water_cn[imolnum]['solute'] = []

        for j in coordinated:
            shift = i + j + 1
            jat = acceptor_atoms[ shift ]
            jat_inwat = atomsDAtype[jat][4]
            # if iat and jat are in a solute. skip
            if ( (not iat_inwat) and (not jat_inwat) ):
                continue

            # if iat and jat are both outside of the grid, skip
            # !!! NO PBC CHECK
            jat_ingrid = True

            if ( acoords[shift][0] < grid_minx or
                 acoords[shift][1] < grid_miny or
                 acoords[shift][2] < grid_minz or
                 acoords[shift][0] > grid_maxx or
                 acoords[shift][1] > grid_maxy or
                 acoords[shift][2] > grid_maxz ):
                # jat is also outside of the grid...
                jat_ingrid = False

            if ( ( not iat_ingrid ) and ( not jat_ingrid ) ):
                continue

            # if iat in water and jat in a solute or
            # if jat in water and iat in a solute
            #     Funny rule is...if acceptor is not in a water molecule
            #     AND if is bonded to a polarH that is h-bonded to iat --> skip

            if ( iat_inwat and not jat_inwat):
                jat_connected_polarHs = atomsDAtype[jat][2]
                for jat_polarH in jat_connected_polarHs:
                    if h_bonds[jat_polarH] == iat:
                        continue

            if ( jat_inwat and not iat_inwat):
                iat_connected_polarHs = atomsDAtype[jat][2]
                for iat_polarH in iat_connected_polarHs:
                    if h_bonds[iat_polarH] == jat:
                        continue

            jmolnum = atomsDAtype[jat][3]

            if iat_inwat and iat_ingrid:
                if jat_inwat:
                    #if firstshellwaters.has_key(jmolnum):
                    if jmolnum in firstshellwaters:
                        cntype = 'fsw'
                    else:
                        cntype = 'bulk'
                else:
                    cntype = 'solute'

                #if not water_cn.has_key(imolnum):
                if not imolnum in water_cn:
                    water_cn[imolnum] = {}
                    water_cn[imolnum]['fsw'] = []
                    water_cn[imolnum]['bulk'] = []
                    water_cn[imolnum]['solute'] = []
                water_cn[imolnum][cntype].append( jat )

            if jat_inwat and jat_ingrid:
                if iat_inwat:
                    #if firstshellwaters.has_key(imolnum):
                    if imolnum in firstshellwaters:
                        cntype = 'fsw'
                    else:
                        cntype = 'bulk'
                else:
                    cntype = 'solute'


                #if not water_cn.has_key(jmolnum):
                if not jmolnum in water_cn:
                    water_cn[jmolnum] = {}
                    water_cn[jmolnum]['fsw'] = []
                    water_cn[jmolnum]['bulk'] = []
                    water_cn[jmolnum]['solute'] = []
                water_cn[jmolnum][cntype].append( iat )
        #print (water_cn)
        #print ("STOP")
        #sys.exit(-1)
    #print (len(water_cn))
    #sys.exit(-1)
    return water_cn

def _orientationalNumbers( water_cn, firstshellwaters, h_bonds, mol_to_atoms, atomsDAtype, BULKCN ):
    water_ori = {}

    # For "normal" water
    #   just sum neighbors, and then convert to ori
    # For each first shell water
    # We know number of solutes number of first shell and number of bulk from water_cn...
    # We need to know if we are donating a h-bond to a solute atom, to another first shell water, to a bulk water (pwl_HB)
    # Then we can compute the effective number and then the orientational number

    waters = list(water_cn.keys())
    #print ("waters", waters)
    for wat in waters:
        #if firstshellwaters.has_key(wat):
        if wat in firstshellwaters:
            # First shell water treatment
            molatoms = mol_to_atoms[wat]

            o_atomnum = None
            h_atomnums = []

            for atnum in molatoms:
                if atomsDAtype[atnum][0] == POLARHFLAG:
                    h_atomnums.append(atnum)
                elif atomsDAtype[atnum][0] == ACCEPTORFLAG:
                    o_atomnum = atnum

            hbonded = []
            # the atomic numbers of the atoms that donate a h-bond to wat oxygen
            h_accepteds = h_bonds[o_atomnum]
            # now find to which heavy atom they are bonded
            for h in h_accepteds:
                hbonded += atomsDAtype[h][2]
            # the atomic numbers of the atoms that accepte a h-bond from wat
            hbonded +=  h_bonds[ h_atomnums[0] ] + h_bonds[ h_atomnums[1] ]

            # For each type of cn member, find out how many are currently h-bonded to that water
            nw_solute = len(water_cn[wat]['solute'])
            nw_solute_hb = 0
            phb_solute = 0
            # for each solute, check if accepting a h-bond from that wat
            for sol in water_cn[wat]['solute']:
                if sol in hbonded:
                    nw_solute_hb += 1
            if nw_solute > 0:
                phb_solute = float(nw_solute_hb) / nw_solute

            # now the fsw
            nw_fsw = len(water_cn[wat]['fsw'])
            nw_fsw_hb = 0
            phb_fsw = 0
            for fsw in water_cn[wat]['fsw']:
                if fsw in hbonded:
                    nw_fsw_hb += 1
            if nw_fsw > 0:
                phb_fsw = float(nw_fsw_hb)/nw_fsw
            # now the bulk water
            nw_bulk = len(water_cn[wat]['bulk'])
            nw_bulk_hb = 0
            phb_bulk = 0
            for bulk in water_cn[wat]['bulk']:
                if bulk in hbonded:
                    #print "hbonded with bulk atom in coordination sphere"
                    nw_bulk_hb += 1
            if nw_bulk > 0:
                phb_bulk = float(nw_bulk_hb)/nw_bulk

            phb_max = max ( max ( phb_solute, phb_fsw ) , phb_bulk )
            if phb_max < small:
                neff = 0
            else:
                neff = nw_solute_hb / phb_max + nw_fsw_hb / phb_max + nw_bulk_hb / phb_max
            # And ori
            nw = BULKCN
            ori = ( ( neff * ( neff-1. ) ) / 2. ) * ( ( nw - 2. ) / nw )**(2-phb_solute)
        else:
            # bulk treatment
            nw = len(water_cn[wat]['bulk']) + len(water_cn[wat]['solute']) + len(water_cn[wat]['fsw'])
            #print nw
            if nw < 1:
                ori = 0.0
            else:
                ori = ( (nw * ( nw-1.) ) / 2. ) * ( ( nw - 2. ) / nw ) ** 2
        water_ori[wat] = ori
    #print water_ori
    #sys.exit(-1)

    return water_ori

def _npdistance2(x0, x1, dimensions):
    delta = np.abs(x0 - x1)
    delta = np.where(delta > 0.5 * dimensions, dimensions - delta, delta)
    return (delta ** 2).sum(axis=-1)

def _writeCellDataBin(cell_dir, data, frame_number, volume):

    keys = list(data.keys())

    keys.sort()

    bincell = struct.pack('<1i1f',frame_number,volume)

    #### making a string cell file for frame0.dat
    if frame_number == 0:
        file = open("cell0.dat", "w")

    for k in keys:
        bincell +=struct.pack('<2i12f', k, data[k]['waterflag'],data[k]['energies'][0], data[k]['solenergies'][0], data[k]['forces'][0],data[k]['forces'][1],data[k]['forces'][2],data[k]['torques'][0],data[k]['torques'][1],data[k]['torques'][2], data[k]['o_coords'][0],data[k]['o_coords'][1],data[k]['o_coords'][2], data[k]['ori'])
        if frame_number == 0:
            strcell = "%8d %8d %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n" % (k, data[k]['waterflag'],data[k]['energies'][0], data[k]['forces'][0],data[k]['forces'][1],data[k]['forces'][2],data[k]['torques'][0],data[k]['torques'][1],data[k]['torques'][2], data[k]['o_coords'][0],data[k]['o_coords'][1],data[k]['o_coords'][2], data[k]['ori'])
        # strcell = "%8d %8.5f %8.5f %s %8.5f " % ( k, data[k]['o_coords'][1], len(data[k]['cn']['bulk'])+len(data[k]['cn']['fsw'])+len(data[k]['cn']['solute']), data[k]['cn'], data[k]['ori'] )
            file.write(strcell)
            #print (strcell)
    #sys.exit(-1)
    stream = io.open("%s/cell_results_%016d.bin" % (cell_dir,frame_number) ,"wb")
    stream.write(bincell)
    stream.close()
    if frame_number == 0:
        file.close()

#### functions specific for cell2grid

def _initGrid_c2g(grid_step):
    grid = {}
    grid['center'] = ( grid_center_x.val, grid_center_y.val, grid_center_z.val )
    grid['size'] = ( grid_plus_x.val, grid_min_x.val, grid_plus_y.val, grid_min_y.val, grid_plus_z.val, grid_min_z.val )
    grid['min'] = (
        grid['center'][0] - grid['size'][1], grid['center'][1] - grid['size'][3], grid['center'][2] - grid['size'][5] )
    grid['max'] = (
        grid['center'][0] + grid['size'][0], grid['center'][1] + grid['size'][2], grid['center'][2] + grid['size'][4] )

    print("# The grid coordinates range from (%8.3f, %8.3f, %8.3f) to (%8.3f, %8.3f, %8.3f) " % (
        grid['min'][0], grid['min'][1], grid['min'][2], grid['max'][0], grid['max'][1], grid['max'][2]))

    nx = int(round(( grid['max'][0] - grid['min'][0] ) / grid_step))
    ny = int(round(( grid['max'][1] - grid['min'][1] ) / grid_step))
    nz = int(round(( grid['max'][2] - grid['min'][2] ) / grid_step))

    origin_x = grid['min'][0] + grid_step / 2.0
    origin_y = grid['min'][1] + grid_step / 2.0
    origin_z = grid['min'][2] + grid_step / 2.0

    grid = {}
    grid['origin'] = [origin_x, origin_y, origin_z]
    grid['dimensions'] = ( nx, ny, nz )
    grid['step'] = grid_step
    grid['points'] = {}
    grid['npoints'] = 0
    grid['nsnapshots'] = 0

    grid_counter = 0

    for x in range(0, nx):
        grid['points'][x] = {}
        for y in range(0, ny):
            grid['points'][x][y] = {}
            for z in range(0, nz):
                grid['points'][x][y][z] = {}

                grid['points'][x][y][z]['avgnrg'] = [0.0, 0.0, 0.0]
                grid['points'][x][y][z]['avgsolnrg'] = [0.0, 0.0, 0.0]
                grid['points'][x][y][z]['avgforces'] = [0.0, 0.0, 0.0]
                grid['points'][x][y][z]['avgtorques'] = [0.0, 0.0, 0.0]
                grid['points'][x][y][z]['avgori'] = 0
                grid['points'][x][y][z]['count'] = 0

                grid_counter += 1

    print("# There are %s grid points ( dimensions %s %s %s step %s ) " % ( grid_counter, nx, ny, nz, grid_step))

    grid['npoints'] = grid_counter

    return grid


def _getCellData(cell_dir, cell_interval):
    cell_files = []

    files = os.listdir(cell_dir)

    for f in files:
        if f.startswith("cell") and f.endswith(".bin"):
            elems = f.split("_")
            num = elems[2].strip(".bin")
            num = float(num)
            if num >= cell_interval[0] and num <= cell_interval[1]:
                cell_files.append(os.path.join(cell_dir, f))

    cell_files.sort()

    return cell_files


def _loadCell(cell_file):
    waters_energies = [0.0, 0.0, 0.0]
    waters_solenergies = [0.0, 0.0, 0.0]
    waters_aforces = [0.0, 0.0, 0.0]
    waters_atorques = [0.0, 0.0, 0.0]
    waters_aori = 0

    nwaters = 0

    cell = {}
    cell['water'] = {}

    stream = open(cell_file, "rb")

    s = stream.read(8)
    framenumber, volume = struct.unpack("if", s)

    cell['frame'] = framenumber
    cell['volume'] = volume

    totnrg = 0.0
    s = stream.read(56)
    while ( len(s) == 56):
        resnum, watflag, cljnrg, cljsolnrg, fx, fy, fz, tx, ty, tz, ox, oy, oz, ori = struct.unpack('<2i12f', s)
        #print "<<<<", resnum, watflag, cljnrg, cljsolnrg, fx, fy, fz, tx, ty, tz, ox, oy, oz, ori, ">>>>"
        #sys.exit(-1)
        #if ori < 1 and watflag:
        #    #print "# low ori ( %s ) for water %s in cell file %s " % ( ori, resnum, cell_file)
        #    #ori = 1
        #    pass

        if watflag:
            cell['water'][resnum] = {}
            cell['water'][resnum]['energies'] = ( cljnrg, 0, 0 )
            cell['water'][resnum]['solenergies'] = ( cljsolnrg, 0, 0)
            cell['water'][resnum]['forces'] = ( fx, fy, fz )
            cell['water'][resnum]['torques'] = ( tx, ty, tz )
            cell['water'][resnum]['o_coords'] = ( ox, oy, oz )
            cell['water'][resnum]['ori'] = ori
            nwaters += 1

            waters_energies[0] += cljnrg
            waters_solenergies[0] += cljsolnrg

            waters_aforces[0] += fx
            waters_aforces[1] += fy
            waters_aforces[2] += fz

            waters_atorques[0] += tx
            waters_atorques[1] += ty
            waters_atorques[2] += tz

            waters_aori += ori

        s = stream.read(56)

    s = stream.read(56)
    if len(s) != 0:
        print("Cell file %s seems corrupted. Aborting." % cell_file)
        sys.exit(-1)

    if nwaters != 0:

        for i in range(0, 3):
            waters_energies[i] /= nwaters
            waters_solenergies[i] /= nwaters
            waters_aforces[i] /= nwaters
            waters_atorques[i] /= nwaters

        waters_aori /= float(nwaters)

    #print (waters_energies, waters_solenergies, waters_aforces)
    #sys.exit(-1)

    cell['waters_energies'] = waters_energies
    cell['waters_solenergies'] = waters_solenergies
    cell['waters_aforces'] = waters_aforces
    cell['waters_atorques'] = waters_atorques
    cell['waters_aori'] = waters_aori
    cell['nwaters'] = nwaters

    #print (cell['waters_energies'][0])
    #sys.exit(-1)

    return cell


def _updateGridPropsCumulative(grid, cell):
    step = grid['step']

    nx, ny, nz = grid['dimensions']
    originx, originy, originz = grid['origin']

    waters = list(cell['water'].keys())

    for water in waters:
        ox = cell['water'][water]['o_coords'][0]
        oy = cell['water'][water]['o_coords'][1]
        oz = cell['water'][water]['o_coords'][2]

        # Convert cartesian into grid coordinates
        gridx = int(round(( ox - originx ) / step))
        gridy = int(round(( oy - originy ) / step))
        gridz = int(round(( oz - originz ) / step))

        # Ignore water molecules outside of grid
        if ( gridx >= nx or gridy >= ny or gridz >= nz or gridx < 0 or gridy < 0 or gridz < 0):
            #print "WARNING water %s outside of grid, it will be ignored " % (water)
            continue

        # We accumulate arithmetic mean of forces/torques/coordination numbers, and arithmetic mean of energies
        cumweight = grid['points'][gridx][gridy][gridz]['count'] + 1.0
        grid['points'][gridx][gridy][gridz]['count'] = cumweight
        newweight = 1.0

        bigratio = ( cumweight - newweight ) / cumweight
        smallratio = 1.0 - bigratio

        # Energies
        oldavgnrg = grid['points'][gridx][gridy][gridz]['avgnrg'][0]
        newnrg = cell['water'][water]['energies'][0]
        grid['points'][gridx][gridy][gridz]['avgnrg'][0] = oldavgnrg * bigratio + newnrg * smallratio

        oldavgsolnrg = grid['points'][gridx][gridy][gridz]['avgsolnrg'][0]
        newsolnrg = cell['water'][water]['solenergies'][0]
        grid['points'][gridx][gridy][gridz]['avgsolnrg'][0] = oldavgsolnrg * bigratio + newsolnrg * smallratio

        # Forces
        oldavgforces = grid['points'][gridx][gridy][gridz]['avgforces']
        newforces = [0.0, 0.0, 0.0]
        newavgforces = [0.0, 0.0, 0.0]
        for m in range(0, 3):
            newforces[m] = cell['water'][water]['forces'][m]
            newavgforces[m] = oldavgforces[m] * bigratio + newforces[m] * smallratio
        grid['points'][gridx][gridy][gridz]['avgforces'] = newavgforces

        # Torques
        oldavgtorques = grid['points'][gridx][gridy][gridz]['avgtorques']
        newtorques = [0.0, 0.0, 0.0]
        newavgtorques = [0.0, 0.0, 0.0]
        for m in range(0, 3):
            newtorques[m] = cell['water'][water]['torques'][m]
            newavgtorques[m] = oldavgtorques[m] * bigratio + newtorques[m] * smallratio
        grid['points'][gridx][gridy][gridz]['avgtorques'] = newavgtorques

        # Orientational numbers
        oldavgori = grid['points'][gridx][gridy][gridz]['avgori']
        newori = cell['water'][water]['ori']
        newavgori = oldavgori * bigratio + newori * smallratio
        grid['points'][gridx][gridy][gridz]['avgori'] = newavgori

    grid['nsnapshots'] += 1


def _updatePropsCumulative(water_props, cell):
    #
    # We accumulate the arithmetic average of the forces, torques
    #

    if len(water_props) == 0:
        water_props['density'] = 0.0
        water_props['forces'] = [0.0, 0.0, 0.0]
        water_props['torques'] = [0.0, 0.0, 0.0]
        water_props['ori'] = 0.0

        water_props['nrg'] = 0.0
        water_props['solnrg'] = 0.0
        water_props['cum_weight'] = 0.0
        water_props['avg_w'] = 0.0

    newweight = cell['nwaters']

    water_props['cum_weight'] += newweight

    cumweight = water_props['cum_weight']

    bigratio = ( cumweight - newweight ) / cumweight
    smallratio = 1.0 - bigratio

    oldavgw = water_props['avg_w']
    neww = cell['nwaters']
    newavgw = oldavgw * bigratio + neww * smallratio
    water_props['avg_w'] = newavgw

    oldavgnrg = water_props['nrg']
    newnrg = cell['waters_energies'][0]
    newavgnrg = oldavgnrg * bigratio + newnrg * smallratio
    water_props['nrg'] = newavgnrg

    oldavgsolnrg = water_props['solnrg']
    newsolnrg = cell['waters_solenergies'][0]
    newavgsolnrg = oldavgsolnrg * bigratio + newsolnrg * smallratio
    water_props['solnrg'] = newavgsolnrg


    # Arithmetic mean forces/torques/ori/density
    oldavgforces = water_props['forces']
    newavgforces = [0.0, 0.0, 0.0]

    # average of forces and torques
    newforces = [0.0, 0.0, 0.0]
    # print water_num
    for m in range(0, 3):
        newforces[m] = cell['waters_aforces'][m]

    oldavgtorques = water_props['torques']
    newavgtorques = [0.0, 0.0, 0.0]
    newtorques = [0.0, 0.0, 0.0]
    for m in range(0, 3):
        newtorques[m] = cell['waters_atorques'][m]

    for m in range(0, 3):
        newavgforces[m] = oldavgforces[m] * bigratio + newforces[m] * smallratio
        newavgtorques[m] = oldavgtorques[m] * bigratio + newtorques[m] * smallratio

    water_props['forces'] = newavgforces
    water_props['torques'] = newavgtorques

    oldavgori = water_props['ori']
    newori = cell['waters_aori']
    newavgori = oldavgori * bigratio + newori * smallratio
    water_props['ori'] = newavgori

    # density..
    oldavgdensity = water_props['density']
    newdensity = _calcDensity(cell['volume'], cell['nwaters'])
    newavgdensity = oldavgdensity * bigratio + newdensity * smallratio
    water_props['density'] = newavgdensity
    #print (newavgdensity)
    #sys.exit(-1)

    return water_props


def _calcAverageProps(water_props, BULKORI, BULKFORCES, BULKTORQUES, BULKNRG, BULKDENSITY, BULKWATMOLPER_VOXELVOL):
    vibratio = 1.
    libratio = 1.
    oriratio = 1.

    DHWW = 0.0

    NW = water_props['avg_w']

    # average over all waters
    meanF = [0.0, 0.0, 0.0]
    meanT = [0.0, 0.0, 0.0]
    meanori = 0.0

    for m in range(0, 3):
        meanF[m] = water_props['forces'][m]
        meanT[m] = water_props['torques'][m]

        # !!!! The forces from the cell files are in kcal/mol/A
    # However, for consistency with the literature, the BULK forces are
    # in N/m x 1e-10 so we need to convert
    conv_factor = ( KCAL_TO_J / NAVOGADRO ) * 1e20
    for m in range(0, 3):
        meanF[m] = meanF[m] * conv_factor
        meanT[m] = meanT[m] * conv_factor

    for m in range(0, 3):
        if meanF[m] == 0:
            vibratio = 1.0
            break
        if meanT[m] == 0:
            libratio = 1.0
            break
        vibratio *= ( BULKFORCES[m] / meanF[m] )
        libratio *= ( BULKTORQUES[m] / meanT[m] )

    ori = water_props['ori']
    if ori < SMALL:
        ori = SMALL
    oriratio *= ( ori / BULKORI )

    #print "vibratio is %s " % vibratio
    #print "libratio is %s " % libratio
    #print "oriratio is %s " % oriratio

    DHSOL = NW * ( water_props['solnrg'] )

    DHSLV = NW * ( water_props['nrg'] - BULKNRG )

    meanwwnrg = water_props['nrg']

    meanori = ori

    DS_vib = -T * KB * NW * math.log(vibratio)

    DS_lib = -T * KB * NW * math.log(libratio)

    DS_ori = -T * KB * NW * math.log(oriratio)

    DSW = DS_vib + DS_lib + DS_ori

    DGW = DHSOL + DHSLV + DSW

    meandensity = water_props['density']
    #print (meandensity, "bah")
    #sys.exit(-1)

    return DS_vib, DS_lib, DS_ori, DSW, DHSOL, DHSLV, DGW, meanF, meanT, meanwwnrg, meanori, meandensity


def _calcPointProps(point, T, BULKORI, BULKFORCES, BULKTORQUES, BULKNRG, BULKDENSITY, BULKWATMOLPER_VOXELVOL):
    vibratio = 1.
    libratio = 1.
    oriratio = 1.

    NW = 1.0

    meanF = [0.0, 0.0, 0.0]
    meanT = [0.0, 0.0, 0.0]

    for m in range(0, 3):
        meanF[m] = point['avgforces'][m]
        meanT[m] = point['avgtorques'][m]

    # !!!! The forces from the cell files are in kcal/mol/A
    # However, for consistency with the literature, the BULK forces are
    # in N/m x 1e-10 so we need to convert
    conv_factor = ( KCAL_TO_J / NAVOGADRO ) * 1e20
    for m in range(0, 3):
        meanF[m] = meanF[m] * conv_factor
        meanT[m] = meanT[m] * conv_factor

    for m in range(0, 3):
        if meanF[m] == 0:
            vibratio = 1
            break
        if meanT[m] == 0:
            libratio = 1
            break
        vibratio *= ( BULKFORCES[m] / meanF[m] )
        libratio *= ( BULKTORQUES[m] / meanT[m] )

    ori = point['avgori']
    if (ori < SMALL):
        ori = SMALL
    oriratio *= ( ori / BULKORI )

    DHSOL = ( point['avgsolnrg'][0] ) * NW

    DHSLV = ( point['avgnrg'][0] - BULKNRG ) * NW

    DS_vib = -T * KB * math.log(vibratio) * NW
    DS_lib = -T * KB * math.log(libratio) * NW
    DS_ori = -T * KB * math.log(oriratio) * NW

    #density = calcDensity( volume,  point['count'] / nsnapshots )
    #
    #relative_density = density / BULKDENSITY

    #return relative_density, DS_vib, DS_lib, DS_ori, DH
    return DS_vib, DS_lib, DS_ori, DHSOL, DHSLV


def _calcDensity(volume, nwats):
    if nwats != 0:
        water_per_cubicangstrom = nwats / volume
        cubicdm_per_mol_water = (1 / water_per_cubicangstrom) * NAVOGADRO * 1e-27
        water_density = (1 / cubicdm_per_mol_water) * MASS_WATER  # g/L
    else:
        water_density = 0
    #print (water_density)
    #sys.exit(-1)
    return water_density


def _outputGrid(grid, grid_dir, grid_count_cutoff, BULKORI, BULKFORCES, BULKTORQUES, BULKNRG, BULKDENSITY,
                BULKWATMOLPER_VOXELVOL):
    #
    # Format a PDB style file with the cartesian coordinate of the grid points. Use columns like
    # occupancy or B-factors to write the entropies enthalpies etc...
    #

    titlepdbDG = '%s/grid_DG.pdb' % (grid_dir)
    filinpdbDG = open(titlepdbDG, 'w')

    titlepdbDH = '%s/grid_DH.pdb' % (grid_dir)
    filinpdbDH = open(titlepdbDH, 'w')

    titlepdbDHsol = '%s/grid_DHsolute.pdb' % (grid_dir)
    filinpdbDHsol = open(titlepdbDHsol, 'w')

    titlepdbDHslv = '%s/grid_DHsolvent.pdb' % (grid_dir)
    filinpdbDHslv = open(titlepdbDHslv, 'w')

    titlepdbDS = '%s/grid_minTDS.pdb' % (grid_dir)
    filinpdbDS = open(titlepdbDS, 'w')

    titlepdbDSvib = '%s/grid_minTDSvib.pdb' % (grid_dir)
    filinpdbDSvib = open(titlepdbDSvib, 'w')

    titlepdbDSlib = '%s/grid_minTDSlib.pdb' % (grid_dir)
    filinpdbDSlib = open(titlepdbDSlib, 'w')

    titlepdbDSori = '%s/grid_minTDSori.pdb' % (grid_dir)
    filinpdbDSori = open(titlepdbDSori, 'w')

    titlepdboccupancy = '%s/grid_density.pdb' % (grid_dir)
    filinpdboccupancy = open(titlepdboccupancy, 'w')

    nx, ny, nz = grid['dimensions']
    originx, originy, originz = grid['origin']
    step = grid['step']
    point_volume = step * step * step

    nsnapshots = grid['nsnapshots']

    index = 0

    for x in range(0, nx):
        for y in range(0, ny):
            for z in range(0, nz):
                count = grid['points'][x][y][z]['count']
                if count > grid_count_cutoff:
                    abs_density = _calcDensity(point_volume, count / nsnapshots)
                    density = abs_density / BULKDENSITY
                    minTDSvib, minTDSlib, minTDSori, DHSOL, DHSLV = _calcPointProps(grid['points'][x][y][z], T, BULKORI,
                                                                                    BULKFORCES, BULKTORQUES, BULKNRG,
                                                                                    BULKDENSITY, BULKWATMOLPER_VOXELVOL)
                else:
                    density = 0.0
                    DH = 0.0
                    minTDSvib = 0.0
                    minTDSlib = 0.0
                    minTDSori = 0.0
                    DHSOL = 0.0
                    DHSLV = 0.0

                index += 1

                minTDS = minTDSvib + minTDSlib + minTDSori
                DH = DHSOL + DHSLV
                DG = DH + minTDS

                grid['points'][x][y][z]['delta_G'] = DG
                grid['points'][x][y][z]['delta_H'] = DH
                grid['points'][x][y][z]['delta_H_sol'] = DHSOL
                grid['points'][x][y][z]['delta_H_slv'] = DHSLV
                grid['points'][x][y][z]['minT_delta_S'] = minTDS
                grid['points'][x][y][z]['minT_delta_S_vib'] = minTDSvib
                grid['points'][x][y][z]['minT_delta_S_lib'] = minTDSlib
                grid['points'][x][y][z]['minT_delta_S_ori'] = minTDSori
                grid['points'][x][y][z]['density'] = density

                cartx = originx + step * x
                carty = originy + step * y
                cartz = originz + step * z

                textDG = 'ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (
                    index, 'H', cartx, carty, cartz, '1.00', DG, density, 'H')
                textDH = 'ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (
                    index, 'H', cartx, carty, cartz, '1.00', DH, density, 'H')
                textDHsol = 'ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (
                    index, 'H', cartx, carty, cartz, '1.00', DHSOL, density, 'H')
                textDHslv = 'ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (
                    index, 'H', cartx, carty, cartz, '1.00', DHSLV, density, 'H')
                textDS = 'ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (
                    index, 'H', cartx, carty, cartz, '1.00', minTDS, density, 'H')
                textDSvib = 'ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (
                    index, 'H', cartx, carty, cartz, '1.00', minTDSvib, density, 'H')
                textDSlib = 'ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (
                    index, 'H', cartx, carty, cartz, '1.00', minTDSlib, density, 'H')
                textDSori = 'ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (
                    index, 'H', cartx, carty, cartz, '1.00', minTDSori, density, 'H')

                textocc = 'ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (
                    index, 'H', cartx, carty, cartz, '1.00', density, count, 'H')

                filinpdbDG.write(textDG + "\n")
                filinpdbDH.write(textDH + "\n")
                filinpdbDHsol.write(textDHsol + "\n")
                filinpdbDHslv.write(textDHslv + "\n")
                filinpdbDS.write(textDS + "\n")
                filinpdbDSvib.write(textDSvib + "\n")
                filinpdbDSlib.write(textDSlib + "\n")
                filinpdbDSori.write(textDSori + "\n")

                filinpdboccupancy.write(textocc + "\n")

    #(-1)

    filinpdbDG.close()
    filinpdbDH.close()
    filinpdbDHsol.close()
    filinpdbDHslv.close()
    filinpdbDS.close()
    filinpdbDSvib.close()
    filinpdbDSlib.close()
    filinpdbDSori.close()
    filinpdboccupancy.close()

    # Format a DX style file with the cartesian coordinate of the grid points. Use columns like
    # occupancy or B-factors to write the entropies enthalpies etc...
    #

    ### keep list of each data file for moe
    dg_moe = ""
    dh_moe = ""
    dhsolute_moe = ""
    dhsolvent_moe = ""
    mintds_moe = ""
    mintdsvib_moe = ""
    mintdslib_moe = ""
    mintdsori_moe = ""
    occupancy_moe = ""

    titledxDG = '%s/grid_DG.dx' % (grid_dir)
    filindxDG = open(titledxDG, 'w')

    titledxDH = '%s/grid_DH.dx' % (grid_dir)
    filindxDH = open(titledxDH, 'w')

    titledxDHsol = '%s/grid_DHsolute.dx' % (grid_dir)
    filindxDHsol = open(titledxDHsol, 'w')

    titledxDHslv = '%s/grid_DHsolvent.dx' % (grid_dir)
    filindxDHslv = open(titledxDHslv, 'w')

    titledxDS = '%s/grid_minTDS.dx' % (grid_dir)
    filindxDS = open(titledxDS, 'w')

    titledxDSvib = '%s/grid_minTDSvib.dx' % (grid_dir)
    filindxDSvib = open(titledxDSvib, 'w')

    titledxDSlib = '%s/grid_minTDSlib.dx' % (grid_dir)
    filindxDSlib = open(titledxDSlib, 'w')

    titledxDSori = '%s/grid_minTDSori.dx' % (grid_dir)
    filindxDSori = open(titledxDSori, 'w')

    titledxoccupancy = '%s/grid_density.dx' % (grid_dir)
    filindxoccupancy = open(titledxoccupancy, 'w')

    nx, ny, nz = grid['dimensions']
    originx, originy, originz = grid['origin']
    step = grid['step']

    index = 0

    ######FORMAT an appropriate DX file for both DH and DS ###############
    #####header of DX file ##############
    #object 1 class gridpositions counts 16 17 15
    line1 = 'object 1 class gridpositions counts %i %i %i\n' % (
        grid['dimensions'][0], grid['dimensions'][1], grid['dimensions'][2])
    #origin -1.0 -1.0  0.0
    line2 = 'origin %4.2f %4.2f %4.2f\n' % (grid['origin'][0], grid['origin'][1], grid['origin'][2])
    #delta 1 0 0
    line3 = 'delta %4.2f 0 0\n' % (grid['step'])
    #delta 0 1 0
    line4 = 'delta 0 %4.2f 0\n' % (grid['step'])
    #delta 0 0 1
    line5 = 'delta 0 0 %4.2f\n' % (grid['step'])
    #object 2 class gridconnections counts 16 17 15
    line6 = 'object 2 class gridconnections counts %4.2f %4.2f %4.2f\n' % (
        grid['dimensions'][0], grid['dimensions'][1], grid['dimensions'][2])
    #object 3 class array type double rank 0 items 4080 data follows
    line7 = 'object 3 class array type double rank 0 items %i data follows\n' % (
        grid['dimensions'][0] * grid['dimensions'][1] * grid['dimensions'][2])

    filindxDS.write(line1 + line2 + line3 + line4 + line5 + line6 + line7)
    filindxDH.write(line1 + line2 + line3 + line4 + line5 + line6 + line7)
    filindxDHsol.write(line1 + line2 + line3 + line4 + line5 + line6 + line7)
    filindxDHslv.write(line1 + line2 + line3 + line4 + line5 + line6 + line7)
    filindxDG.write(line1 + line2 + line3 + line4 + line5 + line6 + line7)
    filindxDSvib.write(line1 + line2 + line3 + line4 + line5 + line6 + line7)
    filindxDSlib.write(line1 + line2 + line3 + line4 + line5 + line6 + line7)
    filindxDSori.write(line1 + line2 + line3 + line4 + line5 + line6 + line7)
    filindxoccupancy.write(line1 + line2 + line3 + line4 + line5 + line6 + line7)

    strdh = ""
    strdhsol = ""
    strdhslv = ""
    strds = ""
    strdg = ""
    strocc = ""
    strvib = ""
    strlib = ""
    strori = ""

    entries = 0
    cond = 0

    for x in range(0, nx):
        for y in range(0, ny):
            for z in range(0, nz):
                count = grid['points'][x][y][z]['count']
                index += 1
                if count > grid_count_cutoff:
                    density = grid['points'][x][y][z]['density']
                    DG = grid['points'][x][y][z]['delta_G']
                    DH = grid['points'][x][y][z]['delta_H']
                    DHSOL = grid['points'][x][y][z]['delta_H_sol']
                    DHSLV = grid['points'][x][y][z]['delta_H_slv']
                    minTDS = grid['points'][x][y][z]['minT_delta_S']
                    minTDSvib = grid['points'][x][y][z]['minT_delta_S_vib']
                    minTDSlib = grid['points'][x][y][z]['minT_delta_S_lib']
                    minTDSori = grid['points'][x][y][z]['minT_delta_S_ori']
                else:
                    density = 0.0
                    DH = 0.0
                    DHSOL = 0.0
                    DHSLV = 0.0
                    minTDSvib = 0.0
                    minTDSlib = 0.0
                    minTDSori = 0.0
                    DH = 0.0

                minTDS = minTDSvib + minTDSlib + minTDSori
                DG = DH + minTDS

                strdh += " %8.4f" % DH
                strdhsol += " %8.4f" % DHSOL
                strdhslv += " %8.4f" % DHSLV
                strds += " %8.4f" % minTDS
                strdg += " %8.4f" % DG
                strocc += " %8.4f" % density
                strvib += " %8.4f" % minTDSvib
                strlib += " %8.4f" % minTDSlib
                strori += " %8.4f" % minTDSori

                if index == 1:
                    dg_moe += "[%8.4f," % DG
                    dh_moe += "[%8.4f," % DH
                    dhsolute_moe += "[%8.4f," % DHSOL
                    dhsolvent_moe += "[%8.4f," % DHSLV
                    mintds_moe += "[%8.4f," % minTDS
                    mintdsvib_moe += "[%8.4f," % minTDSvib
                    mintdslib_moe += "[%8.4f," % minTDSlib
                    mintdsori_moe += "[%8.4f," % minTDSori
                    occupancy_moe += "[%8.4f," % density
                if index > 1:
                    dg_moe += "%8.4f," % DG
                    dh_moe += "%8.4f," % DH
                    dhsolute_moe += "%8.4f," % DHSOL
                    dhsolvent_moe += "%8.4f," % DHSLV
                    mintds_moe += "%8.4f," % minTDS
                    mintdsvib_moe += "%8.4f," % minTDSvib
                    mintdslib_moe += "%8.4f," % minTDSlib
                    mintdsori_moe += "%8.4f," % minTDSori
                    occupancy_moe += "%8.4f," % density

                entries += 1
                if (entries % 3) == 0:
                    strdh += "\n"
                    strdhsol += "\n"
                    strdhslv += "\n"
                    strds += "\n"
                    strdg += "\n"
                    strocc += "\n"
                    strvib += "\n"
                    strlib += "\n"
                    strori += "\n"

    filindxDG.write(strdg + "\n")
    filindxDH.write(strdh + "\n")
    filindxDHsol.write(strdhsol + "\n")
    filindxDHslv.write(strdhslv + "\n")
    filindxDS.write(strds + "\n")
    filindxDSvib.write(strvib + "\n")
    filindxDSlib.write(strlib + "\n")
    filindxDSori.write(strori + "\n")

    filindxoccupancy.write(strocc + "\n")

    filindxDG.close()
    filindxDH.close()
    filindxDHsol.close()
    filindxDHslv.close()
    filindxDS.close()
    filindxDSvib.close()
    filindxDSlib.close()
    filindxDSori.close()
    filindxoccupancy.close()

    #
    # Region file for entire grid
    #
    region = open(os.path.join(grid_dir, "all.region"), 'w')
    regstr = ""
    for x in range(1, entries + 1):
        regstr += " %d " % x
        if (x % 10) == 0:
            regstr += "\n"
    regstr += "\n"
    region.write(regstr)
    region.close()


    ##### ALSO output MOE grid files.  Which has two files one shape file for both weighted and unweighted MOE grid files
    ##### this contains the [x1,x2,....xn],[y1,y2,.....yn], [z1,z2....zn] of all the grid points
    ##### All data points will be labeled each with simply a file with [d_x1_y1_z1, d_x1_y1_z2....d_xn_yn_zn]
    moe_shape = open(os.path.join(grid_dir, "shape"), 'w')
    regstr = "[%s,%s,%s]" % (nx, ny, nz)
    moe_shape.write(regstr)
    moe_shape.close()

    moe_outputs = {"dg_moe": dg_moe, "dh_moe": dh_moe, "dhsolute_moe": dhsolute_moe, "dhsolvent_moe": dhsolvent_moe,
                   "mintds_moe": mintds_moe, "mintdsvib_moe": mintdsvib_moe, "mintdslib_moe": mintdslib_moe,
                   "mintdsori_moe": mintdsori_moe, "occupancy_moe": occupancy_moe}

    for i in moe_outputs.keys():
        output = open(os.path.join(grid_dir, i), 'w')
        moestr = moe_outputs[i] + "]\n"
        output.write(moestr)
        output.close()


def _outputGridParams(grid, grid_dir, grid_count_cutoff, BULKORI, BULKFORCES, BULKTORQUES, BULKNRG, BULKDENSITY,
                      BULKWATMOLPER_VOXELVOL):
    nx, ny, nz = grid['dimensions']
    originx, originy, originz = grid['origin']
    step = grid['step']
    point_volume = step * step * step
    nsnapshots = grid['nsnapshots']

    kcalmolA_to_Nm = ( KCAL_TO_J / NAVOGADRO ) * 1e20

    index = 0

    outfile = os.path.join(grid_dir, "grid.forces")
    gridforces = open(outfile, 'w')

    # grid details
    gridforces.write("#Grid Dimensions Origin Step \n")
    line = "%8d\t%8d\t%8d\t%10.5f\t%10.5f\t%10.5f\t%10.5f\n" % (nx, ny, nz, originx, originy, originz, step)
    gridforces.write(line)
    # Could add standard deviation of each parameter
    gridforces.write("#Point\tX\tY\tZ\tFx\tFy\tFz\tTx\tTy\tTz\tOri\tDensity\tVolume\tEnergy\tSolEnergy\n")

    for x in range(0, nx):
        for y in range(0, ny):
            for z in range(0, nz):
                index += 1
                count = grid['points'][x][y][z]['count']
                if count > grid_count_cutoff:
                    density = _calcDensity(point_volume, count / nsnapshots)
                else:
                    density = 0
                meanF = grid['points'][x][y][z]['avgforces']
                meanT = grid['points'][x][y][z]['avgtorques']
                volume = point_volume
                # !!!! The forces from the cell files are in kcal/mol/A
                # However, for consistency with the literature, the BULK forces are
                # in N/m x  1e-10 so we need to convert
                for m in range(0, 3):
                    meanF[m] = meanF[m] * kcalmolA_to_Nm
                    meanT[m] = meanT[m] * kcalmolA_to_Nm
                meanori = grid['points'][x][y][z]['avgori']
                meannrg = grid['points'][x][y][z]['avgnrg'][0]
                meansolnrg = grid['points'][x][y][z]['avgsolnrg'][0]
                #cartx = originx + step * x
                #carty = originy + step * y
                #cartz = originz + step * z

                line = "%9d %8d %8d %8d %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %12.3f %10.5f %10.5f %10.5f\n" % (
                    index, x, y, z, meanF[0], meanF[1], meanF[2], meanT[0], meanT[1], meanT[2], meanori, density,
                    volume,
                    meannrg, meansolnrg)
                gridforces.write(line)
    gridforces.close()


def _outputWeightedGrid(grid, grid_dir, grid_count_cutoff, BULKORI, BULKFORCES, BULKTORQUES, BULKNRG, BULKDENSITY,
                        BULKWATMOLPER_VOXELVOL):
    #
    # Format a PDB style file with the cartesian coordinate of the grid points. Use columns like
    # occupancy or B-factors to write the entropies enthalpies etc...
    #

    titlepdbDG = '%s/grid_DG_weighted.pdb' % (grid_dir)
    filinpdbDG = open(titlepdbDG, 'w')

    titlepdbDH = '%s/grid_DH_weighted.pdb' % (grid_dir)
    filinpdbDH = open(titlepdbDH, 'w')

    titlepdbDHSOL = '%s/grid_DHsolute_weighted.pdb' % (grid_dir)
    filinpdbDHSOL = open(titlepdbDHSOL, 'w')

    titlepdbDHSLV = '%s/grid_DHsolvent_weighted.pdb' % (grid_dir)
    filinpdbDHSLV = open(titlepdbDHSLV, 'w')

    titlepdbDS = '%s/grid_minTDS_weighted.pdb' % (grid_dir)
    filinpdbDS = open(titlepdbDS, 'w')

    titlepdbDSvib = '%s/grid_minTDSvib_weighted.pdb' % (grid_dir)
    filinpdbDSvib = open(titlepdbDSvib, 'w')

    titlepdbDSlib = '%s/grid_minTDSlib_weighted.pdb' % (grid_dir)
    filinpdbDSlib = open(titlepdbDSlib, 'w')

    titlepdbDSori = '%s/grid_minTDSori_weighted.pdb' % (grid_dir)
    filinpdbDSori = open(titlepdbDSori, 'w')

    #titlepdboccupancy='%s/grid_density.pdb' % (grid_dir)
    #filinpdboccupancy=open( titlepdboccupancy, 'w')

    nx, ny, nz = grid['dimensions']
    originx, originy, originz = grid['origin']
    step = grid['step']
    point_volume = step * step * step

    nsnapshots = grid['nsnapshots']

    index = 0

    for x in range(0, nx):
        for y in range(0, ny):
            for z in range(0, nz):
                count = grid['points'][x][y][z]['count']
                if count > grid_count_cutoff:
                    abs_density = _calcDensity(point_volume, count / nsnapshots)
                    density = abs_density / BULKDENSITY
                    molpervoxel = density * BULKWATMOLPER_VOXELVOL
                    minTDSvib, minTDSlib, minTDSori, DHSOL, DHSLV = _calcPointProps(grid['points'][x][y][z], T, BULKORI,
                                                                                    BULKFORCES, BULKTORQUES, BULKNRG,
                                                                                    BULKDENSITY, BULKWATMOLPER_VOXELVOL)
                else:
                    density = 0.0
                    molpervoxel = 0.0
                    DHSOL = 0.0
                    DHSLV = 0.0
                    minTDSvib = 0.0
                    minTDSlib = 0.0
                    minTDSori = 0.0

                minTDS = minTDSvib + minTDSlib + minTDSori
                DH = DHSOL + DHSLV
                DG = DH + minTDS

                index += 1

                grid['points'][x][y][z]['delta_Gw'] = DG * molpervoxel
                grid['points'][x][y][z]['delta_Hw'] = DH * molpervoxel
                grid['points'][x][y][z]['delta_H_solw'] = DHSOL * molpervoxel
                grid['points'][x][y][z]['delta_H_slvw'] = DHSLV * molpervoxel
                grid['points'][x][y][z]['minT_delta_Sw'] = minTDS * molpervoxel
                grid['points'][x][y][z]['minT_delta_S_vibw'] = minTDSvib * molpervoxel
                grid['points'][x][y][z]['minT_delta_S_libw'] = minTDSlib * molpervoxel
                grid['points'][x][y][z]['minT_delta_S_oriw'] = minTDSori * molpervoxel
                grid['points'][x][y][z]['densityw'] = density

                cartx = originx + step * x
                carty = originy + step * y
                cartz = originz + step * z

                textDG = 'ATOM%7s%2s%24.2f%8.2f%8.2f%5s %10.4f %10.4f%12s' % (
                    index, 'H', cartx, carty, cartz, '1.00', DG * molpervoxel, density, 'H')
                textDH = 'ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (
                    index, 'H', cartx, carty, cartz, '1.00', DH * molpervoxel, density, 'H')
                textDHSOL = 'ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (
                    index, 'H', cartx, carty, cartz, '1.00', DHSOL * molpervoxel, density, 'H')
                textDHSLV = 'ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (
                    index, 'H', cartx, carty, cartz, '1.00', DHSLV * molpervoxel, density, 'H')
                textDS = 'ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (
                    index, 'H', cartx, carty, cartz, '1.00', minTDS * molpervoxel, density, 'H')
                textDSvib = 'ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (
                    index, 'H', cartx, carty, cartz, '1.00', minTDSvib * molpervoxel, density, 'H')
                textDSlib = 'ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (
                    index, 'H', cartx, carty, cartz, '1.00', minTDSlib * molpervoxel, density, 'H')
                textDSori = 'ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s' % (
                    index, 'H', cartx, carty, cartz, '1.00', minTDSori * molpervoxel, density, 'H')

                filinpdbDG.write(textDG + "\n")
                filinpdbDH.write(textDH + "\n")
                filinpdbDHSOL.write(textDHSOL + "\n")
                filinpdbDHSLV.write(textDHSLV + "\n")
                filinpdbDS.write(textDS + "\n")
                filinpdbDSvib.write(textDSvib + "\n")
                filinpdbDSlib.write(textDSlib + "\n")
                filinpdbDSori.write(textDSori + "\n")

                #sys.exit(-1)

    filinpdbDG.close()
    filinpdbDH.close()
    filinpdbDHSOL.close()
    filinpdbDHSLV.close()
    filinpdbDS.close()
    filinpdbDSvib.close()
    filinpdbDSlib.close()
    filinpdbDSori.close()

    # Format a DX style file with the cartesian coordinate of the grid points. Use columns like
    # occupancy or B-factors to write the entropies enthalpies etc...
    #
    titledxDG = '%s/grid_DG_weighted.dx' % (grid_dir)
    filindxDG = open(titledxDG, 'w')

    titledxDH = '%s/grid_DH_weighted.dx' % (grid_dir)
    filindxDH = open(titledxDH, 'w')

    titledxDHSOL = '%s/grid_DHsolute_weighted.dx' % (grid_dir)
    filindxDHSOL = open(titledxDHSOL, 'w')

    titledxDHSLV = '%s/grid_DHsolvent_weighted.dx' % (grid_dir)
    filindxDHSLV = open(titledxDHSLV, 'w')

    titledxDS = '%s/grid_minTDS_weighted.dx' % (grid_dir)
    filindxDS = open(titledxDS, 'w')

    titledxDSvib = '%s/grid_minTDSvib_weighted.dx' % (grid_dir)
    filindxDSvib = open(titledxDSvib, 'w')

    titledxDSlib = '%s/grid_minTDSlib_weighted.dx' % (grid_dir)
    filindxDSlib = open(titledxDSlib, 'w')

    titledxDSori = '%s/grid_minTDSori_weighted.dx' % (grid_dir)
    filindxDSori = open(titledxDSori, 'w')

    nx, ny, nz = grid['dimensions']
    originx, originy, originz = grid['origin']
    step = grid['step']

    index = 0

    ######FORMAT an appropriate DX file for both DH and DS ###############
    #####header of DX file ##############
    #object 1 class gridpositions counts 16 17 15
    line1 = 'object 1 class gridpositions counts %i %i %i\n' % (
        grid['dimensions'][0], grid['dimensions'][1], grid['dimensions'][2])
    #origin -1.0 -1.0  0.0
    line2 = 'origin %4.2f %4.2f %4.2f\n' % (grid['origin'][0], grid['origin'][1], grid['origin'][2])
    #delta 1 0 0
    line3 = 'delta %4.2f 0 0\n' % (grid['step'])
    #delta 0 1 0
    line4 = 'delta 0 %4.2f 0\n' % (grid['step'])
    #delta 0 0 1
    line5 = 'delta 0 0 %4.2f\n' % (grid['step'])
    #object 2 class gridconnections counts 16 17 15
    line6 = 'object 2 class gridconnections counts %i %i %i\n' % (
        grid['dimensions'][0], grid['dimensions'][1], grid['dimensions'][2])
    #object 3 class array type double rank 0 items 4080 data follows
    line7 = 'object 3 class array type double rank 0 items %i data follows\n' % (
        grid['dimensions'][0] * grid['dimensions'][1] * grid['dimensions'][2])

    filindxDS.write(line1 + line2 + line3 + line4 + line5 + line6 + line7)
    filindxDH.write(line1 + line2 + line3 + line4 + line5 + line6 + line7)
    filindxDHSOL.write(line1 + line2 + line3 + line4 + line5 + line6 + line7)
    filindxDHSLV.write(line1 + line2 + line3 + line4 + line5 + line6 + line7)
    filindxDG.write(line1 + line2 + line3 + line4 + line5 + line6 + line7)
    filindxDSvib.write(line1 + line2 + line3 + line4 + line5 + line6 + line7)
    filindxDSlib.write(line1 + line2 + line3 + line4 + line5 + line6 + line7)
    filindxDSori.write(line1 + line2 + line3 + line4 + line5 + line6 + line7)

    ### keep list of each data file for moe
    dg_moe = ""
    dh_moe = ""
    dhsolute_moe = ""
    dhsolvent_moe = ""
    mintds_moe = ""
    mintdsvib_moe = ""
    mintdslib_moe = ""
    mintdsori_moe = ""
    occupancy_moe = ""

    strdh = ""
    strdhsol = ""
    strdhslv = ""
    strds = ""
    strdg = ""
    strvib = ""
    strlib = ""
    strori = ""

    entries = 0
    cond = 0

    for x in range(0, nx):
        for y in range(0, ny):
            for z in range(0, nz):
                index += 1
                DG = grid['points'][x][y][z]['delta_Gw']
                DH = grid['points'][x][y][z]['delta_Hw']
                DHSOL = grid['points'][x][y][z]['delta_H_solw']
                DHSLV = grid['points'][x][y][z]['delta_H_slvw']

                minTDS = grid['points'][x][y][z]['minT_delta_Sw']
                minTDSvib = grid['points'][x][y][z]['minT_delta_S_vibw']
                minTDSlib = grid['points'][x][y][z]['minT_delta_S_libw']
                minTDSori = grid['points'][x][y][z]['minT_delta_S_oriw']
                density = grid['points'][x][y][z]['densityw']

                minTDS = minTDSvib + minTDSlib + minTDSori

                strdh += " %8.4f" % (DH)
                strdhsol += " %8.4f" % (DHSOL)
                strdhslv += " %8.4f" % (DHSLV)
                strds += " %8.4f" % (minTDS)
                strdg += " %8.4f" % (DG)
                strvib += " %8.4f" % (minTDSvib)
                strlib += " %8.4f" % (minTDSlib)
                strori += " %8.4f" % (minTDSori)

                if index == 1:
                    dg_moe += "[%8.4f, " % DG
                    dh_moe += "[%8.4f, " % DH
                    dhsolute_moe += "[%8.4f, " % DHSOL
                    dhsolvent_moe += "[%8.4f, " % DHSLV
                    mintds_moe += "[%8.4f, " % minTDS
                    mintdsvib_moe += "[%8.4f, " % minTDSvib
                    mintdslib_moe += "[%8.4f, " % minTDSlib
                    mintdsori_moe += "[%8.4f, " % minTDSori
                    occupancy_moe += "[%8.4f, " % density
                if index > 1:
                    dg_moe += "%8.4f, " % DG
                    dh_moe += "%8.4f, " % DH
                    dhsolute_moe += "%8.4f, " % DHSOL
                    dhsolvent_moe += "%8.4f, " % DHSLV
                    mintds_moe += "%8.4f, " % minTDS
                    mintdsvib_moe += "%8.4f, " % minTDSvib
                    mintdslib_moe += "%8.4f, " % minTDSlib
                    mintdsori_moe += "%8.4f, " % minTDSori
                    occupancy_moe += "%8.4f, " % density

                entries += 1
                if (entries % 3) == 0:
                    strdh += "\n"
                    strdhsol += "\n"
                    strdhslv += "\n"
                    strds += "\n"
                    strdg += "\n"
                    strvib += "\n"
                    strlib += "\n"
                    strori += "\n"

    filindxDG.write(strdg + "\n")
    filindxDH.write(strdh + "\n")
    filindxDHSOL.write(strdhsol + "\n")
    filindxDHSLV.write(strdhslv + "\n")
    filindxDS.write(strds + "\n")
    filindxDSvib.write(strvib + "\n")
    filindxDSlib.write(strlib + "\n")
    filindxDSori.write(strori + "\n")

    filindxDG.close()
    filindxDH.close()
    filindxDHSOL.close()
    filindxDHSLV.close()
    filindxDS.close()
    filindxDSvib.close()
    filindxDSlib.close()
    filindxDSori.close()

    moe_outputs = {"dg_moe_weighted": dg_moe, "dh_moe_weighted": dh_moe, "dhsolute_moe_weighted": dhsolute_moe,
                   "dhsolvent_moe_weighted": dhsolvent_moe, "mintds_moe_weighted": mintds_moe,
                   "mintdsvib_moe_weighted": mintdsvib_moe, "mintdslib_moe_weighted": mintdslib_moe,
                   "mintdsori_moe_weighted": mintdsori_moe}

    for i in moe_outputs.keys():
        output = open(os.path.join(grid_dir, i), 'w')
        moestr = moe_outputs[i] + "]\n"
        output.write(moestr)
        output.close()


###### region_properties functions

def _loadGrid(gridforces):
    grid = {}

    stream = open(gridforces, 'r')
    buffer = stream.readlines()
    stream.close()

    # the grid dimensions/origin/step
    nx, ny, nz, ox, oy, oz, step = buffer[1].split()

    grid['params'] = [int(nx), int(ny), int(nz), float(ox), float(oy), float(oz), float(step)]
    grid['step'] = float(step)
    grid['dimensions'] = [int(nx), int(ny), int(nz)]
    grid['points'] = {}

    for line in buffer[2:]:
        if line.startswith("#"):
            continue
        elems = line.split()
        index = int(elems[0])
        x = int(elems[1])
        y = int(elems[2])
        z = int(elems[3])
        fx = float(elems[4])
        fy = float(elems[5])
        fz = float(elems[6])
        tx = float(elems[7])
        ty = float(elems[8])
        tz = float(elems[9])
        ori = float(elems[10])
        density = float(elems[11])
        volume = float(elems[12])
        nrg = float(elems[13])
        solnrg = float(elems[14])

        grid['points'][index] = {}
        grid['points'][index]['coords'] = (x, y, z)
        grid['points'][index]['forces'] = (fx, fy, fz)
        grid['points'][index]['torques'] = (tx, ty, tz)
        grid['points'][index]['ori'] = ori
        grid['points'][index]['density'] = density
        grid['points'][index]['volume'] = volume
        grid['points'][index]['nrg'] = nrg
        grid['points'][index]['solnrg'] = solnrg

    return grid


def _loadRegion(regionfile):
    region = []

    stream = open(regionfile, 'r')
    buffer = stream.readlines()
    stream.close()

    for line in buffer:
        if line.startswith("#"):
            continue
        elems = line.split()
        for value in elems:
            region.append(int(value))

    return region


def _getProperties(region, grid, BULKORI, BULKFORCES, BULKTORQUES, BULKNRG, BULKDENSITY):
    #
    # Evaluate thermodynamic properties over the grid region
    #
    #

    region_DHW = 0.0
    region_DHX = 0.0
    density = 0.0
    region_f = [0.0, 0.0, 0.0]
    region_t = [0.0, 0.0, 0.0]
    region_ori = 0.0
    region_volume = 0.0
    region_nw = 0.0

    for p in region:
        if p not in grid:
            print("Problem ! Point %s in your region file is not present in the grid file !. Abort. ")
            sys.exit(-1)
        density = grid[p]['density']
        volume = grid[p]['volume']
        region_volume += volume
        if density < SMALL:
            # No data on this point, so no contribution
            continue
        density_prefactor = ( volume * NAVOGADRO * 1e-27) / ( MASS_WATER )
        # Each point contains a variable number of water molecules, so to recover properties over whole volume from grid, we
        # must weight the contribution of each point by local density (water/A**3)
        density_conv = density * density_prefactor
        #print density, density_conv

        region_DHW += (grid[p]['nrg'] - BULKNRG) * density_conv
        region_DHX += (grid[p]['solnrg']) * density_conv
        region_ori += grid[p]['ori'] * density_conv
        for m in range(0, 3):
            region_f[m] += grid[p]['forces'][m] * density_conv
            region_t[m] += grid[p]['torques'][m] * density_conv

        region_nw += density_conv
    #

    if region_nw > 0:
        region_ori /= region_nw
        for m in range(0, 3):
            region_f[m] /= region_nw
            region_t[m] /= region_nw

    oriratio = 1.0
    vibratio = 1.0
    libratio = 1.0

    if region_nw > 0:
        if region_ori < 1.0:
            region_ori = 1.0
        oriratio *= ( region_ori / BULKORI )
        for m in range(0, 3):
            if region_f[m] == 0:
                continue
            if region_t[m] == 0:
                continue
            vibratio *= ( BULKFORCES[m] / (region_f[m]) )
            libratio *= ( BULKTORQUES[m] / (region_t[m]) )
    #print "vibratio", vibratio
    region_minTDSvib = -KB * T * region_nw * math.log(vibratio)
    #print region_minTDSvib
    #print libratio
    region_minTDSlib = -KB * T * region_nw * math.log(libratio)
    #print oriratio
    region_minTDSori = -KB * T * region_nw * math.log(oriratio)

    #sys.exit(-1)

    region_DH = region_DHW + region_DHX
    region_minTDS = region_minTDSvib + region_minTDSlib + region_minTDSori
    region_DG = region_DH + region_minTDS
    #print region_nw
    #print region_volume
    if region_nw > 0:
        #print region_nw
        #print region_volume
        #print density_prefactor
        region_density = (region_nw / region_volume) / ( ( 1 * NAVOGADRO * 1e-27) / ( MASS_WATER ) )
        #print region_density
    else:
        region_density = 0.0
    #print region_density
    region_density /= BULKDENSITY
    #print region_density

    return (
    region_DG, region_DH, region_minTDS, region_DHW, region_DHX, region_minTDSvib, region_minTDSlib, region_minTDSori,
    region_density, region_nw, region_volume)


def _writeRegionProperties(regionprops, root):
    region_stats_file = root + "_properties.dat"

    stream = open(region_stats_file, 'w')

    #print (regionprops)

    line = "#Excess Thermodynamic properties of water within region\n"
    stream.write(line)
    line = """
Delta_G = %8.4f kcal/mol
Delta_H = %8.4f kcal/mol
-T*Delta_S = %8.4f kcal/mol
Delta_Hw = %8.4f kcal/mol
Delta_Hx = %8.4f kcal/mol
-T*Delta_Svib = %8.4f kcal/mol
-T*Delta_Slib = %8.4f kcal/mol
-T*Delta_Sori = %8.4f kcal/mol
Relative_Density = %8.4f
Average_number_of_water_molecules = %8.4f
Volume = %8.4f A^3
""" % (regionprops[0], regionprops[1], regionprops[2], regionprops[3], regionprops[4], regionprops[5], regionprops[6],
       regionprops[7], regionprops[8], regionprops[9], regionprops[10])

    stream.write(line)

    stream.close()


def _writeRegionPDB(region, grid, root, BULKDENSITY):
    region_pdb_file = root + ".pdb"

    stream = open(region_pdb_file, 'w')

    originx = grid['params'][3]
    originy = grid['params'][4]
    originz = grid['params'][5]
    step = grid['params'][6]

    for p in region:
        if p not in grid['points']:
            print("Problem ! Point %s in your region file is not present in the grid file !. Abort. ")
            sys.exit(-1)
        coords = grid['points'][p]['coords']
        cartx = originx + coords[0] * step
        carty = originy + coords[1] * step
        cartz = originz + coords[2] * step
        density = grid['points'][p]['density']
        #volume = grid[p]['volume']
        #density_prefactor = ( volume * NAVOGADRO * 1e-27) / ( MASS_WATER )
        #density_conv = density * density_prefactor
        #print (density)

        line = 'ATOM%7s%3s   SIT %5d    %8.3f%8.3f%8.3f%6.2f%6.2f%12s\n' % (
        p, 'C', 1, cartx, carty, cartz, 1.00, density / BULKDENSITY, 'C')
        #line='ATOM%7s%2s%24.2f%8.2f%8.2f%5s %8.4f %10.4f%12s\n' % (p, 'C', cartx, carty, cartz, '1.00', 0.00, 0.000, 'H')
        stream.write(line)
    stream.close()


#### clustergrid functions
def _define_neighbor_list(step, neighcut):
    #
    # Precompute list of neighboring points that will be within cutoff distance of a point.
    #
    nlist = []

    gridcut2 = neighcut ** 2
    step2 = step ** 2
    cut2 = gridcut2 / step2

    max = int(math.sqrt(cut2)) + 1

    for x in range(-max, max + 1 ):
        for y in range(-max, max + 1):
            for z in range(-max, max + 1):
                if x == 0 and y == 0 and z == 0:
                    continue
                #print x,y,z
                d2 = x ** 2 + y ** 2 + z ** 2
                #print d2
                if d2 < cut2:
                    nlist.append(( x, y, z))

    #print nlist
    #print len(nlist)
    #sys.exit(-1)

    return nlist


def _cluster(points, grid, lowt, neighcut, BULKDENSITY):
    sites_count = 0
    sites = {}

    #print (grid['step'])
    nlist = _define_neighbor_list(grid['step'], neighcut)

    maxx = grid['params'][0] - 1
    maxy = grid['params'][1] - 1
    maxz = grid['params'][2] - 1

    assigngrid = {}
    aidx = 0
    for x in range(0,maxx + 1):
        assigngrid[x] = {}
        for y in range(0,maxy + 1 ):
            assigngrid[x][y] = {}
            for z in range(0,maxz + 1 ):
                aidx += 1
                assigngrid[x][y][z] = {}
                assigngrid[x][y][z]['index'] = aidx
                assigngrid[x][y][z]['site'] = 0

    cutval = lowt*BULKDENSITY

    for point in points:
        index = point[0]
        val = point[1]
        px = grid['points'][index]['coords'][0]
        py = grid['points'][index]['coords'][1]
        pz = grid['points'][index]['coords'][2]
        #print val, index, px, py, pz
        #print "Doing point %s %s " % (index,val)
        # Skip if already assigned
        if assigngrid[px][py][pz]['site'] > 0:
            #print "Already assigned"
            continue
        if val < cutval:
            #print "Density dropped below threshold, break."
            break
        # Assign to a new site
        sites_count += 1
        sites[sites_count] = {}
        sites[sites_count]['points'] = []
        sites[sites_count]['center'] = index
        sites[sites_count]['density'] = val
        sites[sites_count]['points'].append(index)
        sites[sites_count]['npoints'] = 1
        assigngrid[px][py][pz]['index'] = sites_count
        # Assign all unassigned neighboring points to that site
        #print "Finding neighbors of %s %s %s " % (px,py,pz)
        for n in nlist:
            nx = px + n[0]
            ny = py + n[1]
            nz = pz + n[2]
            if ( ( nx > maxx or nx < 0 ) or
                 ( ny > maxy or ny < 0 ) or
                 ( nz > maxz or nz < 0 )
                 ):
                continue
            #print nx, ny, nz
            if assigngrid[nx][ny][nz]['site'] > 0:
                #print "Already assigned"
                continue
            idx = assigngrid[nx][ny][nz]['index']
            #print idx
            sites[sites_count]['points'].append(idx)
            sites[sites_count]['density'] += grid['points'][idx]['density']
            sites[sites_count]['npoints'] += 1
            assigngrid[nx][ny][nz]['site'] = sites_count
        #sys.exit(-1)

    #print sites
    #print sites_count
    # sys.exit(-1)
    return sites


def _outputSites(sites, sites_dir):
    if os.path.exists(sites_dir):
        cmd = "rm -rf %s " % sites_dir
        os.system(cmd)

    cmd = "mkdir -p %s " % sites_dir
    os.system(cmd)

    siteindices = list(sites.keys())
    siteindices.sort()

    for siteidx in siteindices:
        center = sites[siteidx]['center']
        points = sites[siteidx]['points']
        filename = "site%000006d.region" % siteidx
        #print filename
        path = os.path.join(sites_dir, filename)
        stream = open(path, 'w')
        stream.write("#Center %d\n" % center)
        thestr = ""
        x = 0
        for p in points:
            x += 1
            if (x % 10 ) == 0:
                thestr += "\n"
            thestr += " %d " % p
        thestr += "\n"
        stream.write(thestr)
        stream.close()
        #sys.exit(-1)

    ### get a representation of all clustered centers!
    filecenter = 'site_allcenters.region'
    writeinto = os.path.join(sites_dir, filecenter)
    centerdata = open(writeinto, 'w')
    thestr = ""
    x = 0
    for siteidx in siteindices:
        center = sites[siteidx]['center']
        x += 1
        if (x % 10 ) == 0:
            thestr += "\n"
        thestr += " %d " % center
    thestr += "\n"
    centerdata.write(thestr)
    centerdata.close()


def _getSiteProperties(sites, grids):
    sitesindices = sites.keys()
    sitesindices.sort()

    for siteidx in sitesindices:
        site = sites[siteidx]
        #print site
        sum_density = 0.0
        site_DG = 0.0
        site_DH = 0.0
        site_minTDS = 0.0
        site_minTDSvib = 0.0
        site_minTDSlib = 0.0
        site_minTDSori = 0.0
        for point in site['points']:
            density = grids["density"]['index'][point][0]
            #print density
            DG = grids["DG"]['index'][point][0]
            DH = grids["DH"]['index'][point][0]
            minTDS = grids["minTDS"]['index'][point][0]
            minTDSvib = grids["minTDSvib"]['index'][point][0]
            minTDSlib = grids["minTDSlib"]['index'][point][0]
            minTDSori = grids["minTDSori"]['index'][point][0]
            #print density, minTDSori

            sum_density += density
            site_DG += DG * density
            site_DH += DH * density
            site_minTDS += minTDS * density
            site_minTDSvib += minTDSvib * density
            site_minTDSlib += minTDSlib * density
            site_minTDSori += minTDSori * density
        site_DG /= sum_density
        site_DH /= sum_density
        site_minTDS /= sum_density
        site_minTDSvib /= sum_density
        site_minTDSlib /= sum_density
        site_minTDSori /= sum_density

        #print sum_density


        site['properties'] = {"DG": site_DG, "DH": site_DH, "minTDS": site_minTDS,
                              "minTDSvib": site_minTDSvib, "minTDSlib": site_minTDSlib,
                              "minTDSori": site_minTDSori}
        #print site['properties']
        #sys.exit(-1)

#### choose parameters when running clustering script....does not interface with command line
def _regionproperties_clu(regionfile, gridforces, BULKORI, BULKFORCES, BULKTORQUES, BULKNRG, BULKDENSITY):
    """
    regionproperties works out the excess thermodynamic properties of water molecules that were contained in
    a user defined region of space."""

    name = os.path.split(regionfile)[-1]
    root = name.replace(".region", "")
    #print name, root
    #sys.exit(-1)

    #tb = time.time()

    # Load the grid
    grid = _loadGrid(gridforces)
    # Load the region
    region = _loadRegion(regionfile)
    # Get region properties
    region_props = _getProperties(region, grid['points'], BULKORI, BULKFORCES, BULKTORQUES, BULKNRG, BULKDENSITY)
    # Write region properties
    _writeRegionProperties(region_props, root)
    # Write a pdb to visualize the region
    _writeRegionPDB(region, grid, root, BULKDENSITY)

    #ta = time.time()
    #if benchmark:
    #    print("time to generate region properties %5.3f seconds " % (ta - tb))

##### sub-dx and average-dx functions

def _loadheader(grid_file):
    stream = open(grid_file,"r")
    buffer = stream.readlines()
    stream.close()
    header = buffer[0:7]
    return header

def _loadgrid_dx(grid_file):
    stream = open(grid_file,"r")
    buffer = stream.readlines()
    stream.close()
    data = []
    for line in buffer[7:]:
        elems = line.split()
        for e in elems:
            data.append(float(e))

    return data

def _subtractgrids(grid1, grid2):
    grid_sub = copy.deepcopy(grid1)
    for x in range(0,len(grid2)):
        grid_sub[x] -= grid2[x]
    return grid_sub

def _averagegrids(grids):
    avg_grid = copy.deepcopy(grids[0])
    n = len(grids)
    #print n
    for grid in grids[1:]:
        for x in range(0,len(grid)):
            avg_grid[x] += grid[x]
    for x in range(0,len(avg_grid)):
        avg_grid[x] /= n
    return avg_grid

####### API #######################
@resolveParameters
def traj2cell():
    """
    traj2cell converts a trajectory file into a collection of cell files that stores cell theory
    parameters for the water molecules. If a grid region has been defined, ony water molecules
    within this region are analysed. """

    #print (grid_center_x.val, grid_center_y.val, grid_center_z.val, grid_plus_x.val, grid_min_x.val, grid_plus_y.val, grid_min_y.val, grid_plus_z.val, grid_min_z.val)
    #sys.exit(-1)

    top_file = topfile.val
    crd_file = crdfile.val
    dcd_file = trajfile.val

    start_frame = startframe.val
    end_frame = endframe.val
    BENCHMARK = benchmark.val

    # Now set the bulk coordination number of water to use
    BULKCN = bulk[watermodel.val]

    if not os.path.exists(cell_dir.val):
        cmd = "mkdir %s" % cell_dir.val
        os.system(cmd)

    # First create a Sire system using input top/crd files
    amber = Amber()
    molecules, space = amber.readCrdTop(crd_file, top_file)
    system = _createSystem(molecules)

    # Assign h-bond donors/acceptors
    atomsDAtype, mol_to_atoms = _assignDonorsAcceptors( system )

    # Init the grid boundaries

    grid = _initGrid( grid_center_x.val, grid_center_y.val, grid_center_z.val, grid_plus_x.val, grid_min_x.val, grid_plus_y.val, grid_min_y.val, grid_plus_z.val, grid_min_z.val )
    print (BENCHMARK)
    print (grid)
    #sys.exit(-1)

    e_gridwater = Symbol("E_{gridwater}")
    e_gridwater_solutes = Symbol("E_{gridwater:solutes}")

    # Now load the trajectory, selecting only the frames that have been assigned
    # to this instance
    mdtraj_top = mdtraj.load_prmtop(top_file)
    mdtraj_dcdfile = mdtraj.open(dcd_file,'r')
    nframes = len(mdtraj_dcdfile)

    if end_frame > (nframes - 1):
        end_frame = nframes - 1

    mdtraj_dcdfile.seek(start_frame)

    current_frame = start_frame

    while (current_frame <= end_frame):
        frames_xyz, cell_lengths, cell_angles = mdtraj_dcdfile.read(n_frames=1)
        print ("Processing frame %s " % current_frame)
        tb = time.time()
        system, polarH_coords, acceptor_coords, alternate_acceptor_coords = _updateSystemfromDCDgrid(system, frames_xyz, cell_lengths, cell_angles, atomsDAtype, grid)
        # Dump a PDB of the system if this is the first frame, to facilitate visualisation of cell/grid/sites
        #if current_frame == 0:
        #    PDB().write(system[MGName("all")], "frame-0.pdb" )
        ta = time.time()
        if BENCHMARK:
            print ("time to update coords %5.3f seconds" % (ta - tb) )
        tb = ta
        # Initialise energies, this appears necessary to initialise properly
        # LJPair database for force calculations
        system.energy()

        # Use a custom RBWorkspace to compute forces and torques
        all = system[MGName("all")]
        wspace =  RBWorkspaceJM(all)
        wspace.setSystem(system)

        #print "About to call calculate forces"
        wspace.calculateForces(e_gridwater)

        bead_forces = wspace.beadForcesArray()
        bead_torques = wspace.beadTorquesArray()
        bead_energies = wspace.beadEnergiesArray()

        #print (bead_forces)
        #print ("####")
        #print (bead_torques)

        # Now the solute energy terms
        wspace.calculateForces(e_gridwater_solutes)
        beadsolute_energies = wspace.beadEnergiesArray()

        ta = time.time()
        if BENCHMARK:
            print ("time to get energies and forces/torques %5.3f seconds " % (ta - tb))
        tb = ta

        space = system.property("space")

        # detect waters that are first shell to solute acceptors
        #
        #
        # CHECK NOT BUGGY ! (NEED H-BOND SOLUTE)
        firstshellwaters = _defineFSW(space, acceptor_coords, atomsDAtype, grid)

        ta = time.time()
        if BENCHMARK:
            print ("time to detect first shell waters %5.3f seconds " % (ta - tb))
        tb = ta

        # h-bonds analysis
        h_bonds = _assignHbonds( space, polarH_coords, acceptor_coords, alternate_acceptor_coords, atomsDAtype, grid)

        ta = time.time()
        if BENCHMARK:
            print ("time to assign hbonds %5.3f seconds " % (ta - tb))
        tb = ta

        # a dictionnary of coordination numbers for each water's oxygens
        water_cn = _coordinationNumbers( space, acceptor_coords, h_bonds, atomsDAtype, firstshellwaters, grid)

        ta = time.time()
        if BENCHMARK:
            print ("time to do fastCN %5.3f seconds " % (ta - tb))
        tb = ta

        water_ori = _orientationalNumbers( water_cn, firstshellwaters, h_bonds, mol_to_atoms, atomsDAtype, BULKCN )

        #print (water_ori)

        ta = time.time()
        if BENCHMARK:
            print ("time to do orientational number %5.3f seconds " % (ta - tb))
        tb = ta

        # A complication is that the bead values are not in the same order
        # between different snapshots. So let's store the forces/torques/energies
        # values in dictionnaries keyed by the first residue number of
        # each mol (unique)
        molecules = system[MGName("all")].molecules()
        gridwater = system[MGName("gridwater")]
        molnums = molecules.molNums()
        molnums.sort()

        water_results = {}

        for x in range(0,len(molnums)):
            molnum = molnums[x]
            #import pdb; pdb.set_trace()
            molecule = system.molecule(molnum).molecule()
            # Ignore if not gridwater
            if not gridwater.contains(molecule):
                continue

            # Ignore non water molecules
            nresidues = molecule.nResidues()
            for y in range(nresidues):
                residue = molecule.residue(ResIdx(y))
                #print residue
                resname = residue.name().value()
                # This is a solute residue. Skip it.
                if not resname in Water_residue_names:
                    continue
                # This is a water molecule
                else:
                    # Skip if not within grid volume
                    #if not water_ori.has_key(molnum.value()):
                    #print (water_ori, molnum.value())
                    #sys.exit(-1)
                    if not molnum.value() in water_ori:
                        continue
                    waterflag = 1

                res_num = residue.number().value()
                mol_energies = bead_energies[x]
                mol_solenergies = beadsolute_energies[x]
                mol_forces = [ 0.0, 0.0, 0.0 ]
                mol_torques = [ 0.0, 0.0, 0.0 ]

                for i in range(0,3):
                    mol_forces[i] = 0.5 * abs( bead_forces[x][i] )
                    mol_torques[i] = 0.5 * abs( bead_torques[x][i] )

                water_results[res_num] = {}
                water_results[res_num]['waterflag'] = waterflag
                water_results[res_num]['forces'] = mol_forces
                water_results[res_num]['torques'] = mol_torques
                water_results[res_num]['energies'] = mol_energies
                water_results[res_num]['solenergies'] = mol_solenergies

                oatom = molecule.select( AtomName("O") )
                mol_o_coords = oatom.property("coordinates")
                water_results[res_num]['o_coords'] = mol_o_coords
                water_results[res_num]['ori'] = water_ori[ molnum.value() ]
                water_results[res_num]['cn'] = water_cn[ molnum.value() ]

        # JM May 13 - store frame number rather than time
        _writeCellDataBin(cell_dir.val, water_results, current_frame, grid['volume'] )

        ta = time.time()
        if BENCHMARK:
            print ("time to write output %5.3f seconds " % (ta - tb))
        tb = ta

        # End of loop over trajectory frame
        current_frame += 1
        #sys.exit(-1)
        #import pdb; pdb.set_trace()


@resolveParameters
def cell2grid():
    """
    cell2grid averages the contents of a collection of cell files and output a set of grid files in
    dx, MOE and pdb format that spatially resolve thermodynamic components of the excess free energy of water
    molecules."""

    BULKORI = waterProps[watermodel.val]['BULKORI']
    BULKFORCES = waterProps[watermodel.val]['BULKFORCES']
    BULKTORQUES = waterProps[watermodel.val]['BULKTORQUES']
    BULKNRG = waterProps[watermodel.val]['BULKNRG']
    BULKDENSITY = waterProps[watermodel.val]['BULKDENSITY']
    #################### change so that the bulk wat mol per voxel depends on the volume which
    #################### depends on the grid step size
    BULKWATMOLPER_VOXELVOL = (BULKDENSITY * (grid_step.val * 1e-9) ** 3 / MASS_WATER) * NAVOGADRO

    BENCHMARK = benchmark.val

    if os.path.exists(grid_dir.val):
        if Overwrite:
            cmd = "rm -rf %s " % grid_dir.val
            os.system(cmd)
        else:
            print("Output folder %s exists, but overwrite mode has not been set. Aborting.")
            sys.exit(-1)

    if not os.path.exists(grid_dir.val):
        cmd = "mkdir %s " % grid_dir.val

        #print cmd
        os.system(cmd)

    tb = time.time()

    # Initialise grid
    grid = _initGrid_c2g(grid_step.val)
    cell_interval = (startframe.val, endframe.val)
    # Get list of cell files
    cell_files = _getCellData(cell_dir.val, cell_interval)

    # Keep track of mean forces/torques/energies and orientational numbers for each water
    water_props = {}

    # Update grid with cell data
    ### make a new file with grid data
    stream = open("grid.dat", 'w')
    a = "#Frame       -TDS_vib -TDS_lib -TDS_ori  -TDSW      DHW    DHX    DG (kcal/mol) Forces  (N x m x 1e-10)  Torques (N x m x 1e-20)   MeanWatNrg    OriN   Density (g/L)\n"
    stream.write(a)

    print(
        "#Frame       -TDS_vib -TDS_lib -TDS_ori  -TDSW      DHW    DHX    DG (kcal/mol) Forces  (N x m x 1e-10)  Torques (N x m x 1e-20)   MeanWatNrg    OriN   Density (g/L)")

    framecount = 0

    for cell_file in cell_files:
        cell = _loadCell(cell_file)
        # sys.exit(-1)
        water_props = _updatePropsCumulative(water_props, cell)
        _updateGridPropsCumulative(grid, cell)
        #sys.exit(-1)
        framecount += 1

        if (framecount % frequencyupdate.val) == 0:
            framenum = framecount + cell_interval[0]
            DS_vib, DS_lib, DS_ori, DSW, DHSOL, DHSLV, DGW, meanF, meanT, meannrg, meanori, density = _calcAverageProps(
                water_props, BULKORI, BULKFORCES, BULKTORQUES, BULKNRG, BULKDENSITY, BULKWATMOLPER_VOXELVOL)
            outstr = "%12d %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f" % (
                framenum, DS_vib, DS_lib, DS_ori, DSW, DHSLV, DHSOL, DGW)
            outstr += "%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f " % (
                meanF[0], meanF[1], meanF[2], meanT[0], meanT[1], meanT[2] )
            outstr += "%8.5f %8.5f %8.5f\n" % (meannrg, meanori, density)
            print(outstr)
            stream.write(outstr)
            #sys.exit(-1)
    stream.close()

    _outputGrid(grid, grid_dir.val, grid_count_cutoff.val, BULKORI, BULKFORCES, BULKTORQUES, BULKNRG, BULKDENSITY,
                BULKWATMOLPER_VOXELVOL)
    _outputWeightedGrid(grid, grid_dir.val, grid_count_cutoff.val, BULKORI, BULKFORCES, BULKTORQUES, BULKNRG,
                        BULKDENSITY, BULKWATMOLPER_VOXELVOL)
    _outputGridParams(grid, grid_dir.val, grid_count_cutoff.val, BULKORI, BULKFORCES, BULKTORQUES, BULKNRG, BULKDENSITY,
                      BULKWATMOLPER_VOXELVOL)

    ta = time.time()
    if BENCHMARK:
        print("time to generate cell files %5.3f seconds " % (ta - tb))

@resolveParameters
def regionproperties():
    """
    regionproperties works out the excess thermodynamic properties of water molecules that were contained in
    a user defined region of space."""

    BULKORI = waterProps[watermodel.val]['BULKORI']
    BULKFORCES = waterProps[watermodel.val]['BULKFORCES']
    BULKTORQUES = waterProps[watermodel.val]['BULKTORQUES']
    BULKNRG = waterProps[watermodel.val]['BULKNRG']
    BULKDENSITY = waterProps[watermodel.val]['BULKDENSITY']

    print (regionfile.val)
    name = os.path.split(regionfile.val)[-1]
    root = name.replace(".region", "")
    #print (name, root)
    #sys.exit(-1)

    tb = time.time()

    # Load the grid
    grid = _loadGrid(gridforces.val)
    # Load the region
    region = _loadRegion(regionfile.val)
    # Get region properties
    region_props = _getProperties(region, grid['points'], BULKORI, BULKFORCES, BULKTORQUES, BULKNRG, BULKDENSITY)
    # Write region properties
    _writeRegionProperties(region_props, root)
    # Write a pdb to visualize the region
    _writeRegionPDB(region, grid, root, BULKDENSITY)

    ta = time.time()
    if benchmark.val:
        print("time to generate region properties %5.3f seconds " % (ta - tb))


@resolveParameters
def clustergrids():
    """
    clustergrids clusters grid according to distance cutoff for neighbours and density threshold as a
    multiple of bulk density"""

    BULKORI = waterProps[watermodel.val]['BULKORI']
    BULKFORCES = waterProps[watermodel.val]['BULKFORCES']
    BULKTORQUES = waterProps[watermodel.val]['BULKTORQUES']
    BULKNRG = waterProps[watermodel.val]['BULKNRG']
    BULKDENSITY = waterProps[watermodel.val]['BULKDENSITY']

    tb = time.time()
    # Load grid
    grid = _loadGrid(gridforces.val)

    #### get origin and grid step
    # the grid dimensions/origin/step
    #nx, ny, nz, ox, oy, oz, step = buffer[1].split()
    #grid['params'] = [int(nx), int(ny), int(nz), float(ox), float(oy), float(oz), float(step)]

    ## origin
    origin = [grid['params'][3], grid['params'][4], grid['params'][5]]
    step = grid['params'][6]

    # Sort grid points according to their density
    densities = {}
    for p in grid['points'].keys():
        densities[p] = grid['points'][p]['density']
    sorted_points = sorted(densities.items(), key=operator.itemgetter(1), reverse=True)
    # Group grid points into sites
    sites = _cluster(sorted_points, grid, lowt.val, neighcut.val, BULKDENSITY)
    print ("Clustering analysis generated %s sites " % len(sites.keys()))
    # Output sites in region/pdb format.
    _outputSites(sites, sites_dir.val)

    ta = time.time()
    if benchmark.val:
        print("time to cluster %5.3f seconds " % (ta - tb))

    ####  move into the clustered sites folder.
    os.chdir(sites_dir.val)

    #generate data on all sites
    for files in os.listdir("."):
        if files.endswith(".region") and files.startswith("site"):
            _regionproperties_clu("%s" % files, "../%s" % gridforces.val, BULKORI, BULKFORCES, BULKTORQUES, BULKNRG, BULKDENSITY)

    ### Typical output of regionproperties
    '''
    #Excess Thermodynamic properties of water within region

    Delta_G =  -9.1684 kcal/mol
    Delta_H = -10.4131 kcal/mol
    -T*Delta_S =   1.2447 kcal/mol
    Delta_Hw =   0.6946 kcal/mol
    Delta_Hx = -11.1077 kcal/mol
    -T*Delta_Svib =   0.4485 kcal/mol
    -T*Delta_Slib =   0.4161 kcal/mol
    -T*Delta_Sori =   0.3800 kcal/mol
    Relative_Density =   9.1534
    Average_number_of_water_molecules =   0.9899
    Volume =   3.2500 A^3


    ### generate CSV file with each name grid point at the right followed by components found in region
    ### Dictionary for CSV file conversion
### script to generate csv for the regionproperties

### Typical output of regionproperties
    '''

    ### generate CSV file with each name grid point at the right followed by components found in region
    ### Dictionary for CSV file conversion
    dict_centroid = {}
    parmlist = []
    for file in os.listdir("."):
        if file.endswith("_properties.dat") and file.startswith("site") and file != "site_allcenters_properties.dat":
            stream_cen = open("%s" % file, "r")
            buffer = stream_cen.readlines()
            stream_cen.close()
            site_no = file.lstrip("site").rstrip("_properties.dat")
            #print site_no
            dict_centroid[site_no] = {}
            for line in buffer:
                #print line
                if line.startswith("#"):
                    continue
                ### get parameters into dictionary
                if line.strip():
                    elems = line.split()
                    #print elems
                    dict_centroid[site_no][str(elems[0])] = float(elems[2])
                    #print dict_centroid
                    #sys.exit(-1)

                ### using grid indice get cartesian coordinates of the grid points
                stream = open("site%s.region" % site_no, 'r')
                buffer = stream.readlines()
                stream.close()
                for line in buffer:
                    if line.startswith("#"):
                        elems = line.split()
                        centerindex = int(elems[1])


                #cartx = originx + coords[0] * step
                #carty = originy + coords[1] * step
                #cartz = originz + coords[2] * step

                dict_centroid[site_no]['x'] = origin[0] + float( grid['points'][centerindex]['coords'][0] ) * step
                dict_centroid[site_no]['y'] = origin[1] +  float( grid['points'][centerindex]['coords'][1] ) * step
                dict_centroid[site_no]['z'] = origin[2] + float( grid['points'][centerindex]['coords'][2] ) * step


    writer = csv.writer(open('centroids.csv', 'w'))
    ### First make header of the csv file
    parmlist = ["cluster_no", "x", "y", "z", "Delta_G", "Delta_H", "-T*Delta_S", "Delta_Hw", "Delta_Hx", "-T*Delta_Svib", "-T*Delta_Slib", "-T*Delta_Sori", "Relative_Density", "Average_number_of_water_molecules", "Volume"]
    #print parmlist
    writer.writerow(parmlist)

    ### now populate each centroid(cluster_no) with each of its parameters
    for centroid in sorted(list(dict_centroid.keys()), key=float):
        rowdata = ["%s" % centroid]
        for parm in parmlist[1:]:
            rowdata += [float(dict_centroid[centroid][parm])]
        writer.writerow(rowdata)

@resolveParameters
def subgrids():
    """
    subgrids subtract grid in dx format to get differences where grid1.dx-grid2.dx=diff.dx.  Grids must be of
    identical dimensions and grid density"""

    tb = time.time()

    header = _loadheader(gridf.val)

    ## create file
    stream = open(diffdx.val, 'w')

    grid1 = _loadgrid_dx(gridf.val)
    grid2 = _loadgrid_dx(gridl.val)

    sub_grid = _subtractgrids(grid1, grid2)

    for line in header:
        stream.write("%s \n" % line.rstrip())
    for x in range(0,len(sub_grid),3):
        val0 = sub_grid[x]
        #print x, len(avg_grid)
        if (x + 1 ) >= len(sub_grid):
           str = "%8.5f" % val0
        elif ( x + 2 ) >= len(sub_grid):
            val1 = sub_grid[x+1]
            str = "%8.5f %8.5f" % (val0, val1)
        else:
            val1 = sub_grid[x+1]
            val2 = sub_grid[x+2]
            str = "%8.5f %8.5f %8.5f \n" % (val0, val1, val2)
        stream.write("%s" % str)
    stream.close()
    ta = time.time()
    if benchmark.val:
        print("time to subtract grids %5.3f seconds " % (ta - tb))

@resolveParameters
def avggrids():
    """
    avggrids averages grids in dx format.  Grids must be of identical dimensions and grid density """
    grid_files = avggridfiles.val

    header = _loadheader(grid_files[0])

    ## create file
    stream = open(avgdx.val, 'w')


    grids = []
    for grid_file in grid_files:
        grid = _loadgrid_dx(grid_file)
        grids.append(grid)

    avg_grid = _averagegrids(grids)

    for line in header:
        stream.write("%s \n" % line.rstrip())
    for x in range(0,len(avg_grid),3):
        val0 = avg_grid[x]
    #print x, len(avg_grid)
        if (x + 1 ) >= len(avg_grid):
            str = "%8.5f" % val0
        elif ( x + 2 ) >= len(avg_grid):
            val1 = avg_grid[x+1]
            str = "%8.5f %8.5f" % (val0, val1)
        else:
            val1 = avg_grid[x+1]
            val2 = avg_grid[x+2]
            str = "%8.5f %8.5f %8.5f \n" % (val0, val1, val2)
        stream.write("%s" % str)

    ta = time.time()
    if benchmark.val:
        print("time to average grids %5.3f seconds " % (ta - tb))
