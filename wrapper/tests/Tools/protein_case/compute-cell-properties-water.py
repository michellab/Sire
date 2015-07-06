#!/usr/bin/python
#
# Cell theory analysis of an MD frame
#
# Julien Michel and George Gerogiokas Dec 2012
#
 
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
import os, sys

import struct

import time
import numpy as np

from MDAnalysis import Universe

######  CALCULATION PARAMETERS ##
combining_rules = "arithmetic"
temperature = 298 * kelvin
cutoff = 10.0 * angstrom # Using a reaction field
## the dielectric for the reaction field
rfdielectric=78.3
# Name of water residues in input
Water_residue_names = ["T4P","T3P","Wat","HOH","WAT"]
# The water model used in the simulation
waterModel = "TIP4PEW-SireOpenMM"
# the name of the gridin file
gridin_file = "grid.in"
# Name of output folder for cell data
cell_dir = "./cell"
# Whether or not to overwrite existing output
#Overwrite = False
# Amber/GAFF atom types classified as donors or acceptors
polarH = ["HW","H","HO", "C3","C5"]
acceptors = ["OW","O", "OH", "Cl", "N", "N3", "N2", "Na", "C2", "C6" ]

####### FUNCTIONS  ###############

bulk = { "TIP4PEW-SireOpenMM": 4.88,
	"TIP4PEW-RH" : 4.67 }

# The bulk coordination number of water
BULKCN = bulk[waterModel]
# The cutoff to compute coordination numbers of water to other waters
CUTOFF_CN_WAT = 3.4 
# The curoff to use for coordination number of water to solute acceptors
CUTOFF_CN_SOL = 3.4 # This is for Cl, but do we need different numbers for each solute acceptor?
# GRID BUFFER. When parsing input, keep track of position of donors/acceptors that are within grid dimension + GRIDBUF
# This will cause an incorrect calculation of hydrogen-bonds and coordination numbers for waters that are at the edge of the extended 
# grid, since they will miss some neighbors. It is important that GRIDBUF is "big enough" then. 
GRIDBUFF = 5.0 

NONPOLARFLAG = 0
POLARHFLAG = 1
ACCEPTORFLAG = 2

# Atom with masses ( in g.mol-1) above this are considered "HEAVY" 
LIGHT = 1.10


CUTOFF_CN_WAT2 = CUTOFF_CN_WAT**2
CUTOFF_CN_SOL2 = CUTOFF_CN_SOL**2

BENCHMARK=False

def createSystem(molecules):

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

def setupForcefields(system, space):

    # print "Creating force fields... "

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
                        CHARMMSwitchingFunction(cutoff) )
    system.setProperty( "useReactionField", VariantProperty(True) )
    system.setProperty( "reactionFieldDielectric", VariantProperty(rfdielectric) )
    system.setProperty( "combiningRules", VariantProperty(combining_rules) )

    #total_nrg = gridwaterff.components().total() + gridwater_otherff.components().total() #+\
    #    # otherff.components().total()

    gridwater_nrg =  gridwaterff.components().total() + gridwater_otherwaterff.components().total() +\
                     gridwater_solutesff.components().total()

    # Maybe include otherwater_solutesff.components.total() ?? 

    gridwater_solutes_nrg = gridwater_solutesff.components().total()
    #e_total = system.totalComponent()

    e_gridwater = Symbol("E_{gridwater}")
 
    e_gridwater_solutes = Symbol("E_{gridwater:solutes}")

    #system.setComponent( e_total, total_nrg )

    system.setComponent( e_gridwater, gridwater_nrg )

    system.setComponent( e_gridwater_solutes, gridwater_solutes_nrg )

    return system

def assignDonorsAcceptors( system ):
    
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

def npdistance2(x0, x1, dimensions):
    delta = np.abs(x0 - x1)
    delta = np.where(delta > 0.5 * dimensions, dimensions - delta, delta)
    return (delta ** 2).sum(axis=-1)

def fastassignHbonds( space, polarH_coords, acceptor_coords, alternate_acceptor_coords, atomsDAtype, grid):
    # the key-value mark a h-bond between a polarH and an acceptor
    h_bonds = {}

    polarH_atoms = polarH_coords.keys()
    polarH_atoms.sort()
    ND = len(polarH_atoms)
    
    acceptor_atoms = acceptor_coords.keys()
    acceptor_atoms.sort()
    NA = len(acceptor_atoms)

    for polarH_at in polarH_atoms:
        h_bonds[polarH_at] = []
        
    acoords = np.zeros([NA,3])
    ac = np.zeros(NA)

    dimensions = np.array( [space.dimensions().x(), space.dimensions().y(), space.dimensions().z()] )

    for i in xrange(NA):
        acceptor = acceptor_atoms[i]
        h_bonds[acceptor] = []
        try:
            A_coords = alternate_acceptor_coords[acceptor][1]
            A_charge = atomsDAtype[ alternate_acceptor_coords[acceptor][0] ][1]
            # print "Using alternative coords / charge %s %s " % (A_coords, A_charge)
        except KeyError:
            A_coords = acceptor_coords[acceptor]
            A_charge = atomsDAtype[acceptor][1]
            # print "Using normal coords / charge %s %s " % (A_coords, A_charge)        
        acoords[i] = [ A_coords[0], A_coords[1], A_coords[2] ]
        ac[i] = A_charge

    for polarH_at in polarH_atoms:
        H_coords = polarH_coords[polarH_at]
        point = np.array( [ H_coords[0], H_coords[1], H_coords[2] ] )
    
        connected_acceptors = atomsDAtype[polarH_at][2]
        #print polarH_at
        #print connected_acceptors
     
        hforce = np.empty(NA)
        d2 = npdistance2(point, acoords, dimensions)
       
        #print d2
        hforce = ac / d2
        # Have to mask acceptor(s) bonded to polarH_at
        for acc in connected_acceptors:
            try:
                acceptor_index = acceptor_atoms.index(acc)
                hforce[acceptor_index] = +99
                #print acceptor_index
            except ValueError:
                # This could happen if acceptor was not in the buffered grid
                continue
       
        # Then get index from hforce.argmax()
        hmax = hforce.argmin()
        #print hmax
        hbonded = acceptor_atoms[hmax]
        #print polarH_at, hbonded
        #return h_bonds
        #print hbonded
        #sys.exit(-1)
        h_bonds[polarH_at].append(hbonded)
        h_bonds[hbonded].append(polarH_at)

    return h_bonds

def writeCellDataBin(cell_dir, data, frame_number, volume):

    keys = data.keys()

    keys.sort()

    bincell = struct.pack('<1i1f',frame_number,volume)
    for k in keys:
        print k, data[k]['waterflag'],data[k]['energies'][0], data[k]['solenergies'][0]
        sys.exit(-1)

        bincell +=struct.pack('<2i12f', k, data[k]['waterflag'],data[k]['energies'][0], data[k]['solenergies'][0], data[k]['forces'][0],data[k]['forces'][1],data[k]['forces'][2],data[k]['torques'][0],data[k]['torques'][1],data[k]['torques'][2], data[k]['o_coords'][0],data[k]['o_coords'][1],data[k]['o_coords'][2], data[k]['ori'])
        #strcell = "%8d %8d %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f" % (k, data[k]['waterflag'],data[k]['energies'][0], data[k]['forces'][0],data[k]['forces'][1],data[k]['forces'][2],data[k]['torques'][0],data[k]['torques'][1],data[k]['torques'][2], data[k]['o_coords'][0],data[k]['o_coords'][1],data[k]['o_coords'][2], data[k]['ori'])
        #strcell = "%8d %8.5f %8.5f %s %8.5f " % ( k, data[k]['o_coords'][1], len(data[k]['cn']['bulk'])+len(data[k]['cn']['fsw'])+len(data[k]['cn']['solute']), data[k]['cn'], data[k]['ori'] )
        #print strcell
    stream = open("%s/cell_results_%016d.bin" % (cell_dir,frame_number) ,"wd")
    stream.write(bincell)
    stream.close()

def updateSystemfromDCDgrid(system, dcd_traj, atomsDAtype, grid):

    dcd_coordinates = dcd_traj.atoms.coordinates()
    dcd_dimensions = dcd_traj.dimensions

    dcd_box_x = dcd_dimensions[0].tolist()
    dcd_box_y = dcd_dimensions[1].tolist()
    dcd_box_z = dcd_dimensions[2].tolist()

    dcd_natoms = dcd_traj.atoms.numberOfAtoms()
    
    #print "natoms", dcd_natoms, dcd_box_x, dcd_box_y, dcd_box_z
    
    newmols_coords = {}

    dcd_index = 0
    mol_index = 0

    molnums = system.molNums()
    molnums.sort()

    #print dcd_coordinates[0]
    
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

    #ta = time.time()
    #if BENCHMARK:
    #    print "Time to create new coordinate sets %5.3f s"  % (ta-tb)
    #tb = ta

    if dcd_natoms != dcd_index:
        print "The number of atoms in the system is not equal to the number of atoms in the DCD file ! Aborting."
        sys.exit(-1)

    #grid_box = AABox( Vector(grid['center']), Vector( (grid['max'][0]-grid['min'][0])/2.,
    #                                                  (grid['max'][1]-grid['min'][1])/2., 
    #                                                  (grid['max'][2]-grid['min'][2])/2. ) 

    #print grid['min']
    #print grid['max']
    #grid_box = PeriodicBox( Vector( grid['min'][0], grid['min'][1], grid['min'][2] ),
    #                        Vector( grid['max'][0], grid['max'][1], grid['max'][2] ) ) 

    #grid_box_min =  grid_box.minCoords( Vector( grid['center'] ) )
    #grid_box_max = grid_box.maxCoords( Vector( grid['center'] ) )
   
    grid_box_min = Vector(grid['min'][0], grid['min'][1], grid['min'][2] )
    grid_box_max = Vector(grid['max'][0], grid['max'][1], grid['max'][2] )

    #print dcd_dimensions

    #print grid_box_min
    #print grid_box_max
    #sys.exit(-1)
    #grid_buffer_box =  AABox( Vector(grid['center']), Vector( 2*GRIDBUFF+(grid['max'][0]-grid['min'][0])/2.,
    #                                                          2*GRIDBUFF+(grid['max'][1]-grid['min'][1])/2., 
    #                                                          2*GRIDBUFF+(grid['max'][2]-grid['min'][2])/2. ) 
    #grid_buffer_box_min = grid_buffer_box.minCoords()
    #grid_buffer_box_max = grid_buffer_box.maxCoords()

    grid_buffer_box_min = Vector(grid['min'][0]-GRIDBUFF, grid['min'][1]-GRIDBUFF, grid['min'][2]-GRIDBUFF )
    grid_buffer_box_max = Vector(grid['max'][0]+GRIDBUFF, grid['max'][1]+GRIDBUFF, grid['max'][2]+GRIDBUFF )

    #print grid_buffer_box_min
    #print grid_buffer_box_max

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
        #system.update(mol)

        if mol.residues()[0].name() == ResName('WAT'):
            # If water. check if oxygen atom within grid
            #wat_box = mol.property("coordinates")[CGIdx(0)].aaBox()
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
                #new_other.add(mol)
            else:
                new_otherwater.add(mol)
                #new_gridwater.add(mol)
        else:
            new_solutes.add(mol)

        # Update dictionnaries of donors/acceptors/alternate_acceptors and gridwater group based on new cartesian coordinates 
        molatoms = mol.atoms()
        molnatoms = mol.nAtoms()

        for x in xrange(molnatoms):
            atom = molatoms[x]
            at_num = atom.number().value()
            atom_coord = atom.property("coordinates")

#            if ( atom_coord.x() < grid_buffer_box_min.x() or
#                 atom_coord.y() < grid_buffer_box_min.y() or
#                 atom_coord.z() < grid_buffer_box_min.z() or
#                 atom_coord.x() > grid_buffer_box_max.x() or
#                 atom_coord.y() > grid_buffer_box_max.y() or
#                 atom_coord.z() > grid_buffer_box_max.z() ):
#                #print "Skip atom not in grid_buffer"
#                #print atom_coord
#                continue

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

            #print atom_coord, within_x, within_y, within_z

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

    #print dcd_dimensions

    #print grid_box_min
    #print grid_box_max

    #print grid_buffer_box_min
    #print grid_buffer_box_max

    #print acceptor_coords.keys()
    system.update(changedmols)

    #print system.groups()
    # print "Doing removeAllMoleculeGroups"
    system.removeAllMoleculeGroups()
    # print "Doing removeAllForceFields"
    system.removeAllForceFields()
    #print "Going to remove gridwater"
    #system.remove(MGName("gridwater"))
    #print "going to remove other"
    #system.remove(MGName("other")) 

    #print system.groups()

    new_all = MoleculeGroup("all", changedmols)
    # print "Adding new_all"
    system.add(new_all)
    # print "Adding new_gridwater"
    system.add(new_gridwater)
    # print "Adding new_other"
    system.add(new_otherwater)
    system.add(new_solutes)
    #print system.groups()

    #sys.exit(-1)

    if  dcd_traj.trajectory.periodic: 
        space = PeriodicBox(Vector( dcd_dimensions[0].tolist(),
                                     dcd_dimensions[1].tolist(), 
                                     dcd_dimensions[2].tolist() ) )
    else:
        space = Cartesian()

    # And a new space object
    #space = system.property("space")
    
    #if space.isPeriodic() and dcd_traj.trajectory.periodic: 
    #    space.setDimensions( Vector( dcd_dimensions[0].tolist(),
    #                                 dcd_dimensions[1].tolist(), 
    #                                 dcd_dimensions[2].tolist() ) )
    #    system.setProperty("space",space)
    #elif (not space.isPeriodic and not dcd_traj.trajectory.periodic):
    #    # Do nothing as no dimensions
    #    pass
    #else:
    #    print "Problem !! the Sire system and dcd files do not have the same periodicity! Aborting."
    #    sys.exit(-1)
    system = setupForcefields(system, space)
    

    #ta = time.time()
    #if BENCHMARK:
    #    print "Time to update sys with new coordinates %5.3f s"  % (ta-tb)
    #tb = ta

    #PDB().write(system[MGName("gridwater")], "debug.pdb" )

    return system, polarH_coords, acceptor_coords, alternate_acceptor_coords

def fastercoordinationNumbers( space, acceptor_coords, h_bonds, atomsDAtype, firstshellwaters, grid):
    
    grid_minx = grid['min'][0]
    grid_miny = grid['min'][1]
    grid_minz = grid['min'][2]
    grid_maxx = grid['max'][0]
    grid_maxy = grid['max'][1]
    grid_maxz = grid['max'][2] 

    water_cn = {}
    acceptor_atoms = acceptor_coords.keys()
    nacceptors = len(acceptor_atoms)
    
    acoords = np.zeros([nacceptors,3])
    for i in xrange(nacceptors):
        iat = acceptor_atoms[i]
        acoords[i] = [ acceptor_coords[iat][0], acceptor_coords[iat][1], acceptor_coords[iat][2] ]
    
    dimensions = np.array( [space.dimensions().x(), space.dimensions().y(), space.dimensions().z()] )

    for i in xrange(nacceptors-1):
        point = acoords[i]
        iat = acceptor_atoms[i]
        iat_inwat = atomsDAtype[iat][4]

        iat_ingrid = True
        if ( acoords[i][0] < grid_minx or 
             acoords[i][1] < grid_miny or
             acoords[i][2] < grid_minz or 
             acoords[i][0] > grid_maxx or
             acoords[i][1] > grid_maxy or 
             acoords[i][2] > grid_maxz ):
            # iat is outside of the grid...
            iat_ingrid = False

        #print "Looking at %s ingrid %s "% (iat, iat_ingrid)
        d2 = npdistance2( point, acoords[i+1:], dimensions)
        # Select acceptors within cutoff
        check = np.where(d2 < CUTOFF_CN_WAT2, 1, 0)
        #print len(check)
        coordinated = check.nonzero()[0]
        #print coordinated
   
        # Need to consider case where there are no coordinated particles 
        # but particle is in grid then water_cn has to be initialised

        imolnum = atomsDAtype[iat][3]

        if iat_ingrid and not water_cn.has_key(imolnum):
            water_cn[imolnum] = {}
            water_cn[imolnum]['fsw'] = []
            water_cn[imolnum]['bulk'] = []
            water_cn[imolnum]['solute'] = []

        for j in coordinated:
            shift = i + j + 1
            jat = acceptor_atoms[ shift ]
            #print jat
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

            #imolnum = atomsDAtype[iat][3]
            jmolnum = atomsDAtype[jat][3]

            #if (not iat_inwat) or (not jat_inwat):
            #    print "PAIR involving solute %s %s " % (imolnum, jmolnum)

            if iat_inwat and iat_ingrid:
                if jat_inwat:
                    if firstshellwaters.has_key(jmolnum):
                        cntype = 'fsw'
                    else:
                        cntype = 'bulk'
                else:
                    cntype = 'solute'

                if not water_cn.has_key(imolnum):
                    water_cn[imolnum] = {}
                    water_cn[imolnum]['fsw'] = []
                    water_cn[imolnum]['bulk'] = []
                    water_cn[imolnum]['solute'] = []
                water_cn[imolnum][cntype].append( jat )

            if jat_inwat and jat_ingrid:
                if iat_inwat:
                    if firstshellwaters.has_key(imolnum):
                        cntype = 'fsw'
                    else:
                        cntype = 'bulk'
                else:
                    cntype = 'solute'

                #print soluteinc, fswinc, bulkinc

                if not water_cn.has_key(jmolnum):
                    water_cn[jmolnum] = {}
                    water_cn[jmolnum]['fsw'] = []
                    water_cn[jmolnum]['bulk'] = []
                    water_cn[jmolnum]['solute'] = []
                water_cn[jmolnum][cntype].append( iat )
        #sys.exit(-1)

    return water_cn

def orientationalNumbers( water_cn, firstshellwaters, h_bonds, mol_to_atoms, atomsDAtype ):
    water_ori = {}

    # For "normal" water
    #   just sum neighbors, and then convert to ori
    # For each first shell water
    # We know number of solutes number of first shell and number of bulk from water_cn...
    # We need to know if we are donating a h-bond to a solute atom, to another first shell water, to a bulk water (pwl_HB)
    # Then we can compute the effective number and then the orientational number

    waters = water_cn.keys()
 
    for wat in waters:
        #print wat
        if firstshellwaters.has_key(wat):
            # First shell water treatment
            molatoms = mol_to_atoms[wat]
     
            #print water_cn[wat]

            o_atomnum = None
            h_atomnums = []

            for atnum in molatoms:
                if atomsDAtype[atnum][0] == POLARHFLAG:
                    h_atomnums.append(atnum)
                elif atomsDAtype[atnum][0] == ACCEPTORFLAG:
                    o_atomnum = atnum 
            #print o_atomnum
            #print h_atomnums
 
            hbonded = []
            # the atomic numbers of the atoms that donate a h-bond to wat oxygen
            h_accepteds = h_bonds[o_atomnum]
            # now find to which heavy atom they are bonded
            for h in h_accepteds:
                hbonded += atomsDAtype[h][2]
            # the atomic numbers of the atoms that accepte a h-bond from wat  
            hbonded +=  h_bonds[ h_atomnums[0] ] + h_bonds[ h_atomnums[1] ] 

            #print hbonded
            # For each type of cn member, find out how many are currently h-bonded to that water
            nw_solute = len(water_cn[wat]['solute'])
            nw_solute_hb = 0
            phb_solute = 0
            # for each solute, check if accepting a h-bond from that wat
            for sol in water_cn[wat]['solute']:
                if sol in hbonded:
                    #print "hbonded with solute atom in coordination sphere"
                    nw_solute_hb += 1    
            if nw_solute > 0:
                phb_solute = float(nw_solute_hb) / nw_solute
            
            # now the fsw
            nw_fsw = len(water_cn[wat]['fsw'])
            nw_fsw_hb = 0
            phb_fsw = 0
            for fsw in water_cn[wat]['fsw']:
                if fsw in hbonded:
                    #print "hbonded with fsw atom in coordination sphere"
                    nw_fsw_hb += 1
                    #sys.exit(-1)
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

            #print nw_solute, nw_solute_hb
            #print nw_fsw, nw_fsw_hb
            #print nw_bulk, nw_bulk_hb
            #print phb_solute, phb_fsw, phb_bulk
            phb_max = max ( max ( phb_solute, phb_fsw ) , phb_bulk )
            if phb_max < small:
                neff = 0
            else:
                neff = nw_solute_hb / phb_max + nw_fsw_hb / phb_max + nw_bulk_hb / phb_max
            # And ori
            nw = BULKCN
            #nw = nw_solute + nw_fsw + nw_bulk
            ori = ( ( neff * ( neff-1. ) ) / 2. ) * ( ( nw - 2. ) / nw )**(2-phb_solute)
            #print nw
            #print neff
            #print "fsw ", nw, neff, ori
            #if ori < 1 : 
            #    ori = 1.0
        else:
            # bulk treatment
            nw = len(water_cn[wat]['bulk']) + len(water_cn[wat]['solute']) + len(water_cn[wat]['fsw'])
            #print nw
            if nw < 1:
                ori = 0.0
            else:
                ori = ( (nw * ( nw-1.) ) / 2. ) * ( ( nw - 2. ) / nw ) ** 2
            #print "bulk ", nw, ori
            #if ori < 1.0:
            #    ori = 1.0
            #print ori
        water_ori[wat] = ori
        #sys.exit(-1)

    #sys.exit(-1)
    return water_ori

def fastdefineFSW(space, acceptor_coords, atomsDAtype, grid ):
    
    cutoff_cn2 = CUTOFF_CN_SOL2 # Angstroms
    
    acceptor_atoms = acceptor_coords.keys()

    nacceptors = len(acceptor_atoms)

    waters_firstshell = {}

    # Loop over acceptor_atoms and split in solute/waters
    water_acceptors = []
    solute_acceptors = []
    
    nwat = 0
    nsol = 0
    for i in xrange(nacceptors):
        iat = acceptor_atoms[i] 	
        iat_inwat = atomsDAtype[iat][4]
        if iat_inwat:
            water_acceptors.append(iat)
            nwat += 1
        else:
            solute_acceptors.append(iat)
            nsol += 1
     
    water_coords = np.zeros([nwat,3])
    for i in xrange(nwat):
        iat = water_acceptors[i]
        water_coords[i] = [ acceptor_coords[iat][0], acceptor_coords[iat][1], acceptor_coords[iat][2] ]

    solute_coords = np.zeros([nsol,3])
    for i in xrange(nsol):
        iat = solute_acceptors[i]
        solute_coords[i] = [ acceptor_coords[iat][0], acceptor_coords[iat][1], acceptor_coords[iat][2] ]

    #print nwat
    #print nsol
    #print len(water_acceptors)
    #print len(solute_acceptors)

    dimensions = np.array( [space.dimensions().x(), space.dimensions().y(), space.dimensions().z()] )

    fsw = np.zeros(nwat)

    for i in xrange(nsol):
        point = solute_coords[i]
        
        d2 = npdistance2(point, water_coords, dimensions)
        
        contact = np.where(d2 < cutoff_cn2, 1, 0)
        fsw += contact
        #print d2
        #sys.exit(-1)

    #print fsw
    fsw_indices = fsw.nonzero()[0]
    for idx in fsw_indices:
        iat = water_acceptors[idx]
        iatmolnum = atomsDAtype[iat][3]
        #print iatmolnum
        waters_firstshell[iatmolnum] = True

    return waters_firstshell

def readGrid( grid_file ):
    
    stream = open( grid_file , 'r' )
    buffer = stream.readlines()
    elems = buffer[1].split()
    cx = float( elems[0] )
    cy = float( elems[1] )
    cz = float( elems[2] )
    plusdx = float( elems[3] )
    mindx = float( elems[4] )
    plusdy = float( elems[5] ) 
    mindy = float( elems[6] )
    plusdz = float( elems[7] )
    mindz = float( elems[8] )
    
    grid = {}
    grid['center'] = ( cx, cy, cz )
    grid['step'] = ( plusdx, mindx, plusdy, mindy, plusdz, mindz )
    grid['min'] = ( cx - mindx, cy - mindy, cz - mindz )
    grid['max'] = ( cx + plusdx, cy + plusdy, cz + plusdz )

    vol = ( plusdx + mindx ) * ( plusdy + mindy ) * ( plusdz + mindz )
    grid['volume'] = vol

    return grid

######## MAIN SCRIPT #############

if __name__ == "__main__":
    try:
        top_file = sys.argv[1]
        crd_file = sys.argv[2]
        dcd_file = sys.argv[3]
        grid_file = sys.argv[4]
    except IndexError:
        print "Usage is %s top_file crd_file dcd_file grid_file [start_frame/end_frame]" % ( sys.argv[0] )
        sys.exit(-1)
    #print j

    # By default we process every frame
    start_frame = -1
    end_frame = 1e10

    try:
        start_frame = sys.argv[5]
        end_frame = sys.argv[6]
        start_frame = int(start_frame)
        end_frame = int(end_frame)
    except IndexError:
        pass
    
    if not os.path.exists(cell_dir):
        cmd = "mkdir %s" % cell_dir
        os.system(cmd)

    # First create a Sire system using input top/crd files
    
    amber = Amber()
    molecules, space = amber.readCrdTop(crd_file, top_file)
    system = createSystem(molecules)
    #system = setupForcefields(system, space)

    # Assign h-bond donors/acceptors
    atomsDAtype, mol_to_atoms = assignDonorsAcceptors( system )

    # Read the grid boundaries
    grid = readGrid(grid_file)   

    #sys.exit(-1)

    # Now load the dcd trajectory
    dcd_traj = Universe(top_file, dcd_file)
    nframes = len(dcd_traj.trajectory)

    #e_total = Symbol("E_{total}")
    e_gridwater = Symbol("E_{gridwater}")
    e_gridwater_solutes = Symbol("E_{gridwater:solutes}")

    # For each frame in the dcd trajectory...
    for framenumber in range(0, nframes):
        # Time in picoseconds
        frametime = round(dcd_traj.trajectory.time,3)
        #print dcd_traj.trajectory.time
        if framenumber < start_frame: 
            if (framenumber+1) < nframes:
                dcd_traj.trajectory.next()
                #print "Below start chunk. Continue. "
                continue
        elif framenumber > end_frame:
            #print "Above end chunk. Break."
            break
        #sys.exit(-1)
        print "Frame number %s start_f %sf end_f %s " % (framenumber, start_frame, end_frame) 

        tb = time.time()
        print "Processing frame %s at time %16.3f ..." % (framenumber, frametime)
        system, polarH_coords, acceptor_coords, alternate_acceptor_coords = updateSystemfromDCDgrid(system, dcd_traj, atomsDAtype, grid)

        # Dump a PDB of the system if this is the first frame, to facilitate visualisation of cell/grid/sites
        if framenumber == 0:
            PDB().write(system[MGName("all")], "frame-0.pdb" )

        ta = time.time()
        if BENCHMARK:
            print "time to update coords %5.3f seconds" % (ta - tb) 
        tb = ta
        #sys.exit(-1)

        # Initialise energies, this appears necessary to initialise properly 
        # LJPair database for force calculations
        system.energy()
        #print system.energy()
        #print system.energies()
        #sys.exit(-1)

        # Use a custom RBWorkspace to compute forces and torques
        all = system[MGName("all")]
        #gridwater = system[MGName("gridwater")]
        #wspace = RBWorkspaceJM(gridwater)
        wspace =  RBWorkspaceJM(all)
        wspace.setSystem(system)
 
        #print "About to call calculate forces"
        #wspace.calculateForces(e_total)
        wspace.calculateForces(e_gridwater)
        #sys.exit(-1)

        bead_forces = wspace.beadForcesArray()
        bead_torques = wspace.beadTorquesArray()
        bead_energies = wspace.beadEnergiesArray()       

        #print len(bead_forces)
        #print bead_forces
        #sys.exit(-1)
        #print bead_torques[-1]
        #print bead_energies[-1]

        # Now the solute energy terms
        wspace.calculateForces(e_gridwater_solutes)
        beadsolute_energies = wspace.beadEnergiesArray()

	ta = time.time()
        if BENCHMARK:
            print "time to get energies and forces/torques %5.3f seconds " % (ta - tb)
        tb = ta

        #sys.exit(-1)

        space = system.property("space")
        
        # detect waters that are first shell to solute acceptors
        firstshellwaters = fastdefineFSW(space, acceptor_coords, atomsDAtype, grid)

        ta = time.time()
        if BENCHMARK:
            print "time to detect first shell waters %5.3f seconds " % (ta - tb)
        tb = ta

        #sys.exit(-1)

        # h-bonds analysis
        h_bonds = fastassignHbonds( space, polarH_coords, acceptor_coords, alternate_acceptor_coords, atomsDAtype, grid)

        ta = time.time()
        if BENCHMARK:
            print "time to assign hbonds %5.3f seconds " % (ta - tb)
        tb = ta

        #sys.exit(-1)

        # a dictionnary of coordination numbers for each water's oxygens
        water_cn = fastercoordinationNumbers( space, acceptor_coords, h_bonds, atomsDAtype, firstshellwaters, grid)

        ta = time.time()
        if BENCHMARK:
            print "time to do fastCN %5.3f seconds " % (ta - tb)
        tb = ta

        #sys.exit(-1)
        #print water_cn

        water_ori = orientationalNumbers( water_cn, firstshellwaters, h_bonds, mol_to_atoms, atomsDAtype )
        #print len(water_cn)

        ta = time.time()
        if BENCHMARK:
            print "time to do orientational number %5.3f seconds " % (ta - tb)
        tb = ta

        # A complication is that the bead values are not in the same order 
        # between different snapshots. So let's store the forces/torques/energies 
        # values in dictionnaries keyed by the first residue number of 
        # each mol (unique)
        molecules = system[MGName("all")].molecules()
        gridwater = system[MGName("gridwater")]
        molnums = molecules.molNums()
        
        water_results = {}

        for x in range(0,len(molnums)):
            molnum = molnums[x]
            molecule = system.molecule(molnum).molecule()
            # Ignore if not gridwater
            if not gridwater.contains(molecule):
                #print "not in gridwater"
                continue

            # Ignore non water molecules
            nresidues = molecule.nResidues()
            for y in range(nresidues):   
                residue = molecule.residue(ResIdx(y))
                #print residue
                resname = residue.name().value()
                # This is a solute residue. Skip it. 
                if not resname in Water_residue_names:
                    #print "Solute residue - continue -"
                    #waterflag = 0
                    continue
                # This is a water molecule
                else:
                    #print "Water"
                    # Skip if not within grid volume
                    if not water_ori.has_key(molnum.value()):
                        #print "Not in grid"
                        continue
                    waterflag = 1
                    #print "Should be in grid"
                    
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

                #if waterflag:
                # Also get coordinates for O to track on grid
                oatom = molecule.select( AtomName("O") )
                mol_o_coords = oatom.property("coordinates")
                water_results[res_num]['o_coords'] = mol_o_coords
                water_results[res_num]['ori'] = water_ori[ molnum.value() ]
                water_results[res_num]['cn'] = water_cn[ molnum.value() ]
                #else:
                #    water_results[res_num]['o_coords'] = ( 0.0, 0.0 ,0.0 )
                #    water_results[res_num]['ori'] = 1.0

        #sys.exit(-1)
        # JM May 13 - store frame number rather than time
        #writeCellDataBin(cell_dir, water_results, frametime, grid['volume'] )
        print water_results[753]
        sys.exit(-1)
        writeCellDataBin(cell_dir, water_results, framenumber, grid['volume'] )

        ta = time.time()
        if BENCHMARK:
            print "time to write output %5.3f seconds " % (ta - tb)
        tb = ta

        # PDB().write(system[MGName("all")], "test-%s.pdb" % framenumber)

        # Advance to the next frame if not at last one

        if (framenumber+1) < nframes:
            dcd_traj.trajectory.next()
        # End of loop over DCD frame
         #sys.exit(-1)


