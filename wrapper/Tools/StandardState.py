#Script to calculate the free energy change
# upon removing a set of host-guest distance
# restraints and applying standard state conditions
# @authors: Stefano Bosisio and Julien Michel

import os,sys, random
import math
from math import pi, cos, sin
from Sire.Tools.OpenMMMD import *
from Sire.Tools import Parameter, resolveParameters
# Python dependencies

try:
    import mdtraj
except ImportError:
    print ("StandarState.py depends on a working install of the python module mdtraj. Please install mdtraj in your sire python.")
    sys.exit(-1)

try:
    import ast
except ImportError:
    print ("StandarState.py depends on a working install of the python module ast. Please install ast in your sire python.")
    sys.exit(-1)

try:
    import shutil
except ImportError:
    print ("StandarState.py depends on a working install of the python module shutil. Please install shutil in your sire python.")
    sys.exit(-1)

try:
    import numpy
except ImportError:
    print ("StandarState.py depends on a working install of the numpy module. Please install mdtraj in your sire python.")
    sys.exit(-1)

# Constant for conversion
NM_TO_ANG = 10.0


trajfile = Parameter("trajfile", "traj000000001.dcd",
                    """File name of the trajectory to process.""")

stepframe = Parameter("step_frame",1,"""The number of frames to step to between two succcessive evaluations.""")
simfile  = Parameter("simfile", "sim.cfg", """ Configuration file with distance restraints dictionary""")
topfile = Parameter("topfile", "SYSTEM.top",
                    """File name of the topology file containing the system to be simulated.""")

buff = Parameter("buff", 5.0, """Buffer to be added to coordinates to create domain of integration""")

verbose = Parameter("verbose", False, """Print debug output""")


def averageCoordinates(restr_dict):
    r"""This function computes the average x,y and z coordinates for host atoms
    Parameters
    ----------
    restr_dict : dictionary
                 restr_dict[lig_idx,host_idx] = ([req,D,K],[x0,y0,z0],[x1,y1,z1]...)
                 where x0 is the x coordinate at frame 0 and so on
    Returns
    ----------
    restr_dict : dictionary
                 restr_dict[lig_idx,host_idx] = ([req,D,K],[avgx,avgy,avgz])
                 where avgx,avgy and avgz are the average coordinates
    """
    #restr_dict[idx]=([req,K,D],[coords])
    #Calculation of the mean coordinate for every atoms
    for pairs in restr_dict:
        coords1 = restr_dict[pairs][1] #here the list of all the coords
        coords2 = restr_dict[pairs][2]
        coord1_avg = [ 0.0, 0.0, 0.0 ]
        coord2_avg = [ 0.0, 0.0, 0.0 ]
        for x in range(0,len(coords1)):
            val1 = coords1[x]
            val2 = coords2[x]
            coord1_avg[0] += (val1[0] - coord1_avg[0])/(x+1)
            coord1_avg[1] += (val1[1] - coord1_avg[1])/(x+1)
            coord1_avg[2] += (val1[2] - coord1_avg[2])/(x+1)
            coord2_avg[0] += (val2[0] - coord2_avg[0])/(x+1)
            coord2_avg[1] += (val2[1] - coord2_avg[1])/(x+1)
            coord2_avg[2] += (val2[2] - coord2_avg[2])/(x+1)

        #Substitution of values and reset each avg coords
        # Note that mdtraj coordinates are in nm, but potential parameters
        # in angstrom
        restr_dict[pairs][1] = [ coord1_avg[0] * NM_TO_ANG,
                                 coord1_avg[1] * NM_TO_ANG,
                                 coord1_avg[2] * NM_TO_ANG ]
        restr_dict[pairs][2] = [ coord2_avg[0] * NM_TO_ANG,
                                 coord2_avg[1] * NM_TO_ANG,
                                 coord2_avg[2] * NM_TO_ANG ]

    return restr_dict

def defineIntegrationDomain(restr_dict):
    r"""Definition of the integration domain starting from host coordiantes
    Parameters
    ----------
    restr_dict : dictionary
                 restr_dict[lig_idx,host_idx] = ([req,D,K],[avg_lig_x, avg_lig_y, avg_lig_z], [avg_hostx,avg_hosty,avg_hostz])
                 where avgx,avgy and avgz are the average coordinates of the ligand and host atoms defined by lig_idx and host_idx
    Returns
    ----------
    space : 3D array
            space = [(-x,-y,-z)(x,y,z)]
            where the coordinates define the minimum and maximum coordinates
            of the integration domain

    """
    max_x = -999999
    max_y = -999999
    max_z = -999999
    min_x = +999999
    min_y = +999999
    min_z = +999999

    for pairs in restr_dict:
        coord = restr_dict[pairs][2]
        for val in coords:
            if val[0] > max_x:
                max_x = val[0]
            if val[0] < min_x:
                min_x = val[0]
            if val[1] > max_y:
                max_y = val[1]
            if val[1] < min_y :
                min_y = val[1]
            if val[2] > max_z:
                max_z = val[2]
            if val[2] < min_z :
                min_z = val[2]
    #print(max_x,max_y,max_z,min_x,min_y,min_z)
    # Adding a buffer region
    max_x += buff.val
    max_y += buff.val
    max_z += buff.val
    min_x -= buff.val
    min_y -= buff.val
    min_z -= buff.val
    space = [(min_x,min_y, min_z), (max_x,max_y,max_z)]
    return space

def genOrientations(restr_dict, norientations=5):
    r"""Generates a set of orientations for guest atoms
    The coordinates are in a frame of reference centered on the COG
    of the guest atoms.
    Norient orientations are generated by multiplying the input coordinates
    by a series of rotation matrices, each corresponding to a rotation along
    Euler angles theta, phi and psi.

    Parameters
    ----------
    restr_dict : dictionary
                 restr_dict[lig_idx,host_idx] = ([req,D,K],[avg_lig_x, avg_lig_y, avg_lig_z], [avg_hostx,avg_hosty,avg_hostz])
                 where avgx,avgy and avgz are the average coordinates of the ligand and host atoms defined by lig_idx and host_idx
    Returns
    ----------
    orientations : array of norient**3 * n_guest atom 3D coordinates
    """
    # if only one guest atom then return [ [0.0, 0.0, 0.0] ]

    # First pass: work out COG of guest atoms
    guest_cog = [0.0, 0.0, 0.0]
    for pairs in restr_dict:
        guest_cog[0] += restr_dict[pairs][1][0]*(1/len(restr_dict))
        guest_cog[1] += restr_dict[pairs][1][1]*(1/len(restr_dict))
        guest_cog[2] += restr_dict[pairs][1][2]*(1/len(restr_dict))
    # Second pass: Subtract COG to get COG centered coordinates
    body = []
    for pairs in restr_dict:
        new_x = restr_dict[pairs][1][0] - guest_cog[0]
        new_y = restr_dict[pairs][1][1] - guest_cog[1]
        new_z = restr_dict[pairs][1][2] - guest_cog[2]
        body.append( [new_x, new_y, new_z] )

    import pdb; pdb.set_trace()
    # Now work out set of rotations along Euler Angles
    TWOPI = 2*math.pi
    for x in range(0,norientations):
        theta = (x*TWOPI)/norientations
        for y in range(0,norientations):
            phi = (y*TWOPI)/norientations
            for z in range(0,norientations):
                psi = (z*TWOPI)/norientations
                rot00 = cos(phi)*cos(psi)-cos(theta)*sin(phi)*sin(psi)
                rot10 = sin(phi)*cos(psi)+cos(theta)*cos(phi)*sin(psi)
                rot20 = sin(theta)*sin(psi)
                rot10 = -cos(phi)*sin(psi)-cos(theta)*sin(phi)*cos(psi)
                rot11 = -sin(phi)*sin(psi)+cos(theta)*cos(phi)*cos(psi)
                rot12 = sin(theta)*cos(phi)
                rot20 = sin(theta)*sin(phi)
                rot21 = -sin(theta)*cos(phi)
                rot22 = cos(theta)


@resolveParameters
def run():
    try:
        host = os.environ['HOSTNAME']
    except KeyError:
        host = "unknown"
    print("### Running Standard state correction calculation on %s ###" % host)

    if verbose.val:
        print("###================= Simulation Parameters=====================###")
        Parameter.printAll()
        print ("###===========================================================###\n")

    #Constants
    delta = 0.10
    delta_over_two = delta/2.0
    deltavol = delta*delta*delta
    ROT = 8 * pi**2
    NORIENT = 5
    deltarot = ROT/NORIENT
    kb = 0.001987# GET FROM SIRE
    T = 298 # SHOULD READ THIS FROM CFG FILE
    kbT = kb*T # GET FROM SIRE
    beta = 1/kbT # GET FROM SIRE

    Ztrans = 0.0
    Uavg = 0

    sim_dictionary = distance_restraints_dict.val
    if sim_dictionary == {}:
        print ("Error, no distance restraints dictionary was found in the supplied config file. Abort.")
        sys.exit(-1)
    #now create a dictionary in this way:
    #dict[pairs] = {[Req,D,K] [coords] [coords]...}
    restr_dict = {}
    #create a list of host indexes to be used with mdtraj for alignment
    hosts = []
    for pairs in sim_dictionary:
        req = sim_dictionary[pairs][0]
        K   = sim_dictionary[pairs][1]
        D   = sim_dictionary[pairs][2]
        # First entry are parameters, second and third entries for coordinates of first and second atom
        restr_dict[pairs]=[[req,K,D],[],[]]
        #idx = max(pairs)
        #hosts.append(idx)
    #FIXME: code assumes guest atoms have lower indices than host atoms
    # a more reliable algorithm could work out whether the atoms belong
    # to a guest residue?

    #load the trajectory
    start_frame = 1
    end_frame = 1000000000
    step_frame = stepframe.val
    #Check the extension of the topology file
    if ".top" in topfile.val:
        shutil.copy(topfile.val,"SYSTEM.prmtop")
        #top doesn't work with mdtraj
        top_file = "SYSTEM.prmtop"
    else:
        top_file = topfile.val

    print("Loading trajectory and topology files")
    #
    mdtraj_trajfile = mdtraj.load(trajfile.val,top=top_file)
    nframes = len(mdtraj_trajfile)
    if end_frame > (nframes - 1):
        end_frame = nframes - 1
    current_frame = start_frame

    #Aligning everything along the first frame
    # FIXME: Code fails if too few host restrained atoms.
    # Use all heavy atoms for alignment?
    selection = "not water and not resname 'Na+' and not resname 'Cl-' and mass > 1"
    align_indices = mdtraj_trajfile.topology.select(selection)


    print("Aligning frames along first frame of trajectory")
    #aligned_traj = mdtraj_trajfile.superpose(mdtraj_trajfile,0, atom_indices=hosts)
    aligned_traj = mdtraj_trajfile.superpose(mdtraj_trajfile,0, atom_indices=align_indices)

    print("Processing frames")
    while (current_frame <= end_frame):
        # FIXME: Save restrained guest atoms as well
        #now for each lig:host pairs append the host coordinates
        for pairs in restr_dict:
            #host_idx = max(pairs)
            idx1 = pairs[0]
            idx2 = pairs[1]
            coord1 = aligned_traj.xyz[current_frame,idx1,:].tolist()
            coord2 = aligned_traj.xyz[current_frame,idx2,:].tolist()
            restr_dict[pairs][1].append(coord1)
            restr_dict[pairs][2].append(coord2)

        current_frame += step_frame

    #now restr_dict has:
    #restr_dict[lig,host]=[ [req,K,D], [ [coords]...] ,[ [coords],...] ]
    print("Calculating average coordinates for restrained atoms")
    restr_dict = averageCoordinates(restr_dict)
    #now the restr_dict is:
    #restr_dict[pairs]=[[req,K,D],[avgx,avgy,avgz]]

    guest_orientations = genOrientations(restr_dict, norientations=NORIENT)
    import pdb;pdb.set_trace()
    # FIXME: Create N orientations of restrained guest atoms by
    # rigib body rotations around COM
    # See Sire::Maths::rotate, Sire::Maths::matrix

    # FIXME: Only consider coordinates of host atoms to define the
    # integration domain
    space = defineIntegrationDomain(restr_dict)
    if verbose.val:
        print("Integration space")
        print(space)

    #Grid creation
    Nx = int ( round ( ( space[1][0] - space[0][0] ) / delta ) )
    Ny = int ( round ( ( space[1][0] - space[0][0] ) / delta ) )
    Nz = int ( round ( ( space[1][0] - space[0][0] ) / delta ) )
    #Constants
    print("Number of elements to be evaluated %d" %(Nx*Ny*Nz))
    print("Evaluation...")
    for i in range(0,Nx):
        for j in range(0,Ny):
            for k in range(0,Nz):
                x = space[0][0] + delta*i + delta_over_two
                y = space[0][1] + delta*j + delta_over_two
                z = space[0][2] + delta*k + delta_over_two
                # FIXME: Update all orientations centering COM on grid point
                counter= 0
                # FIXME: Compute restraint energy
                for pairs in restr_dict:
                    #restr_dict[pairs]=[[req,K,D],[coords]]
                    x_dict = float(restr_dict[pairs][1][0])
                    y_dict = float(restr_dict[pairs][1][1])
                    z_dict = float(restr_dict[pairs][1][2])

                    distance = ((x-x_dict)**2 + (y-y_dict)**2 + (z-z_dict)**2)

                    upper_bound = (restr_dict[pairs][0][0]+ restr_dict[pairs][0][2])**2
                    intmd_bound = (restr_dict[pairs][0][0])**2
                    lower_bound = max(0,(restr_dict[pairs][0][0]- restr_dict[pairs][0][2])**2 )

                    if distance <= upper_bound and distance >= intmd_bound:
                        if counter == len(restr_dict)-1:
                            U = 0.0
                        else:
                            counter+=1

                    elif distance <= intmd_bound and distance >= lower_bound:
                        if counter== len(restr_dict)-1:
                            U = 0.0
                        else:
                            counter+=1
                    else:

                        dist = (math.sqrt(x**2 + y**2 + z**2))
                        K = (restr_dict[pairs][0][1])
                        D = (restr_dict[pairs][0][2])
                        U =(K*(dist-D)**2)
                        break

                #FIXME: Add rotational volume element
                Boltz = math.exp(-beta*U)*deltavol
                Ztrans += (Boltz)
                Uavg += U*Boltz*ROT
    #Calculation of Ztot, Uavg, S, Frestraint:
    #FIXME: Remove ROT
    Ztot = Ztrans*ROT
    Uavg /= (Ztot)

    Zideal = 1661.*ROT
    Delta_F = -kbT*math.log(Ztot/Zideal)
    minTDelta_S = -T*(kb*math.log(Ztot/Zideal)+Uavg/T)


    print ("Ztrans  = %8.5f Angstrom^3" % Ztrans)
    print ("Free energy Cost of removing the restraint = %8.5f kcal/mol" % -Delta_F)

    #tidy up the folder by removing prmtop
    cmd = "rm SYSTEM.prmtop"
    os.system(cmd)
