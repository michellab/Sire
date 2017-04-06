#Script to calculate the free energy change
# upon removing a set of host-guest distance
# restraints and applying standard state conditions
# @authors: Stefano Bosisio and Julien Michel

import os,sys, random
import math
from math import pi, cos, sin
from Sire.Tools.OpenMMMD import *
from Sire.Tools import Parameter, resolveParameters
from Sire.Maths import *
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
buffer = Parameter("buffer", 5.0, """Buffer to be added to coordinates to create domain of integration""")
delta_trans = Parameter("dtrans", 0.20, """Step size for translational volume elements, in Angstrom""")
norient = Parameter("norient", 4, """Number of orientations per [0,2pi] Euler Angles interval""")
temperature = Parameter("temperature",298 * kelvin,"""The temperature of the system""")

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
        val = restr_dict[pairs][2]
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
    max_x += buffer.val
    max_y += buffer.val
    max_z += buffer.val
    min_x -= buffer.val
    min_y -= buffer.val
    min_z -= buffer.val
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

    # First pass: work out COG of guest atoms
    guest_cog = [0.0, 0.0, 0.0]
    guest_indices = []
    for pairs in restr_dict:
        if not (pairs[0] in guest_indices):
            guest_indices.append(pairs[0])
        guest_cog[0] += restr_dict[pairs][1][0]*(1/len(restr_dict))
        guest_cog[1] += restr_dict[pairs][1][1]*(1/len(restr_dict))
        guest_cog[2] += restr_dict[pairs][1][2]*(1/len(restr_dict))
    # if only one guest atom then return null coordinates
    body = []
    if len(guest_indices) < 2:
        print ("Restraints apply to a single guest atom, isotropic case, collapsing orientations.")
        coords = []
        for pairs in restr_dict:
            coords.append( [0.0, 0.0, 0.0] )
        body.append(coords)
        return body
    # Second pass: Subtract COG to get COG centered coordinates
    body = []
    for pairs in restr_dict:
        new_x = restr_dict[pairs][1][0] - guest_cog[0]
        new_y = restr_dict[pairs][1][1] - guest_cog[1]
        new_z = restr_dict[pairs][1][2] - guest_cog[2]
        body.append( Vector(new_x, new_y, new_z) )

    # Now work out set of rotations along Euler Angles
    PI = math.pi
    TWOPI = 2*PI
    orientations = []
    for x in range(0,norientations):
        phi = (x*TWOPI)/norientations
        for y in range(0,int(norientations/2)):
            theta = (y*PI)/(norientations/2)
            for z in range(0,norientations):
                psi = (z*TWOPI)/norientations
                rot00 = cos(phi)*cos(psi)-cos(theta)*sin(phi)*sin(psi)
                rot10 = sin(phi)*cos(psi)+cos(theta)*cos(phi)*sin(psi)
                rot20 = sin(theta)*sin(psi)
                rot01 = -cos(phi)*sin(psi)-cos(theta)*sin(phi)*cos(psi)
                rot11 = -sin(phi)*sin(psi)+cos(theta)*cos(phi)*cos(psi)
                rot21 = sin(theta)*cos(phi)
                rot02 = sin(theta)*sin(phi)
                rot12 = -sin(theta)*cos(phi)
                rot22 = cos(theta)
                rotmat = Matrix(rot00,rot01,rot02,\
                                rot10,rot11,rot12,\
                                rot20,rot21,rot22)
                rotvecs = []
                #import pdb; pdb.set_trace()
                for vec in body:
                    rotvec = rotmat*vec
                    rotvecs.append(rotvec)
                orientations.append(rotvecs)
    return orientations

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
    delta_over_two = delta_trans.val/2.0
    deltavol = delta_trans.val*delta_trans.val*delta_trans.val
    ROT = 8 * pi**2
    deltarot = ROT# To adjsut after number of orientations has been determined
    kb = Sire.Units.k_boltz
    T = temperature.val.value()
    kbT = kb*T
    beta = 1/kbT

    Ztot = 0.0
    Uavg = 0

    sim_dictionary = distance_restraints_dict.val
    if sim_dictionary == {}:
        print ("Error, no distance restraints dictionary was found in the supplied config file. Abort.")
        sys.exit(-1)
    #now create a dictionary in this way:
    #dict[pairs] = {[Req,D,K] [coords] [coords]...}
    restr_dict = {}
    #create a list of host/guest indexes to be used with mdtraj for alignment
    host_indices = []
    guest_indices = []
    for pairs in sim_dictionary:
        req = sim_dictionary[pairs][0]
        K   = sim_dictionary[pairs][1]
        D   = sim_dictionary[pairs][2]
        # First entry are parameters, second and third entries for coordinates of first and second atom
        restr_dict[pairs]=[[req,K,D],[],[]]
        if not pairs[1] in host_indices:
            host_indices.append(pairs[1])
        if not pairs[0] in guest_indices:
            guest_indices.append(pairs[0])
        #idx = max(pairs)
        #hosts.append(idx)
    #FIXME: code assumes guest atoms have lower indices than host atoms
    # a more reliable algorithm could work out whether the atoms belong
    # to a guest residue?

    # FIXME:
    # Support loading a crd file only to define restraint parameters

    #load the trajectory
    start_frame = 0
    end_frame = 1000000000
    step_frame = stepframe.val
    #Check the extension of the topology file
    if ".top" in topfile.val:
        shutil.copy(topfile.val,"SYSTEM.prmtop")
        #top doesn't work with mdtraj
        top_file = "SYSTEM.prmtop"
    else:
        top_file = topfile.val
    #Check the extension of the traj file
    if ".crd" in trajfile.val:
        shutil.copy(trajfile.val,"SYSTEM.rst7")
        #crd doesn't work with mdtraj
        traj_file = "SYSTEM.rst7"
    else:
        traj_file = trajfile.val


    print("Loading trajectory and topology files")
    #
    mdtraj_trajfile = mdtraj.load(traj_file,top=top_file)
    nframes = len(mdtraj_trajfile)
    if end_frame > (nframes - 1):
        end_frame = nframes - 1
    current_frame = start_frame

    #Aligning everything along the first frame
    # Either use restrained host atoms, or a default selection if less than 3 unique atoms
    selection_default = "not water and not resname 'Na+' and not resname 'Cl-' and mass > 1"
    if len(host_indices) >2:
        selection = "index %s " % host_indices[0]
        for idx in host_indices[1:]:
            selection += " or index %s" % idx
    else:
        selection = selection_default

    print (selection)
    align_indices = mdtraj_trajfile.topology.select(selection)
    #print (align_indices)

    # FIXME: check whether alignment tolerates PBC artefacts
    print("Host: Aligning frames along first frame of trajectory")
    aligned_traj = mdtraj_trajfile.superpose(mdtraj_trajfile,0, atom_indices=align_indices)

    # First pass, collect aligned host coordinates
    while (current_frame <= end_frame):
        for pairs in restr_dict:
            idx1 = pairs[0]
            idx2 = pairs[1]
            #coord1 = aligned_traj.xyz[current_frame,idx1,:].tolist()
            coord2 = aligned_traj.xyz[current_frame,idx2,:].tolist()
            #restr_dict[pairs][1].append(coord1)
            restr_dict[pairs][2].append(coord2)
        current_frame += step_frame

    # Second pass, align against guest coordinates
    # and collect aligned guest coordinates
    if len(guest_indices) >2:
        selection = "index %s" % (guest_indices[0])
        for idx in guest_indices[1:]:
            selection += " or index %s" % idx
    else:
        selection = selection_default

    print (selection)
    align_indices = mdtraj_trajfile.topology.select(selection)
    #print (align_indices)
    # FIXME: check whether alignment tolerates PBC artefacts
    print("Guest: Aligning frames along first frame of trajectory")
    aligned_traj = mdtraj_trajfile.superpose(mdtraj_trajfile,0, atom_indices=align_indices)

    current_frame = start_frame
    while (current_frame <= end_frame):
        for pairs in restr_dict:
            idx1 = pairs[0]
            idx2 = pairs[1]
            coord1 = aligned_traj.xyz[current_frame,idx1,:].tolist()
            restr_dict[pairs][1].append(coord1)
        current_frame += step_frame

    #now restr_dict has:
    #restr_dict[lig,host]=[ [req,K,D], [ [coords]...] ,[ [coords],...] ]
    print("Calculating average coordinates for restrained atoms")
    #import pdb; pdb.set_trace()
    restr_dict = averageCoordinates(restr_dict)
    #now the restr_dict is:
    #restr_dict[pairs]=[[req,K,D],[avgx,avgy,avgz]]

    # Create N orientations of restrained guest atoms by
    # rigib body rotations around COM
    guest_orientations = genOrientations(restr_dict, norientations=norient.val)
    deltarot /= len(guest_orientations)
    # FIXME: make sure space extends well into regions where restraint energy
    # is high
    space = defineIntegrationDomain(restr_dict)
    if verbose.val:
        print("Integration space")
        print(space)

    #Grid creation
    Nx = int ( round ( ( space[1][0] - space[0][0] ) / delta_trans.val ) )
    Ny = int ( round ( ( space[1][0] - space[0][0] ) / delta_trans.val ) )
    Nz = int ( round ( ( space[1][0] - space[0][0] ) / delta_trans.val ) )
    print("Number of grid points to be evaluated %d (%d orientations per point)" %(Nx*Ny*Nz,norient.val*(norient.val/2)*norient.val))
    print("Evaluation...")
    count = 0
    free = 0
    loweight = 0
    for i in range(0,Nx):
        # FIXME: OpenMP ?
        for j in range(0,Ny):
            for k in range(0,Nz):
                count += 1
                if ( (count % 100000) == 0):
                    print ("Done %s grid points..." % (count))
                xgrid = space[0][0] + delta_trans.val*i + delta_over_two
                ygrid = space[0][1] + delta_trans.val*j + delta_over_two
                zgrid = space[0][2] + delta_trans.val*k + delta_over_two
                for orientation in guest_orientations:
                    pos = 0
                    U = 0.0
                    for pairs in restr_dict:
                        req = restr_dict[pairs][0][0]
                        k = restr_dict[pairs][0][1]
                        dtol = restr_dict[pairs][0][2]
                        host_coord = restr_dict[pairs][2]
                        guest_coord = orientation[pos]
                        # Accumulate energy
                        d2 = ((guest_coord[0]+xgrid) - host_coord[0])**2+\
                             ((guest_coord[1]+ygrid) - host_coord[1])**2+\
                             ((guest_coord[2]+zgrid) - host_coord[2])**2
                        d = math.sqrt(d2)
                        if (d > (req+dtol)):
                            U += k*(d-req-dtol)**2
                        elif (d < (req-dtol)):
                            U += k*(d-req+dtol)**2
                        else:
                            U += 0.0
                        pos += 1
                        #print ("d %s U %s " % (d,U))
                    Boltz = math.exp(-beta*U)*deltavol*deltarot
                    Uavg += U*Boltz
                    Ztot += Boltz
                    if (U < 0.000001):
                        free += 1/(norient.val*(norient.val/2)*norient.val)
                    if ( U*beta < 10):
                        loweight += 1/(norient.val*(norient.val/2)*norient.val)
                #import pdb;pdb.set_trace()
    #Calculation of Ztot, Uavg, S, Frestraint:
    free_vol = free*deltavol
    loweight_frac = loweight/float(count)
    print ("Volume where restraint is null %8.2f Angstrom^3" % (free_vol))
    print ("Fraction of points considered were restraint is under 10kbT %8.2f" % loweight_frac)
    if (loweight_frac > 0.25):
        print ("WARNING !!! The integration domain does not contain a significant number of datapoints with high restraint energy. It is possible that the domain does not cover all low restraint energy regions. Please check and if necessary increase the buffer keyword.")
        #sys.exit(-1)

    Uavg /= (Ztot)

    Zideal = 1661.*ROT
    Delta_F = -kbT*math.log(Zideal/Ztot)
    minTDelta_S = -T*(kb*math.log(Zideal/Ztot)-Uavg/T)
    print ("Ztot  = %8.2f Angstrom^3" % Ztot)
    print ("WARNING !!! Have you checked that Ztot does not increase significantly when the value of the command line argument -b/--buffer size is increased?")
    print ("WARNING !!! This calculation was done with the argument -b %s Angstrom" % buffer.val)
    print ("Free energy change upon removing the restraint and applying standard state conditions = %8.2f kcal/mol" % Delta_F)

    #tidy up the folder by removing prmtop
    cmd = "rm -f SYSTEM.prmtop SYSTEM.rst7"
    os.system(cmd)
    #import pdb; pdb.set_trace()
