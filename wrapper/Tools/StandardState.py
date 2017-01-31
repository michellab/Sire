#Script to calculate the Standard state correction
#It works also for different values of Req and D


import os,sys, random
import math
import numpy
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
    counter = 0
    x_avg = 0.0
    y_avg = 0.0
    z_avg = 0.0
    for pairs in restr_dict:
        coords = restr_dict[pairs][1:] #here the list of all the coords
        for val in coords:
            if counter == 0 :
                x_avg = val[0]
                y_avg = val[1]
                z_avg = val[2]
                counter+=1

            else:
                x_avg = x_avg + (val[0] - x_avg)/counter
                y_avg = y_avg + (val[1] - y_avg)/counter
                z_avg = z_avg + (val[2] - z_avg)/counter
                counter+=1

        #Substitution of values and reset each avg coords
        restr_dict[pairs][1:]=[[x_avg*10,y_avg*10,z_avg*10]]
        counter=0
        x_avg = 0.0
        y_avg = 0.0
        z_avg = 0.0

    return restr_dict

def defineIntegrationDomain(restr_dict):
    r"""Definition of the integration domain starting from host coordiantes
    Parameters
    ----------
    restr_dict : dictionary
                 restr_dict[lig_idx,host_idx] = ([req,D,K],[avgx,avgy,avgz])
                 where avgx,avgy and avgz are the average coordinates
    Returns
    ----------
    space : 3D array
            space = [(-x,-y,-z)(x,y,z)]
            where -x,-y and -z are the min values of the space and x,y and z the
            max values of the space
    """
    counter = 0
    for pairs in restr_dict:
        coords = restr_dict[pairs][1:]
        for val in coords:
            if counter==0:
                max_x = val[0]
                min_x = val[0]
                max_y = val[1]
                min_y = val[1]
                max_z = val[2]
                min_z = val[2]
                counter+=1
            else:
                if val[0]> max_x:
                    max_x = val[0]
                if val[0]<min_x:
                    min_x = val[0]
                if val[1]> max_y:
                    max_y = val[1]
                if val[1]< min_y :
                    min_y = val[1]
                if val[2]>max_z:
                    max_z = val[2]
                if val[2]<min_z :
                    min_z = val[2]

    #print(max_x,max_y,max_z,min_x,min_y,min_z)
    max_x += buff.val
    max_y += buff.val
    max_z += buff.val
    min_x -= buff.val
    min_y -= buff.val
    min_z -= buff.val
    #Adding a buffer region
    space = [(min_x,min_y, min_z), (max_x,max_y,max_z)]
    return space


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
    kb = 0.001987
    T = 298
    kbT = kb*T
    beta = 1/kbT
    ROT = 8 * pi**2
    Ztrans = 0.0
    Uavg = 0
    #import pdb; pdb.set_trace()
    #open the simfile and extract the reestraint dictionary
    #sim = open(simfile.val,"r")
    #for line in sim.readlines():
    #    if "dictionary" in line:
    #        sim_dictionary = line
    #        print("Found restraint dictionary: %s" % sim_dictionary)
    #    else:
    #        continue
    #
    #Sanity Check
    #if the dictionary was not found in the simfile exit from the script
    #if line==None:
    #    print("Error! Impossible to find dictionary restraint in sim.cfg")
    #    sys.exit(-1)
    #ast automatically transform a string into a dictionary
    #sim_dictionary=ast.literal_eval(sim_dictionary.strip("distance restraint dictionary="))

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
        restr_dict[pairs]=[[req,K,D]]
        idx = max(pairs)
        hosts.append(idx)

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
    print("Aligning frames along first frame of trajectory")
    aligned_traj = mdtraj_trajfile.superpose(mdtraj_trajfile,0, atom_indices=hosts)

    print("Processing frames")
    while (current_frame <= end_frame):
        #now for each lig:host pairs append the host coordinates
        for pairs in restr_dict:
            host_idx = max(pairs)
            restr_dict[pairs].append(aligned_traj.xyz[current_frame,host_idx,:].tolist())

        current_frame += step_frame

    #now restr_dict has:
    #restr_dict[lig,host]=[ [req,K,D], [coords],[coords],...]
    print("Calculating average coordinates for restrained atoms")
    restr_dict = averageCoordinates(restr_dict)
    #now the restr_dict is:
    #restr_dict[pairs]=[[req,K,D],[avgx,avgy,avgz]]
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

                counter= 0
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

                Boltz = math.exp(-beta*U)*deltavol
                Ztrans += (Boltz)
                Uavg += U*Boltz*ROT
    #Calculation of Ztot, Uavg, S, Frestraint:
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
