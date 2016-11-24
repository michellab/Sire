#
# Evaluates the free energy cost to remove a set of distance
# restraints under standard state conditions

import os,sys, random
import math
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

    #Calculation of the mean coordinate for every atoms
    for idx in restr_dict.keys():
        values = restr_dict[idx]
        #Moltiplication by 10 since mdtraj use nm
        x_average  = (sum(val[0] for val in values)/len(values))*10
        y_average  = (sum(val[1] for val in values)/len(values))*10
        z_average  = (sum(val[2] for val in values)/len(values))*10
        #Substitution of values
        restr_dict[idx]=(x_average,y_average,z_average)
    return restr_dict

def defineIntegrationDomain(restr_dict):

    #Find the max x,y,z:
    coord = list(restr_dict.values())
    max_x = max(coordx[0] for coordx in coord)
    max_y = max(coordy[1] for coordy in coord)
    max_z = max(coordz[2] for coordz in coord)

    min_x = min(coordx[0] for coordx in coord)
    min_y = min(coordy[1] for coordy in coord)
    min_z = min(coordz[2] for coordz in coord)

    #Creation of the greatest domain plus a buffer
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
    #Check if the configuration file contains distance restraints, otherwise stop
    sim = open(simfile.val,"r")
    for line in sim.readlines():
        if "dictionary" in line:
            restraint_dictionary = line
            print("Found restraint dictionary %s" % restraint_dictionary)
        else:
            continue
    if restraint_dictionary is None:
        print("No restraint dictionary found in %s" % simfile.val)
        sys.exit(-1)

    #ast automatically transform a string into a dictionary
    rest_dict=ast.literal_eval(restraint_dictionary.strip("distance restraint dictionary="))
    #Dictionary of restrained atoms
    restrained_atoms={}
    for key in rest_dict.keys():
        for atom_idx in key:
            if atom_idx in restrained_atoms.keys():
                continue
            else:
                restrained_atoms[atom_idx] = []

    #FIXME: for the moment  K,req,D must be equal for all the restrained atoms
    req = list(rest_dict.values())[0][0]
    K=list(rest_dict.values())[0][1]
    D = list(rest_dict.values())[0][2]

    #FIXME:supposing ligand_idx < atoms_idx, remove the ligand
    ligand_idx = min(restrained_atoms.keys())
    del restrained_atoms[ligand_idx]
    print("Restrained atoms are:")
    print(list(restrained_atoms.keys()))

    start_frame = 1
    end_frame = 1000000000
    step_frame = stepframe.val
    #Check the extension of the topology file
    #.top does not work with mdtraj
    if ".top" in topfile.val:
        shutil.copy(topfile.val,"SYSTEM.prmtop")
        top_file = "SYSTEM.prmtop"
    else:
        top_file = topfile.val
    print("Loading trajectory and topology files")
    mdtraj_trajfile = mdtraj.load(trajfile.val,top=top_file)
    nframes = len(mdtraj_trajfile)
    if end_frame > (nframes - 1):
        end_frame = nframes - 1

    current_frame = start_frame

    #Aligning everything along the first frame
    print("Aligning frames along first frame of trajectory")

    aligned_traj = mdtraj_trajfile.superpose(mdtraj_trajfile,0, atom_indices=list(restrained_atoms.keys()))

    print("Processing frames")
    while (current_frame <= end_frame):

        for idx in restrained_atoms.keys():

            restrained_atoms[idx].append(aligned_traj.xyz[current_frame,idx,:].tolist())

        current_frame += step_frame


    print("Calculating average coordinates for restrained atoms")
    restrained_atoms = averageCoordinates(restrained_atoms)
    space = defineIntegrationDomain(restrained_atoms)

    if verbose.val:
        print("Integration space")
        print(space)

    #Grid creation
    Nx = int ( round ( ( space[1][0] - space[0][0] ) / delta ) )
    Ny = int ( round ( ( space[1][0] - space[0][0] ) / delta ) )
    Nz = int ( round ( ( space[1][0] - space[0][0] ) / delta ) )
    #Constants
    print("Number of elements to be evaluated %d" %(Nx*Ny*Nz))
    #Integration
    #print(restrained_atoms)
    upper_bound=(req+D)**2
    lower_bound = max(req-D,0)**2
    for i in range(0,Nx):
        for j in range(0,Ny):
            for k in range(0,Nz):
                x = space[0][0] + delta*i + delta_over_two
                y = space[0][1] + delta*j + delta_over_two
                z = space[0][2] + delta*k + delta_over_two


                for idx in range(0,len(restrained_atoms)):
                    values_x = list(restrained_atoms.values())[idx][0]
                    values_y = list(restrained_atoms.values())[idx][1]
                    values_z = list(restrained_atoms.values())[idx][2]
                    #is the point within the spheres of restraint?
                    d = (x-values_x)**2 + (y-values_y)**2 + (z-values_z)**2
                    if d <= upper_bound and d >= lower_bound:
                        # Checked all the atoms:
                        if idx == len(restrained_atoms) - 1:
                            U = 0.0
                        else:
                            idx+=1
                    else:
                        d = sqrt(x**2+y**2+z**2)
                        U = K*(d-D)**2
                        break

                Boltz = exp(-beta*U)*deltavol
                Ztrans += (Boltz)
                Uavg += U*Boltz*ROT
    #Calculation of Ztot, Uavg, S, Frestraint:
    Ztot = Ztrans*ROT
    Uavg /= (Ztot)

    Zideal = 1661.*ROT
    Delta_F = -kbT*log(Ztot/Zideal)
    minTDelta_S = -T*(kb*log(Ztot/Zideal)+Uavg/T)


    print ("Ztrans  = %8.5f Angstrom^3" % Ztrans)
    print ("Free energy Cost of removing the restraint = %8.5f kcal/mol" % -Delta_F)
