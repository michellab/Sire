#
# Evaluates the free energy cost to remove a set of distance
# restraints under standard state conditions

import os,sys, random
import math
from Sire.Tools.OpenMMMD import *
from Sire.Tools import Parameter, resolveParameters

# Python dependencies
#
try:
    import mdtraj
except ImportError:
    print ("StandarState.py depends on a working install of the python module mdtraj. Please install mdtraj in your sire python.")
    sys.exit(-1)

trajfile = Parameter("trajfile", "traj000000001.dcd",
                    """File name of the trajectory to process.""")

stepframe = Parameter("step_frame",1,"""The number of frames to step to between two succcessive evaluations.""")

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
 
    # FIXME: Check that the config files contained distance restraints, otherwise stop
    print (distance_restraints_dict.val)

    #sys.exit(-1)

    # Now loop over snapshots in dcd
    start_frame = 1
    end_frame = 1000000000
    step_frame = stepframe.val

    mdtraj_trajfile = mdtraj.open(trajfile.val,'r')
    # FIXME: Load topology via mdtraj
    #mdtraj_topfile = Load topology file
    nframes = len(mdtraj_trajfile)
    if end_frame > (nframes - 1):
        end_frame = nframes - 1
    mdtraj_trajfile.seek(start_frame)
    current_frame = start_frame

    # FIXME Store coordinates of first_frame

    restrained_atoms_coordinates = []

    while (current_frame <= end_frame):
        print ("Processing frame %s " % current_frame)
        frames_xyz, cell_lengths, cell_angles = mdtraj_trajfile.read(n_frames=1)
        ## FIXME Align coordinates to first frame of the trajectory
        #aligned_frame = alignment(frames_xyz, first_frame)
        ## FIXEME Store coordinates of atoms involved in restraints using mdtraj
        ## to select atoms by restraint indices
        #restrained_atoms_coordinates_frame = selectCoordinates(aligned_frame)
        #restrained_atoms_coordinates.append(restrained_atom_coordinates_frame)
        current_frame += step_frame
        mdtraj_trajfile.seek(current_frame)

    ## FIXME Compute average coordinate of atoms involved in restraints
    #averaged_restrained_atoms_coordinates = averageCoordinates(restrained_atoms_coordinates)
    ## FIXME Define bounding cube that encompass all coordinates
    #domain = defineIntegrationDomain(average_restrained_atoms_coordinates)
    ## FIXME Perform integration of configuration integral
    #Ztrans = Integration(domain, distance_restraints_dict )
    ## FIXME Convert to free energy, assuming standard state conditions
    #DG_restraint = restraintFreeEnergy(Ztrans)
    #print ("DG_restraint = %8.5f kcal/mol")
