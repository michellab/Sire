from Sire.Tools import Nautilus
from Sire.Tools import readParams

import Sire.Config

import argparse
import os
import sys

parser = argparse.ArgumentParser(description="Generate cell files from a passed trajectory",
                                 epilog="nautilus-traj2cell is built using Sire, Numpy and mdtraj and is distributed "
                                        "under the GPL. For more information please visit "
                                        "http://siremol.org/nautilus",
                                 prog="nautilus")

parser.add_argument('-C', '--config', nargs="?", 
                    help='Supply an optional Nautilus CONFIG file to control the calculation.')

parser.add_argument('--author', action="store_true",
                    help="Get information about the authors of this script.")

parser.add_argument('--version', action="store_true",
                    help="Get version information about this script.")

parser.add_argument('-t', '--topology_file', nargs="?",
                    help="The Amber topology file containing the system.")

parser.add_argument('-c', '--coordinate_file', nargs="?",
                    help="The Amber coordinate file giving the coordinates "
                         "of all of the atoms in the passed topology file.")

parser.add_argument('-d', '--data_file', nargs="?",
                    help="The simulation trajectory file containing coordinates and box information.")

parser.add_argument('-s', '--start_frame', nargs="?",
                    help="The frame number of the first frame to analyse.")

parser.add_argument('-e', '--end_frame', nargs="?",
                    help="The frame number of the last frame to analyse.")

parser.add_argument('-b', '--benchmark', action='store_true',
                    help="Benchmark the Nautilus subroutines.")

sys.stdout.write("\n")

args = parser.parse_args()

must_exit = False

if args.author:
    print("\n nautilus-traj2cell was written by Georgios Gerogiokas and Julien Michel (C) 2014")
    print("It is based on the Nautilus Sire module.")
    must_exit = True

if args.version:
    print("\nautilus-traj2cell version 1.0")
    #print(Sire.Config.versionString())
    must_exit = True

if must_exit:
    sys.exit(0)

# If we have been given a CONFIG file, read it now
params = {}

if args.config:
    print("Loading configuration information from file %s" % args.config)
    params = readParams(args.config)

if args.coordinate_file:
    coord_file = args.coordinate_file
    params["crdfile"] = coord_file
elif "crdfile" in params:
    coord_file = params["crdfile"]
else:
    coord_file = "system.crd"
    params["crdfile"] = coord_file

if args.topology_file:
    top_file = args.topology_file
    params["topfile"] = top_file
elif "topfile" in params:
    top_file = params["topfile"]
else:
    top_file = "system.top"
    params["topfile"] = top_file
    
if args.data_file:
    traj_file = args.data_file
    params["trajfile"] = traj_file
elif "trajfile" in params:
    traj_file = params["trajfile"]
else:
    traj_file = "traj000000001.dcd"
    params["trajfile"] = traj_file

if args.start_frame:
    start_frame = int(args.start_frame)
    params["start_frame"] = start_frame

if args.end_frame:
    end_frame = int(args.end_frame)
    params["end_frame"] = end_frame

if args.benchmark:
    params["benchmark"] = True

if not (os.path.exists(coord_file) and os.path.exists(top_file)):
    parser.print_help()
    print("\nPlease supply the name of an existing topology and coordinate file.")
    if not os.path.exists(coord_file):
        print("(cannot find coordinate file %s)" % coord_file)
    if not os.path.exists(top_file):
        print("(cannot find topology file %s)" % top_file)
    sys.exit(-1)

if not (os.path.exists(traj_file)): 
    parser.print_help()
    print("\nPlease supply the name of an existing trajectory file.")
    sys.exit(-1)

print("\nRunning nautilus-traj2cell.py using files %s, %s and %s " % (top_file, coord_file, traj_file) )
Nautilus.traj2cell(params)

