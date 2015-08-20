from Sire.Tools import Nautilus
from Sire.Tools import readParams

import Sire.Config

import argparse
import os
import sys

parser = argparse.ArgumentParser(description="Generate grid files from cell files containing water parameters over a whole trajectory in a particular volume defined by a grid",
                                 epilog="nautilus-cell2grid is built using Sire, Numpy and mdtraj and is distributed "
                                        "under the GPL. For more information please visit "
                                        "http://siremol.org/nautilus",
                                 prog="nautilus")

parser.add_argument('-C', '--config', nargs="?", 
                    help='Supply an optional Nautilus CONFIG file to control the calculation.')

parser.add_argument('--author', action="store_true",
                    help="Get information about the authors of this script.")

parser.add_argument('--version', action="store_true",
                    help="Get version information about this script.")

parser.add_argument('-c', '--cell_dir', nargs="?",
                    help="The Amber topology file containing the system.")

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
    print("\n nautilus-cell2grid was written by Georgios Gerogiokas and Julien Michel (C) 2014")
    print("It is based on the Nautilus Sire module.")
    must_exit = True

if args.version:
    print("nautilus-cell2grid -- from Sire release version <%s>" %Sire.__version__)
    print("This particular release can be downloaded here: "
          "https://github.com/michellab/Sire/releases/tag/v%s" %Sire.__version__)
    must_exit = True

if must_exit:
    sys.exit(0)

# If we have been given a CONFIG file, read it now
params = {}

if args.config:
    print("Loading configuration information from file %s" % args.config)
    params = readParams(args.config)

if args.cell_dir:
    cell_dir = args.cell_dir
    params["cell_dir"] = cell_dir
elif "cell_dir" in params:
    cell_dir = params["cell_dir"]
else:
    cell_dir= "cell"
    params["cell_dir"] = cell_dir

if args.start_frame:
    start_frame = int(args.start_frame)
    params["start_frame"] = start_frame

if args.end_frame:
    end_frame = int(args.end_frame)
    params["end_frame"] = end_frame

if args.benchmark:
    params["benchmark"] = True

if not (os.path.exists(cell_dir)): 
    parser.print_help()
    print("\nPlease supply the cell file directory, it is missing")
    sys.exit(-1)

print("\nRunning nautilus-cell2grid.py using files %s " % (cell_dir) )
Nautilus.cell2grid(params)

