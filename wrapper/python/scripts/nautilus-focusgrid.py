from Sire.Tools import Nautilus
from Sire.Tools import readParams

import Sire.Config

import argparse
import os
import sys

import ast

parser = argparse.ArgumentParser(description="Focusgrid, a protocol where selected coordinates are used to define "
                                             "regions with a spherical cutoff.  The thermodynamics of these regions are computed."
                                             "These coordinates are usually selected "
                                             "by using clustered densities of interest",
                                 epilog="nautilus-focusgrid is built using Sire, Numpy and mdtraj and is distributed "
                                        "under the GPL. For more information please visit "
                                        "http://siremol.org/nautilus",
                                 prog="nautilus")

parser.add_argument('-C', '--config', nargs="?",
                    help='Supply an optional Nautilus CONFIG file to control the calculation.')

parser.add_argument('--author', action="store_true",
                    help="Get information about the authors of this script.")

parser.add_argument('--version', action="store_true",
                    help="Get version information about this script.")

parser.add_argument('-g', '--gridforces', nargs="?",
                    help="Grid forces file containing averages of parameters in the grid")

parser.add_argument('-R', '--R', nargs="?",
                    help="Distance from the coordinate defining a sphere")

parser.add_argument('-c', '--coords', nargs="?",
                    help="Coordinate list defining where spheres begin")

parser.add_argument('-b', '--benchmark', action='store_true',
                    help="Benchmark the Nautilus subroutines.")
sys.stdout.write("\n")

args = parser.parse_args()

must_exit = False

if args.author:
    print("\n nautilus-focusgrids was written by George Gerogiokas and Julien Michel (C) 2014")
    print("It is based on the Nautilus Sire module.")
    must_exit = True

if args.version:
    print("nautilus-focusgrids -- from Sire release version <%s>" %Sire.__version__)
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

if args.gridforces:
    gridforces = args.gridforces
    params["gridforces"] = gridforces
elif "gridforces" in params:
    gridforces = params["gridforces"]
else:
    gridforces = "grid.forces"
    params["gridforces"] = gridforces

if args.R:
    R = float(args.R)
    params["R"] = R
else:
    R = 6
    params["R"] = R

if args.benchmark:
    params["benchmark"] = True

if not (os.path.exists(gridforces)):
    parser.print_help()
    print("\nPlease supply the name of an existing grid force file and coordinates")
    if not os.path.exists(gridforces):
        print("(cannot find grid file %s)" % gridforces)
    sys.exit(-1)

print("\nRunning nautilus-subgrids.py using file %s" % (gridforces) )
Nautilus.subgrids(params)

