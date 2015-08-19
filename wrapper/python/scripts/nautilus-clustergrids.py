from Sire.Tools import Nautilus
from Sire.Tools import readParams

import Sire.Config

import argparse
import os
import sys

parser = argparse.ArgumentParser(description="Cluster grid points into centroids using densities from the grid and distance cutoffs for neighbour lists",
                                 epilog="nautilus-clustergrids is built using Sire, Numpy and mdtraj and is distributed "
                                        "under the GPL. For more information please visit "
                                        "http://siremol.org/nautilus",
                                 prog="nautilus")

parser.add_argument('--author', action="store_true",
                    help="Get information about the authors of this script.")

parser.add_argument('--version', action="store_true",
                    help="Get version information about this script.")

parser.add_argument('-C', '--config', nargs="?",
                    help='Supply an optional Nautilus CONFIG file to control the calculation.')

parser.add_argument('-g', '--gridforces', nargs="?",
                    help="Grid.forces file which specifies average parameters of each grid point.")

parser.add_argument('-n', '--neighcut', nargs="?",
                    help="The maximum distance (in Angstroms) between grid points to consider them neighbors, recommended value of 1.5")

parser.add_argument('-lt', '--lowt', nargs="?",
                    help="The density threshold to terminate clustering, recommended value of 1.5X greater than bulk")

parser.add_argument('-b', '--benchmark', action='store_true',
                    help="Benchmark the Nautilus subroutines.")

sys.stdout.write("\n")

args = parser.parse_args()

must_exit = False

if args.author:
    print("\n nautilus-clustergrids was written by Georgios Gerogiokas and Julien Michel (C) 2014")
    print("It is based on the Nautilus Sire module.")
    must_exit = True

if args.version:
    print("\n nautilus-clustergrids version 1.0")
    #print(Sire.Config.versionString())
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

if args.neighcut:
    neighcut = float(args.neighcut)
    params["neighcut"] = neighcut

if args.lowt:
    lowt = float(args.lowt)
    params["lowt"] = lowt

if args.benchmark:
    params["benchmark"] = True

if not (os.path.exists(gridforces)):
    parser.print_help()
    print("\nPlease supply the name of an existing grid force file.")
    if not os.path.exists(gridforces):
        print("(cannot find grid file %s)" % gridforces)
    sys.exit(-1)

print("\nRunning nautilus-clustergrids.py using file %s" % (gridforces) )
Nautilus.clustergrids(params)

