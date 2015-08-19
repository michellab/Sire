from Sire.Tools import Nautilus
from Sire.Tools import readParams

import Sire.Config

import argparse
import os
import sys

import ast

parser = argparse.ArgumentParser(description=   " With given center and length of an edge (in Ang) look at grid points of interest "
                                                " Defined by a cube of grid points ",
                                 epilog="nautilus-boxgrid is built using Sire, Numpy and mdtraj and is distributed "
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
                    help='Supply an optional Nautilus CONFIG file to control the calculation.')

parser.add_argument('-c', '--center', nargs="?",
                    help="")

parser.add_argument('-R', '--R', nargs="?",
                    help="Distance from the center, half edge length used to form a cubic area")


parser.add_argument('-b', '--benchmark', action='store_true',
                    help="Benchmark the Nautilus subroutines.")
sys.stdout.write("\n")

args = parser.parse_args()

must_exit = False

if args.author:
    print("\n nautilus-boxgrid was written by George Gerogiokas and Julien Michel (C) 2014")
    print("It is based on the Nautilus Sire module.")
    must_exit = True

if args.version:
    print("\n nautilus-boxgrid version 0.1")
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

if args.R:
    R = float(args.R)
    params["R"] = R
else:
    R = 10
    params["R"] = R

if args.center:
    center = ast.literal_eval(center)
    params["center"] = center
elif "center" in params:
    center = parms["center"]
else:
    print ("Using grid center specified in config file")

if args.benchmark:
    params["benchmark"] = True

if not (os.path.exists(gridforces)):
    parser.print_help()
    print("\nPlease supply the name of an existing grid force file and a center coordinate and the distance from it")
    if not os.path.exists(gridforces):
        print("(cannot find grid file %s)" % gridforces)
    sys.exit(-1)

print("\nRunning nautilus-boxgrid.py using file %s" % (gridforces) )
Nautilus.boxgrid(params)


