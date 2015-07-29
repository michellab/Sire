from Sire.Tools import Nautilus
from Sire.Tools import readParams

import Sire.Config

import argparse
import os
import sys

parser = argparse.ArgumentParser(description="Subgrids, subtracts grids in dx format to get differences where gridf.dx-"
                                             "gridl.dx=diff.dx.  Grids must be of identical dimensions and grid density ",
                                 epilog="nautilus-subgrids is built using Sire, Numpy and mdtraj and is distributed "
                                        "under the GPL. For more information please visit "
                                        "http://siremol.org/nautilus",
                                 prog="nautilus")

parser.add_argument('--author', action="store_true",
                    help="Get information about the authors of this script.")

parser.add_argument('--version', action="store_true",
                    help="Get version information about this script.")

parser.add_argument('-gf', '--gridf', nargs="?",
                    help="Grid dx file f to be subtracted")

parser.add_argument('-gl', '--gridl', nargs="?",
                    help="Grid dx file l used to subtract")

parser.add_argument('-d', '--diffdx', nargs="?",
                    help="Name of difference dx file")

parser.add_argument('-b', '--benchmark', action='store_true',
                    help="Benchmark the Nautilus subroutines.")
sys.stdout.write("\n")

args = parser.parse_args()

must_exit = False

if args.author:
    print("\n nautilus-subgrids was written by Georgios Gerogiokas and Julien Michel (C) 2014")
    print("It is based on the Nautilus Sire module.")
    must_exit = True

if args.version:
    print("\n nautilus-subgrids version 1.0")
    #print(Sire.Config.versionString())
    must_exit = True

if must_exit:
    sys.exit(0)

# If we have been given a CONFIG file, read it now
params = {}

if args.gridf:
    gridf = args.gridf
    params["gridf"] = gridf
elif "gridf" in params:
    gridf = params["gridf"]
else:
    gridf = "gridf.dx"
    params["gridf"] = gridf

if args.gridl:
    gridl = args.gridl
    params["gridl"] = gridl
elif "gridl" in params:
    gridl = params["gridl"]
else:
    gridl = "gridl.dx"
    params["gridl"] = gridl

if args.diffdx:
    diffdx = args.diffdx
    params["diffdx"] = diffdx
elif "diffdx" in params:
    diffdx= params["diffdx"]
else:
    diffdx = "diff.dx"
    params["diffdx"] = diffdx

if args.benchmark:
    params["benchmark"] = True

#print (params)

if not (os.path.exists(gridf) and os.path.exists(gridl)):
    parser.print_help()
    print("\nPlease supply the names of an dx files to be subtracted.")
    if not os.path.exists(os.path.exists(gridf) and os.path.exists(gridl)):
        print("(cannot find dx files %s %s)" % (gridf, gridl))
    sys.exit(-1)

print("\nRunning nautilus-subgrids.py using files %s %s" % (gridf, gridl) )
Nautilus.subgrids(params)

