try:
    import sire
    sire.use_old_api()
except ImportError:
    pass

from Sire.Tools import Nautilus
from Sire.Tools import readParams

import Sire.Config

import argparse
import os
import sys

parser = argparse.ArgumentParser(description="Generate cell files from a passed trajectory",
                                 epilog="nautilus-regionproperties is built using Sire, Numpy and mdtraj and is distributed "
                                        "under the GPL. For more information please visit "
                                        "http://siremol.org/nautilus",
                                 prog="nautilus")

parser.add_argument('--author', action="store_true",
                    help="Get information about the authors of this script.")

parser.add_argument('--version', action="store_true",
                    help="Get version information about this script.")

parser.add_argument('-g', '--gridforces', nargs="?",
                    help="Grid.forces file which specifies average parameters of each grid point.")

parser.add_argument('-r', '--regionfile', nargs="?",
                    help="Region file which specifies grid points to be averaged.  Text file containing only grid indices as seen in gridforces file")

parser.add_argument('-b', '--benchmark', action='store_true',
                    help="Benchmark the Nautilus subroutines.")

sys.stdout.write("\n")

args = parser.parse_args()

must_exit = False

if args.author:
    print("\n nautilus-regionproperties was written by Georgios Gerogiokas and Julien Michel (C) 2014")
    print("It is based on the Nautilus Sire module.")
    must_exit = True

if args.version:
    print("nautilus-regionproperties -- from Sire release version <%s>" %Sire.__version__)
    print("This particular release can be downloaded here: "
          "https://github.com/michellab/Sire/releases/tag/v%s" %Sire.__version__)
    must_exit = True

if must_exit:
    sys.exit(0)

# If we have been given a CONFIG file, read it now
params = {}

if args.gridforces:
    gridforces = args.gridforces
    params["gridforces"] = gridforces
elif "gridforces" in params:
    gridforces = params["gridforces"]
else:
    gridforces = "grid.forces"
    params["gridforces"] = gridforces

if args.regionfile:
    regionfile = args.regionfile
    params["regionfile"] = regionfile
elif "regionfile" in params:
    regionfile = params["regionfile"]
else:
    regionfile = "all.region"
    params["regionfile"] = regionfile

if args.benchmark:
    params["benchmark"] = True

if not (os.path.exists(gridforces) and os.path.exists(regionfile)):
    parser.print_help()
    print("\nPlease supply the name of an existing grid and region file.")
    if not os.path.exists(regionfile):
        print("cannot find region file %s" % regionfile)
    if not os.path.exists(gridforces):
        print("cannot find grid file %s" % gridforces)
    sys.exit(-1)

#print (params,args)
print("\nRunning nautilus-regionproperties.py using files %s and %s" % (gridforces, regionfile) )
Nautilus.regionproperties(params)

