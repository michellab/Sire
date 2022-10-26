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

import ast

parser = argparse.ArgumentParser(description="avggrids averages grids in dx format.  Grids must be of identical dimensions and grid density",
                                 epilog="nautilus-avggrids is built using Sire, Numpy and mdtraj and is distributed "
                                        "under the GPL. For more information please visit "
                                        "http://siremol.org/nautilus",
                                 prog="nautilus")

parser.add_argument('--author', action="store_true",
                    help="Get information about the authors of this script.")

parser.add_argument('--version', action="store_true",
                    help="Get version information about this script.")

parser.add_argument('-gs', '--avggridfiles', nargs="?",
                    help='''List of grid dx files to be averaged, "['file1','file2'..'filen']" ''')

parser.add_argument('-a', '--avgdx', nargs="?",
                    help="Name of average dx file")

parser.add_argument('-b', '--benchmark', action='store_true',
                    help="Benchmark the Nautilus subroutines.")
sys.stdout.write("\n")

args = parser.parse_args()

must_exit = False

if args.author:
    print("\n nautilus-avggrids was written by Georgios Gerogiokas and Julien Michel (C) 2014")
    print("It is based on the Nautilus Sire module.")
    must_exit = True

if args.version:
    print("nautilus-avggrids -- from Sire release version <%s>" %Sire.__version__)
    print("This particular release can be downloaded here: "
          "https://github.com/michellab/Sire/releases/tag/v%s" %Sire.__version__)
    must_exit = True

if must_exit:
    sys.exit(0)

params = {}

if args.avggridfiles:
    print (args.avggridfiles)
    avggridfiles = ast.literal_eval(str(args.avggridfiles))
    params["avggridfiles"] = avggridfiles
else:
    parser.print_help()
    print ("\nNo grid files to average specified")
    sys.exit(-1)

if args.avgdx:
    grid2 = args.avgdx
    params["avgdx"] = avgdx
elif "avgdx" in params:
    avgdx = params["avgdx"]
else:
    avgdx = "avg.dx"
    params["avgdx"] = avgdx

if args.benchmark:
    params["benchmark"] = True

if not (os.path.exists(avggridfiles[0])):
    parser.print_help()
    print("\nPlease supply the names of relevant dx files.")
    if not os.path.exists(avggridfiles[0]):
        print("(cannot find dx files %s)" % (avggridfiles[0]))
    sys.exit(-1)

print("\nRunning nautilus-avggrids.py using list of files %s" % (avggridfiles) )
Nautilus.avggrids(params)

