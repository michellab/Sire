try:
    import sire
    sire.use_old_api()
except ImportError:
    pass

import Sire.Config

import argparse
import os
import sys

try:
    numpy = Sire.try_import("numpy")
except:
    pass
from Sire.Tools import OpenMMMD
from Sire.Tools import readParams


parser = argparse.ArgumentParser(description="Perform molecular dynamics single topology free energy calculations "
                                             "using OpenMM",
                                 epilog="somd-freenrg is built using Sire and OpenMM and is distributed "
                                        "under the GPL. For more information please visit "
                                        "http://siremol.org",
                                 prog="somd-freenrg")

parser.add_argument('-C', '--config', nargs="?",
                    help='Supply an optional CONFIG file to control the calculation.')

parser.add_argument('-H', '--help-config', action="store_true",
                    help="Get additional help regarding all of the parameters "
                         "(and their default values) that can be "
                         "set in the optionally-supplied CONFIG file")

parser.add_argument('--author', action="store_true",
                    help="Get information about the authors of this script.")

parser.add_argument('--version', action="store_true",
                    help="Get version information about this script.")

parser.add_argument('-t', '--topology_file', nargs="?",
                    help="The Amber topology file containing the system.")

parser.add_argument('-c', '--coordinate_file', nargs="?",
                    help="The Amber coordinate file giving the coordinates "
                         "of all of the atoms in the passed topology file.")

parser.add_argument('-m', '--morph_file', nargs="?",
                    help="The morph file describing the single topology "
                         "calculation to be performed.")

parser.add_argument('-d', '--device', nargs="?",
                    help="The device ID of the GPU on which you want to run the simulation.")

parser.add_argument('-n', '--nmoves', nargs="?",
                    help="The number of Molecular Dynamics moves you want to run.")

parser.add_argument('-p', '--platform', nargs="?",
                    help="The OpenMM platform on which you want to run the simulation.")

parser.add_argument('-l', '--lambda_val', nargs="?",
                    help="The lambda value at which you want to run the simulation.")

sys.stdout.write("\n")
args = parser.parse_args()

must_exit = False

if args.author:
    print("\nsomd-freenrg was written by Gaetano Calabro, Julien Michel and Christopher Woods (C) 2014")
    print("It is based on the OpenMMMD module distributed in Sire.")
    must_exit = True

if args.version:
    print("somd-freenrg -- from Sire release version <%s>" %Sire.__version__)
    print("This particular release can be downloaded here: "
          "https://github.com/michellab/Sire/releases/tag/v%s" %Sire.__version__)
    must_exit = True

if args.help_config:
    OpenMMMD.Parameter.printAll(True)
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

if args.morph_file:
    morph_file = args.morph_file
    params["morphfile"] = morph_file
elif "morphfile" in params:
    morph_file = params["morphfile"]
else:
    morph_file = "system.morph"
    params["morphfile"] = morph_file

if args.device:
    device = args.device
    params["gpu"] = device

if args.platform:
    platform = args.platform
    params["platform"] = platform

if args.nmoves:
    nmoves = int(args.nmoves)
    params["nmoves"] = nmoves

if args.lambda_val:
    lambda_val = float(args.lambda_val)
    params["lambda_val"] = lambda_val

if not (os.path.exists(coord_file) and os.path.exists(top_file) and os.path.exists(morph_file)):
    parser.print_help()
    print("\nPlease supply the name of an existing topology, coordinate and morph file.")
    if not os.path.exists(coord_file):
        print("(cannot find coordinate file %s)" % coord_file)
    if not os.path.exists(top_file):
        print("(cannot find topology file %s)" % top_file)
    if not os.path.exists(morph_file):
        print("(cannot find morph file %s)" % morph_file)
    sys.exit(-1)

print("\nRunning a somd-freenrg calculation using files %s, %s and %s." % (top_file, coord_file, morph_file))

# Now lets run the OpenMMMD free energy calculation
OpenMMMD.runFreeNrg(params)
