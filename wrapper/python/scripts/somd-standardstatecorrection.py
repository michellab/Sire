description="""
somd-standardstatecorrection is a trajectory post-processing app that computes a the free energy 
cost for removing a set of distance restraints."""

from Sire.Tools import StandardState
from Sire.Tools import readParams
from Sire.Units import *

import Sire.Config

import argparse
import os
import sys

parser = argparse.ArgumentParser(description="Evaluates the free energy cost for removing a restraint and setting standard state concentration",
                                epilog="somd-standardstatecorrection is built using Sire and is distributed "
                                        "under the GPL. For more information please visit "
                                        "http://siremol.org",
                                 prog="somd-standardstatecorrection")

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

parser.add_argument('-r', '--traj_file', nargs="?",
                    help="The trajectory file to process.")

parser.add_argument('-s', '--step', nargs="?",
                    help="The number of frames to skip between two snapshot evaluations.")

sys.stdout.write("\n")
args = parser.parse_args()

must_exit = False

if args.author:
    print("\nsomd-standardstatecorrection was written by Julien Michel (C) 2015")
    must_exit = True

if args.version:
    print("somd-standardstatecorrection -- from Sire release version <%s>" %Sire.__version__)
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

if args.topology_file:
    top_file = args.topology_file
    params["topfile"] = top_file
elif "topfile" in params:
    top_file = params["topfile"]
else:
    top_file = "system.top"
    params["topfile"] = top_file

if args.traj_file:
    traj_file = args.traj_file
    params["trajfile"] = traj_file
elif "trajfile" in params:
    traj_file = params["trajfile"]
else:
    traj_file = "traj000000001.dcd"
    params["trajfile"] = traj_file

if args.step:
    step_frame = int(args.step)
    params["step_frame"] = step_frame

if not (os.path.exists(top_file) \
        and os.path.exists(traj_file)):
    parser.print_help()
    print("\nPlease supply the name of an existing topology and trajectory file.")
    if not os.path.exists(top_file):
        print("(cannot find topology file %s)" % top_file)
    if not os.path.exists(traj_file):
        print("(cannot find traj file %s)" % traj_file)
    sys.exit(-1)

print("\nRunning a somd-standardstatcorrection calculation using files %s and %s." % (top_file, traj_file))

print (args)

# Now lets run the OpenMMMD free energy calculation
StandardState.run(params)
