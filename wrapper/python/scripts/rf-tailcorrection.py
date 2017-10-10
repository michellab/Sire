description="""
rf-tailcorrection is a trajectory post-processing app that computes a correction 
to computed free energy changes. This app evaluates the contributions of electrostatic interactions 
from beyond the cutoff used in the original simulation according to a reaction field formalism. 
"""
from Sire.Tools import RFcutoff
from Sire.Tools import readParams
from Sire.Units import *

import Sire.Config

import argparse
import os
import sys

parser = argparse.ArgumentParser(description="Evaluates contribution of missing electrostatic interactions"
                                 "to a free energy change",
                                 epilog="rf-tailcorrection is built using Sire and is distributed "
                                        "under the GPL. For more information please visit "
                                        "http://siremol.org",
                                 prog="lj-tailcorrection")

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
parser.add_argument('-l', '--lambda_val', nargs="?",
                    help="The lambda value at which you want to run the simulation.")

parser.add_argument('-r', '--traj_file', nargs="?",
                    help="The trajectory file to process.")

parser.add_argument('-s', '--step', nargs="?",
                    help="The number of frames to skip between two snapshot evaluations.")

parser.add_argument('-i','--use_solute_inter_nrg',action="store_true",
                    help="Wheter to use the solute intermolecular energy instead of the total electrostatic potential energy.")

parser.add_argument('-f','--disable_crf',action="store_true",
                    help="Whether to disable the reaction field crf term.")

sys.stdout.write("\n")
args = parser.parse_args()

must_exit = False

if args.author:
    print("\nrf-tailcorrection was written by Stefano Bosisio and Julien Michel (C) 2017")
    print("It is based on the OpenMMMD module distributed in Sire.")
    must_exit = True

if args.version:
    print("rf-tailcorrection -- from Sire release version <%s>" %Sire.__version__)
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

if args.lambda_val:
    lambda_val = float(args.lambda_val)
    params["lambda_val"] = lambda_val

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

if args.use_solute_inter_nrg:
    params["use_solute_inter_nrg"]=True
elif "use_solute_inter_nrg" in params:
    use_solute_inter_nrg = params["use_solute_inter_nrg"]

if args.disable_crf:
    params["disable_crf"]=True
elif "disable_crf" in params:
    disable_crf = params["disable_crf"]


if not (os.path.exists(coord_file) and os.path.exists(top_file) \
        and os.path.exists(morph_file) and os.path.exists(traj_file)):
    parser.print_help()
    print("\nPlease supply the name of an existing topology, coordinate, morph and trajectory file.")
    if not os.path.exists(coord_file):
        print("(cannot find coordinate file %s)" % coord_file)
    if not os.path.exists(top_file):
        print("(cannot find topology file %s)" % top_file)
    if not os.path.exists(morph_file):
        print("(cannot find morph file %s)" % morph_file)
    if not os.path.exists(traj_file):
        print("(cannot find traj file %s)" % traj_file)
    sys.exit(-1)

print("\nRunning a rf-tailcorrection calculation using files %s, %s, %s and %s." % (top_file, coord_file, morph_file, traj_file))

print (args)

#Now lets run the LJCutoff correction code
RFcutoff.runLambda(params)
