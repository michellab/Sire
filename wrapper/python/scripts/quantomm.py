
from Sire.Tools import QuantumToMM
from Sire.Tools import readParams

import Sire.Config

import argparse
import os
import sys

parser = argparse.ArgumentParser(description="Calculate MM to QM/MM correction free "
                                             "energies using Metropolis-Hastings and "
                                             "the Warshel free energy cycle",
                                 epilog="quantomm is built using Sire and is distributed "
                                        "under the GPL. For more information please visit "
                                        "http://siremol.org/waterswap",
                                 prog="quantomm")

parser.add_argument('-C', '--config', nargs="?", 
                    help='Supply an optional CONFIG file to control the calculation.')

parser.add_argument('-H', '--help-config', action="store_true",
                    help="Get additional help regarding all of the parameters "
                         "(and their default values) that can be "
                         "set in the optionally-supplied CONFIG file")

parser.add_argument('-l', '--ligand', nargs="?",
                    help="Supply the name of one of the residues in the ligand whose "
                         "MM to QM/MM correction free energy is to be calculated. By default, "
                         "the ligand will be the first non-protein, non-solvent molecule in the "
                         "input topology file.")

parser.add_argument('--author', action="store_true",
                    help="Get information about the authors of this script.")

parser.add_argument('--version', action="store_true",
                    help="Get version information about this script.")

parser.add_argument('-t', '--topology_file', nargs="?",
                    help="The Amber topology file containing the solvated complex.")

parser.add_argument('-c', '--coordinate_file', nargs="?",
                    help="The Amber coordinate file (with periodic box) giving the coordinates "
                         "of all of the atoms in the passed topology file.")

parser.add_argument('--lambda_values', type=float, nargs='+',
                    help='Lambda values for the windows used in the free energy calculation')

parser.add_argument('-q', '--qm_method', nargs="?",
                    help="The quantum method to use when modelling the ligand. This should "
                         "be a string that is understood by Molpro.")

parser.add_argument('-b', '--basis_set', nargs="?",
                    help="The basis set to use when modelling the ligand. This should be "
                         "a string that is understood by Molpro.")

sys.stdout.write("\n")
args = parser.parse_args()

must_exit = False

if args.author:
    print("\nquantomm was written by Christopher Woods (C) 2013")
    print("It is based on the QuantumToMM module distributed in Sire.")
    must_exist = True

if args.version:
    print("\quantomm version 0.1")
    print(Sire.Config.versionString())
    must_exit = True

if args.help_config:
    QuantumToMM.Parameter.printAll(True)
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

if args.qm_method:
    qm_method = args.qm_method
    params["qm method"] = qm_method

if args.basis_set:
    basis_set = args.basis_set
    params["basis set"] = basis_set

if not (os.path.exists(coord_file) and os.path.exists(top_file)):
    parser.print_help()
    print("\nPlease supply the name of an existing topology and coordinate file.")
    if not os.path.exists(coord_file):
        print("(cannot find coordinate file %s)" % coord_file)
    if not os.path.exists(top_file):
        print("(cannot find topology file %s)" % top_file)

    sys.exit(-1)

print("\nRunning a quantomm calculation using files %s and %s." % (top_file,coord_file))

ligand = None
if args.ligand:
    ligand = args.ligand
    params["ligand name"] = ligand
elif "ligand name" in params:
    ligand = params["ligand name"]

if ligand:
    print("The MM to QM/MM correction free energy of the molecule containing "
          "residue %s will be calculated.\n" % (ligand))
    
else:
    print("The MM to QM/MM correction free energy of the first non-protein, non-solvent "
          "molecule will be calculated.\n")

lambda_values = args.lambda_values

if lambda_values:
    print("Using lambda windows %s\n" % lambda_values)
    params["lambda values"] = lambda_values

# Now lets run the QuantumToMM calculation
QuantumToMM.run(params)
