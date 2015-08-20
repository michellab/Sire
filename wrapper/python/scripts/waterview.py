description = """
waterview is a simple app that can be used to visualise how water binds to a protein. The input for the app is a snapshot from a dynamics simulation of a solvated protein-ligand complex (solvated by TIP3P in a cubic periodic box). waterview deletes the ligand from the complex, replaces it with an equal volume of water, and then runs Monte Carlo on all of the waters (the ones in the binding site and the ones swapped in in-place of the ligand). The output from the app is a series of PDB files showing the dynamics of the water, together with a grid DX file that contains the average water occupancy at each point on a 3D grid in the binding site.

Assuming that your protein ligand complex is contained in the Amber topology and coordinate files “system.top” and “system.crd”, and the ligand contains a residue called “LIG”, then the command;

sire.app/bin/waterview -c system.crd -t system.top -l LIG

will run the calculation. The calculation is performed using two stages; stage 1 looks at the water molecules when the ligand is present in the binding site, while stage 2 looks at the water molecules when the ligand is substituted by water. The simulation will write the output to PDB files called “stage1_mobile_XXXX.pdb” and “stage2_mobile_XXXX.pdb” (where XXXX is an iteration number). The average grid occupancy files are “stage1_vol.dx” and “stage2_vol.dx”.
waterview is a relatively quick calculation, taking minutes to hours to complete (depending on your machine and the number of iterations you want to perform). There are some more advanced options that can be set in a config file (see “--help-config”). As with the other apps, you can specify a config file using the “-C” option.

If you need more help understanding or interpreting the results of a ligandswap calculation then please feel free to get in touch via the Sire users mailing list.
"""

from Sire.Tools import WaterView
from Sire.Tools import readParams

import Sire.Config

import argparse
import os
import sys

parser = argparse.ArgumentParser(description="Visualise the distribution of water molecules in "
                                             "a binding site, both with and without a ligand present.",
                                 epilog="waterview is built using Sire and is distributed "
                                        "under the GPL. For more information please visit "
                                        "http://siremol.org/waterview",
                                 prog="waterview")

parser.add_argument('--description', action="store_true",
                    help="Print a complete description of this program.")

parser.add_argument('--author', action="store_true",
                    help="Get information about the authors of this script.")

parser.add_argument('--version', action="store_true",
                    help="Get version information about this script.")

parser.add_argument('-H', '--help-config', action="store_true",
                    help="Get additional help regarding all of the parameters "
                         "(and their default values) that can be "
                         "set in the optionally-supplied CONFIG file")

parser.add_argument('-l', '--ligand', nargs="?",
                    help="Supply the name of one of the residues in the ligand around "
                         "which the water distribution will be viewed. By default, the ligand "
                         "will be the first non-protein, non-solvent molecule in the "
                         "input topology file.")

parser.add_argument('-t', '--topology_file', nargs="?",
                    help="The Amber topology file containing the solvated complex.")

parser.add_argument('-c', '--coordinate_file', nargs="?",
                    help="The Amber coordinate file (with periodic box) giving the coordinates "
                         "of all of the atoms in the passed topology file.")

parser.add_argument('-C', '--config', nargs="?", 
                    help='Supply an optional CONFIG file to control the calculation.')

sys.stdout.write("\n")
args = parser.parse_args()

must_exit = False

if args.description:
    print("%s\n" % description)
    must_exit = True

if args.author:
    print("\nwaterview was written by Christopher Woods (C) 2013")
    print("It is based on the WaterView module distributed in Sire.")
    must_exit = True

if args.version:
    print("waterview -- from Sire release version <%s>" %Sire.__version__)
    print("This particular release can be downloaded here: "
          "https://github.com/michellab/Sire/releases/tag/v%s" %Sire.__version__)
    must_exit = True

if args.help_config:
    WaterView.Parameter.printAll(True)
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
elif "protein crdfile" in params:
    coord_file = params["crdfile"]
else:
    coord_file = "complex.crd"
    params["crdfile"] = coord_file

if args.topology_file:
    top_file = args.topology_file
    params["topfile"] = top_file
elif "topfile" in params:
    top_file = params["topfile"]
else:
    top_file = "complex.top"
    params["topfile"] = top_file

if not (os.path.exists(coord_file) and os.path.exists(top_file)):
    parser.print_help()
    print("\nPlease supply the name of an existing topology and coordinate file.")
    if not os.path.exists(coord_file):
        print("(cannot find coordinate file %s)" % coord_file)
    if not os.path.exists(top_file):
        print("(cannot find topology file %s)" % top_file)

    sys.exit(-1)

print("\nRunning a waterview calculation using files %s and %s." % (top_file,coord_file))

ligand = None
if args.ligand:
    ligand = args.ligand
    params["ligand name"] = ligand
elif "ligand name" in params:
    ligand = params["ligand name"]

if ligand:
    print("The water molecules around the ligand containing residue %s "
          "will be visualised.\n" % (ligand))
    
else:
    print("The water molecules around the first non-protein, non-solvent "
          "molecule will be visualised.\n")

# Now lets run the WaterView calculation
WaterView.run(params)
