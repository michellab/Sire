description = """
waterswap is a method developed and implemented using Sire that allows absolute protein-ligand binding free energies to be calculated from first-principles, condensed-phase simulations. The method is described in;

Woods, C. J., Malaisree, M., Hannongbua, S., Mulholland, A.J., “A water-swap reaction coordinate for the calculation of absolute protein-ligand binding free energies”, J. Chem. Phys. 134, 054114, 2011, DOI:10.1063/1.3519057

and

Woods, C. J., Malaisree, M., Michel, J., Long, B., McIntosh-Smith, S., Mulholland, A. J., "Rapid Decomposition and Visualisation of Protein-Ligand Binding Free Energies by Residue and by Water", Faraday Discussions 169: Molecular Simulation and Visualisation, 2014, DOI:10.1039/C3FD00125C

The method works by constructing a reaction coordinate that swaps the ligand bound to the protein with an equivalent volume of water molecules. The affect is to move the from being bound to the protein, to being free in solution, while simultaneously transferring an equivalent volume of water from being free in solution to being bound to the protein.

The input files for a waterswap calculation are the Amber format coordinate and topology files that represent the solvated protein-ligand complex (solvated using TIP3P water in a periodic, orthorhombic/cubic box).

Assuming that these files are called “complex.crd” and “complex.top”, and the ligand to be swapped with water has residue name “LIG”, then the command to run a waterswap simulation is;

sire.app/bin/waterswap -c complex.crd -t complex.top -l LIG

Sire will run the calculation using a default configuration that should be sufficient for most use cases. If you want to change any of the configuration parameters, then you can do so by writing a configuration file, the help for which can be found by running

sire.app/bin/waterswap --help-config

Once you have written a configuration file, e.g. called “CONFIG”, then you can use it via

sire.app/bin/waterswap -c complex.crd -t complex.top -C CONFIG -l LIG

Sire will automatically use all of the processor cores available on your compute node. The calculation is not fast, and the free energy averages (collected simultaneously via thermodynamic integration (TI), free energy perturbation (FEP) and Bennetts Acceptance Ratio (BAR) method) will take 1-4 days of compute time to converge.

Sire performs the calculation as a series of iterations (1000 by default), with the binding free energy (and binding free energy components) written to an output results file at the end of each iteration. This file, called output/results_????.log (where ???? is the iteration number) can be monitored during the simulation to check for convergence. At the end of the simulation, you can analyse the results using the Sire app sire.app/bin/analyse_freenrg, e.g.

sire.app/bin/analyse_freenrg -i output/freenrgs.s3 -o results.txt

This will calculate the potentials of mean force (PMFs) from the FEP, TI and BAR averages and will write them all to the file 'results.txt'. At the bottom of the results will be four estimates of the unbinding free energy (unbinding as waterswap pulls the ligand out of the protein). These four estimates are; estimate from analytic integration of TI, estimate from quadrature based integration of TI, estimate from FEP and estimate from BAR. The absolute binding free energy is the negative of the average of these four estimates, while an error can be approximated by looking at the spread of these values (e.g. by a standard deviation). If the simulation is well-converged, then these four estimates should be roughly equal.

In addition to the output/results_????.log files, Sire will also write a restart file (wsrc_restart.s3) and will write all of the free energies into freenrgs_????.s3 files (for the total free energy, and also for the residue/water components and bound and free parts). These .s3 files contain the streamed versions of the Sire objects, and can be used to restart the simulation, or to inspect the free energies or perform statistical analysis (e.g. recalculating the free energies using different integration methods, examining convergence of thermodynamic integration compared to free energy perturbation etc.). Also, Sire can be instructed to write out PDB coordinate files of the intermediates in the calculation, e.g. so you can see how the protein changes conformation as the ligand is exchanged with water. The PDB output files show only the atoms that move during the simulation, so do not worry if you only see a small cutout of your protein.

While waterswap aims to calculate the absolute binding free energy, it is not magic, and cannot overcome errors in the parameters or model of the complex. Waterswap will only be as accurate as the underlying model of the complex, and as it neglects terms such as polarisability, ionic effects and concentration effects, the results cannot be compared directly with experiment (like all absolute binding methods). Because waterswap ignores many effects, including the entropy cost of putting the ligand into the binding site, the waterswap binding free energy is an overestimate of the true value. Typically, we find that waterswap overestimates the binding free energy by 15-20 kcal mol-1 (e.g. we get values in the range -25 to -35 kcal mol-1. As such, waterswap is best used to rank binding of different ligands, or to compare binding free energies of different binding modes of the same ligand. Care should be taken when interpreting the results of waterswap, and, ideally, repeat calculations should be performed, e.g. by running on snapshots taken from an equilibrated molecular dynamics simulation.

If you need more help understanding or interpreting the results of a waterswap calculation then please feel free to get in touch via the Sire users mailing list.
"""

try:
    import sire
    sire.use_old_api()
except ImportError:
    pass

from Sire.Tools import WSRC
from Sire.Tools import readParams

import Sire.Config

import argparse
import os
import sys

parser = argparse.ArgumentParser(description="Calculate absolute binding free "
                                             "energies using waterswap",
                                 epilog="waterswap is built using Sire and is distributed "
                                        "under the GPL. For more information please visit "
                                        "http://siremol.org/waterswap.html, or type "
                                        "'waterswap --description'",
                                 prog="waterswap")

parser.add_argument('--description', action="store_true",
                    help="Print a complete description of this program.")

parser.add_argument('-H', '--help-config', action="store_true",
                    help="Get additional help regarding all of the parameters "
                         "(and their default values) that can be "
                         "set in the optionally-supplied CONFIG file")

parser.add_argument('--author', action="store_true",
                    help="Get information about the authors of this script.")

parser.add_argument('--version', action="store_true",
                    help="Get version information about this script.")

parser.add_argument('-l', '--ligand', nargs="?",
                    help="Supply the name of one of the residues in the ligand whose "
                         "binding free energy is to be calculated. By default, the ligand "
                         "will be the first non-protein, non-solvent molecule in the "
                         "input topology file.")

parser.add_argument('-t', '--topology_file', nargs="?",
                    help="The Amber topology file containing the solvated complex.")

parser.add_argument('-c', '--coordinate_file', nargs="?",
                    help="The Amber coordinate file (with periodic box) giving the coordinates "
                         "of all of the atoms in the passed topology file.")

parser.add_argument('-C', '--config', nargs="?",
                    help='Supply an optional CONFIG file to control the calculation.')

parser.add_argument('--lambda_values', type=float, nargs='+',
                    help='Lambda values for the windows used in the free energy calculation')

parser.add_argument('-n', '--num_iterations', type=int, nargs="?",
                    help='The number of waterswap iterations to perform (default 1000)')

parser.add_argument('-f', '--fast', action="store_true",
                    help="Activate the 'fast' version of waterswap used for fast, yet "
                         "potentially inaccurate calculations that are useful for equilibration.")

sys.stdout.write("\n")
args = parser.parse_args()

must_exit = False

if args.description:
    print("%s\n" % description)
    must_exit = True

if args.author:
    print("\nwaterswap was written by Christopher Woods (C) 2013-2014")
    print("It is based on the WSRC module distributed in Sire.")
    must_exit = True

if args.version:
    print("waterswap -- from Sire release version <%s>" % Sire.__version__)
    print("This particular release can be downloaded here: "
          "https://github.com/michellab/Sire/releases/tag/v%s" % Sire.__version__)
    must_exit = True

if args.help_config:
    WSRC.Parameter.printAll(True)
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
    params["protein crdfile"] = coord_file
elif "protein crdfile" in params:
    coord_file = params["protein crdfile"]
else:
    coord_file = "proteinbox.crd"
    params["protein crdfile"] = coord_file

if args.topology_file:
    top_file = args.topology_file
    params["protein topfile"] = top_file
elif "protein topfile" in params:
    top_file = params["protein topfile"]
else:
    top_file = "proteinbox.top"
    params["protein topfile"] = top_file

if not (os.path.exists(coord_file) and os.path.exists(top_file)):
    parser.print_help()
    print("\nPlease supply the name of an existing topology and coordinate file.")
    if not os.path.exists(coord_file):
        print("(cannot find coordinate file %s)" % coord_file)
    if not os.path.exists(top_file):
        print("(cannot find topology file %s)" % top_file)

    sys.exit(-1)

print("\nRunning a waterswap calculation using files %s and %s." % (top_file, coord_file))

ligand = None
if args.ligand:
    ligand = args.ligand
    params["ligand name"] = ligand
elif "ligand name" in params:
    ligand = params["ligand name"]

if ligand:
    print("The absolute binding free energy of the molecule containing residue %s "
          "will be calculated.\n" %ligand)

else:
    print("The absolute binding free energy of the first non-protein, non-solvent "
          "molecule will be calculated.\n")

lambda_values = args.lambda_values

if lambda_values:
    print("Using lambda windows %s\n" % lambda_values)
    params["lambda values"] = lambda_values

nits = args.num_iterations

if args.fast:
    print("Activating 'fast' mode. Results will be less accurate.")
    params["fast simulation"] = True
    if nits:
        nits = min(nits, 100)
    else:
        nits = 100

if nits:
    nits = int(nits)
    print("Number of iterations to perform == %d\n" % nits)
    params["nmoves"] = nits

#  Now lets run the WSRC calculation
WSRC.run(params)
