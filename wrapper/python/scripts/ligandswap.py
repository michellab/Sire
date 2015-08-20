description = """
ligandswap is a method developed and implemented using Sire that allows relative protein-ligand binding free energies to be calculated from first-principles, condensed-phase simulations. The method is currently under development and is not yet published (so use at your own risk!). The method is similar to waterswap, but instead of swapping a ligand with water, ligandswap swaps one ligand with another.

The method works by constructing a reaction coordinate that swaps one ligand (ligand 0) bound to the protein with another ligand (ligand 1) solvated in water. The affect is to unbind ligand 0 while simultaneously binding ligand 1. To ensure there is no bias towards ligand 0, a second stage calculation is performed which swaps ligand 1 bound to the protein with ligand 0 free in water. This should give the negative of the stage 1 calculation.

The input files for a ligandswap calculation are the Amber format coordinate and topology files that represent the solvated protein-ligand 0 and protein-ligand 1 complexes (solvated using TIP3P water in a periodic, orthorhombic/cubic box).

Assuming that these files are called “complex0.crd” and “complex0.top”, and that ligand 0 is called “LIG0”, and “complex1.crd” and “complex1.top”, and that ligand 1 is called “LIG1”, then the command to run a ligandswap simulation is;

sire.app/bin/ligandswap -c0 complex0.crd -t0 complex0.top -l0 LIG0 -c1 complex1.crd -t1 complex1.top -l1 LIG1

Sire will run the calculation using a default configuration that should be sufficient for most use cases. If you want to change any of the configuration parameters, then you can do so by writing a configuration file, the help for which can be found by running

sire.app/bin/ligandswap --help-config

Once you have written a configuration file, e.g. called “CONFIG”, then you can use it via

sire.app/bin/ligandswap -c0 complex0.crd -t0 complex0.top -l0 LIG0 -c1 complex1.crd -t1 complex1.top -l1 LIG1 -C config

As well as calculating relative binding free energies, ligandswap can calculate relative hydration free energies by swapping two ligands between a solvent box and vacuum. To swap to vacuum, use the ‘--vacuum’ option (see ‘--help’).

Sire will automatically use all of the processor cores available on your compute node. The calculation is not fast, and the free energy averages (collected simultaneously via thermodynamic integration (TI), free energy perturbation (FEP) and Bennetts Acceptance Ratio (BAR) method) will take 1-4 days of compute time to converge. 

Sire performs the calculation as a series of iterations (1000 by default), with the binding free energy (and binding free energy components) written to an output results file at the end of each iteration. This is performed for each stage (stage 1 being ligand 0 bound to ligand 1 bound, and stage 2 being ligand 1 bound to ligand 0 bound). These files, called output/stageX/results_????.log (where X is the stage number and ???? is the iteration number) can be monitored during the simulation to check for convergence. At the end of the simulation, you can analyse the results of each stage using the Sire app analyse_freenrg, e.g.

sire.app/bin/analyse_freenrg -i output/stage1/freenrgs.s3 -o stage1_results.txt
sire.app/bin/analyse_freenrg -i output/stage2/freenrgs.s3 -o stage2_results.txt

This will calculate the potentials of mean force (PMFs) from the FEP, TI and BAR averages and will write them all to the files 'stage1_results.txt' and ‘stage2_results.txt’. At the bottom of the results will be four estimates of the relative binding free energy. These four estimates are; estimate from analytic integration of TI, estimate from quadrature based integration of TI, estimate from FEP and estimate from BAR. The relative binding free energy is the average of the four estimates from stage 1, and the negative of the four estimates from stage 2 (stage 2 is the reverse process of stage 1). An error can be approximated by looking at the spread of these values (e.g. by a standard deviation). If the simulation is well-converged, and there is no bias between the different binding modes of ligand 0 and ligand 1, then these eight estimates should be roughly equal.

In addition to the output/results_????.log files, Sire will also write a restart file (lsrc_restart.s3) and will write all of the free energies into freenrgs_????.s3 files (for the total free energy, and also for the residue/water components and bound and free parts). These .s3 files contain the streamed versions of the Sire objects, and can be used to restart the simulation, or to inspect the free energies or perform statistical analysis (e.g. recalculating the free energies using different integration methods, examining convergence of thermodynamic integration compared to free energy perturbation etc.). Also, Sire can be instructed to write out PDB coordinate files of the intermediates in the calculation, e.g. so you can see how the protein changes conformation as the ligand is exchanged with water. The PDB output files show only the atoms that move during the simulation, so do not worry if you only see a small cutout of your protein.

While ligandswap aims to calculate the relative binding free energy, it is not magic, and cannot overcome errors in the parameters or model of the complex. ligandswap will only be as accurate as the underlying model of the complex, and as it neglects terms such as polarisability, ionic effects and concentration effects. Care should be taken when interpreting the results of ligandswap, and, ideally, repeat calculations should be performed, e.g. by running on snapshots taken from an equilibrated molecular dynamics simulation.

If you need more help understanding or interpreting the results of a ligandswap calculation then please feel free to get in touch via the Sire users mailing list.
"""

from Sire.Tools import LSRC
from Sire.Tools import readParams

import Sire.Config

import argparse
import os
import sys

parser = argparse.ArgumentParser(description="Calculate relative binding free "
                                             "energies using ligandswap",
                                 epilog="ligandswap is built using Sire and is distributed "
                                        "under the GPL. For more information please visit "
                                        "http://siremol.org/ligandswap",
                                 prog="ligandswap")

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

parser.add_argument('-l0', '--ligand0', nargs="?",
                    help="Supply the name of one of the residues in ligand 0 whose "
                         "binding free energy is to be calculated. By default, the ligand "
                         "will be the first non-protein, non-solvent molecule in the "
                         "input topology file. ligandswap calculates the relative binding "
                         "free energy of ligand0 and ligand1 (dG{ligand1} - dG{ligand0}).")

parser.add_argument('-l1', '--ligand1', nargs="?",
                    help="Supply the name of one of the residues in ligand 1 whose "
                         "binding free energy is to be calculated. By default, the ligand "
                         "will be the first non-protein, non-solvent molecule in the "
                         "input topology file. ligandswap calculates the relative binding "
                         "free energy of ligand0 and ligand1 (dG{ligand1} - dG{ligand0}).")

parser.add_argument('-t0', '--topology_file0', nargs="?",
                    help="The Amber topology file containing the solvated, equilbrated protein-ligand0 complex.")

parser.add_argument('-t1', '--topology_file1', nargs="?",
                    help="The Amber topology file containing the solvated, equilbrated protein-ligand1 complex.")

parser.add_argument('-c0', '--coordinate_file0', nargs="?",
                    help="The Amber coordinate file (with periodic box) giving the coordinates "
                         "of all of the atoms in the passed topology file of the protein-ligand0 complex.")

parser.add_argument('-c1', '--coordinate_file1', nargs="?",
                    help="The Amber coordinate file (with periodic box) giving the coordinates "
                         "of all of the atoms in the passed topology file of the protein-ligand1 complex.")

parser.add_argument('-C', '--config', nargs="?", 
                    help='Supply an optional CONFIG file to control the calculation.')

parser.add_argument('--lambda_values', type=float, nargs='+',
                    help='Lambda values for the windows used in the free energy calculation')

parser.add_argument('-n', '--num_iterations', type=int, nargs="?",
                    help='The number of waterswap iterations to perform (default 1000)')

parser.add_argument('--vacuum', action="store_true",
                    help="Swap ligands into a vacuum box rather than a water box (for relative hydration calculations).") 

sys.stdout.write("\n")
args = parser.parse_args()

must_exit = False

if args.description:
    print("%s\n" % description)
    must_exit = True

if args.author:
    print("\nligandswap was written by Christopher Woods (C) 2013-2014")
    print("It is based on the LSRC module distributed in Sire.")
    must_exit = True

if args.version:
    print("ligandswap -- from Sire release version <%s>" %Sire.__version__)
    print("This particular release can be downloaded here: "
          "https://github.com/michellab/Sire/releases/tag/v%s" %Sire.__version__)
    must_exit = True

if args.help_config:
    LSRC.Parameter.printAll(True)
    must_exit = True

if must_exit:
    sys.exit(0)

# If we have been given a CONFIG file, read it now
params = {}

if args.config:
    print("Loading configuration information from file %s" % args.config)
    params = readParams(args.config)

if args.coordinate_file0:
    coord_file0 = args.coordinate_file0
    params["crdfile0"] = coord_file0
elif "protein crdfile0" in params:
    coord_file0 = params["crdfile0"]
else:
    coord_file0 = "complex0.crd"
    params["crdfile0"] = coord_file0

if args.coordinate_file1:
    coord_file1 = args.coordinate_file1
    params["crdfile1"] = coord_file1
elif "protein crdfile1" in params:
    coord_file1 = params["crdfile1"]
else:
    coord_file1 = "complex1.crd"
    params["crdfile1"] = coord_file1

if args.topology_file0:
    top_file0 = args.topology_file0
    params["topfile0"] = top_file0
elif "topfile0" in params:
    top_file0 = params["topfile0"]
else:
    top_file0 = "complex0.top"
    params["topfile0"] = top_file0

if args.topology_file1:
    top_file1 = args.topology_file1
    params["topfile1"] = top_file1
elif "topfile1" in params:
    top_file1 = params["topfile1"]
else:
    top_file1 = "complex1.top"
    params["topfile1"] = top_file1

if args.vacuum:
    params["vacuum calculation"] = True
    print("\nPerforming a relative hydration free energy calculation.\n")

if not (os.path.exists(coord_file0) and os.path.exists(top_file0) and
        os.path.exists(coord_file1) and os.path.exists(top_file1)):
    parser.print_help()
    print("\nPlease supply the name of an existing topology and coordinate files.")
    if not os.path.exists(coord_file0):
        print("(cannot find coordinate file %s)" % coord_file0)
    if not os.path.exists(top_file0):
        print("(cannot find topology file %s)" % top_file0)
    if not os.path.exists(coord_file1):
        print("(cannot find coordinate file %s)" % coord_file1)
    if not os.path.exists(top_file1):
        print("(cannot find topology file %s)" % top_file1)

    sys.exit(-1)

ligand0 = None
if args.ligand0:
    ligand0 = args.ligand0
    params["ligand0"] = ligand0
elif "ligand0" in params:
    ligand0 = params["ligand0"]

if ligand0:
    print("Ligand 0 will be located by finding the first molecule containing residue %s" % ligand0)
    
else:
    print("Ligand 0 will be the first non-protein, non-solvent molecule found in system.")

ligand1 = None
if args.ligand1:
    ligand1 = args.ligand1
    params["ligand1"] = ligand1
elif "ligand1" in params:
    ligand1 = params["ligand1"]

if ligand1:
    print("Ligand 1 will be located by finding the first molecule containing residue %s" % ligand1)
    
else:
    print("Ligand 1 will be the first non-protein, non-solvent molecule found in system.")

print("\nRunning a ligandswap calculation calculating the difference in free energy between"
      "\nligands 0 and 1 using files %s|%s and %s|%s." % (top_file0,coord_file0,top_file1,coord_file1))

lambda_values = args.lambda_values

if lambda_values:
    print("Using lambda windows %s\n" % lambda_values)
    params["lambda values"] = lambda_values

nits = args.num_iterations

if nits:
    nits = int(nits)
    print("Number of iterations to perform == %d\n" % nits)
    params["nmoves"] = nits

# Now lets run the LSRC calculation
LSRC.run(params)
