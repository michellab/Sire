description = """
analyse_freenrg_mbar is an analysis app that has been designed to analyse the output of all free energy calculations in Sire.
 analyse_freenrg_mbar reads simdata.dat files and compute free energy differences for TI and MBAR.

sire.app/bin/analyse_freenrg_mbar -i simdata.dat -o results.txt -efficiency

The option efficiency will take the statistical inefficiency of the data into account and subsample them accordingly.

If you need more help understanding or interpreting the results of an analyse_freenrg analysis then please feel free to
get in touch via the Sire users mailing list, or by creating a github issue.
"""

import argparse
import sys
import os
import numpy as np
from Sire.Tools.FreeEnergyAnalysis import SubSample
from Sire.Tools.FreeEnergyAnalysis import FreeEnergies

parser = argparse.ArgumentParser(description="Analyse free energy files to calculate "
                                             "free energies, PMFs and to view convergence.",
                                 epilog="analyse_freenrg_mbar is built using Sire and is distributed "
                                        "under GPL. For more information please visit ",
                                 prog="analyse_freenrg_mbar")

parser.add_argument('--description', action="store_true",
                    help="Print a complete description of this program.")

parser.add_argument('--author', action="store_true",
                    help="Get information about the authors of this script.")

parser.add_argument('--version', action="store_true",
                    help="Get version information about this script.")

parser.add_argument('-i', '--input', nargs='*',
                    help="Supply the name of the Sire simulation file/s containing gradients and perturbed energies to "
                         "be analysed. ")

parser.add_argument('-o', '--output', nargs=1,
                    help="""Supply the name of the file in which to write the output.""")

parser.add_argument('--efficiency', action="store_true",
                    help="Use statistical inefficiency to subsample the data")

parser.add_argument('--lam', nargs='*', type=float,
                    help="The values of lambda at which a PMF should be evaluated.")

parser.add_argument('--kT', type=float, help= 'kT in either kJ/mol or kcal/mol at which the simulation was '
                                                       'generated')

sys.stdout.write("\n")
args = parser.parse_args()

must_exit = False

if args.description:
    print("%s\n" % description)
    must_exit = True

if args.author:
    print("\nanalyse_freenrg_mbar was written by Antonia Mey (C) 2015")
    print("It is based on the pymbar analysis scripts")
    must_exit = True

if args.version:
    print("analyse_freenrg_mbar -- from Sire release version <%s>" % Sire.__version__)
    print("This particular release can be downloaded here: "
          "https://github.com/michellab/Sire/releases/tag/v%s" % Sire.__version__)
    must_exit = True

if must_exit:
    sys.exit(0)

if args.input:
    input_file = args.input
else:
    input_file = None

if args.output:
    output_file = args.output[0]
else:
    output_file = None

if args.efficiency:
    efficiency = True
else:
    efficiency = False

if args.kT:
    kT = args.kT
else:
    kT = None

if not args.lam is None:
    lam = np.array(args.lam)
else:
    lam = None

if not input_file:
    parser.print_help()
    print("\nPlease supply the name of the simulation file/s containing the reduced perturbed energies/gradients to be "
          "analysed.")
    sys.exit(-1)

if output_file:
    print("# Writing all output to file %s" % output_file)
    FILE = open(output_file, "w")
else:
    print("# Writing all output to stdout")
    FILE = sys.stdout


FILE.write("# Analysing data contained in file(s) \"%s\"\n" % input_file)

num_inputfiles = len(input_file)
if len(lam) != num_inputfiles:
    print (len(lam))
    print (num_inputfiles)
    print ("The lambda array you have provided does not match the number of simulation files provided please revise!")
    sys.exit(-1)

data = []
for i in range(num_inputfiles):
    data.append(np.loadtxt(input_file[i]))

if lam is None:
    print ("Lambda array was not giving, trying to infer lambda values from simulation files, this however is not yet "
           "implemented")
    sys.exit(-1)

if kT is None:
    print("simulation kT was not giving trying to infer kT value from simulation files, this however is not yet "
          "implemented")
    sys.exit(-1)

#Now we do the data estimation
data = np.array(data)
grad_kn = data[:,:,2]
energies_kn = data[:,:,1]
N_k = np.zeros(shape=lam.shape[0])
for k in range(0, lam.shape[0]):
    N_k[k] = data[k].shape[0]

u_kln = []
for k in range(0, len(lam)):
    u_kln.append(data[k][:,5:].transpose())
u_kln=np.array(u_kln)

subsample_obj = SubSample(grad_kn,energies_kn,u_kln,N_k)
subsample_obj.subsample_energies()
subsample_obj.subsample_gradients()

free_energy_obj = FreeEnergies(subsample_obj.u_kln, subsample_obj.N_k_energies, lam, subsample_obj.gradients_kn)
free_energy_obj.run_mbar()
free_energy_obj.run_ti()

pmf_mbar = free_energy_obj.pmf_mbar
if kT != None:
    pmf_mbar[:,1] = pmf_mbar[:,1]*kT
np.savetxt(output_file, pmf_mbar)

pmf_ti = free_energy_obj.pmf_ti


ti_out = os.path.join(os.path.dirname(output_file),'TI_'+os.path.basename(output_file))
np.savetxt(ti_out, pmf_ti, fmt=['%d', '%f'])

df_mbar_kcal = free_energy_obj.deltaF_mbar
if kT != None:
    df_mbar_kcal = df_mbar_kcal*kT

print("Free energy difference estimated with mbar is: %f kcal/mol" %df_mbar_kcal)
print("Free energy differences estimated with TI is: %f kcal/mol" %free_energy_obj.deltaF_ti)