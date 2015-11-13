description = """
analyse_freenrg_mbar is an analysis app that has been designed to analyse the output of all free energy calculations in Sire.
 analyse_freenrg_mbar reads simdata.dat files and compute free energy differences for TI and MBAR.

sire.app/bin/analyse_freenrg_mbar -i simdata.dat -o results.txt --subsampling timeseries

The option subsampling will take the statistical inefficiency of the data into account and subsample them accordingly.

If you need more help understanding or interpreting the results of an analyse_freenrg analysis then please feel free to
get in touch via the Sire users mailing list, or by creating a github issue.
"""

import argparse
import sys
import os
import numpy as np
from Sire.Tools.FreeEnergyAnalysis import SubSample
from Sire.Tools.FreeEnergyAnalysis import FreeEnergies
from Sire.Units import *

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

parser.add_argument('--subsampling',
                    help="subsample according to either [timeseries] or [percentage] ")

parser.add_argument('--percentage', type=float,
                    help="Percentage of the data to be discarded towards equilibration.")

parser.add_argument('--lam', nargs='*', type=float,
                    help="The values of lambda at which a PMF should be evaluated.")

parser.add_argument('--temperature', type=float, help= 'temperature in [Kelvin] at which the simulation was generated.')

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

if args.subsampling:
    subsampling = args.subsampling
else:
    subsampling = 'timeseries'

if args.percentage:
    percentage = args.percentage
else:
    percentage = 100

if args.temperature:
    T = args.temperature
else:
    T = None

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
    FILE = open(output_file, "wb")

else:
    print("# Writing all output to stdout")
    FILE = sys.stdout


FILE.write(bytes("# Analysing data contained in file(s) \"%s\"\n" % input_file,  "UTF-8"))

#Todo: add header sanity checks from simfiles

#Now we do some sanity checking
num_inputfiles = len(input_file)
if len(lam) != num_inputfiles:
    print (len(lam))
    print (num_inputfiles)
    print ("The lambda array you have provided does not match the number of simulation files provided please revise!")
    sys.exit(-1)

#We will load all the data now
data = []
for i in range(num_inputfiles):
    data.append(np.loadtxt(input_file[i]))

#array of generating lambdas to check that they agree with provided lambdas
lam_sanity = []
if lam is None:
    print ("Lambda array was not given, trying to infer lambda values from simulation files, this however is not yet "
           "implemented")
    sys.exit(-1)

#sanity check array of generating temperatures to make sure they are all the same as well as match the  given
# temperature
T_sanity = []
if T is None:
    print("simulation T was not giving trying to infer T value from simulation files, this however is not yet "
          "implemented")
    sys.exit(-1)

k_boltz_J = 0.0083144621
data = np.array(data)
grad_kn = data[:,:,2] #extract the reduced gradients
energies_kn = data[:,:,1] #extract the potential energies

#N_k is the number of samples at generating thermodynamic state (lambda) k
N_k = np.zeros(shape=lam.shape[0])
for k in range(0, lam.shape[0]):
    N_k[k] = data[k].shape[0]

#Are the reduced perturbed potential energies generated at thermodynamic state k evaluated at state l, over all n
# samples. This information is contained as is in the simulation file.
u_kln = []
for k in range(0, len(lam)):
    u_kln.append(data[k][:,5:].transpose())
u_kln=np.array(u_kln)

#now we use the subsampling information to subsample the data.
subsample_obj = SubSample(grad_kn, energies_kn, u_kln, N_k, percentage=percentage, subsample=subsampling)
subsample_obj.subsample_energies()
subsample_obj.subsample_gradients()

free_energy_obj = FreeEnergies(subsample_obj.u_kln, subsample_obj.N_k_energies, lam, subsample_obj.gradients_kn)
free_energy_obj.run_mbar()
free_energy_obj.run_ti()


pmf_mbar = free_energy_obj.pmf_mbar
if T != None:
    pmf_mbar[:,1] = pmf_mbar[:,1]*T*k_boltz
np.savetxt(FILE, pmf_mbar, fmt=['%f.2', '%f'])
FILE.close()

pmf_ti = free_energy_obj.pmf_ti
if T != None:
    pmf_ti[:,1] = pmf_ti[:,1]*T*k_boltz_J


ti_out = os.path.join(os.path.dirname(output_file),'TI_'+os.path.basename(output_file))
np.savetxt(ti_out, pmf_ti, fmt=['%f.2', '%f'])

#writing out free energy differences and errors
#TODO: rewrite this with a decorator!
df_mbar_kcal = free_energy_obj.deltaF_mbar
df_mbar_kJ = free_energy_obj.deltaF_mbar
df_ti_kcal = free_energy_obj.deltaF_ti
dDf_mbar_kcal = free_energy_obj.errorF_mbar
dDf_mbar_kJ = free_energy_obj.errorF_mbar
if T != None:
    df_mbar_kcal = df_mbar_kcal*T*k_boltz
    df_mbar_kJ = df_mbar_kcal*T*k_boltz_J
    #df_ti_kcal = df_ti_kcal*T*k_boltz
    dDf_mbar_kcal = dDf_mbar_kcal*T*k_boltz
    dDf_mbar_kJ = dDf_mbar_kJ*T*k_boltz_J

print("Free energy difference estimated with MBAR is: %f kcal/mol +/- %f kcal/mol" %(df_mbar_kcal,dDf_mbar_kcal))
print("Free energy difference estimated with MBAR is: %f kJ/mol +/- %f kJ/mol" %(df_mbar_kJ,dDf_mbar_kJ))
print("Free energy difference estimated with TI is: %f kcal/mol" %df_ti_kcal)