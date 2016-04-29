description = """
analyse_freenrg is an analysis app that has been designed to analyse the output of all free energy calculations in Sire.
 analyse_freenrg reads in a Sire Saved Stream (.s3) file that contains a list of Sire.Analysis free energy objects (e.g.
 FEP, TI, Bennetts). analyse_freenrg will average and analyse these free energies according the to options you supply,
 e.g. assuming that the free energies are stored in freenrgs.s3, and you want to average iterations 100-200 from the
 simulation, and write the results to ‘results.txt’, type;

sire.app/bin/analyse_freenrg -i freenrgs.s3 -r 100 200 -o results.txt

Alternatively, if you just want to average over the last 60% of iterations, type;

sire.app/bin/analyse_freenrg -i freenrgs.s3 -o results.txt

(you can specify the percentage to average using the ‘--percent’ option)

analyse_freenrg automatically knows how many free energies are contained in the s3 file, what their types are and what
should be done to analyse the results. For example, the waterswap, ligandswap and quantomm apps all output s3 files that
contain FEP, Bennetts and TI free energy data, so analyse_freenrg knows automatically to perform FEP, Bennetts and TI
analysis on that data and to report all of the results. analyse_freenrg also knows whether or not finite difference
approximations have been used, whether forwards and backwards windows were evaluated, the temperature and conditions of
the simulation etc. The aim is that it should handle everything for you, so you can concentrate on looking at the
potential of mean force (PMF) or final result.

For FEP data, analyse_freenrg will return the FEP PMF across lambda, together with errors based on statistical
convergence (95% standard error) and the difference between forwards and backwards free energies (if available).

For Bennetts data, analyse_freenrg will return the Bennetts Acceptance Ratio PMF across lambda, with errors based on
statistical convergence (95% standard error).

For TI data, analyse_free energy will return the PMF across lambda based on polynomial fitting of the gradients and
analytic integration of the resulting function. It will also return the integral across lambda using simple quadrature.
Errors are based on statistical convergence (95% standard error) and on the difference between the forwards and
backwards finite difference gradients (if available, and if finite-difference TI was used).

If you need more help understanding or interpreting the results of an analyse_freenrg analysis then please feel free to
get in touch via the Sire users mailing list, or by creating a github issue.
"""

from Sire.Analysis import *
import Sire.Stream

import argparse
import sys
import os

parser = argparse.ArgumentParser(description="Analyse free energy files to calculate "
                                             "free energies, PMFs and to view convergence.",
                                 epilog="analyse_freenrg is built using Sire and is distributed "
                                        "under the GPL. For more information please visit "
                                        "http://siremol.org/analyse_freenrg",
                                 prog="analyse_freenrg")

parser.add_argument('--description', action="store_true",
                    help="Print a complete description of this program.")

parser.add_argument('--author', action="store_true",
                    help="Get information about the authors of this script.")

parser.add_argument('--version', action="store_true",
                    help="Get version information about this script.")

parser.add_argument('-i', '--input', nargs=1,
                    help="Supply the name of the Sire Streamed Save (.s3) file containing the "
                         "free energies to be analysed.")

parser.add_argument('-g', '--gradients', nargs='*',
                    help="Supply the name of the Sire Streamed gradients (.s3) files containing the "
                         "gradients to be analysed.")

parser.add_argument('-o', '--output', nargs=1,
                    help="""Supply the name of the file in which to write the output.""")

parser.add_argument('-r', '--range', nargs=2,
                    help="Supply the range of iterations over which to average. "
                         "By default, this will be over the last 60 percent of iterations.")

parser.add_argument('-p', '--percent', nargs=1,
                    help="Supply the percentage of iterations over which to average. By default "
                         "the average will be over the last 60 percent of iterations.")

sys.stdout.write("\n")
args = parser.parse_args()

must_exit = False

if args.description:
    print("%s\n" % description)
    must_exit = True

if args.author:
    print("\nanalyse_freenrg was written by Christopher Woods (C) 2014")
    print("It is based on the analysis tools in Sire.Analysis")
    must_exit = True

if args.version:
    print("analyse_freenrg -- from Sire release version <%s>" % Sire.__version__)
    print("This particular release can be downloaded here: "
          "https://github.com/michellab/Sire/releases/tag/v%s" % Sire.__version__)
    must_exit = True

if must_exit:
    sys.exit(0)

if args.input:
    input_file = args.input
else:
    input_file = None

if args.gradients:
    gradient_files = args.gradients
else:
    gradient_files = None

if args.output:
    output_file = args.output[0]
else:
    output_file = None

if args.range:
    range_start = int(args.range[0])
    range_end = int(args.range[1])
    if range_end < 1:
        range_end = 1
    if range_start < 1:
        range_start = 1
    if range_start > range_end:
        tmp = range_start
        range_start = range_end
        range_end = tmp
else:
    range_start = None
    range_end = None

if args.percent:
    percent = float(args.percent[0])
else:
    percent = 60.0

if not input_file:
    if gradient_files:
        input_file = gradient_files

if not input_file:
    parser.print_help()
    print("\nPlease supply the name of the .s3 file(s) containing the free energies/gradients to be analysed.")
    sys.exit(-1)

# elif not os.path.exists(input_file):
#    parser.print_help()
#    print("\nPlease supply the name of the .s3 file containing the free energies to be analysed.")
#    print("(cannot find file %s)" % input_file)
#    sys.exit(-1)

if output_file:
    print("# Writing all output to file %s" % output_file)
    FILE = open(output_file, "w")
else:
    print("# Writing all output to stdout")
    FILE = sys.stdout

#input_file = os.path.realpath(input_file)

FILE.write("# Analysing free energies contained in file(s) \"%s\"\n" % input_file)

num_inputfiles = len(input_file)

if gradient_files:
    # Multiple input files provided. Assume we have several gradients files that must be combined
    grads = {}
    fwds_grads = {}
    bwds_grads = {}

    delta_lambda = None

    for i in range(0, num_inputfiles):
        grad = Sire.Stream.load(input_file[i])

        #print(grad)
        analytic_data = grad.analyticData()
        fwds_data = grad.forwardsData()
        bwds_data = grad.backwardsData()

        #print(analytic_data)
        #print(fwds_data)
        #print(bwds_data)

        if len(analytic_data) > 0:
            # analytic gradients
            #print(analytic_data.keys())
            lamval = list(analytic_data.keys())[0]
            grads[lamval] = analytic_data[lamval]
        else:
            # finite difference gradients 
            lamval = list(fwds_data.keys())[0]
            fwds_grads[lamval] = fwds_data[lamval]
            bwds_grads[lamval] = bwds_data[lamval]
            delta_lambda = grad.deltaLambda()

    ti = None

    if len(grads) > 0:
        ti = TI(Gradients(grads))
    else:
        ti = TI(Gradients(fwds_grads, bwds_grads, delta_lambda))

    input_file = "freenrgs.s3"
    Sire.Stream.save(ti, input_file)
    #freenrgs = ti

# Only one input file provided, assumes it contains freenrgs
freenrgs = Sire.Stream.load(input_file)

results = []


def processFreeEnergies(nrgs, FILE):
    # try to merge the free enegies - this will raise an exception
    # if this object is not a free energy collection
    nrgs.merge(0, 0)

    FILE.write("# Processing object %s\n" % nrgs)

    name = nrgs.typeName().split("::")[-1]

    nits = nrgs.count()

    # get the convergence of the free energy
    convergence = {}

    for i in range(1, nits):
        try:
            convergence[i] = nrgs[i].sum().values()[-1].y()
        except:
            try:
                convergence[i] = nrgs[i].integrate().values()[-1].y()
            except:
                pass

    # now get the averaged PMF
    if range_start:
        if range_start > nits - 1:
            start = nits - 1
        else:
            start = range_start

        if range_end > nits - 1:
            end = nits - 1
        else:
            end = range_end

    else:
        end = nits - 1
        start = end - int(percent * end / 100.0)

    FILE.write("# Averaging over iterations %s to %s\n" % (start, end))

    nrg = nrgs.merge(start, end)

    try:
        pmf = nrg.sum()
    except:
        pmf = nrg.integrate()

    return (name, convergence, pmf)


try:
    results.append(processFreeEnergies(freenrgs, FILE))
except:
    for freenrg in freenrgs:
        results.append(processFreeEnergies(freenrg, FILE))

FILE.write("# Convergence\n")
FILE.write("# Iteration \n")
for result in results:
    FILE.write("# %s " % result[0])
FILE.write("\n")

i = 1
has_value = True
while has_value:
    values = []
    has_value = False
    for result in results:
        if i in result[1]:
            has_value = True
            values.append(result[1][i])
        else:
            values.append(0.0)

    if has_value:
        FILE.write("%s " % i)

        for value in values:
            FILE.write("%s " % value)

        FILE.write("\n")
        i += 1

FILE.write("# PMFs\n")

for result in results:
    FILE.write("# %s\n" % result[0])
    FILE.write("# Lambda  PMF  Maximum  Minimum \n")

    for value in result[2].values():
        FILE.write(
            "%s  %s  %s  %s\n" % (value.x(), value.y(), value.y() + value.yMaxError(), value.y() - value.yMaxError()))

FILE.write("# Free energies \n")

for result in results:
    FILE.write("# %s = %s +/- %s kcal mol-1" % (result[0], result[2].deltaG(), result[2].values()[-1].yMaxError()))

    try:
        FILE.write(" (quadrature = %s kcal mol-1)" % result[2].quadrature())
    except:
        pass

    FILE.write("#\n")

if gradient_files:
    cmd = "rm freenrgs.s3"
    os.system(cmd)
