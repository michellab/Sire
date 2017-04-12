# coding=utf-8
import bz2
import glob
import os
import sys
import warnings

import Sire.Stream
import argparse
from Sire.Analysis import *
from Sire.Tools.FreeEnergyAnalysis import FreeEnergies
from Sire.Tools.FreeEnergyAnalysis import SubSample
from Sire.Units import *

# Import the asciiplot (ap) plotting library - this needs numpy
try:
    numpy = Sire.try_import("numpy")
    from Sire.Tools import ap
except:
    pass

description_mbar = """
analyse_freenrg mbar is an analysis app that has been designed to analyse the output of all free energy calculations in
Sire. analyse_freenrg_mbar reads simdata.dat files and compute free energy differences for TI and MBAR.
sire.app/bin/analyse_freenrg mbar -i simdata.dat -o results.txt --subsampling timeseries
The option subsampling will take the statistical inefficiency of the data into account and subsample them accordingly.
If you need more help understanding or interpreting the results of an analyse_freenrg analysis then please feel free to
get in touch via the Sire users mailing list, or by creating a github issue.
"""

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


def processFreeEnergies(nrgs, FILE):
    # try to merge the free enegies - this will raise an exception
    # if this object is not a free energy collection
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

    return name, convergence, pmf


if __name__ == '__main__':

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

    subparser = parser.add_subparsers(description="Analyse free energy files to calculate "
                                                  "free energies, PMFs and to view convergence.",
                                      prog="analyse_freenrg")
    parser_mbar = subparser.add_parser('mbar')

    parser_mbar.add_argument('--description', action="store_true",
                             help="Print a complete description of this program.")

    parser_mbar.add_argument('--author', action="store_true",
                             help="Get information about the authors of this script.")

    parser_mbar.add_argument('--version', action="store_true",
                             help="Get version information about this script.")

    parser_mbar.add_argument('-i', '--input', nargs='*',
                             help="Supply the name of the Sire simulation.dat files containing "
                                  "gradients and perturbed energies to be analysed. Valid options are: "
                                  "'simfile1.dat, simfile2.dat', or 'simfile1.dat simfile2.dat', or "
                                  "wildcard arguments 'lambda*/simfile.dat'")

    parser_mbar.add_argument('-o', '--output',
                             help="""Supply the name of the file in which to write the output.""")

    parser_mbar.add_argument('-p', '--percent', type=float,
                             help="Supply the percentage of iterations over which to average. By default "
                                  "the average will be over the last 60 percent of iterations.")
    parser_mbar.add_argument('--lam', type=float, nargs='*',
                             help="lambda array")

    parser_mbar.add_argument('--temperature', type=float,
                             help="lambda array")

    parser_mbar.add_argument('--subsampling', action="store_true",
                             help="Subsampling of data according to statistical inefficiency should be done.")

    parser_mbar.add_argument('--overlap', action="store_true",
                             help="Compute the overlap matrix")

    args, unknown = parser.parse_known_args()
    if len(sys.argv) < 2:
        parser.print_help()
        print("\nPlease supply the name of the .s3 file(s) containing the free energies/gradients to be analysed.")
        sys.exit(1)

    ##########################################
    #         Free energy with s3 files      #
    ##########################################
    if sys.argv[1] != 'mbar':
        for u in unknown:
            mbar_list = ['--temperature', '--subsampling', '--lam', '--overlap']
            if u in mbar_list:
                parser_mbar.print_help()
                print ('analyse_freenrg: error: mbar command argument: %s ' % u)
                print ('Try adding the mabr command option')
                sys.exit(1)
            else:
                parser_mbar.print_help()
                print ('analyse_freenrg: error: unrecognized arguments: %s ' % u)
                sys.exit(1)

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
            input_file = args.input[0]
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
                print (input_file)

        if not input_file:
            parser.print_help()
            print("\nPlease supply the name of the .s3 file(s) containing the free energies/gradients to be analysed.")
            sys.exit(-1)
        if not input_file[0].endswith('.s3'):
            parser.print_help()
            print(
                "\nInvalid file format. Please supply the name of the .s3 file(s) containing the free energies/gradients "
                "to be analysed.")
            sys.exit(-1)

        # elif not os.path.exists(input_file):
        #    parser.print_help()
        #    print("\nPlease supply the name of the .s3 file containing the free energies to be analysed.")
        #    print("(cannot find file %s)" % input_file)
        #    sys.exit(-1)

        if output_file:
            print("# Writing all output to file %s" % output_file)
            FILE = open(output_file, "wb")
        else:
            print("# Writing all output to stdout")
            FILE = sys.stdout

        # input_file = os.path.realpath(input_file)

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

                # print(grad)
                analytic_data = grad.analyticData()
                fwds_data = grad.forwardsData()
                bwds_data = grad.backwardsData()

                # print(analytic_data)
                # print(fwds_data)
                # print(bwds_data)

                if len(analytic_data) > 0:
                    # analytic gradients
                    # print(analytic_data.keys())
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
            # freenrgs = ti

        # Only one input file provided, assumes it contains freenrgs
        freenrgs = Sire.Stream.load(input_file)

        results = []

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

        x = []
        y = []

        i = 1
        has_value = True
        while has_value:
            values = []
            has_value = False
            for result in results:
                if i in result[1]:
                    has_value = True
                    values.append(result[1][i])
                    x.append(i)
                    y.append(result[1][i])
                else:
                    values.append(0.0)

            if has_value:
                FILE.write("%s " % i)

                for value in values:
                    FILE.write("%s " % value)

                FILE.write("\n")
                i += 1

        # now plot the graph of the convergence (just the last 90% of iterations)
        try:
            p = ap.AFigure()
            strt = int(len(x) / 10)
            conv_plot = p.plot(numpy.array(x[strt:]), numpy.array(y[strt:]), ".")
            print("\nPlot of free energy versus iteration")
            print(conv_plot + "\n")
        except Exception as e:
            print("Error plotting graph: %s" % e)

        FILE.write("# PMFs\n")

        for result in results:
            x = []
            y = []
            FILE.write("# %s\n" % result[0])
            FILE.write("# Lambda  PMF  Maximum  Minimum \n")

            for value in result[2].values():
                FILE.write(
                    "%s  %s  %s  %s\n" % (
                        value.x(), value.y(), value.y() + value.yMaxError(), value.y() - value.yMaxError()))
                x.append(value.x())
                y.append(value.y())

            try:
                p = ap.AFigure()
                pmf_plot = p.plot(numpy.array(x), numpy.array(y), '_o')
                print("\nPMF Plot of free energy versus lambda")
                print(pmf_plot + "\n")
            except:
                pass

        FILE.write("# Free energies \n")

        for result in results:
            FILE.write(
                "# %s = %s +/- %s kcal mol-1" % (result[0], result[2].deltaG(), result[2].values()[-1].yMaxError()))

            try:
                FILE.write(" (quadrature = %s kcal mol-1)" % result[2].quadrature())
            except:
                pass

            FILE.write("#\n")

        if gradient_files:
            cmd = "rm freenrgs.s3"
            os.system(cmd)


    else:
        must_exit = False
        # MBAR is true and we run and MBAR analysis
        print('Simulation data is analysed using the python module pymbar')
        print('----------------------------------------------------------')

        if args.author:
            print("\nanalyse_freenrg mbar was written by Antonia Mey (C) 2015-2017")
            print("It is based on the pymbar analysis scripts (https://github.com/choderalab/pymbar)")
            must_exit = True
        if args.version:
            print("analyse_freenrg mbar -- from Sire release version <%s>" % Sire.__version__)
            print("This particular release can be downloaded here: "
                  "https://github.com/michellab/Sire/releases/tag/v%s" % Sire.__version__)
            must_exit = True
        if args.description:
            print("%s\n" % description_mbar)
            must_exit = True
        if must_exit:
            sys.exit(0)

        input_files = args.input
        output_file = args.output
        percentage = args.percent
        T = args.temperature

        # boolean arguemtns
        test_overlap = args.overlap
        subsampling = args.subsampling

        if not input_files:
            parser_mbar.print_help()
            print(
                "\nPlease supply the name of the simulation file/s containing the reduced perturbed energies/gradients "
                "to be analysed.")
            sys.exit(-1)
        if input_files[0].endswith('.s3'):
            parser_mbar.print_help()
            print(
                "\nInvalid file format for MBAR analysis. Please supply the names of the simulation.dat file/s "
                "containing the reduced perturbed energies/gradients to be analysed.")
            sys.exit(-1)
        elif '*' in input_files[0]:
            input_files = glob.glob(input_files[0])
            input_files.sort()

        if output_file:
            print("# Writing all output to file %s" % output_file)
            FILE = open(output_file, "wb")
        else:
            print("# Writing all output to stdout")
            FILE = sys.stdout.buffer

        FILE.write(bytes("# Analysing data contained in file(s) %s\n" % input_files, "UTF-8"))

        # sanity checking of other input parameters:
        # processing lambda values
        lamvals = None
        if not args.lam:
            print ("#Lambda array was not given, trying to infer lambda values from simulation files...")
        else:
            lamvals = args.lam
        lam = None
        for f in input_files:
            print ('working on input file %s' % f)
            # Compressed file
            if f.endswith('.bz2'):
                bz_file = bz2.BZ2File(f)
                for line in bz_file:
                    if line.startswith(b'#Alchemical'):
                        if lamvals is None:
                            lamvals = (line.split(b'(')[-1].split(b')')[0].split(b','))
                            lamvals = numpy.array([float(i) for i in lamvals])
                            lam = lamvals

                        else:
                            lam = (line.split(b'(')[-1].split(b')')[0].split(b','))
                            lam = numpy.array([float(i) for i in lam])
                        break
                        # Normal file
            else:
                for line in open(f):
                    if line.startswith('#Alchemical'):
                        if lamvals is None:
                            lamvals = (line.split('(')[-1].split(')')[0].split(','))
                            lamvals = numpy.array([float(i) for i in lamvals])
                            lam = lamvals

                        else:
                            lam = (line.split('(')[-1].split(')')[0].split(','))
                            lam = numpy.array([float(i) for i in lam])
                        break
            if not numpy.array_equal(lam, lamvals):
                print ("Lambda arrays do not match! Make sure your input data is consistent")
                print (lam)
                print (lamvals)
                sys.exit(-1)

        # We will load all the data now
        data = []
        for f in input_files:
            data.append(numpy.loadtxt(f))

        # N_k is the number of samples at generating thermodynamic state (lambda) k
        N_k = numpy.zeros(shape=lamvals.shape[0], dtype='int32')
        for k in range(0, lamvals.shape[0]):
            N_k[k] = data[k].shape[0]

        max_sample = int(max(N_k))
        grad_kn = numpy.zeros(shape=(lamvals.shape[0], max_sample))
        energies_kn = numpy.zeros(shape=(lamvals.shape[0], max_sample))

        for k in range(0, N_k.shape[0]):
            grad_kn[k, 0:N_k[k]] = data[k][:, 2]  # get the gradient information
            energies_kn[k, 0:N_k[k]] = data[k][:, 1]  # get the potential energies from the file.

        # Are the reduced perturbed potential energies generated at thermodynamic state k evaluated at state l, over
        # all n samples. This information is contained as is in the simulation file.
        u_kln = numpy.zeros(shape=(lamvals.shape[0], lamvals.shape[0], max_sample))
        for k in range(0, lamvals.shape[0]):
            u_kln[k, :, 0:N_k[k]] = data[k][:, 5:].transpose()

        # now we use the subsampling information to subsample the data.
        subsample_obj = SubSample(grad_kn, energies_kn, u_kln, N_k, percentage=percentage, subsample=subsampling)
        subsample_obj.subsample_energies()
        subsample_obj.subsample_gradients()

        free_energy_obj = FreeEnergies(subsample_obj.u_kln, subsample_obj.N_k_energies, lamvals,
                                       subsample_obj.gradients_kn)
        print ('#running mbar ====================================================')
        free_energy_obj.run_mbar(test_overlap)
        print ('#running mbar done ===============================================')
        free_energy_obj.run_ti()

        ti_warn_msg = ''
        mbar_warn_msg = ''
        if subsample_obj.gradients_kn.shape[1]<50:
            ti_warn_msg = ' #WARNING SUBSAMLING GRADIENTS RESULTED IN LESS THAN 50 SAMPLES, CONSIDER RERUN WITHOUT SUBSAMPLE OPTION'

        if subsample_obj.u_kln.shape[2]<50:
            mbar_warn_msg = ' #WARNING SUBSAMLING ENERGIES RESULTED IN LESS THAN 50 SAMPLES, CONSIDER RERUN WITHOUT SUBSAMPLE OPTION'

        if test_overlap:
            M = free_energy_obj.overlap_matrix
            diag_elements = numpy.array([numpy.diag(M, k=1), numpy.diag(M, k=-1)])
            if numpy.min(diag_elements) < 0.03:
                warnings.warn(
                    'Off diagonal elements of the overlap matrix are smaller than 0.03! Your free energy estiamte is '
                    'not reliable!')
            FILE.write(bytes('#Overlap matrix\n', "UTF-8"))
            numpy.savetxt(FILE, M, fmt='%.4f')

        # mbar DG for neighbourting lambda
        pairwise_F = free_energy_obj.pairwise_F
        if T is not None:
            pairwise_F[:, 2] = pairwise_F[:, 2] * T * k_boltz
            pairwise_F[:, 3] = pairwise_F[:, 3] * T * k_boltz
            FILE.write(bytes('#PMF from MBAR in kcal/mol\n', "UTF-8"))
        else:
            FILE.write(bytes('#PMF from MBAR in reduced units\n', "UTF-8"))
        numpy.savetxt(FILE, pairwise_F, fmt=['%.4f', '%.4f', '%.4f', '%.4f'])

        # mbar pmf
        pmf_mbar = free_energy_obj.pmf_mbar
        if T is not None:
            pmf_mbar[:, 1] = pmf_mbar[:, 1] * T * k_boltz
            pmf_mbar[:, 2] = pmf_mbar[:, 2] * T * k_boltz
            FILE.write(bytes('#PMF from MBAR in kcal/mol\n', "UTF-8"))
        else:
            FILE.write(bytes('#PMF from MBAR in reduced units\n', "UTF-8"))
        numpy.savetxt(FILE, pmf_mbar, fmt=['%.4f', '%.4f', '%.4f'])

        # ti mean gradients and std.
        grad_mean = numpy.mean(subsample_obj.gradients_kn, axis=1)
        grad_std = numpy.std(subsample_obj.gradients_kn, axis=1)
        if T is not None:
            grad_mean = grad_mean * T * k_boltz
            grad_std = grad_std * T * k_boltz
            FILE.write(bytes('#TI average gradients and standard deviation in kcal/mol\n', "UTF-8"))
        else:
            FILE.write(bytes('#TI average gradients and standard deviation in reduced units\n', "UTF-8"))
        grad_data = numpy.column_stack((lamvals, grad_mean, grad_std))
        numpy.savetxt(FILE, grad_data, fmt=['%.4f', '%.4f', '%.4f'])

        # TI pmf
        pmf_ti = free_energy_obj.pmf_ti
        if T is not None:
            pmf_ti[:, 1] = pmf_ti[:, 1] * T * k_boltz
            FILE.write(bytes('#PMF from TI in kcal/mol\n', "UTF-8"))
        else:
            FILE.write(bytes('#PMF from TI in reduced units\n', "UTF-8"))
        numpy.savetxt(FILE, pmf_ti, fmt=['%.4f', '%.4f'])

        df_mbar_kcal = free_energy_obj.deltaF_mbar
        df_mbar_kJ = free_energy_obj.deltaF_mbar
        df_ti_kcal = free_energy_obj.deltaF_ti
        dDf_mbar_kcal = free_energy_obj.errorF_mbar
        dDf_mbar_kJ = free_energy_obj.errorF_mbar
        if T is not None:
            df_mbar_kcal = df_mbar_kcal * T * k_boltz
            dDf_mbar_kcal = dDf_mbar_kcal * T * k_boltz
            df_ti_kcal = free_energy_obj.deltaF_ti * T * k_boltz
            FILE.write(
                bytes("#MBAR free energy difference in kcal/mol: \n%f, %f %s\n" % (df_mbar_kcal, dDf_mbar_kcal, mbar_warn_msg), "UTF-8"))
            FILE.write(bytes("#TI free energy difference in kcal/mol: \n%f %s \n" % (df_ti_kcal, ti_warn_msg), "UTF-8"))

        else:
            FILE.write(
                bytes("#MBAR free energy difference in reduced units: \n%f, %f %s\n" % (df_mbar_kcal, dDf_mbar_kcal, mbar_warn_msg),
                      "UTF-8"))
            FILE.write(bytes("#TI free energy difference in reduced units: \n%f %s\n" % (df_ti_kcal, ti_warn_msg), "UTF-8"))

        FILE.close()
