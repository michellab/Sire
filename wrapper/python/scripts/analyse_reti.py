description = """
analyse_reti is an analysis app that has been designed to analyse the replica exchange trajectory of all RETI
simulations in Sire. analyse_reti reads in a Sire Saved Stream (.s3) file that contains a Sire simulation restart file
(typically called ???_restart.s3). analyse_reti will extract the replica exchange trajectory and print it out in a
format that will allow easy graphing. If you want to analyse the replica exchange moves from iterations 100 to 200 from
the files restart.s3 and to write the results out to 'results.txt' then type;

sire.app/bin/analyse_reti -i restart.s3 -r 100 200 -o results.txt

Alternatively, if you just want to analyse the last 60% of iterations, type;

sire.app/bin/analyse_reti -i restart.s3 --percent 60 -o results.txt

You can also analyse all replica exchange iterations using

sire.app/bin/analyse_reti -i restart.s3 -o results.txt

If you don't supply the name of the output file, then the analysis is printed to the screen.

If you need more help understanding or interpreting the results of an analyse_reti analysis then please feel free to get
in touch via the Sire users mailing list, or by creating a github issue.
"""

import Sire.Stream

import argparse
import sys
import os

parser = argparse.ArgumentParser(description="Analyse Sire restart files to get information "
                                             "about the replica exchange moves performed during the simulation.",
                                 epilog="analyse_reti is built using Sire and is distributed "
                                        "under the GPL. For more information please visit "
                                        "http://siremol.org/analyse_reti",
                                 prog="analyse_reti")

parser.add_argument('--description', action="store_true",
                    help="Print a complete description of this program.")

parser.add_argument('--author', action="store_true",
                    help="Get information about the authors of this script.")

parser.add_argument('--version', action="store_true",
                    help="Get version information about this script.")

parser.add_argument('-i', '--input', nargs=1,
                    help="Supply the name of the Sire Streamed Save (.s3) file containing the "
                         "free energies to be analysed.")

parser.add_argument('-o', '--output', nargs=1,
                    help="""Supply the name of the file in which to write the output.""")

parser.add_argument('-r', '--range', nargs=2,
                    help="Supply the range of iterations over which to analyse. "
                         "By default, this will be over the last 60 percent of iterations.")

parser.add_argument('-p', '--percent', nargs=1,
                    help="Supply the percentage of iterations over which to analyse. By default "
                         "the analysis will be over the last 100 percent of iterations.")

sys.stdout.write("\n")
args = parser.parse_args()

must_exit = False

if args.description:
    print("%s\n" % description)
    must_exit = True

if args.author:
    print("\nanalyse_reti was written by Christopher Woods (C) 2014")
    must_exit = True

if args.version:
    print("analyse_reti -- from Sire release version <%s>" %Sire.__version__)
    print("This particular release can be downloaded here: "
          "https://github.com/michellab/Sire/releases/tag/v%s" %Sire.__version__)
    must_exit = True

if must_exit:
    sys.exit(0)

if args.input:
    input_file = args.input[0]
else:
    input_file = None

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
    percent = 100.0

if not input_file:
    parser.print_help()
    print("\nPlease supply the name of the .s3 file containing the Sire restart file.")
    sys.exit(-1)

elif not os.path.exists(input_file):
    parser.print_help()
    print("\nPlease supply the name of the .s3 file containing the Sire restart file.")
    print("(cannot find file %s)" % input_file)
    sys.exit(-1)

if output_file:
    print("Writing all output to file %s\n" % output_file)
    FILE = open(output_file, "w")
else:
    print("Writing all output to stdout\n")
    FILE = sys.stdout

input_file = os.path.realpath(input_file)

FILE.write("Analysing the replica exchange moves contained in file \"%s\"\n" % input_file)

(system, moves) = Sire.Stream.load(input_file)  

def getRepExMoves(moves):
    if moves.what() == "SireMove::RepExMove":
        return moves
    else:
        try:
            for i in range(0, moves.nSubMoveTypes()):
                if moves[i].what() == "SireMove::RepExMove":
                    return moves[i]
        except:
            pass

        FILE.write("Cannot find any Replica Exchange moves (SireMove::RepExMove) objects to analyse!\n")
        sys.exit(0)

# get the replica exchange move object and print out the acceptance ratio 
repex = getRepExMoves(moves)
FILE.write("\nReplica exchange moves: %s accepted, %s attempted, acceptance ratio = %.1f %%\n" % \
              (repex.nAccepted(), repex.nAttempted(), 100 * repex.acceptanceRatio()) )

# Now get the replica exchange history
try:
    hist = system.lambdaTrajectoryHistory()
except:
    FILE.write("Cannot extract the lambda trajectory history from the RETI system file!\n")
    sys.exit(0)

nits = len(hist)

if range_start:
    if range_start > nits-1:
        start = nits-1
    else:
        start = range_start

    if range_end > nits-1:
        end = nits-1
    else:
        end = range_end
else:
    end = nits-1
    start = end - int(percent * end / 100.0)

FILE.write("\nPrinting lambda trajectory history from iterations %s to %s\n\n" % (start+1, end+1))

for i in range(start, end+1):
    FILE.write("%d " % (i+1))
    for val in hist[i]:
        FILE.write("%s " % val)
    FILE.write("\n")
