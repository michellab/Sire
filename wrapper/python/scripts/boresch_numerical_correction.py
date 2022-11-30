description="""
An app to calculate the correction for releasing Boresch restraints
using numerical integration - see Equation 12 in
J. Phys. Chem. B 2003, 107, 35, 9535–9551. This is more robust than
the analytical correction, which can result in substantial errors in
certain regimes."""

try:
    import sire
    sire.use_old_api()
except ImportError:
    pass

from Sire.Tools import BoreschNumericalCorrection
from Sire.Tools import readParams
import Sire.Config

import argparse
import sys
import os

parser = argparse.ArgumentParser(description="An app to calculate the correction for releasing Boresch restraints "
                                             "using numerical integration - see Equation 12 in J. Phys. Chem. B 2003, 107,"
                                              "35, 9535–9551. This is more robust than the analytical correction, which"
                                              "can result in substantial errors in certain regimes.",
                                epilog="boresch_numerical_correction is built using Sire and is distributed "
                                        "under the GPL. For more information please visit "
                                        "http://siremol.org",
                                 prog="boresch_numerical_correction")

parser.add_argument('-C', '--config', nargs="?",
                    help='A config file used for the simulation (only the restraints '
                         'dictionary and temperature are used to calculate the correction).')

parser.add_argument('--author', action="store_true",
                    help="Get information about the authors of this script.")

parser.add_argument('--version', action="store_true",
                    help="Get version information about this script.")

parser.add_argument('--verbose', action="store_true",
                    help="Print debug output.")

sys.stdout.write("\n")
args = parser.parse_args()

must_exit = False

if args.author:
    print("\nboresch_numerical_correction was written by Finlay Clark (C) 2022, but is very closely"
          "\nbased on standardstatecorrection, written by Julien Michel and Stefano Bosisio (C) 2017")
    must_exit = True

if args.version:
    print("boresch_numerical_correction -- from Sire release version <%s>" %Sire.__version__)
    print("This particular release can be downloaded here: "
          "https://github.com/michellab/Sire/releases/tag/v%s" %Sire.__version__)
    must_exit = True

if must_exit:
    sys.exit(0)

if args.config:
    config_file = args.config
    if not os.path.exists(config_file):
        print("(cannot find configuration file %s)" % config_file)
        sys.exit(-1)
    else:
        print("Loading temperature and restraints information from file %s" % args.config)
        params = readParams(args.config)
else:
    print("Please supply a configuration file")
    sys.exit(-1)

if args.verbose:
    params["verbose"] = True

print("\nCalculating numerical correction for Boresch restraints using temperature and "
      "restraint information from %s." % args.config)

BoreschNumericalCorrection.run(params)
