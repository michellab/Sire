from Sire.Tools import Nautilus
from Sire.Tools import readParams

import Sire.Config

import argparse
import os
import sys

import ast

parser = argparse.ArgumentParser(description="Nautilus protocol is run with default settings.  Default settings are controlled within the nautilus-protocol.cfg file",
                                 epilog="nautilus-protocol is built using Sire, Numpy and mdtraj and is distributed "
                                        "under the GPL. For more information please visit "
                                        "http://siremol.org/nautilus",
                                 prog="nautilus")

parser.add_argument('--author', action="store_true",
                    help="Get information about the authors of this script.")

parser.add_argument('--version', action="store_true",
                    help="Get version information about this script.")

parser.add_argument('-C', '--config', nargs="?",
                    help='Supply a Nautilus CONFIG file to control the calculation.')

parser.add_argument('-b', '--benchmark', action='store_true',
                    help="Benchmark the Nautilus subroutines.")
sys.stdout.write("\n")

args = parser.parse_args()

must_exit = False

if args.author:
    print("\n nautilus-protocol was written by Georgios Gerogiokas and Julien Michel (C) 2014")
    print("It is based on the Nautilus Sire module.")
    must_exit = True

if args.version:
    print("\n nautilus-protocol version 1.0")
    #print(Sire.Config.versionString())
    must_exit = True

if must_exit:
    sys.exit(0)

params = {}

if args.config:
    print("Loading configuration information from file %s" % args.config)
    params = readParams(args.config)

if args.benchmark:
    params["benchmark"] = True

print ( args.config)
if not (os.path.exists("%s" % args.config)):
    parser.print_help()
    print("\nPlease supply the configuration file, %s." % args.config)
    if not os.path.exists("%s" % args.config):
        print("(cannot find configuration file, %s" % args.config )
    sys.exit(-1)

print("\nRunning nautilus-protocol.py using configuration file, %s" % args.config)

print ("running nautilus-traj2cell")
# traj2cell
cmd = "%s/bin/nautilus-traj2cell -C %s" % (params["sire"], args.config)
print (cmd)
os.system(cmd)

print ("running nautilus-cell2grid")
# cell2grid
cmd = "%s/bin/nautilus-cell2grid -C %s " % (params["sire"], args.config)
print (cmd)
os.system(cmd)

print ("running nautilus-clustergrids")
# clustergrid
cmd = ''' %s/bin/nautilus-clustergrids -C %s ''' % (params["sire"], args.config)
print (cmd)
os.system(cmd)
