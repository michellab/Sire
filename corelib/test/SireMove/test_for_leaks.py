import Sire.Stream
from Sire.Units import *

import sys

print("Loading %s..." % sys.argv[1])
(replicas, moves) = Sire.Stream.load( sys.argv[1] )

replica0 = replicas[0]

print("Unpacking the replica...")
replica0.unpack()

system = replica0.subSystem()

print("Calculating the old energies...")
old_nrgs = system.energies()

system.mustNowRecalculateFromScratch()

print("Calculating the new energies...")
new_nrgs = system.energies()

keys = list(old_nrgs.keys())
keys.sort()

all_ok = True

print("Checking the energies are unchanged...")

for key in keys:
    old_nrg = old_nrgs[key]
    new_nrg = new_nrgs[key]

    diff_nrg = new_nrg - old_nrg

    if ( (diff_nrg < -0.001) or (diff_nrg > 0.001) ):
        print("ENERGY DISCREPANCY IN %s: %s vs. %s (difference = %s)" % \
              (key, old_nrg, new_nrg, diff_nrg))

        all_ok = False

if all_ok:
    print("ALL ENERGIES OK - CODE IS WORKING")
