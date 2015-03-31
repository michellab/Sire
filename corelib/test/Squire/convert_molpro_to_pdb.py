
## Converts a Molpro input file into a PDB so 
## that I can check that the cutoffs are implemented
## correctly

import sys
import os
import re

molpro = open(sys.argv[1], "r")

line = molpro.readline()

while (line):
      m = re.match(r"\s*(\d+)\s*!\s*number of atoms", line)
      
      if (m):
          natoms = int(m.group(1))
          print("natoms = %d" % natoms)
          break
          
      line = molpro.readline()
      
npdb = 0
      
#read in the next natoms lines
while (line and natoms > 0):
      line = molpro.readline()
      
      m = re.match(r"\s*(\w+),\s*([\-0-9\.]+)\s*,\s*([\-0-9\.]+)\s*,\s*([\-0-9\.]+)",
                   line)
                   
      if m:
            natoms = natoms - 1
            npdb = npdb + 1

            print("ATOM %6d  %3s %4s %4d   %8.3f%8.3f%8.3f" % \
                     (npdb, m.group(1), "QM", 1, float(m.group(2)),
                                                 float(m.group(3)),
                                                 float(m.group(4))))

print("TER")

while line:
    line = molpro.readline()

    if re.match("BEGIN_DATA",line):
        break
        
line = molpro.readline()

while line:
    if re.match("END",line):
        break
        
      
    m = re.match(r"\s*([\-0-9\.]+)\s*,\s*([\-0-9\.]+)\s*,\s*([\-0-9\.]+)\s*,\s*([\-0-9\.]+)",
                 line)
                 
    if m:
        #must convert MM coordinates from bohr radii to angstrom
        npdb = npdb + 1
        print("ATOM %6d  %3s %4s %4d   %8.3f%8.3f%8.3f" % \
                     (npdb, "DUM", "MM", 2, float(m.group(1)) * 0.529177249,
                                            float(m.group(2)) * 0.529177249,
                                            float(m.group(3)) * 0.529177249))

    line = molpro.readline()

print("END")
