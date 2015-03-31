#!/usr/bin/env python

import sys

for filename in sys.argv:
   lines = open(filename,"r").readlines()

   found = False

   for line in lines:
      if line.find("SIRE_END_HEADER") != -1:
          found = True
          break

   if not found:
       print "%s does not contain SIRE_END_HEADER" % filename

