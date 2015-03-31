#!/bin/env python

import sys
import os

args = sys.argv[1:]

files = args[0:-1]
newdir = args[-1]

for file in files:
    cmd = "svn mv %s %s/" % (file,newdir)

    print cmd
    os.system(cmd)
    
