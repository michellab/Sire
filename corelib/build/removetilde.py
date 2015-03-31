#!/bin/env python

import os
import string

def getTildeFiles(arg, dirname, fnames):
    for file in fnames:
        if (file.endswith("~")):
            arg.append(dirname + "/" + file)
            
files = []
os.path.walk(".", getTildeFiles, files)

for file in files:
    print "Removing %s..." % file
    os.remove(file)
