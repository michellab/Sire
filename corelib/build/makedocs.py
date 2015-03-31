#!/bin/env python

import sys
import os

runmodules = []

if len(sys.argv) > 1:
    runmodules = sys.argv[1:]

# Relative path to the api doc directory
docdir = "docs/api"

# get the top-level directory (directory to build docs from)
topdir = os.getcwd()

#get the template doxyfile
doxyfile = "%s/Doxyfile" % topdir

#definition of the Module class
class Module:
    def __init__(self,name,indir,outdir):
        self._name = name
        self._id = outdir
        self._indir = indir
        self._outdir = "../../../%s/%s" % (docdir,outdir)
        self._tagfile = "%s/%s.tag" % (self._outdir,self._name)
        self._html = "../../../%s/html" % outdir
        
    def name(self):
        return self._name

    def ID(self):
        return self._id

    def indir(self):
        return self._indir
        
    def outdir(self):
        return self._outdir

    def tagfile(self):
        return self._tagfile

    def html(self):
        return self._html

# List the modules whose documentation should
# be compiled, together with the location of the 
# source code, and where to place the docs
modules = [
            Module("SireStream", "src/libs/SireStream", "libs/SireStream"),
            Module("SireError", "src/libs/SireError","libs/SireError"),
            Module("SireUnits", "src/libs/SireUnits","libs/SireUnits"),
            Module("SireMaths", "src/libs/SireMaths","libs/SireMaths"),
            Module("SireCluster", "src/libs/SireCluster","libs/SireCluster"),
            Module("SireCAS", "src/libs/SireCAS","libs/SireCAS"),
            Module("SireBase", "src/libs/SireBase", "libs/SireBase"),
            Module("SireMol", "src/libs/SireMol", "libs/SireMol"),
            Module("SireVol", "src/libs/SireVol", "libs/SireVol"),
            Module("SireDB", "src/libs/SireDB", "libs/SireDB"),
            Module("SireFF", "src/libs/SireFF", "libs/SireFF"),
            Module("SireMM", "src/libs/SireMM", "libs/SireMM"),
            Module("SireIO", "src/libs/SireIO", "libs/SireIO"),
            Module("SireSystem", "src/libs/SireSystem", "libs/SireSystem"),
            Module("SireMove", "src/libs/SireMove", "libs/SireMove"),
            Module("SireSim", "src/libs/SireSim", "libs/SireSim"),
            Module("Spier", "src/libs/Spier", "libs/Spier"),
            Module("Squire", "src/libs/Squire", "libs/Squire"),
            Module("SireTest", "src/libs/SireTest", "libs/SireTest"),
            Module("SireUnitTest", "src/libs/SireUnitTest", "libs/SireUnitTest"),
            Module("SirePy", "src/libs/SirePy", "libs/SirePy"),
           
            Module("Spier", "src/apps/spier", "apps/spier"),
            Module("siretest", "src/apps/siretest", "apps/siretest")
          ]

#index all of the modules
modulehash = {}

for module in modules:
    modulehash[module.ID()] = module

#first, parse the Doxyfile and save each key-value pair
template = file(doxyfile,"r").readlines()
          
#variable used to save the paths and names of all of the
#generated tag files
tagfiles = {}

#save the mapping from a modules tagfile to the location
#of the html docs
for module in modules:
    tagfiles[module.tagfile()] = module.html()

def removeDir(directory):
    
    try:
      os.remove(directory)
    except OSError:
      for filename in os.listdir(directory):
        removeDir("%s/%s" % (directory,filename))
        
      os.rmdir(directory)

def processModule(module):
    #change to the top-level directory
    os.chdir(topdir)
    
    print "Making documentation for %s" % (module.name())
    
    #change into the input directory
    os.chdir(module.indir())
    
    #create the output directory
    try:
        os.makedirs(module.outdir())
    except OSError:
        print "Directory %s already exists!" % module.outdir()
        #remove the existing directory
        removeDir(module.outdir())
    
    print "input from %s, output to %s" % (os.getcwd(),module.outdir())

    pipe = os.popen("doxygen -","w")
    #pipe = file("test.txt","w")
    
    #write the template configuration commands
    for line in template:
       pipe.write(line)
    
    #now write the local-template configuration commands
    try:
       localtmpl = file("Doxyfile","r").readlines()
       
       for line in localtmpl:
           pipe.write(line)

    except IOError:
       #there are no local options for this module
       nolocal = True
    
    #now write the local configuration for this module
    pipe.write("PROJECT_NAME = %s\n" % module.name())
    pipe.write("INPUT = .\n")
    pipe.write("OUTPUT_DIRECTORY = %s\n" % module.outdir())
    pipe.write("HTML_HEADER = ../../../build/apidocs/header.html\n")
    pipe.write("HTML_FOOTER = ../../../build/apidocs/footer.html\n")
   
    #what is the name of this modules tag file?
    pipe.write("GENERATE_TAGFILE = %s\n" % module.tagfile())
    
    #link to existing tag files
    pipe.write("TAGFILES = ")
       
    for taglink in tagfiles:
       pipe.write(" \\ \n      %s=%s" % (taglink,tagfiles[taglink]))
           
    pipe.write("\n")
    
    pipe.close()

if len(runmodules) > 0:
    for module in runmodules:
        processModule(modulehash[module])
        
else:
    for module in modules:
        processModule(module)
