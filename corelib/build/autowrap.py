#!/bin/env python

import os
import sys

def writeWrapper(classname, name, returns, args, isconst):

    if (name.startswith("&")):
        name = name[1:]
        returns = returns + "&"
        
    if (isconst):
        print "        .def( \"%s\", (const %s (%s::*)(%s))" \
                         % (name,returns,classname,args)
        print "                        &%s::%s," %(classname,name)
        print "                        return_value_policy<copy_const_reference>() )"
    else:
        print "        .def( \"%s\", (%s (%s::*)(%s))" \
                         % (name,returns,classname,args)
        print "                        &%s::%s )" % (classname,name)

def process(classname, line):
    """Turn a function declaration into a boost python wrapper"""
    
    if (words[0] == "const"):
        returns = words[1]
        name = words[2].split("(")[0]
        rest = line[ line.find(name) + len(name) +1 : -3]

        writeWrapper(classname,name,returns,rest,True)
    else:
        returns = words[0]
        name = words[1].split("(")[0]
        rest = line[ line.find(name) + len(name) +1 : -3]
        
        writeWrapper(classname,name,returns,rest,False)

lines = file(sys.argv[1],"r").readlines()

readheader = True

classname = "Unknown"

for line in lines:

    if (readheader):
        sys.stdout.write(line)

        if line.find("class_") != -1:
             classname = line.split("<")[1].split(",")[0].split(">")[0]
             readheader = False

    else:        
        #now we are reading the raw functions to be converted, 
        #one per line
        words = line.split()

        if (len(words) <= 1 or words[0].startswith(".def")):
            sys.stdout.write(line)
        else:
            process(classname, line)
        
