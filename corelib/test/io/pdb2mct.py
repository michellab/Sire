#!/bin/env python

def convertPDB2MCT(filename):
    """Convert the file 'filename' from PDB format to MCT format"""
    
    f = file(filename,"r")
    lines = f.readlines()
    f.close()
    
    print("%MCT-version_1.0")
    print("INFO produced_by pdb2mct.py, and will be broken!")
    
    oldmol = 0
    newmol = True
    oldres = 0
    
    for line in lines:
        if (line[0:4] == "ATOM" or line[0:6] == "HETATM"):
            #get the info for this atom
            words = line.split()
        
            #make sure we are working on a molecule
            if (newmol):
                oldmol = oldmol + 1
                newmol = False
                print("molecule %d filename" % oldmol)

            #see if the residue number has changed
            if (oldres != int(words[4])):
                oldres = int(words[4])
                print("residue %d %s" % (oldres,words[3]))
            
            #get the element name
            elsym = "?"
            try:
                elsym = int(words[2][0:1])
                elsym = words[2][1:2]
            except:
                elsym = words[2][0:1]
            
            #now print out the atom info
            print("atom %4s %1s %8.3f %8.3f %8.3f" % (words[2],elsym, \
                                       float(words[5]),float(words[6]),float(words[7])))
        elif (line[0:3] == "TER"):
            newmol = True
        elif (line[0:3] == "END"):
            return


if (__name__ == "__main__"):

    import sys

    args = sys.argv

    if (len(args) < 2):
        print("USAGE: pdb2mct.py file.pdb > file.mct")
        sys.exit(0)
    
    convertPDB2MCT(args[1])
