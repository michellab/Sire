
import os
import sys

from multiprocessing import Pool

dirs = [ "Analysis", \
         "Base", \
         "CAS", \
         "Cluster", \
         "FF", \
         "ID", \
         "IO", \
         "MM", \
         "Maths", \
         "Mol", \
         "Move", \
         "Squire", \
         "Stream", \
         "System", \
         "Units", \
         "Vol" ]

def create_wrappers(dir): 
    os.chdir(dir)
    os.system("%s ../AutoGenerate/create_wrappers.py" % sys.executable)
    os.chdir("..")

if __name__ == "__main__":
    pool = Pool()

    pool.map( create_wrappers, dirs )

    #os.chdir("Qt")
    #os.system("%s create_qt_wrappers.py" % sys.executable)
    #os.chdir("..")

    # restore this function, as it doesn't change, and
    # rewrapping causes an obsolete gamma function to be exposed
    os.system("git checkout Maths/_Maths_free_functions.pypp.cpp")
