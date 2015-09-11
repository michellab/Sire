
import os
import sys

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

for dir in dirs:
    create_wrappers(dir)

#os.chdir("Qt")
#os.system("%s create_qt_wrappers.py" % sys.executable)
#os.chdir("..")
