
import os

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
    os.system("python ../AutoGenerate/create_wrappers.py")
    os.chdir("..")

for dir in dirs:
    create_wrappers(dir)

#os.chdir("Qt")
#os.system("python create_qt_wrappers.py")
#os.chdir("..")
