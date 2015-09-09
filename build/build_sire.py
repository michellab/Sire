## This script is used to finish the compile of Sire corelib and wrappers
## This script must be called using the python from a miniconda or
## anaconda distribution

import sys
import os

conda_base = os.path.abspath( "%s/.." % os.path.dirname(sys.executable))
install_script = sys.argv[0]
build_dir = os.path.abspath( os.path.dirname(install_script) )

# Get the number of cores to use for any compiling - if this is 0, then 
# we use all of the cores
if "NCORES" in os.environ:
    NCORES = int(os.environ["NCORES"])
else:
    import multiprocessing
    NCORES = multiprocessing.cpu_count()

print("Number of cores used for compilation = %d" % NCORES)

if not os.path.exists("%s/bin/conda" % conda_base):
    print("Cannot find a 'conda' binary in directory '%s/bin'. "
          "Are you running this script using the python executable "
          "from a valid miniconda or anaconda installation?" % conda_base) 
    sys.exit(-1)

print("Continuing the Sire install using %s/bin/python %s" \
          % (conda_base,sys.argv[0]))

# now go through all of the python modules that need to be available
# into this conda installation, and make sure they have been installed

# first, pip
try:
    import pip
    print("pip is already installed...")
except:
    print("Installing pip using '%s/bin/conda install pip'" % conda_base)
    os.system("%s/bin/conda install --yes pip" % conda_base)

# ipython
try:
    import IPython
    print("ipython is already installed...")
except:
    print("Installing ipython using '%s/bin/conda install ipython'" % conda_base)
    os.system("%s/bin/conda install --yes ipython" % conda_base)

# nose
try:
    import nose
    print("nose is already installed...")
except:
    print("Installing nose using '%s/bin/conda install nose'" % conda_base)
    os.system("%s/bin/conda install --yes nose" % conda_base)

# openmm
try:
    import simtk.openmm
    print("openmm is already installed...")
except:
    print("Installing openmm from the omnia repository...")
    os.system("%s/bin/conda config --add channels http://conda.binstar.org/omnia" % conda_base)
    os.system("%s/bin/conda install --yes openmm" % conda_base)

# Now that the miniconda distribution is ok, the next step
# is to use cmake to build the corelib and wrapper in the build/corelib
# and build/wrapper directories

# change into the build/corelib directory
OLDPWD = os.path.abspath(os.curdir)

coredir = "%s/corelib" % build_dir

if not os.path.exists(coredir):
    os.makedirs(coredir)

if not os.path.isdir(coredir):
    print("SOMETHING IS WRONG. %s is not a directory?" % coredir)
    sys.exit(-1)

os.chdir(coredir)

os.system("which cmake")
os.system("locate cmake")

if os.path.exists("CMakeCache.txt"):
    # we have run cmake in this directory before. Run it again.
    status = os.system("cmake .")
else:
    # this is the first time we are running cmake
    sourcedir = os.path.abspath("../../corelib")

    if not os.path.exists("%s/CMakeLists.txt" % sourcedir):
        print("SOMETHING IS WRONG. There is no file %s/CMakeLists.txt" % coredir)
        sys.exit(-1)

    status = os.system("cmake -D ANACONDA_BUILD=ON -D ANACONDA_BASE=%s -D BUILD_NCORES=%s %s" \
                     % (conda_base,NCORES,sourcedir) )

if status != 0:
    print("SOMETHING WENT WRONG WHEN USING CMAKE ON CORELIB!")
    sys.exit(-1)

# Now that cmake has run, we can compile corelib
status = os.system("make -j %s" % NCORES)

if status != 0:
    print("SOMETHING WENT WRONG WHEN COMPILING CORELIB!")
    sys.exit(-1)

# Now the compilation has finished, install corelib
status = os.system("make -j %s install/strip" % NCORES)

if status != 0:
    print("SOMETHING WENT WRONG WHEN INSTALLING CORELIB!")
    sys.exit(-1)

# Ok, that is all complete. Next we must work on the
# python wrappers
os.chdir(OLDPWD)

wrapperdir = "%s/wrapper" % build_dir
    
if not os.path.exists(wrapperdir):
    os.makedirs(wrapperdir)

if not os.path.isdir(wrapperdir):
    print("SOMETHING IS WRONG. %s is not a directory?" % wrapperdir)
    sys.exit(-1)

os.chdir(wrapperdir)

if os.path.exists("CMakeCache.txt"):
    # we have run cmake in this directory before. Run it again.
    status = os.system("cmake .")
else:
    # this is the first time we are running cmake
    sourcedir = os.path.abspath("../../wrapper")   

    if not os.path.exists("%s/CMakeLists.txt" % sourcedir):
        print("SOMETHING IS WRONG. There is no file %s/CMakeLists.txt" % wrapperdir)
        sys.exit(-1)
    
    status = os.system("cmake -D ANACONDA_BUILD=ON -D ANACONDA_BASE=%s -D BUILD_NCORES=%s %s" \
                     % (conda_base,NCORES,sourcedir) )

if status != 0: 
    print("SOMETHING WENT WRONG WHEN USING CMAKE ON WRAPPER!")
    sys.exit(-1)

# Now that cmake has run, we can compile wrapper
status = os.system("make -j %s" % NCORES)

if status != 0:
    print("SOMETHING WENT WRONG WHEN COMPILING WRAPPER!")
    sys.exit(-1)

# Now the compilation has finished, install wrapper
status = os.system("make -j %s install/strip" % NCORES)

if status != 0:
    print("SOMETHING WENT WRONG WHEN INSTALLING WRAPPER!")
    sys.exit(-1)

# Now that everything has been installed, we should be able 
# to import Sire
try:
    import Sire.CAS
    x = Sire.CAS.Symbol("x")
    f = x**2 + 5*x - 10
    g = f.differentiate(x)
except:
    print("Something went wrong when trying to test the Sire installation.")
    print("Please check things manually yourself...")
    sys.exit(-1)

print("=================================")
print("Congratulations. Everything has installed :-)")
