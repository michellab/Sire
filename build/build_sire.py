## This script is used to finish the compile of Sire corelib and wrappers
## This script must be called using the python from a miniconda or
## anaconda distribution

import sys
import os
import time

from multiprocessing import Pool

def print_progress():
    while True:
        print("Build is in progress...")
        time.sleep(300)

if __name__ == "__main__":

    pool = Pool()

    # run a background process that prints to the screen so that people know that we are still alive...
    result = pool.apply_async(print_progress)

    conda_base = os.path.abspath( os.path.dirname(sys.executable) )

    if conda_base.endswith("bin"):
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

    # Get the number of cores to use for compiling the python wrappers - this 
    # defaults to NCORES
    if "NPYCORES" in os.environ:
        NPYCORES = int(os.environ["NPYCORES"])
    else:
        NPYCORES = NCORES

    print("Number of cores used for compilation = %d" % NCORES)

    if NPYCORES != NCORES:
        print("Number of cores used for wrapper compilation = %d" % NPYCORES)

    python_exe = None
    conda_exe = None

    if os.path.exists("%s/bin/conda" % conda_base):
        python_exe = "%s/bin/python" % conda_base
        conda_exe = "%s/bin/conda" % conda_base
    elif os.path.exists("%s\python.exe" % conda_base):
        python_exe = "%s/python.exe" % conda_base
        conda_exe = "%s/Scripts/conda.exe" % conda_base
    else:
        print("Cannot find a 'conda' binary in directory '%s'. "
              "Are you running this script using the python executable "
              "from a valid miniconda or anaconda installation?" % conda_base) 
        sys.exit(-1)

    print("Continuing the Sire install using %s %s" \
              % (python_exe,sys.argv[0]))

    # now go through all of the python modules that need to be available
    # into this conda installation, and make sure they have been installed

    # first, pip
    try:
        import pip
        print("pip is already installed...")
    except:
        print("Installing pip using '%s install pip'" % conda_exe)
        os.system("%s install --yes pip" % conda_exe)

    # ipython
    try:
        import IPython
        print("ipython is already installed...")
    except:
        print("Installing ipython using %s install ipython" % conda_exe)
        os.system("%s install --yes ipython" % conda_exe)

    # pytest
    try:
        import pytest
        print("pytest is already installed...")
    except:
        print("Installing pytest using %s install pytest" % conda_exe)
        os.system("%s install --yes pytest" % conda_exe)

    # nose
    try:
        import nose
        print("nose is already installed...")
    except:
        print("Installing nose using '%s install nose'" % conda_exe)
        os.system("%s install --yes nose" % conda_exe)

    # openmm
    try:
        import simtk.openmm
        print("openmm is already installed...")
    except:
        print("Installing openmm from the conda-forge repository...")
        os.system("%s conda install -c omnia -c conda-forge openmm" % conda_exe)
        #os.system("%s install --yes openmm=7.1" % conda_exe)

    # libnetcdf
    try:
        import netCDF4
        print("netCDF4 is already installed...")
    except:
        print("Installing netCDF4 using '%s install netcdf4'" % conda_exe)
        os.system("%s install --yes netcdf4" % conda_exe)

    # Now that the miniconda distribution is ok, the next step
    # is to use cmake to build the corelib and wrapper in the build/corelib
    # and build/wrapper directories

    # first, get the value of the CXX environment variable
    cxx = os.getenv("CXX")
    compiler_ext = None

    if cxx:
        if cxx.find("icpc") != -1:
            compiler_ext = "intel"

    # change into the build/corelib directory
    OLDPWD = os.path.abspath(os.curdir)

    if compiler_ext:
        coredir = "%s/corelib_%s" % (build_dir,compiler_ext)
    else:
        coredir = "%s/corelib" % build_dir

    if not os.path.exists(coredir):
        os.makedirs(coredir)

    if not os.path.isdir(coredir):
        print("SOMETHING IS WRONG. %s is not a directory?" % coredir)
        sys.exit(-1)

    os.chdir(coredir)

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
    status = os.system("make -j %s install" % NCORES)

    if status != 0:
        print("SOMETHING WENT WRONG WHEN INSTALLING CORELIB!")
        sys.exit(-1)

    # Ok, that is all complete. Next we must work on the
    # python wrappers
    os.chdir(OLDPWD)

    if compiler_ext:
        wrapperdir = "%s/wrapper_%s" % (build_dir,compiler_ext)
    else:
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
    status = os.system("make -j %s" % NPYCORES)

    if status != 0:
        print("SOMETHING WENT WRONG WHEN COMPILING WRAPPER!")
        sys.exit(-1)

    # Now the compilation has finished, install wrapper
    status = os.system("make -j %s install" % NPYCORES)

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
    except Exception as e:
        print("Something went wrong when trying to test the Sire installation.")
        print(e)
        print("Please check things manually yourself...")

    print("\n\n=================================")
    print("Congratulations. Everything has installed :-)")

    sys.exit(0)
