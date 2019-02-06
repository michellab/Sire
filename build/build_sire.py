## This script is used to finish the compile of Sire corelib and wrappers
## This script must be called using the python from a miniconda or
## anaconda distribution

import sys
import os
import glob
import time
import platform

if __name__ == "__main__":

    is_linux = False
    is_windows = False
    is_osx = False

    if platform.system() == "Linux":
        is_linux = True
        print("Compiling on Linux")
    elif platform.system() == "Darwin":
        is_osx = True
        print("Compiling on OS X")
    elif platform.system() == "Windows":
        print("Sorry - compiling into miniconda on Windows is not supported yet")
        is_windows = True
        sys.exit(-1)
    else:
        print("Unrecognised build platform: %s" % platform.system())
        sys.exit(-1)

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
        # default to half the number to save memory
        NPYCORES = NCORES / 2

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

    conda_pkgs = []

    # first, pip
    try:
        import pip
        print("pip is already installed...")
    except:
        conda_pkgs.append("pip")

    # ipython
    try:
        import IPython
        print("ipython is already installed...")
    except:
        conda_pkgs.append("ipython")

    # pytest
    try:
        import pytest
        print("pytest is already installed...")
    except:
        conda_pkgs.append("pytest")

    # nose
    try:
        import nose
        print("nose is already installed...")
    except:
        conda_pkgs.append("nose")

    # boost
    if os.path.exists("%s/include/boost/python.hpp" % conda_base):
        print("boost is already installed...")
    else:
        conda_pkgs.append("boost")

    # gsl
    if os.path.exists("%s/include/gsl/gsl_version.h" % conda_base):
        print("gsl is already installed...")
    else:
        conda_pkgs.append("gsl")

    # tbb
    if os.path.exists("%s/include/tbb/tbb.h" % conda_base):
        print("TBB is already installed...")
    else:
        conda_pkgs.append("tbb")
        conda_pkgs.append("tbb-devel")

    # Qt5
    try:
        import PyQt5
        print("Qt5 is already installed...")
    except:
        conda_pkgs.append("pyqt")

    # libnetcdf
    try:
        import netCDF4
        print("netCDF4 is already installed...")
    except:
        conda_pkgs.append("netcdf4")

    # compilers (so we keep binary compatibility
    if is_osx:
        try:
            CXX = glob.glob("%s/bin/clang++" % conda_base)[0]
            CC = glob.glob("%s/bin/clang" % conda_base)[0]
            print("clang++ is already installed...")
        except:
            conda_pkgs.append("clang_osx-64")
            conda_pkgs.append("clangxx_osx-64")
    elif is_linux:
        try:
            CXX = glob.glob("%s/bin/*-g++" % conda_base)[0]
            CC = glob.glob("%s/bin/*-gcc" % conda_base)[0]
        except:
            conda_pkgs.append("gcc_linux-64")
            conda_pkgs.append("gxx_linux-64")

    if os.path.exists("%s/bin/make" % conda_base):
        print("make is already installed...")
    else:
        conda_pkgs.append("make")

    make = "%s/bin/make" % conda_base

    if os.path.exists("%s/bin/cmake" % conda_base):
        print("cmake is already installed...")
    else:
        conda_pkgs.append("cmake")

    cmake = "%s/bin/cmake" % conda_base

    installed_something = False

    if len(conda_pkgs) > 0:
        cmd = "%s install --yes %s" % (conda_exe, " ".join(conda_pkgs))
        print("Installing packages using '%s'" % cmd)
        status = os.system(cmd)
        installed_something = True
        if status != 0:
            print("Something went wrong installing dependencies!")
            sys.exit(-1)

    # openmm last as different repo
    try:
        import simtk.openmm
        print("openmm is already installed...")
    except:
        print("Installing openmm from the conda-forge repository...")
        os.system("%s install --yes -c omnia -c conda-forge openmm" % conda_exe)
        #os.system("%s install --yes openmm=7.1" % conda_exe)
        installed_something = True

    if installed_something:
        # need to fix numpy - breaks because of MKL blas!
        # see https://github.com/numpy/numpy/issues/11481
        os.system("%s uninstall --yes --force blas" % conda_exe)
        os.system("%s install --yes \"blas=*=openblas\"" % conda_exe)
        os.system("%s update --yes --force numpy" % conda_exe)

    # make sure we really have found the compilers
    if is_osx:
        try:
            CXX = glob.glob("%s/bin/clang++" % conda_base)[0]
            CC = glob.glob("%s/bin/clang" % conda_base)[0]
            print("clang++ is already installed...")
        except:
            print("Cannot find the conda clang++ binaries!")
            sys.exit(-1)
    elif is_linux:
        try:
            CXX = glob.glob("%s/bin/*-g++" % conda_base)[0]
            CC = glob.glob("%s/bin/*-gcc" % conda_base)[0]
        except:
            print("Cannot find the conda g++ binaries!")
            sys.exit(-1)

    print("Using compilers %s | %s" % (CC, CXX))

    # Make sure all of the above output is printed to the screen
    # before we start running any actual compilation
    sys.stdout.flush()

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

    if os.path.exists("CMakeCache.txt"):
        # we have run cmake in this directory before. Run it again.
        status = os.system("%s ." % cmake)
    else:
        # this is the first time we are running cmake
        sourcedir = os.path.abspath("../../corelib")

        if not os.path.exists("%s/CMakeLists.txt" % sourcedir):
            print("SOMETHING IS WRONG. There is no file %s/CMakeLists.txt" % coredir)
            sys.exit(-1)

        cmd = "CC=%s CXX=%s %s -D ANACONDA_BUILD=ON -D ANACONDA_BASE=%s -D BUILD_NCORES=%s %s" \
                         % (CC,CXX,cmake,conda_base,NCORES,sourcedir) 
        print(cmd)
        status = os.system(cmd)

    if status != 0:
        print("SOMETHING WENT WRONG WHEN USING CMAKE ON CORELIB!")
        sys.exit(-1)

    # Now that cmake has run, we can compile and install corelib
    status = os.system("%s -j %s install" % (make,NCORES))

    if status != 0:
        print("SOMETHING WENT WRONG WHEN COMPILING CORELIB!")
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
        status = os.system("%s ." % cmake)
    else:
        # this is the first time we are running cmake
        sourcedir = os.path.abspath("../../wrapper")   

        if not os.path.exists("%s/CMakeLists.txt" % sourcedir):
            print("SOMETHING IS WRONG. There is no file %s/CMakeLists.txt" % wrapperdir)
            sys.exit(-1)
    
        cmd = "CC=%s CXX=%s %s -D ANACONDA_BUILD=ON -D ANACONDA_BASE=%s -D BUILD_NCORES=%s %s" \
                         % (CC,CXX,cmake,conda_base,NCORES,sourcedir)
        print(cmd)
        status = os.system(cmd)

    if status != 0: 
        print("SOMETHING WENT WRONG WHEN USING CMAKE ON WRAPPER!")
        sys.exit(-1)

    # Now that cmake has run, we can compile wrapper
    status = os.system("%s -j %s install" % (make,NPYCORES))

    if status != 0:
        print("SOMETHING WENT WRONG WHEN COMPILING WRAPPER!")
        sys.exit(-1)

    print("\n\n=================================")
    print("Congratulations. Everything has installed :-)")
