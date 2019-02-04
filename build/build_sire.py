## This script is used to finish the compile of Sire corelib and wrappers
## This script must be called using the python from a miniconda or
## anaconda distribution

import sys
import os
import time
import argparse
import multiprocessing


def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("-C", "--corelib", action="append", nargs=1,
        metavar=("PARAMETER=VALUE",), default=[],
        help="pass CMake definitions for corelib")
    parser.add_argument("-W", "--wrapper", action="append", nargs=1,
        metavar=("PARAMETER=VALUE",), default=[],
        help="pass CMake definitions for wrapper")
    parser.add_argument("-G", dest="generator", action="store", nargs=1,
        metavar=("GENERATOR",), default="",
        help="pass CMake generator")
    parser.add_argument("-n", "--ncores", action="store", type=int, nargs=1,
        metavar=("N_CORES",), default=multiprocessing.cpu_count(),
        help="Number of CPU cores used for compiling corelib "
        "(defaults to all available on the current system)")
    parser.add_argument("-N", "--npycores", action="store", type=int, nargs=1,
        metavar=("N_PYTHON_CORES",), default=-1,
        help="Number of CPU cores used for compiling Python wrappers "
        "(defaults to the number of CPU cores used for compiling corelib)")
    parser.add_argument("--pip", action="store_true",
        help="Use pip rather than conda to install dependencies")
    return parser.parse_args()


def make_cmd(ncores, install = False):
    try:
        make = os.environ["MAKE"]
    except KeyError:
        make = None
        with open("CMakeCache.txt", "r") as hnd:
            line = True
            while line:
                line = hnd.readline()
                if (line and line.startswith("CMAKE_GENERATOR:INTERNAL")
                    and line.split("=")[-1].startswith("Visual Studio")):
                    make = "MSBuild.exe"
                    action = "INSTALL" if install else "BUILD_ALL"
                    make_args = "/m:%s /p:Configuration=Release /p:Platform=x64 %s.vcxproj" % (ncores, action)
                    break
        if (make is None):
            make = "make"
    if ("make" in make):
        action = " install" if install else ""
        make_args = "-j %s%s" % (ncores, action)
    return make, make_args

if __name__ == "__main__":
    args = parse_args()
    
    install_script = sys.argv[0]
    build_dir = os.path.abspath( os.path.dirname(install_script) )

    # Get the number of cores to use for any compiling - if this is 0, then 
    # we use all of the cores
    try:
        NCORES = int(os.environ["NCORES"])
    except KeyError:
        NCORES = args.ncores

    # Get the number of cores to use for compiling the python wrappers - this 
    # defaults to NCORES
    try:
        NPYCORES = int(os.environ["NPYCORES"])
    except KeyError:
        NPYCORES = args.npycores if args.npycores > 0 else NCORES

    print("Number of cores used for compilation = %d" % NCORES)

    if NPYCORES != NCORES:
        print("Number of cores used for wrapper compilation = %d" % NPYCORES)

    if args.pip:
        python_exe = sys.executable
        py_module_install = "%s -m pip install" % python_exe
        conda_base = ""
    else:
        conda_base = os.path.abspath(os.path.dirname(sys.executable))

        if os.path.basename(conda_base) == "bin":
            conda_base = os.path.dirname(conda_base)

        python_exe = None
        conda_exe = None

        if os.path.exists(os.path.join(conda_base, "bin", "conda")):
            python_exe = os.path.join(conda_base, "bin", "python")
            conda_exe = os.path.join(conda_base, "bin", "conda")
        elif os.path.exists(os.path.join(conda_base, "python.exe")):
            python_exe = os.path.join(conda_base, "python.exe")
            conda_exe = os.path.join(conda_base, "Scripts", "conda.exe")
        else:
            print("Cannot find a 'conda' binary in directory '%s'. "
                  "Are you running this script using the python executable "
                  "from a valid miniconda or anaconda installation?" % conda_base) 
            sys.exit(-1)
        py_module_install = "%s install --yes" % conda_exe

    print("Continuing the Sire install using %s %s" \
              % (python_exe, sys.argv[0]))

    # now go through all of the python modules that need to be available
    # into this conda installation, and make sure they have been installed

    # first, pip
    try:
        import pip
        print("pip is already installed...")
    except ImportError:
        if args.pip:
            print("Could not import pip - please make sure "
                "pip is available within your current Python environment "
                "then re-run this script.")
            sys.exit(-1)
        else:
            print("Installing pip using '%s install pip'" % conda_exe)
            os.system("%s pip" % py_module_install)

    # ipython
    try:
        import IPython
        print("ipython is already installed...")
    except ImportError:
        print("Installing ipython using '%s ipython'" % py_module_install)
        os.system("%s ipython" % py_module_install)

    # pytest
    try:
        import pytest
        print("pytest is already installed...")
    except ImportError:
        print("Installing pytest using '%s pytest'" % py_module_install)
        os.system("%s pytest" % py_module_install)

    # nose
    try:
        import nose
        print("nose is already installed...")
    except ImportError:
        print("Installing nose using '%s nose'" % py_module_install)
        os.system("%s nose" % py_module_install)

    # openmm
    try:
        import simtk.openmm
        print("openmm is already installed...")
    except ImportError:
        if args.pip:
            print("It looks like the openmm Python modules are not "
                "available - please check your openmm instllation")
            sys.exit(-1)
        else:
            print("Installing openmm from the conda-forge repository...")
            os.system("%s install --yes -c omnia -c conda-forge openmm" % conda_exe)
            #os.system("%s install --yes openmm=7.1" % conda_exe)

    # libnetcdf
    try:
        import netCDF4
        print("netCDF4 is already installed...")
    except ImportError:
        print("Installing netCDF4 using '%s netcdf4'" % py_module_install)
        os.system("%s netcdf4" % py_module_install)

    # Make sure all of the above output is printed to the screen
    # before we start running any actual compilation
    sys.stdout.flush()

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
    OLDPWD = os.getcwd()

    if compiler_ext:
        coredir = os.path.join(build_dir, "corelib_%s" % build_dir)
    else:
        coredir = os.path.join(build_dir, "corelib")

    if not os.path.exists(coredir):
        os.makedirs(coredir)

    if not os.path.isdir(coredir):
        print("SOMETHING IS WRONG. %s is not a directory?" % coredir)
        sys.exit(-1)

    os.chdir(coredir)
    sys.stderr.write("1) I am now here: %s\n" % os.getcwd())
    sys.stderr.flush()

    def add_default_cmake_defs(cmake_defs):
        for a in ("ANACONDA_BUILD=ON", "ANACONDA_BASE=%s" % conda_base, "BUILD_NCORES=%s" % NCORES):
            if (args.pip and a.startswith("ANACONDA")):
                continue
            found = False
            for d in cmake_defs:
                if (a in d[0]):
                    found = True
                    break
            if (not found):
                cmake_defs.append([a])

    if os.path.exists("CMakeCache.txt"):
        # we have run cmake in this directory before. Run it again.
        status = os.system("cmake .")
    else:
        # this is the first time we are running cmake
        sourcedir = os.path.join(os.path.dirname(os.path.dirname(
            os.getcwd())), "corelib")

        if not os.path.exists(os.path.join(sourcedir, "CMakeLists.txt")):
            print("SOMETHING IS WRONG. There is no file %s" % os.path.join(sourcedir, "CMakeLists.txt"))
            sys.exit(-1)

        #status = os.system("cmake -D ANACONDA_BUILD=ON -D ANACONDA_BASE=%s -D BUILD_NCORES=%s %s" \
        #                 % (conda_base,NCORES,sourcedir) )
        for a in ("NetCDF_ROOT_DIR", "OPENMM_ROOT_DIR"):
            for i, d in enumerate(args.corelib):
                if (a in d[0]):
                    v = args.corelib.pop(i)[0]
                    if (not a in os.environ):
                        os.environ[a] = v.split("=")[-1]
        add_default_cmake_defs(args.corelib)
        sys.stderr.write("2) args.corelib: %s\n" % str(args.corelib))
        sys.stderr.flush()
        cmake_cmd = "cmake %s %s" % (" ".join(["-D \"%s\"" % d[0] for d in args.corelib]), sourcedir)
        sys.stderr.write("1) cmake_cmd: %s\n" % cmake_cmd)
        sys.stderr.flush()
        status = os.system(cmake_cmd)

    if status != 0:
        print("SOMETHING WENT WRONG WHEN USING CMAKE ON CORELIB!")
        sys.exit(-1)

    # Now that cmake has run, we can compile corelib
    make, make_args = make_cmd(NCORES)
                
    status = os.system("%s %s" % (make, make_args))

    if status != 0:
        print("SOMETHING WENT WRONG WHEN COMPILING CORELIB!")
        sys.exit(-1)

    # Now the compilation has finished, install corelib
    make, make_args = make_cmd(NCORES, True)

    status = os.system("%s %s" % (make, make_args))

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
        sourcedir = os.path.join(os.path.dirname(os.path.dirname(
            os.getcwd())), "wrapper")

        if not os.path.exists(os.path.join(sourcedir, "CMakeLists.txt")):
            print("SOMETHING IS WRONG. There is no file %s" % os.path.join(sourcedir, "CMakeLists.txt"))
            sys.exit(-1)

        #status = os.system("cmake -D ANACONDA_BUILD=ON -D ANACONDA_BASE=%s -D BUILD_NCORES=%s %s" \
        #                 % (conda_base,NCORES,sourcedir) )
        add_default_cmake_defs(args.wrapper)
        sys.stderr.write("2) args.wrapper: %s\n" % str(args.wrapper))
        sys.stderr.flush()
        cmake_cmd = "cmake %s %s" % (" ".join(["-D \"%s\"" % d[0] for d in args.wrapper]), sourcedir)
        sys.stderr.write("2) cmake_cmd: %s\n" % cmake_cmd)
        sys.stderr.flush()
        status = os.system(cmake_cmd)

    if status != 0: 
        print("SOMETHING WENT WRONG WHEN USING CMAKE ON WRAPPER!")
        sys.exit(-1)

    #Now that cmake has run, we can compile wrapper
    #status = os.system("make -j %s" % NPYCORES)
    make, make_args = make_cmd(NPYCORES)

    status = os.system("%s %s" % (make, make_args))

    if status != 0:
        print("SOMETHING WENT WRONG WHEN COMPILING WRAPPER!")
        sys.exit(-1)

    # Now the compilation has finished, install wrapper
    make, make_args = make_cmd(NPYCORES, True)

    status = os.system("%s %s" % (make, make_args))

    if status != 0:
        print("SOMETHING WENT WRONG WHEN INSTALLING WRAPPER!")
        sys.exit(-1)

    print("\n\n=================================")
    print("Congratulations. Everything has installed :-)")
