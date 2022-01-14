## This script is used to finish the compile of Sire corelib and wrappers
## This script must be called using the python from a miniconda or
## anaconda distribution

import sys
import sysconfig
import os
import glob
import time
import argparse
import subprocess
import multiprocessing
import platform


def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("-C", "--corelib", action="append", nargs=1,
        metavar=("PARAMETER=VALUE",), default=[],
        help="pass CMake definitions for corelib")
    parser.add_argument("-W", "--wrapper", action="append", nargs=1,
        metavar=("PARAMETER=VALUE",), default=[],
        help="pass CMake definitions for wrapper")
    parser.add_argument("-G", "--generator", action="append", nargs=1,
        metavar=("GENERATOR",), default=[],
        help="pass CMake generator")
    parser.add_argument("-n", "--ncores", action="store", type=int, nargs=1,
        metavar=("N_CORES",), default=multiprocessing.cpu_count(),
        help="Number of CPU cores used for compiling corelib "
        "(defaults to all available on the current system)")
    parser.add_argument("-N", "--npycores", action="store", type=int, nargs=1,
        metavar=("N_PYTHON_CORES",), default=-1,
        help="Number of CPU cores used for compiling Python wrappers "
        "(defaults to the number of CPU cores used for compiling corelib)")
    parser.add_argument("--noconda", action="store_true",
        help="Use pip rather than conda to install dependencies")
    parser.add_argument("--no-openmm", action="store_true",
        help="Compile Sire without OpenMM support")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    is_linux = False
    is_windows = False
    is_mingw = False
    is_osx = False

    exe_suffix = ""
    if platform.system() == "Linux":
        is_linux = True
        print("Compiling on Linux")
    elif platform.system() == "Darwin":
        is_osx = True
        print("Compiling on OS X")
    elif platform.system() == "Windows":
        exe_suffix = ".exe"
        #print("Sorry - compiling into miniconda on Windows is not supported yet")
        #args.noconda = True
        is_mingw = (sysconfig.get_platform() == "mingw")
        is_windows = (sysconfig.get_platform().startswith("win"))
    else:
        print("Unrecognised build platform: %s" % platform.system())
        sys.exit(-1)

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
        # default to half the number to save memory
        NPYCORES = args.npycores if args.npycores > 0 else NCORES

    print("Number of cores used for compilation = %d" % NCORES)

    if NPYCORES != NCORES:
        print("Number of cores used for wrapper compilation = %d" % NPYCORES)

    if args.noconda:
        python_exe = sys.executable
        py_module_install = [python_exe, "-m", "pip", "install"]
        conda_base = ""
    else:
        conda_base = os.path.abspath(os.path.dirname(sys.executable))

        if os.path.basename(conda_base) == "bin":
            conda_base = os.path.dirname(conda_base)

        python_exe = None
        conda_exe = None

        if os.path.exists(os.path.join(conda_base, "bin", "conda")):
            conda_bin = os.path.join(conda_base, "bin")
            python_exe = os.path.join(conda_bin, "python")
            conda_exe = os.path.join(conda_bin, "conda")
        elif os.path.exists(os.path.join(conda_base, "python.exe")):
            conda_bin = os.path.join(conda_base, "Library", "bin")
            python_exe = os.path.join(conda_base, "python.exe")
            conda_exe = os.path.join(conda_base, "Scripts", "conda.exe")
        else:
            print("Cannot find a 'conda' binary in directory '%s'. "
                  "Are you running this script using the python executable "
                  "from a valid miniconda or anaconda installation?" % conda_base)
            sys.exit(-1)
        py_module_install = [conda_exe, "install", "--yes"]

    print("Continuing the Sire install using %s %s" \
              % (python_exe, sys.argv[0]))

    # now go through all of the python modules that need to be available
    # into this conda installation, and make sure they have been installed

    conda_pkgs = []

    # first, pip
    try:
        import pip
        print("pip is already installed...")
    except ImportError:
        if args.noconda:
            print("Could not import pip - please make sure "
                "pip is available within your current Python environment "
                "then re-run this script.")
            sys.exit(-1)
        else:
            print("Installing pip using '%s install pip'" % conda_exe)
            conda_pkgs.append("pip")

    # ipython
    try:
        import IPython
        print("ipython is already installed...")
    except ImportError:
        conda_pkgs.append("ipython")

    # pytest
    try:
        import pytest
        print("pytest is already installed...")
    except ImportError:
        conda_pkgs.append("pytest")

    # nose
    try:
        import nose
        print("nose is already installed...")
    except ImportError:
        conda_pkgs.append("nose")

    # libnetcdf
    try:
        import netCDF4
        print("netCDF4 is already installed...")
    except ImportError:
        # had to step back from 1.5.8 as this isn't available
        # yet on Linux (works on MacOS)
        conda_pkgs.append("netcdf4=1.5.7")

    CC = None
    CXX = None
    cmake = "cmake%s" % exe_suffix
    if (not args.noconda):
        # We pin Conda packages to the highest version that is available on
        # Linux, macOS, and Windows-x64.

        # boost
        if os.path.exists(os.path.join(conda_base, "include", "boost", "python.hpp")):
            print("boost is already installed...")
        else:
            conda_pkgs.append("boost=1.74.0")

        # gsl
        if os.path.exists(os.path.join(conda_base, "include", "gsl", "gsl_version.h")):
            print("gsl is already installed...")
        else:
            conda_pkgs.append("gsl=2.7")

        # tbb
        if os.path.exists(os.path.join(conda_base, "include", "tbb", "tbb.h")):
            print("TBB is already installed...")
        else:
            conda_pkgs.append("tbb=2021.5.0")
            conda_pkgs.append("tbb-devel=2021.5.0")

        # Qt5
        try:
            import PyQt5
            print("Qt5 is already installed...")
        except ImportError:
            conda_pkgs.append("pyqt=5.12.3")
            # This is the version for Apple M1
            #conda_pkgs.append("pyqt=5.15.2")

        # compilers (so we keep binary compatibility)
        if is_osx:
            try:
                CXX = glob.glob(os.path.join(conda_bin, "clang++"))[0]
                CC = glob.glob(os.path.join(conda_bin, "clang"))[0]
                print("clang++ is already installed...")
            except:
                conda_pkgs.append("clang_osx-64")
                conda_pkgs.append("clangxx_osx-64")
        elif is_linux:
            try:
                CXX = glob.glob(os.path.join(conda_bin, "*-g++"))[0]
                CC = glob.glob(os.path.join(conda_bin, "*-gcc"))[0]
            except:
                conda_pkgs.append("gcc_linux-64")
                conda_pkgs.append("gxx_linux-64")

                # trying not to fix these, as these cause issues later...
                conda_pkgs.append("libgfortran5")
                conda_pkgs.append("libgcc-ng")

        if (not is_windows):
            if os.path.exists(os.path.join(conda_bin, "make")):
                print("make is already installed...")
            else:
                conda_pkgs.append("make")
            make = os.path.join(conda_bin, "make")
            if os.path.exists(os.path.join(conda_bin, "libtool")):
                print("libtool is already installed...")
            else:
                conda_pkgs.append("libtool")
            if os.path.exists(os.path.join(conda_bin, "autoreconf")):
                print("autoconf is already installed...")
            else:
                conda_pkgs.append("autoconf")
            if os.path.exists(os.path.join(conda_bin, "aclocal")):
                print("automake is already installed...")
            else:
                conda_pkgs.append("automake")

        if os.path.exists(os.path.join(conda_bin, "cmake%s" % exe_suffix)):
            print("cmake is already installed...")
        else:
            conda_pkgs.append("cmake")

        cmake = os.path.join(conda_bin, "cmake%s" % exe_suffix)

    installed_something = False

    if (not args.noconda) and conda_pkgs:
        # Do we still need to do this - my miniconda looks in conda-forge
        # already. Maybe this is what is slowing the environment processing?
        #cmd = "%s config --prepend channels conda-forge" % conda_exe
        #print("Activating conda-forge channel using: '%s'" % cmd)
        #status = subprocess.run(cmd.split())
        #if status.returncode != 0:
        #    print("Failed to add conda-forge channel!")
        #    sys.exit(-1)

        # Need to run this command to prevent conda errors on
        # some platforms - see
        # https://github.com/ContinuumIO/anaconda-issues/issues/11246
        #cmd = "%s config --set channel_priority false" % conda_exe
        #print("Setting channel priority to false using: '%s'" % cmd)
        #status = subprocess.run(cmd.split())
        #if status.returncode != 0:
        #    print("Failed to set channel priority!")
        #    sys.exit(-1)

        cmd = [*py_module_install, *conda_pkgs]
        print("Installing packages using: '%s'" % " ".join(cmd))
        status = subprocess.run(cmd)
        installed_something = True
        if status.returncode != 0:
            print("Something went wrong installing dependencies!")
            sys.exit(-1)

    # check if the user wants to link against openmm
    use_openmm = True

    if args.no_openmm:
        print("DISABLING OPENMM SUPPORT")
        use_openmm = False

    for d in args.corelib:
        if ("SIRE_USE_OPENMM=OFF" in d[0]):
            use_openmm = False
            break
    # openmm last as different repo
    if use_openmm:
        try:
            import simtk.openmm
            print("openmm is already installed...")
        except ImportError:
            if args.noconda:
                print("It looks like the openmm Python modules are not "
                    "available - please check your openmm installation")
                sys.exit(-1)
            else:
                print("Installing openmm from the Omnia channel...")
                subprocess.run(("%s install --yes openmm=7.7.0" % conda_exe).split())
                installed_something = True

    # make sure we really have found the compilers
    if (not args.noconda):
        if is_osx:
            try:
                CXX = glob.glob(os.path.join(conda_bin, "clang++"))[0]
                CC = glob.glob(os.path.join(conda_bin, "clang"))[0]
                print("clang++ is already installed...")
            except:
                print("Cannot find the conda clang++ binaries!")
                sys.exit(-1)
        elif is_linux:
            try:
                CXX = glob.glob(os.path.join(conda_bin, "*-g++"))[0]
                CC = glob.glob(os.path.join(conda_bin, "*-gcc"))[0]
            except:
                print("Cannot find the conda g++ binaries!")
                sys.exit(-1)

    print("Using compilers %s | %s" % (CC, CXX))

    # Make sure all of the above output is printed to the screen
    # before we start running any actual compilation
    sys.stdout.flush()

    # Now that the miniconda distribution is ok, the next step
    # is to use cmake to build the corelib and wrapper in the build/corelib
    # and build/wrapper directories

    # change into the build/corelib directory
    OLDPWD = os.getcwd()

    coredir = os.path.join(build_dir, "corelib")

    if not os.path.exists(coredir):
        os.makedirs(coredir)

    if not os.path.isdir(coredir):
        print("SOMETHING IS WRONG. %s is not a directory?" % coredir)
        sys.exit(-1)

    os.chdir(coredir)

    def add_default_cmake_defs(cmake_defs):
        for a in ("ANACONDA_BUILD=ON", "ANACONDA_BASE=%s" % conda_base.replace("\\", "/"), "BUILD_NCORES=%s" % NCORES):
            if (args.noconda and a.startswith("ANACONDA")):
                continue
            found = False
            for d in cmake_defs:
                if (a in d[0]):
                    found = True
                    break
            if (not found):
                cmake_defs.append([a])

    def make_cmd(ncores, install = False):
        if is_windows:
            action = "INSTALL" if install else "ALL_BUILD"
            make_args = "%s -- /m:%s /p:Configuration=Release /p:Platform=x64" % (action, ncores)
        else:
            action = "install" if install else ""
            make_args = "%s -- VERBOSE=1 -j %s" % (action, ncores)
        return make_args.split()

    if os.path.exists("CMakeCache.txt"):
        # we have run cmake in this directory before. Run it again.
        status = subprocess.run([cmake, "."])
    else:
        # this is the first time we are running cmake
        sourcedir = os.path.join(os.path.dirname(os.path.dirname(
            os.getcwd())), "corelib")

        if not os.path.exists(os.path.join(sourcedir, "CMakeLists.txt")):
            print("SOMETHING IS WRONG. There is no file %s" % os.path.join(sourcedir, "CMakeLists.txt"))
            sys.exit(-1)

        for a in ("NetCDF_ROOT_DIR", "OPENMM_ROOT_DIR"):
            for i, d in enumerate(args.corelib):
                if (a in d[0]):
                    v = args.corelib.pop(i)[0]
                    if (not a in os.environ):
                        os.environ[a] = v.split("=")[-1]
        add_default_cmake_defs(args.corelib)
        cmake_cmd = [cmake, *sum([["-D", d[0]] for d in args.corelib], []),
                     *sum([["-G", g[0]] for g in args.generator], []),
                     sourcedir]
        if (CC):
            os.environ["CC"] = CC
        if (CXX):
            os.environ["CXX"] = CXX
        print(" ".join(cmake_cmd))
        sys.stdout.flush()
        status = subprocess.run(cmake_cmd)

    if status.returncode != 0:
        print("SOMETHING WENT WRONG WHEN USING CMAKE ON CORELIB!")
        sys.exit(-1)

    # Now that cmake has run, we can compile and install corelib
    make_args = make_cmd(NCORES, True)

    print("NOW RUNNING \"%s\" --build . --target %s" % (cmake, " ".join(make_args)))
    sys.stdout.flush()
    status = subprocess.run([cmake, "--build", ".", "--target", *make_args])

    if status.returncode != 0:
        print("SOMETHING WENT WRONG WHEN COMPILING CORELIB!")
        sys.exit(-1)

    # Ok, that is all complete. Next we must work on the
    # python wrappers
    os.chdir(OLDPWD)

    wrapperdir = os.path.join(build_dir, "wrapper")

    if not os.path.exists(wrapperdir):
        os.makedirs(wrapperdir)

    if not os.path.isdir(wrapperdir):
        print("SOMETHING IS WRONG. %s is not a directory?" % wrapperdir)
        sys.exit(-1)

    os.chdir(wrapperdir)

    if os.path.exists("CMakeCache.txt"):
        # we have run cmake in this directory before. Run it again.
        status = subprocess.run([cmake, "."])
    else:
        # this is the first time we are running cmake
        sourcedir = os.path.join(os.path.dirname(os.path.dirname(
            os.getcwd())), "wrapper")

        if not os.path.exists(os.path.join(sourcedir, "CMakeLists.txt")):
            print("SOMETHING IS WRONG. There is no file %s" % os.path.join(sourcedir, "CMakeLists.txt"))
            sys.exit(-1)

        add_default_cmake_defs(args.wrapper)
        cmake_cmd = [cmake, *sum([["-D", d[0]] for d in args.wrapper], []),
                     *sum([["-G", g[0]] for g in args.generator], []),
                     sourcedir]
        print(" ".join(cmake_cmd))
        sys.stdout.flush()
        status = subprocess.run(cmake_cmd)

    if status.returncode != 0:
        print("SOMETHING WENT WRONG WHEN USING CMAKE ON WRAPPER!")
        sys.exit(-1)

    make_args = make_cmd(NPYCORES, True)

    # Now that cmake has run, we can compile and install wrapper
    print("NOW RUNNING \"%s\" --build . --target %s" % (cmake, " ".join(make_args)))
    sys.stdout.flush()
    status = subprocess.run([cmake, "--build", ".", "--target", *make_args])

    if status.returncode != 0:
        print("SOMETHING WENT WRONG WHEN COMPILING WRAPPER!")
        sys.exit(-1)

    print("\n\n=================================")
    print("Congratulations. Everything has installed :-)")
