"""
Installation script for Sire

This assumes that the python that is used to execute this script
is the conda / miniconda / miniforge environment into
which you want to install Sire.

USAGE:

    python setup.py install_requires   : Will install all of the dependencies

    python setup.py build              : Will install requires and will then
                                         compile sire (takes a long time!)

    python setup.py install            : Will build sire and will then install

    python setup.py install_module     : Will only install the Python module

You can use `--skip-deps` to skip the installation of the conda dependencies
You can use `--skip-build` to skip the building of the corelib and wrappers
"""

import sys
import os
import platform
import subprocess
import shutil
import glob

# Debug - we need to print out all of the environment variables
# for key, value in os.environ.items():
#     print(f"{key}\n{value}\n")

# We can only run this script from the sire directory
curdir = os.path.abspath(".")

if os.path.abspath(os.path.dirname(sys.argv[0])) != curdir:
    print("You can only run this script from the sire directory")
    sys.exit(-1)

# We need to import the 'parse_requirements' function to get the list
# of requirements
sys.path.insert(0, os.path.join(curdir, "actions"))
from parse_requirements import parse_requirements

# Next we need to verify that this is a Python that is part of a
# conda installation

# Find the path to the conda or mamba executable
conda_base = os.path.abspath(os.path.dirname(sys.executable))

if os.path.basename(conda_base) == "bin":
    conda_base = os.path.dirname(conda_base)

python_exe = None
conda = None

if "CONDA_EXE" in os.environ:
    conda = os.environ["CONDA_EXE"]
else:
    conda = None

if "CONDA_DEFAULT_ENV" in os.environ:
    conda_env = os.environ["CONDA_DEFAULT_ENV"]
else:
    conda_env = None

if os.path.exists(os.path.join(conda_base, "python.exe")):
    # Windows
    conda_bin = os.path.join(conda_base, "Library", "bin")
    python_exe = os.path.join(conda_base, "python.exe")

    if conda is None:
        conda = os.path.join(conda_base, "Scripts", "conda.exe")
elif os.path.exists(os.path.join(conda_base, "bin", "python")):
    # MacOS and Linux
    conda_bin = os.path.join(conda_base, "bin")
    python_exe = os.path.join(conda_bin, "python")

    if conda is None:
        conda = os.path.join(conda_bin, "conda")
else:
    print("Cannot find a 'python' binary in directory '%s'. "
          "Are you running this script using the python executable "
          "from a valid miniconda or anaconda installation?" % conda_base)
    sys.exit(-1)


def find_mamba():
    """Find mamba"""
    if conda.endswith(".exe"):
        m = os.path.join(os.path.dirname(conda), "mamba.exe")
    else:
        m = os.path.join(os.path.dirname(conda), "mamba")

    if os.path.exists(m):
        return m
    else:
        return None


mamba = find_mamba()

# Get the build operating system and processor
is_linux = False
is_windows = False
is_macos = False
platform_name = platform.system()
machine = platform.machine()

if machine in ["aarch64", "arm64"]:
    machine = "arm64"
elif machine in ["x86_64", "AMD64"]:
    machine = "x86_64"
else:
    print(f"Unrecognised architecture ({machine}). Compile at your own risk.")

exe_suffix = ""
if platform_name == "Linux":
    platform_name = "linux"
    is_linux = True
elif platform_name == "Darwin":
    platform_name = "macos"
    is_macos = True
elif platform_name == "Windows":
    platform_name = "windows"
    exe_suffix = ".exe"
    is_windows = True
else:
    print("Unrecognised build platform: %s" % platform_name)
    print("We cannot compile sire.")
    sys.exit(-1)

platform_string = f"{platform_name}_{machine}"


def parse_args():
    import argparse
    import multiprocessing
    parser=argparse.ArgumentParser()

    ncores = multiprocessing.cpu_count()

    if ncores % 2 == 0:
        npycores = int(ncores / 2)
    else:
        npycores = int((ncores + 1) / 2)

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
        metavar=("N_CORES",), default=[ncores],
        help="Number of CPU cores used for compiling corelib "
        "(defaults to all available on the current system)")
    parser.add_argument("-N", "--npycores", action="store", type=int, nargs=1,
        metavar=("N_PYTHON_CORES",), default=[npycores],
        help="Number of CPU cores used for compiling Python wrappers "
        "(defaults to the number of CPU cores used for compiling corelib)")
    parser.add_argument("--install-bss-deps", action="store_true", default=False,
        help="Install BioSimSpace's dependencies too. This helps ensure "
             "compatibility between Sire's and BioSimSpace's dependencies.")
    parser.add_argument("--skip-deps", action="store_true", default=False,
        help="Skip the installation of the dependencies (only use if you know "
             "that they are already installed)")
    parser.add_argument("--skip-build", action="store_true", default=False,
        help="Skip the build of the C++ code (only use if you know that "
             "the C++ code is already built)")
    parser.add_argument("action", nargs="*",
        help="Should be one of 'install_requires', 'build', 'install' or 'install_module'.\n"
             "\n [install_requires] : Just install the conda dependencies.\n"
             " [build] : 'install_requires' plus compile and install corelib, and just compile the wrappers.\n"
             " [install] : 'build' plus install the wrappers and install the module.\n"
             " [install_module] : Just install the module (no compilation or conda dependencies)."
             )
    return parser.parse_args()


_installed_deps = None

def _get_installed(conda: str):
    """Return the list of installed conda dependencies"""
    global _installed_deps

    if _installed_deps is None:
        p = subprocess.Popen([conda, "list"], stdout=subprocess.PIPE)
        _installed_deps = str(p.stdout.read())

    return _installed_deps


def is_installed(dep: str, conda: str) -> bool:
    """Return whether or not the passed dependency is installed"""
    installed = _get_installed(conda=conda)
    return installed.find(dep) != -1


def _add_to_dependencies(dependencies, lines):
    import re

    for line in lines:
        line = line.lstrip().rstrip()

        words = re.split("[<>]*=", line)

        if len(words) > 0:
            package = words[0]
            dependencies[package] = line


_is_conda_prepped = False

def conda_install(dependencies, install_bss_reqs=False):
    """Install the passed list of dependencies using conda"""

    conda_exe = conda

    global _is_conda_prepped

    if not _is_conda_prepped:
        if install_bss_reqs:
            cmd = "%s config --prepend channels michellab/label/dev" % conda_exe
            print("Activating michellab channel channel using: '%s'" % cmd)
            status = subprocess.run(cmd.split())
            if status.returncode != 0:
                print("Failed to add michellab channel!")
                sys.exit(-1)

        print("\nSetting channel priorities to favour conda-forge")
        cmd = "%s config --prepend channels conda-forge" % conda_exe
        print("Activating conda-forge channel using: '%s'" % cmd)
        status = subprocess.run(cmd.split())
        if status.returncode != 0:
            print("Failed to add conda-forge channel!")
            sys.exit(-1)

        cmd = "%s config --set channel_priority strict" % conda_exe
        print("Setting channel priority to strict using: '%s'" % cmd)
        status = subprocess.run(cmd.split())
        if status.returncode != 0:
            print("Failed to set channel priority!")
            sys.exit(-1)

        _is_conda_prepped = True

    if mamba is None:
        conda_install = [conda, "install", "--yes"]
    else:
        conda_install = [mamba, "install", "--yes"]

    cmd = [*conda_install, *dependencies]
    print("\nInstalling packages using: '%s'" % " ".join(cmd))
    status = subprocess.run(cmd)

    if status.returncode != 0:
        if mamba is not None:
            # try with conda, as mamba was broken?
            conda_install = [conda, "install", "--yes"]
            cmd = [*conda_install, *dependencies]
            print("\nTrying again using: '%s'" % " ".join(cmd))
            status = subprocess.run(cmd)

        if status.returncode != 0:
            print("Something went wrong installing dependencies!")
            sys.exit(-1)


def install_requires(install_bss_reqs=False):
    """Installs all of the dependencies. This can safely be called
       multiple times, as it will cache the result to prevent future
       installs taking too long
    """
    print(f"Installing requirements for {platform_string}")

    if not os.path.exists(conda):
        print("\nSire can only be installed into a conda or miniconda environment.")
        print("Please install conda, miniconda, miniforge or similar, then "
              "activate the conda environment, then rerun this installation "
              "script.")
        sys.exit(-1)

    reqs = parse_requirements("requirements.txt")
    build_reqs = parse_requirements("requirements_build.txt")

    if install_bss_reqs:
        bss_reqs = parse_requirements("requirements_bss.txt")
        reqs = reqs + bss_reqs

    print("\nUsing dependencies:")
    dependencies = build_reqs + reqs

    for dep in dependencies:
        print(dep)

    global mamba

    if mamba is None:
        # install mamba first!
        conda_install(["mamba"])
        mamba = find_mamba()

    conda_install(dependencies, install_bss_reqs)


def add_default_cmake_defs(cmake_defs, ncores):
    for a in ("ANACONDA_BASE=%s" % conda_base.replace("\\", "/"),
              "BUILD_NCORES=%s" % ncores):
        found = False
        for d in cmake_defs:
            if (a in d[0]):
                found = True
                break
        if (not found):
            cmake_defs.append([a])

    if is_macos:
        # don't compile with AVX as the resulting binaries won't
        # work on M1 macs
        cmake_defs.append(["SIRE_DISABLE_AVX=ON"])
        cmake_defs.append(["SIRE_DISABLE_AVX512F=ON"])


def make_cmd(ncores, install = False):
    if is_windows:
        action = "INSTALL" if install else "ALL_BUILD"
        make_args = "%s -- /m:%s /p:Configuration=Release /p:Platform=x64" % (action, ncores)
    else:
        action = "install" if install else "all"
        make_args = "%s -- VERBOSE=1 -j %s" % (action, ncores)

    return make_args.split()


def _get_build_ext():
    if "CONDA_BUILD" in os.environ and os.environ["CONDA_BUILD"] == "1":
        return "conda_build"
    else:
        if conda_env is not None:
            ext = "_" + conda_env.replace(" ", "_").replace(".", "_")
        else:
            ext = ""

        return os.path.basename(conda_base.replace(" ", "_").replace(".", "_")) + ext


def _get_bin_dir():
    if "CONDA_BUILD" in os.environ and os.environ["CONDA_BUILD"] == "1":
        bindir = os.environ["BUILD_PREFIX"]

        if is_windows:
            return os.path.join(bindir, "Library", "bin")
        else:
            return os.path.join(bindir, "bin")
    else:
        return conda_bin


def build(ncores: int = 1, npycores: int = 1,
          coredefs=[], pydefs=[]):
    print("\nCompiling the C++ code")

    CC=None
    CXX=None

    cmake = "cmake"

    if is_windows:
        cmake = f"{cmake}.exe"

    bindir = _get_bin_dir()
    cmake = os.path.join(bindir, cmake)

    if "CONDA_BUILD" in os.environ and os.environ["CONDA_BUILD"] == "1":
        conda_build = True
    else:
        os.environ["CONDA_BUILD"] = "0"
        conda_build = False

    # get the compilers
    if conda_build:
        print("This is a conda build")
        CXX = os.environ["CXX"]
        CC = os.environ["CC"]

        # make sure that these compilers are in the path
        CXX_bin = shutil.which(CXX)
        CC_bin = shutil.which(CC)

        print(f"{CXX} => {CXX_bin}")
        print(f"{CC} => {CC_bin}")

        if CXX_bin is None or CC_bin is None:
            print("Cannot find the compilers requested by conda-build in the PATH")
            print("Please check that the compilers are installed and available.")
            sys.exit(-1)

        # use the full paths, in case CMake struggles
        CXX = CXX_bin
        CC = CC_bin

    elif is_macos:
        try:
            CXX = glob.glob(os.path.join(bindir, "clang++"))[0]
            CC = glob.glob(os.path.join(bindir, "clang"))[0]
        except:
            print("Cannot find the conda clang++ binaries!")
            print("Please install these, e.g. via")
            print("conda install clang clangxx")
            sys.exit(-1)
    elif is_linux:
        try:
            CXX = glob.glob(os.path.join(bindir, "*-g++"))[0]
            CC = glob.glob(os.path.join(bindir, "*-gcc"))[0]
        except:
            print("Cannot find the conda g++ binaries!")
            print("Please install these, e.g. via")
            print("conda install gcc gxx")
            sys.exit(-1)

    print(f"Using compilers {CC} | {CXX}")

    # Make sure all of the above output is printed to the screen
    # before we start running any actual compilation
    sys.stdout.flush()

    # Now that the dependencies are installed, the next step
    # is to use cmake to build the corelib and wrapper in the build/corelib
    # and build/wrapper directories

    # change into the build/corelib directory
    OLDPWD = os.getcwd()

    build_ext = _get_build_ext()

    build_dir = os.path.abspath("build")

    coredir = os.path.join(build_dir, f"{build_ext}_corelib")

    if not os.path.exists(coredir):
        os.makedirs(coredir)

    if not os.path.isdir(coredir):
        print("SOMETHING IS WRONG. %s is not a directory?" % coredir)
        sys.exit(-1)

    os.chdir(coredir)

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

        add_default_cmake_defs(args.corelib, ncores)

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
        print("\n== OUTPUT LOG ==")
        try:
            for line in open("CMakeFiles/CMakeOutput.log").readlines():
                print(line.strip())
        except Exception as e:
            print(e)

        print("\n== ERROR LOG ==")
        try:
            for line in open("CMakeFiles/CMakeError.log").readlines():
                print(line.strip())
        except Exception as e:
            print(e)

        sys.exit(-1)

    # Now that cmake has run, we can compile and install corelib
    ######
    ###### Compiling and installing corelib
    ######
    # Compile and install, as need to install to compile the wrappers
    make_args = make_cmd(ncores, True)

    print("NOW RUNNING \"%s\" --build . --target %s" % (cmake, " ".join(make_args)))
    sys.stdout.flush()
    status = subprocess.run([cmake, "--build", ".", "--target", *make_args])

    if status.returncode != 0:
        print("SOMETHING WENT WRONG WHEN COMPILING CORELIB!")
        sys.exit(-1)

    ######
    ###### Compiling wrapper
    ######

    # Make sure that the Python in conda is used
    pydefs.append([f"PYTHON_EXECUTABLE={python_exe}"])

    os.chdir(OLDPWD)

    wrapperdir = os.path.join(build_dir, f"{build_ext}_wrapper")

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

        add_default_cmake_defs(args.wrapper, npycores)

        cmake_cmd = [cmake, *sum([["-D", d[0]] for d in args.wrapper], []),
                     *sum([["-G", g[0]] for g in args.generator], []),
                     sourcedir]

        print(" ".join(cmake_cmd))
        sys.stdout.flush()
        status = subprocess.run(cmake_cmd)

    if status.returncode != 0:
        print("SOMETHING WENT WRONG WHEN USING CMAKE ON WRAPPER!")
        sys.exit(-1)

    # Just compile the wrappers
    make_args = make_cmd(npycores, False)

    print("NOW RUNNING \"%s\" --build . --target %s" % (cmake, " ".join(make_args)))
    sys.stdout.flush()
    status = subprocess.run([cmake, "--build", ".", "--target", *make_args])

    if status.returncode != 0:
        print("SOMETHING WENT WRONG WHEN COMPILING THE WRAPPERS!")
        sys.exit(-1)

    os.chdir(OLDPWD)


def install_module(ncores: int = 1):
    print("\nInstalling the module")

    OLDPWD = os.getcwd()

    build_ext = _get_build_ext()

    build_dir = os.path.abspath("build")

    cmake = "cmake"

    if is_windows:
        cmake = f"{cmake}.exe"

    bindir = _get_bin_dir()
    cmake = os.path.join(bindir, cmake)

    moduledir = os.path.join(build_dir, f"{build_ext}_module")

    if not os.path.exists(moduledir):
        os.makedirs(moduledir)

    if not os.path.isdir(moduledir):
        print("SOMETHING IS WRONG. %s is not a directory?" % moduledir)
        sys.exit(-1)

    os.chdir(moduledir)

    if os.path.exists("CMakeCache.txt"):
        # we have run cmake in this directory before. Run it again.
        status = subprocess.run([cmake, "."])
    else:
        # this is the first time we are running cmake
        sourcedir = os.path.join(os.path.dirname(os.path.dirname(
            os.getcwd())), "src", "sire")

        if not os.path.exists(os.path.join(sourcedir, "CMakeLists.txt")):
            print("SOMETHING IS WRONG. There is no file %s" % os.path.join(sourcedir, "CMakeLists.txt"))
            sys.exit(-1)

        add_default_cmake_defs(args.wrapper, ncores)
        cmake_cmd = [cmake, *sum([["-D", d[0]] for d in args.wrapper], []),
                     *sum([["-G", g[0]] for g in args.generator], []),
                     sourcedir]
        print(" ".join(cmake_cmd))
        sys.stdout.flush()
        status = subprocess.run(cmake_cmd)

    if status.returncode != 0:
        print("SOMETHING WENT WRONG WHEN USING CMAKE ON MODULE!")
        sys.exit(-1)

    make_args = make_cmd(ncores, True)

    # Now that cmake has run, we can compile and install wrapper
    print("NOW RUNNING \"%s\" --build . --target %s" % (cmake, " ".join(make_args)))
    sys.stdout.flush()
    status = subprocess.run([cmake, "--build", ".", "--target", *make_args])

    if status.returncode != 0:
        print("SOMETHING WENT WRONG WHEN COMPILING WRAPPER!")
        sys.exit(-1)


def install(ncores: int = 1, npycores: int = 1):
    print("\nInstalling sire")

    OLDPWD = os.getcwd()

    build_ext = _get_build_ext()

    build_dir = os.path.abspath("build")

    wrapperdir = os.path.join(build_dir, f"{build_ext}_wrapper")

    cmake = "cmake"

    if is_windows:
        cmake = f"{cmake}.exe"

    bindir = _get_bin_dir()
    cmake = os.path.join(bindir, cmake)

    if not os.path.exists(wrapperdir):
        os.makedirs(wrapperdir)

    if not os.path.isdir(wrapperdir):
        print("SOMETHING IS WRONG. %s is not a directory?" % wrapperdir)
        sys.exit(-1)

    os.chdir(wrapperdir)

    # Now install the wrappers
    make_args = make_cmd(npycores, True)

    print("NOW RUNNING \"%s\" --build . --target %s" % (cmake, " ".join(make_args)))
    sys.stdout.flush()
    status = subprocess.run([cmake, "--build", ".", "--target", *make_args])

    if status.returncode != 0:
        print("SOMETHING WENT WRONG WHEN COMPILING THE WRAPPERS!")
        sys.exit(-1)

    os.chdir(OLDPWD)

    install_module(ncores=ncores)


if __name__ == "__main__":
    args = parse_args()

    if len(args.action) != 1:
        print("Please use either 'install_requires', 'build' or 'install'")
        sys.exit(-1)

    install_bss = args.install_bss_deps

    action = args.action[0]

    if is_windows and (args.generator is None or len(args.generator) == 0):
        #args.generator = [["Visual Studio 17 2022"]]
        args.generator = [["Visual Studio 15 2017 Win64"]]  # preferred VC version for conda

    if action == "install":
        if not (args.skip_deps or args.skip_build):
            install_requires(install_bss_reqs=install_bss)

        if not args.skip_build:
            build(ncores=args.ncores[0], npycores=args.npycores[0],
                  coredefs=args.corelib, pydefs=args.wrapper)

        install(ncores=args.ncores[0], npycores=args.npycores[0])

    elif action == "build":
        if not args.skip_deps:
            install_requires(install_bss_reqs=install_bss)

        build(ncores=args.ncores[0], npycores=args.npycores[0],
              coredefs=args.corelib, pydefs=args.wrapper)

    elif action == "install_requires":
        install_requires(install_bss_reqs=install_bss)

    elif action == "install_module":
        install_module(ncores=args.ncores[0])

    else:
        print(f"Unrecognised action '{action}'. Please use 'install_requires', "
              "'build', 'install' or 'install_module'")
