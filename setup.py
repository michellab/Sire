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
"""

import sys
import os
import shutil
import platform
import subprocess

# We can only run this script from the sire directory
curdir = os.path.abspath(".")

if os.path.abspath(os.path.dirname(sys.argv[0])) != curdir:
    print("You can only run this script from the sire directory")
    sys.exit(-1)

# Next we need to verify that this is a Python that is part of a
# conda installation

# Find the path to the conda or mamba executable
conda_base = os.path.abspath(os.path.dirname(sys.executable))

if os.path.basename(conda_base) == "bin":
    conda_base = os.path.dirname(conda_base)

python_exe = None
conda = None

if os.path.exists(os.path.join(conda_base, "bin", "conda")):
    conda_bin = os.path.join(conda_base, "bin")
    python_exe = os.path.join(conda_bin, "python")
    conda = os.path.join(conda_bin, "conda")
elif os.path.exists(os.path.join(conda_base, "python.exe")):
    conda_bin = os.path.join(conda_base, "Library", "bin")
    python_exe = os.path.join(conda_base, "python.exe")
    conda = os.path.join(conda_base, "Scripts", "conda.exe")
else:
    print("Cannot find a 'conda' binary in directory '%s'. "
          "Are you running this script using the python executable "
          "from a valid miniconda or anaconda installation?" % conda_base)
    sys.exit(-1)

if not os.path.exists(conda):
    print("\nSire can only be installed into a conda or miniconda environment.")
    print("Please install conda, miniconda, miniforge or similar, then "
          "activate the conda environment, then rerun this installation "
          "script.")
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
        metavar=("N_PYTHON_CORES",), default=multiprocessing.cpu_count(),
        help="Number of CPU cores used for compiling Python wrappers "
        "(defaults to the number of CPU cores used for compiling corelib)")
    parser.add_argument("--skip-deps", action="store_true", default=False,
        help="Skip the installation of the dependencies (only use if you know "
             "that they are already installed)")
    parser.add_argument("--skip-build", action="store_true", default=False,
        help="Skip the build of the C++ code (only use if you know that "
             "the C++ code is already built)")
    parser.add_argument("action", nargs="*",
        help="Should be one of 'install_requires', 'build' or 'install'")
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

def conda_install(dependencies):
    """Install the passed list of dependencies using conda"""
    conda_exe = conda

    global _is_conda_prepped

    if not _is_conda_prepped:
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


def install_requires():
    """Installs all of the dependencies. This can safely be called
       multiple times, as it will cache the result to prevent future
       installs taking too long
    """
    print(f"Installing requirements for {platform_string}")

    dependencies = {}

    lines = open("requirements.txt", "r").readlines()
    _add_to_dependencies(dependencies, lines)

    # now read the platform-specific requirements (these may overload
    # some of the above requirements)
    try:
        lines = open(f"requirements_{platform_name}.txt", "r").readlines()
        _add_to_dependencies(dependencies, lines)
    except IOError:
        pass

    try:
        lines = open(f"requirements_{machine}.txt", "r").readlines()
        _add_to_dependencies(dependencies, lines)
    except IOError:
        pass

    try:
        lines = open(f"requirements_{platform_string}.txt", "r").readlines()
        _add_to_dependencies(dependencies, lines)
    except IOError:
        pass

    deps = list(dependencies.keys())
    deps.sort()

    print("\nUsing dependencies:")
    d = []
    for dep in deps:
        d.append(dependencies[dep])
        print(f"{dependencies[dep]}")

    dependencies = d

    global conda
    global mamba

    if mamba is None:
        # install mamba first!
        conda_install(["mamba"])
        mamba = find_mamba()

    conda_install(dependencies)


def build():
    print("\nCompiling the C++ code")


def install():
    print("\nInstalling sire")


if __name__ == "__main__":
    args = parse_args()

    if len(args.action) != 1:
        print("Please use either 'install_requires', 'build' or 'install'")
        sys.exit(-1)

    action = args.action[0]

    if action == "install":
        if not (args.skip_deps or args.skip_build):
            install_requires()

        if not args.skip_build:
            build()

        install()

    elif action == "build":
        if not args.skip_deps:
            install_requires()

        build()

    elif action == "install_requires":
        install_requires()
