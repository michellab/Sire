import simtk.openmm as mm
import subprocess
import re
import os
import sys
import importlib

import Sire.Base

# Check if the Cuda platform is recognised in platforms
platforms = [ mm.Platform.getPlatform(index).getName() for index in range(mm.Platform.getNumPlatforms()) ]

if 'CUDA' in platforms:
    print('The CUDA platform is successfully recognised by openMM! Nothing needs to be done and you can use SOMD')
    print(platforms)
    sys.exit(0)
else:
    print('CUDA platform is not recognised by OpenMM!')
    print('available platforms are: ')
    print(platforms)
    print("Let's see if we can do something about this....\n")

    # Find OpenMM version
    mm_version = mm.__version__

    # Figure out where the Conda executable is.

    # 1) First assume that this is a Sire Mininconda. The executable will be located
    #    at something like $HOME/sire.app/bin/conda
    conda_bin = Sire.Base.getBinDir()
    conda_exe = None
    if os.path.exists("%s/conda" % conda_bin):
        conda_exe = "%s/conda" % conda_bin

    # 2) If no excutable was found, then search the system path for the conda binary.
    #    If this is a Conda install, and the environment is activated, then the
    #    correct binary will be in the path.
    else:
        try:
            conda_exe = Sire.Base.findExe("conda").absoluteFilePath()
        except:
            raise Exception("Cannot find a 'conda' binary in the system. Are you "
                            "running this script using the Python interpreter from "
                            "a valid Sire application or Conda installation of Sire?") from None

    # does nvcc exists in the system?
    nvcc_out = None
    nvcc_release = None
    nvcc = subprocess.run(["nvcc --version"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if nvcc.returncode == 0:
        nvcc_out = nvcc.stdout.decode("ascii").strip()
        match = re.search('release [+-]?([0-9]*[.]?[0-9])+', nvcc_out)
        nvcc_release = str(match.group().split(' ')[-1])
        print('Found a CUDA toolkit release version: '+nvcc_release)
    else:
        nvcc_out = nvcc.stderr.decode("ascii").strip()
        print('Are you sure CUDA toolkit is installed? nvcc does not seem to exist!')
        print('If you want to use SOMD with CUDA plese double check your CUDA installation')
        sys.exit(-1)

    # Found CUDA version and openMM version, now trying to update this with conda!
    print('Trying to update OpenMM to match your CUDA version %s for your OpenMM version %s' %(nvcc_release,mm_version))
    print('This may take a little while. Please hold tight!')
    print('................................................')
    version = nvcc_release.replace('.','')
    conda_cmd = conda_exe +'  install --yes --no-deps -c omnia/label/cuda%s openmm=%s' %(version,mm_version)
    conda_proc = subprocess.run(conda_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    importlib.reload(mm)
    if conda_proc.returncode == 0:
        platform = [ mm.Platform.getPlatform(index).getName() for index in range(mm.Platform.getNumPlatforms()) ]
        print(platform)
        if 'CUDA' in platform:
            print('You have now successfully configured OpenMM for use with Sire')
            sys.exit(0)
        else:
            print("Something didn't work out with the update of OpenMM via conda, have a look at the output.")
            print(conda_proc.stdout.decode("ascii".strip()))
            print(conda_proc.stderr.decode("ascii".strip()))
    else:
        print('Something went wrong with the conda installation. Please take a look at the ouptut')
        print('Executing command %s failed' %conda_cmd)
        print(conda_proc.stderr.decode("ascii".strip()))

