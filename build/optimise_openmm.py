import simtk.openmm as mm
import subprocess
import re
import os
import sys

# Check if the Cuda platform is recognised in platforms
platforms = [ mm.Platform.getPlatform(index).getName() for index in range(mm.Platform.getNumPlatforms()) ]
print (platforms)

if 'CUDA' in platforms:
    print('The CUDA platform is successfully recognised by openMM! Nothing needs to be done and you can use SOMD')
    sys.exit(0)
else:
    # Find OpenMM version
    mm_version = mm.__version__
    print(mm_version)

    # Let's find out which version of conda you are using!
    conda_base = os.path.abspath(os.path.dirname(sys.executable))
    conda_exe = None
    if os.path.exists("%s/conda" % conda_base):
        conda_exe = "%s/conda" % conda_base
    else:
        print("Cannot find a 'conda' binary in directory '%s'. "
              "Are you running this script using the python executable "
              "from a valid miniconda or anaconda installation within Sire?" % conda_base)
        sys.exit(-1)

    # does nvcc exists in the system?
    nvcc_out = None
    nvcc_release = None
    nvcc = subprocess.run(["nvcc --version"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
if nvcc.returncode == 0:
    nvcc_out = nvcc.stdout.decode("ascii").strip()
    match = re.search('release [+-]?([0-9]*[.]?[0-9])+', nvcc_out)
    print(match)
    nvcc_release = str(match.group().split(' ')[-1])
    print(type(nvcc_release))
    print('Found a CUDA toolkit release version: '+nvcc_release)
else:
    nvcc_out = nvcc.stderr.decode("ascii").strip()  
    print('Are you sure CUDA toolkit is installed? nvcc does not seem to exist!')
    print('If you want to use SOMD with CUDA plese double check your CUDA installation')
    sys.exit(-1)

if 'CUDA' in platforms:
    print('The CUDA platform is successfully recognised by openMM! You can now start using SOMD')
    sys.exit(0)
else:
    print('CUDA platform is not recognised by OpenMM!')
    print('Trying to update OpenMM to match your CUDA version %s for your OpenMM version %s' %(nvcc_release,mm_version))
    version = nvcc_release.replace('.','')
    conda_cmd = conda_exe +'  install --yes --no-deps -c omnia/label/cuda%s openmm=%s' %(version,mm_version)
    conda_proc = subprocess.run(conda_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if conda_proc.returncode == 0:
        platforms = [ mm.Platform.getPlatform(index).getName() for index in range(mm.Platform.getNumPlatforms()) ]
        print(platforms)
        if 'CUDA' in platforms:
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

