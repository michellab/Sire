
from Sire.Base import *

import shutil
import sys
import time
import os
import tempfile
import glob

from os.path import getsize, join

siredir = os.path.abspath( getInstallDir() )

print("\n#################################################################")
print(  "# This script will create a sire.run install file for the copy")
print(  "# of Sire running in directory %s" % siredir)
print(  "#################################################################\n")

sys.stdout.write("Please provide a name for the installer (sire.run) : ",)
sys.stdout.flush()

sire_run = sys.stdin.readline()

try:
    sire_run = sire_run.lstrip().rstrip()
except:
    sire_run = ''

if len(sire_run) == 0:
    sire_run = "sire.run"

print( "\nCreating Sire installer \"%s\"..." % (sire_run) )

python_exe = join(getInstallDir(), "bin", "python")
makeself =  join(getShareDir(), "build", "makeself.sh")
install_sire = join(getShareDir(), "build", "install_sire.sh")
remove_path = join(getShareDir(), "build", "remove_path.py")

# create a directory that will contain the files to be packages
with tempfile.TemporaryDirectory() as tempdir:
    tmp_sire = join(tempdir,"tmp_sire.app")

    print("Copying %s to %s" % (getInstallDir(),tmp_sire))
    shutil.copytree(getInstallDir(), "%s" % tmp_sire, symlinks=True)

    # Remove all of the .pyc and .pyo files as these can be regenerated and
    # take up too much space
    py_saved_space = 0
    print("Removing cached .pyc and .pyo files...")
    for root, dirs, files in os.walk(tmp_sire):
        for file in files:
            if file.endswith(".pyc") or file.endswith(".pyo"):
                filename = join(root,file)
                py_saved_space += getsize(filename)
                os.remove(filename)

    print("...files removed. Saved %d MB of space" % (py_saved_space/(1024*1024)))

    # Now use the share/Sire/remove_root.py script to replace all
    # absolute paths with {[{ROOT}]}
    print("Removing absolute directory references... (may take a while)...")
    os.system("%s %s %s %s" % (python_exe, remove_path, tmp_sire, siredir))

    print("Using 'makeself' to create the self-extracting installer...")
    share_dir = Sire.Base.getShareDir().replace(Sire.Base.getInstallDir(),"./tmp_sire.app")
    os.system("%s --current %s %s \"Sire Molecular Simulation Framework\" %s/build/install_sire.sh" \
                   % (makeself, tempdir, sire_run, share_dir))

print( "\nAll done :-). Just type %s to install Sire." % sire_run )
