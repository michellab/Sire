
from Sire.Base import *

import shutil
import sys
import os
import tempfile
import glob

from os.path import getsize, join

siredir = getInstallDir()

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

makeself =  join(getShareDir(), "build", "makeself.sh")
install_sire = join(getShareDir(), "build", "install_sire.sh")

# create a directory that will contain the files to be packages
with tempfile.TemporaryDirectory() as tempdir:
    print("Copying %s to %s" % (getInstallDir(),join(tempdir,"tmp_sire.app")))
    shutil.copytree(getInstallDir(), "%s" % (join(tempdir,"tmp_sire.app")), symlinks=True)
    #print("Copying %s to %s" % (install_sire, join(tempdir,"install_sire.sh")))
    #shutil.copyfile(install_sire, "%s" % join(tempdir,"install_sire.sh"))

    # Remove all of the .pyc and .pyo files as these can be regenerated and
    # take up too much space
    py_saved_space = 0
    print("Removing cached .pyc and .pyo files...")
    for root, dirs, files in os.walk(getBundledLibDir()):
        for file in files:
            if file.endswith(".pyc") or file.endswith(".pyo"):
                filename = join(root,file)
                py_saved_space += getsize(filename)
                os.remove(filename)
    print("...files removed. Saved %d MB of space" % (py_saved_space/(1024*1024)))

    # Extract the development files and package them up
    print("Using 'makeself' to create a self-extracting installer for the development files...")
    develdir = join(tempdir,"devel")
    os.makedirs(develdir)
    os.makedirs(join(develdir,"bundled"))
    os.makedirs(join(develdir,"bundled","lib"))

    print("Moving development files to %s" % develdir)
    shutil.move(join(tempdir,"tmp_sire.app","include"), develdir)
    shutil.move(join(tempdir,"tmp_sire.app","bundled","include"), join(develdir,"bundled"))
    shutil.move(join(tempdir,"tmp_sire.app","bundled","mkspecs"), join(develdir,"bundled"))
    shutil.move(join(tempdir,"tmp_sire.app","bundled","lib","pkgconfig"), join(develdir,"bundled","lib"))
    
    for lib in glob.glob(join(tempdir,"tmp_sire.app","bundled","lib","*debug*")):
        shutil.move(lib, join(develdir,"bundled","lib"))

    devel_saved_space = 0
    for root, dirs, files in os.walk(develdir):
        for file in files:
            devel_saved_space += getsize(join(root,file))

    os.system("%s --current %s %s \"Sire Development Files\" echo \"Unpacked files\"" \
                % (makeself, develdir, join(tempdir,"tmp_sire.app","sire_devel.run")))

    print("...package made. Saved %d MB of space" % (devel_saved_space/(1024*1024)))

    # remove the devel directory
    shutil.rmtree(develdir, ignore_errors=True)

    print("Using 'makeself' to create the self-extracting installer...")
    os.system("%s --current %s %s \"Sire Molecular Simulation Framework\" ./tmp_sire.app/share/Sire/build/install_sire.sh" \
                   % (makeself, tempdir, sire_run))

print( "\nAll done :-). Just type %s to install Sire." % sire_run )
