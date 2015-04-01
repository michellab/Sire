
from Sire.Config import *
from Sire.Base import *

import os
import sys

siredir = getInstallDir()

print("\n#################################################################")
print(  "# This script will download the source code for the copy")
print(  "# of Sire running in directory %s" % siredir)
print(  "# (note that you must have subversion installed for this to work)")
print(  "#################################################################\n")

print("\n...looking for subversion executable...")
svn = findExe("svn").absoluteFilePath()
print("...found! (%s)" % svn)


print("\n...checking the version of the source code that was used to compile this binary...")
coreversion = sire_repository_version
pythonversion = sire_python_repository_version

try:
    coreversion = int(coreversion)
    pythonversion = int(pythonversion)
except:
    print("Either the version of the corelib (%s) or python wrappers (%s)" % (coreversion,pythonversion))
    print("is not integer, meaning that the version of the source used to ")
    print("compile this binary of Sire cannot be determined.")
    print("Unfortunately, this means the source cannot be downloaded.")
    sys.exit(-1)

print("...corelib version %d, python wrapper version %d.\n" % \
              (coreversion,pythonversion))

sys.stdout.write("Please provide the directory in which to download the source (current directory) : ",)
sys.stdout.flush()

prefix = sys.stdin.readline()

try:
    prefix = prefix.lstrip().rstrip()

    if len(prefix) == 0:
        prefix = '.'
except:
    prefix = '.'

srcdir = "%s/SireSource" % prefix

print("\nDownloading the Sire source code into %s" % srcdir)

try:
    os.makedirs(srcdir)
except:
    print("Cannot create the directory %s." % srcdir)
    print("Please make sure that this path is accessible, ")
    print("and that the directory doesn't already exist.")
    sys.exit(-1)

cmd = "%s co -r %s %s %s/corelib" % (svn,coreversion,sire_repository_url,srcdir)
print("\nDownloading corelib source using command\n%s" % cmd)

if os.system(cmd) != 0:
    print("WARNING: Something went wrong with the download!")
else:
    print("\n...COMPLETE")

cmd = "%s co -r %s %s %s/python" % (svn,pythonversion,sire_python_repository_url,srcdir)
print("\nDownloading python wrapper source using command\n%s" % cmd)

if os.system(cmd) != 0:
    print("WARNING: Something went wrong with the download!")
else:
    print("\n...COMPLETE")

print("\n...creating build directories and placing README, INSTALL and COPYING files...")
try:
    os.makedirs("%s/build" % srcdir)
    os.makedirs("%s/build/corelib" % srcdir)
    os.makedirs("%s/build/python" % srcdir)
except:
    pass

try:
    shutil.copyfile("%s/corelib/build/INSTALL" % srcdir, "%s/INSTALL" % srcdir)
    shutil.copyfile("%s/corelib/build/README" % srcdir, "%s/README" % srcdir)
    shutil.copyfile("%s/corelib/build/COPYING" % srcdir, "%s/COPYING" % srcdir)
except:
    pass

print("\nFinished :-)")

