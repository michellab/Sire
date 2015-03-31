
import os
import sys
import shutil

rebuild = False

try:
    if sys.argv[1] == "--rebuild":
        rebuild = True
except:
    pass

home = os.getenv("HOME")
testdir = "%s/sire_test" % home

branch = "branches/devel"

corelib = "https://sire.googlecode.com/svn/corelib/%s" % branch
python = "https://sire.googlecode.com/svn/python/%s" % branch

print("Running tests in directory %s" % testdir)
print("Testing corelib from %s" % corelib)
print("Testing python from %s" % python)

if not os.path.exists(testdir):
    os.mkdir(testdir)

print("Changing into directory %s" % testdir)
os.chdir(testdir)

if not os.path.exists("corelib"):
    print("Downloading corelib from %s" % corelib)
    if os.system("svn co %s ./corelib" % corelib) != 0:
        print("Failed to checkout the source for corelib")
        sys.exit(-1)
else:
    print("Updating corelib...")
    if os.system("svn update ./corelib") != 0:
        print("Failed to update the source for corelib")
        sys.exit(-1)

if not os.path.exists("python"):
    print("Downloading python wrappers from %s" % python)
    if os.system("svn co %s ./python" % python) != 0:
        print("Failed to checkout the source for the python wrappers")
        sys.exit(-1)
else:
    print("Updating the python wrappers...")
    if os.system("svn update ./python") != 0:
        print("Failed to update the source for the python wrappers")
        sys.exit(-1)

if rebuild:
    if os.path.exists("build"):
        print("Removing existing build directory...")
        shutil.rmtree("./build")

if not os.path.exists("build"):
    os.mkdir("build")

if not os.path.exists("build/corelib"):
    os.mkdir("build/corelib")

if not os.path.exists("build/python"):
    os.mkdir("build/python")

sire_app = "%s/sire.app" % testdir
print("Compiling Sire to install in %s" % sire_app)

print("Configuring corelib...")
os.chdir("build/corelib")

if os.system("nice cmake ../../corelib -DSIRE_INSTALL_PREFIX=%s" % sire_app) != 0:
    print("Could not successfully run cmake on corelib")
    sys.exit(-1)

if os.system("nice make -j 4") != 0:
    print("Could not successfully compile corelib")
    sys.exit(-1)

if os.system("nice make -j 4 install/strip") != 0:
    print("Could not successfully install corelib")
    sys.exit(-1)

print("Configuring the python wrappers...")
os.chdir("../python")

if os.system("nice cmake ../../python -DSIRE_APP=%s" % sire_app) != 0:
    print("Could not successfully run cmake on python")
    sys.exit(-1)

if os.system("nice make -j 4") != 0:
    print("Could not successfully compile the python wrappers")
    sys.exit(-1)

if os.system("nice make -j 4 install/strip") != 0:
    print("Could not successfully install the python wrappers")
    sys.exit(-1)
