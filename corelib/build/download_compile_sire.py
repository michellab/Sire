
import os
import sys
import shutil

for arg in sys.argv[1:]:
    if arg == "-h" or arg == "--help":
        print("python download_compile_sire.py OPTIONS")
        print("\nScript to download and (optionally) compile and install Sire.")
        print("\nOptions:")
        print("    -r / --rebuild    Rebuild Sire from scratch every time you run this script.")
        print("    --download-only   Download Sire only. Don't compile or install.")
        print("    -b / --branch     Choose which branch of the source to download and install.")
        print("                      By default the 'trunk' (most up-to-date, stable version)")
        print("                      will be downloaded.")
        print("    -d / --directory  Supply the directory into which the source will be downloaded")
        print("                      and Sire will be compiled and installed. By default, the")
        print("                      current directory will be used.")
        print("    --no-execute      Only show the commands that will be run. Don't actually run anything.")
        print("\nSire is released under the GPL. For more information see http://siremol.org")
        sys.exit(0)


must_exit = False

sys.stdout.write("\n")

download_only = False
rebuild = False
no_execute = False
rundir = os.getcwd()
branch = "trunk"

for i in range(1,len(sys.argv)):
    arg = sys.argv[i]

    if arg == "--download-only":
        print("\nSire will only be downloaded. It will not be compiled or installed.")
        download_only = True

    elif arg == "-r" or arg == "--rebuild":
        print("\nSire will be rebuilt from scratch after download. This will be quite slow...")
        rebuild = True

    elif arg == "--no-execute":
        print("\nThis script will only print the commands to be run. It won't actually run anything...")
        no_execute = True

    elif arg == "-d" or arg == "--directory":
        print("\nDownloading, compiling and installing Sire in directory %s" % sys.argv[i+1])
        rundir = sys.argv[i+1]

    elif arg == "-b" or arg == "--branch":
        print("\nUsing the %s branch of Sire" % sys.argv[i+1])
        branch = sys.argv[i+1]

if must_exit:
    sys.exit(0)

if branch == "trunk":
    corelib = "https://sire.googlecode.com/svn/trunk/corelib"
    python = "https://sire.googlecode.com/svn/trunk/python"
else:
    corelib = "https://sire.googlecode.com/svn/corelib/%s" % branch
    python = "https://sire.googlecode.com/svn/python/%s" % branch

source_dir = "%s/Sire" % rundir
sire_app = "%s/sire.app" % rundir

if no_execute:
    print("mkdir %s" % source_dir)
    print("mkdir %s" % sire_app)
else:
    if not os.path.exists(source_dir):
        os.mkdir(source_dir)

    if not os.path.exists(sire_app):
        os.mkdir(sire_app)

print("Changing into directory %s" % source_dir)

print("cd %s" % source_dir)

if not no_execute:
    os.chdir(source_dir)

if not os.path.exists("corelib"):
    print("Downloading corelib from %s" % corelib)

    if no_execute:
        print("svn co %s ./corelib" % corelib)
    elif os.system("svn co %s ./corelib" % corelib) != 0:
        print("Failed to checkout the source for corelib")
        sys.exit(-1)
else:
    print("Updating corelib...")

    if no_execute:
        print("svn update ./corelib")
    elif os.system("svn update ./corelib") != 0:
        print("Failed to update the source for corelib")
        sys.exit(-1)

if not os.path.exists("python"):
    print("Downloading python wrappers from %s" % python)

    if no_execute:
        print("svn co %s ./python" % python)
    elif os.system("svn co %s ./python" % python) != 0:
        if os.system("svn co %s ./python" % python.replace("python","python2")) != 0:
            print("Failed to checkout the source for the python wrappers")
            sys.exit(-1)
else:
    print("Updating the python wrappers...")
    
    if no_execute:
        print("svn update ./python")
    elif os.system("svn update ./python") != 0:
        print("Failed to update the source for the python wrappers")
        sys.exit(-1)

if download_only:
    print("Everything downloaded or updated...")
    sys.exit(0)

if rebuild:
    if os.path.exists("build"):
        print("Rebuilding from scratch so removing existing build directory...")
        
        if no_execute:
            print("rm -rf ./build")
        else:
            shutil.rmtree("./build")

if not os.path.exists("build"):
    if no_execute:
        print("mkdir build")
    else:
        os.mkdir("build")

if not os.path.exists("build/corelib"):
    if no_execute:
        print("mkdir build/corelib")
    else:
        os.mkdir("build/corelib")

if not os.path.exists("build/python"):
    if no_execute:
        print("mkdir build/python")
    else:
        os.mkdir("build/python")

print("Compiling Sire to install in %s" % sire_app)

print("Configuring corelib...")

print("cd build/corelib")

if not no_execute:
    os.chdir("build/corelib")

if no_execute:
    print("nice cmake ../../corelib -DSIRE_INSTALL_PREFIX=%s" % sire_app)
elif os.system("nice cmake ../../corelib -DSIRE_INSTALL_PREFIX=%s" % sire_app) != 0:
    print("Could not successfully run cmake on corelib")
    sys.exit(-1)

if no_execute:
    print("nice make -j 4")
elif os.system("nice make -j 4") != 0:
    print("Could not successfully compile corelib")
    sys.exit(-1)

if no_execute:
    print("nice make -j 4 install/strip")
elif os.system("nice make -j 4 install/strip") != 0:
    print("Could not successfully install corelib")
    sys.exit(-1)

print("Configuring the python wrappers...")
print("cd ../python")

if not no_execute:
    os.chdir("../python")

if no_execute:
    print("nice cmake ../../python -DSIRE_APP=%s" % sire_app)
elif os.system("nice cmake ../../python -DSIRE_APP=%s" % sire_app) != 0:
    print("Could not successfully run cmake on python")
    sys.exit(-1)

if no_execute:
    print("nice make -j 4")
elif os.system("nice make -j 4") != 0:
    print("Could not successfully compile the python wrappers")
    sys.exit(-1)

if no_execute:
    print("nice make -j 4 install/strip")
elif os.system("nice make -j 4 install/strip") != 0:
    print("Could not successfully install the python wrappers")
    sys.exit(-1)
