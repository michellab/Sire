
import os
import sys

from Sire import try_import

import Sire.Config

nose = try_import("nose")

import Sire.Base

branch = Sire.Base.getRepositoryBranch()
is_clean = Sire.Base.getRepositoryVersionIsClean()

testdir_exists = False

testdir = "%s/test" % Sire.Base.getShareDir()

if not os.path.exists(testdir):
    os.makedirs(testdir)

unittestdir = "%s/SireUnitTests/unittests" % testdir

old_cwd = os.getcwd()

def downloadTestsFromWebsite():
    """Function downloads the test suite from the website"""
    try:
        pycurl = try_import("pycurl", package_registry={"pycurl":"PyCurl"})
        test_package = "unittests_%s.tar.bz2" % Sire.Base.getReleaseVersion()
        test_package_file = "%s/unittests.tar.bz2" % testdir
        with open(test_package_file, "wb") as f:
            c = pycurl.Curl()
            c.setopt(c.URL, 
              "http://siremol.org/largefiles/sire_releases/download.php?name=%s" \
                      % test_package )
            c.setopt(c.WRITEDATA, f)
            c.perform()
            c.close()

            os.chdir(testdir)
            os.system("tar -jxvf unittests.tar.bz2")

            if not os.path.exists("SireUnitTests/README.md"):
                print("Could not find SireUnitTests/README.md in download - everything ok?")
                raise IOError()
            
            os.chdir(old_cwd)
            return True
    except:
        os.chdir(old_cwd)
        return False

if os.path.exists(unittestdir):
    testdir_exists

    try:
        gitexe = Sire.Base.findExe("git").absoluteFilePath()
        os.chdir("%s/.." % unittestdir)
        os.system("%s fetch" % gitexe)
        os.system("%s checkout %s" % (gitexe,"devel"))
        os.system("%s pull" % gitexe)
        os.chdir(old_cwd)
    except:
        pass

else:
    #if is_clean and branch == "master":
    #    testdir_exists = downloadTestsFromWebsite()
        
    if not testdir_exists:
        #things will be cloned
        try:
            gitexe = Sire.Base.findExe("git").absoluteFilePath()
            os.chdir(testdir)
            gitcmd = "%s clone https://github.com/michellab/SireUnitTests.git -b %s" % (gitexe,"devel")
            print("Cloning unittests from git repository - %s" % gitcmd)
            os.system(gitcmd)
            os.chdir(old_cwd)
        except:
            testdir_exists = downloadTestsFromWebsite()

print("You should find the unit tests in %s" % unittestdir)

print("\nRunning C++ unit tests...\n")
try:
    Sire.Base.UnitTest.runAll()
except:
    pass

print("\nNow running Python-based unit tests...\n")

testdir = unittestdir

if len(sys.argv) > 1:
    testdirs = sys.argv[1:]
    sys.argv = [sys.argv[0]]

    drop = []

    for dir in testdirs:
        if dir.startswith("-"):
            drop.append(dir)

    if len(drop) > 0:
        testdirs = []

        dirs = os.listdir(testdir)

        for dir in dirs:
            for d in drop:
                if d[1:] != dir:
                    testdirs.append(dir)

else:
    testdirs = os.listdir(testdir)

failures = []

for file in testdirs:
    subdir = os.path.join(testdir, file)

    if not file.startswith(".") and os.path.isdir(subdir):
        print("\n\nRunning tests in directory %s..." % subdir)
        print("############################################")
        os.chdir(subdir)
        success = nose.run()

        if not success:
            failures.append(file)

        os.chdir(old_cwd)

if len(failures) > 0:
    print("\n\nWARNING: SOME OF THE TEST JOBS FAILED!!!")
    print("#########################################")

    for failure in failures:
        print("One of more jobs in %s failed!" % failure)
else:
    print("\n\nHOORAY - ALL OF THE UNIT TESTS PASSED!!!")
    print("\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/")
