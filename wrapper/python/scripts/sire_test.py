
import os
import sys

from Sire import try_import

import Sire.Config

try_import("pycurl", package_registry={"pycurl":"PyCurl"})
try_import("nose")

sys.exit(-1)

testdir = Sire.Config.test_directory

old_cwd = os.getcwd()

if len(sys.argv) > 1:
    testdirs = sys.argv[1:]
    sys.argv = [sys.argv[0]]
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
