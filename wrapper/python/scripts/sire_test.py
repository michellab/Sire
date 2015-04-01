
import os
import sys

import Sire.Config

try:
    import nose
except:
    print("Cannot import nose. Must install nose via sire.app/bin/install_package install nose")
    sys.exit(-1)

testdir = Sire.Config.test_directory

old_cwd = os.getcwd()

for file in os.listdir(testdir):
    subdir = os.path.join(testdir, file)

    if not file.startswith(".") and os.path.isdir(subdir):
        print("Running tests in directory %s..." % subdir)
        os.chdir(subdir)
        nose.run()
        os.chdir(old_cwd)
