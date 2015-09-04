
import os
import sys
import shutil

# This script copies the Sire files from its pkg
# directory into the appropriate anaconda directory.
# The "copy" is really hard linking, to save space,
# and the files are only copied if they don't already
# exist (nothing is overwritten)

try:
    # The directory to copy from
    fromdir = sys.argv[1]

    # The directroy to copy to
    todir = sys.argv[2]
except:
    print("USAGE: copy_into_conda.py from_dir to_dir")
    sys.exit(-1)

try:
    force = (sys.argv[3] == "force")
except:
    force = False

if force:
    print("force copying files from %s to %s" % (fromdir,todir))
else:
    print("copying files from %s to %s" % (fromdir,todir))

for file in os.listdir(fromdir):
    fromfile = "%s/%s" % (fromdir,file)
    tofile = "%s/%s" % (todir,file)

    # ignore all subdirectories
    if os.path.isdir(fromfile):
        continue

    if os.path.exists(tofile):
        if force:
            os.remove(tofile)
        else:
            continue

    # now copy symbolic links directly
    if os.path.islink(fromfile):
        shutil.copy(fromfile, tofile, follow_symlinks=False)
    else:
        os.link(fromfile, tofile)

