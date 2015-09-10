
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

if not os.path.exists(todir):
    os.makedirs(todir)

for file in os.listdir(fromdir):
    fromfile = "%s/%s" % (fromdir,file)
    tofile = "%s/%s" % (todir,file)

    # ignore all subdirectories
    if os.path.isdir(fromfile):
        # recurse into this directory
        if not os.path.exists(tofile):
            os.makedirs(tofile)
        elif not os.path.isdir(tofile):
            # cannot copy a directory into a file
            continue

        if force:
            os.system("%s %s %s/%s %s/%s force" % \
               (sys.executable,sys.argv[0],sys.argv[1],file,sys.argv[2],file))
        else:
            os.system("%s %s %s/%s %s/%s" % \
               (sys.executable,sys.argv[0],sys.argv[1],file,sys.argv[2],file))

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

