
## This script scans through the sire.app directory and
## finds all textual references to {[{ROOT}]} and replaces
## that with the the base directory that is supplied.
## This should allow absolute paths to be updated
## after a file move.

import os
import sys
import shutil

special = "{[{ROOT}]}"

root_dir = sys.argv[1]

print("Scanning through %s" % (root_dir))

def is_binary(filename):
    """
    Return true if the given filename appears to be binary.
    File is considered to be binary if it contains a NULL byte.
    FIXME: This approach incorrectly reports UTF-16 as binary.
    """
    with open(filename, 'rb') as f:
        for block in f:
            if '\0' in block:
                return True
    return False

def contains_special(filename):
    if filename.endswith("restore_path.py"):
        return False

    FILE = open(filename, "r")

    line = FILE.readline()

    while line:
        if line.find(special) != -1:
            FILE.close()
            return True

        line = FILE.readline()

    return False

def restore_root(filename, root):
    FILE = open(filename, "r")
    FILE2 = open("%s.tmpcopy" % filename, "w")

    line = FILE.readline()

    while line:
        if line.find(special) != -1:
            line = line.replace(special, root)

        FILE2.write(line)

        line = FILE.readline()

    FILE.close()
    FILE2.close()

    shutil.copystat(filename, "%s.tmpcopy" % filename)
    shutil.move("%s.tmpcopy" % filename, filename)    

def scanDir(root_dir, top_root_dir):
    for file in os.listdir(root_dir):
        fullfile = "%s/%s" % (root_dir,file)

        if os.path.isdir(fullfile):
            scanDir(fullfile, top_root_dir)
            continue

        if os.path.islink(fullfile):
            continue

        if is_binary(fullfile):
            continue

        # this is a text file and not a symbolic link
        # see if it contains the text '{[{ROOT}]}'
        if contains_special(fullfile):
            restore_root(fullfile, top_root_dir)

scanDir(root_dir, root_dir)
