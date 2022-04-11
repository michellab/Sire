##
## Simple script that reinstalls all of the __init__.py
## files in Sire, without having to go through a full
## compile. Use it from this directory, via
## {SIRE_APP}/bin/python quick_wrapper.py
## This will copy all of the files, avoiding
## the need for a slower install.
##
## Only use this if you know what you are doing ;-)
##

import Sire

import glob
import os
import shutil

inits = glob.glob("../wrapper/*/__init__.py")

inits.append("../wrapper/__init__.py")
path = os.path.dirname(Sire.__file__)

for init in inits:
    dest = f"{path}/{init.replace('../wrapper/', '')}"

    if os.path.exists(dest):
        cmd = f"cp {init} {dest}"
        print(cmd)
        shutil.copy(init, dest)
