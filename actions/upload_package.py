
import os
import sys
import glob

script = os.path.abspath(sys.argv[0])

# go up one directories to get the source directory
# (this script is in Sire/actions/)
srcdir = os.path.dirname(os.path.dirname(script))

print(f"Sire source is in {srcdir}")

# Get the anaconda token to authorise uploads
conda_token = os.environ["ANACONDA_TOKEN"]

# get the root conda directory
conda = os.environ["CONDA"]

# Set the path to the conda-bld directory.
conda_bld = os.path.join(conda, "envs", "sire_build", "conda-bld")

# Find the packages to upload
sire_pkg = glob.glob(os.path.join(conda_bld, "*", "sire-*.tar.bz2"))
fkcombu_pkg = glob.glob(os.path.join(conda_bld, "*", "fkcombu-*.tar.bz2"))

print(sire_pkg)
print(fkcombu_pkg)

def run_cmd(cmd):
    import subprocess
    p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    return str(p.stdout.read().decode("utf-8")).lstrip().rstrip()

gitdir = os.path.join(srcdir, ".git")

tag = run_cmd(f"git --git-dir={gitdir} --work-tree={stddir} tag --contains")

print(tag)

# If the tag is not empty, then set the label to main (this is a release)
if tag is not None:
    label = "main"
else:
    # this is a development release
    label = "devel"

print(label)


