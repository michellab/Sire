from re import I


import sys
import os
import subprocess

script = os.path.abspath(sys.argv[0])

# go up one directories to get the source directory
# (this script is in Sire/actions/)
srcdir = os.path.dirname(os.path.dirname(script))

print(f"Sire source is in {srcdir}")

condadir = os.path.join(srcdir, "recipes", "sire")

print(f"conda recipe in {condadir}")

# Store the name of the recipe and template YAML files.
recipe = os.path.join(condadir, "meta.yaml")
template = os.path.join(condadir, "template.yaml")

def run_cmd(cmd):
    p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    return str(p.stdout.read().decode("utf-8")).lstrip().rstrip()

gitdir = os.path.join(srcdir, ".git")

# Get the Sire version. (Latest tag.)
sire_version = run_cmd(f"git --git-dir={gitdir} --work-tree={srcdir} describe --tags --abbrev=0")
print(sire_version)

# Get the build number. (Number of commits since last tag.)
sire_build = len(run_cmd(f"git --git-dir={gitdir} --work-tree={srcdir} log --oneline {sire_version}..").split("\n"))
print(sire_build)

# Get the Sire branch.
sire_branch = run_cmd(f"git --git-dir={gitdir} --work-tree={srcdir} rev-parse --abbrev-ref HEAD")
print(sire_branch)

lines = open(template, "r").readlines()

with open(recipe, "w") as FILE:
    for line in lines:
        line = line.replace("SIRE_VERSION", sire_version)
        line = line.replace("SIRE_BUILD", str(sire_build))
        line = line.replace("SIRE_BRANCH", sire_branch)

        FILE.write(line)
