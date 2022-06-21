
import sys
import os
import subprocess

script = os.path.abspath(sys.argv[0])

# we want to import the 'get_requirements' package from this directory
sys.path.insert(0, os.path.dirname(script))

from parse_requirements import parse_requirements

# go up one directories to get the source directory
# (this script is in Sire/actions/)
srcdir = os.path.dirname(os.path.dirname(script))

print(f"Sire source is in {srcdir}")

condadir = os.path.join(srcdir, "recipes", "sire")

print(f"conda recipe in {condadir}")

# Store the name of the recipe and template YAML files.
recipe = os.path.join(condadir, "meta.yaml")
template = os.path.join(condadir, "template.yaml")

# Now parse all of the requirements
run_reqs = parse_requirements(os.path.join(srcdir, "requirements.txt"))
print(run_reqs)
build_reqs = parse_requirements(os.path.join(srcdir, "requirements_build.txt"))
print(build_reqs)
bss_reqs = parse_requirements(os.path.join(srcdir, "requirements_bss.txt"))
print(bss_reqs)


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


def dep_lines(deps):
    lines = []

    for dep in deps:
        lines.append(f"    - {dep}\n")

    return "".join(lines)


run_reqs = dep_lines(run_reqs)
build_reqs = dep_lines(build_reqs)
bss_reqs = dep_lines(bss_reqs)


with open(recipe, "w") as FILE:
    for line in lines:
        if line.find("SIRE_BUILD_REQUIREMENTS") != -1:
            line = build_reqs
        elif line.find("SIRE_RUN_REQUIREMENTS") != -1:
            line = run_reqs
        elif line.find("SIRE_BSS_REQUIREMENTS") != -1:
            line = bss_reqs
        else:
            line = line.replace("SIRE_VERSION", sire_version)
            line = line.replace("SIRE_BUILD", str(sire_build))
            line = line.replace("SIRE_BRANCH", sire_branch)

        FILE.write(line)
