
import os
import sys
import glob

script = os.path.abspath(sys.argv[0])

# go up one directories to get the source directory
# (this script is in Sire/actions/)
srcdir = os.path.dirname(os.path.dirname(script))

print(f"sire source is in {srcdir}\n")

# Get the anaconda token to authorise uploads
if "ANACONDA_TOKEN" in os.environ:
    conda_token = os.environ["ANACONDA_TOKEN"]
else:
    conda_token = "TEST"

# get the root conda directory
conda = os.environ["CONDA"]

# Set the path to the conda-bld directory.
conda_bld = os.path.join(conda, "envs", "sire_build", "conda-bld")

print(f"conda_bld = {conda_bld}")

# Find the packages to upload
sire_pkg = glob.glob(os.path.join(conda_bld, "*-*", "sire-*.tar.bz2"))
kcombu_pkg = glob.glob(os.path.join(conda_bld, "*-*", "kcombu-*.tar.bz2"))

if len(sire_pkg) == 0:
    print("No sire packages to upload?")
    sys.exit(-1)

packages = sire_pkg

if len(kcombu_pkg) == 0:
    packages = packages + kcombu_pkg

print(f"Uploading packages:")
print(" * ", "\n *  ".join(packages))

packages = " ".join(packages)


def run_cmd(cmd):
    import subprocess
    p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    return str(p.stdout.read().decode("utf-8")).lstrip().rstrip()

gitdir = os.path.join(srcdir, ".git")

tag = run_cmd(f"git --git-dir={gitdir} --work-tree={srcdir} tag --contains")

# If the tag is not empty, then set the label to main (this is a release)
if tag is not None and tag.lstrip().rstrip() != "":
    print(f"\nTag {tag} is set. This is a 'main' release.")
    label = "--label main --label dev"
else:
    # this is a development release
    print("\nNo tag is set. This is a 'devel' release.")
    label = "--label dev"

# Upload the packages to the michellab channel on Anaconda Cloud.
cmd = f"anaconda --token {conda_token} upload --user michellab {label} --force {packages}"

print(f"\nUpload command:\n\n{cmd}\n")

# Label release packages with main and dev so that dev is at least as new as
# main. Only need to uncomment the libcpuid and kcombu package uploads when
# there new versions are released.
if conda_token == "TEST":
    print("Not uploading as the ANACONDA_TOKEN is not set!")
    sys.exit(-1)

output = run_cmd(cmd)

print(output)

print("Package uploaded!")
