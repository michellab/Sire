import os

for file in os.popen("git ls-files --deleted", "r").readlines():
    cmd = "git checkout %s" % file.lstrip().rstrip()
    print(cmd)
    os.system(cmd)
