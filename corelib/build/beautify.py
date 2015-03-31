
import os
import sys

for file in sys.argv[1:]:
    os.system("astyle --style=ansi --indent=spaces=4 --max-code-length=100 %s" % file)

