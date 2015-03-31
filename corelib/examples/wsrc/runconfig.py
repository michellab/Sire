
# Import the WSRC module that provides all of the code
# to run a waterswap simulation
from Sire.Tools import WSRC

# Import the function used to read the config file
from Sire.Tools import readParams

import sys

params = readParams(sys.argv[1])

WSRC.run(params)
