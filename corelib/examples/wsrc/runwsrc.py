
# Import the WSRC module that provides all of the code
# to run a waterswap simulation
from Sire.Tools import WSRC

# Import the Units module so that we can specify the dimensions
# (units) of input parameters
from Sire.Units import *

# Import the Stream module so that we can save restart files
import Sire.Stream

# Create a dictionary of parameters for waterswap
params = {}

# Give the name of the input amber topology and coordinate files
# that contain the solvated, bound protein-ligand system
params["protein topfile"] = "proteinbox.top"
params["protein crdfile"] = "proteinbox.crd"

## Specify the name of one of the residues in the ligand to be swapped.
# The WSRC module will find the ligand molecule by looking for the first
# molecule loaded that contains a residue with this name
params["ligand name"] = "ose"

# Specify the radius of the reflection sphere
params["reflection radius"] = 10 * angstroms

# Specify the temperature of the simulation
params["temperature"] = 25 * celsius

# Specify the lambda values to use
params["lambda values"] = [ 0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.49,
                            0.51, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99 ]
params["lambda values"] = [ 0.01, 0.3, 0.7, 0.99 ]

# Specify the number of RETI iterations
params["nmoves"] = 10

# Specify the number of Monte Carlo moves to perform per RETI iteration
params["nsubmoves"] = 100

params["waterbox only"] = True

# Specify the number of equilibration moves to perform when setting up the system
params["nequilmoves"] = 50

# Switch to use the TIP4P water model
# params["water model"] = "tip4p"

# Specify the frequency of updating the residue- and water-based energy monitors
params["energy monitor frequency"] = 5

# I want to save the PDBs for debugging...
params["save pdb"] = True

params["pdb frequency"] = 5

params["restart	frequency"] = 6

# Actually run the simulation :-)
WSRC.run(params)

# Note that you can resume the simulation by rerunning this script. The script
# records its progress so will always resume from where it last finished. This
# means that you can increase the number of RETI iterations by increasing
# params["nsubmoves"] and rerunning the script.

# If you want the simulation to be restarted from the beginning, you need to
# remove all restart and output files, e.g. by running "rm -rf output *.s3"
