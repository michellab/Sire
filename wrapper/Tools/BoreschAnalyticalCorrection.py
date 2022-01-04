# Script to calculate the analytical correction for releasing Boresch
# restraints - see Equation 32 in J. Phys. Chem. B 2003, 107, 35, 9535–9551
# This includes the standard state correction and will be reliable when
# restraints are sufficiently strong, r is sufficiently far from zero,
# and r2, r1, l1, and r1, l1, l2 are sufficiently far from collinear. If
# this is not the case, the semi-analytical expression should be used.
# @author: Finlay Clark
# Based very closely on StandardState.py by Stefano Bosisio and Julien Michel

import numpy as np
import os
from math import pi, sin, log

from Sire.Tools import Parameter, resolveParameters
#from Sire.Units import *
from Sire.Tools.OpenMMMD import *

# Constants
v0 = 1660.53907 # A^3, the standard state volume
R = 0.00198720425864083 # kcal mol-1, the molar gas constant

#temperature = Parameter("temperature", 25 * celsius, """Simulation temperature""")

@resolveParameters
def run():
    try:
        host = os.environ['HOSTNAME']
    except KeyError:
        host = "unknown"
    print("### Running Standard state correction calculation on %s ###" % host)

    if verbose.val:
        print("###====================Utilised Parameters=====================###")
        print(temperature)
        print(boresch_restraints_dict)
        print ("###===========================================================###\n")

    # Get Boresch restraint dict in dict form
    boresch_dict = dict(boresch_restraints_dict.val)

    # Params
    T = temperature.val.value() # K
    r0 = boresch_dict['equilibrium_values']['r0'] # A
    thetaA0 = boresch_dict["equilibrium_values"]["thetaA0"] # rad
    thetaB0 = boresch_dict["equilibrium_values"]["thetaB0"] # rad

    force_constants = list(boresch_dict["force_constants"].values()) # kcal mol-1 A-2 or rad-2
    prod_force_constants = np.prod(force_constants)

    # Calculation
    numerator = 8*(pi**2)*v0*np.sqrt(prod_force_constants)
    denominator = (r0**2)*sin(thetaA0)*sin(thetaB0)*(2*pi*R*T)**3

    dg = -R*T*log(numerator/denominator)

    print(f"Analytical correction for releasing Boresch restraints = {dg:.2f} kcal mol-1")
    print("WARNING !!! The analytical correction is only reliable when restraints are "
          "sufficiently strong, r is sufficiently far from zero, and r2, r1, l1, and" 
          "r1, l1, l2 are sufficiently far from collinear")