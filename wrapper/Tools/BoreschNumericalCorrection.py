# Script to calculate the correction for releasing Boresch
# restraints with numerical integration- see Equation 12 in
# J. Phys. Chem. B 2003, 107, 35, 9535–9551. This is more robust
# than calculating the correction analytically, which requires the
# assumption of equilibrium values and is subject to instabilities
# under certain regimes.
# @author: Finlay Clark
# Overall layout based very closely on StandardState.py by Stefano Bosisio
# and Julien Michel, and calculation copied from restraints.py as in
# Yank https://github.com/choderalab/yank

import numpy as np
import os
from math import pi, sin, log
import scipy.integrate
from Sire import Units

from Sire.Tools import Parameter, resolveParameters
from Sire.Tools.OpenMMMD import *

# Constants
v0 = ((Units.meter3/1000)/Units.mole.value()).value() # A^3, the standard state volume
R = Units.gasr # kcal mol-1, the molar gas constant

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


    def numerical_distance_integrand(r, r0, kr):
        """Integrand for harmonic distance restraint. Domain is on [0, infinity],
        but this will be truncated to [0, 8 RT] for practicality.

        Args:
            r (float): Distance to be integrated, in Angstrom
            r0 (float): Equilibrium distance, in Angstrom
            kr (float): Force constant, in kcal mol-1 A-2

        Returns:
            float: Value of integrand
        """
        return (r**2)*np.exp(-(kr*(r-r0)**2)/(2*R*T))


    def numerical_angle_integrand(theta, theta0, ktheta):
        """Integrand for harmonic angle restraints. Domain is on [0,pi].

        Args:
            theta (float): Angle to be integrated, in radians
            theta0 (float): Equilibrium angle, in radians
            ktheta (float): Force constant, in kcal mol-1 rad-2

        Returns:
            float: Value of integrand
        """
        return sin(theta)*np.exp(-(ktheta*(theta-theta0)**2)/(2*R*T))


    def numerical_dihedral_integrand(phi, phi0, kphi):
        """Integrand for the harmonic dihedral restraints. Domain is on [-pi,pi].

        Args:
            phi (float): Angle to be integrated, in radians
            phi0 (float): Equilibrium angle, in radians
            kphi (float): Force constant, in kcal mol-1 rad-2

        Returns:
            float: Value of integrand
        """
        d_phi = abs(phi - phi0)
        d_phi_corrected = min(d_phi, 2*pi-d_phi) # correct for periodic boundaries
        return np.exp(-(kphi*d_phi_corrected**2)/(2*R*T))


    # Get Boresch restraint dict in dict form
    boresch_dict = dict(boresch_restraints_dict.val)
    T = temperature.val.value() # K

    # Force constants defined as E = k*x**2, so need to multiply all force constants
    # by 2 to correct for original definition (E= 0.5*k*x**2)

    # Radial
    r0 = boresch_dict['equilibrium_values']['r0'] # A
    kr = boresch_dict['force_constants']['kr']*2 # kcal mol-1 A-1
    dist_at_8RT = 4*np.sqrt((R*T)/kr) # Dist. which gives restraint energy = 8 RT
    r_min = max(0, r0-dist_at_8RT)
    r_max = r0 + dist_at_8RT
    integrand = lambda r: numerical_distance_integrand(r, r0, kr)
    z_r = scipy.integrate.quad(integrand, r_min, r_max)[0]

    # Angular
    for angle in ["thetaA", "thetaB"]:
        theta0 = boresch_dict["equilibrium_values"][f"{angle}0"] # rad
        ktheta = boresch_dict["force_constants"][f"k{angle}"]*2 # kcal mol-1 rad-2
        integrand = lambda theta: numerical_angle_integrand(theta, theta0, ktheta)
        z_r *= scipy.integrate.quad(integrand, 0, pi)[0]

    # Dihedral
    for dihedral in ["phiA", "phiB", "phiC"]:
        phi0 = boresch_dict["equilibrium_values"][f"{dihedral}0"] # rad
        kphi = boresch_dict["force_constants"][f"k{dihedral}"]*2 # kcal mol-1 rad-2
        integrand = lambda phi: numerical_dihedral_integrand(phi, phi0, kphi)
        z_r *= scipy.integrate.quad(integrand, -pi, pi)[0]

    dg = -R*T*log(8*pi**2*v0/z_r)

    print(f"Numerical correction for releasing Boresch restraints = {dg:.2f} kcal mol-1")
