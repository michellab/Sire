
from Sire.IO import *
from Sire.Mol import *
from Sire.CAS import *
from Sire.System import *
from Sire.Move import *
from Sire.MM import *
from Sire.FF import *
from Sire.Units import *
from Sire.Vol import *
from Sire.Maths import *
from Sire.Base import *
from Sire.Qt import *

import os,re,sys
import shutil


combining_rules = "arithmetic"
temperature = 25 * celsius
pressure = 1 * atm
coulomb_cutoff = 1000 * angstrom
coulomb_feather = 999.5 * angstrom
lj_cutoff = 1000 * angstrom
lj_feather = 999.5 * angstrom
#############################################################

top_file = "test/io/SYSTEM.top"
crd_file = "test/io/SYSTEM.crd"
print("Solvated protein/ligand complex: Loading a top file...")
timer = QTime()
timer.start()

amber = Amber()
molecules, space = amber.readCrdTop(crd_file, top_file)
ms = timer.elapsed()
print("...took %d s" % (ms/1000.))
# Overload, we want to calc the energy in a non periodic box for comparison with Sander
space = Cartesian()

moleculeNumbers = molecules.molNums()
moleculeList = []

for moleculeNumber in moleculeNumbers:
    molecule = molecules.molecule(moleculeNumber).molecule()
    moleculeList.append(molecule)

system = System()

solute = MoleculeGroup("solute", moleculeList[0])
protein = MoleculeGroup("protein", moleculeList[1])
solvent = MoleculeGroup("solvent")
for molecule in moleculeList[2:]:
    solvent.add(molecule)

all = MoleculeGroup("all")
all.add(solute)
all.add(solvent)
all.add(protein)
# Add these groups to the System
system.add(solute)
system.add(solvent)
system.add(protein)
system.add(all)

# Now create all of the forcefields
# - first solvent-solvent coulomb/LJ (CLJ) energy
solventff = InterCLJFF("solvent:solvent")
solventff.add(solvent)

# Now solute bond, angle, dihedral energy
solute_intraff = InternalFF("solute_intraff")
solute_intraff.add(solute)

# Now solute intramolecular CLJ energy
solute_intraclj = IntraCLJFF("solute_intraclj")
solute_intraclj.add(solute)

# The protein bond, angle, dihedral energy
protein_intraff = InternalFF("protein_intraff")
protein_intraff.add(protein)

# The protein intramolecular CLJ energy
protein_intraclj = IntraCLJFF("protein_intraclj")
protein_intraclj.add(protein)

# Now the solute-solvent CLJ energy
solute_solventff = InterGroupCLJFF("solute:solvent")
solute_solventff.add(solute, MGIdx(0))
solute_solventff.add(solvent, MGIdx(1))

# Now the solute-protein CLJ energy
solute_proteinff = InterGroupCLJFF("solute:protein")
solute_proteinff.add(solute, MGIdx(0))
solute_proteinff.add(protein, MGIdx(1))
  
# The protein-solvent energy 
protein_solventff = InterGroupCLJFF("protein:solvent")
protein_solventff.add(protein, MGIdx(0))
protein_solventff.add(solvent, MGIdx(1))

# Here is the list of all forcefields
forcefields = [ solute_intraff, solute_intraclj,
                solventff, solute_solventff,
                protein_intraff, protein_intraclj,
                solute_proteinff, protein_solventff ] 
# Add these forcefields to the system
for forcefield in forcefields:
    system.add(forcefield)

system.setProperty( "space", space )
system.setProperty( "switchingFunction", 
                    HarmonicSwitchingFunction(coulomb_cutoff, coulomb_feather,
                                              lj_cutoff, lj_feather) )
system.setProperty( "combiningRules", VariantProperty(combining_rules) )

total_nrg = solute_intraclj.components().total() + solute_intraff.components().total() +\
    solventff.components().total() + solute_solventff.components().total() +\
    protein_intraclj.components().total() + protein_intraff.components().total() + \
    solute_proteinff.components().total() + protein_solventff.components().total() 

e_total = system.totalComponent()
system.setComponent( e_total, total_nrg )

# Add a space wrapper that wraps all molecules into the box centered at (0,0,0)
#system.add( SpaceWrapper(Vector(0,0,0), all) )

print("\nTotal energy ")
print(system.energy())

print("Components energies ")
for component in list(system.energyComponents().keys()):
    print(component, system.energyComponents().value(component) * kcal_per_mol)

# Note that tip3p water are likely to have bonds between hydrogen atoms.
PDB().write(all, "out.pdb")

print("The AMBER11/sander energies for this system are ") 
print("""
# NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
# Etot   =    -47010.2216  EKtot   =         0.0000  EPtot      =    -47010.2216
# BOND   =       898.1982  ANGLE   =      5310.2620  DIHED      =      2922.5644
# 1-4 NB =       790.8755  1-4 EEL =      7702.0145  VDWAALS    =      7345.0484
# EELEC  =    -71979.1846  EHBOND  =         0.0000  RESTRAINT  =         0.0000
# EKCMT  =         0.0000  VIRIAL  =         0.0000  VOLUME     =   1856243.3813
""")


print("The total energies differ by %12.8f kcal/mol" % ( system.energy().value() - -47010.2216 )) 
print("""Some difference is expected in the coulombic energies due to the 
use of different constants by Sander ( %12.8f ) and Sire ( %12.8f ). 
The other terms should be in excellent agreement. """ % (18.2223*18.2223, 
                                                         one_over_four_pi_eps0))  

#############################################################


#top_file = "/home/julien/projects/cyps/sire/setup/leap/LI8/LIG_SOLV.top"
#crd_file = "/home/julien/projects/cyps/sire/setup/leap/LI8/LIG_SOLV.crd"

#top_file = "/home/julien/projects/cyps/sire/setup/leap/LI8/LIG_GAS.top"
#crd_file = "/home/julien/projects/cyps/sire/setup/leap/LI8/LIG_GAS.crd"

#top_file = "/home/julien/projects/cyps/sire/setup/leap/methanol/LIG_GAS.top"
#crd_file = "/home/julien/projects/cyps/sire/setup/leap/methanol/LIG_GAS.crd"

#top_file = "/home/julien/projects/cyps/sire/setup/leap/ethane/LIG_GAS.top"
#crd_file = "/home/julien/projects/cyps/sire/setup/leap/ethane/LIG_GAS.crd"

#top_file = "/home/julien/projects/cyps/sire/setup/leap/ethane-edit/LIG_GAS.top"
#crd_file = "/home/julien/projects/cyps/sire/setup/leap/ethane-edit/LIG_GAS.crd"

#top_file = "/home/julien/projects/cyps/sire/setup/leap/benzene/LIG_GAS.top"
#crd_file = "/home/julien/projects/cyps/sire/setup/leap/benzene/LIG_GAS.crd"

#top_file = "/home/julien/projects/cyps/sire/setup/leap/benzene-distorted/LIG_GAS.top"
#crd_file = "/home/julien/projects/cyps/sire/setup/leap/benzene-distorted/LIG_GAS.crd"

#top_file = "/home/julien/projects/cyps/sire/setup/leap/pentane/LIG_GAS.top"
#crd_file = "/home/julien/projects/cyps/sire/setup/leap/pentane/LIG_GAS.crd"

#top_file = "/home/julien/projects/cyps/sire/setup/leap/dimedone/LIG_GAS.top"
#crd_file = "/home/julien/projects/cyps/sire/setup/leap/dimedone/LIG_GAS.crd"

#top_file = "/home/julien/projects/cyps/sire/setup/leap/cyclopropane/LIG_GAS.top"
#crd_file = "/home/julien/projects/cyps/sire/setup/leap/cyclopropane/LIG_GAS.crd"

#top_file = "/home/julien/projects/cyps/sire/setup/leap/cyclohexane/LIG_GAS.top"
#crd_file = "/home/julien/projects/cyps/sire/setup/leap/cyclohexane/LIG_GAS.crd"

#top_file = "/home/julien/projects/cyps/sire/setup/leap/cyclohexane-edit/LIG_GAS.top"
#crd_file = "/home/julien/projects/cyps/sire/setup/leap/cyclohexane-edit/LIG_GAS.crd"

#top_file = "/home/julien/projects/cyps/sire/setup/leap/LI8/COMPLEX_GAS.top"
#crd_file = "/home/julien/projects/cyps/sire/setup/leap/LI8/COMPLEX_GAS.crd"
#
# SOLUTE solvated, non periodic
#
#Etot   =     -4270.5482  EKtot   =         0.0000  EPtot      =     -4270.5482
# BOND   =       246.3819  ANGLE   =        39.9317  DIHED      =        18.2875
#1-4 NB =        17.8213  1-4 EEL =      -100.7708  VDWAALS    =      1484.6639
# EELEC  =     -5976.8636  EHBOND  =        0.0000  RESTRAINT  =         0.0000
#KCMT  =         0.0000  VIRIAL  =         0.0000  VOLUME     =    704972.9702

# SOLUTE in gas
# NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
# Etot   =      1035.2593  EKtot   =         0.0000  EPtot      =      1035.2593
# BOND   =       246.3483  ANGLE   =        39.9317  DIHED      =        18.2875
# 1-4 NB =        17.8213  1-4 EEL =      -100.7804  VDWAALS    =       800.6806
# EELEC  =        12.9704  EHBOND  =         0.0000  RESTRAINT  =         0.0000
# EKCMT  =         0.0000  VIRIAL  =         0.0000  VOLUME     =    291200.3229
#                                                    Density    =         0.0012
# Methanol gas phase
# NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
# Etot   =         4.5332  EKtot   =         0.0000  EPtot      =         4.5332
# BOND   =         0.1647  ANGLE   =         0.3909  DIHED      =         0.1363
# 1-4 NB =         0.0000  1-4 EEL =         3.8414  VDWAALS    =         0.0000
# EELEC  =         0.0000  EHBOND  =         0.0000  RESTRAINT  =         0.0000
# EKCMT  =         0.0000  VIRIAL  =         0.0000  VOLUME     =    239214.0118
#                                                    Density    =         0.0002
# Ethane
# NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
# Etot   =         1.8211  EKtot   =         0.0000  EPtot      =         1.8211
# BOND   =         0.0352  ANGLE   =         0.3727  DIHED      =         0.3701
# 1-4 NB =         0.1184  1-4 EEL =         0.9247  VDWAALS    =         0.0000
# EELEC  =         0.0000  EHBOND  =         0.0000  RESTRAINT  =         0.0000
# EKCMT  =         0.0000  VIRIAL  =         0.0000  VOLUME     =    248073.6160
#                                                    Density    =         0.0002
# Benzene
# NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
# Etot   =         6.9367  EKtot   =         0.0000  EPtot      =         6.9367
# BOND   =         1.4488  ANGLE   =         0.0003  DIHED      =         0.0000
# 1-4 NB =         3.3090  1-4 EEL =        -0.1645  VDWAALS    =        -0.1587
# EELEC  =         2.5019  EHBOND  =         0.0000  RESTRAINT  =         0.0000
# EKCMT  =         0.0000  VIRIAL  =         0.0000  VOLUME     =    251003.3702
#                                                    Density    =         0.0005
# ------------------------------------------------------------------------------
# benzene-distorted
# NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
# Etot   =        43.3166  EKtot   =         0.0000  EPtot      =        43.3166
# BOND   =        12.1533  ANGLE   =         5.1175  DIHED      =        20.4993
# 1-4 NB =         3.3749  1-4 EEL =        -0.1636  VDWAALS    =        -0.1592
# EELEC  =         2.4944  EHBOND  =         0.0000  RESTRAINT  =         0.0000
#
# Pentane
# NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
# Etot   =         6.2335  EKtot   =         0.0000  EPtot      =         6.2335
# BOND   =         3.0295  ANGLE   =         0.2986  DIHED      =         0.8000
# 1-4 NB =         1.5475  1-4 EEL =        -0.3847  VDWAALS    =        -0.1756
# EELEC  =         1.1182  EHBOND  =         0.0000  RESTRAINT  =         0.0000
#
# DIMEDONE
# NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
# Etot   =       -57.1783  EKtot   =         0.0000  EPtot      =       -57.1783
# BOND   =         1.0758  ANGLE   =         1.8962  DIHED      =        10.5221
# 1-4 NB =         4.1964  1-4 EEL =      -150.4665  VDWAALS    =        -1.4352
# EELEC  =        77.0329  EHBOND  =         0.0000  RESTRAINT  =         0.0000
#
# Cyclopropane
#
# NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
# Etot   =        15.1811  EKtot   =         0.0000  EPtot      =        15.1811
# BOND   =         1.8459  ANGLE   =         1.4451  DIHED      =         6.2051
# 1-4 NB =         0.0893  1-4 EEL =         5.5957  VDWAALS    =         0.0000
# EELEC  =         0.0000  EHBOND  =         0.0000  RESTRAINT  =         0.0000
#
# Cyclohexane
# NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
# Etot   =         8.2912  EKtot   =         0.0000  EPtot      =         8.2912
# BOND   =         2.5535  ANGLE   =         0.2117  DIHED      =         2.7506
# 1-4 NB =         2.7541  1-4 EEL =        -0.8982  VDWAALS    =        -0.5223
# EELEC  =         1.4417  EHBOND  =         0.0000  RESTRAINT  =         0.0000
#
# LI8 + PROTEIN + 168 water HOH
# NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
# Etot   =      6474.5265  EKtot   =         0.0000  EPtot      =      6474.5265
# BOND   =       897.8876  ANGLE   =      5310.2620  DIHED      =      2922.5644
# 1-4 NB =       790.8755  1-4 EEL =      7702.0144  VDWAALS    =       859.0611
# EELEC  =    -12008.1385  EHBOND  =         0.0000  RESTRAINT  =         0.0000
#
# Li8 + PROTEIN + WATERBOX
# NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
# Etot   =    -47010.2216  EKtot   =         0.0000  EPtot      =    -47010.2216
# BOND   =       898.1982  ANGLE   =      5310.2620  DIHED      =      2922.5644
# 1-4 NB =       790.8755  1-4 EEL =      7702.0145  VDWAALS    =      7345.0484
# EELEC  =    -71979.1846  EHBOND  =         0.0000  RESTRAINT  =         0.0000


