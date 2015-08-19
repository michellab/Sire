
from Sire.IO import *
from Sire.Move import *
from Sire.Mol import *
from Sire.MM import *
from Sire.System import *
from Sire.Units import *

from nose.tools import assert_almost_equal

(mols, space) = Amber().readCrdTop("../io/ala.crd", "../io/ala.top")

coul_cutoff = 10 * angstrom
lj_cutoff = 10 * angstrom

rf_diel = 78.3

switchfunc = HarmonicSwitchingFunction(coul_cutoff,coul_cutoff,lj_cutoff,lj_cutoff)

interff = InterCLJFF("interclj")
interff.setProperty("space", space)
interff.setProperty("switchingFunction", switchfunc)
interff.setUseReactionField(True)
interff.setReactionFieldDielectric(rf_diel)
interff.add(mols)

intraclj = IntraCLJFF("intraclj")
intraclj.setProperty("space", space)
intraclj.setProperty("switchingFunction", switchfunc )
intraclj.setUseReactionField(True)
intraclj.setReactionFieldDielectric(rf_diel)
intraclj.add(mols)

intraff = InternalFF("intraff")
intraff.add(mols)

system = System()
system.add(interff)
system.add(intraff)
system.add(intraclj)
system.add(mols)

system.setProperty("space", space)

def _pvt_createOpenMM(mols, temperature, pressure):
    openmm = OpenMMMDIntegrator(mols)
    openmm.setPlatform("CPU")
    openmm.setConstraintType("none")
    openmm.setCutoffType("cutoffperiodic")
    openmm.setIntegrator("leapfrogverlet")
    openmm.setFriction(0.1 * picosecond)
    openmm.setPrecision("double")
    openmm.setTimetoSkip(0*picosecond)
    openmm.setDeviceIndex("0")
    openmm.setLJDispersion(False)
    openmm.setFieldDielectric(rf_diel)
    openmm.setCMMremovalFrequency(0)
    openmm.setBufferFrequency(0)
    openmm.setRestraint(False)
    openmm.setTemperature(temperature)
    openmm.setAndersen(True)
    openmm.setAndersenFrequency(10)
    openmm.setPressure(pressure)
    openmm.setMCBarostat(True)
    openmm.setMCBarostatFrequency(25)
    openmm.initialise()
    return openmm

def test_setup(verbose = False):

    sire_nrg = system.energy().value()

    if verbose:
        print("\nInitial Sire energy = %s kcal mol-1" % sire_nrg)

    openmm = _pvt_createOpenMM(system[mols.number()], 25*celsius, 1*atm)

    openmm_nrg = openmm.getPotentialEnergy(system).value()

    if verbose:
        print("\nInitial OpenMM energy = %s kcal mol-1" % openmm_nrg)

    assert_almost_equal( sire_nrg, openmm_nrg, 1 )

def test_nve(verbose = False):

    sire_nrg = system.energy().value()

    # build OpenMM NVE integrator




if __name__ == "__main__":
    test_setup(True)
