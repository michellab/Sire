/********************************************   \
 *
 *  Sire - Molecular Simulation Framework
 *
 *  Copyright (C) 2009  Christopher Woods
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *  For full details of the license please see the COPYING file
 *  that should have come with this distribution.
 *
 *  You can contact the authors via the developer's mailing list
 *  at http://siremol.org
 *
\*********************************************/

#include <iostream>
#include <cmath>

#include "openmmpmefep.h"
#include "ensemble.h"

#include "SireMol/moleculegroup.h"
#include "SireMol/partialmolecule.h"
#include "SireMol/molecule.h"
#include "SireMol/atommasses.h"
#include "SireMol/atomcoords.h"
#include "SireMol/moleditor.h"

#include "SireMol/amberparameters.h"

#include "SireSystem/system.h"

#include "SireBase/variantproperty.h"

#include "SireFF/forcetable.h"

#include "SireMaths/rangenerator.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "SireUnits/units.h"
#include "SireUnits/temperature.h"
#include "SireUnits/convert.h"

// ADDED BY JM
#include "SireMol/connectivity.h"
#include "SireMol/bondid.h"
#include "SireMol/atomcharges.h"
#include "SireMM/internalff.h"
#include "SireIO/amber.h"
#include "SireMM/atomljs.h"

#include "SireVol/periodicbox.h"
#include "SireVol/triclinicbox.h"

#include "SireMove/flexibility.h"

#include "SireMaths/constants.h"

//ADDED BY GAC
#include "SireMaths/vector.h"
#include "SireMol/mgname.h"
#include "SireMol/perturbation.h"
#include "SireMM/internalperturbation.h"

#include <QDebug>
#include <QTime>
#include <boost/tuple/tuple.hpp>


using namespace SireMove;
using namespace SireSystem;
using namespace SireMol;
using namespace SireFF;
using namespace SireCAS;
using namespace SireVol;
using namespace SireBase;
using namespace SireStream;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireMM;
using namespace SireIO;
using namespace std;
using boost::tuples::tuple;


enum CutoffTypes {
    NOCUTOFF = 0,
    CUTOFFNONPERIODIC = 1,
    CUTOFFPERIODIC = 2,
    PME = 3         // only PME used here
};

enum ConstraintTypes {
    NONE = 0,
    HBONDS = 1,
    ALLBONDS = 2,
    HANGLES = 3

};

enum CombinationRules {
    ARITHMETIC = 0,
    GEOMETRIC = 1
};


enum ForceGroups {
    NONBONDED_FCG = 0,
    RECIP_FCG = 1,
    DIRECT_FCG = 2,
    CORR_FCG = 3,
    BOND_FCG = 4
};

// force group selection mask
const int group_mask = 1 << NONBONDED_FCG | 1 << RECIP_FCG
    | 1 << DIRECT_FCG | 1 << CORR_FCG;

static const RegisterMetaType<OpenMMPMEFEP> r_openmmint;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const OpenMMPMEFEP &velver)
{
    writeHeader(ds, r_openmmint, 1);

    SharedDataStream sds(ds);

    sds << velver.frequent_save_velocities << velver.molgroup << velver.solute
        << velver.solutehard << velver.solutetodummy << velver.solutefromdummy
        << velver.combiningRules
        << velver.CutoffType << velver.cutoff_distance << velver.field_dielectric
        << velver.Andersen_flag << velver.Andersen_frequency
        << velver.MCBarostat_flag << velver.MCBarostat_frequency
        << velver.ConstraintType << velver.Pressure << velver.Temperature
        << velver.platform_type << velver.Restraint_flag
        << velver.CMMremoval_frequency << velver.buffer_frequency
        << velver.energy_frequency
        << velver.device_index << velver.precision << velver.current_lambda
        << velver.coulomb_power << velver.shift_delta << velver.delta_alchemical
        << velver.alchemical_array
        << velver.Integrator_type
        << velver.friction << velver.integration_tol
        << velver.timeskip << velver.reinitialise_context
        << static_cast<const Integrator&> (velver);
    // Free OpenMM pointers??

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, OpenMMPMEFEP &velver)
{
    VersionID v = readHeader(ds, r_openmmint);

    if (v == 1) {
        SharedDataStream sds(ds);
        sds >> velver.frequent_save_velocities >> velver.molgroup
            >> velver.solute >> velver.solutehard >> velver.solutetodummy
            >> velver.solutefromdummy >> velver.combiningRules
            >> velver.CutoffType >> velver.cutoff_distance
            >> velver.field_dielectric >> velver.Andersen_flag
            >> velver.Andersen_frequency >> velver.MCBarostat_flag
            >> velver.MCBarostat_frequency >> velver.ConstraintType
            >> velver.Pressure >> velver.Temperature >> velver.platform_type
            >> velver.Restraint_flag >> velver.CMMremoval_frequency
            >> velver.buffer_frequency >> velver.energy_frequency
            >> velver.device_index >> velver.precision >> velver.current_lambda
            >> velver.coulomb_power >> velver.shift_delta >> velver.delta_alchemical
            >> velver.alchemical_array
            >> velver.Integrator_type >> velver.friction >> velver.integration_tol
            >> velver.timeskip >> velver.reinitialise_context
            >> static_cast<Integrator&> (velver);

        // Maybe....need to reinitialise from molgroup because openmm system was not serialised...
        velver.isSystemInitialised = false;
        velver.isContextInitialised = false;

        //qDebug() << " Re-initialisation of OpenMMPMEFEP from datastream";

        velver.initialise();
    }
    else
        throw version_error(v, "1", r_openmmint, CODELOC);

    return ds;
}

/** Constructor*/
OpenMMPMEFEP::OpenMMPMEFEP(bool frequent_save)
    : ConcreteProperty<OpenMMPMEFEP, Integrator>(),
      frequent_save_velocities(frequent_save),
      molgroup(MoleculeGroup()), solute(MoleculeGroup()), solutehard(MoleculeGroup()),
      solutetodummy(MoleculeGroup()), solutefromdummy(MoleculeGroup()),
      openmm_system(0), openmm_context(0), isSystemInitialised(false),
      isContextInitialised(false),
      combiningRules("arithmetic"),
      CutoffType("PME"), cutoff_distance(1.0 * nanometer), field_dielectric(78.3),
      Andersen_flag(false), Andersen_frequency(90.0), MCBarostat_flag(false),
      MCBarostat_frequency(25), ConstraintType("none"),
      Pressure(1.0 * bar), Temperature(300.0 * kelvin), platform_type("Reference"),
      Restraint_flag(false),
      CMMremoval_frequency(0), buffer_frequency(0), energy_frequency(100),
      device_index("0"), precision("single"), current_lambda(0.5), coulomb_power(0),
      shift_delta(2.0), delta_alchemical(0.001), alchemical_array(),
      finite_diff_gradients(), pot_energies(), perturbed_energies(),
      reduced_perturbed_energies(),
      forward_Metropolis(), backward_Metropolis(),
      Integrator_type("leapfrogverlet"), friction(1.0 / picosecond),
      integration_tol(0.001), timeskip(0.0 * picosecond),
      reinitialise_context(false), Debug(false)
{
}

/** Constructor using the passed molecule groups */
OpenMMPMEFEP::OpenMMPMEFEP(const MoleculeGroup &molecule_group,
                           const MoleculeGroup &solute_group, const MoleculeGroup &solute_hard,
                           const MoleculeGroup &solute_todummy, const MoleculeGroup &solute_fromdummy,
                           bool frequent_save)
    : ConcreteProperty<OpenMMPMEFEP, Integrator>(),
      frequent_save_velocities(frequent_save),
      molgroup(molecule_group), solute(solute_group), solutehard(solute_hard),
      solutetodummy(solute_todummy), solutefromdummy(solute_fromdummy),
      openmm_system(0), openmm_context(0), isSystemInitialised(false),
      isContextInitialised(false),
      combiningRules("arithmetic"),
      CutoffType("PME"), cutoff_distance(1.0 * nanometer), field_dielectric(78.3),
      Andersen_flag(false), Andersen_frequency(90.0), MCBarostat_flag(false),
      MCBarostat_frequency(25), ConstraintType("none"),
      Pressure(1.0 * bar), Temperature(300.0 * kelvin), platform_type("Reference"),
      Restraint_flag(false),
      CMMremoval_frequency(0), buffer_frequency(0), energy_frequency(100),
      device_index("0"), precision("single"), current_lambda(0.5), coulomb_power(0),
      shift_delta(2.0), delta_alchemical(0.001), alchemical_array(),
      finite_diff_gradients(), pot_energies(), perturbed_energies(),
      reduced_perturbed_energies(), forward_Metropolis(), backward_Metropolis(),
      Integrator_type("leapfrogverlet"), friction(1.0 / picosecond),
      integration_tol(0.001), timeskip(0.0 * picosecond),
      reinitialise_context(false), Debug(false)
{
}

/** Copy constructor */
OpenMMPMEFEP::OpenMMPMEFEP(const OpenMMPMEFEP &other)
    : ConcreteProperty<OpenMMPMEFEP, Integrator>(other),
      frequent_save_velocities(other.frequent_save_velocities),
      molgroup(other.molgroup), solute(other.solute), solutehard(other.solutehard),
      solutetodummy(other.solutetodummy), solutefromdummy(other.solutefromdummy),
      openmm_system(other.openmm_system), openmm_context(other.openmm_context),
      isSystemInitialised(other.isSystemInitialised),
      isContextInitialised(other.isContextInitialised),
      combiningRules(other.combiningRules),
      CutoffType(other.CutoffType), cutoff_distance(other.cutoff_distance),
      field_dielectric(other.field_dielectric), Andersen_flag(other.Andersen_flag),
      Andersen_frequency(other.Andersen_frequency),
      MCBarostat_flag(other.MCBarostat_flag),
      MCBarostat_frequency(other.MCBarostat_frequency),
      ConstraintType(other.ConstraintType),
      Pressure(other.Pressure), Temperature(other.Temperature),
      platform_type(other.platform_type),
      Restraint_flag(other.Restraint_flag),
      CMMremoval_frequency(other.CMMremoval_frequency),
      buffer_frequency(other.buffer_frequency),
      energy_frequency(other.energy_frequency), device_index(other.device_index),
      precision(other.precision), current_lambda(other.current_lambda),
      coulomb_power(other.coulomb_power), shift_delta(other.shift_delta),
      delta_alchemical(other.delta_alchemical),
      alchemical_array(other.alchemical_array),
      finite_diff_gradients(other.finite_diff_gradients),
      pot_energies(other.pot_energies),
      perturbed_energies(other.perturbed_energies),
      reduced_perturbed_energies(other.reduced_perturbed_energies),
      forward_Metropolis(other.forward_Metropolis),
      backward_Metropolis(other.backward_Metropolis),
      Integrator_type(other.Integrator_type), friction(other.friction),
      integration_tol(other.integration_tol), timeskip(other.timeskip),
      reinitialise_context(other.reinitialise_context), Debug(other.Debug)
{
}

/** Destructor */
OpenMMPMEFEP::~OpenMMPMEFEP()
{
    //delete openmm_system;
}

/** Copy assignment operator */
OpenMMPMEFEP& OpenMMPMEFEP::operator=(const OpenMMPMEFEP &other)
{
    Integrator::operator=(other);
    frequent_save_velocities = other.frequent_save_velocities;
    molgroup = other.molgroup;
    solute = other.solute;
    solutehard = other.solutehard;
    solutetodummy = other.solutetodummy;
    solutefromdummy = other.solutefromdummy;
    openmm_system = other.openmm_system;
    openmm_context = other.openmm_context;
    isSystemInitialised = other.isSystemInitialised;
    isContextInitialised = other.isContextInitialised;
    combiningRules = other.combiningRules;
    CutoffType = other.CutoffType;
    cutoff_distance = other.cutoff_distance;
    field_dielectric = other.field_dielectric;
    Andersen_flag = other.Andersen_flag;
    Andersen_frequency = other.Andersen_frequency;
    MCBarostat_flag = other.MCBarostat_flag;
    MCBarostat_frequency = other.MCBarostat_frequency;
    ConstraintType = other.ConstraintType;
    Pressure = other.Pressure;
    Temperature = other.Temperature;
    platform_type = other.platform_type;
    Restraint_flag = other.Restraint_flag;
    CMMremoval_frequency = other.CMMremoval_frequency;
    buffer_frequency = other.buffer_frequency;
    energy_frequency = other.energy_frequency;
    device_index = other.device_index;
    precision = other.precision;
    current_lambda = other.current_lambda;
    coulomb_power = other.coulomb_power;
    shift_delta = other.shift_delta;
    delta_alchemical = other.delta_alchemical;
    alchemical_array = other.alchemical_array;
    finite_diff_gradients = other.finite_diff_gradients;
    pot_energies = other.pot_energies;
    reduced_perturbed_energies = other. reduced_perturbed_energies;
    forward_Metropolis = other.forward_Metropolis;
    backward_Metropolis = other.backward_Metropolis;
    perturbed_energies = other.perturbed_energies;
    Integrator_type = other.Integrator_type;
    friction = other.friction;
    integration_tol = other.integration_tol;
    timeskip = other.timeskip;
    reinitialise_context = other.reinitialise_context;
    Debug = other.Debug;
    return *this;
}

/**
 * <Compariton operator>
 * @param other
 * @return boolean
 */
bool OpenMMPMEFEP::operator==(const OpenMMPMEFEP &other) const
{
    return frequent_save_velocities == other.frequent_save_velocities
           and isSystemInitialised == other.isSystemInitialised
           and isContextInitialised == other.isContextInitialised
           and combiningRules == other.combiningRules
           and CutoffType == other.CutoffType
           and cutoff_distance == other.cutoff_distance
           and field_dielectric == other.field_dielectric
           and Andersen_flag == other.Andersen_flag
           and Andersen_frequency == other.Andersen_frequency
           and MCBarostat_flag == other.MCBarostat_flag
           and MCBarostat_frequency == other.MCBarostat_frequency
           and ConstraintType == other.ConstraintType
           and Pressure == other.Pressure
           and Temperature == other.Temperature
           and platform_type == other.platform_type
           and Restraint_flag == other.Restraint_flag
           and CMMremoval_frequency == other.CMMremoval_frequency
           and buffer_frequency == other.buffer_frequency
           and energy_frequency == other.energy_frequency
           and device_index == other.device_index
           and precision == other.precision
           and current_lambda == other.current_lambda
           and coulomb_power == other.coulomb_power
           and shift_delta == other.shift_delta
           and delta_alchemical == other.delta_alchemical
           and alchemical_array == other.alchemical_array
           and finite_diff_gradients == other.finite_diff_gradients
           and pot_energies == other.pot_energies
           and reduced_perturbed_energies == other.reduced_perturbed_energies
           and forward_Metropolis == other.forward_Metropolis
           and backward_Metropolis == other. backward_Metropolis
           and perturbed_energies == other.perturbed_energies
           and Integrator_type == other.Integrator_type
           and friction == other.friction
           and integration_tol == other.integration_tol
           and timeskip == other.timeskip
           and reinitialise_context == other.reinitialise_context
           and Debug == other.Debug
           and Integrator::operator==(other);
}

/** Comparison operator */
bool OpenMMPMEFEP::operator!=(const OpenMMPMEFEP &other) const
{
    return not OpenMMPMEFEP::operator==(other);
}

/** Return a string representation of this integrator */
QString OpenMMPMEFEP::toString() const
{
    return QObject::tr("OpenMMPMEFEP()");
}

static void addPerParticleParameters(OpenMM::CustomNonbondedForce &force,
                                     std::vector<std::string> params)
{
    for (auto const &param : params)
        force.addPerParticleParameter(param);
}

static void addPerBondParameters(OpenMM::CustomBondForce &force,
                                 std::vector<std::string> params)
{
    for (auto const &param : params)
        force.addPerBondParameter(param);
}

static void addPerAngleParameters(OpenMM::CustomAngleForce &force,
                                  std::vector<std::string> params)
{
    for (auto const &param : params)
        force.addPerAngleParameter(param);
}

void OpenMMPMEFEP::addAndersenThermostat(OpenMM::System &system)
{
    const double converted_Temperature =
	convertTo(Temperature.value(), kelvin);

    auto thermostat = new OpenMM::AndersenThermostat(converted_Temperature,
						     Andersen_frequency);

    thermostat->setRandomNumberSeed(random_seed);
    system.addForce(thermostat);

    if (Debug) {
        qDebug() << "Andersen Thermostat set";
        qDebug() << "Temperature =" << converted_Temperature << "K";
        qDebug() << "Frequency collisions =" << Andersen_frequency << "1/ps";
    }
}

void OpenMMPMEFEP::addMCBarostat(OpenMM::System &system)
{
    const double converted_Temperature =
	convertTo(Temperature.value(), kelvin);
    const double converted_Pressure = convertTo(Pressure.value(), bar);

    auto barostat = new OpenMM::MonteCarloBarostat(converted_Pressure,
						   converted_Temperature,
						   MCBarostat_frequency);

    barostat->setRandomNumberSeed(random_seed);
    system.addForce(barostat);

    if (Debug) {
        qDebug() << "Monte Carlo Barostat set";
        qDebug() << "Temperature =" << converted_Temperature << " K";
        qDebug() << "Pressure =" << converted_Pressure << " bar";
        qDebug() << "Frequency every" << MCBarostat_frequency << " steps";
    }
}


// NOTE: only for debugging with simple non-dummy systems like ions
const bool useOffset = true; // use true for production

/*
  A simple system of an ion in a water box.  This is really for debugging
  purposes only and assumes that the parm7 and rst7 file are describing
  such a system exactly as in the code below.
 */

#include "posvel.h"
void OpenMMPMEFEP::initialise_ion(bool fullPME, bool doCharge)
{
    if (Debug) {
        qDebug() << "Initialising OpenMMPMEFEP";
        const std::string version = OpenMM::Platform::getOpenMMVersion();
        qDebug() << "OpenMM Version:" << QString::fromUtf8(version.data(),
                 version.size());
	qDebug() << "Running single ion debug system";
	qDebug() << "fullPME =" << fullPME;
	qDebug() << "useOffset =" << useOffset;
	qDebug() << "doCharge =" << doCharge;
    }

    // Create a workspace using the stored molgroup
    const MoleculeGroup moleculegroup = this->molgroup.read();

    if (moleculegroup.isEmpty())
        throw SireError::program_bug(
            QObject::tr("Cannot initialise OpenMMPMEFEP because molgroup "
			"has not been defined"),
            CODELOC);

    AtomicVelocityWorkspace ws =
        this->createWorkspace(moleculegroup).read().asA<AtomicVelocityWorkspace>();

    const int nmols = ws.nMolecules();

    int nats = 0;

    for (int i = 0; i < nmols; ++i) {
        nats = nats + ws.nAtoms(i);
    }

    if (Debug)
        qDebug() << "There are" << nats << "atoms. " << "There are" << nmols
                 << "molecules";

    OpenMM::Platform::loadPluginsFromDirectory(
        OpenMM::Platform::getDefaultPluginsDirectory());

    auto system = new OpenMM::System();

    if (Andersen_flag)
        addAndersenThermostat(*system);

    if (MCBarostat_flag)
        addMCBarostat(*system);

    if (CMMremoval_frequency > 0) {
        auto cmmotionremover =
	    new OpenMM::CMMotionRemover(CMMremoval_frequency);

        system->addForce(cmmotionremover);

        if (Debug)
            qDebug() << "\nWill remove Center of Mass motion every" <<
                CMMremoval_frequency << "steps\n";
    }

    auto recip_space = new OpenMM::NonbondedForce();
    std::vector<std::pair<int,int> > pairs;

    const double cutoff = convertTo(cutoff_distance.value(), nanometer);

    system->addForce(recip_space);
    system->setDefaultPeriodicBoxVectors(OpenMM::Vec3(boxl_ion, 0.0, 0.0),
					 OpenMM::Vec3(0.0, boxl_ion, 0.0),
					 OpenMM::Vec3(0.0, 0.0, boxl_ion));

    recip_space->setForceGroup(RECIP_FCG);
    recip_space->setNonbondedMethod(OpenMM::NonbondedForce::PME);
    recip_space->setCutoffDistance(cutoff);
    recip_space->setIncludeDirectSpace(fullPME);
    recip_space->setUseDispersionCorrection(false);

    // sodium assumed to come first
    system->addParticle(atom_Na.mass);
    const unsigned int idx =
        recip_space->addParticle(atom_Na.charge, atom_Na.sigma,
                                 atom_Na.epsilon);

    double charge_scale = 0.0;
    double sigma_scale = 0.0;
    double epsilon_scale = 0.0;

    if (doCharge) {            // discharge
        charge_scale = -atom_Na.charge;
    } else {                    // shrink
        sigma_scale = -atom_Na.sigma;
        epsilon_scale = -atom_Na.epsilon;

        // make sure the charge is zero before attempting to shrink
        recip_space->setParticleParameters(idx, 0.0,
                                           atom_Na.sigma, atom_Na.epsilon);
    }

    // linear scaling of the charges
    if (useOffset) {
	recip_space->addGlobalParameter("lambda_offset", current_lambda);
	recip_space->addParticleParameterOffset("lambda_offset", idx,
						charge_scale,
						sigma_scale, epsilon_scale);
    }

    // NOTE: can't be *new for some reason
    auto direct_space = new OpenMM::CustomNonbondedForce(GENERAL_ION);
    std::vector<double> params(10);
    auto corr_recip = new OpenMM::CustomBondForce(CORR_ION);
    QVector<bool> perturbed_energies_tmp{false, false, false, false, false,
	false, false, false, false};

    /*
     * Direct space and reciprocal correction setup through explict energy
     * expressions
     */

    if (!fullPME) {
        double tolerance_PME = recip_space->getEwaldErrorTolerance();
        double alpha_PME = (1.0 / cutoff)
            * std::sqrt(-log(2.0 * tolerance_PME));

        system->addForce(direct_space);
	perturbed_energies_tmp[0] = true;

        direct_space->setForceGroup(DIRECT_FCG);
        direct_space->setNonbondedMethod
            (OpenMM::CustomNonbondedForce::CutoffPeriodic);
        direct_space->setCutoffDistance(cutoff);

        direct_space->addGlobalParameter("lam", current_lambda);
        direct_space->addGlobalParameter("alpha_pme", alpha_PME);
	direct_space->addGlobalParameter("SPOnOff", 0.0); // not needed

	addPerParticleParameters(*direct_space,
                             {"qstart", "qend", "sigmastart", "sigmaend",
			      "epstart", "epend", "isHD", "isTD",
                              "isFD", "isSolvent"});

	if (doCharge) {
            params = {
                atom_Na.charge, 0.0,
                atom_Na.sigma, atom_Na.sigma,
                atom_Na.epsilon, atom_Na.epsilon,
                1.0, 0.0, 0.0, 0.0
            };
        } else {
            params = {
                0.0, 0.0,
                atom_Na.sigma, 0.0,
                atom_Na.epsilon, 0.0,
                1.0, 0.0, 0.0, 0.0
            };
        }

        direct_space->addParticle(params);

        system->addForce(corr_recip);
        corr_recip->setForceGroup(CORR_FCG);
	perturbed_energies_tmp[8] = true;

        corr_recip->addGlobalParameter("lam_corr", current_lambda);
        corr_recip->addGlobalParameter("alpha_pme", alpha_PME);

	addPerBondParameters(*corr_recip, {"qcstart", "qcend"});
    }

    for (int i = 0; i < nwater_ion; i++) {
        int idx_O = system->addParticle(water[0].mass);
        int idx_H1 = system->addParticle(water[1].mass);
        int idx_H2 = system->addParticle(water[2].mass);

        recip_space->addParticle(water[0].charge, water[0].sigma,
                                water[0].epsilon);
        recip_space->addParticle(water[1].charge, water[1].sigma,
                                water[1].epsilon);
        recip_space->addParticle(water[2].charge, water[2].sigma,
                                water[2].epsilon);

        if (!fullPME) {
            direct_space->addParticle({water[0].charge, water[0].charge,
                    water[0].sigma, water[0].sigma,
                    water[0].epsilon, water[0].epsilon,
                    1.0, 0.0, 0.0, 1.0});
            direct_space->addParticle({water[1].charge, water[1].charge,
                    water[1].sigma, water[1].sigma,
                    water[1].epsilon, water[1].epsilon,
                    1.0, 0.0, 0.0, 1.0});
            direct_space->addParticle({water[2].charge, water[2].charge,
                    water[2].sigma, water[2].sigma,
                    water[2].epsilon, water[2].epsilon,
                    1.0, 0.0, 0.0, 1.0});
        }

        pairs.push_back(std::make_pair(idx_O, idx_H1));
        pairs.push_back(std::make_pair(idx_O, idx_H2));
	pairs.push_back(std::make_pair(idx_H1, idx_H2));

        system->addConstraint(idx_O, idx_H1, r_OH);
        system->addConstraint(idx_O, idx_H2, r_OH);
        system->addConstraint(idx_H1, idx_H2, r_HH);
    }

    recip_space->createExceptionsFromBonds(pairs, Coul14_scale_ion,
					   LJ14_scale_ion);

    /*
     * Setup of exclusions for direct space and reciprocal space corrections
     * which are custom bond forces
     */

    if (!fullPME) {
        int p1, p2;
        std::vector<double> p1_params(10);
        std::vector<double> p2_params(10);

        double charge_prod, sigma_avg, epsilon_avg;
        double qprod_start, qprod_end;
        double Qstart_p1, Qend_p1, Qstart_p2, Qend_p2;

	unsigned int num_exceptions = recip_space->getNumExceptions();

	if (Debug)
	    qDebug() << "Number of exceptions =" << num_exceptions;

        for (unsigned int i = 0; i < num_exceptions; i++) {
            recip_space->getExceptionParameters
                (i, p1, p2, charge_prod, sigma_avg, epsilon_avg);

            direct_space->addExclusion(p1, p2);

            direct_space->getParticleParameters(p1, p1_params);
            direct_space->getParticleParameters(p2, p2_params);

            // NOTE: these are just the TIP3P atom charges
            Qstart_p1 = p1_params[0];
            Qend_p1 = p1_params[1];
            Qstart_p2 = p2_params[0];
            Qend_p2 = p2_params[1];

            qprod_start = Qstart_p1 * Qstart_p2;
            qprod_end = Qend_p1 * Qend_p2;

            corr_recip->addBond(p1, p2, {qprod_start, qprod_end});
        }
    }

    perturbed_energies = perturbed_energies_tmp;
    this->openmm_system = system;
    this->isSystemInitialised = true;
} // OpenMMPMEFEP::initialise_ion END

/**
 * initialises the openMM Free energy single topology calculation
 * Initialise must be called before anything else happens.
 */
void OpenMMPMEFEP::initialise(bool fullPME)
{
    if (Debug) {
        qDebug() << "Initialising OpenMMPMEFEP";
        const std::string version = OpenMM::Platform::getOpenMMVersion();
        qDebug() << "OpenMM Version:" << QString::fromUtf8(version.data(),
                 version.size());
        qDebug() << "fullPME =" << fullPME;
	qDebug() << "useOffset =" << useOffset;
    }

    // Create a workspace using the stored molgroup
    const MoleculeGroup moleculegroup = this->molgroup.read();

    if (moleculegroup.isEmpty())
        throw SireError::program_bug(
            QObject::tr("Cannot initialise OpenMMPMEFEP because molgroup has not been defined"),
            CODELOC);

    const MoleculeGroup solute = this->solute.read();
    const MoleculeGroup solutehard = this->solutehard.read();
    const MoleculeGroup solutetodummy = this->solutetodummy.read();
    const MoleculeGroup solutefromdummy = this->solutefromdummy.read();

    AtomicVelocityWorkspace ws =
	this->createWorkspace(moleculegroup).read().asA<AtomicVelocityWorkspace>();

    const int nmols = ws.nMolecules();

    int nats = 0;

    for (int i = 0; i < nmols; ++i) {
        nats = nats + ws.nAtoms(i);
    }

    if (Debug)
        qDebug() << "There are" << nats << "atoms. " << "There are" << nmols
		 << "molecules";

    int flag_combRules;

    if (combiningRules == "arithmetic")
        flag_combRules = ARITHMETIC;
    else if (combiningRules == "geometric")
        flag_combRules = GEOMETRIC;
    else
        throw SireError::program_bug(
            QObject::tr("The combining rules have not been specified. Possible choises: arithmetic, geometric"),
            CODELOC);

    bool flag_noperturbedconstraints = false;
    int flag_constraint;
    bool flag_constraint_water = false;

    if (ConstraintType == "none")
        flag_constraint = NONE;
    else if (ConstraintType == "hbonds")
        flag_constraint = HBONDS;
    else if (ConstraintType == "allbonds")
        flag_constraint = ALLBONDS;
    else if (ConstraintType == "hangles")
        flag_constraint = HANGLES;
    else if (ConstraintType == "hbonds-notperturbed") {
        flag_constraint = HBONDS;
        flag_noperturbedconstraints = true;
    }
    else if (ConstraintType == "none-notwater") {
        flag_constraint = NONE;
        flag_constraint_water = true;
    }
    else
        throw SireError::program_bug(
            QObject::tr("The Constraints method has not been specified."
                        "Possible choises: none, hbonds, allbonds, hangles, hbonds-notperturbed, none-notwater"),
            CODELOC);

    // Load Plugins from the OpenMM standard Plugin Directory
    OpenMM::Platform::loadPluginsFromDirectory(
        OpenMM::Platform::getDefaultPluginsDirectory());

    // the system will hold all atoms
    auto system_openmm = new OpenMM::System();

    // Andersen thermostat
    if (Andersen_flag)
        addAndersenThermostat(*system_openmm);

    // Monte Carlo Barostat
    if (MCBarostat_flag)
        addMCBarostat(*system_openmm);

    if (CMMremoval_frequency > 0) {
        auto cmmotionremover = new OpenMM::CMMotionRemover(CMMremoval_frequency);

        system_openmm->addForce(cmmotionremover);
    }

    system_openmm->setDefaultPeriodicBoxVectors(OpenMM::Vec3(6, 0, 0),
						OpenMM::Vec3(0, 6, 0),
						OpenMM::Vec3(0, 0, 6));

    const double converted_cutoff_distance = convertTo(cutoff_distance.value(),
						       nanometer);

    // HHL
    // Use NonbondedForce to compute Ewald reciprocal and self terms
    // Direct space and LJ need to be implemented via expressions to
    // custom forces, see above
    auto recip_space = new OpenMM::NonbondedForce();
    recip_space->setNonbondedMethod(OpenMM::NonbondedForce::PME);
    recip_space->setCutoffDistance(converted_cutoff_distance);
    recip_space->setIncludeDirectSpace(fullPME);
    recip_space->setUseDispersionCorrection(false);

    // scale the charges in the reciprocal space
    if (useOffset) {
	recip_space->addGlobalParameter("lambda_offset", current_lambda);

	if (Debug)
	    qDebug() << "Adding lambda offset to reciprocal space";
    }

    // use default tolerance for the moment
    double tolerance_PME = recip_space->getEwaldErrorTolerance();

    // from NonbondedForceImpl.cpp
    double alpha_PME = (1.0 / converted_cutoff_distance)
                       * std::sqrt(-log(2.0 * tolerance_PME));


    /*** NON-BONDED FORCE FIELDS ***/

    /* general force field: direct space Coulomb and LJ */
    auto direct_space =
	new OpenMM::CustomNonbondedForce(GENERAL_ION);

    /* correction term for 1-2, 1-3, 1-4 exceptions in reciprocal space */
    auto custom_corr_recip =
	new OpenMM::CustomBondForce(CORR_ION);

    if (!fullPME) {
	// This ensures that also the direct space is subject to PBC
	direct_space->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffPeriodic);
	direct_space->setCutoffDistance(converted_cutoff_distance);

	direct_space->addGlobalParameter("lam", current_lambda);
	direct_space->addGlobalParameter("delta", shift_delta);
	direct_space->addGlobalParameter("n", coulomb_power);
	direct_space->addGlobalParameter("SPOnOff", 0.0);
	direct_space->addGlobalParameter("alpha_pme", alpha_PME);

	addPerParticleParameters(*direct_space,
				 {"qstart", "qend", "epstart", "epend",
				  "sigmastart", "sigmaend", "isHD", "isTD",
				  "isFD", "isSolvent"});

	custom_corr_recip->addGlobalParameter("lam_corr", current_lambda);
	custom_corr_recip->addGlobalParameter("n_corr", coulomb_power);
	custom_corr_recip->addGlobalParameter("alpha_pme", alpha_PME);
	custom_corr_recip->addGlobalParameter("cutoff", converted_cutoff_distance);

	addPerBondParameters(*custom_corr_recip, {"qcstart", "qcend"});
    }

    /*** BUILD OpenMM SYSTEM ***/

    /* add all atoms to the system */
    for (int i = 0; i < nmols; ++i) {
        const int nats_mol = ws.nAtoms(i);
        const double *m = ws.massArray(i);

        for (int j = 0; j < nats_mol; ++j) {
            system_openmm->addParticle(m[j]);
        }
    }

    int num_atoms_till_i = 0;

    QVector<bool> perturbed_energies_tmp{false, false, false, false, false,
	false, false, false, false};

    // the default AMBER 1-4 scaling factors
    double const Coulomb14Scale = 1.0 / 1.2;
    double const LennardJones14Scale = 1.0 / 2.0;

    std::vector<std::pair<int, int> > bondPairs;

    // A list of 1,4 atom pairs with non default scale factors
    // for each entry, first pair has pair of indices, second has pair of
    // scale factors
    QHash< QPair<int, int>, QPair<double, double> > custom14pairs;

    for (int i = 0; i < nmols; i++) {
        const Vector *c = ws.coordsArray(i);

        Molecule molecule = moleculegroup.moleculeAt(i).molecule();

        int num_atoms_molecule = molecule.nAtoms();

        std::vector<double> custom_non_bonded_params(10);

        if (Debug)
            qDebug() << "Molecule number = " << i << "name =" << molecule.name();

        // NON BONDED TERMS

        // The atomic parameters
        AtomLJs atomvdws = molecule.property("LJ").asA<AtomLJs>();
        AtomCharges atomcharges = molecule.property("charge").asA<AtomCharges>();
        QVector<SireMM::LJParameter> ljparameters = atomvdws.toVector();
        QVector<SireUnits::Dimension::Charge> charges = atomcharges.toVector();

        QVector<SireUnits::Dimension::Charge> start_charges;
        QVector<SireUnits::Dimension::Charge> final_charges;
        QVector<SireMM::LJParameter> start_LJs;
        QVector<SireMM::LJParameter> final_LJs;

        if (molecule.hasProperty("perturbations")) {
            if (Debug)
                qDebug() << "  ... molecule is perturbed";

            AtomCharges atomcharges_start =
                molecule.property("initial_charge").asA<AtomCharges>();
            AtomCharges atomcharges_final =
                molecule.property("final_charge").asA<AtomCharges>();

            start_charges = atomcharges_start.toVector();
            final_charges = atomcharges_final.toVector();

            AtomLJs atomvdws_start = molecule.property("initial_LJ").asA<AtomLJs>();
            AtomLJs atomvdws_final = molecule.property("final_LJ").asA<AtomLJs>();

            start_LJs = atomvdws_start.toVector();
            final_LJs = atomvdws_final.toVector();
        }

        int nonbond_idx = 0;
        double charge_start = 0.0, charge_final = 0.0;

	/* add non-bonded parameters to direct and reciprocal spaces
	   and restraints if applicable */

        // Iterate over all atoms in the molecules:
        // ljparameters.size() is used here as the number of atoms
        for (int j = 0; j < ljparameters.size(); j++) {
            double charge = charges[j].value();
            double sigma = ljparameters[j].sigma() * OpenMM::NmPerAngstrom;
            double epsilon = ljparameters[j].epsilon() * OpenMM::KJPerKcal;
            double charge_diff = 0.0;

            nonbond_idx = recip_space->addParticle(charge, sigma, epsilon);

            Atom atom = molecule.molecule().atoms()(j);

            if (molecule.hasProperty("perturbations")) {
                // Is atom a hard (changing charge and LJ), from dummy or to dummy type?
                bool ishard = false;
                bool istodummy = false;
                bool isfromdummy = false;

                charge_start = start_charges[j].value();
                charge_final = final_charges[j].value();

                charge_diff = charge_final - charge_start;

                if (abs(charge_diff) < 0.00001)
                    charge_diff = 0.0;

                if (useOffset && charge_diff != 0.0) {
		    // charge = charge_start + lambda_offset * charge_diff
		    recip_space->addParticleParameterOffset("lambda_offset",
							    nonbond_idx,
							    charge_diff,
							    0.0, 0.0); // sigma, epsilon not needed

                    if (Debug)
                        qDebug() << "Adding offset to atom idx" << nonbond_idx
                                 << "; charge_diff =" << charge_diff;
                }

                double sigma_start = start_LJs[j].sigma() * OpenMM::NmPerAngstrom;
                double sigma_final = final_LJs[j].sigma() * OpenMM::NmPerAngstrom;
                double epsilon_start = start_LJs[j].epsilon() * OpenMM::KJPerKcal;
                double epsilon_final = final_LJs[j].epsilon() * OpenMM::KJPerKcal;

                for (int l = 0; l < solutehard.nViews(); l++) {
                    Selector<Atom> view_atoms = solutehard.viewAt(l).atoms();

                    for (int m = 0; m < view_atoms.count(); m++) {
                        Atom view_atom = view_atoms(m);

                        if (atom == view_atom) {
                            ishard = true;
                            break;
                        }
                    }

                    if (ishard)
                        break;
                }

                // if not hard check if to_dummy
                if (!ishard) {
                    for (int l = 0; l < solutetodummy.nViews(); l++) {
                        Selector<Atom> view_atoms = solutetodummy.viewAt(l).atoms();

                        for (int m = 0; m < view_atoms.count(); m++) {
                            Atom view_atom = view_atoms(m);

                            if (atom == view_atom) {
                                istodummy = true;
                                break;
                            }
                        }//end for
                        if (istodummy)
                            break;
                    }//end for
                }

                // if not todummy, check if fromdummy
                if (!istodummy && !ishard) {
                    for (int l = 0; l < solutefromdummy.nViews(); l++) {
                        Selector<Atom> view_atoms = solutefromdummy.viewAt(l).atoms();

                        for (int m = 0; m < view_atoms.count(); m++) {
                            Atom view_atom = view_atoms(m);

                            if (atom == view_atom) {
                                isfromdummy = true;
                                break;
                            }
                        }//end for
                        if (isfromdummy)
                            break;
                    }//end for
                }

                if (ishard || istodummy || isfromdummy) {
                    custom_non_bonded_params[0] = charge_start;
                    custom_non_bonded_params[1] = charge_final;
                    custom_non_bonded_params[2] = epsilon_start;
                    custom_non_bonded_params[3] = epsilon_final;
                    custom_non_bonded_params[4] = sigma_start;
                    custom_non_bonded_params[5] = sigma_final;
                }
                else {
                    // unperturbed atoms
                    custom_non_bonded_params[0] = charge;
                    custom_non_bonded_params[1] = charge;
                    custom_non_bonded_params[2] = epsilon;
                    custom_non_bonded_params[3] = epsilon;
                    custom_non_bonded_params[4] = sigma;
                    custom_non_bonded_params[5] = sigma;
                }

                custom_non_bonded_params[6] = 0.0; // isHard
                custom_non_bonded_params[7] = 0.0; // isTodummy
                custom_non_bonded_params[8] = 0.0; // isFromdummy
                custom_non_bonded_params[9] = 0.0; // isSolventProtein

                if (ishard) {
                    custom_non_bonded_params[6] = 1.0;

                    if (Debug)
                        qDebug() << "hard solute = " << atom.index();
                }
                // JM July 13 THIS NEEDS FIXING TO DEAL WITH GROUPS THAT CONTAIN MORE THAN ONE MOLECULE
                else if (istodummy) {
                    custom_non_bonded_params[7] = 1.0;

                    if (Debug)
                        qDebug() << "to dummy solute = " << atom.index();
                }
                else if (isfromdummy) {
                    custom_non_bonded_params[8] = 1.0;

                    if (Debug)
                        qDebug() << "from dummy solute = " << atom.index();
                }
                else {
                    custom_non_bonded_params[6] = 1.0; // isHard

                    if (Debug)
                        qDebug() << " unperturbed solute atom " << atom.index();
                }
            }
            else {
                // solvent atom like hard
                custom_non_bonded_params[0] = charge;
                custom_non_bonded_params[1] = charge;
                custom_non_bonded_params[2] = epsilon;
                custom_non_bonded_params[3] = epsilon;
                custom_non_bonded_params[4] = sigma;
                custom_non_bonded_params[5] = sigma;
                custom_non_bonded_params[6] = 1.0; // isHard
                custom_non_bonded_params[7] = 0.0; // isTodummy
                custom_non_bonded_params[8] = 0.0; // isFromdummy
                custom_non_bonded_params[9] = 1.0; // isSolventProtein

                if (Debug)
                    qDebug() << "Solvent = " << atom.index();
            } // end if perturbations

            if (Debug) {
                qDebug() << "Charge start = " << custom_non_bonded_params[0];
                qDebug() << "Charge end = " << custom_non_bonded_params[1];
                qDebug() << "Eps start = " << custom_non_bonded_params[2];
                qDebug() << "Eps end = " << custom_non_bonded_params[3];
                qDebug() << "Sig start = " << custom_non_bonded_params[4];
                qDebug() << "Sig end = " << custom_non_bonded_params[5];
                qDebug() << "is Hard = " << custom_non_bonded_params[6];
                qDebug() << "is To dummy = " << custom_non_bonded_params[7];
                qDebug() << "is From dummy = " << custom_non_bonded_params[8];
                qDebug() << "is Solvent = " << custom_non_bonded_params[9] << "\n";
            }

            // Adds the custom parmaters to _all_ atoms
            // Must be in the same order as in the System
	    if (!fullPME)
		direct_space->addParticle(custom_non_bonded_params);
        }

        // single atoms like ions
        if (!molecule.hasProperty("connectivity")) {
            num_atoms_till_i = num_atoms_till_i + num_atoms_molecule;
            continue;
        }

        // The bonded parameters are stored in "amberparameters"
        AmberParameters amber_params =
            molecule.property("amberparameters").asA<AmberParameters>();
        QList<BondID> bonds_ff = amber_params.getAllBonds();
        QVector<BondID> bonds = bonds_ff.toVector();
        ResName molfirstresname = molecule.residues()(0).name();

        // Bonds
        for (int j = 0; j < bonds_ff.length(); j++) {
            BondID bond_ff = bonds_ff[j];
            QList<double> bond_params = amber_params.getParams(bond_ff);
            double r0 = bond_params[1];

            int idx0 = bonds[j].atom0().asA<AtomIdx>().value();
            int idx1 = bonds[j].atom1().asA<AtomIdx>().value();

            //Select the atom type
            QString atom0 = molecule.atom(AtomIdx(idx0)).toString();
            QString atom1 = molecule.atom(AtomIdx(idx1)).toString();

            idx0 = idx0 + num_atoms_till_i;
            idx1 = idx1 + num_atoms_till_i;

	    system_openmm->addConstraint(idx0, idx1, r0 * OpenMM::NmPerAngstrom);

	    if (Debug)
		qDebug() << "Adding constraint to atoms" << idx0 << "and"
			 << idx1 << "with distance" << r0;

            // Bond exclusion List
            bondPairs.push_back(std::make_pair(idx0, idx1));
        }

        num_atoms_till_i = num_atoms_till_i + num_atoms_molecule;
    } // end of loop over molecules

    /*** EXCEPTION HANDLING ***/

    // Exclude the 1-2, 1-3 bonded atoms from nonbonded forces, and scale
    // down 1-4 bonded atoms
    recip_space->createExceptionsFromBonds(bondPairs, Coulomb14Scale,
                                           LennardJones14Scale);

    int num_exceptions = recip_space->getNumExceptions();

    if (Debug)
        qDebug() << "Number of exceptions =" << num_exceptions;

    for (int i = 0; i < num_exceptions; i++) {
        int p1, p2;

	double qprod_start, qprod_end;
	double Qstart_p1, Qend_p1, Qstart_p2, Qend_p2;
        double charge_prod, sigma_avg, epsilon_avg;

	std::vector<double> p1_params(10);
	std::vector<double> p2_params(10);

        recip_space->getExceptionParameters(i, p1, p2, charge_prod, sigma_avg,
                                            epsilon_avg);

	if (!fullPME) {
	    direct_space->addExclusion(p1, p2);

	    direct_space->getParticleParameters(p1, p1_params);
	    direct_space->getParticleParameters(p2, p2_params);
	}

        Qstart_p1 = p1_params[0];
        Qend_p1 = p1_params[1];
        Qstart_p2 = p2_params[0];
        Qend_p2 = p2_params[1];

        qprod_start = Qstart_p1 * Qstart_p2;
        qprod_end = Qend_p1 * Qend_p2;

        if (Debug)
            qDebug() << "Exception =" << i << ", p1 =" << p1 << ", p2 ="
                     << p2 << ", charge prod =" << charge_prod
                     << ", sigma avg =" << sigma_avg << ", epsilon_avg ="
                     << epsilon_avg;

	if (!fullPME)
	    custom_corr_recip->addBond(p1, p2, {qprod_start, qprod_end});

	if (Debug)
	    qDebug() << "Adding qprod_start =" << qprod_start
		     << "and qprod_end =" << qprod_end
		     << "to correction term for bond"
		     << p1 << "-" << p2;
    } // end of loop over exceptions

    /*** add non-bonded force fields to System ***/

    unsigned int nAtoms = recip_space->getNumParticles();
    unsigned int npairs = (nAtoms * (nAtoms - 1)) / 2;

    if (Debug)
        qDebug() << "Num pairs  = " << npairs;

    system_openmm->addForce(recip_space);
    recip_space->setForceGroup(RECIP_FCG);

    if (!fullPME) {
        if (npairs != num_exceptions) {
            direct_space->setForceGroup(DIRECT_FCG);
            system_openmm->addForce(direct_space);
            perturbed_energies_tmp[0] = true;

	    qDebug() << "Number of atoms in GENERAL" << direct_space->getNumParticles();
        }

        if (custom_corr_recip->getNumBonds() != 0) {
            custom_corr_recip->setForceGroup(CORR_FCG);
            system_openmm->addForce(custom_corr_recip);
            perturbed_energies_tmp[8] = true;

	    qDebug() << "Number of bonds in CORR" << custom_corr_recip->getNumBonds();
        }
    } // if (!fullPME)

    perturbed_energies = perturbed_energies_tmp;

    this->openmm_system = system_openmm;
    this->isSystemInitialised = true;
} // OpenMMPMEFEP::initialise END

/**
 *
 * @param workspace
 * @param timestep
 */

void OpenMMPMEFEP::createContext(IntegratorWorkspace &workspace,
                                 SireUnits::Dimension::Time timestep)
{
    if (Debug) {
        qDebug() << "OpenMMPMEFEP::createContext(): isContextInitialised ="
                 << isContextInitialised << ", reinitialise_context ="
                 << reinitialise_context;
    }

    // Check that the openmm system has been initialised
    // !! Should check that the workspace is compatible with molgroup
    if (not this->isSystemInitialised) {
        qDebug() << "Not initialised ! ";
        throw SireError::program_bug(
            QObject::tr("OpenMMPMEFEP should have been initialised before calling integrate."),
            CODELOC);
    }

    OpenMM::System *system_openmm = openmm_system;

    int nats = system_openmm->getNumParticles();

    if (Debug)
        qDebug() << "openmm nats " << nats;

    // Integrator

    const double dt = convertTo(timestep.value(), picosecond);
    const double converted_Temperature = convertTo(Temperature.value(), kelvin);
    const double converted_friction = convertTo(friction.value(), picosecond);
    OpenMM::Integrator *integrator_openmm = NULL;

    if (!isContextInitialised || (isContextInitialised && reinitialise_context)) {
        if (Integrator_type == "leapfrogverlet")
            integrator_openmm = new OpenMM::VerletIntegrator(dt);
        else if (Integrator_type == "variableleapfrogverlet")
            integrator_openmm = new OpenMM::VariableVerletIntegrator
            (integration_tol); // integration tolerance error unitless
        else if (Integrator_type == "langevin")
            integrator_openmm = new OpenMM::LangevinIntegrator
            (converted_Temperature, converted_friction, dt);
        else if (Integrator_type == "langevinmiddle")
            integrator_openmm = new OpenMM::LangevinMiddleIntegrator
            (converted_Temperature, converted_friction, dt);
        else if (Integrator_type == "variablelangevin")
            integrator_openmm = new OpenMM::VariableLangevinIntegrator
            (converted_Temperature, converted_friction, integration_tol);
        else if (Integrator_type == "brownian")
            integrator_openmm = new OpenMM::BrownianIntegrator
            (converted_Temperature, converted_friction, dt);
        else
            throw SireError::program_bug
            (QObject::tr("The user defined Integrator type is not "
                         "supported. Available types are leapfrogverlet, "
                         "variableleapfrogverlet, langevin, langevinmiddle"
                         "variablelangevin, brownian"), CODELOC);

        if (Debug) {
            qDebug() << "Using Integrator:" << Integrator_type;
            qDebug() << "Integration step =" << dt << " ps";

            if (Integrator_type == "variablelangevin"
                    || Integrator_type == "variableleapfrogverlet") {
                qDebug() << "Integration Tol = " << integration_tol;
            }

            if (Integrator_type == "langevin" || Integrator_type == "variablelangevin"
                    || Integrator_type == "brownian") {
                qDebug() << "Converted Friction = " << converted_friction << "1/ps";
            }
        }

        OpenMM::Platform& platform_openmm = OpenMM::Platform::getPlatformByName(
                                                platform_type.toStdString());

        if (platform_type == "OpenCL") {
            const std::string prop = std::string("OpenCLDeviceIndex");
            const std::string prec = std::string("OpenCLPrecision");

            platform_openmm.setPropertyDefaultValue(prop, device_index.toStdString());
            platform_openmm.setPropertyDefaultValue(prec, precision.toStdString());

            if (Debug) {
                qDebug() << "Setting up OpenCL default Index to " << device_index;
                qDebug() << "Setting up OpenCL precision to" << precision;
            }
        }
        else if (platform_type == "CUDA") {
            const std::string prop = std::string("CudaDeviceIndex");
            const std::string prec = std::string("CudaPrecision");

            platform_openmm.setPropertyDefaultValue(prop, device_index.toStdString());
            platform_openmm.setPropertyDefaultValue(prec, precision.toStdString());

            if (Debug) {
                qDebug() << "Setting up CUDA default Index to " << device_index;
                qDebug() << "Setting up CUDA precision to" << precision;
            }
        }

        delete openmm_context;

        if (Debug) {
            qDebug() << "Deleted openmm_context";
        }

        openmm_context =
            new OpenMM::Context(*system_openmm, *integrator_openmm,
                                platform_openmm);
        this->isContextInitialised = true;

        if (Debug) {
            qDebug() << "New OpenMM Context created";
        }
    }

    if (Debug)
        qDebug() << "\nUsing OpenMM platform = " <<
                 openmm_context->getPlatform().getName().c_str() << "\n";

    // Now update coordinates / velocities / dimensions with sire data
    AtomicVelocityWorkspace &ws = workspace.asA<AtomicVelocityWorkspace>();

    const System & ptr_sys = ws.system();
    const PropertyName &space_property = PropertyName("space");

    // PeriodicBox
    if (ptr_sys.property(space_property).isA<PeriodicBox>()) {
        const PeriodicBox &space = ptr_sys.property(space_property).asA<PeriodicBox>();

        const double Box_x_Edge_Length = space.dimensions()[0] *
                                         OpenMM::NmPerAngstrom; //units in nm
        const double Box_y_Edge_Length = space.dimensions()[1] *
                                         OpenMM::NmPerAngstrom; //units in nm
        const double Box_z_Edge_Length = space.dimensions()[2] *
                                         OpenMM::NmPerAngstrom; //units in nm

        if (Debug)
            qDebug() << "\nBOX SIZE [A] = (" << space.dimensions()[0] << " , " <<
                     space.dimensions()[1] << " ,  " << space.dimensions()[2] << ")\n\n";

        //Set Periodic Box Condition

        system_openmm->setDefaultPeriodicBoxVectors(
            OpenMM::Vec3(Box_x_Edge_Length, 0, 0),
            OpenMM::Vec3(0, Box_y_Edge_Length, 0),
            OpenMM::Vec3(0, 0, Box_z_Edge_Length));

        openmm_context->setPeriodicBoxVectors(OpenMM::Vec3(Box_x_Edge_Length, 0, 0),
                                              OpenMM::Vec3(0, Box_y_Edge_Length, 0),
                                              OpenMM::Vec3(0, 0, Box_z_Edge_Length));
    }
    // TriclinicBox
    else if (ptr_sys.property(space_property).isA<TriclinicBox>()) {
        const TriclinicBox &space = ptr_sys.property(
                                        space_property).asA<TriclinicBox>();

        // Get the three triclinic box vectors.
        // FIXME: not good practice of using auto here
        const auto v0 = space.vector0();
        const auto v1 = space.vector1();
        const auto v2 = space.vector2();

        // Get cell matrix components in nm.
        const double xx = v0.x() * OpenMM::NmPerAngstrom;
        const double xy = v0.y() * OpenMM::NmPerAngstrom;
        const double xz = v0.z() * OpenMM::NmPerAngstrom;
        const double yx = v1.x() * OpenMM::NmPerAngstrom;
        const double yy = v1.y() * OpenMM::NmPerAngstrom;
        const double yz = v1.z() * OpenMM::NmPerAngstrom;
        const double zx = v2.x() * OpenMM::NmPerAngstrom;
        const double zy = v2.y() * OpenMM::NmPerAngstrom;
        const double zz = v2.z() * OpenMM::NmPerAngstrom;

        system_openmm->setDefaultPeriodicBoxVectors(
            OpenMM::Vec3(xx, xy, xz),
            OpenMM::Vec3(yx, yy, yz),
            OpenMM::Vec3(zx, zy, zz));

        openmm_context->setPeriodicBoxVectors(OpenMM::Vec3(xx, xy, xz),
                                              OpenMM::Vec3(yx, yy, yz),
                                              OpenMM::Vec3(zx, zy, zz));
    }

    openmm_context->reinitialize();

    //OpenMM vector coordinate
    std::vector<OpenMM::Vec3> positions_openmm(nats);

    //OpenMM vector momenta
    std::vector<OpenMM::Vec3> velocities_openmm(nats);

    // Conversion factor because sire units of time are in AKMA, whereas OpenMM uses picoseconds

    double AKMAPerPs = 0.04888821;
    double PsPerAKMA = 1.0 / AKMAPerPs;

    const int nmols = ws.nMolecules();

    int system_index = 0;

    for (int i = 0; i < nmols; ++i) {
        const int nats_mol = ws.nAtoms(i);

        Vector *c = ws.coordsArray(i);
        Vector *p = ws.momentaArray(i);
        const double *m = ws.massArray(i);

        for (int j = 0; j < nats_mol; ++j) {
            positions_openmm[system_index] = OpenMM::Vec3(c[j].x() *
                                             (OpenMM::NmPerAngstrom), c[j].y() * (OpenMM::NmPerAngstrom),
                                             c[j].z() * (OpenMM::NmPerAngstrom));

            if (m[j] == 0.0)
                qDebug() << "\nWARNING - THE MASS OF PARTICLE " << system_index << " is ZERO\n";

            if (m[j] > SireMaths::small) {
                velocities_openmm[system_index] = OpenMM::Vec3(p[j].x() / m[j] *
                                                  (OpenMM::NmPerAngstrom) * PsPerAKMA,
                                                  p[j].y() / m[j] * (OpenMM::NmPerAngstrom) * PsPerAKMA,
                                                  p[j].z() / m[j] * (OpenMM::NmPerAngstrom) * PsPerAKMA);
            }
            else {
                velocities_openmm[system_index] = OpenMM::Vec3(0.0, 0.0, 0.0);
            }

            if (Debug) {
                qDebug() << "Particle num = " << system_index;
                qDebug() << "Particle mass = " << m[j];
                qDebug() << "X = " << positions_openmm[system_index][0] * OpenMM::AngstromsPerNm
                         << " A" <<
                         " Y = " << positions_openmm[system_index][1] * OpenMM::AngstromsPerNm << " A" <<
                         " Z = " << positions_openmm[system_index][2] * OpenMM::AngstromsPerNm << " A";
                qDebug() << "Vx = " << velocities_openmm[system_index][0] << " Vy = " <<
                         velocities_openmm[system_index][1] << " Vz = " <<
                         velocities_openmm[system_index][2] << "\n";
            }
            system_index++;
        }
    }

    if (system_index != nats) {
        if (Debug)
            qDebug() << " system_index " << system_index << " nats " << nats;
        throw SireError::program_bug(
            QObject::tr("The number of atoms in the openmm system does not match the number of atoms in the sire workspace"),
            CODELOC);
    }

    openmm_context->setPositions(positions_openmm);
    openmm_context->setVelocities(velocities_openmm);
}

void OpenMMPMEFEP::destroyContext()
{
    if (this->isContextInitialised) {
        delete openmm_context;
        openmm_context = 0;
        this->isContextInitialised = false;
    }
}

MolarEnergy OpenMMPMEFEP::getPotentialEnergy(const System &system)
{
    IntegratorWorkspacePtr ws = this->createWorkspace(molgroup);
    ws.edit().setSystem(system);

    createContext(ws.edit(), 2 * femtosecond);

    OpenMM::State state_openmm = openmm_context->getState(OpenMM::State::Energy);
    MolarEnergy nrg = state_openmm.getPotentialEnergy() *
                      kJ_per_mol; // convert to kcal/mol

    this->destroyContext();

    return nrg;
}

/**
 * <Runs an energy Minimization on the current system.>
 * minimizeEnergy will find the nearest local potential energy minimum,
 * given the current Sire::System. It calls the
 * LocalEnergyMinimizer :: minimize() function of OpenMM.
 * @param system                Sire System including molegroup, forcefield
 *                              positions etc
 * @param tolerance             Default = 1. This specifies how precisely the
 * energy minimum must be located. Minimisation will be halted once the
 * root-mean-square value of all force components reaches this tolerance.
 * @param max_iteration         Default = 1000. this specifies the number of
 * iterations are run for the minimisation. If max_iteration = 0, the
 * iteration will run until convergence.
 *
 * @return                      Sire System, with the updated energy
 * minimised coordinates.
 */
System OpenMMPMEFEP::minimiseEnergy(System &system, double tolerance = 1.0e-10,
                                    int max_iteration = 1)
{
    const MoleculeGroup moleculegroup = this->molgroup.read();
    IntegratorWorkspacePtr workspace = this->createWorkspace(moleculegroup);

    if (system.nMolecules() != moleculegroup.nMolecules()) {
        std::cerr << "Number of molecules do not agree!";
        exit(1);
    }

    workspace.edit().setSystem(system);

    // Use helper function to create a Context
    SireUnits::Dimension::Time timestep = 0.0 * picosecond;
    createContext(workspace.edit(), timestep);

    int stateTypes = OpenMM::State::Positions;

    if (Debug)
        stateTypes |= OpenMM::State::Energy;

    OpenMM::State state_openmm = openmm_context->getState(stateTypes);
    std::vector<OpenMM::Vec3> old_positions_openmm = state_openmm.getPositions();

    if (Debug) {
	MolarEnergy Epot = state_openmm.getPotentialEnergy() * kJ_per_mol;

        qDebug() << "Total energy before minimisation:" << Epot
                 << "kcal/mol at lambda =" << current_lambda;

	const OpenMM::State state1 = openmm_context->getState
	    (OpenMM::State::Energy, false, 1 << RECIP_FCG);
	const OpenMM::State state2 = openmm_context->getState
	    (OpenMM::State::Energy, false, 1 << DIRECT_FCG);
	const OpenMM::State state4 = openmm_context->getState
	    (OpenMM::State::Energy, false, 1 << CORR_FCG);

	qDebug() << "Reciprocal energy ="
                 << state1.getPotentialEnergy() * kJ_per_mol;
	qDebug() << "Direct energy ="
                 << state2.getPotentialEnergy() * kJ_per_mol;
	qDebug() << "Correction energy ="
                 << state4.getPotentialEnergy() * kJ_per_mol;

    }

    // Step 2 minimise
    OpenMM::LocalEnergyMinimizer::minimize(*openmm_context, tolerance,
                                           max_iteration);

    // Step 3 update the positions in the system
    state_openmm = openmm_context->getState(stateTypes);
    std::vector<OpenMM::Vec3> positions_openmm = state_openmm.getPositions();

    // Recast to atomicvelocityworkspace because want to use commitCoordinates() method to update system
    AtomicVelocityWorkspace &ws = workspace.edit().asA<AtomicVelocityWorkspace>();
    const int nmols = ws.nMolecules();
    int k = 0;

    for (int i = 0; i < nmols; i++) {
        Vector *sire_coords = ws.coordsArray(i);
        for (int j = 0; j < ws.nAtoms(i); j++) {
            sire_coords[j] = Vector(positions_openmm[j + k][0] * (OpenMM::AngstromsPerNm),
                                    positions_openmm[j + k][1] * (OpenMM::AngstromsPerNm),
                                    positions_openmm[j + k][2] * (OpenMM::AngstromsPerNm));
            if (Debug) {
                qDebug() << "oX =" << old_positions_openmm[j + k][0] * OpenMM::AngstromsPerNm <<
                         " A" <<
                         " oY =" << old_positions_openmm[j + k][1] * OpenMM::AngstromsPerNm << " A" <<
                         " oZ =" << old_positions_openmm[j + k][2] * OpenMM::AngstromsPerNm << " A";
                qDebug() << "nX =" << positions_openmm[j + k][0] * OpenMM::AngstromsPerNm <<
                         " A" <<
                         " nY =" << positions_openmm[j + k][1] * OpenMM::AngstromsPerNm << " A" <<
                         " nZ =" << positions_openmm[j + k][2] * OpenMM::AngstromsPerNm << " A";
            }
        }
        k = k + ws.nAtoms(i);
    }

    // This causes the workspace to update the system coordinates with the
    // contents of *sire_coords. Note that velocities aren't touched.
    ws.commitCoordinates();

    if (Debug) {
        MolarEnergy Epot = state_openmm.getPotentialEnergy() * kJ_per_mol;

        qDebug() << "Total energy after minimisation:" << Epot
                 << "kcal/mol at lambda =" << current_lambda;
    }

    // Step 4 delete the context
    // JM 04/15 FIXME: See comment above at step 1
    this->destroyContext();
    // Step 5. Return pointer to the workspace's system
    const System & ptr_sys = ws.system();
    return ptr_sys;
}

/**
 * annealSystemToLambda will anneal the system to the current alchemical lambda
 * value of the system
 * @param system                Sire System including molegroup, forcefield
 *                              positions etc
 * @param timestep              Default = 0.005. Time step used of the
 * equilibration to the desired lambda
 * @param annealingSteps        Default = 1000. Number of steps used for the
 * annealing
 * @return                      Sire system with updated coordinates and
 * velocities.
 */

System OpenMMPMEFEP::annealSystemToLambda(System &system,
        SireUnits::Dimension::Time anneal_step_size,
        int annealing_steps)
{
    const double AKMAPerPs = 0.04888821;

    const MoleculeGroup moleculegroup = this->molgroup.read();
    IntegratorWorkspacePtr workspace = this->createWorkspace(moleculegroup);
    //TODO: Add some sanity checks here.
    if (system.nMolecules() != moleculegroup.nMolecules()) {
        std::cerr << "Number of molecules in do not agree!";
        exit(1);
    }

    workspace.edit().setSystem(system);
    //SireUnits::Dimension::Time timestep = stepSize * picosecond;
    createContext(workspace.edit(), anneal_step_size);

    int max = ceil(current_lambda / 0.1);

    double lam = 0.0;

    for (int i = 0; i < max + 1; i++) {
        updateOpenMMContextLambda(lam);
        (openmm_context->getIntegrator()).step(annealing_steps);

        if (i == max - 1)
            lam = current_lambda;
        else
            lam = lam + 0.1;
    }
    int stateTypes = OpenMM::State::Positions | OpenMM::State::Velocities;
    OpenMM::State state_openmm = openmm_context->getState(stateTypes);

    std::vector<OpenMM::Vec3> positions_openmm = state_openmm.getPositions();
    std::vector<OpenMM::Vec3> velocities_openmm = state_openmm.getVelocities();

    // Recast to atomicvelocityworkspace because want to use commitCoordinates() method to update system
    AtomicVelocityWorkspace &ws = workspace.edit().asA<AtomicVelocityWorkspace>();

    const int nmols = ws.nMolecules();
    int k = 0;

    for (int i = 0; i < nmols; i++) {
        Vector *sire_coords = ws.coordsArray(i);
        Vector *sire_momenta = ws.momentaArray(i);
        const double *m = ws.massArray(i);

        for (int j = 0; j < ws.nAtoms(i); j++) {
            sire_coords[j] = Vector(positions_openmm[j + k][0] * (OpenMM::AngstromsPerNm),
                                    positions_openmm[j + k][1] * (OpenMM::AngstromsPerNm),
                                    positions_openmm[j + k][2] * (OpenMM::AngstromsPerNm));

            sire_momenta[j] = Vector(velocities_openmm[j + k][0] * m[j] *
                                     (OpenMM::AngstromsPerNm) * AKMAPerPs,
                                     velocities_openmm[j + k][1] * m[j] * (OpenMM::AngstromsPerNm) * AKMAPerPs,
                                     velocities_openmm[j + k][2] * m[j] * (OpenMM::AngstromsPerNm) * AKMAPerPs);
        }
        k = k + ws.nAtoms(i);
    }


    ws.commitCoordinatesAndVelocities();
    //Now we also want to update the systems box in case it is very different!
    if(MCBarostat_flag) {
        // dummy buffered dimensions vector, maybe there is better solution
        //to this than just passing an empty vector
        QVector<QVector<Vector>> dimensions;
        updateBoxDimensions(state_openmm, dimensions, ws);
    }
    this->destroyContext();
    // Step 5. Return pointer to the workspace's system
    const System & ptr_sys = ws.system();
    return ptr_sys;
}

/**
 * Main integration methods for advancing dynamics
 * @param workspace             Sire Integrator workspace
 * @param nrg_component
 * @param timestep              Default = 0.002. Integration timestep
 * @param nmoves                Number of moves
 * @param record_stats          boolean that tracks recording.
 */
void OpenMMPMEFEP::integrate(IntegratorWorkspace &workspace,
                             const Symbol &nrg_component,
                             SireUnits::Dimension::Time timestep,
                             int nmoves, bool record_stats)
{
    createContext(workspace, timestep);
    const int nats = openmm_system->getNumParticles();

    AtomicVelocityWorkspace &ws = workspace.asA<AtomicVelocityWorkspace>();

    const int nmols = ws.nMolecules();

    const double AKMAPerPs = 0.04888821;

    const double dt = convertTo(timestep.value(), picosecond);

    if (Debug) {
        qDebug() << " Doing " << nmoves << " steps of dynamics ";
    }

    int n_samples = nmoves / energy_frequency;

    if (nmoves < energy_frequency)
        throw SireError::program_bug(QObject::tr("You are requesting to save "
                                     "energy every  %1 steps, which is above the total number of "
                                     "%2 steps.").arg(energy_frequency, nmoves), CODELOC);

    if (Debug)
        qDebug() << "Number Energy Samples = " << n_samples << "\n\n";

    int coord_freq = buffer_frequency;
    int nframes = 0;
    int MAXFRAMES = 1000;

    if (coord_freq > 0) {
        nframes = (nmoves / coord_freq);
        // Check that we are saving snapshots modulo frequency
        int remainder = buffer_frequency % energy_frequency;
        if (buffer_frequency < energy_frequency or remainder != 0)
            throw SireError::program_bug(QObject::tr("You are requesting to "
                                         "buffer snapshots every %1 steps, but this must number must"
                                         " be a positive integer multiple of the frequency of saving "
                                         "energies %2").arg(buffer_frequency, energy_frequency), CODELOC);
    }

    if (buffer_frequency > nmoves)
        throw SireError::program_bug(QObject::tr("You are requesting to buffer "
                                     "snapshots every  %1 steps, which is above the total number of "
                                     "%2 steps.").arg(buffer_frequency, nmoves), CODELOC);

    // Limit excessive internal buffering
    if (coord_freq > 0) {
        if (nframes > MAXFRAMES) {
            throw SireError::program_bug(QObject::tr("You are requesting to "
                                         "buffer %1 frames, which is above the hardcoded limit "
                                         "of %2.").arg(n_samples, MAXFRAMES), CODELOC);
        }
    }
    else {
        nframes = 0;
    }

    QVector< std::vector<OpenMM::Vec3> > buffered_positions;
    QVector<QVector<Vector>> buffered_dimensions;

    OpenMM::Vec3 a;
    OpenMM::Vec3 b;
    OpenMM::Vec3 c;

    const double beta = 1.0 / (0.0083144621 * convertTo(Temperature.value(),
                               kelvin)); //mol/kJ

    int stateTypes = OpenMM::State::Energy | OpenMM::State::Positions
	| OpenMM::State::Velocities;

    OpenMM::State state_openmm;
    OpenMM::Integrator &integrator = openmm_context->getIntegrator();

    int sample_count = 1;

    if (Debug && coord_freq > 0)
        qDebug() << "Saving atom coordinates every " << coord_freq << "\n";

    if (Debug) {
        qDebug() << "kT (beta) =" << beta;
        for (int i = 0; i < perturbed_energies.size(); i++)
            qDebug() << "Perturbed energy flag index" << i << " Value = " <<
                     perturbed_energies[i];
    }

    //OpenMM vector coordinate
    std::vector<OpenMM::Vec3> positions_openmm(nats);
    //OpenMM vector momenta
    std::vector<OpenMM::Vec3> velocities_openmm(nats);
    //Time skipping
    const double time_skip = convertTo(timeskip.value(), picosecond);

    if (time_skip != 0.0) {
        if (Debug)
            qDebug() << "Time to Skip = " << time_skip << "ps";

        int new_nmoves = time_skip / dt;

        if (new_nmoves >= nmoves) {
            throw SireError::program_bug(
                QObject::tr("Time to Skip is greater than the simulation time"), CODELOC);
            exit(-1);
        }

        integrator.step(new_nmoves);

        n_samples = (nmoves - new_nmoves) / energy_frequency;

        if (coord_freq > 0)
            nframes = (nmoves - new_nmoves) / coord_freq;
    }

    bool IsFiniteNumber = true;
    double increment = delta_alchemical;
    double incr_plus = current_lambda + increment;
    double incr_minus = current_lambda - increment;
    double actual_gradient = 0.0;

    emptyContainers();

    while (sample_count <= n_samples) {
	//*********************MD STEPS****************************
        integrator.step(energy_frequency);

        state_openmm = openmm_context->getState(stateTypes, false, group_mask);
        double p_energy_lambda = state_openmm.getPotentialEnergy();

        if (Debug) {
            qDebug() << "Lambda =" << current_lambda << "Potential energy ="
                     << p_energy_lambda * OpenMM::KcalPerKJ << "kcal/mol";
        }

        IsFiniteNumber = (p_energy_lambda <= DBL_MAX && p_energy_lambda >= -DBL_MAX);

        if (!IsFiniteNumber) {
            qDebug() << "NaN or Inf has been generated along the simulation";
            exit(-1);
        }

        pot_energies.append(p_energy_lambda * OpenMM::KcalPerKJ);

        if (perturbed_energies[0]) {
            // Solvent-Solvent and Protein Protein Non Bonded OFF
            // NOTE: this can dramatically change the potential energy and so the
            //       biases (reduced energies)
            openmm_context->setParameter("SPOnOff", 1.0);
        }

        // get new state as SPOnOff was set above
        state_openmm = openmm_context->getState(stateTypes, false, group_mask);

        if (Debug)
            qDebug() << "Total Time = " << state_openmm.getTime() << " ps";

        // Because looping from 1 to n_samples
        int modulo = 1;

        if (coord_freq > 0)
            modulo = sample_count % ((coord_freq / energy_frequency));

        if (Debug)
            qDebug() << "modulo is " << modulo;

        if (coord_freq > 0 and modulo == 0) {
            if (Debug)
                qDebug() << "buffering coordinates and dimensions";

            positions_openmm = state_openmm.getPositions();
            buffered_positions.append(positions_openmm);

            if (MCBarostat_flag == true) {
                state_openmm.getPeriodicBoxVectors(a, b, c);

                Vector v0 = Vector(a[0] * OpenMM::AngstromsPerNm, a[1] * OpenMM::AngstromsPerNm,
                                   a[2] * OpenMM::AngstromsPerNm);
                Vector v1 = Vector(b[0] * OpenMM::AngstromsPerNm, b[1] * OpenMM::AngstromsPerNm,
                                   b[2] * OpenMM::AngstromsPerNm);
                Vector v2 = Vector(c[0] * OpenMM::AngstromsPerNm, c[1] * OpenMM::AngstromsPerNm,
                                   c[2] * OpenMM::AngstromsPerNm);

                QVector<Vector> lattice_vectors{v0, v1, v2};

                buffered_dimensions.append(lattice_vectors);
            }
        }

        //Computing the potential energies and gradients
        p_energy_lambda = state_openmm.getPotentialEnergy();

        //Let's calculate the gradients
        double m_forward, m_backward;
        boost::tuples::tie(actual_gradient, m_forward,
                           m_backward) = calculateGradient(incr_plus,
                                         incr_minus, p_energy_lambda, beta);

        if (alchemical_array.size() > 1) {
            //Let's calculate the biased energies
            reduced_perturbed_energies.append(computeReducedPerturbedEnergies(beta));
        }

        //Now we append all the calculated information to the useful accumulation arrays
        finite_diff_gradients.append(actual_gradient * beta);
        forward_Metropolis.append(m_forward);
        backward_Metropolis.append(m_backward);

        //RESET coupling parameter to its original value
        if (perturbed_energies[0]) {
            // Solvent-Solvent and Protein Protein Non Bonded ON
            openmm_context->setParameter("SPOnOff", 0.0);
        }

        updateOpenMMContextLambda(current_lambda);

        sample_count = sample_count + 1.0;
    } // end while (sample_count <= n_samples)

    if (time_skip != 0) {
        timeskip = SireUnits::Dimension::Time(0.0);
    }

    state_openmm = openmm_context->getState(stateTypes);
    positions_openmm = state_openmm.getPositions();
    velocities_openmm = state_openmm.getVelocities();

    // Vector of Vector of molecules that are vector of atomic coordinates...
    QVector< QVector< QVector< Vector > > > buffered_workspace(nframes);
    for (int i = 0; i < buffered_workspace.size(); i++) {
        buffered_workspace[i].resize(nmols);

        for (int j = 0; j < nmols; j++) {
            int nats = ws.nAtoms(j);
            buffered_workspace[i][j].resize(nats);
        }
    }

    int k = 0;

    for (int i = 0; i < nmols; i++) {
        Vector *sire_coords = ws.coordsArray(i);
        Vector *sire_momenta = ws.momentaArray(i);
        const double *m = ws.massArray(i);

        for (int j = 0; j < ws.nAtoms(i); j++) {
            sire_coords[j] = Vector(positions_openmm[j + k][0] * (OpenMM::AngstromsPerNm),
                                    positions_openmm[j + k][1] * (OpenMM::AngstromsPerNm),
                                    positions_openmm[j + k][2] * (OpenMM::AngstromsPerNm));

            if (Debug)
                qDebug() << "X = " << positions_openmm[j + k][0] * OpenMM::AngstromsPerNm <<
                         " A" <<
                         " Y = " << positions_openmm[j + k][1] * OpenMM::AngstromsPerNm << " A" <<
                         " Z = " << positions_openmm[j + k][2] * OpenMM::AngstromsPerNm << " A";

            for (int l = 0; l < nframes; l++) {
                //qDebug() << " i " << i << " j " << j << " k " << k << " l " << l;
                Vector buffered_atcoord = Vector(buffered_positions[l][j + k][0] *
                                                 (OpenMM::AngstromsPerNm),
                                                 buffered_positions[l][j + k][1] * (OpenMM::AngstromsPerNm),
                                                 buffered_positions[l][j + k][2] * (OpenMM::AngstromsPerNm));

                buffered_workspace[l][i][j] = buffered_atcoord;
            }

            sire_momenta[j] = Vector(velocities_openmm[j + k][0] * m[j] *
                                     (OpenMM::AngstromsPerNm) * AKMAPerPs,
                                     velocities_openmm[j + k][1] * m[j] * (OpenMM::AngstromsPerNm) * AKMAPerPs,
                                     velocities_openmm[j + k][2] * m[j] * (OpenMM::AngstromsPerNm) * AKMAPerPs);

        }
        k = k + ws.nAtoms(i);
    }

    if (nframes <= 0)
        ws.commitCoordinatesAndVelocities();
    else
        ws.commitBufferedCoordinatesAndVelocities(buffered_workspace);

    //Now the box dimensions
    if (MCBarostat_flag == true) {
        updateBoxDimensions(state_openmm, buffered_dimensions, ws);
    }
    // Clear all buffers

    buffered_workspace.clear();
    buffered_dimensions.clear();
    System & ptr_sys = ws.nonConstsystem();
    ptr_sys.mustNowRecalculateFromScratch();

    return;
}

double OpenMMPMEFEP::getPotentialEnergyAtLambda(double lambda)
{
    double curr_potential_energy = 0.0;

    updateOpenMMContextLambda(lambda);

    OpenMM::State state_openmm =
        openmm_context->getState(OpenMM::State::Energy, false, group_mask);
    curr_potential_energy = state_openmm.getPotentialEnergy();

    return curr_potential_energy;
}

void OpenMMPMEFEP::updateOpenMMContextLambda(double lambda)
{
    // nonbonded terms
    if (perturbed_energies[0]) {
        openmm_context->setParameter("lam", lambda); // 1-5 HD

	if (Debug)                                                                          qDebug() << "Updating direct space lambda tp" << lambda;
    }

    // reciprocal space corrections for 1-2, 1-3 and scaled 1-4
    if (perturbed_energies[8]) {
        openmm_context->setParameter("lam_corr", lambda);

	if (Debug)
	    qDebug() << "Updating correction lambda to" << lambda;
    }

    // 1-4 Interactions
    if (perturbed_energies[1])
        openmm_context->setParameter("lamhd", lambda); // 1-4 HD
    if (perturbed_energies[2])
        openmm_context->setParameter("lamtd", 1.0 - lambda); // 1-4 To Dummy
    if (perturbed_energies[3])
        openmm_context->setParameter("lamfd", lambda); // 1-4 From Dummy
    if (perturbed_energies[4])
        openmm_context->setParameter("lamftd", lambda); // 1-4 From Dummy to Dummy

    // bonded perturbed terms
    if (perturbed_energies[5])
        openmm_context->setParameter("lambond", lambda); // Bonds
    if (perturbed_energies[6])
        openmm_context->setParameter("lamangle", lambda); // Angles
    if (perturbed_energies[7])
        openmm_context->setParameter("lamdih", lambda); // Torsions

    // lambda for the offsets (linear scaling) of the charges in
    // reciprocal space
    if (useOffset) {
	openmm_context->setParameter("lambda_offset", lambda);

	if (Debug)
	    qDebug() << "Updating lambda_offset to" << lambda;
    }
}

boost::tuples::tuple<double, double, double> OpenMMPMEFEP::calculateGradient(
    double incr_plus, double incr_minus, double p_energy_lambda, double beta)
{
    double double_increment = incr_plus - incr_minus;
    double gradient = 0;
    double potential_energy_lambda_plus_delta;
    double potential_energy_lambda_minus_delta;
    double forward_m;
    double backward_m;
    if (incr_plus < 1.0) {
        potential_energy_lambda_plus_delta = getPotentialEnergyAtLambda(incr_plus);
    }
    if (incr_minus > 0.0) {
        potential_energy_lambda_minus_delta = getPotentialEnergyAtLambda(incr_minus);
    }
    if (incr_minus < 0.0) {
        gradient = (potential_energy_lambda_plus_delta-p_energy_lambda)
                   *2/double_increment;
        backward_m = exp(beta * (potential_energy_lambda_plus_delta - p_energy_lambda));
        forward_m = exp(-beta * (potential_energy_lambda_plus_delta - p_energy_lambda));
    }
    else if(incr_plus > 1.0) {
        gradient = -(potential_energy_lambda_minus_delta-p_energy_lambda)
                   *2/double_increment;
        backward_m = exp(-beta * (potential_energy_lambda_minus_delta -
                                  p_energy_lambda));
        forward_m = exp(beta * (potential_energy_lambda_minus_delta - p_energy_lambda));
    }
    else {
        gradient = (potential_energy_lambda_plus_delta
                    -potential_energy_lambda_minus_delta)/double_increment;

        backward_m = exp(-beta * (potential_energy_lambda_minus_delta -
                                  p_energy_lambda));
        forward_m = exp(-beta * (potential_energy_lambda_plus_delta - p_energy_lambda));
    }
    return boost::tuples::make_tuple(gradient, forward_m, backward_m);
}

QVector<double> OpenMMPMEFEP::computeReducedPerturbedEnergies(double beta)
{
    QVector<double> perturbed;

    for (auto const &lam: alchemical_array)
        perturbed.append(getPotentialEnergyAtLambda(lam) * beta);

    if (Debug) {
        for (auto const &red_en: perturbed)
            qDebug() << "bias is: " << red_en;
    }

    return perturbed;
}

void OpenMMPMEFEP::emptyContainers()
{
    finite_diff_gradients.clear();
    pot_energies.clear();
    forward_Metropolis.clear();
    backward_Metropolis.clear();
    reduced_perturbed_energies.clear();
}

void OpenMMPMEFEP::updateBoxDimensions(
    OpenMM::State &state_openmm, QVector<QVector<Vector>> &buffered_dimensions,
    AtomicVelocityWorkspace &ws)
{
    OpenMM::Vec3 a;
    OpenMM::Vec3 b;
    OpenMM::Vec3 c;

    System & ptr_sys = ws.nonConstsystem();

    // TriclinicBox.
    if (ptr_sys.property("space").isA<TriclinicBox>()) {
        state_openmm.getPeriodicBoxVectors(a, b, c);
        Vector v0 = Vector(a[0] * OpenMM::AngstromsPerNm, a[1] * OpenMM::AngstromsPerNm,
                           a[2] * OpenMM::AngstromsPerNm);
        Vector v1 = Vector(b[0] * OpenMM::AngstromsPerNm, b[1] * OpenMM::AngstromsPerNm,
                           b[2] * OpenMM::AngstromsPerNm);
        Vector v2 = Vector(c[0] * OpenMM::AngstromsPerNm, c[1] * OpenMM::AngstromsPerNm,
                           c[2] * OpenMM::AngstromsPerNm);

        System & ptr_sys = ws.nonConstsystem();
        TriclinicBox sp(v0, v1, v2);

        const QString string = "space";
        ptr_sys.setProperty(string, sp);

        /** Buffer dimensions if necessary */
        for (int k = 0; k < buffered_dimensions.size(); k++) {
            const QString buffered_space = "buffered_space_" + QString::number(k);
            TriclinicBox buff_space = TriclinicBox(buffered_dimensions[k][0],
                                                   buffered_dimensions[k][1],
                                                   buffered_dimensions[k][2]);
            ptr_sys.setProperty(buffered_space, buff_space);
        }
    }
    // PeriodicBox.
    else {
        state_openmm.getPeriodicBoxVectors(a, b, c);
        Vector new_dims = Vector(a[0] * OpenMM::AngstromsPerNm,
                                 b[1] * OpenMM::AngstromsPerNm, c[2] * OpenMM::AngstromsPerNm);
        if (Debug)
            qDebug() << " a " << a[0] << " b " << b[1] << " c " << c[2];

        System & ptr_sys = ws.nonConstsystem();
        PeriodicBox sp = ptr_sys.property("space").asA<PeriodicBox>();

        sp.setDimensions(new_dims);
        const QString string = "space";
        ptr_sys.setProperty(string, sp);

        /** Buffer dimensions if necessary */
        for (int k = 0; k < buffered_dimensions.size(); k++) {
            const QString buffered_space = "buffered_space_" + QString::number(k);
            Vector dims(buffered_dimensions[k][0].x(),
                        buffered_dimensions[k][1].y(),
                        buffered_dimensions[k][2].z());
            PeriodicBox buff_space = PeriodicBox(dims);
            ptr_sys.setProperty(buffered_space, buff_space);
        }
    }
}

/** Get the combining rules type: arithmetic, geometric*/
QString OpenMMPMEFEP::getCombiningRules(void)
{
    return combiningRules;
}

/** Set the combining rules type: arithmetic, geometric*/
void OpenMMPMEFEP::setCombiningRules(QString combining_rules)
{
    combiningRules = combining_rules;
}

/** Get the cutoff type: nocutoff, cutoffnonperiodic, cutoffperiodic */
QString OpenMMPMEFEP::getCutoffType(void)
{
    return CutoffType;
}

/** Get the cutoff distance in A */
SireUnits::Dimension::Length OpenMMPMEFEP::getCutoffDistance(void)
{
    return cutoff_distance;
}

/** Set the cutoff distance in A */
void OpenMMPMEFEP::setCutoffDistance(SireUnits::Dimension::Length distance)
{
    cutoff_distance = distance;
}

/** Get the dielectric constant */
double OpenMMPMEFEP::getFieldDielectric(void)
{
    return field_dielectric;
}

/** Set the dielectric constant */
void OpenMMPMEFEP::setFieldDielectric(double dielectric)
{
    field_dielectric = dielectric;
}

/** Set Andersen thermostat */

void OpenMMPMEFEP::setAndersen(bool andersen)
{
    Andersen_flag = andersen;
}

/** Get Andersen thermostat status on/off */
bool OpenMMPMEFEP::getAndersen(void)
{
    return Andersen_flag;
}

/** Get the Andersen Thermostat frequency collision */
double OpenMMPMEFEP::getAndersenFrequency(void)
{
    return Andersen_frequency;
}

/** Set the Andersen Thermostat frequency collision */
void OpenMMPMEFEP::setAndersenFrequency(double freq)
{
    Andersen_frequency = freq;
}

/** Get the Integrator random seed */
int OpenMMPMEFEP::getRandomSeed(void)
{
    return random_seed;
}

/** Set the Integrator random seed */
void OpenMMPMEFEP::setRandomSeed(int seed)
{
    random_seed = seed;
}

/** Get the bath Temperature */
SireUnits::Dimension::Temperature OpenMMPMEFEP::getTemperature(void)
{
    return Temperature;
}

/** Set the Temperature */
void OpenMMPMEFEP::setTemperature(SireUnits::Dimension::Temperature temperature)
{
    Temperature = temperature;
}

/** Set Monte Carlo Barostat on/off */

void OpenMMPMEFEP::setMCBarostat(bool MCBarostat)
{
    MCBarostat_flag = MCBarostat;
}

/** Get Andersen thermostat status on/off */
bool OpenMMPMEFEP::getMCBarostat(void)
{
    return MCBarostat_flag;
}

/** Get the Monte Carlo Barostat frequency in time speps */
int OpenMMPMEFEP::getMCBarostatFrequency(void)
{
    return MCBarostat_frequency;
}

/** Set the Monte Carlo Barostat frequency in time speps */
void OpenMMPMEFEP::setMCBarostatFrequency(int freq)
{
    MCBarostat_frequency = freq;

}

/** Get the Presure */
SireUnits::Dimension::Pressure OpenMMPMEFEP::getPressure(void)
{
    return Pressure;
}

/** Set the Pressure */
void OpenMMPMEFEP::setPressure(SireUnits::Dimension::Pressure pressure)
{
    Pressure = pressure;
}

/** Get the Constraint type: none, hbonds, allbonds, hangles */
QString OpenMMPMEFEP::getConstraintType(void)
{
    return ConstraintType;
}

/** Set the Constraint type: none, hbonds, allbonds, hangles */
void OpenMMPMEFEP::setConstraintType(QString constrain)
{
    ConstraintType = constrain;
}

/** Get the OpenMMMD Platform: CUDA, OpenCL, CPU */
QString OpenMMPMEFEP::getPlatform(void)
{
    return platform_type;
}

/** Set the OpenMM Platform: CUDA, OpenCL, CPU */
void OpenMMPMEFEP::setPlatform(QString platform)
{
    platform_type = platform;
}

/** Get the OpenMMMD Platform: CUDA, OpenCL, CPU */
QString OpenMMPMEFEP::getDeviceIndex(void)
{
    return device_index;
}

/** Set the OpenMM Precision */
void OpenMMPMEFEP::setPrecision(QString prec)
{
    precision = prec;
}

/** Get the OpenMMMD Precision */
QString OpenMMPMEFEP::getPrecision(void)
{
    return precision;
}

/** Set the OpenMM Platform: CUDA, OpenCL, CPU */
void OpenMMPMEFEP::setDeviceIndex(QString deviceidx)
{
    device_index = deviceidx;
}

/** Get the Restaint mode*/
bool OpenMMPMEFEP::getRestraint(void)
{
    return Restraint_flag;
}

/** Set the Retraint mode */
void OpenMMPMEFEP::setRestraint(bool Restraint)
{
    Restraint_flag = Restraint;
}

/** Get the Center of Mass motion removal frequency */
int OpenMMPMEFEP::getCMMremovalFrequency(void)
{
    return CMMremoval_frequency;
}

/** Set the Center of Mass motion removal frequency */
void OpenMMPMEFEP::setCMMremovalFrequency(int frequency)
{
    CMMremoval_frequency = frequency;
}

/** Get the frequency of buffering coordinates */
int OpenMMPMEFEP::getBufferFrequency()
{
    return buffer_frequency;
}

/** Set the Center of Mass motion removal frequency */
void OpenMMPMEFEP::setBufferFrequency(int frequency)
{
    buffer_frequency = frequency;
}

/** Get the frequency of buffering coordinates */
int OpenMMPMEFEP::getEnergyFrequency()
{
    return energy_frequency;
}

/** Set the Center of Mass motion removal frequency */
void OpenMMPMEFEP::setEnergyFrequency(int frequency)
{
    energy_frequency = frequency;
}

/** Get the alchemical value used to calculate the free energy change via TI method*/
double OpenMMPMEFEP::getAlchemicalValue(void)
{
    return current_lambda;
}

/** Set the alchemical value used to calculate the free energy change via TI method*/
void OpenMMPMEFEP::setAlchemicalValue(double lambda_value)
{
    current_lambda = max(0.0, min(1.0, lambda_value));
}

void OpenMMPMEFEP::setAlchemicalArray(QVector<double> lambda_array)
{
    for (int i =0; i< lambda_array.size(); i++) {
        alchemical_array.append(max(0.0, min(1.0, lambda_array[i])));
    }
}

/** Get the coulomb power used in the soft core potential*/
//int OpenMMPMEFEP::getCoulomb_power(void)

float OpenMMPMEFEP::getCoulombPower(void)
{
    return coulomb_power;
}

/** Set the coulomb power used in the soft core potential*/
//void OpenMMPMEFEP::setCoulomb_power(int coulomb)

void OpenMMPMEFEP::setCoulombPower(float coulomb)
{
    coulomb_power = coulomb;
}

/** Get the shift used in the soft core potential*/
double OpenMMPMEFEP::getShiftDelta(void)
{
    return shift_delta;
}

/**
 * <Set the shift used in the soft core potential>
 * @param shiftdelta
 */
void OpenMMPMEFEP::setShiftDelta(double shiftdelta)
{
    shift_delta = shiftdelta;
}

/** Get the delta alchemical used in the FEP method*/
double OpenMMPMEFEP::getDeltaAlchemical(void)
{
    return delta_alchemical;
}

/**
 * Set the delta alchemical used in the FEP method
 * @param deltaalchemical
 */
void OpenMMPMEFEP::setDeltatAlchemical(double deltaalchemical)
{
    delta_alchemical = deltaalchemical;
}

/** Calculated Gradients*/
QVector<double> OpenMMPMEFEP::getGradients(void)
{
    return finite_diff_gradients;
}

/** Average energies*/
QVector<double> OpenMMPMEFEP::getEnergies(void)
{
    return pot_energies;
}
/** Average energies*/
QVector<double> OpenMMPMEFEP::getForwardMetropolis(void)
{
    return forward_Metropolis;
}
/** Average energies*/
QVector<double> OpenMMPMEFEP::getBackwardMetropolis(void)
{
    return backward_Metropolis;
}


QVector<QVector <double> > OpenMMPMEFEP::getReducedPerturbedEnergies(void)
{
    return reduced_perturbed_energies;
}

/** Get the Integrator type*/
QString OpenMMPMEFEP::getIntegrator(void)
{
    return Integrator_type;
}

/** Set the Integrator type*/
void OpenMMPMEFEP::setIntegrator(QString intgrator)
{
    Integrator_type = intgrator;
}

/** Get the friction used in specific Integrator type*/
SireUnits::Dimension::Time OpenMMPMEFEP::getFriction(void)
{
    return friction;
}

/** Set the friction used in specific Integrator type*/
void OpenMMPMEFEP::setFriction(SireUnits::Dimension::Time thefriction)
{
    friction = thefriction;
}

/** Get the integration tolerance */
double OpenMMPMEFEP::getIntegrationTolerance(void)
{
    return integration_tol;
}

/** Set the integration tolerance*/
void OpenMMPMEFEP::setIntegrationTolerance(double tolerance)
{
    integration_tol = tolerance;
}

/** Get total time to skip*/
SireUnits::Dimension::Time OpenMMPMEFEP::getTimetoSkip(void)
{
    return timeskip;
}

/** Get total time to skip*/
void OpenMMPMEFEP::setTimetoSkip(SireUnits::Dimension::Time skip)
{
    timeskip = skip;
}

/** Set the flag to reinitialise the context*/
void OpenMMPMEFEP::setReinitialiseContext(bool reinitialise)
{
    reinitialise_context = reinitialise;
}

/** Create an empty workspace */
IntegratorWorkspacePtr OpenMMPMEFEP::createWorkspace(const PropertyMap &map)
const
{
    return IntegratorWorkspacePtr(new AtomicVelocityWorkspace(map));
}

/** Return the ensemble of this integrator */
Ensemble OpenMMPMEFEP::ensemble() const
{
    return Ensemble::NVE();
}

/** Return whether or not this integrator is time-reversible */
bool OpenMMPMEFEP::isTimeReversible() const
{
    return true;
}

/** Create a workspace for this integrator for the molecule group 'molgroup' */
IntegratorWorkspacePtr OpenMMPMEFEP::createWorkspace(const MoleculeGroup
        &molgroup, const PropertyMap &map) const
{
    return IntegratorWorkspacePtr(new AtomicVelocityWorkspace(molgroup, map));
}

const char* OpenMMPMEFEP::typeName()
{
    return QMetaType::typeName(qMetaTypeId<OpenMMPMEFEP>());
}

void OpenMMPMEFEP::setDebug(bool debug)
{
    Debug = debug;
}
