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


// General force field
// HHL
// FIXME: disable SPOnOff and see if it works with PME
//
// NOTE: There is one single namespace for global parameters but each parameter
//       must added to the force it is used in.  This is relevant for all
//       global lambdas below because they need to be changed during MD.
//       Cutoff, delta and n could use a single name each as they are constant
//       throughout the simulation.

#define COULOMB_SHIFT "rCoul = lam_diff + r;" // can we shift?
//#define COULOMB_SHIFT "rCoul = r;"

// Direct space PME and LJ
tmpl_str OpenMMPMEFEP::GENERAL =
    "(1.0 - isSolvent1 * isSolvent2 * SPOnOff) * (U_direct + U_LJ);"

    // need to subtract scaled 1-4 interactions with erf() because those are
    // computed in reciprocal space, the same for 1-2 and 1-3
    "U_direct = %1 138.935456 * q_prod * erfc(alpha_pme*rCoul) / rCoul;"
    COULOMB_SHIFT

    "U_LJ = 4.0 * eps_avg * (TWSIX3*TWSIX3 - TWSIX3);"
    "TWSIX3 = ((sigma_avg * sigma_avg) / rLJ)^3;"
    "rLJ = delta*sigma_avg*lam_diff + r*r;"

    "lam_diff = (1.0 - lambda) * 0.1;"  // 0.1 to convert to nm
    "lambda = Logic_lam * lam + Logic_om_lam * (1.0-lam) + Logic_hard;"

    "Logic_hard = isHD1 * isHD2 * (1.0-isTD1) * (1.0-isTD2) * (1.0-isFD1) * (1.0-isFD2);"
    "Logic_om_lam = max((1.0-isHD1)*(1.0-isHD2)*isTD1*isTD2*(1.0-isFD1)*(1.0-isFD2), B_om_lam);"
    "B_om_lam = max(isHD1*(1.0-isHD2)*isTD1*(1.0-isTD2)*(1.0-isFD1)*(1.0-isFD2), C_om_lam);"
    "C_om_lam = max(isHD1*(1.0-isHD2)*(1.0-isTD1)*isTD2*(1.0-isFD1)*(1.0-isFD2), D_om_lam);"
    "D_om_lam = max((1.0-isHD1)*isHD2*isTD1*(1.0-isTD2)*(1.0-isFD1)*(1.0-isFD2), E_om_lam);"
    "E_om_lam = (1.0-isHD1)*isHD2*(1.0-isTD1)*isTD2*(1.0-isFD1)*(1.0-isFD2);"
    "Logic_lam = max((1.0-isHD1)*(1.0-isHD2)*(1.0-isTD1)*(1.0-isTD2)*isFD1*isFD2, B_lam);"
    "B_lam = max(isHD1*(1.0-isHD2)*(1.0-isTD1)*(1.0-isTD2)*isFD1*(1.0-isFD2), C_lam);"
    "C_lam = max(isHD1*(1.0-isHD2)*(1.0-isTD1)*(1.0-isTD2)*(1.0-isFD1)*isFD2, D_lam);"
    "D_lam = max((1.0-isHD1)*isHD2*(1.0-isTD1)*(1.0-isTD2)*isFD1*(1.0-isFD2), E_lam);"
    "E_lam = (1.0-isHD1)*isHD2*(1.0-isTD1)*(1.0-isTD2)*(1.0-isFD1)*isFD2;"
    "B_mix = max((1.0-isHD1)*(1.0-isHD2)*isTD1*(1.0-isTD2)*(1.0-isFD1)*isFD2, C_mix);"
    "C_mix = max((1.0-isHD1)*(1.0-isHD2)*(1.0-isTD1)*isTD2*isFD1*(1.0-isFD2), D_mix);"
    "D_mix= (1.0-isHD1)*(1.0-isHD2)*(1.0-isTD1)*isTD2*(1.0-isFD1)*isFD2;"

    "q_prod = (qend1*lam + (1.0-lam)*qstart1) * (qend2*lam + (1.0-lam)*qstart2);"
    "eps_avg = sqrt((epend1*lam+(1.0-lam)*epstart1)*(epend2*lam+(1.0-lam)*epstart2));"
    "sigma_avg=";
tmpl_str OpenMMPMEFEP::GENERAL_SIGMA[2] = {
    "0.5*((sigmaend1*lam+(1.0-lam)*sigmastart1)+(sigmaend2*lam+(1.0-lam)*sigmastart2));",
    "sqrt((sigmaend1*lam+(1.0-lam)*sigmastart1)*(sigmaend2*lam+(1.0-lam)*sigmastart2));"
};

// subtract 1-2, 1-3 and 1-4 interactions that have been calculated in reciprocal space
tmpl_str OpenMMPMEFEP::CORR_RECIP =
    // cutoff shouldn't be needed because 1-4 should be shorter than cutoff
    "-U_corr * withinCutoff;"
    "withinCutoff = step(cutoff - r);"

    // erf() instead of erfc(), see PME implementation in OpenMM
    "U_corr = 138.935456 * q_prod * erf(alpha_pme*r) / r;"

    "q_prod = lam_corr*qcend + (1-lam_corr)*qcstart;";   // this is symmetrical

// 1-4 term for shrinking atoms
// NOTE: passed-in lambda (lamtd) is 1-lambda
tmpl_str OpenMMPMEFEP::TODUMMY =
    "withinCutoff*(U_direct + U_LJ);"
    "withinCutoff = step(cutofftd - r);"

    "U_direct = %1 138.935456 * q_prod * erfc(alpha_pme*rCoul) / rCoul;"
    COULOMB_SHIFT

    "U_LJ = 4.0 * eps_avg * (TWSIX3*TWSIX3 - TWSIX3);"
    "TWSIX3 = ((sigma_avg * sigma_avg) / rLJ)^3;"
    "rLJ = deltatd*sigma_avg*lam_diff + r*r;"

    "lam_diff = (1.0 - lamtd) * 0.1;"
    "eps_avg = sqrt((1-lamtd)*(1-lamtd)*eaend + lamtd*lamtd*eastart + lamtd*(1-lamtd)*emix);"
    "q_prod = (1-lamtd)*(1-lamtd)*qpend + lamtd*lamtd*qpstart + lamtd*(1-lamtd)*qmix;"
    "sigma_avg=";
tmpl_str OpenMMPMEFEP::TODUMMY_SIGMA[2] = {
    "(1-lamtd)*saend + lamtd*sastart;",
    "sqrt((1-lamtd)*(1-lamtd)*saend + lamtd*lamtd*sastart + lamtd*(1-lamtd)*samix);"
};

// 1-4 term for growing atoms
tmpl_str OpenMMPMEFEP::FROMDUMMY =
    "withinCutoff*(U_direct + U_LJ);"
    "withinCutoff=step(cutofffd - r);"

    "U_direct = %1 138.935456 * q_prod * erfc(alpha_pme*rCoul) / rCoul;"
    COULOMB_SHIFT

    "U_LJ = 4.0 * eps_avg * (TWSIX3*TWSIX3 - TWSIX3);"
    "TWSIX3 = ((sigma_avg * sigma_avg) / rLJ)^3;"
    "rLJ = deltafd*sigma_avg*lam_diff + r*r;"

    "lam_diff = (1.0 - lamfd) * 0.1;"
    "eps_avg = sqrt(lamfd*lamfd*eaend + (1-lamfd)*(1-lamfd)*eastart + lamfd*(1-lamfd)*emix);"
    "q_prod = lamfd*lamfd*qpend + (1-lamfd)*(1-lamfd)*qpstart + lamfd*(1-lamfd)*qmix;"
    "sigma_avg=";
tmpl_str OpenMMPMEFEP::FROMDUMMY_SIGMA[2] = {
    "lamfd*saend + (1-lamfd)*sastart;",
    "sqrt(lamfd*lamfd*saend + (1-lamfd)*(1-lamfd)*sastart + lamfd*(1-lamfd)*samix);"
};

// 1-4 term between shrinking and growing atoms
tmpl_str OpenMMPMEFEP::FROMTODUMMY =
    "withinCutoff*(U_direct + U_LJ);"
    "withinCutoff = step(cutoffftd - r);"

    "U_direct = 138.935456 * q_prod * erfc(alpha_pme*rCoul) / rCoul;"
    COULOMB_SHIFT

    "U_LJ = 4.0 * eps_avg * (TWSIX3*TWSIX3 - TWSIX3);"
    "TWSIX3 = ((sigma_avg * sigma_avg) / rLJ)^3;"
    "rLJ = deltaftd*sigma_avg*0.1 + r*r;"

    "eps_avg = sqrt(lamftd*lamftd*eaend + (1-lamftd)*(1-lamftd)*eastart + lamftd*(1-lamftd)*emix);"
    "q_prod = lamftd*lamftd*qpend + (1-lamftd)*(1-lamftd)*qpstart + lamftd*(1-lamftd)*qmix;"
    "sigma_avg=";
tmpl_str OpenMMPMEFEP::FROMTODUMMY_SIGMA[2] = {
    "lamftd*saend + (1-lamftd)*sastart;",
    "sqrt(lamftd*lamftd*saend + (1-lamftd)*(1-lamftd)*sastart + lamftd*(1-lamftd)*samix);"
};

// standard 1-4 Coulomb and LJ terms, scaling is done per each exception below
tmpl_str OpenMMPMEFEP::INTRA_14_CLJ =
    "withinCutoff*(U_direct + U_LJ);"
    "withinCutoff = step(cutoffhd - r);"

    "U_direct = 138.935456 * q_prod * erfc(alpha_pme*r) / r;" // no lambda^n, no shift

    "U_LJ = 4.0 * eps_avg * ((sigma_avg/r)^12 - (sigma_avg/r)^6);" // no shift

    "eps_avg = sqrt(lamhd*lamhd*eaend + (1-lamhd)*(1-lamhd)*eastart + lamhd*(1-lamhd)*emix);"
    "q_prod = lamhd*lamhd*qpend + (1-lamhd)*(1-lamhd)*qpstart + lamhd*(1-lamhd)*qmix;"
    "sigma_avg=";
tmpl_str OpenMMPMEFEP::INTRA_14_CLJ_SIGMA[2] = {
    "lamhd*saend + (1-lamhd)*sastart;",
    "sqrt(lamhd*lamhd*saend + (1-lamhd)*(1-lamhd)*sastart + lamhd*(1-lamhd)*samix);"
};


// NOTE: only for debugging with simple non-dummy systems like ions
const bool fullPME = false;   // use false for production
const bool useOffset = true; // use true for production

/**
 * initialises the openMM Free energy single topology calculation
 * Initialise must be called before anything else happens.
 */
void OpenMMPMEFEP::initialise()
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

    if (Debug)
        qDebug() << "combiningRules =" << combiningRules;

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

    if (Debug)
        qDebug() << "Constraint Type =" << ConstraintType;

    // Load Plugins from the OpenMM standard Plugin Directory
    OpenMM::Platform::loadPluginsFromDirectory(
        OpenMM::Platform::getDefaultPluginsDirectory());

    // the system will hold all
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

        if (Debug)
            qDebug() << "\nWill remove Center of Mass motion every" <<
		CMMremoval_frequency << "steps\n";
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

    if (Debug) {
        qDebug() << "PME alpha =" << alpha_PME
                 << "computed from PME error tolerance =" << tolerance_PME;
    }

    /*** NON-BONDED FORCE FIELDS ***/

    QString lam_pre = "";

    /* general force field: direct space Coulomb and LJ */

    // This check is necessary to avoid nan errors on the GPU platform caused
    // by the calculation of 0^0
    if (coulomb_power > 0)
        lam_pre = "(lambda^n) *";

    QString general_ff = GENERAL.arg(lam_pre);
    general_ff.append(GENERAL_SIGMA[flag_combRules]);

    auto direct_space =
	new OpenMM::CustomNonbondedForce(general_ff.toStdString());

    /* correction term for 1-2, 1-3, 1-4 exceptions in reciprocal space */
    auto custom_corr_recip =
	new OpenMM::CustomBondForce(CORR_RECIP.toStdString());

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

    if (coulomb_power > 0)
        lam_pre = "(lamtd^ntd) *";

    /* to dummy force field */

    QString intra_14_todummy = TODUMMY.arg(lam_pre);
    intra_14_todummy.append(TODUMMY_SIGMA[flag_combRules]);

    auto custom_intra_14_todummy =
        new OpenMM::CustomBondForce(intra_14_todummy.toStdString());
    custom_intra_14_todummy->addGlobalParameter("lamtd", 1.0-current_lambda);
    custom_intra_14_todummy->addGlobalParameter("deltatd", shift_delta);
    custom_intra_14_todummy->addGlobalParameter("ntd", coulomb_power);
    custom_intra_14_todummy->addGlobalParameter("cutofftd",
						converted_cutoff_distance);
    custom_intra_14_todummy->addGlobalParameter("alpha_pme", alpha_PME);

    if (coulomb_power > 0)
        lam_pre = "(lamfd^nfd) *";

    /* from dummy force field */

    QString intra_14_fromdummy = FROMDUMMY.arg(lam_pre);
    intra_14_fromdummy.append(FROMDUMMY_SIGMA[flag_combRules]);

    auto custom_intra_14_fromdummy =
        new OpenMM::CustomBondForce(intra_14_fromdummy.toStdString());
    custom_intra_14_fromdummy->addGlobalParameter("lamfd", current_lambda);
    custom_intra_14_fromdummy->addGlobalParameter("deltafd", shift_delta);
    custom_intra_14_fromdummy->addGlobalParameter("nfd", coulomb_power);
    custom_intra_14_fromdummy->addGlobalParameter("cutofffd",
						  converted_cutoff_distance);
    custom_intra_14_fromdummy->addGlobalParameter("alpha_pme", alpha_PME);

    /* from and to dummy force field */

    auto intra_14_fromdummy_todummy = QString(FROMTODUMMY);
    intra_14_fromdummy_todummy.append(FROMTODUMMY_SIGMA[flag_combRules]);

    auto custom_intra_14_fromdummy_todummy =
        new OpenMM::CustomBondForce(intra_14_fromdummy_todummy.toStdString());

    custom_intra_14_fromdummy_todummy->addGlobalParameter("lamftd",
            current_lambda);
    custom_intra_14_fromdummy_todummy->addGlobalParameter("deltaftd", shift_delta);
    custom_intra_14_fromdummy_todummy->addGlobalParameter("nftd", coulomb_power);
    custom_intra_14_fromdummy_todummy->addGlobalParameter("cutoffftd",
            converted_cutoff_distance);
    custom_intra_14_fromdummy_todummy->addGlobalParameter("alpha_pme", alpha_PME);

    /* hard, perturbed atom force field */

    QString intra_14_clj(INTRA_14_CLJ);
    intra_14_clj.append(INTRA_14_CLJ_SIGMA[flag_combRules]);

    auto custom_intra_14_clj =
        new OpenMM::CustomBondForce(intra_14_clj.toStdString());

    custom_intra_14_clj->addGlobalParameter("lamhd", current_lambda);
    custom_intra_14_clj->addGlobalParameter("cutoffhd", converted_cutoff_distance);
    custom_intra_14_clj->addGlobalParameter("alpha_pme", alpha_PME);

    std::vector<std::string> paramList = {
	"qpstart", "qpend", "qmix", "eastart", "eaend", "emix",
	"sastart", "saend", "samix"
    };
    addPerBondParameters(*custom_intra_14_todummy, paramList);
    addPerBondParameters(*custom_intra_14_fromdummy, paramList);
    addPerBondParameters(*custom_intra_14_fromdummy_todummy, paramList);
    addPerBondParameters(*custom_intra_14_clj, paramList);

    if (Debug) {
        qDebug() << "\nCutoff type = " << CutoffType;
        qDebug() << "CutOff distance = " << converted_cutoff_distance << " Nm";
        qDebug() << "Dielectric constant = " << field_dielectric;
        qDebug() << "Lambda = " << current_lambda << " Coulomb Power = " <<
                 coulomb_power << " Delta Shift = " << shift_delta;
    }


    /*** BONDED FORCE FIELDS ***/

    auto bondStretch_openmm = new OpenMM::HarmonicBondForce();
    auto bondBend_openmm = new OpenMM::HarmonicAngleForce();
    auto bondTorsion_openmm = new OpenMM::PeriodicTorsionForce();

    auto solute_bond_perturbation = new OpenMM::CustomBondForce(
        "0.5*B*(r-req)^2;"
        "B=bend*lambond+(1.0-lambond)*bstart;"
        "req=rend*lambond+(1.0-lambond)*rstart");

    solute_bond_perturbation->addGlobalParameter("lambond", current_lambda);
    addPerBondParameters(*solute_bond_perturbation,
			 {"bstart", "bend", "rstart", "rend"});


    auto solute_angle_perturbation = new OpenMM::CustomAngleForce(
        "0.5*A*(theta-thetaeq)^2;"
        "A=aend*lamangle+(1.0-lamangle)*astart;"
        "thetaeq=thetaend*lamangle+(1.0-lamangle)*thetastart");

    solute_angle_perturbation->addGlobalParameter("lamangle", current_lambda);
    addPerAngleParameters(*solute_angle_perturbation,
			  {"astart", "aend", "thetastart", "thetaend"});


    /*** RESTRAINTS ***/

    OpenMM::CustomExternalForce *positionalRestraints_openmm = NULL;

    if (Restraint_flag) {
        positionalRestraints_openmm = new OpenMM::CustomExternalForce
	    (
	     "k*d2;"
	     "d2 = max(0.0, d1 - d^2);"
	     "d1 = (x-xref)^2 + (y-yref)^2  + (z-zref)^2"
	    );
        positionalRestraints_openmm->addPerParticleParameter("xref");
        positionalRestraints_openmm->addPerParticleParameter("yref");
        positionalRestraints_openmm->addPerParticleParameter("zref");
        positionalRestraints_openmm->addPerParticleParameter("k");
        positionalRestraints_openmm->addPerParticleParameter("d");

        system_openmm->addForce(positionalRestraints_openmm);

        if (Debug)
            qDebug() << "\n\nRestraint is ON\n\n";
    }

    /*** BOND LINK FORCE FIELD ***/
    /* NOTE: CustomBondForce does not (OpenMM 6.2) apply PBC checks so code will be buggy if
       restraints involve one atom that diffuses out of the box. */

    auto custom_link_bond =
        new OpenMM::CustomBondForce("kl * max(0, d - dl*dl);"
                                    "d = (r-reql) * (r-reql)");
    custom_link_bond->addPerBondParameter("reql");
    custom_link_bond->addPerBondParameter("kl");
    custom_link_bond->addPerBondParameter("dl");


    /*** BUILD OpenMM SYSTEM ***/

    int system_index = 0;

    // To avoid possible mismatch between the index in which atoms are added to the openmm system arrays and
    // their atomic numbers in sire, one array is populated while filling up the openmm global arrays
    QHash<int, int> AtomNumToOpenMMIndex;

    /* add all atoms to the system */

    for (int i = 0; i < nmols; ++i) {
        const int nats_mol = ws.nAtoms(i);

        const double *m = ws.massArray(i);

        MolNum molnum = moleculegroup.molNumAt(i);

        const ViewsOfMol &molview = moleculegroup[molnum].data();

        const Molecule &mol = molview.molecule();

        Selector<Atom> molatoms = mol.atoms();

        for (int j = 0; j < nats_mol; ++j) {
            // This adds each atom to the system via its mass
            // JM 10/16 make sure that perturbed atoms have mass of heaviest end-state
            system_openmm->addParticle(m[j]);

            Atom at = molatoms(j);
            AtomNum atnum = at.number();

            if (Debug)
                qDebug() << " openMM_index " << system_index << " Sire Atom Number "
                         << atnum.toString() << " Mass particle = " << m[j];

            AtomNumToOpenMMIndex[atnum.value()] = system_index;

            // JM Nov 12
            // The code below implements a ThreeParticleAverageSite for virtual
            // sites for EPW atoms present in a WAT residue
            // This is very AMBER specific.

            AtomName atname = at.name();

            if (Debug)
                qDebug() << " atname " << atname.value() << " mol " << i;

            if (atname == AtomName("EPW")) {
                ResName resname = at.residue().name();

                if (resname == ResName("WAT")) {
                    Atom oatom = molatoms.select(AtomName("O"));
                    Atom h1atom = molatoms.select(AtomName("H1"));
                    Atom h2atom = molatoms.select(AtomName("H2"));

                    AmberParameters amber_params =
                        mol.property("amberparameters").asA<AmberParameters>();
                    QList<BondID> bonds_ff = amber_params.getAllBonds();

                    double distoh = -1.0;
                    double disthh = -1.0;
                    double distoe = -1.0;

                    for (int k = 0; k < bonds_ff.length(); k++) {
                        BondID bond_ff = bonds_ff[k];
                        QList<double> bond_params = amber_params.getParams(bond_ff);

                        double r0 = bond_params[1];

                        AtomName at0name = mol.select(bond_ff.atom0()).name();
                        AtomName at1name = mol.select(bond_ff.atom1()).name();

                        if ((at0name == AtomName("O") and at1name == AtomName("H1"))
                                or ( at0name == AtomName("H1") and at1name == AtomName("O"))) {
                            distoh = r0;
                        }
                        else if ((at0name == AtomName("H1") and at1name == AtomName("H2"))
                                 or ( at0name == AtomName("H2") and at1name == AtomName("H1"))) {
                            disthh = r0;
                        }
                        else if ((at0name == AtomName("EPW") and at1name == AtomName("O"))
                                 or ( at0name == AtomName("O") and at1name == AtomName("EPW"))) {
                            distoe = r0;
                        }
                    }

                    if (distoh < 0.0 or disthh < 0.0 or distoe < 0.0)
                        throw SireError::program_bug(
                            QObject::tr("Could not find expected atoms in TIP4P water molecule."), CODELOC);

                    double weightH = distoe / sqrt((distoh * distoh) - (0.25 * disthh * disthh));

                    int o_index = AtomNumToOpenMMIndex[oatom.number().value()];
                    int h1_index = AtomNumToOpenMMIndex[h1atom.number().value()];
                    int h2_index = AtomNumToOpenMMIndex[h2atom.number().value()];

                    if (Debug)
                        qDebug() << "virtual site " << system_index << " o "
                                 << o_index << " h1 " << h1_index << " h2 "
                                 << h2_index << " 1 - weightH " << 1 - weightH
                                 << " weightH/2 " << weightH / 2;

                    auto vsite =
                        new OpenMM::ThreeParticleAverageSite(o_index, h1_index,
                                h2_index, 1 - weightH,
                                weightH / 2, weightH / 2);

                    system_openmm->setVirtualSite(system_index, vsite);
                }
            }

            system_index = system_index + 1;

        } // end of loop on atoms in molecule

    } // end of loop on molecules in workspace

    int num_atoms_till_i = 0;

    // JM July 13. This also needs to be changed because there could be more than one perturbed molecule
    // Molecule solutemol = solute.moleculeAt(0).molecule();
    int nions = 0;

    QVector<bool> perturbed_energies_tmp{false, false, false, false, false,
	false, false, false, false};

    // the default AMBER 1-4 scaling factors
    double const Coulomb14Scale = 1.0 / 1.2;
    double const LennardJones14Scale = 1.0 / 2.0;

    std::vector<std::pair<int, int> > bondPairs;

    // A list of 1,4 atom pairs with non default scale factors
    // for each entry, first pair has pair of indices, second has pair of scale factors
    QHash< QPair<int, int>, QPair<double, double> > custom14pairs;

    bool special_14 = false;

    for (int i = 0; i < nmols; i++) {
        const Vector *c = ws.coordsArray(i);

        Molecule molecule = moleculegroup.moleculeAt(i).molecule();

        int num_atoms_molecule = molecule.nAtoms();

        std::vector<double> custom_non_bonded_params(10);

        if (Debug)
            qDebug() << "Molecule number = " << i;

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
                qDebug() << "Molecule Perturbed number = " << i;

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

            // This really only adds the nonbonded parameters
            // The parameters need to be added in the same order as they
            // appear in the System
            nonbond_idx = recip_space->addParticle(charge, sigma, epsilon);

            Atom atom = molecule.molecule().atoms()(j);

            if (molecule.hasProperty("perturbations")) {
                // Is atom a hard (changing charge and LJ), from dummy or to dummy type?
                bool ishard = false;
                bool istodummy = false;
                bool isfromdummy = false;

                charge_start = start_charges[j].value();
                charge_final = final_charges[j].value();

                // HHL
                // Lambda scaling for 1-5+ (see exceptions below) in reciprocal
                // space complimentary to scaling in direct space
                // need to provide the parameter (lambda) and the chargeScale for
                // reciprocal PME
                charge_diff = charge_final - charge_start;

                // FIXME: really needed? const for small value
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

        if (Restraint_flag) {
            bool hasRestrainedAtoms = molecule.hasProperty("restrainedatoms");

            if (hasRestrainedAtoms) {
                Properties restrainedAtoms =
                    molecule.property("restrainedatoms").asA<Properties>();

                int nrestrainedatoms = restrainedAtoms.property(
                                           QString("nrestrainedatoms")).asA<VariantProperty>().toInt();

                if (Debug)
                    qDebug() << "nrestrainedatoms = " << nrestrainedatoms;

                for (int i = 0; i < nrestrainedatoms; i++) {
                    int atomnum = restrainedAtoms.property(QString("AtomNum(%1)").arg(
                            i)).asA<VariantProperty>().toInt();
                    double xref = restrainedAtoms.property(QString("x(%1)").arg(
                            i)).asA<VariantProperty>().toDouble();
                    double yref = restrainedAtoms.property(QString("y(%1)").arg(
                            i)).asA<VariantProperty>().toDouble();
                    double zref = restrainedAtoms.property(QString("z(%1)").arg(
                            i)).asA<VariantProperty>().toDouble();
                    double k = restrainedAtoms.property(QString("k(%1)").arg(
                                                            i)).asA<VariantProperty>().toDouble();
                    double d = restrainedAtoms.property(QString("d(%1)").arg(
                                                            i)).asA<VariantProperty>().toDouble();

                    int openmmindex = AtomNumToOpenMMIndex[atomnum];

                    if (Debug) {
                        qDebug() << "atomnum " << atomnum << " openmmindex " << openmmindex << " x " <<
                                 xref << " y " << yref << " z " << zref << " k " << k << " d " << d;
                    }

                    int posrestrdim = 5;
                    std::vector<double> params(posrestrdim);

                    params[0] = xref * OpenMM::NmPerAngstrom;
                    params[1] = yref * OpenMM::NmPerAngstrom;
                    params[2] = zref * OpenMM::NmPerAngstrom;
                    params[3] = k * (OpenMM::KJPerKcal * OpenMM::AngstromsPerNm *
                                     OpenMM::AngstromsPerNm);
                    params[4] = d * OpenMM::NmPerAngstrom;

                    positionalRestraints_openmm->addParticle(openmmindex, params);
                }
            }
        } // end of restraint flag

        // single atoms like ions

        bool hasConnectivity = molecule.hasProperty("connectivity");

        if (!hasConnectivity) {
            num_atoms_till_i = num_atoms_till_i + num_atoms_molecule;

            if (Debug) {
                qDebug() << "\nAtoms = " << num_atoms_molecule << " Num atoms till i =" <<
                         num_atoms_till_i;
                qDebug() <<
                         "*********************MONOATOMIC MOLECULE DETECTED**************************\n";
            }

            nions = nions + 1;
            continue;
        }

        // bonded terms

        QList< BondID > bond_pert_list;
        QList< BondID > bond_pert_swap_list;
        QList< AngleID > angle_pert_list;
        QList< AngleID > angle_pert_swap_list;
        QList< DihedralID > dihedral_pert_list;
        QList< DihedralID > dihedral_pert_swap_list;
        QList< ImproperID > improper_pert_list;
        QList< ImproperID > improper_pert_swap_list;

        if (solute.contains(molecule)) {
            Perturbations pert_params =
                molecule.property("perturbations").asA<Perturbations>();

            QList< PropPtr<Perturbation> > perturbation_list = pert_params.perturbations();

            std::vector<double> solute_bond_perturbation_params(4);
            std::vector<double> solute_angle_perturbation_params(4);
            std::vector<double> solute_torsion_perturbation_params(1);

            QHash<BondID, double> bond_pert_eq_list;

            for (QList< PropPtr<Perturbation> >::const_iterator it =
                        perturbation_list.constBegin(); it != perturbation_list.constEnd(); ++it) {
                const Perturbation &pert = *it;

                if (pert.isA<InternalPerturbation>()) {
                    QString str = pert.what();

                    if (str == "SireMM::TwoAtomPerturbation") {
                        const TwoAtomPerturbation &two = pert.asA<TwoAtomPerturbation>();
                        int idx0 = two.atom0().asA<AtomIdx>().value() + num_atoms_till_i;
                        int idx1 = two.atom1().asA<AtomIdx>().value() + num_atoms_till_i;
                        double rstart = two.initialForms()[Symbol("r0")].toString().toDouble();
                        double bstart = two.initialForms()[Symbol("k")].toString().toDouble();
                        double rend = two.finalForms()[Symbol("r0")].toString().toDouble();
                        double bend = two.finalForms()[Symbol("k")].toString().toDouble();

                        solute_bond_perturbation_params[0] = bstart * 2.0 * OpenMM::KJPerKcal *
                                                             OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm;
                        solute_bond_perturbation_params[1] = bend * 2.0 * OpenMM::KJPerKcal *
                                                             OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm;
                        solute_bond_perturbation_params[2] = rstart * OpenMM::NmPerAngstrom;
                        solute_bond_perturbation_params[3] = rend * OpenMM::NmPerAngstrom;

                        /* JM 10/16 Also apply this if 'no solute constraints' flag is on*/
                        if (flag_constraint == NONE) {
                            solute_bond_perturbation->addBond(idx0, idx1, solute_bond_perturbation_params);
                        }
                        else if (flag_constraint == ALLBONDS || flag_constraint == HANGLES) {
                            /* JM 10/16 ALLBONDS and HANGLES may be unwise with current free energy implementation !*/
                            double pert_eq_distance = solute_bond_perturbation_params[3] * current_lambda
                                                      + (1.0 - current_lambda) * solute_bond_perturbation_params[2];
                            system_openmm->addConstraint(idx0, idx1, pert_eq_distance);
                            bond_pert_eq_list.insert(BondID(two.atom0(), two.atom1()),
                                                     pert_eq_distance * OpenMM::AngstromsPerNm);
                            if (Debug) {
                                qDebug() << "bond start distance = " << solute_bond_perturbation_params[2] <<
                                         " Nm";
                                qDebug() << "bond end distance = " << solute_bond_perturbation_params[3] <<
                                         " Nm";
                                qDebug() << "Perturbation bond equilibrium distance = " << pert_eq_distance <<
                                         " Nm";
                            }
                        }
                        /* JM 10/16 */
                        /*  Here add code to constraint hbonds only if initial and final parameters are unperturbed*/
                        /*  check also what is the mass of the atoms in that case */
                        else if (flag_constraint == HBONDS and flag_noperturbedconstraints) {
                            const SireMol::Atom atom0 = molecule.select(two.atom0());
                            //double m0 = atom0.property("mass").value();
                            double m0 = system_openmm->getParticleMass(idx0);
                            const SireMol::Atom atom1 = molecule.select(two.atom1());
                            //double m1 = atom1.property("mass").value();
                            double m1 = system_openmm->getParticleMass(idx1);
                            double deltar = abs(rend-rstart);
                            double deltak = abs(bend-bstart);
                            // only constraint if m0 < 1.1 g.mol-1 or m1 < 1.1 g.mol-1
                            // AND the initial and final parameters differ
                            if (Debug) {
                                qDebug() << " m0 " << m0 << " m1 " << m1 << "\n";
                                qDebug() << " deltar " << deltar << " " << " deltak " << deltak;
                            }
                            /* bonds that do not change parameters are constrained*/
                            double pert_eq_distance = solute_bond_perturbation_params[3] * current_lambda
                                                      + (1.0 - current_lambda) * solute_bond_perturbation_params[2];
                            if (deltar < SMALL and deltak < SMALL) {
                                system_openmm->addConstraint(idx0, idx1, pert_eq_distance);
                                if (Debug) {
                                    qDebug() << "perturbed bond but no parameter changes so constrained " <<
                                             atom0.name().toString()
                                             << "-" << atom1.name().toString() << "\n";
                                }
                            }
                            /* bonds that change parameters and have one of the atoms with a mass < HMASS are constrained*/
                            else if (m0 < HMASS or m1 < HMASS) {
                                system_openmm->addConstraint(idx0, idx1, pert_eq_distance);
                                if (Debug) {
                                    qDebug() << "perturbed bond parameter changes but involving "
                                             << " light mass so constrained " << atom0.name().toString()
                                             << "- " << atom1.name().toString() << "\n";
                                }
                            }
                            /* other bonds are flexible */
                            else {
                                solute_bond_perturbation->addBond(idx0, idx1, solute_bond_perturbation_params);
                                if (Debug) {
                                    qDebug() << "perturbed bond flexible " << atom0.name().toString()
                                             << "- " << atom1.name().toString() << "\n";
                                }
                            }
                        }
                        else if (flag_constraint == HBONDS) {
                            const SireMol::Atom atom0 = molecule.select(two.atom0());
                            QString initial_type_atom0 = atom0.property<QString>("initial_ambertype");
                            QString final_type_atom0 = atom0.property<QString>("final_ambertype");

                            const SireMol::Atom atom1 = molecule.select(two.atom1());
                            QString initial_type_atom1 = atom1.property<QString>("initial_ambertype");
                            QString final_type_atom1 = atom1.property<QString>("final_ambertype");

                            if (initial_type_atom0.startsWith("h", Qt::CaseInsensitive)
                                    || final_type_atom0.startsWith("h", Qt::CaseInsensitive) ||
                                    initial_type_atom1.startsWith("h", Qt::CaseInsensitive)
                                    || final_type_atom1.startsWith("h", Qt::CaseInsensitive)) {
                                double pert_eq_distance = solute_bond_perturbation_params[3] * current_lambda
                                                          + (1.0 - current_lambda) * solute_bond_perturbation_params[2];
                                system_openmm->addConstraint(idx0, idx1, pert_eq_distance);

                                if (Debug) {
                                    qDebug() << "Two/one bond atom(s) start(s) or end(s) with h/H";
                                    qDebug() << "bond start distance = " << solute_bond_perturbation_params[2] <<
                                             " Nm";
                                    qDebug() << "bond end distance = " << solute_bond_perturbation_params[3] <<
                                             " Nm";
                                    qDebug() << "Perturbation bond equilibrium distance = " << pert_eq_distance <<
                                             " Nm";
                                }
                            }
                            else {
                                solute_bond_perturbation->addBond(idx0, idx1, solute_bond_perturbation_params);
                            }

                            if (Debug) {
                                qDebug() << "Atom0 initil type = " << initial_type_atom0;
                                qDebug() << "Atom0 final type = " << final_type_atom0;
                                qDebug() << "Atom1 initil type = " << initial_type_atom1;
                                qDebug() << "Atom1 final type = " << final_type_atom1;
                            }

                        }

                        bond_pert_list.append(BondID(two.atom0(), two.atom1()));
                        bond_pert_swap_list.append(BondID(two.atom1(), two.atom0()));

                        bondPairs.push_back(std::make_pair(idx0, idx1));

                        if (Debug) {
                            qDebug() << "Atom0 = " << two.atom0().asA<AtomIdx>().value() <<
                                     "Atom1 = " << two.atom1().asA<AtomIdx>().value();
                            qDebug() << "IDX0 = " << idx0 << "IDX1 = " << idx1 << "\n";
                            qDebug() << "rstart = " << rstart << " A" <<
                                     "rend = " << rend << " A";
                            qDebug() << "bstart = " << bstart << " kcal/A A" <<
                                     "bend = " << bend << " kcal/A A" << "\n";
                        }
                    }
                    if (str == "SireMM::ThreeAtomPerturbation") {
                        const ThreeAtomPerturbation &three = pert.asA<ThreeAtomPerturbation>();
                        int idx0 = three.atom0().asA<AtomIdx>().value() + num_atoms_till_i;
                        int idx1 = three.atom1().asA<AtomIdx>().value() + num_atoms_till_i;
                        int idx2 = three.atom2().asA<AtomIdx>().value() + num_atoms_till_i;
                        double astart = three.initialForms()[Symbol("k")].toString().toDouble();
                        double thetastart =
                            three.initialForms()[Symbol("theta0")].toString().toDouble();
                        double aend = three.finalForms()[Symbol("k")].toString().toDouble();
                        double thetaend = three.finalForms()[Symbol("theta0")].toString().toDouble();

                        solute_angle_perturbation_params[0] = astart * 2.0 * OpenMM::KJPerKcal;
                        solute_angle_perturbation_params[1] = aend * 2.0 * OpenMM::KJPerKcal;
                        solute_angle_perturbation_params[2] = thetastart;
                        solute_angle_perturbation_params[3] = thetaend;

                        if (Debug) {
                            qDebug() << "astart = " << solute_angle_perturbation_params[0] << " kJ/rad rad"
                                     <<
                                     " aend = " << solute_angle_perturbation_params[1] << " kJ/rad rad";
                            qDebug() << "thetastart = " << solute_angle_perturbation_params[2] << " rad" <<
                                     "thetaend = " << solute_angle_perturbation_params[3] << " rad";
                        }

                        if (flag_constraint == HANGLES) {
                            const SireMol::Atom atom0 = molecule.select(three.atom0());
                            QString initial_type_atom0 = atom0.property<QString>("initial_ambertype");
                            QString final_type_atom0 = atom0.property<QString>("final_ambertype");

                            const SireMol::Atom atom1 = molecule.select(three.atom1());
                            QString initial_type_atom1 = atom1.property<QString>("initial_ambertype");
                            QString final_type_atom1 = atom1.property<QString>("final_ambertype");

                            const SireMol::Atom atom2 = molecule.select(three.atom2());
                            QString initial_type_atom2 = atom2.property<QString>("initial_ambertype");
                            QString final_type_atom2 = atom2.property<QString>("final_ambertype");

                            bool H_X_H = (initial_type_atom0.startsWith("h", Qt::CaseInsensitive)
                                          || final_type_atom0.startsWith("h", Qt::CaseInsensitive)) &&
                                         (initial_type_atom2.startsWith("h", Qt::CaseInsensitive)
                                          || final_type_atom2.startsWith("h", Qt::CaseInsensitive));

                            bool H_O_X = (initial_type_atom0.startsWith("h", Qt::CaseInsensitive)
                                          || final_type_atom0.startsWith("h", Qt::CaseInsensitive)) &&
                                         (initial_type_atom1.startsWith("o", Qt::CaseInsensitive)
                                          || final_type_atom1.startsWith("o", Qt::CaseInsensitive));

                            bool X_O_H = (initial_type_atom1.startsWith("o", Qt::CaseInsensitive)
                                          || final_type_atom1.startsWith("o", Qt::CaseInsensitive)) &&
                                         (initial_type_atom2.startsWith("h", Qt::CaseInsensitive)
                                          || final_type_atom2.startsWith("h", Qt::CaseInsensitive));

                            if (Debug) {
                                if (H_X_H)
                                    qDebug() << "type =  H_X_H";
                                if (H_O_X)
                                    qDebug() << "type =  H_O_X";
                                if (X_O_H)
                                    qDebug() << "type =  X_O_H";
                            }

                            if (H_X_H || H_O_X || X_O_H) {
                                const BondID * first_alchemical_bond = NULL;
                                const BondID * second_alchemical_bond = NULL;

                                double first_alchemical_distance = -1.0;
                                double second_alchemical_distance = -1.0;

                                if (bond_pert_eq_list.contains(BondID(three.atom0(), three.atom1()))) {
                                    first_alchemical_bond = &((bond_pert_eq_list.find(BondID(three.atom0(),
                                                               three.atom1()))).key());
                                    first_alchemical_distance = bond_pert_eq_list.value(*first_alchemical_bond);
                                    if (Debug)
                                        qDebug() << "Atom0 - Atom1";
                                }
                                else if (bond_pert_eq_list.contains(BondID(three.atom1(), three.atom0()))) {
                                    first_alchemical_bond = &((bond_pert_eq_list.find(BondID(three.atom1(),
                                                               three.atom0()))).key());
                                    first_alchemical_distance = bond_pert_eq_list.value(*first_alchemical_bond);
                                    if (Debug)
                                        qDebug() << "Atom1 - Atom0";
                                }
                                else {
                                    if (Debug)
                                        qDebug() << "First perturbed bond was not foud in the perturned list";
                                    first_alchemical_bond = new BondID(three.atom0(), three.atom1());
                                }


                                if (bond_pert_eq_list.contains(BondID(three.atom1(), three.atom2()))) {
                                    second_alchemical_bond = &((bond_pert_eq_list.find(BondID(three.atom1(),
                                                                three.atom2()))).key());
                                    second_alchemical_distance = bond_pert_eq_list.value(*second_alchemical_bond);
                                    if (Debug)
                                        qDebug() << "Atom1 - Atom2";
                                }
                                else if (bond_pert_eq_list.contains(BondID(three.atom2(), three.atom1()))) {
                                    second_alchemical_bond = &((bond_pert_eq_list.find(BondID(three.atom2(),
                                                                three.atom1()))).key());
                                    second_alchemical_distance = bond_pert_eq_list.value(*second_alchemical_bond);
                                    if (Debug)
                                        qDebug() << "Atom2 - Atom1";
                                }
                                else {
                                    if (Debug)
                                        qDebug() << "Second perturbed bond was not foud in the perturned list";
                                }


                                if (Debug)
                                    qDebug() << "First Alchemical distance = " << first_alchemical_distance
                                             << "Second Alchemical distance = " << second_alchemical_distance;

                                SireMaths::Vector bond1_vec;
                                SireMaths::Vector bond2_vec;

                                if (first_alchemical_bond->atom0() == second_alchemical_bond->atom0()) {
                                    SireMaths::Vector tmp1 = (molecule.atom(
                                                                  first_alchemical_bond->atom1())).property<SireMaths::Vector>("coordinates");
                                    SireMaths::Vector tmp2 = (molecule.atom(
                                                                  second_alchemical_bond->atom1())).property<SireMaths::Vector>("coordinates");
                                    SireMaths::Vector tmpc = (molecule.atom(
                                                                  first_alchemical_bond->atom0())).property<SireMaths::Vector>("coordinates");
                                    bond1_vec = tmp1 - tmpc;
                                    bond2_vec = tmp2 - tmpc;
                                    if (Debug)
                                        qDebug() << "Bond1 Atom0 = Bond2 Atom0";
                                }
                                else if (first_alchemical_bond->atom0() == second_alchemical_bond->atom1()) {
                                    SireMaths::Vector tmp1 = (molecule.atom(
                                                                  first_alchemical_bond->atom1())).property<SireMaths::Vector>("coordinates");
                                    SireMaths::Vector tmp2 = (molecule.atom(
                                                                  second_alchemical_bond->atom0())).property<SireMaths::Vector>("coordinates");
                                    SireMaths::Vector tmpc = (molecule.atom(
                                                                  first_alchemical_bond->atom0())).property<SireMaths::Vector>("coordinates");
                                    bond1_vec = tmp1 - tmpc;
                                    bond2_vec = tmp2 - tmpc;
                                    if (Debug)
                                        qDebug() << "Bond1 Atom0 = Bond2 Atom1";
                                }
                                else if (first_alchemical_bond->atom1() == second_alchemical_bond->atom0()) {
                                    SireMaths::Vector tmp1 = (molecule.atom(
                                                                  first_alchemical_bond->atom0())).property<SireMaths::Vector>("coordinates");
                                    SireMaths::Vector tmp2 = (molecule.atom(
                                                                  second_alchemical_bond->atom1())).property<SireMaths::Vector>("coordinates");
                                    SireMaths::Vector tmpc = (molecule.atom(
                                                                  first_alchemical_bond->atom1())).property<SireMaths::Vector>("coordinates");
                                    bond1_vec = tmp1 - tmpc;
                                    bond2_vec = tmp2 - tmpc;
                                    if (Debug)
                                        qDebug() << "Bond1 Atom1 = Bond2 Atom0";
                                }
                                else if (first_alchemical_bond->atom1() == second_alchemical_bond->atom1()) {
                                    SireMaths::Vector tmp1 = (molecule.atom(
                                                                  first_alchemical_bond->atom0())).property<SireMaths::Vector>("coordinates");
                                    SireMaths::Vector tmp2 = (molecule.atom(
                                                                  second_alchemical_bond->atom0())).property<SireMaths::Vector>("coordinates");
                                    SireMaths::Vector tmpc = (molecule.atom(
                                                                  first_alchemical_bond->atom1())).property<SireMaths::Vector>("coordinates");
                                    bond1_vec = tmp1 - tmpc;
                                    bond2_vec = tmp2 - tmpc;
                                    if (Debug)
                                        qDebug() << "Bond1 Atom1 = Bond2 Atom1";
                                }
                                else
                                    throw SireError::program_bug(QObject::tr("No coorner bond"), CODELOC);

                                if (Debug) {
                                    if (first_alchemical_distance != -1.0) {
                                        qDebug() << "First vector X = " << (bond1_vec.normalise() *
                                                                            first_alchemical_distance).x();
                                        qDebug() << "First vector Y = " << (bond1_vec.normalise() *
                                                                            first_alchemical_distance).y();
                                        qDebug() << "First vector Z = " << (bond1_vec.normalise() *
                                                                            first_alchemical_distance).z();
                                    }
                                    else {
                                        qDebug() << "First vector X = " << (bond1_vec).x();
                                        qDebug() << "First vector Y = " << (bond1_vec).y();
                                        qDebug() << "First vector Z = " << (bond1_vec).z();
                                    }

                                    if (second_alchemical_distance != -1.0) {
                                        qDebug() << "Second vector X = " << (bond2_vec.normalise() *
                                                                             second_alchemical_distance).x();
                                        qDebug() << "Second vector Y = " << (bond2_vec.normalise() *
                                                                             second_alchemical_distance).y();
                                        qDebug() << "Second vector Z = " << (bond2_vec.normalise() *
                                                                             second_alchemical_distance).z();
                                    }
                                    else {
                                        qDebug() << "Second vector X = " << (bond2_vec).x();
                                        qDebug() << "Second vector Y = " << (bond2_vec).y();
                                        qDebug() << "Second vector Z = " << (bond2_vec).z();
                                    }
                                }

                                double constraint_distance;

                                double eq_angle = solute_angle_perturbation_params[3] * current_lambda +
                                                  (1.0 - current_lambda) * solute_angle_perturbation_params[2];

                                if (first_alchemical_distance == -1.0 && second_alchemical_distance != -1.0) {
                                    //Carnot's theorem a^2 = c^2 + b^2 - a*b*c*cos(bc)
                                    double sq = bond1_vec.length() * bond1_vec.length() +
                                                (bond2_vec.normalise() * second_alchemical_distance).length() *
                                                (bond2_vec.normalise() * second_alchemical_distance).length();

                                    double dp = 2.0 * bond1_vec.length() * (bond2_vec.normalise() *
                                                                            second_alchemical_distance).length() * cos(eq_angle);

                                    constraint_distance = sqrt(sq - dp);

                                    //constraint_distance = (bond1_vec - bond2_vec.normalise() * second_alchemical_distance).length();
                                }
                                else if (first_alchemical_distance != -1.0
                                         && second_alchemical_distance == -1.0) {
                                    //Carnot theorem a^2 = c^2 + b^2 - a*b*c*cos(bc)
                                    double sq = bond2_vec.length() * bond2_vec.length() +
                                                (bond1_vec.normalise() * first_alchemical_distance).length() *
                                                (bond1_vec.normalise() * first_alchemical_distance).length();

                                    double dp = 2.0 * bond2_vec.length() * (bond2_vec.normalise() *
                                                                            first_alchemical_distance).length() * cos(eq_angle);

                                    constraint_distance = sqrt(sq - dp);

                                    //constraint_distance = (bond1_vec.normalise() * first_alchemical_distance - bond2_vec).length();
                                }
                                else if (first_alchemical_distance != -1.0
                                         && second_alchemical_distance != -1.0) {
                                    //Carnot's theorem a^2 = c^2 + b^2 - a*b*c*cos(bc)
                                    double sq = (bond1_vec.normalise() * first_alchemical_distance).length() *
                                                (bond1_vec.normalise() * first_alchemical_distance).length() +
                                                (bond2_vec.normalise() * second_alchemical_distance).length() *
                                                (bond2_vec.normalise() * second_alchemical_distance).length();

                                    double dp = 2.0 * (bond1_vec.normalise() * first_alchemical_distance).length()
                                                * (bond2_vec.normalise() * second_alchemical_distance).length() * cos(eq_angle);

                                    constraint_distance = sqrt(sq - dp);

                                    //constraint_distance = (bond1_vec.normalise() * first_alchemical_distance - bond2_vec.normalise() * second_alchemical_distance).length();
                                }
                                else
                                    throw SireError::program_bug(
                                        QObject::tr("The angle does not contain perturbed bond"), CODELOC);


                                system_openmm->addConstraint(idx0, idx2,
                                                             constraint_distance * OpenMM::NmPerAngstrom);

                                if (Debug)
                                    qDebug() << "CONSTRAINT DISTANCE = " << constraint_distance << " A";
                            }

                        }//end if HANGLES
                        else {
                            solute_angle_perturbation->addAngle(idx0, idx1, idx2,
                                                                solute_angle_perturbation_params);
                            if (Debug)
                                qDebug() << "Added perturbed angle";
                        }

                        angle_pert_list.append(AngleID(three.atom0(), three.atom1(), three.atom2()));
                        angle_pert_swap_list.append(AngleID(three.atom2(), three.atom1(),
                                                            three.atom0()));

                        if (Debug) {
                            qDebug() << "Atom0 = " << three.atom0().asA<AtomIdx>().value() <<
                                     "Atom1 = " << three.atom1().asA<AtomIdx>().value() <<
                                     "Atom2 = " << three.atom2().asA<AtomIdx>().value();

                            qDebug() << "IDX0 = " << idx0 << "IDX1 = " << idx1 << "IDX2 = " << idx2 << "\n";

                            qDebug() << "astart = " << astart << " kcal/rad rad" <<
                                     "aend = " << " kcal/rad rad" << aend;

                            qDebug() << "thetastart = " << thetastart << " rad" <<
                                     "thetaend = " << thetaend << " rad" << "\n";
                        }
                    }
                    if (str == "SireMM::FourAtomPerturbation") {

                        const FourAtomPerturbation &four = pert.asA<FourAtomPerturbation>();
                        int idx0 = four.atom0().asA<AtomIdx>().value() + num_atoms_till_i;
                        int idx1 = four.atom1().asA<AtomIdx>().value() + num_atoms_till_i;
                        int idx2 = four.atom2().asA<AtomIdx>().value() + num_atoms_till_i;
                        int idx3 = four.atom3().asA<AtomIdx>().value() + num_atoms_till_i;

                        QString tmp = four.perturbExpression().toOpenMMString();
                        tmp.replace(QString("phi"), QString("theta"));
                        tmp.replace(QString("lambda"), QString("lamdih"));
                        tmp = "( " + tmp + " ) * " + "KJPerKcal";

                        std::string openmm_str = tmp.toStdString();

                        if (Debug) {
                            qDebug() << "IDX0 = " << idx0 << "IDX1 = " << idx1 << "IDX2 = " << idx2 <<
                                     "IDX3 = " << idx3;
                            qDebug() << "Dihedral String = " << openmm_str.c_str();
                            qDebug() << "Dihedral Normal String = " << four.perturbExpression().toString()
                                     << "\n";
                        }

                        auto solute_torsion_perturbation = new OpenMM::CustomTorsionForce(openmm_str);
                        solute_torsion_perturbation->addPerTorsionParameter("KJPerKcal");
                        solute_torsion_perturbation_params[0] = 4.184;
                        solute_torsion_perturbation->addGlobalParameter("lamdih", current_lambda);
                        solute_torsion_perturbation->addTorsion(idx0, idx1, idx2, idx3,
                                                                solute_torsion_perturbation_params);

                        //********************************BONDED ENERGY TORSIONS ARE ADDED TO THE SYSTEM*****************************
                        solute_torsion_perturbation->setForceGroup(
                            NONBONDED_FCG); // FIXME: why in this force group?

                        if (!fullPME)
                            system_openmm->addForce(solute_torsion_perturbation);

                        perturbed_energies_tmp[7] = true; //Torsions are added to the system

                        dihedral_pert_list.append(DihedralID(four.atom0(), four.atom1(), four.atom2(),
                                                             four.atom3()));
                        dihedral_pert_swap_list.append(DihedralID(four.atom3(), four.atom1(),
                                                       four.atom2(), four.atom0()));

                        improper_pert_list.append(ImproperID(four.atom0(), four.atom1(), four.atom2(),
                                                             four.atom3()));
                        improper_pert_swap_list.append(ImproperID(four.atom0(), four.atom1(),
                                                       four.atom3(), four.atom2()));

                        if (Debug) {
                            qDebug() << "Atom0 = " << four.atom0().asA<AtomIdx>().value() <<
                                     "Atom1 = " << four.atom1().asA<AtomIdx>().value() <<
                                     "Atom2 = " << four.atom2().asA<AtomIdx>().value() <<
                                     "Atom3 = " << four.atom3().asA<AtomIdx>().value() << "\n";
                        }
                    }

                }
            } // end for perturbations

        } // end solute molecule perturbation

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
            double k = bond_params[0];
            double r0 = bond_params[1];

            int idx0 = bonds[j].atom0().asA<AtomIdx>().value();
            int idx1 = bonds[j].atom1().asA<AtomIdx>().value();

            if (solute.contains(molecule)) {
                if (bond_pert_list.indexOf(bond_ff) != -1
                        || bond_pert_swap_list.indexOf(bond_ff) != -1) {
                    //Solute molecule. Check if the current solute bond is in the perturbed bond list
                    // JM July 13 --> Note, we should still have the ability to constrain the bond to its r(lambda) equilibrium distance
                    if (Debug)
                        qDebug() << "Found Perturbed Bond\n";
                    continue;
                }
            }

            //Select the atom type
            QString atom0 = molecule.atom(AtomIdx(idx0)).toString();
            QString atom1 = molecule.atom(AtomIdx(idx1)).toString();

            idx0 = idx0 + num_atoms_till_i;
            idx1 = idx1 + num_atoms_till_i;

            // JM July 13. The constraint code has to be revised to handle massless particles so TIP4P can be dealt with.
            // Should check if a constraint involves an atom with a null mass, and if so skip constraint.

            if (flag_constraint == NONE) {
                //JM 10/16 If constraint water flag is on and if molecule is a water molecule then apply constraint
                if (flag_constraint_water and molfirstresname == ResName("WAT"))
                    system_openmm->addConstraint(idx0, idx1, r0 * OpenMM::NmPerAngstrom);
                else
                    bondStretch_openmm->addBond(idx0, idx1, r0 * OpenMM::NmPerAngstrom,
                                                k * 2.0 * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
            }
            else if (flag_constraint == ALLBONDS || flag_constraint == HANGLES) {
                system_openmm->addConstraint(idx0, idx1, r0 * OpenMM::NmPerAngstrom);
            }
            else if (flag_constraint == HBONDS) {
                if ((atom0[6] == 'H') || (atom1[6] == 'H')) {
                    system_openmm->addConstraint(idx0, idx1, r0 * OpenMM::NmPerAngstrom);
                }
                else {
                    bondStretch_openmm->addBond(idx0, idx1, r0 * OpenMM::NmPerAngstrom,
                                                k * 2.0 * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
                }
            }

            // Bond exclusion List
            bondPairs.push_back(std::make_pair(idx0, idx1));
        }

        // Angles

        QList<AngleID> angles_ff = amber_params.getAllAngles();
        QVector<AngleID> angles = angles_ff.toVector();

        for (int j = 0; j < angles_ff.length(); j++) {
            AngleID angle_ff = angles_ff[j];
            QList<double> angle_params = amber_params.getParams(angle_ff);

            double k = angle_params[0];
            double theta0 = angle_params[1]; // It is already in radiant

            int idx0 = angles[j].atom0().asA<AtomIdx>().value();
            int idx1 = angles[j].atom1().asA<AtomIdx>().value();
            int idx2 = angles[j].atom2().asA<AtomIdx>().value();

            if (solute.contains(molecule)) {
                if (angle_pert_list.indexOf(angle_ff) != -1
                        || angle_pert_swap_list.indexOf(angle_ff) != -1) {
                    //Solute molecule. Check if the current solute angle is in the perturbed angle list
                    if (Debug)
                        qDebug() << "Found Perturbed Angle\n";
                    continue;
                }
                else {
                    if (Debug)
                        qDebug() << "Solute normal Angle - Atom0 = " << idx0 << "Atom1 = " << idx1 <<
                                 "Atom2 = " << idx2 << "theta0 = " << theta0 << " k = " << k << "\n";

                    idx0 = idx0 + num_atoms_till_i;
                    idx1 = idx1 + num_atoms_till_i;
                    idx2 = idx2 + num_atoms_till_i;
                    bondBend_openmm->addAngle(idx0, idx1, idx2, theta0,
                                              k * 2.0 * OpenMM::KJPerKcal);
                    continue;
                }
            }
            if (Debug)
                qDebug() << "Angle - Atom0 = " << idx0 << "Atom1 = " << idx1 << "Atom2 = " <<
                         idx2 << "\n";

            QString atom0 = molecule.atom(AtomIdx(idx0)).toString();
            QString atom1 = molecule.atom(AtomIdx(idx1)).toString();
            QString atom2 = molecule.atom(AtomIdx(idx2)).toString();

            Vector diff = c[idx2] - c[idx0];

            idx0 = idx0 + num_atoms_till_i;
            idx1 = idx1 + num_atoms_till_i;
            idx2 = idx2 + num_atoms_till_i;

            if (flag_constraint == HANGLES) {
                if (((atom0[6] == 'H') && (atom2[6] == 'H'))) {
                    system_openmm->addConstraint(idx0, idx2, diff.length() * OpenMM::NmPerAngstrom);
                }
                else if (((atom0[6] == 'H') && (atom1[6] == 'O')) || ((atom1[6] == 'O')
                         && (atom2[6] == 'H'))) {
                    system_openmm->addConstraint(idx0, idx2, diff.length() * OpenMM::NmPerAngstrom);
                }
                else {
                    bondBend_openmm->addAngle(idx0, idx1, idx2, theta0,
                                              k * 2.0 * OpenMM::KJPerKcal);
                }
            }
            else {
                bondBend_openmm->addAngle(idx0, idx1, idx2, theta0,
                                          k * 2.0 * OpenMM::KJPerKcal);
            }
        } // end of angles

        // Dihedrals

        QList<DihedralID> dihedrals_ff = amber_params.getAllDihedrals();
        QVector<DihedralID> dihedrals = dihedrals_ff.toVector();

        for (int j = 0; j < dihedrals_ff.length(); j++) {
            DihedralID dihedral_ff = dihedrals_ff[j];
            QList<double> dihedral_params = amber_params.getParams(dihedral_ff);

            int idx0 = dihedrals[j].atom0().asA<AtomIdx>().value() + num_atoms_till_i;
            int idx1 = dihedrals[j].atom1().asA<AtomIdx>().value() + num_atoms_till_i;
            int idx2 = dihedrals[j].atom2().asA<AtomIdx>().value() + num_atoms_till_i;
            int idx3 = dihedrals[j].atom3().asA<AtomIdx>().value() + num_atoms_till_i;

            if (Debug) {
                qDebug() << "TOTAL Dihedral between atom global index " << idx0 -
                         num_atoms_till_i <<
                         " and " << idx1 - num_atoms_till_i <<
                         " and " << idx2 - num_atoms_till_i <<
                         " and " << idx3 - num_atoms_till_i << "\n";
            }

            if (solute.contains(molecule)) {
                if (dihedral_pert_list.indexOf(dihedral_ff) != -1
                        || dihedral_pert_swap_list.indexOf(dihedral_ff) != -1) {
                    //Solute molecule. Check if the current solute dihedral is in the perturbed dihedral list
                    if (Debug)
                        qDebug() << "Found Perturbed Dihedral\n";
                    continue;
                }
            }

            // Variable number of parameters
            for (int k = 0; k < dihedral_params.length(); k = k + 3) {
                double v = dihedral_params[ k ];
                int periodicity = dihedral_params[ k + 1 ];
                double phase = dihedral_params[ k + 2 ];
                bondTorsion_openmm->addTorsion(idx0, idx1, idx2, idx3, periodicity, phase,
                                               v * OpenMM::KJPerKcal);
                if (Debug) {
                    qDebug() << "Dihedral between atom global index " << idx0 - num_atoms_till_i <<
                             " and " << idx1 - num_atoms_till_i <<
                             " and " << idx2 - num_atoms_till_i <<
                             " and " << idx3 - num_atoms_till_i << "\n";
                    qDebug() << "Amplitude_dih = " << v << " periodicity " << periodicity <<
                             " phase " << phase << "\n";
                }
            }
        } // end of dihedrals

        // Improper Dihedrals

        QList<ImproperID> impropers_ff = amber_params.getAllImpropers();
        QVector<ImproperID> impropers = impropers_ff.toVector();

        for (int j = 0; j < impropers_ff.length(); j++) {
            ImproperID improper_ff = impropers_ff[j];
            QList<double> improper_params = amber_params.getParams(improper_ff);

            int idx0 = impropers[j].atom0().asA<AtomIdx>().value() + num_atoms_till_i;
            int idx1 = impropers[j].atom1().asA<AtomIdx>().value() + num_atoms_till_i;
            int idx2 = impropers[j].atom2().asA<AtomIdx>().value() + num_atoms_till_i;
            int idx3 = impropers[j].atom3().asA<AtomIdx>().value() + num_atoms_till_i;

            if (Debug) {
                qDebug() << "TOTAL Improper between atom global index " << idx0 -
                         num_atoms_till_i <<
                         " and " << idx1 - num_atoms_till_i <<
                         " and " << idx2 - num_atoms_till_i <<
                         " and " << idx3 - num_atoms_till_i << "\n";
            }

            if (solute.contains(molecule)) {
                //Solute molecule. Check if the current solute dihedral is in the perturbed improper list
                if (improper_pert_list.indexOf(improper_ff) != -1
                        || improper_pert_swap_list.indexOf(improper_ff) != -1) {
                    if (Debug)
                        qDebug() << "Found Perturbed Improper\n";
                    continue;
                }
            }

            // Variable number of parameters
            for (int k = 0; k < improper_params.length(); k = k + 3) {
                double v = improper_params[ k ];
                int periodicity = improper_params[ k + 1 ];
                double phase = improper_params[ k + 2 ];

                bondTorsion_openmm->addTorsion(idx0, idx1, idx2, idx3, periodicity, phase,
                                               v * OpenMM::KJPerKcal);
                if (Debug) {
                    qDebug() << "Improper between atom global index " << idx0 - num_atoms_till_i <<
                             " and " << idx1 - num_atoms_till_i <<
                             " and " << idx2 - num_atoms_till_i <<
                             " and " << idx3 - num_atoms_till_i << "\n";
                    qDebug() << "Amplitude_imp = " << v << " periodicity " << periodicity <<
                             " phase " << phase << "\n";
                }
            }
        } // end of impropers

        // Variable 1,4 scaling factors
        QList<BondID> pairs14_ff = amber_params.getAll14Pairs();
        QVector<BondID> pairs14 = pairs14_ff.toVector();

        for (int j = 0; j < pairs14_ff.length(); j++) {
            BondID pair14_ff = pairs14_ff[j];

            QList<double> pair14_params = amber_params.get14PairParams(pair14_ff);

            double cscl = pair14_params[0];
            double ljscl = pair14_params[1];

            if (Debug)
                qDebug() << " cscl@ " << cscl << " ljscl " << ljscl;

            // Add to custom pairs if scale factor differs from default
            if (abs(cscl - Coulomb14Scale) > 0.0001
                    or abs(ljscl - LennardJones14Scale) > 0.0001) {
                int idx0 = pair14_ff.atom0().asA<AtomIdx>().value() + num_atoms_till_i;
                int idx1 = pair14_ff.atom1().asA<AtomIdx>().value() + num_atoms_till_i;

                QPair<int, int> indices_pair(idx0, idx1);
                QPair<double, double> scl_pair(cscl, ljscl);
                custom14pairs.insert(indices_pair, scl_pair);

                special_14 = true;

                if (Debug)
                    qDebug() << "IDX0 = " << idx0 << " IDX1 =" << idx1 << "14 OpenMM Index";
            }
        } // end of variable 1,4 scaling factors

        num_atoms_till_i = num_atoms_till_i + num_atoms_molecule;

    } // end of loop over molecules

    if (Debug) {
        if (nions != 0)
            qDebug() << "\nNumber of ions = " << nions << "\n";
    }


    /*** EXCEPTION HANDLING ***/

    // Exclude the 1-2, 1-3 bonded atoms from nonbonded forces, and scale
    // down 1-4 bonded atoms
    recip_space->createExceptionsFromBonds(bondPairs, Coulomb14Scale,
                                           LennardJones14Scale);

    int num_exceptions = recip_space->getNumExceptions();

    if (Debug)
        qDebug() << "num exceptions =" << num_exceptions << "\n";

    for (int i = 0; i < num_exceptions; i++) {
        int p1, p2;

	double qprod_diff, qprod_start, qprod_end, qprod_mix;
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
        qprod_mix = Qend_p1 * Qstart_p2 + Qstart_p1 * Qend_p2;

        if (Debug)
            qDebug() << "Exception =" << i << ", p1 =" << p1 << ", p2 ="
                     << p2 << ", charge prod =" << charge_prod
                     << ", sigma avg =" << sigma_avg << ", epsilon_avg ="
                     << epsilon_avg;

        // run over the 1-4 exceptions
        if (!(charge_prod == 0 && sigma_avg == 1 && epsilon_avg == 0)) {
            QVector<double> perturbed_14_tmp(13);

            double Epstart_p1 = p1_params[2];
            double Epend_p1 = p1_params[3];
            double Sigstart_p1 = p1_params[4];
            double Sigend_p1 = p1_params[5];
            double isHard_p1 = p1_params[6];
            double isTodummy_p1 = p1_params[7];
            double isFromdummy_p1 = p1_params[8];

            double Epstart_p2 = p2_params[2];
            double Epend_p2 = p2_params[3];
            double Sigstart_p2 = p2_params[4];
            double Sigend_p2 = p2_params[5];
            double isHard_p2 = p2_params[6];
            double isTodummy_p2 = p2_params[7];
            double isFromdummy_p2 = p2_params[8];

            double sigma_avg_start, sigma_avg_end, sigma_avg_mix;
            double epsilon_avg_start, epsilon_avg_end, epsilon_avg_mix;

            double Coulomb14Scale_tmp = Coulomb14Scale;
            double LennardJones14Scale_tmp = LennardJones14Scale;

            if (special_14) {
                QPair<double, double> sc_factors;

                QPair<int, int> indices_pair(p1, p2);
                QHash< QPair<int, int>, QPair<double, double> >::const_iterator i_pair =
                    custom14pairs.find(indices_pair);

                if (i_pair != custom14pairs.end()) {
                    sc_factors = i_pair.value();
                    Coulomb14Scale_tmp = sc_factors.first;
                    LennardJones14Scale_tmp = sc_factors.second;

                    if (Debug)
                        qDebug() << "The pair (" << p1 << ", " << p2
				 << ") is 1-4 special no swap pair";
                }
                else {
                    QPair<int, int> indices_swap_pair(p2, p1);
                    QHash< QPair<int, int>, QPair<double, double> >::const_iterator i_swap_pair =
                        custom14pairs.find(indices_swap_pair);

                    if (i_swap_pair != custom14pairs.end()) {
                        sc_factors = i_swap_pair.value();
                        Coulomb14Scale_tmp = sc_factors.first;
                        LennardJones14Scale_tmp = sc_factors.second;

                        if (Debug)
                            qDebug() << "The pair (" << p2 << ", " << p1
				     << ") is 1-4 special swap pair";
                    }
                }
            }

            qprod_start *= Coulomb14Scale_tmp;
            qprod_end *= Coulomb14Scale_tmp;
            qprod_mix *= Coulomb14Scale_tmp;

            if (flag_combRules == ARITHMETIC) {
                sigma_avg_start = (Sigstart_p1 + Sigstart_p2) / 2.0;
                sigma_avg_end = (Sigend_p1 + Sigend_p2) / 2.0;
                sigma_avg_mix = (Sigend_p1 * Sigstart_p2 + Sigstart_p1 * Sigend_p2) / 2.0;
            }
            else if (flag_combRules == GEOMETRIC) {
                sigma_avg_start = Sigstart_p1 * Sigstart_p2 ;
                sigma_avg_end = Sigend_p1 * Sigend_p2 ;
                sigma_avg_mix = Sigend_p1 * Sigstart_p2 + Sigstart_p1 * Sigend_p2 ;
            }

            epsilon_avg_start = Epstart_p1 * Epstart_p2 * LennardJones14Scale_tmp *
                                LennardJones14Scale_tmp;
            epsilon_avg_end = Epend_p1 * Epend_p2 * LennardJones14Scale_tmp *
                              LennardJones14Scale_tmp;
            epsilon_avg_mix = (Epend_p1 * Epstart_p2 + Epstart_p1 * Epend_p2) *
                              LennardJones14Scale_tmp * LennardJones14Scale_tmp;

            // ["qpstart", "qpend", "qmix", "eastart", "eaend", "emix", "sastart", "saend", "samix"]
            // see setPerParicleParameters and expressions above
            std::vector<double> params(9);

            params[0] = qprod_start;
            params[1] = qprod_end;
            params[2] = qprod_mix;
            params[3] = epsilon_avg_start;
            params[4] = epsilon_avg_end;
            params[5] = epsilon_avg_mix;
            params[6] = sigma_avg_start;
            params[7] = sigma_avg_end;
            params[8] = sigma_avg_mix;

            if (Debug) {
                qDebug() << "Particle p1 = " << p1 << "\nQstart = " << Qstart_p1 << "\nQend = "
                         << Qend_p1
                         << "\nEpstart = " << Epstart_p1 << "\nEpend = " << Epend_p1
                         << "\nSgstart = " << Sigstart_p1 << "\nSgend = " << Sigend_p1
                         << "\nisHard = " << isHard_p1 << "\nisTodummy = " << isTodummy_p1 <<
                         "\nisFromdummy = " << isFromdummy_p1 << "\n";
                qDebug() << "Particle p2 = " << p2 << "\nQstart = " << Qstart_p2 << "\nQend = "
                         << Qend_p2
                         << "\nEpstart = " << Epstart_p2 << "\nEpend = " << Epend_p2
                         << "\nSgstart = " << Sigstart_p2 << "\nSgend = " << Sigend_p2
                         << "\nisHard = " << isHard_p2 << "\nisTodummy = " << isTodummy_p2 <<
                         "\nisFromdummy = " << isFromdummy_p2 << "\n";

                qDebug() << "Product Charge start = " << qprod_start <<
                         "\nProduct Charge end = " << qprod_end << "\nProduct Charge mixed = " <<
                         qprod_mix
                         << "\nEpsilon average start = " << epsilon_avg_start <<
                         "\nEpsilon average end = " << epsilon_avg_end << "\nEpsilon average mixed = " <<
                         qprod_mix
                         << "\nSigma average start = " << sigma_avg_start << "\nSigma average end = " <<
                         sigma_avg_end;
                qDebug() << "Coulombic Scale Factor = " << Coulomb14Scale_tmp <<
                         " Lennard-Jones Scale Factor = " << LennardJones14Scale_tmp << "\n";
            }

            if ((isHard_p1 == 1.0 && isHard_p2 == 1.0)) {
                custom_intra_14_clj->addBond(p1, p2, params);

                if (Debug)
                    qDebug() << "Added clj Hard 1-4 exception";
            }
            else if ((isTodummy_p1 == 1.0 && isTodummy_p2 == 1.0) || (isHard_p1 == 1.0
                     && isTodummy_p2 == 1.0) || (isHard_p2 == 1.0 && isTodummy_p1 == 1.0)) {
                custom_intra_14_todummy->addBond(p1, p2, params);

                if (Debug)
                    qDebug() << "Added soft TO dummy 1-4 exception";
            }

            else if ((isFromdummy_p1 == 1.0 && isFromdummy_p2 == 1.0) || (isHard_p1 == 1.0
                     && isFromdummy_p2 == 1.0) || (isHard_p2 == 1.0 && isFromdummy_p1 == 1.0)) {
                custom_intra_14_fromdummy->addBond(p1, p2, params);

                if (Debug)
                    qDebug() << "Added soft FROM dummy 1-4 exception";
            }

            else if ((isFromdummy_p1 == 1.0 && isTodummy_p2 == 1.0)
                     || (isFromdummy_p2 == 1.0 && isTodummy_p1 == 1.0)) {
                custom_intra_14_fromdummy_todummy->addBond(p1, p2, params);

                if (Debug)
                    qDebug() << "Added soft FROM dummy TO dummy 1-4 exception";
            }
        } // 1-4 exceptions

        qprod_diff = qprod_end - qprod_start;

        if (useOffset && qprod_diff != 0.0) {
            recip_space->addExceptionParameterOffset("lambda_offset", i,
                    qprod_diff, 0.0, 0.0);

            if (Debug)
                qDebug() << "Adding exception offset for atom idx" << i
                         << "; qprod_diff =" << qprod_diff;

        }

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

    if (Debug) {
        qDebug() << "Num pairs  = " << npairs;
        qDebug() << "Num bonds 1-4 Hard = " << custom_intra_14_clj->getNumBonds();
        qDebug() << "Num bonds 1-4 To Dummy = " <<
                 custom_intra_14_todummy->getNumBonds();
        qDebug() << "Num bonds 1-4 From Dummy = " <<
                 custom_intra_14_fromdummy->getNumBonds();
        qDebug() << "Num bonds 1-4 From Dummy To Dummy = " <<
                 custom_intra_14_fromdummy_todummy->getNumBonds();
    }

    system_openmm->addForce(recip_space);
    recip_space->setForceGroup(RECIP_FCG);

    if (!fullPME) {
        if (npairs != num_exceptions) {
            direct_space->setForceGroup(DIRECT_FCG);
            system_openmm->addForce(direct_space);
            perturbed_energies_tmp[0] = true;

            if (Debug)
                qDebug() << "Added 1-5 direct space (PME, LJ):" << general_ff;
        }

        if (custom_corr_recip->getNumBonds() != 0) {
            custom_corr_recip->setForceGroup(CORR_FCG);
            system_openmm->addForce(custom_corr_recip);
            perturbed_energies_tmp[8] = true;

            if (Debug)
                qDebug() << "Added reciprocal correction term:" << CORR_RECIP;
        }

        if (custom_intra_14_clj->getNumBonds() != 0) {
            custom_intra_14_clj->setForceGroup(NONBONDED_FCG);
            system_openmm->addForce(custom_intra_14_clj);
            perturbed_energies_tmp[1] = true;

	    if (Debug)
                qDebug() << "Added 1-4 CLJ:" << intra_14_clj;
        }

        if (custom_intra_14_todummy->getNumBonds() != 0) {
            custom_intra_14_todummy->setForceGroup(NONBONDED_FCG);
            system_openmm->addForce(custom_intra_14_todummy);
            perturbed_energies_tmp[2] = true;

	    if (Debug)
                qDebug() << "Added 1-4 To Dummy:" << intra_14_todummy;
        }

        if (custom_intra_14_fromdummy->getNumBonds() != 0) {
            custom_intra_14_fromdummy->setForceGroup(NONBONDED_FCG);
            system_openmm->addForce(custom_intra_14_fromdummy);
            perturbed_energies_tmp[3] = true;

	    if (Debug)
                qDebug() << "Added 1-4 From Dummy:" << intra_14_fromdummy;
        }

        if (custom_intra_14_fromdummy_todummy->getNumBonds() != 0) {
            custom_intra_14_fromdummy_todummy->setForceGroup(NONBONDED_FCG);
            system_openmm->addForce(custom_intra_14_fromdummy_todummy);
            perturbed_energies_tmp[4] = true;

	    if (Debug)
                qDebug() << "Added 1-4 From Dummy To Dummy:"
			 << intra_14_fromdummy_todummy;
        }
    } // if (!fullPME)


    /*** add bonded force fields to System ***/

    if (bondStretch_openmm->getNumBonds() != 0) {
        bondStretch_openmm->setForceGroup(BOND_FCG);
        system_openmm->addForce(bondStretch_openmm);
        if (Debug)
            qDebug() << "Added Internal Bond energy term";
    }

    if (bondBend_openmm->getNumAngles() != 0) {
        bondBend_openmm->setForceGroup(BOND_FCG);
        system_openmm->addForce(bondBend_openmm);
        if (Debug)
            qDebug() << "Added Internal Angle energy term";
    }

    if (bondTorsion_openmm->getNumTorsions() != 0) {
        bondTorsion_openmm->setForceGroup(BOND_FCG);
        system_openmm->addForce(bondTorsion_openmm);
        if (Debug)
            qDebug() << "Added Internal Torsion energy term";
    }

    if (!fullPME) {
        if (solute_bond_perturbation->getNumBonds() != 0) {
            solute_bond_perturbation->setForceGroup(BOND_FCG);
            system_openmm->addForce(solute_bond_perturbation);
            perturbed_energies_tmp[5] = true; //Custom bonded is added to the system
            if (Debug)
                qDebug() << "Added Perturbed Internal Bond energy term";
        }

        if (solute_angle_perturbation->getNumAngles() != 0) {
            solute_angle_perturbation->setForceGroup(BOND_FCG);
            system_openmm->addForce(solute_angle_perturbation);
            perturbed_energies_tmp[6] = true; //Custom bonded is added to the system
            if (Debug)
                qDebug() << "Added Perturbed Internal Angle energy term";
        }
    } // if (!fullPME)

    perturbed_energies = perturbed_energies_tmp;

    // IMPORTANT: PERTURBED ENERGY TORSIONS ARE ADDED ABOVE
    bool UseLink_flag = true;

    // Distance Restraint. All the information are stored in the first molecule only.

    if (UseLink_flag) {
        Molecule molecule = moleculegroup.moleculeAt(0).molecule();

        bool haslinkinfo = molecule.hasProperty("linkbonds");

        if (haslinkinfo) {
            std::vector<double> custom_bond_link_par(3);

            Properties linkprop = molecule.property("linkbonds").asA<Properties>();

            int nlinks = linkprop.property(
                             QString("nbondlinks")).asA<VariantProperty>().toInt();

            if (Debug)
                qDebug() << "Number of constraint links = " << nlinks;

            for (int i = 0; i < nlinks; i++) {
                int atomnum0 = linkprop.property(QString("AtomNum0(%1)").arg(
                                                     i)).asA<VariantProperty>().toInt();
                int atomnum1 = linkprop.property(QString("AtomNum1(%1)").arg(
                                                     i)).asA<VariantProperty>().toInt();
                double reql = linkprop.property(QString("reql(%1)").arg(
                                                    i)).asA<VariantProperty>().toDouble();
                double kl = linkprop.property(QString("kl(%1)").arg(
                                                  i)).asA<VariantProperty>().toDouble();
                double dl = linkprop.property(QString("dl(%1)").arg(
                                                  i)).asA<VariantProperty>().toDouble();

                int openmmindex0 = AtomNumToOpenMMIndex[atomnum0];
                int openmmindex1 = AtomNumToOpenMMIndex[atomnum1];

                custom_bond_link_par[0] = reql * OpenMM::NmPerAngstrom; //req
                custom_bond_link_par[1] = kl * (OpenMM::KJPerKcal * OpenMM::AngstromsPerNm *
                                                OpenMM::AngstromsPerNm); //k
                custom_bond_link_par[2] = dl * OpenMM::NmPerAngstrom; //dl

                if (Debug) {
                    qDebug() << "atomnum0 = " << atomnum0 << " openmmindex0 =" << openmmindex0;
                    qDebug() << "atomnum1 = " << atomnum1 << " openmmindex1 =" << openmmindex1;
                    qDebug() << "Req = " << reql << " kl = " << kl << " dl = " << dl;
                }

                custom_link_bond->addBond(openmmindex0, openmmindex1, custom_bond_link_par);
            }

            system_openmm->addForce(custom_link_bond);
        }

    } // end of bond link flag

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
