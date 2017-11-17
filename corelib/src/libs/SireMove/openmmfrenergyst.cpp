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

#define _GLIBCXX_USE_CXX11_ABI 0

#include "openmmfrenergyst.h"
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

#include "SireMove/flexibility.h"

#include "SireMaths/constants.h"

//ADDED BY GAC
#include "SireMaths/vector.h"
#include "SireMol/mgname.h"
#include "SireMol/perturbation.h"
#include "SireMM/internalperturbation.h"
#include <iostream>
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

enum
{
    NOCUTOFF = 0,
    CUTOFFNONPERIODIC = 1,
    CUTOFFPERIODIC = 2
};

enum
{
    NONE = 0,
    HBONDS = 1,
    ALLBONDS = 2,
    HANGLES = 3

};

static const RegisterMetaType<OpenMMFrEnergyST> r_openmmint;

/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds, const OpenMMFrEnergyST &velver)
{
    writeHeader(ds, r_openmmint, 1);

    SharedDataStream sds(ds);

    sds << velver.frequent_save_velocities << velver.molgroup << velver.solute
        << velver.solutehard << velver.solutetodummy << velver.solutefromdummy
        << velver.CutoffType << velver.cutoff_distance << velver.field_dielectric
        << velver.Andersen_flag << velver.Andersen_frequency
        << velver.MCBarostat_flag << velver.MCBarostat_frequency
        << velver.ConstraintType << velver.Pressure << velver.Temperature
        << velver.platform_type << velver.Restraint_flag
        << velver.CMMremoval_frequency << velver.buffer_frequency
        << velver.energy_frequency
        << velver.device_index << velver.precision << velver.Alchemical_value
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
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, OpenMMFrEnergyST &velver)
{
    VersionID v = readHeader(ds, r_openmmint);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> velver.frequent_save_velocities >> velver.molgroup
            >> velver.solute >> velver.solutehard >> velver.solutetodummy
            >> velver.solutefromdummy >> velver.CutoffType >> velver.cutoff_distance
            >> velver.field_dielectric >> velver.Andersen_flag
            >> velver.Andersen_frequency >> velver.MCBarostat_flag
            >> velver.MCBarostat_frequency >> velver.ConstraintType
            >> velver.Pressure >> velver.Temperature >> velver.platform_type
            >> velver.Restraint_flag >> velver.CMMremoval_frequency
            >> velver.buffer_frequency >> velver.energy_frequency
            >> velver.device_index >> velver.precision >> velver.Alchemical_value
            >> velver.coulomb_power >> velver.shift_delta >> velver.delta_alchemical
            >> velver.alchemical_array
            >> velver.Integrator_type >> velver.friction >> velver.integration_tol
            >> velver.timeskip >> velver.reinitialise_context
            >> static_cast<Integrator&> (velver);

        // Maybe....need to reinitialise from molgroup because openmm system was not serialised...
        velver.isSystemInitialised = false;
        velver.isContextInitialised = false;

        qDebug() << " Re-initialisation of OpenMMFrEnergyST from datastream";

        velver.initialise();
    }
    else
        throw version_error(v, "1", r_openmmint, CODELOC);

    return ds;
}

/** Constructor*/
OpenMMFrEnergyST::OpenMMFrEnergyST(bool frequent_save)
: ConcreteProperty<OpenMMFrEnergyST, Integrator>(),
frequent_save_velocities(frequent_save),
molgroup(MoleculeGroup()), solute(MoleculeGroup()), solutehard(MoleculeGroup()), solutetodummy(MoleculeGroup()), solutefromdummy(MoleculeGroup()),
openmm_system(0), openmm_context(0), isSystemInitialised(false), isContextInitialised(false),
CutoffType("nocutoff"), cutoff_distance(1.0 * nanometer), field_dielectric(78.3),
Andersen_flag(false), Andersen_frequency(90.0), MCBarostat_flag(false),
MCBarostat_frequency(25), ConstraintType("none"),
Pressure(1.0 * bar), Temperature(300.0 * kelvin), platform_type("Reference"), Restraint_flag(false),
CMMremoval_frequency(0), buffer_frequency(0), energy_frequency(100), device_index("0"), precision("single"), Alchemical_value(0.5), coulomb_power(0),
shift_delta(2.0), delta_alchemical(0.001), alchemical_array(),
    finite_diff_gradients(), pot_energies(), perturbed_energies(), reduced_perturbed_energies(),
    forward_Metropolis(), backward_Metropolis(),
Integrator_type("leapfrogverlet"), friction(1.0 / picosecond), integration_tol(0.001), timeskip(0.0 * picosecond),
reinitialise_context(false), Debug(false)
{
}

/** Constructor using the passed molecule groups */
OpenMMFrEnergyST::OpenMMFrEnergyST(const MoleculeGroup &molecule_group, const MoleculeGroup &solute_group, const MoleculeGroup &solute_hard, const MoleculeGroup &solute_todummy, const MoleculeGroup &solute_fromdummy, bool frequent_save)
: ConcreteProperty<OpenMMFrEnergyST, Integrator>(),
frequent_save_velocities(frequent_save),
molgroup(molecule_group), solute(solute_group), solutehard(solute_hard), solutetodummy(solute_todummy), solutefromdummy(solute_fromdummy),
openmm_system(0), openmm_context(0), isSystemInitialised(false), isContextInitialised(false),
CutoffType("nocutoff"), cutoff_distance(1.0 * nanometer), field_dielectric(78.3),
Andersen_flag(false), Andersen_frequency(90.0), MCBarostat_flag(false),
MCBarostat_frequency(25), ConstraintType("none"),
Pressure(1.0 * bar), Temperature(300.0 * kelvin), platform_type("Reference"), Restraint_flag(false),
CMMremoval_frequency(0), buffer_frequency(0), energy_frequency(100), device_index("0"), precision("single"), Alchemical_value(0.5), coulomb_power(0),
shift_delta(2.0), delta_alchemical(0.001), alchemical_array(), finite_diff_gradients(), pot_energies(), perturbed_energies(),
    reduced_perturbed_energies(), forward_Metropolis(), backward_Metropolis(),
Integrator_type("leapfrogverlet"), friction(1.0 / picosecond), integration_tol(0.001), timeskip(0.0 * picosecond),
reinitialise_context(false), Debug(false)
{
}

/** Copy constructor */
OpenMMFrEnergyST::OpenMMFrEnergyST(const OpenMMFrEnergyST &other)
: ConcreteProperty<OpenMMFrEnergyST, Integrator>(other),
frequent_save_velocities(other.frequent_save_velocities),
molgroup(other.molgroup), solute(other.solute), solutehard(other.solutehard),
solutetodummy(other.solutetodummy), solutefromdummy(other.solutefromdummy),
openmm_system(other.openmm_system), openmm_context(other.openmm_context), isSystemInitialised(other.isSystemInitialised),
isContextInitialised(other.isContextInitialised),
CutoffType(other.CutoffType), cutoff_distance(other.cutoff_distance),
field_dielectric(other.field_dielectric), Andersen_flag(other.Andersen_flag),
Andersen_frequency(other.Andersen_frequency), MCBarostat_flag(other.MCBarostat_flag),
MCBarostat_frequency(other.MCBarostat_frequency), ConstraintType(other.ConstraintType),
Pressure(other.Pressure), Temperature(other.Temperature), platform_type(other.platform_type),
Restraint_flag(other.Restraint_flag), CMMremoval_frequency(other.CMMremoval_frequency),
buffer_frequency(other.buffer_frequency), energy_frequency(other.energy_frequency), device_index(other.device_index),
precision(other.precision), Alchemical_value(other.Alchemical_value),
coulomb_power(other.coulomb_power), shift_delta(other.shift_delta),
delta_alchemical(other.delta_alchemical), alchemical_array(other.alchemical_array), finite_diff_gradients(other.finite_diff_gradients), pot_energies(other.pot_energies),
perturbed_energies(other.perturbed_energies),reduced_perturbed_energies(other.reduced_perturbed_energies),
    forward_Metropolis(other.forward_Metropolis), backward_Metropolis(other.backward_Metropolis),
Integrator_type(other.Integrator_type), friction(other.friction), integration_tol(other.integration_tol), timeskip(other.timeskip),
reinitialise_context(other.reinitialise_context), Debug(other.Debug)
{
}

/** Destructor */
OpenMMFrEnergyST::~OpenMMFrEnergyST()
{
    //delete openmm_system;
}

/** Copy assignment operator */
OpenMMFrEnergyST& OpenMMFrEnergyST::operator=(const OpenMMFrEnergyST &other)
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
    Alchemical_value = other.Alchemical_value;
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
bool OpenMMFrEnergyST::operator==(const OpenMMFrEnergyST &other) const
{
    return frequent_save_velocities == other.frequent_save_velocities
        and isSystemInitialised == other.isSystemInitialised
        and isContextInitialised == other.isContextInitialised
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
        and Alchemical_value == other.Alchemical_value
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
bool OpenMMFrEnergyST::operator!=(const OpenMMFrEnergyST &other) const
{
    return not OpenMMFrEnergyST::operator==(other);
}

/** Return a string representation of this integrator */
QString OpenMMFrEnergyST::toString() const
{
    return QObject::tr("OpenMMFrEnergyST()");
}

/**
 * initialises the openMM Free energy single topology calculation
 * Initialise must be called before anything else happens.
 */
void OpenMMFrEnergyST::initialise()
{

    bool Debug = false;
    if (Debug)
    {
        qDebug() << "Initialising OpenMMFrEnergyST";
        const std::string version = OpenMM::Platform::getOpenMMVersion();
        qDebug() << "OpenMM Version: " << QString::fromUtf8(version.data(), version.size());
    }

    // Create a workspace using the stored molgroup

    const MoleculeGroup moleculegroup = this->molgroup.read();

    if (moleculegroup.isEmpty())
    {
        throw SireError::program_bug(QObject::tr("Cannot initialise OpenMMFrEnergyST because molgroup has not been defined"), CODELOC);
    }

    const MoleculeGroup solute = this->solute.read();

    //if ( solute.isEmpty() ){
    //    throw SireError::program_bug(QObject::tr("Cannot initialise OpenMMFrEnergyST because solute group has not been defined"), CODELOC);
    //}


    const MoleculeGroup solutehard = this->solutehard.read();

    const MoleculeGroup solutetodummy = this->solutetodummy.read();

    const MoleculeGroup solutefromdummy = this->solutefromdummy.read();

    AtomicVelocityWorkspace ws = this->createWorkspace(moleculegroup).read().asA<AtomicVelocityWorkspace>();

    const int nmols = ws.nMolecules();

    int nats = 0;

    for (int i = 0; i < nmols; ++i)
    {
        nats = nats + ws.nAtoms(i);
    }

    if (Debug)
        qDebug() << "There are " << nats << " atoms " << "There are " << nmols << " molecules" << "\n";

    int flag_cutoff;
    int flag_constraint;

    if (CutoffType == "nocutoff")
        flag_cutoff = NOCUTOFF;
    else if (CutoffType == "cutoffnonperiodic")
        flag_cutoff = CUTOFFNONPERIODIC;
    else if (CutoffType == "cutoffperiodic")
        flag_cutoff = CUTOFFPERIODIC;
    else
        throw SireError::program_bug(QObject::tr("The CutOff method has not been specified. Possible choises: nocutoff, cutoffnonperiodic, cutoffperiodic"), CODELOC);

    if (Debug)
        qDebug() << "\nCutoffType = " << CutoffType << "\n";

    bool flag_noperturbedconstraints = false;
    bool flag_constraint_water = false;
    if (ConstraintType == "none")
        flag_constraint = NONE;
    else if (ConstraintType == "hbonds")
        flag_constraint = HBONDS;
    else if (ConstraintType == "allbonds")
        flag_constraint = ALLBONDS;
    else if (ConstraintType == "hangles")
        flag_constraint = HANGLES;
    else if (ConstraintType == "hbonds-notperturbed")
    {
        flag_constraint = HBONDS;
        flag_noperturbedconstraints = true;
    }
    else if (ConstraintType == "none-notwater")
    {
        flag_constraint = NONE;
        flag_constraint_water = true;
    }
    else
        throw SireError::program_bug(QObject::tr("The Constraints method has not been specified."
        "Possible choises: none, hbonds, allbonds, hangles, hbonds-notperturbed, none-notwater"), CODELOC);

    if (Debug)
        qDebug() << "\nConstraint Type = " << ConstraintType << "\n";

    //Load Plugins from the OpenMM standard Plugin Directory
    OpenMM::Platform::loadPluginsFromDirectory(OpenMM::Platform::getDefaultPluginsDirectory());

    OpenMM::System * system_openmm = new OpenMM::System();

    system_openmm->setDefaultPeriodicBoxVectors(OpenMM::Vec3(6, 0, 0),
                                                OpenMM::Vec3(0, 6, 0),
                                                OpenMM::Vec3(0, 0, 6));

    //The Standard Non Bonded is only defined to extract 1-2,1-3,1-4 pairs from the system
    OpenMM::NonbondedForce * nonbond_openmm = new OpenMM::NonbondedForce();

    nonbond_openmm->setUseDispersionCorrection(false);

    //CUSTOM NON BONDED FORCE FIELD

    OpenMM::CustomNonbondedForce * custom_force_field = NULL;

    //1-4 interactions
    OpenMM::CustomBondForce * custom_intra_14_clj = NULL;
    OpenMM::CustomBondForce * custom_intra_14_todummy = NULL;
    OpenMM::CustomBondForce * custom_intra_14_fromdummy = NULL;
    OpenMM::CustomBondForce * custom_intra_14_fromdummy_todummy = NULL;


    if (flag_cutoff == NOCUTOFF)
    {

        if (coulomb_power > 0)
        { //This is necessary to avoid nan errors on the GPUs platform caused by the calculation of 0^0

            custom_force_field = new OpenMM::CustomNonbondedForce("(1.0 - isSolvent1 * isSolvent2 * SPOnOff) * (Hcs + Hls);"
                                                                  "Hcs = (lambda^n) * 138.935456 * q_prod/sqrt(diff_cl+r^2);"
                                                                  "diff_cl = (1.0-lambda) * 0.01;"
                                                                  "Hls = 4.0 * eps_avg * (LJ*LJ-LJ);"
                                                                  "LJ=((sigma_avg * sigma_avg)/soft)^3;"
                                                                  "soft=(diff_lj*delta*sigma_avg + r*r);"
                                                                  "diff_lj=(1.0-lambda) * 0.1;"
                                                                  "lambda = Logic_lam * lam + Logic_om_lam * (1.0-lam) + Logic_mix_lam * max(lam,1.0-lam) + Logic_hard;"
                                                                  "Logic_hard = isHD1 * isHD2 * (1.0-isTD1) * (1.0-isTD2) * (1.0-isFD1) * (1.0-isFD2);"
                                                                  "Logic_om_lam = max((1.0-isHD1)*(1.0-isHD2)*isTD1*isTD2*(1.0-isFD1)*(1.0-isFD2), B_om_lam);"
                                                                  "B_om_lam = max(isHD1*(1.0-isHD2)*isTD1*(1.0-isTD2)*(1.0-isFD1)*(1.0-isFD2), C_om_lam);"
                                                                  "C_om_lam = max(isHD1*(1.0-isHD2)*(1.0-isTD1)*isTD2*(1.0-isFD1)*(1.0-isFD2) , D_om_lam);"
                                                                  "D_om_lam = max((1.0-isHD1)*isHD2*isTD1*(1.0-isTD2)*(1.0-isFD1)*(1.0-isFD2), E_om_lam);"
                                                                  "E_om_lam = (1.0-isHD1)*isHD2*(1.0-isTD1)*isTD2*(1.0-isFD1)*(1.0-isFD2);"
                                                                  "Logic_lam = max((1.0-isHD1)*(1.0-isHD2)*(1.0-isTD1)*(1.0-isTD2)*isFD1*isFD2, B_lam);"
                                                                  "B_lam = max(isHD1*(1.0-isHD2)*(1.0-isTD1)*(1.0-isTD2)*isFD1*(1.0-isFD2), C_lam);"
                                                                  "C_lam = max(isHD1*(1.0-isHD2)*(1.0-isTD1)*(1.0-isTD2)*(1.0-isFD1)*isFD2 , D_lam);"
                                                                  "D_lam = max((1.0-isHD1)*isHD2*(1.0-isTD1)*(1.0-isTD2)*isFD1*(1.0-isFD2), E_lam);"
                                                                  "E_lam = (1.0-isHD1)*isHD2*(1.0-isTD1)*(1.0-isTD2)*(1.0-isFD1)*isFD2;"
                                                                  "Logic_mix_lam = max((1.0-isHD1)*(1.0-isHD2)*isTD1*(1.0-isTD2)*isFD1*(1.0-isFD2), B_mix);"
                                                                  "B_mix = max((1.0-isHD1)*(1.0-isHD2)*isTD1*(1.0-isTD2)*(1.0-isFD1)*isFD2, C_mix);"
                                                                  "C_mix = max((1.0-isHD1)*(1.0-isHD2)*(1.0-isTD1)*isTD2*isFD1*(1.0-isFD2) , D_mix);"
                                                                  "D_mix= (1.0-isHD1)*(1.0-isHD2)*(1.0-isTD1)*isTD2*(1.0-isFD1)*isFD2;"
                                                                  "q_prod = (qend1 * lam+(1.0-lam) * qstart1) * (qend2 * lam+(1.0-lam) * qstart2);"
                                                                  "eps_avg = sqrt((epend1*lam+(1.0-lam)*epstart1)*(epend2*lam+(1.0-lam)*epstart2));"
                                                                  "sigma_avg = 0.5*((sigmaend1*lam+(1.0-lam)*sigmastart1)+(sigmaend2*lam+(1.0-lam)*sigmastart2))");

            custom_force_field->addGlobalParameter("lam", Alchemical_value);
            custom_force_field->addGlobalParameter("delta", shift_delta);
            custom_force_field->addGlobalParameter("n", coulomb_power);
            custom_force_field->addGlobalParameter("SPOnOff", 0.0);

            custom_force_field->setNonbondedMethod(OpenMM::CustomNonbondedForce::NoCutoff);



            custom_intra_14_todummy = new OpenMM::CustomBondForce("Hcs + Hls;"
                                                                  "Hcs=(lamtd^ntd)*138.935456*q_prod/sqrt(diff_cl+r^2);"
                                                                  "diff_cl=(1.0-lamtd)*0.01;"
                                                                  "Hls=4.0*eps_avg*(LJ*LJ-LJ);"
                                                                  "LJ=((sigma_avg*sigma_avg)/soft)^3;"
                                                                  "soft=(diff_lj*deltatd*sigma_avg+r*r);"
                                                                  "diff_lj=(1.0-lamtd)*0.1;"
                                                                  "eps_avg = sqrt((1-lamtd)*(1-lamtd)*eaend + lamtd*lamtd*eastart + lamtd*(1-lamtd)*emix);"
                                                                  "sigma_avg = (1-lamtd)*saend + lamtd*sastart;"
                                                                  "q_prod = (1-lamtd)*(1-lamtd)*qpend + lamtd*lamtd*qpstart + lamtd*(1-lamtd)*qmix");


            custom_intra_14_todummy->addGlobalParameter("lamtd", 1.0 - Alchemical_value);
            custom_intra_14_todummy->addGlobalParameter("deltatd", shift_delta);
            custom_intra_14_todummy->addGlobalParameter("ntd", coulomb_power);


            custom_intra_14_fromdummy = new OpenMM::CustomBondForce("Hcs + Hls;"
                                                                    "Hcs=(lamfd^nfd)*138.935456*q_prod/sqrt(diff_cl+r^2);"
                                                                    "diff_cl=(1.0-lamfd)*0.01;"
                                                                    "Hls=4.0*eps_avg*(LJ*LJ-LJ);"
                                                                    "LJ=((sigma_avg*sigma_avg)/soft)^3;"
                                                                    "soft=(diff_lj*deltafd*sigma_avg+r*r);"
                                                                    "diff_lj=(1.0-lamfd)*0.1;"
                                                                    "eps_avg = sqrt(lamfd*lamfd*eaend + (1-lamfd)*(1-lamfd)*eastart + lamfd*(1-lamfd)*emix);"
                                                                    "sigma_avg = lamfd*saend + (1-lamfd)*sastart;"
                                                                    "q_prod = lamfd*lamfd*qpend + (1-lamfd)*(1-lamfd)*qpstart + lamfd*(1-lamfd)*qmix");

            custom_intra_14_fromdummy->addGlobalParameter("lamfd", Alchemical_value);
            custom_intra_14_fromdummy->addGlobalParameter("deltafd", shift_delta);
            custom_intra_14_fromdummy->addGlobalParameter("nfd", coulomb_power);


            custom_intra_14_fromdummy_todummy = new OpenMM::CustomBondForce("Hcs + Hls;"
                                                                            "Hcs=(lamFTD^nftd)*138.935456*q_prod/sqrt(diff_cl+r^2);"
                                                                            "diff_cl=(1.0-lamFTD)*0.01;"
                                                                            "Hls=4.0*eps_avg*(LJ*LJ-LJ);"
                                                                            "LJ=((sigma_avg*sigma_avg)/soft)^3;"
                                                                            "soft=(diff_lj*deltaftd*sigma_avg+r*r);"
                                                                            "diff_lj=(1.0-lamFTD)*0.1;"
                                                                            "eps_avg = sqrt(lamftd*lamftd*eaend + (1-lamftd)*(1-lamftd)*eastart + lamftd*(1-lamftd)*emix);"
                                                                            "sigma_avg = lamftd*saend + (1-lamftd)*sastart;"
                                                                            "q_prod = lamftd*lamftd*qpend + (1-lamftd)*(1-lamftd)*qpstart + lamftd*(1-lamftd)*qmix;"
                                                                            "lamFTD = max(lamftd,1-lamftd)");

            custom_intra_14_fromdummy_todummy->addGlobalParameter("lamftd", Alchemical_value);
            custom_intra_14_fromdummy_todummy->addGlobalParameter("deltaftd", shift_delta);
            custom_intra_14_fromdummy_todummy->addGlobalParameter("nftd", coulomb_power);

        }
        else
        {// coulomb_power == 0. //This is necessary to avoid nan errors on the GPUs platform caused by the calculation of 0^0

            custom_force_field = new OpenMM::CustomNonbondedForce("(1.0 - isSolvent1 * isSolvent2 * SPOnOff) * (Hcs + Hls);"
                                                                  "Hcs = 138.935456 * q_prod/sqrt(diff_cl+r^2);"
                                                                  "diff_cl = (1.0-lambda) * 0.01;"
                                                                  "Hls = 4.0 * eps_avg * (LJ*LJ-LJ);"
                                                                  "LJ=((sigma_avg * sigma_avg)/soft)^3;"
                                                                  "soft=(diff_lj*delta*sigma_avg + r*r);"
                                                                  "diff_lj=(1.0-lambda) * 0.1;"
                                                                  "lambda = Logic_lam * lam + Logic_om_lam * (1.0-lam) + Logic_mix_lam * max(lam,1.0-lam) + Logic_hard;"
                                                                  "Logic_hard = isHD1 * isHD2 * (1.0-isTD1) * (1.0-isTD2) * (1.0-isFD1) * (1.0-isFD2);"
                                                                  "Logic_om_lam = max((1.0-isHD1)*(1.0-isHD2)*isTD1*isTD2*(1.0-isFD1)*(1.0-isFD2), B_om_lam);"
                                                                  "B_om_lam = max(isHD1*(1.0-isHD2)*isTD1*(1.0-isTD2)*(1.0-isFD1)*(1.0-isFD2), C_om_lam);"
                                                                  "C_om_lam = max(isHD1*(1.0-isHD2)*(1.0-isTD1)*isTD2*(1.0-isFD1)*(1.0-isFD2) , D_om_lam);"
                                                                  "D_om_lam = max((1.0-isHD1)*isHD2*isTD1*(1.0-isTD2)*(1.0-isFD1)*(1.0-isFD2), E_om_lam);"
                                                                  "E_om_lam = (1.0-isHD1)*isHD2*(1.0-isTD1)*isTD2*(1.0-isFD1)*(1.0-isFD2);"
                                                                  "Logic_lam = max((1.0-isHD1)*(1.0-isHD2)*(1.0-isTD1)*(1.0-isTD2)*isFD1*isFD2, B_lam);"
                                                                  "B_lam = max(isHD1*(1.0-isHD2)*(1.0-isTD1)*(1.0-isTD2)*isFD1*(1.0-isFD2), C_lam);"
                                                                  "C_lam = max(isHD1*(1.0-isHD2)*(1.0-isTD1)*(1.0-isTD2)*(1.0-isFD1)*isFD2 , D_lam);"
                                                                  "D_lam = max((1.0-isHD1)*isHD2*(1.0-isTD1)*(1.0-isTD2)*isFD1*(1.0-isFD2), E_lam);"
                                                                  "E_lam = (1.0-isHD1)*isHD2*(1.0-isTD1)*(1.0-isTD2)*(1.0-isFD1)*isFD2;"
                                                                  "Logic_mix_lam = max((1.0-isHD1)*(1.0-isHD2)*isTD1*(1.0-isTD2)*isFD1*(1.0-isFD2), B_mix);"
                                                                  "B_mix = max((1.0-isHD1)*(1.0-isHD2)*isTD1*(1.0-isTD2)*(1.0-isFD1)*isFD2, C_mix);"
                                                                  "C_mix = max((1.0-isHD1)*(1.0-isHD2)*(1.0-isTD1)*isTD2*isFD1*(1.0-isFD2) , D_mix);"
                                                                  "D_mix= (1.0-isHD1)*(1.0-isHD2)*(1.0-isTD1)*isTD2*(1.0-isFD1)*isFD2;"
                                                                  "q_prod = (qend1 * lam+(1.0-lam) * qstart1) * (qend2 * lam+(1.0-lam) * qstart2);"
                                                                  "eps_avg = sqrt((epend1*lam+(1.0-lam)*epstart1)*(epend2*lam+(1.0-lam)*epstart2));"
                                                                  "sigma_avg = 0.5*((sigmaend1*lam+(1.0-lam)*sigmastart1)+(sigmaend2*lam+(1.0-lam)*sigmastart2))");

            custom_force_field->addGlobalParameter("lam", Alchemical_value);
            custom_force_field->addGlobalParameter("delta", shift_delta);
            custom_force_field->addGlobalParameter("n", coulomb_power);
            custom_force_field->addGlobalParameter("SPOnOff", 0.0);

            custom_force_field->setNonbondedMethod(OpenMM::CustomNonbondedForce::NoCutoff);



            custom_intra_14_todummy = new OpenMM::CustomBondForce("Hcs + Hls;"
                                                                  "Hcs=138.935456*q_prod/sqrt(diff_cl+r^2);"
                                                                  "diff_cl=(1.0-lamtd)*0.01;"
                                                                  "Hls=4.0*eps_avg*(LJ*LJ-LJ);"
                                                                  "LJ=((sigma_avg*sigma_avg)/soft)^3;"
                                                                  "soft=(diff_lj*deltatd*sigma_avg+r*r);"
                                                                  "diff_lj=(1.0-lamtd)*0.1;"
                                                                  "eps_avg = sqrt((1-lamtd)*(1-lamtd)*eaend + lamtd*lamtd*eastart + lamtd*(1-lamtd)*emix);"
                                                                  "sigma_avg = (1-lamtd)*saend + lamtd*sastart;"
                                                                  "q_prod = (1-lamtd)*(1-lamtd)*qpend + lamtd*lamtd*qpstart + lamtd*(1-lamtd)*qmix");


            custom_intra_14_todummy->addGlobalParameter("lamtd", 1.0 - Alchemical_value);
            custom_intra_14_todummy->addGlobalParameter("deltatd", shift_delta);
            custom_intra_14_todummy->addGlobalParameter("ntd", coulomb_power);


            custom_intra_14_fromdummy = new OpenMM::CustomBondForce("Hcs + Hls;"
                                                                    "Hcs=138.935456*q_prod/sqrt(diff_cl+r^2);"
                                                                    "diff_cl=(1.0-lamfd)*0.01;"
                                                                    "Hls=4.0*eps_avg*(LJ*LJ-LJ);"
                                                                    "LJ=((sigma_avg*sigma_avg)/soft)^3;"
                                                                    "soft=(diff_lj*deltafd*sigma_avg+r*r);"
                                                                    "diff_lj=(1.0-lamfd)*0.1;"
                                                                    "eps_avg = sqrt(lamfd*lamfd*eaend + (1-lamfd)*(1-lamfd)*eastart + lamfd*(1-lamfd)*emix);"
                                                                    "sigma_avg = lamfd*saend + (1-lamfd)*sastart;"
                                                                    "q_prod = lamfd*lamfd*qpend + (1-lamfd)*(1-lamfd)*qpstart + lamfd*(1-lamfd)*qmix");

            custom_intra_14_fromdummy->addGlobalParameter("lamfd", Alchemical_value);
            custom_intra_14_fromdummy->addGlobalParameter("deltafd", shift_delta);
            custom_intra_14_fromdummy->addGlobalParameter("nfd", coulomb_power);


            custom_intra_14_fromdummy_todummy = new OpenMM::CustomBondForce("Hcs + Hls;"
                                                                            "Hcs=138.935456*q_prod/sqrt(diff_cl+r^2);"
                                                                            "diff_cl=(1.0-lamFTD)*0.01;"
                                                                            "Hls=4.0*eps_avg*(LJ*LJ-LJ);"
                                                                            "LJ=((sigma_avg*sigma_avg)/soft)^3;"
                                                                            "soft=(diff_lj*deltaftd*sigma_avg+r*r);"
                                                                            "diff_lj=(1.0-lamFTD)*0.1;"
                                                                            "eps_avg = sqrt(lamftd*lamftd*eaend + (1-lamftd)*(1-lamftd)*eastart + lamftd*(1-lamftd)*emix);"
                                                                            "sigma_avg = lamftd*saend + (1-lamftd)*sastart;"
                                                                            "q_prod = lamftd*lamftd*qpend + (1-lamftd)*(1-lamftd)*qpstart + lamftd*(1-lamftd)*qmix;"
                                                                            "lamFTD = max(lamftd,1-lamftd)");

            custom_intra_14_fromdummy_todummy->addGlobalParameter("lamftd", Alchemical_value);
            custom_intra_14_fromdummy_todummy->addGlobalParameter("deltaftd", shift_delta);
            custom_intra_14_fromdummy_todummy->addGlobalParameter("nftd", coulomb_power);


        }



        custom_intra_14_clj = new OpenMM::CustomBondForce("Hl+Hc;"
                                                          "Hl=4*eps_avg*((sigma_avg/r)^12-(sigma_avg/r)^6);"
                                                          "Hc=138.935456*q_prod/r;"
                                                          "eps_avg = sqrt(lamhd*lamhd*eaend + (1-lamhd)*(1-lamhd)*eastart + lamhd*(1-lamhd)*emix);"
                                                          "sigma_avg = lamhd*saend + (1-lamhd)*sastart;"
                                                          "q_prod = lamhd*lamhd*qpend + (1-lamhd)*(1-lamhd)*qpstart + lamhd*(1-lamhd)*qmix");

        custom_intra_14_clj->addGlobalParameter("lamhd", Alchemical_value);

        if (Debug)
        {
            qDebug() << "\nCut off type = " << CutoffType << "\n";
            qDebug() << "Lambda = " << Alchemical_value << " Coulomb Power = " << coulomb_power << " Delta Shift = " << shift_delta << "\n";
        }

    }
    else
    {//CUTOFF PERIODIC OR NON PERIODIC

        const double converted_cutoff_distance = convertTo(cutoff_distance.value(), nanometer);

        double eps2 = (field_dielectric - 1.0) / (2 * field_dielectric + 1.0);
        double kvalue = eps2 / (converted_cutoff_distance * converted_cutoff_distance * converted_cutoff_distance);
        double cvalue = (1.0 / converted_cutoff_distance)*(3.0 * field_dielectric) / (2.0 * field_dielectric + 1.0);

        if (coulomb_power > 0)
        {//This is necessary to avoid nan errors on the GPUs platform caused by the calculation of 0^0

            custom_force_field = new OpenMM::CustomNonbondedForce("(1.0 - isSolvent1 * isSolvent2 * SPOnOff) * (Hls + Hcs);"
                                                                  "Hcs = (lambda^n) * 138.935456 * q_prod*(1/sqrt(diff_cl+r*r) + krflam*(diff_cl+r*r)-crflam);"
                                                                  "crflam = crf * src;"
                                                                  "krflam = krf * src * src * src;"
                                                                  "src = cutoff/sqrt(diff_cl + cutoff*cutoff);"
                                                                  "diff_cl = (1.0-lambda) * 0.01;"
                                                                  "Hls = 4.0 * eps_avg * (LJ*LJ-LJ);"
                                                                  "LJ=((sigma_avg * sigma_avg)/soft)^3;"
                                                                  "soft=(diff_lj*delta*sigma_avg + r*r);"
                                                                  "diff_lj=(1.0-lambda) * 0.1;"
                                                                  "lambda = Logic_lam * lam + Logic_om_lam * (1.0-lam) + Logic_mix_lam * max(lam,1.0-lam) + Logic_hard;"
                                                                  "Logic_hard = isHD1 * isHD2 * (1.0-isTD1) * (1.0-isTD2) * (1.0-isFD1) * (1.0-isFD2);"
                                                                  "Logic_om_lam = max((1.0-isHD1)*(1.0-isHD2)*isTD1*isTD2*(1.0-isFD1)*(1.0-isFD2), B_om_lam);"
                                                                  "B_om_lam = max(isHD1*(1.0-isHD2)*isTD1*(1.0-isTD2)*(1.0-isFD1)*(1.0-isFD2), C_om_lam);"
                                                                  "C_om_lam = max(isHD1*(1.0-isHD2)*(1.0-isTD1)*isTD2*(1.0-isFD1)*(1.0-isFD2) , D_om_lam);"
                                                                  "D_om_lam = max((1.0-isHD1)*isHD2*isTD1*(1.0-isTD2)*(1.0-isFD1)*(1.0-isFD2), E_om_lam);"
                                                                  "E_om_lam = (1.0-isHD1)*isHD2*(1.0-isTD1)*isTD2*(1.0-isFD1)*(1.0-isFD2);"
                                                                  "Logic_lam = max((1.0-isHD1)*(1.0-isHD2)*(1.0-isTD1)*(1.0-isTD2)*isFD1*isFD2, B_lam);"
                                                                  "B_lam = max(isHD1*(1.0-isHD2)*(1.0-isTD1)*(1.0-isTD2)*isFD1*(1.0-isFD2), C_lam);"
                                                                  "C_lam = max(isHD1*(1.0-isHD2)*(1.0-isTD1)*(1.0-isTD2)*(1.0-isFD1)*isFD2 , D_lam);"
                                                                  "D_lam = max((1.0-isHD1)*isHD2*(1.0-isTD1)*(1.0-isTD2)*isFD1*(1.0-isFD2), E_lam);"
                                                                  "E_lam = (1.0-isHD1)*isHD2*(1.0-isTD1)*(1.0-isTD2)*(1.0-isFD1)*isFD2;"
                                                                  "Logic_mix_lam = max((1.0-isHD1)*(1.0-isHD2)*isTD1*(1.0-isTD2)*isFD1*(1.0-isFD2), B_mix);"
                                                                  "B_mix = max((1.0-isHD1)*(1.0-isHD2)*isTD1*(1.0-isTD2)*(1.0-isFD1)*isFD2, C_mix);"
                                                                  "C_mix = max((1.0-isHD1)*(1.0-isHD2)*(1.0-isTD1)*isTD2*isFD1*(1.0-isFD2) , D_mix);"
                                                                  "D_mix= (1.0-isHD1)*(1.0-isHD2)*(1.0-isTD1)*isTD2*(1.0-isFD1)*isFD2;"
                                                                  "q_prod = (qend1 * lam+(1.0-lam) * qstart1) * (qend2 * lam+(1.0-lam) * qstart2);"
                                                                  "eps_avg = sqrt((epend1*lam+(1.0-lam)*epstart1)*(epend2*lam+(1.0-lam)*epstart2));"
                                                                  "sigma_avg = 0.5*((sigmaend1*lam+(1.0-lam)*sigmastart1)+(sigmaend2*lam+(1.0-lam)*sigmastart2))");

            custom_force_field->setCutoffDistance(converted_cutoff_distance);

            custom_force_field->addGlobalParameter("lam", Alchemical_value);
            custom_force_field->addGlobalParameter("delta", shift_delta);
            custom_force_field->addGlobalParameter("n", coulomb_power);
            custom_force_field->addGlobalParameter("krf", kvalue);
            custom_force_field->addGlobalParameter("crf", cvalue);
            custom_force_field->addGlobalParameter("cutoff", converted_cutoff_distance);
            custom_force_field->addGlobalParameter("SPOnOff", 0.0);


            if (flag_cutoff == CUTOFFNONPERIODIC)
            {
                custom_force_field->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);

            }
            else
            {
                custom_force_field->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffPeriodic);
            }

            //NO REACTION FIELD IS APPLIED TO 1-4 INTERACTIONS. If the scaling factor is one (Glycam ff) then the OpenMM potential energy is not equal to he Sire energy. This is caused by the application of the reaction field on the 14 pairs in Sire.


            custom_intra_14_todummy = new OpenMM::CustomBondForce("withinCutoff*(Hcs + Hls);"
                                                                  "withinCutoff=step(cutofftd-r);"
                                                                  "Hcs=(lamtd^ntd)*138.935456*q_prod/sqrt(diff_cl+r^2);"
                                                                  "diff_cl=(1.0-lamtd)*0.01;"
                                                                  "Hls=4.0*eps_avg*(LJ*LJ-LJ);"
                                                                  "LJ=((sigma_avg*sigma_avg)/soft)^3;"
                                                                  "soft=(diff_lj*deltatd*sigma_avg+r*r);"
                                                                  "diff_lj=(1.0-lamtd)*0.1;"
                                                                  "eps_avg = sqrt((1-lamtd)*(1-lamtd)*eaend + lamtd*lamtd*eastart + lamtd*(1-lamtd)*emix);"
                                                                  "sigma_avg = (1-lamtd)*saend + lamtd*sastart;"
                                                                  "q_prod = (1-lamtd)*(1-lamtd)*qpend + lamtd*lamtd*qpstart + lamtd*(1-lamtd)*qmix");

            custom_intra_14_todummy->addGlobalParameter("lamtd", 1.0 - Alchemical_value);
            custom_intra_14_todummy->addGlobalParameter("deltatd", shift_delta);
            custom_intra_14_todummy->addGlobalParameter("ntd", coulomb_power);
            custom_intra_14_todummy->addGlobalParameter("cutofftd", converted_cutoff_distance);

            custom_intra_14_fromdummy = new OpenMM::CustomBondForce("withinCutoff*(Hcs + Hls);"
                                                                    "withinCutoff=step(cutofffd-r);"
                                                                    "Hcs=(lamfd^nfd)*138.935456*q_prod/sqrt(diff_cl+r^2);"
                                                                    "diff_cl=(1.0-lamfd)*0.01;"
                                                                    "Hls=4.0*eps_avg*(LJ*LJ-LJ);"
                                                                    "LJ=((sigma_avg*sigma_avg)/soft)^3;"
                                                                    "soft=(diff_lj*deltafd*sigma_avg+r*r);"
                                                                    "diff_lj=(1.0-lamfd)*0.1;"
                                                                    "eps_avg = sqrt(lamfd*lamfd*eaend + (1-lamfd)*(1-lamfd)*eastart + lamfd*(1-lamfd)*emix);"
                                                                    "sigma_avg = lamfd*saend + (1-lamfd)*sastart;"
                                                                    "q_prod = lamfd*lamfd*qpend + (1-lamfd)*(1-lamfd)*qpstart + lamfd*(1-lamfd)*qmix");

            custom_intra_14_fromdummy->addGlobalParameter("lamfd", Alchemical_value);
            custom_intra_14_fromdummy->addGlobalParameter("deltafd", shift_delta);
            custom_intra_14_fromdummy->addGlobalParameter("nfd", coulomb_power);
            custom_intra_14_fromdummy->addGlobalParameter("cutofffd", converted_cutoff_distance);


            custom_intra_14_fromdummy_todummy = new OpenMM::CustomBondForce("withinCutoff*(Hcs + Hls);"
                                                                            "withinCutoff=step(cutoffftd-r);"
                                                                            "Hcs=(lamFTD^nftd)*138.935456*q_prod/sqrt(diff_cl+r^2);"
                                                                            "diff_cl=(1.0-lamFTD)*0.01;"
                                                                            "Hls=4.0*eps_avg*(LJ*LJ-LJ);"
                                                                            "LJ=((sigma_avg*sigma_avg)/soft)^3;"
                                                                            "soft=(diff_lj*deltaftd*sigma_avg+r*r);"
                                                                            "diff_lj=(1.0-lamFTD)*0.1;"
                                                                            "eps_avg = sqrt(lamftd*lamftd*eaend + (1-lamftd)*(1-lamftd)*eastart + lamftd*(1-lamftd)*emix);"
                                                                            "sigma_avg = lamftd*saend + (1-lamftd)*sastart;"
                                                                            "q_prod = lamftd*lamftd*qpend + (1-lamftd)*(1-lamftd)*qpstart + lamftd*(1-lamftd)*qmix;"
                                                                            "lamFTD = max(lamftd,1-lamftd)");

            custom_intra_14_fromdummy_todummy->addGlobalParameter("lamftd", Alchemical_value);
            custom_intra_14_fromdummy_todummy->addGlobalParameter("deltaftd", shift_delta);
            custom_intra_14_fromdummy_todummy->addGlobalParameter("nftd", coulomb_power);
            custom_intra_14_fromdummy_todummy->addGlobalParameter("cutoffftd", converted_cutoff_distance);


        }
        else
        {//coulomb_power == 0. //This is necessary to avoid nan errors on the GPUs platform caused by the calculation of 0^0

            custom_force_field = new OpenMM::CustomNonbondedForce("(1.0 - isSolvent1 * isSolvent2 * SPOnOff) * (Hls + Hcs);"
                                                                  "Hcs = 138.935456 * q_prod*(1/sqrt(diff_cl+r*r) + krflam*(diff_cl+r*r)-crflam);"
                                                                  "crflam = crf * src;"
                                                                  "krflam = krf * src * src * src;"
                                                                  "src = cutoff/sqrt(diff_cl + cutoff*cutoff);"
                                                                  "diff_cl = (1.0-lambda) * 0.01;"
                                                                  "Hls = 4.0 * eps_avg * (LJ*LJ-LJ);"
                                                                  "LJ=((sigma_avg * sigma_avg)/soft)^3;"
                                                                  "soft=(diff_lj*delta*sigma_avg + r*r);"
                                                                  "diff_lj=(1.0-lambda) * 0.1;"
                                                                  "lambda = Logic_lam * lam + Logic_om_lam * (1.0-lam) + Logic_mix_lam * max(lam,1.0-lam) + Logic_hard;"
                                                                  "Logic_hard = isHD1 * isHD2 * (1.0-isTD1) * (1.0-isTD2) * (1.0-isFD1) * (1.0-isFD2);"
                                                                  "Logic_om_lam = max((1.0-isHD1)*(1.0-isHD2)*isTD1*isTD2*(1.0-isFD1)*(1.0-isFD2), B_om_lam);"
                                                                  "B_om_lam = max(isHD1*(1.0-isHD2)*isTD1*(1.0-isTD2)*(1.0-isFD1)*(1.0-isFD2), C_om_lam);"
                                                                  "C_om_lam = max(isHD1*(1.0-isHD2)*(1.0-isTD1)*isTD2*(1.0-isFD1)*(1.0-isFD2) , D_om_lam);"
                                                                  "D_om_lam = max((1.0-isHD1)*isHD2*isTD1*(1.0-isTD2)*(1.0-isFD1)*(1.0-isFD2), E_om_lam);"
                                                                  "E_om_lam = (1.0-isHD1)*isHD2*(1.0-isTD1)*isTD2*(1.0-isFD1)*(1.0-isFD2);"
                                                                  "Logic_lam = max((1.0-isHD1)*(1.0-isHD2)*(1.0-isTD1)*(1.0-isTD2)*isFD1*isFD2, B_lam);"
                                                                  "B_lam = max(isHD1*(1.0-isHD2)*(1.0-isTD1)*(1.0-isTD2)*isFD1*(1.0-isFD2), C_lam);"
                                                                  "C_lam = max(isHD1*(1.0-isHD2)*(1.0-isTD1)*(1.0-isTD2)*(1.0-isFD1)*isFD2 , D_lam);"
                                                                  "D_lam = max((1.0-isHD1)*isHD2*(1.0-isTD1)*(1.0-isTD2)*isFD1*(1.0-isFD2), E_lam);"
                                                                  "E_lam = (1.0-isHD1)*isHD2*(1.0-isTD1)*(1.0-isTD2)*(1.0-isFD1)*isFD2;"
                                                                  "Logic_mix_lam = max((1.0-isHD1)*(1.0-isHD2)*isTD1*(1.0-isTD2)*isFD1*(1.0-isFD2), B_mix);"
                                                                  "B_mix = max((1.0-isHD1)*(1.0-isHD2)*isTD1*(1.0-isTD2)*(1.0-isFD1)*isFD2, C_mix);"
                                                                  "C_mix = max((1.0-isHD1)*(1.0-isHD2)*(1.0-isTD1)*isTD2*isFD1*(1.0-isFD2) , D_mix);"
                                                                  "D_mix= (1.0-isHD1)*(1.0-isHD2)*(1.0-isTD1)*isTD2*(1.0-isFD1)*isFD2;"
                                                                  "q_prod = (qend1 * lam+(1.0-lam) * qstart1) * (qend2 * lam+(1.0-lam) * qstart2);"
                                                                  "eps_avg = sqrt((epend1*lam+(1.0-lam)*epstart1)*(epend2*lam+(1.0-lam)*epstart2));"
                                                                  "sigma_avg = 0.5*((sigmaend1*lam+(1.0-lam)*sigmastart1)+(sigmaend2*lam+(1.0-lam)*sigmastart2))");

            custom_force_field->setCutoffDistance(converted_cutoff_distance);

            custom_force_field->addGlobalParameter("lam", Alchemical_value);
            custom_force_field->addGlobalParameter("delta", shift_delta);
            custom_force_field->addGlobalParameter("n", coulomb_power);
            custom_force_field->addGlobalParameter("krf", kvalue);
            custom_force_field->addGlobalParameter("crf", cvalue);
            custom_force_field->addGlobalParameter("cutoff", converted_cutoff_distance);
            custom_force_field->addGlobalParameter("SPOnOff", 0.0);


            if (flag_cutoff == CUTOFFNONPERIODIC)
            {
                custom_force_field->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
            }
            else
            {
                custom_force_field->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffPeriodic);
            }



            //NO REACTION FIELD IS APPLIED TO 1-4 INTERACTIONS. If the scaling factor is one (Glycam ff) then the OpenMM potential energy is not equal to he Sire energy. This is caused by the application of the reaction field on the 14 pairs in Sire.



            custom_intra_14_todummy = new OpenMM::CustomBondForce("withinCutoff*(Hcs + Hls);"
                                                                  "withinCutoff=step(cutofftd-r);"
                                                                  "Hcs=138.935456*q_prod/sqrt(diff_cl+r^2);"
                                                                  "diff_cl=(1.0-lamtd)*0.01;"
                                                                  "Hls=4.0*eps_avg*(LJ*LJ-LJ);"
                                                                  "LJ=((sigma_avg*sigma_avg)/soft)^3;"
                                                                  "soft=(diff_lj*deltatd*sigma_avg+r*r);"
                                                                  "diff_lj=(1.0-lamtd)*0.1;"
                                                                  "eps_avg = sqrt((1-lamtd)*(1-lamtd)*eaend + lamtd*lamtd*eastart + lamtd*(1-lamtd)*emix);"
                                                                  "sigma_avg = (1-lamtd)*saend + lamtd*sastart;"
                                                                  "q_prod = (1-lamtd)*(1-lamtd)*qpend + lamtd*lamtd*qpstart + lamtd*(1-lamtd)*qmix");

            custom_intra_14_todummy->addGlobalParameter("lamtd", 1.0 - Alchemical_value);
            custom_intra_14_todummy->addGlobalParameter("deltatd", shift_delta);
            custom_intra_14_todummy->addGlobalParameter("ntd", coulomb_power);
            custom_intra_14_todummy->addGlobalParameter("cutofftd", converted_cutoff_distance);




            custom_intra_14_fromdummy = new OpenMM::CustomBondForce("withinCutoff*(Hcs + Hls);"
                                                                    "withinCutoff=step(cutofffd-r);"
                                                                    "Hcs=138.935456*q_prod/sqrt(diff_cl+r^2);"
                                                                    "diff_cl=(1.0-lamfd)*0.01;"
                                                                    "Hls=4.0*eps_avg*(LJ*LJ-LJ);"
                                                                    "LJ=((sigma_avg*sigma_avg)/soft)^3;"
                                                                    "soft=(diff_lj*deltafd*sigma_avg+r*r);"
                                                                    "diff_lj=(1.0-lamfd)*0.1;"
                                                                    "eps_avg = sqrt(lamfd*lamfd*eaend + (1-lamfd)*(1-lamfd)*eastart + lamfd*(1-lamfd)*emix);"
                                                                    "sigma_avg = lamfd*saend + (1-lamfd)*sastart;"
                                                                    "q_prod = lamfd*lamfd*qpend + (1-lamfd)*(1-lamfd)*qpstart + lamfd*(1-lamfd)*qmix");


            custom_intra_14_fromdummy->addGlobalParameter("lamfd", Alchemical_value);
            custom_intra_14_fromdummy->addGlobalParameter("deltafd", shift_delta);
            custom_intra_14_fromdummy->addGlobalParameter("nfd", coulomb_power);
            custom_intra_14_fromdummy->addGlobalParameter("cutofffd", converted_cutoff_distance);




            custom_intra_14_fromdummy_todummy = new OpenMM::CustomBondForce("withinCutoff*(Hcs + Hls);"
                                                                            "withinCutoff=step(cutoffftd-r);"
                                                                            "Hcs=138.935456*q_prod/sqrt(diff_cl+r^2);"
                                                                            "diff_cl=(1.0-lamFTD)*0.01;"
                                                                            "Hls=4.0*eps_avg*(LJ*LJ-LJ);"
                                                                            "LJ=((sigma_avg*sigma_avg)/soft)^3;"
                                                                            "soft=(diff_lj*deltaftd*sigma_avg+r*r);"
                                                                            "diff_lj=(1.0-lamFTD)*0.1;"
                                                                            "eps_avg = sqrt(lamftd*lamftd*eaend + (1-lamftd)*(1-lamftd)*eastart + lamftd*(1-lamftd)*emix);"
                                                                            "sigma_avg = lamftd*saend + (1-lamftd)*sastart;"
                                                                            "q_prod = lamftd*lamftd*qpend + (1-lamftd)*(1-lamftd)*qpstart + lamftd*(1-lamftd)*qmix;"
                                                                            "lamFTD = max(lamftd,1-lamftd)");

            custom_intra_14_fromdummy_todummy->addGlobalParameter("lamftd", Alchemical_value);
            custom_intra_14_fromdummy_todummy->addGlobalParameter("deltaftd", shift_delta);
            custom_intra_14_fromdummy_todummy->addGlobalParameter("nftd", coulomb_power);
            custom_intra_14_fromdummy_todummy->addGlobalParameter("cutoffftd", converted_cutoff_distance);



        }


        custom_intra_14_clj = new OpenMM::CustomBondForce("withinCutoff*(Hl+Hc);"
                                                          "withinCutoff=step(cutoffhd-r);"
                                                          "Hl=4*eps_avg*((sigma_avg/r)^12-(sigma_avg/r)^6);"
                                                          "Hc=138.935456*q_prod/r;"
                                                          "eps_avg = sqrt(lamhd*lamhd*eaend + (1-lamhd)*(1-lamhd)*eastart + lamhd*(1-lamhd)*emix);"
                                                          "sigma_avg = lamhd*saend + (1-lamhd)*sastart;"
                                                          "q_prod = lamhd*lamhd*qpend + (1-lamhd)*(1-lamhd)*qpstart + lamhd*(1-lamhd)*qmix");


        custom_intra_14_clj->addGlobalParameter("lamhd", Alchemical_value);
        custom_intra_14_clj->addGlobalParameter("cutoffhd", converted_cutoff_distance);


        //REACTION FIELD 14 IMPLEMENTATION FOR FUTURE USE
        /*custom_intra_14_clj = new OpenMM::CustomBondForce("(Hl+Hc);"
                                                          "Hl=4*eps_avg*((sigma_avg/r)^12-(sigma_avg/r)^6);"
                                                          "Hc = 138.935456 * q_prod*(1.0/r + krf*(r*r)-crf);"
                                                          "eps_avg = sqrt(lamhd*lamhd*eaend + (1-lamhd)*(1-lamhd)*eastart + lamhd*(1-lamhd)*emix);"
                                                          "sigma_avg = lamhd*saend + (1-lamhd)*sastart;"
                                                          "q_prod = lamhd*lamhd*qpend + (1-lamhd)*(1-lamhd)*qpstart + lamhd*(1-lamhd)*qmix");


        custom_intra_14_clj->addGlobalParameter("lamhd",Alchemical_value);
        custom_intra_14_clj->addGlobalParameter("krf",kvalue);
        custom_intra_14_clj->addGlobalParameter("crf",cvalue);*/

        if (Debug)
        {
            qDebug() << "\nCut off type = " << CutoffType;
            qDebug() << "CutOff distance = " << converted_cutoff_distance << " Nm";
            qDebug() << "Dielectric constant = " << field_dielectric;
            qDebug() << "Lambda = " << Alchemical_value << " Coulomb Power = " << coulomb_power << " Delta Shift = " << shift_delta;

        }

    }

    // Andersen thermostat
    if (Andersen_flag == true)
    {
        const double converted_Temperature = convertTo(Temperature.value(), kelvin);

        OpenMM::AndersenThermostat * thermostat = new OpenMM::AndersenThermostat(converted_Temperature, Andersen_frequency);

        //Set The random seed
        thermostat->setRandomNumberSeed(random_seed);

        system_openmm->addForce(thermostat);

        if (Debug)
        {
            qDebug() << "\nAndersen Thermostat set\n";
            qDebug() << "Temperature = " << converted_Temperature << " K\n";
            qDebug() << "Frequency collisions = " << Andersen_frequency << " 1/ps\n";
        }
    }


    // Monte Carlo Barostat
    if (MCBarostat_flag == true)
    {
        const double converted_Temperature = convertTo(Temperature.value(), kelvin);
        const double converted_Pressure = convertTo(Pressure.value(), bar);

        OpenMM::MonteCarloBarostat * barostat = new OpenMM::MonteCarloBarostat(converted_Pressure, converted_Temperature, MCBarostat_frequency);

        //Set The random seed
        barostat->setRandomNumberSeed(random_seed);

        system_openmm->addForce(barostat);

        if (Debug)
        {
            qDebug() << "\nMonte Carlo Barostat set\n";
            qDebug() << "Temperature = " << converted_Temperature << " K\n";
            qDebug() << "Pressure = " << converted_Pressure << " bar\n";
            qDebug() << "Frequency every " << MCBarostat_frequency << " steps\n";
        }

    }
    /*******************************************************BONDED INTERACTIONS******************************************************/


    OpenMM::HarmonicBondForce * bondStretch_openmm = new OpenMM::HarmonicBondForce();

    OpenMM::HarmonicAngleForce * bondBend_openmm = new OpenMM::HarmonicAngleForce();

    OpenMM::PeriodicTorsionForce * bondTorsion_openmm = new OpenMM::PeriodicTorsionForce();


    OpenMM::CustomBondForce* solute_bond_perturbation = NULL;

    OpenMM::CustomAngleForce* solute_angle_perturbation = NULL;


    solute_bond_perturbation = new OpenMM::CustomBondForce("0.5*B*(r-req)^2;"
                                                           "B=bend*lambond+(1.0-lambond)*bstart;"
                                                           "req=rend*lambond+(1.0-lambond)*rstart");


    solute_bond_perturbation->addGlobalParameter("lambond", Alchemical_value);

    solute_angle_perturbation = new OpenMM::CustomAngleForce("0.5*A*(theta-thetaeq)^2;"
                                                             "A=aend*lamangle+(1.0-lamangle)*astart;"
                                                             "thetaeq=thetaend*lamangle+(1.0-lamangle)*thetastart");

    solute_angle_perturbation->addGlobalParameter("lamangle", Alchemical_value);



    /************************************************************RESTRAINTS********************************************************/

    OpenMM::CustomExternalForce * positionalRestraints_openmm = NULL;

    if (Restraint_flag == true)
    {

        positionalRestraints_openmm = new OpenMM::CustomExternalForce("k*d2;"
                                                                      "d2 = max(0.0, d1 - d^2);"
                                                                      "d1 = (x-xref)^2 + (y-yref)^2  + (z-zref)^2");

        positionalRestraints_openmm->addPerParticleParameter("xref");
        positionalRestraints_openmm->addPerParticleParameter("yref");
        positionalRestraints_openmm->addPerParticleParameter("zref");
        positionalRestraints_openmm->addPerParticleParameter("k");
        positionalRestraints_openmm->addPerParticleParameter("d");

        system_openmm->addForce(positionalRestraints_openmm);

        if (Debug)
            qDebug() << "\n\nRestraint is ON\n\n";
    }

    /****************************************BOND LINK POTENTIAL*****************************/
    /* !! CustomBondForce does not (OpenMM 6.2) apply PBC checks so code will be buggy is restraints involve one atom that diffuses
       out of the box. */

    OpenMM::CustomBondForce * custom_link_bond = new OpenMM::CustomBondForce("kl*max(0,d-dl*dl);"
                                                                             "d=(r-reql)*(r-reql)");
    custom_link_bond->addPerBondParameter("reql");
    custom_link_bond->addPerBondParameter("kl");
    custom_link_bond->addPerBondParameter("dl");

    //OpenMM vector coordinate
    std::vector<OpenMM::Vec3> positions_openmm(nats);

    //OpenMM vector momenta
    std::vector<OpenMM::Vec3> velocities_openmm(nats);

    std::vector<std::pair<int, int> > bondPairs;

    int system_index = 0;

    // To avoid possible mismatch between the index in which atoms are added to the openmm system arrays and
    // their atomic numbers in sire, one array is populated while filling up the openmm global arrays
    //  AtomNumtoopenmmIndex
    QHash<int, int> AtomNumToOpenMMIndex;


    for (int i = 0; i < nmols; ++i)
    {
        
        const int nats_mol = ws.nAtoms(i);

        const double *m = ws.massArray(i);

        MolNum molnum = moleculegroup.molNumAt(i);

        const ViewsOfMol &molview = moleculegroup[molnum].data();

        const Molecule &mol = molview.molecule();

        Selector<Atom> molatoms = mol.atoms();

        for (int j = 0; j < nats_mol; ++j)
        {
            /*JM 10/16 make sure that perturbed atoms have mass of heaviest end-state */
            system_openmm->addParticle(m[j]);

            Atom at = molatoms(j);
            AtomNum atnum = at.number();

            if (Debug)
                qDebug() << " openMM_index " << system_index << " Sire Atom Number " << atnum.toString() << " Mass particle = " << m[j];

            AtomNumToOpenMMIndex[atnum.value()] = system_index;

            // JM Nov 12
            // The code below implements a ThreeParticleAverageSite for virtual sites for EPW atoms present in a WAT residue
            // This is very AMBER specific.

            AtomName atname = at.name();

            if (Debug)
                qDebug() << " atname " << atname.value() << " mol " << i;

            if (atname == AtomName("EPW"))
            {

                ResName resname = at.residue().name();

                if (resname == ResName("WAT"))
                {

                    Atom oatom = molatoms.select(AtomName("O"));
                    Atom h1atom = molatoms.select(AtomName("H1"));
                    Atom h2atom = molatoms.select(AtomName("H2"));

                    AmberParameters amber_params = mol.property("amberparameters").asA<AmberParameters>();
                    QList<BondID> bonds_ff = amber_params.getAllBonds();

                    double distoh = -1.0;
                    double disthh = -1.0;
                    double distoe = -1.0;

                    for (int k = 0; k < bonds_ff.length(); k++)
                    {

                        BondID bond_ff = bonds_ff[k];
                        QList<double> bond_params = amber_params.getParams(bond_ff);

                        double r0 = bond_params[1];

                        AtomName at0name = mol.select(bond_ff.atom0()).name();
                        AtomName at1name = mol.select(bond_ff.atom1()).name();

                        // qDebug() << " at0name " << at0name.toString() << " at1name " << at1name.toString();

                        if ((at0name == AtomName("O") and at1name == AtomName("H1")) or ( at0name == AtomName("H1") and at1name == AtomName("O")))
                        {
                            distoh = r0;
                        }
                        else if ((at0name == AtomName("H1") and at1name == AtomName("H2")) or ( at0name == AtomName("H2") and at1name == AtomName("H1")))
                        {
                            disthh = r0;
                        }
                        else if ((at0name == AtomName("EPW") and at1name == AtomName("O")) or ( at0name == AtomName("O") and at1name == AtomName("EPW")))
                        {
                            distoe = r0;
                        }
                    }

                    if (distoh < 0 or disthh < 0 or distoe < 0)
                    {
                        throw SireError::program_bug(QObject::tr(
                                                                 "Could not find expected atoms in TIP4P water molecule."), CODELOC);
                    }

                    //qDebug() << " distoe " << distoe << " distoh " << distoh << " disthh " << disthh;

                    double weightH = distoe / sqrt((distoh * distoh) - (0.25 * disthh * disthh));

                    int o_index = AtomNumToOpenMMIndex[oatom.number().value()];
                    int h1_index = AtomNumToOpenMMIndex[h1atom.number().value()];
                    int h2_index = AtomNumToOpenMMIndex[h2atom.number().value()];

                    if (Debug)
                        qDebug() << "virtual site " << system_index << " o " << o_index << " h1 " << h1_index << " h2 " << h2_index << " 1 - weightH " << 1 - weightH << " weightH/2 " << weightH / 2;

                    OpenMM::ThreeParticleAverageSite * vsite = new OpenMM::ThreeParticleAverageSite(o_index, h1_index, h2_index, 1 - weightH, weightH / 2, weightH / 2);

                    system_openmm->setVirtualSite(system_index, vsite);

                }
            }

            system_index = system_index + 1;

        }// end of loop on atoms in molecule

    }//end of loop on molecules in workspace


    int num_atoms_till_i = 0;

    /*NON BONDED PER PARTICLE PARAMETERS*/

    custom_force_field->addPerParticleParameter("qstart");
    custom_force_field->addPerParticleParameter("qend");
    custom_force_field->addPerParticleParameter("epstart");
    custom_force_field->addPerParticleParameter("epend");
    custom_force_field->addPerParticleParameter("sigmastart");
    custom_force_field->addPerParticleParameter("sigmaend");
    custom_force_field->addPerParticleParameter("isHD");
    custom_force_field->addPerParticleParameter("isTD");
    custom_force_field->addPerParticleParameter("isFD");
    custom_force_field->addPerParticleParameter("isSolvent");

    custom_intra_14_clj->addPerBondParameter("qpstart");
    custom_intra_14_clj->addPerBondParameter("qpend");
    custom_intra_14_clj->addPerBondParameter("qmix");
    custom_intra_14_clj->addPerBondParameter("eastart");
    custom_intra_14_clj->addPerBondParameter("eaend");
    custom_intra_14_clj->addPerBondParameter("emix");
    custom_intra_14_clj->addPerBondParameter("sastart");
    custom_intra_14_clj->addPerBondParameter("saend");

    custom_intra_14_todummy->addPerBondParameter("qpstart");
    custom_intra_14_todummy->addPerBondParameter("qpend");
    custom_intra_14_todummy->addPerBondParameter("qmix");
    custom_intra_14_todummy->addPerBondParameter("eastart");
    custom_intra_14_todummy->addPerBondParameter("eaend");
    custom_intra_14_todummy->addPerBondParameter("emix");
    custom_intra_14_todummy->addPerBondParameter("sastart");
    custom_intra_14_todummy->addPerBondParameter("saend");

    custom_intra_14_fromdummy->addPerBondParameter("qpstart");
    custom_intra_14_fromdummy->addPerBondParameter("qpend");
    custom_intra_14_fromdummy->addPerBondParameter("qmix");
    custom_intra_14_fromdummy->addPerBondParameter("eastart");
    custom_intra_14_fromdummy->addPerBondParameter("eaend");
    custom_intra_14_fromdummy->addPerBondParameter("emix");
    custom_intra_14_fromdummy->addPerBondParameter("sastart");
    custom_intra_14_fromdummy->addPerBondParameter("saend");

    custom_intra_14_fromdummy_todummy->addPerBondParameter("qpstart");
    custom_intra_14_fromdummy_todummy->addPerBondParameter("qpend");
    custom_intra_14_fromdummy_todummy->addPerBondParameter("qmix");
    custom_intra_14_fromdummy_todummy->addPerBondParameter("eastart");
    custom_intra_14_fromdummy_todummy->addPerBondParameter("eaend");
    custom_intra_14_fromdummy_todummy->addPerBondParameter("emix");
    custom_intra_14_fromdummy_todummy->addPerBondParameter("sastart");
    custom_intra_14_fromdummy_todummy->addPerBondParameter("saend");


    /*BONDED PER PARTICLE PARAMETERS*/

    solute_bond_perturbation->addPerBondParameter("bstart");
    solute_bond_perturbation->addPerBondParameter("bend");
    solute_bond_perturbation->addPerBondParameter("rstart");
    solute_bond_perturbation->addPerBondParameter("rend");

    solute_angle_perturbation->addPerAngleParameter("astart");
    solute_angle_perturbation->addPerAngleParameter("aend");
    solute_angle_perturbation->addPerAngleParameter("thetastart");
    solute_angle_perturbation->addPerAngleParameter("thetaend");


    // JM July 13. This also needs to be changed because there could be more than one perturbed molecule
    //Molecule solutemol = solute.moleculeAt(0).molecule();


    int nions = 0;

    QVector<bool> perturbed_energies_tmp(8);

    for (int i = 0; i < perturbed_energies_tmp.size(); i++)
        perturbed_energies_tmp[i] = false;


    // The default 1,4 scaling factors
    double const Coulomb14Scale = 1.0 / 1.2;
    double const LennardJones14Scale = 1.0 / 2.0;

    // A list of 1,4 atom pairs with non default scale factors
    // for each entry, first pair has pair of indices, second has pair of scale factors
    //QList< QPair< QPair<int,int>, QPair<double, double > > > custom14pairs;
    QHash< QPair<int, int>, QPair<double, double> > custom14pairs;

    bool special_14 = false;

    for (int i = 0; i < nmols; i++)
    {

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

        if (molecule.hasProperty("perturbations"))
        {

            if (Debug)
                qDebug() << "Molecule Perturbed number = " << i;

            AtomCharges atomcharges_start = molecule.property("initial_charge").asA<AtomCharges>();
            AtomCharges atomcharges_final = molecule.property("final_charge").asA<AtomCharges>();

            start_charges = atomcharges_start.toVector();
            final_charges = atomcharges_final.toVector();

            AtomLJs atomvdws_start = molecule.property("initial_LJ").asA<AtomLJs>();
            AtomLJs atomvdws_final = molecule.property("final_LJ").asA<AtomLJs>();
            start_LJs = atomvdws_start.toVector();
            final_LJs = atomvdws_final.toVector();
        }

        for (int j = 0; j < ljparameters.size(); j++)
        {

            double sigma = ljparameters[j].sigma();
            double epsilon = ljparameters[j].epsilon();
            double charge = charges[j].value();

            nonbond_openmm->addParticle(charge, sigma * OpenMM::NmPerAngstrom, epsilon * OpenMM::KJPerKcal);

            Atom atom = molecule.molecule().atoms()(j);

            if (molecule.hasProperty("perturbations"))
            {

                // Is atom a hard, from dummy or to dummy type?
                bool ishard = false;
                bool istodummy = false;
                bool isfromdummy = false;

                for (int l = 0; l < solutehard.nViews(); l++)
                {

                    Selector<Atom> view_atoms = solutehard.viewAt(l).atoms();

                    for (int m = 0; m < view_atoms.count(); m++)
                    {

                        Atom view_atom = view_atoms(m);

                        if (atom == view_atom)
                        {
                            ishard = true;
                            break;
                        }
                    }//end for

                    if (ishard)
                        break;
                }//end for

                // if not hard check if to_dummy
                if (!ishard)
                {

                    for (int l = 0; l < solutetodummy.nViews(); l++)
                    {

                        Selector<Atom> view_atoms = solutetodummy.viewAt(l).atoms();

                        for (int m = 0; m < view_atoms.count(); m++)
                        {

                            Atom view_atom = view_atoms(m);

                            if (atom == view_atom)
                            {
                                istodummy = true;
                                break;
                            }
                        }//end for
                        if (istodummy)
                            break;
                    }//end for
                }

                // if not todummy, check if fromdummy
                if (!istodummy && !ishard)
                {

                    for (int l = 0; l < solutefromdummy.nViews(); l++)
                    {

                        Selector<Atom> view_atoms = solutefromdummy.viewAt(l).atoms();

                        for (int m = 0; m < view_atoms.count(); m++)
                        {

                            Atom view_atom = view_atoms(m);

                            if (atom == view_atom)
                            {
                                isfromdummy = true;
                                break;
                            }
                        }//end for
                        if (isfromdummy)
                            break;
                    }//end for
                }

                if (ishard)
                {//hard solute atom

                    double charge_start = start_charges[j].value();
                    double charge_final = final_charges[j].value();

                    double epsilon_start = start_LJs[j].epsilon();
                    double epsilon_final = final_LJs[j].epsilon();
                    double sigma_start = start_LJs[j].sigma();
                    double sigma_final = final_LJs[j].sigma();

                    custom_non_bonded_params[0] = charge_start;
                    custom_non_bonded_params[1] = charge_final;
                    custom_non_bonded_params[2] = epsilon_start * OpenMM::KJPerKcal;
                    custom_non_bonded_params[3] = epsilon_final * OpenMM::KJPerKcal;
                    custom_non_bonded_params[4] = sigma_start * OpenMM::NmPerAngstrom;
                    custom_non_bonded_params[5] = sigma_final * OpenMM::NmPerAngstrom;
                    custom_non_bonded_params[6] = 1.0; //isHard
                    custom_non_bonded_params[7] = 0.0; //isTodummy
                    custom_non_bonded_params[8] = 0.0; //isFromdummy
                    custom_non_bonded_params[9] = 0.0; //isSolventProtein

                    if (Debug)
                        qDebug() << "hard solute = " << atom.index();
                }
                    // JM July 13 THIS NEEDS FIXING TO DEAL WITH GROUPS THAT CONTAIN MORE THAN ONE MOLECULE
                else if (istodummy)
                {//to dummy solute atom

                    double charge_start = start_charges[j].value();
                    double charge_final = final_charges[j].value();
                    double epsilon_start = start_LJs[j].epsilon();
                    double epsilon_final = final_LJs[j].epsilon();
                    double sigma_start = start_LJs[j].sigma();
                    double sigma_final = final_LJs[j].sigma();

                    custom_non_bonded_params[0] = charge_start;
                    custom_non_bonded_params[1] = charge_final;
                    custom_non_bonded_params[2] = epsilon_start * OpenMM::KJPerKcal;
                    custom_non_bonded_params[3] = epsilon_final * OpenMM::KJPerKcal;
                    custom_non_bonded_params[4] = sigma_start * OpenMM::NmPerAngstrom;
                    custom_non_bonded_params[5] = sigma_final * OpenMM::NmPerAngstrom;
                    custom_non_bonded_params[6] = 0.0; //isHard
                    custom_non_bonded_params[7] = 1.0; //isTodummy
                    custom_non_bonded_params[8] = 0.0; //isFromdummy
                    custom_non_bonded_params[9] = 0.0; //isSolventProtein

                    if (Debug)
                        qDebug() << "to dummy solute = " << atom.index();
                }
                else if (isfromdummy)
                {//from dummy solute atom

                    double charge_start = start_charges[j].value();
                    double charge_final = final_charges[j].value();
                    double epsilon_start = start_LJs[j].epsilon();
                    double epsilon_final = final_LJs[j].epsilon();
                    double sigma_start = start_LJs[j].sigma();
                    double sigma_final = final_LJs[j].sigma();

                    custom_non_bonded_params[0] = charge_start;
                    custom_non_bonded_params[1] = charge_final;
                    custom_non_bonded_params[2] = epsilon_start * OpenMM::KJPerKcal;
                    custom_non_bonded_params[3] = epsilon_final * OpenMM::KJPerKcal;
                    custom_non_bonded_params[4] = sigma_start * OpenMM::NmPerAngstrom;
                    custom_non_bonded_params[5] = sigma_final * OpenMM::NmPerAngstrom;
                    custom_non_bonded_params[6] = 0.0; //isHard
                    custom_non_bonded_params[7] = 0.0; //isTodummy
                    custom_non_bonded_params[8] = 1.0; //isFromdummy
                    custom_non_bonded_params[9] = 0.0; //isSolventProtein

                    if (Debug)
                        qDebug() << "from dummy solute = " << atom.index();
                }

                else
                {//What if some atoms were not perturbed at all in the pert file? Use default params

                    custom_non_bonded_params[0] = charge;
                    custom_non_bonded_params[1] = charge;
                    custom_non_bonded_params[2] = epsilon * OpenMM::KJPerKcal;
                    custom_non_bonded_params[3] = epsilon * OpenMM::KJPerKcal;
                    custom_non_bonded_params[4] = sigma * OpenMM::NmPerAngstrom;
                    custom_non_bonded_params[5] = sigma * OpenMM::NmPerAngstrom;
                    custom_non_bonded_params[6] = 1.0; //isHard
                    custom_non_bonded_params[7] = 0.0; //isTodummy
                    custom_non_bonded_params[8] = 0.0; //isFromdummy
                    custom_non_bonded_params[9] = 0.0; //isSolventProtein

                    if (Debug)
                        qDebug() << " unperturbed solute atom " << atom.index();
                }
            }//end if perturbation section
            else
            {//solvent atom like hard

                custom_non_bonded_params[0] = charge;
                custom_non_bonded_params[1] = charge;
                custom_non_bonded_params[2] = epsilon * OpenMM::KJPerKcal;
                custom_non_bonded_params[3] = epsilon * OpenMM::KJPerKcal;
                custom_non_bonded_params[4] = sigma * OpenMM::NmPerAngstrom;
                custom_non_bonded_params[5] = sigma * OpenMM::NmPerAngstrom;
                custom_non_bonded_params[6] = 1.0; //isHard
                custom_non_bonded_params[7] = 0.0; //isTodummy
                custom_non_bonded_params[8] = 0.0; //isFromdummy
                custom_non_bonded_params[9] = 1.0; //isSolventProtein

                if (Debug)
                    qDebug() << "Solvent = " << atom.index();

            }

            if (Debug)
            {
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

            custom_force_field->addParticle(custom_non_bonded_params);

        }


        /****************************************************RESTRAINTS*******************************************************/

        if (Restraint_flag == true)
        {

            bool hasRestrainedAtoms = molecule.hasProperty("restrainedatoms");

            if (hasRestrainedAtoms)
            {

                Properties restrainedAtoms = molecule.property("restrainedatoms").asA<Properties>();

                int nrestrainedatoms = restrainedAtoms.property(QString("nrestrainedatoms")).asA<VariantProperty>().toInt();

                if (Debug)
                    qDebug() << "nrestrainedatoms = " << nrestrainedatoms;

                for (int i = 0; i < nrestrainedatoms; i++)
                {

                    int atomnum = restrainedAtoms.property(QString("AtomNum(%1)").arg(i)).asA<VariantProperty>().toInt();
                    double xref = restrainedAtoms.property(QString("x(%1)").arg(i)).asA<VariantProperty>().toDouble();
                    double yref = restrainedAtoms.property(QString("y(%1)").arg(i)).asA<VariantProperty>().toDouble();
                    double zref = restrainedAtoms.property(QString("z(%1)").arg(i)).asA<VariantProperty>().toDouble();
                    double k = restrainedAtoms.property(QString("k(%1)").arg(i)).asA<VariantProperty>().toDouble();
                    double d = restrainedAtoms.property(QString("d(%1)").arg(i)).asA<VariantProperty>().toDouble();

                    int openmmindex = AtomNumToOpenMMIndex[atomnum];

                    if (Debug)
                    {
                        qDebug() << "atomnum " << atomnum << " openmmindex " << openmmindex << " x " << xref << " y " << yref << " z " << zref << " k " << k << " d " << d;
                    }

                    int posrestrdim = 5;
                    std::vector<double> params(posrestrdim);

                    params[0] = xref * OpenMM::NmPerAngstrom;
                    params[1] = yref * OpenMM::NmPerAngstrom;
                    params[2] = zref * OpenMM::NmPerAngstrom;
                    params[3] = k * (OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
                    params[4] = d * OpenMM::NmPerAngstrom;

                    positionalRestraints_openmm->addParticle(openmmindex, params);
                }
            }
        }//end of restraint flag


        // IONS


        bool hasConnectivity = molecule.hasProperty("connectivity");

        if (!hasConnectivity)
        {

            num_atoms_till_i = num_atoms_till_i + num_atoms_molecule;

            if (Debug)
            {
                qDebug() << "\nAtoms = " << num_atoms_molecule << " Num atoms till i =" << num_atoms_till_i << "\n";
                qDebug() << "\n*********************MONOATOMIC MOLECULE DETECTED**************************\n";
            }
            nions = nions + 1;
            continue;
        }

        //BONDED TERMS

        QList< BondID > bond_pert_list;
        QList< BondID > bond_pert_swap_list;
        QList< AngleID > angle_pert_list;
        QList< AngleID > angle_pert_swap_list;
        QList< DihedralID > dihedral_pert_list;
        QList< DihedralID > dihedral_pert_swap_list;
        QList< ImproperID > improper_pert_list;
        QList< ImproperID > improper_pert_swap_list;


        double HMASS = 1.10;/* g per mol-1*/
        //double HEAVYH=12.0;/* g per mol-1*/
        double SMALL = 0.0001;

        if (solute.contains(molecule))
        {
            Perturbations pert_params = molecule.property("perturbations").asA<Perturbations>();

            QList< PropPtr<Perturbation> > perturbation_list = pert_params.perturbations();

            std::vector<double> solute_bond_perturbation_params(4);
            std::vector<double> solute_angle_perturbation_params(4);
            std::vector<double> solute_torsion_perturbation_params(1);

            QHash<BondID, double> bond_pert_eq_list;

            for (QList< PropPtr<Perturbation> >::const_iterator it = perturbation_list.constBegin(); it != perturbation_list.constEnd(); ++it)
            {

                const Perturbation &pert = *it;

                if (pert.isA<InternalPerturbation>())
                {

                    QString str = pert.what();

                    if (str == "SireMM::TwoAtomPerturbation")
                    {
                        const TwoAtomPerturbation &two = pert.asA<TwoAtomPerturbation>();
                        int idx0 = two.atom0().asA<AtomIdx>().value() + num_atoms_till_i;
                        int idx1 = two.atom1().asA<AtomIdx>().value() + num_atoms_till_i;
                        double rstart = two.initialForms()[Symbol("r0")].toString().toDouble();
                        double bstart = two.initialForms()[Symbol("k")].toString().toDouble();
                        double rend = two.finalForms()[Symbol("r0")].toString().toDouble();
                        double bend = two.finalForms()[Symbol("k")].toString().toDouble();

                        solute_bond_perturbation_params[0] = bstart * 2.0 * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm;
                        solute_bond_perturbation_params[1] = bend * 2.0 * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm;
                        solute_bond_perturbation_params[2] = rstart * OpenMM::NmPerAngstrom;
                        solute_bond_perturbation_params[3] = rend * OpenMM::NmPerAngstrom;

                        /* JM 10/16 Also apply this if 'no solute constraints' flag is on*/
                        if (flag_constraint == NONE)
                        {
                            solute_bond_perturbation->addBond(idx0, idx1, solute_bond_perturbation_params);
                        }
                        else if (flag_constraint == ALLBONDS || flag_constraint == HANGLES)
                        {
                            /* JM 10/16 ALLBONDS and HANGLES may be unwise with current free energy implementation !*/
                            double pert_eq_distance = solute_bond_perturbation_params[3] * Alchemical_value + (1.0 - Alchemical_value) * solute_bond_perturbation_params[2];
                            system_openmm->addConstraint(idx0, idx1, pert_eq_distance);
                            bond_pert_eq_list.insert(BondID(two.atom0(), two.atom1()), pert_eq_distance * OpenMM::AngstromsPerNm);
                            if (Debug)
                            {
                                qDebug() << "bond start distance = " << solute_bond_perturbation_params[2] << " Nm";
                                qDebug() << "bond end distance = " << solute_bond_perturbation_params[3] << " Nm";
                                qDebug() << "Perturbation bond equilibrium distance = " << pert_eq_distance << " Nm";
                            }
                        }
                        /* JM 10/16 */
                        /*  Here add code to constraint hbonds only if initial and final parameters are unperturbed*/
                        /*  check also what is the mass of the atoms in that case */
                        else if (flag_constraint == HBONDS and flag_noperturbedconstraints)
                        {
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
                          if (Debug)
                          {
                              qDebug() << " m0 " << m0 << " m1 " << m1 << "\n";
                              qDebug() << " deltar " << deltar << " " << " deltak " << deltak;
                          }
                          /* bonds that do not change parameters are constrained*/
                          double pert_eq_distance = solute_bond_perturbation_params[3] * Alchemical_value + (1.0 - Alchemical_value) * solute_bond_perturbation_params[2];
                          if (deltar < SMALL and deltak < SMALL)
                          {
                              system_openmm->addConstraint(idx0, idx1, pert_eq_distance);
                              if (Debug)
                              {
                                  qDebug() << "perturbed bond but no parameter changes so constrained " << atom0.name().toString()
                                           << "-" << atom1.name().toString() << "\n";
                              }
                          }
                          /* bonds that change parameters and have one of the atoms with a mass < HMASS are constrained*/
                          else if (m0 < HMASS or m1 < HMASS)
                          {
                              system_openmm->addConstraint(idx0, idx1, pert_eq_distance);
                              if (Debug)
                              {
                                  qDebug() << "perturbed bond parameter changes but involving " 
                                           << " light mass so constrained " << atom0.name().toString()
                                           << "- " << atom1.name().toString() << "\n";
                              }
                          }
                          /* other bonds are flexible */
                          else
                          {
                              solute_bond_perturbation->addBond(idx0, idx1, solute_bond_perturbation_params);
                               if (Debug)
                               {
                                   qDebug() << "perturbed bond flexible " << atom0.name().toString()
                                            << "- " << atom1.name().toString() << "\n"; 
                               }
                          }
                        }
                        else if (flag_constraint == HBONDS)
                        {
                            const SireMol::Atom atom0 = molecule.select(two.atom0());
                            QString initial_type_atom0 = atom0.property<QString>("initial_ambertype");
                            QString final_type_atom0 = atom0.property<QString>("final_ambertype");

                            const SireMol::Atom atom1 = molecule.select(two.atom1());
                            QString initial_type_atom1 = atom1.property<QString>("initial_ambertype");
                            QString final_type_atom1 = atom1.property<QString>("final_ambertype");

                            if (initial_type_atom0.startsWith("h", Qt::CaseInsensitive) || final_type_atom0.startsWith("h", Qt::CaseInsensitive) ||
                                initial_type_atom1.startsWith("h", Qt::CaseInsensitive) || final_type_atom1.startsWith("h", Qt::CaseInsensitive))
                            {
                                double pert_eq_distance = solute_bond_perturbation_params[3] * Alchemical_value + (1.0 - Alchemical_value) * solute_bond_perturbation_params[2];
                                system_openmm->addConstraint(idx0, idx1, pert_eq_distance);

                                if (Debug)
                                {
                                    qDebug() << "Two/one bond atom(s) start(s) or end(s) with h/H";
                                    qDebug() << "bond start distance = " << solute_bond_perturbation_params[2] << " Nm";
                                    qDebug() << "bond end distance = " << solute_bond_perturbation_params[3] << " Nm";
                                    qDebug() << "Perturbation bond equilibrium distance = " << pert_eq_distance << " Nm";
                                }
                            }
                            else
                            {
                                solute_bond_perturbation->addBond(idx0, idx1, solute_bond_perturbation_params);
                            }

                            if (Debug)
                            {
                                qDebug() << "Atom0 initil type = " << initial_type_atom0;
                                qDebug() << "Atom0 final type = " << final_type_atom0;
                                qDebug() << "Atom1 initil type = " << initial_type_atom1;
                                qDebug() << "Atom1 final type = " << final_type_atom1;
                            }

                        }

                        bond_pert_list.append(BondID(two.atom0(), two.atom1()));
                        bond_pert_swap_list.append(BondID(two.atom1(), two.atom0()));

                        bondPairs.push_back(std::make_pair(idx0, idx1));

                        if (Debug)
                        {
                            qDebug() << "Atom0 = " << two.atom0().asA<AtomIdx>().value() <<
                                "Atom1 = " << two.atom1().asA<AtomIdx>().value();
                            qDebug() << "IDX0 = " << idx0 << "IDX1 = " << idx1 << "\n";
                            qDebug() << "rstart = " << rstart << " A" <<
                                "rend = " << rend << " A";
                            qDebug() << "bstart = " << bstart << " kcal/A A" <<
                                "bend = " << bend << " kcal/A A" << "\n";
                        }
                    }
                    if (str == "SireMM::ThreeAtomPerturbation")
                    {
                        const ThreeAtomPerturbation &three = pert.asA<ThreeAtomPerturbation>();
                        int idx0 = three.atom0().asA<AtomIdx>().value() + num_atoms_till_i;
                        int idx1 = three.atom1().asA<AtomIdx>().value() + num_atoms_till_i;
                        int idx2 = three.atom2().asA<AtomIdx>().value() + num_atoms_till_i;
                        double astart = three.initialForms()[Symbol("k")].toString().toDouble();
                        double thetastart = three.initialForms()[Symbol("theta0")].toString().toDouble();
                        double aend = three.finalForms()[Symbol("k")].toString().toDouble();
                        double thetaend = three.finalForms()[Symbol("theta0")].toString().toDouble();

                        solute_angle_perturbation_params[0] = astart * 2.0 * OpenMM::KJPerKcal;
                        solute_angle_perturbation_params[1] = aend * 2.0 * OpenMM::KJPerKcal;
                        solute_angle_perturbation_params[2] = thetastart;
                        solute_angle_perturbation_params[3] = thetaend;

                        if (Debug)
                        {
                            qDebug() << "astart = " << solute_angle_perturbation_params[0] << " kJ/rad rad" <<
                                " aend = " << solute_angle_perturbation_params[1] << " kJ/rad rad";
                            qDebug() << "thetastart = " << solute_angle_perturbation_params[2] << " rad" <<
                                "thetaend = " << solute_angle_perturbation_params[3] << " rad";
                        }

                        if (flag_constraint == HANGLES)
                        {
                            const SireMol::Atom atom0 = molecule.select(three.atom0());
                            QString initial_type_atom0 = atom0.property<QString>("initial_ambertype");
                            QString final_type_atom0 = atom0.property<QString>("final_ambertype");

                            const SireMol::Atom atom1 = molecule.select(three.atom1());
                            QString initial_type_atom1 = atom1.property<QString>("initial_ambertype");
                            QString final_type_atom1 = atom1.property<QString>("final_ambertype");

                            const SireMol::Atom atom2 = molecule.select(three.atom2());
                            QString initial_type_atom2 = atom2.property<QString>("initial_ambertype");
                            QString final_type_atom2 = atom2.property<QString>("final_ambertype");

                            bool H_X_H = (initial_type_atom0.startsWith("h", Qt::CaseInsensitive) || final_type_atom0.startsWith("h", Qt::CaseInsensitive)) &&
                                (initial_type_atom2.startsWith("h", Qt::CaseInsensitive) || final_type_atom2.startsWith("h", Qt::CaseInsensitive));

                            bool H_O_X = (initial_type_atom0.startsWith("h", Qt::CaseInsensitive) || final_type_atom0.startsWith("h", Qt::CaseInsensitive)) &&
                                (initial_type_atom1.startsWith("o", Qt::CaseInsensitive) || final_type_atom1.startsWith("o", Qt::CaseInsensitive));

                            bool X_O_H = (initial_type_atom1.startsWith("o", Qt::CaseInsensitive) || final_type_atom1.startsWith("o", Qt::CaseInsensitive)) &&
                                (initial_type_atom2.startsWith("h", Qt::CaseInsensitive) || final_type_atom2.startsWith("h", Qt::CaseInsensitive));

                            if (Debug)
                            {
                                if (H_X_H)
                                    qDebug() << "type =  H_X_H";
                                if (H_O_X)
                                    qDebug() << "type =  H_O_X";
                                if (X_O_H)
                                    qDebug() << "type =  X_O_H";
                            }

                            if (H_X_H || H_O_X || X_O_H)
                            {

                                const BondID * first_alchemical_bond = NULL;
                                const BondID * second_alchemical_bond = NULL;

                                double first_alchemical_distance = -1.0;
                                double second_alchemical_distance = -1.0;

                                if (bond_pert_eq_list.contains(BondID(three.atom0(), three.atom1())))
                                {
                                    first_alchemical_bond = &((bond_pert_eq_list.find(BondID(three.atom0(), three.atom1()))).key());
                                    first_alchemical_distance = bond_pert_eq_list.value(*first_alchemical_bond);
                                    if (Debug)
                                        qDebug() << "Atom0 - Atom1";
                                }
                                else if (bond_pert_eq_list.contains(BondID(three.atom1(), three.atom0())))
                                {
                                    first_alchemical_bond = &((bond_pert_eq_list.find(BondID(three.atom1(), three.atom0()))).key());
                                    first_alchemical_distance = bond_pert_eq_list.value(*first_alchemical_bond);
                                    if (Debug)
                                        qDebug() << "Atom1 - Atom0";
                                }
                                else
                                {

                                    if (Debug)
                                        qDebug() << "First perturbed bond was not foud in the perturned list";
                                    first_alchemical_bond = new BondID(three.atom0(), three.atom1());
                                }


                                if (bond_pert_eq_list.contains(BondID(three.atom1(), three.atom2())))
                                {
                                    second_alchemical_bond = &((bond_pert_eq_list.find(BondID(three.atom1(), three.atom2()))).key());
                                    second_alchemical_distance = bond_pert_eq_list.value(*second_alchemical_bond);
                                    if (Debug)
                                        qDebug() << "Atom1 - Atom2";
                                }
                                else if (bond_pert_eq_list.contains(BondID(three.atom2(), three.atom1())))
                                {
                                    second_alchemical_bond = &((bond_pert_eq_list.find(BondID(three.atom2(), three.atom1()))).key());
                                    second_alchemical_distance = bond_pert_eq_list.value(*second_alchemical_bond);
                                    if (Debug)
                                        qDebug() << "Atom2 - Atom1";
                                }
                                else
                                {
                                    if (Debug)
                                        qDebug() << "Second perturbed bond was not foud in the perturned list";
                                    second_alchemical_bond = new BondID(three.atom2(), three.atom1());
                                }


                                if (Debug)
                                    qDebug() << "First Alchemical distance = " << first_alchemical_distance
                                    << "Second Alchemical distance = " << second_alchemical_distance;

                                SireMaths::Vector bond1_vec;
                                SireMaths::Vector bond2_vec;

                                if (first_alchemical_bond->atom0() == second_alchemical_bond->atom0())
                                {

                                    SireMaths::Vector tmp1 = (molecule.atom(first_alchemical_bond->atom1())).property<SireMaths::Vector>("coordinates");
                                    SireMaths::Vector tmp2 = (molecule.atom(second_alchemical_bond->atom1())).property<SireMaths::Vector>("coordinates");
                                    SireMaths::Vector tmpc = (molecule.atom(first_alchemical_bond->atom0())).property<SireMaths::Vector>("coordinates");
                                    bond1_vec = tmp1 - tmpc;
                                    bond2_vec = tmp2 - tmpc;
                                    if (Debug)
                                        qDebug() << "Bond1 Atom0 = Bond2 Atom0";
                                }
                                else if (first_alchemical_bond->atom0() == second_alchemical_bond->atom1())
                                {
                                    SireMaths::Vector tmp1 = (molecule.atom(first_alchemical_bond->atom1())).property<SireMaths::Vector>("coordinates");
                                    SireMaths::Vector tmp2 = (molecule.atom(second_alchemical_bond->atom0())).property<SireMaths::Vector>("coordinates");
                                    SireMaths::Vector tmpc = (molecule.atom(first_alchemical_bond->atom0())).property<SireMaths::Vector>("coordinates");
                                    bond1_vec = tmp1 - tmpc;
                                    bond2_vec = tmp2 - tmpc;
                                    if (Debug)
                                        qDebug() << "Bond1 Atom0 = Bond2 Atom1";
                                }
                                else if (first_alchemical_bond->atom1() == second_alchemical_bond->atom0())
                                {
                                    SireMaths::Vector tmp1 = (molecule.atom(first_alchemical_bond->atom0())).property<SireMaths::Vector>("coordinates");
                                    SireMaths::Vector tmp2 = (molecule.atom(second_alchemical_bond->atom1())).property<SireMaths::Vector>("coordinates");
                                    SireMaths::Vector tmpc = (molecule.atom(first_alchemical_bond->atom1())).property<SireMaths::Vector>("coordinates");
                                    bond1_vec = tmp1 - tmpc;
                                    bond2_vec = tmp2 - tmpc;
                                    if (Debug)
                                        qDebug() << "Bond1 Atom1 = Bond2 Atom0";
                                }
                                else if (first_alchemical_bond->atom1() == second_alchemical_bond->atom1())
                                {
                                    SireMaths::Vector tmp1 = (molecule.atom(first_alchemical_bond->atom0())).property<SireMaths::Vector>("coordinates");
                                    SireMaths::Vector tmp2 = (molecule.atom(second_alchemical_bond->atom0())).property<SireMaths::Vector>("coordinates");
                                    SireMaths::Vector tmpc = (molecule.atom(first_alchemical_bond->atom1())).property<SireMaths::Vector>("coordinates");
                                    bond1_vec = tmp1 - tmpc;
                                    bond2_vec = tmp2 - tmpc;
                                    if (Debug)
                                        qDebug() << "Bond1 Atom1 = Bond2 Atom1";
                                }
                                else
                                    throw SireError::program_bug(QObject::tr("No coorner bond"), CODELOC);

                                if (Debug)
                                {
                                    if (first_alchemical_distance != -1.0)
                                    {
                                        qDebug() << "First vector X = " << (bond1_vec.normalise() * first_alchemical_distance).x();
                                        qDebug() << "First vector Y = " << (bond1_vec.normalise() * first_alchemical_distance).y();
                                        qDebug() << "First vector Z = " << (bond1_vec.normalise() * first_alchemical_distance).z();
                                    }
                                    else
                                    {
                                        qDebug() << "First vector X = " << (bond1_vec).x();
                                        qDebug() << "First vector Y = " << (bond1_vec).y();
                                        qDebug() << "First vector Z = " << (bond1_vec).z();
                                    }

                                    if (second_alchemical_distance != -1.0)
                                    {
                                        qDebug() << "Second vector X = " << (bond2_vec.normalise() * second_alchemical_distance).x();
                                        qDebug() << "Second vector Y = " << (bond2_vec.normalise() * second_alchemical_distance).y();
                                        qDebug() << "Second vector Z = " << (bond2_vec.normalise() * second_alchemical_distance).z();
                                    }
                                    else
                                    {
                                        qDebug() << "Second vector X = " << (bond2_vec).x();
                                        qDebug() << "Second vector Y = " << (bond2_vec).y();
                                        qDebug() << "Second vector Z = " << (bond2_vec).z();
                                    }
                                }

                                double constraint_distance;

                                double eq_angle = solute_angle_perturbation_params[3] * Alchemical_value + (1.0 - Alchemical_value) * solute_angle_perturbation_params[2];

                                if (first_alchemical_distance == -1.0 && second_alchemical_distance != -1.0)
                                {

                                    //Carnot's theorem a^2 = c^2 + b^2 - a*b*c*cos(bc)
                                    double sq = bond1_vec.length() * bond1_vec.length() +
                                        (bond2_vec.normalise() * second_alchemical_distance).length() * (bond2_vec.normalise() * second_alchemical_distance).length();

                                    double dp = 2.0 * bond1_vec.length() * (bond2_vec.normalise() * second_alchemical_distance).length() * cos(eq_angle);

                                    constraint_distance = sqrt(sq - dp);

                                    //constraint_distance = (bond1_vec - bond2_vec.normalise() * second_alchemical_distance).length();
                                }
                                else if (first_alchemical_distance != -1.0 && second_alchemical_distance == -1.0)
                                {

                                    //Carnot theorem a^2 = c^2 + b^2 - a*b*c*cos(bc)
                                    double sq = bond2_vec.length() * bond2_vec.length() +
                                        (bond1_vec.normalise() * first_alchemical_distance).length() * (bond1_vec.normalise() * first_alchemical_distance).length();

                                    double dp = 2.0 * bond2_vec.length() * (bond2_vec.normalise() * first_alchemical_distance).length() * cos(eq_angle);

                                    constraint_distance = sqrt(sq - dp);

                                    //constraint_distance = (bond1_vec.normalise() * first_alchemical_distance - bond2_vec).length();
                                }
                                else if (first_alchemical_distance != -1.0 && second_alchemical_distance != -1.0)
                                {

                                    //Carnot's theorem a^2 = c^2 + b^2 - a*b*c*cos(bc)
                                    double sq = (bond1_vec.normalise() * first_alchemical_distance).length() * (bond1_vec.normalise() * first_alchemical_distance).length() +
                                        (bond2_vec.normalise() * second_alchemical_distance).length() * (bond2_vec.normalise() * second_alchemical_distance).length();

                                    double dp = 2.0 * (bond1_vec.normalise() * first_alchemical_distance).length() * (bond2_vec.normalise() * second_alchemical_distance).length() * cos(eq_angle);

                                    constraint_distance = sqrt(sq - dp);

                                    //constraint_distance = (bond1_vec.normalise() * first_alchemical_distance - bond2_vec.normalise() * second_alchemical_distance).length();
                                }
                                else
                                    throw SireError::program_bug(QObject::tr("The angle does not contain perturbed bond"), CODELOC);


                                system_openmm->addConstraint(idx0, idx2, constraint_distance * OpenMM::NmPerAngstrom);

                                if (Debug)
                                    qDebug() << "CONSTRAINT DISTANCE = " << constraint_distance << " A";
                            }

                        }//end if HANGLES
                        else
                        {
                            solute_angle_perturbation->addAngle(idx0, idx1, idx2, solute_angle_perturbation_params);
                            if (Debug)
                                qDebug() << "Added perturbed angle";
                        }


                        angle_pert_list.append(AngleID(three.atom0(), three.atom1(), three.atom2()));
                        angle_pert_swap_list.append(AngleID(three.atom2(), three.atom1(), three.atom0()));

                        if (Debug)
                        {
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
                    if (str == "SireMM::FourAtomPerturbation")
                    {

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

                        if (Debug)
                        {
                            qDebug() << "IDX0 = " << idx0 << "IDX1 = " << idx1 << "IDX2 = " << idx2 << "IDX3 = " << idx3;
                            qDebug() << "Dihedral String = " << openmm_str.c_str();
                            qDebug() << "Dihedral Normal String = " << four.perturbExpression().toString() << "\n";
                        }

                        OpenMM::CustomTorsionForce* solute_torsion_perturbation = NULL;

                        solute_torsion_perturbation = new OpenMM::CustomTorsionForce(openmm_str);
                        solute_torsion_perturbation->addPerTorsionParameter("KJPerKcal");
                        solute_torsion_perturbation_params[0] = 4.184;
                        solute_torsion_perturbation->addGlobalParameter("lamdih", Alchemical_value);
                        solute_torsion_perturbation->addTorsion(idx0, idx1, idx2, idx3, solute_torsion_perturbation_params);

                        //********************************BONDED ENERGY TORSIONS ARE ADDED TO THE SYSTEM*****************************
                        solute_torsion_perturbation->setForceGroup(0);
                        system_openmm->addForce(solute_torsion_perturbation);

                        perturbed_energies_tmp[7] = true; //Torsions are added to the system

                        dihedral_pert_list.append(DihedralID(four.atom0(), four.atom1(), four.atom2(), four.atom3()));
                        dihedral_pert_swap_list.append(DihedralID(four.atom3(), four.atom1(), four.atom2(), four.atom0()));

                        improper_pert_list.append(ImproperID(four.atom0(), four.atom1(), four.atom2(), four.atom3()));
                        improper_pert_swap_list.append(ImproperID(four.atom0(), four.atom1(), four.atom3(), four.atom2()));

                        if (Debug)
                        {
                            qDebug() << "Atom0 = " << four.atom0().asA<AtomIdx>().value() <<
                                "Atom1 = " << four.atom1().asA<AtomIdx>().value() <<
                                "Atom2 = " << four.atom2().asA<AtomIdx>().value() <<
                                "Atom3 = " << four.atom3().asA<AtomIdx>().value() << "\n";
                        }
                    }

                }
            }//end for perturbations

        }//end solute molecule perturbation


        // The bonded parameters are stored in "amberparameters"
        AmberParameters amber_params = molecule.property("amberparameters").asA<AmberParameters>();
        QList<BondID> bonds_ff = amber_params.getAllBonds();
        QVector<BondID> bonds = bonds_ff.toVector();
        ResName molfirstresname = molecule.residues()(0).name();
        //BOND

        for (int j = 0; j < bonds_ff.length(); j++)
        {

            BondID bond_ff = bonds_ff[j];
            QList<double> bond_params = amber_params.getParams(bond_ff);
            double k = bond_params[0];
            double r0 = bond_params[1];

            int idx0 = bonds[j].atom0().asA<AtomIdx>().value();
            int idx1 = bonds[j].atom1().asA<AtomIdx>().value();

            if (solute.contains(molecule))
            {

                if (bond_pert_list.indexOf(bond_ff) != -1 || bond_pert_swap_list.indexOf(bond_ff) != -1)
                {//Solute molecule. Check if the current solute bond is in the perturbed bond list
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

            if (flag_constraint == NONE)
            {
                //JM 10/16 If constraint water flag is on and if molecule is a water molecule then apply constraint
                if (flag_constraint_water and molfirstresname == ResName("WAT"))
                    system_openmm->addConstraint(idx0, idx1, r0 * OpenMM::NmPerAngstrom);
                else
                    bondStretch_openmm->addBond(idx0, idx1, r0 * OpenMM::NmPerAngstrom, k * 2.0 * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);

              //cout << "\nBOND ADDED TO "<< atom0.toStdString() << " AND " << atom1.toStdString() << "\n";
            }
            else if (flag_constraint == ALLBONDS || flag_constraint == HANGLES)
            {
                system_openmm->addConstraint(idx0, idx1, r0 * OpenMM::NmPerAngstrom);
                //cout << "\nALLBONDS or HANGLES ADDED BOND CONSTRAINT TO " << atom0.toStdString() << " AND " << atom1.toStdString() << "\n";
            }
            else if (flag_constraint == HBONDS)
            {

                if ((atom0[6] == 'H') || (atom1[6] == 'H'))
                {
                    system_openmm->addConstraint(idx0, idx1, r0 * OpenMM::NmPerAngstrom);
                    //cout << "\nHBONDS ADDED BOND CONSTRAINT TO " << atom0.toStdString() << " AND " << atom1.toStdString() << "\n";
                }
                else
                {
                    bondStretch_openmm->addBond(idx0, idx1, r0 * OpenMM::NmPerAngstrom, k * 2.0 * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
                    //cout << "\nHBONDS ADDED BOND TO " << atom0.toStdString() << " AND " << atom1.toStdString() << "\n";
                }
            }

            //Bond exclusion List
            bondPairs.push_back(std::make_pair(idx0, idx1));
        }


        //Angles

        QList<AngleID> angles_ff = amber_params.getAllAngles();
        QVector<AngleID> angles = angles_ff.toVector();

        for (int j = 0; j < angles_ff.length(); j++)
        {

            AngleID angle_ff = angles_ff[j];
            QList<double> angle_params = amber_params.getParams(angle_ff);


            double k = angle_params[0];
            double theta0 = angle_params[1]; // It is already in radiant

            int idx0 = angles[j].atom0().asA<AtomIdx>().value();
            int idx1 = angles[j].atom1().asA<AtomIdx>().value();
            int idx2 = angles[j].atom2().asA<AtomIdx>().value();

            if (solute.contains(molecule))
            {
                if (angle_pert_list.indexOf(angle_ff) != -1 || angle_pert_swap_list.indexOf(angle_ff) != -1)
                {//Solute molecule. Check if the current solute angle is in the perturbed angle list
                    if (Debug)
                        qDebug() << "Found Perturbed Angle\n";
                    continue;
                }
                else
                {

                    if (Debug)
                        qDebug() << "Solute normal Angle - Atom0 = " << idx0 << "Atom1 = " << idx1 << "Atom2 = " << idx2 << "theta0 = " << theta0 << " k = " << k << "\n";

                    idx0 = idx0 + num_atoms_till_i;
                    idx1 = idx1 + num_atoms_till_i;
                    idx2 = idx2 + num_atoms_till_i;
                    bondBend_openmm->addAngle(idx0, idx1, idx2, theta0, k * 2.0 * OpenMM::KJPerKcal);
                    continue;
                }
            }
            if (Debug)
                qDebug() << "Angle - Atom0 = " << idx0 << "Atom1 = " << idx1 << "Atom2 = " << idx2 << "\n";

            QString atom0 = molecule.atom(AtomIdx(idx0)).toString();
            QString atom1 = molecule.atom(AtomIdx(idx1)).toString();
            QString atom2 = molecule.atom(AtomIdx(idx2)).toString();


            Vector diff = c[idx2] - c[idx0];

            idx0 = idx0 + num_atoms_till_i;
            idx1 = idx1 + num_atoms_till_i;
            idx2 = idx2 + num_atoms_till_i;

            if (flag_constraint == HANGLES)
            {

                if (((atom0[6] == 'H') && (atom2[6] == 'H')))
                {
                    system_openmm->addConstraint(idx0, idx2, diff.length() * OpenMM::NmPerAngstrom);
                }
                else if (((atom0[6] == 'H') && (atom1[6] == 'O')) || ((atom1[6] == 'O') && (atom2[6] == 'H')))
                {
                    system_openmm->addConstraint(idx0, idx2, diff.length() * OpenMM::NmPerAngstrom);
                }
                else
                {
                    bondBend_openmm->addAngle(idx0, idx1, idx2, theta0, k * 2.0 * OpenMM::KJPerKcal);
                }
            }
            else
            {
                bondBend_openmm->addAngle(idx0, idx1, idx2, theta0, k * 2.0 * OpenMM::KJPerKcal);
            }
        }//end of angles

        //Dihedrals

        QList<DihedralID> dihedrals_ff = amber_params.getAllDihedrals();
        QVector<DihedralID> dihedrals = dihedrals_ff.toVector();

        for (int j = 0; j < dihedrals_ff.length(); j++)
        {

            DihedralID dihedral_ff = dihedrals_ff[j];
            QList<double> dihedral_params = amber_params.getParams(dihedral_ff);

            int idx0 = dihedrals[j].atom0().asA<AtomIdx>().value() + num_atoms_till_i;
            int idx1 = dihedrals[j].atom1().asA<AtomIdx>().value() + num_atoms_till_i;
            int idx2 = dihedrals[j].atom2().asA<AtomIdx>().value() + num_atoms_till_i;
            int idx3 = dihedrals[j].atom3().asA<AtomIdx>().value() + num_atoms_till_i;


            if (Debug)
            {
                qDebug() << "TOTAL Dihedral between atom global index " << idx0 - num_atoms_till_i <<
                    " and " << idx1 - num_atoms_till_i <<
                    " and " << idx2 - num_atoms_till_i <<
                    " and " << idx3 - num_atoms_till_i << "\n";
            }

            if (solute.contains(molecule))
            {
                if (dihedral_pert_list.indexOf(dihedral_ff) != -1 || dihedral_pert_swap_list.indexOf(dihedral_ff) != -1)
                {//Solute molecule. Check if the current solute dihedral is in the perturbed dihedral list
                    if (Debug)
                        qDebug() << "Found Perturbed Dihedral\n";
                    continue;
                }
            }

            // Variable number of parameters
            for (int k = 0; k < dihedral_params.length(); k = k + 3)
            {
                double v = dihedral_params[ k ];
                int periodicity = dihedral_params[ k + 1 ];
                double phase = dihedral_params[ k + 2 ];
                bondTorsion_openmm->addTorsion(idx0, idx1, idx2, idx3, periodicity, phase, v * OpenMM::KJPerKcal);
                if (Debug)
                {
                    qDebug() << "Dihedral between atom global index " << idx0 - num_atoms_till_i <<
                        " and " << idx1 - num_atoms_till_i <<
                        " and " << idx2 - num_atoms_till_i <<
                        " and " << idx3 - num_atoms_till_i << "\n";
                    qDebug() << "Amplitude_dih = " << v << " periodicity " << periodicity << " phase " << phase << "\n";
                }
            }
        } // end of dihedrals


        //Improper Dihedrals

        QList<ImproperID> impropers_ff = amber_params.getAllImpropers();
        QVector<ImproperID> impropers = impropers_ff.toVector();

        for (int j = 0; j < impropers_ff.length(); j++)
        {
            ImproperID improper_ff = impropers_ff[j];
            QList<double> improper_params = amber_params.getParams(improper_ff);

            int idx0 = impropers[j].atom0().asA<AtomIdx>().value() + num_atoms_till_i;
            int idx1 = impropers[j].atom1().asA<AtomIdx>().value() + num_atoms_till_i;
            int idx2 = impropers[j].atom2().asA<AtomIdx>().value() + num_atoms_till_i;
            int idx3 = impropers[j].atom3().asA<AtomIdx>().value() + num_atoms_till_i;

            if (Debug)
            {
                qDebug() << "TOTAL Improper between atom global index " << idx0 - num_atoms_till_i <<
                    " and " << idx1 - num_atoms_till_i <<
                    " and " << idx2 - num_atoms_till_i <<
                    " and " << idx3 - num_atoms_till_i << "\n";
            }

            if (solute.contains(molecule))
            {//Solute molecule. Check if the current solute dihedral is in the perturbed improper list
                if (improper_pert_list.indexOf(improper_ff) != -1 || improper_pert_swap_list.indexOf(improper_ff) != -1)
                {
                    if (Debug)
                        qDebug() << "Found Perturbed Improper\n";
                    continue;
                }
            }

            // Variable number of parameters
            for (int k = 0; k < improper_params.length(); k = k + 3)
            {
                double v = improper_params[ k ];
                int periodicity = improper_params[ k + 1 ];
                double phase = improper_params[ k + 2 ];

                bondTorsion_openmm->addTorsion(idx0, idx1, idx2, idx3, periodicity, phase, v * OpenMM::KJPerKcal);
                if (Debug)
                {
                    qDebug() << "Improper between atom global index " << idx0 - num_atoms_till_i <<
                        " and " << idx1 - num_atoms_till_i <<
                        " and " << idx2 - num_atoms_till_i <<
                        " and " << idx3 - num_atoms_till_i << "\n";
                    qDebug() << "Amplitude_imp = " << v << " periodicity " << periodicity << " phase " << phase << "\n";
                }
            }
        }//end of impropers


        // Variable 1,4 scaling factors
        QList<BondID> pairs14_ff = amber_params.getAll14Pairs();
        QVector<BondID> pairs14 = pairs14_ff.toVector();

        for (int j = 0; j < pairs14_ff.length(); j++)
        {

            BondID pair14_ff = pairs14_ff[j];

            QList<double> pair14_params = amber_params.get14PairParams(pair14_ff);

            double cscl = pair14_params[0];
            double ljscl = pair14_params[1];

            if (Debug)
                qDebug() << " cscl@ " << cscl << " ljscl " << ljscl;

            // Add to custom pairs if scale factor differs from default
            if (abs(cscl - Coulomb14Scale) > 0.0001 or abs(ljscl - LennardJones14Scale) > 0.0001)
            {

                int idx0 = pair14_ff.atom0().asA<AtomIdx>().value() + num_atoms_till_i;
                int idx1 = pair14_ff.atom1().asA<AtomIdx>().value() + num_atoms_till_i;

                QPair<int, int> indices_pair(idx0, idx1);
                QPair<double, double> scl_pair(cscl, ljscl);
                custom14pairs.insert(indices_pair, scl_pair);

                special_14 = true;

                if (Debug)
                    qDebug() << "IDX0 = " << idx0 << " IDX1 =" << idx1 << "14 OpenMM Index";
            }
        }// end of variable 1,4 scaling factors

        num_atoms_till_i = num_atoms_till_i + num_atoms_molecule;

    }// end of loop over molecules


    if (Debug)
    {
        if (nions != 0)
            qDebug() << "\n\nNumber of ions = " << nions << "\n\n";
    }

    //Exclude the 1-2, 1-3 bonded atoms from nonbonded forces, and scale down 1-4 bonded atoms

    nonbond_openmm->createExceptionsFromBonds(bondPairs, Coulomb14Scale, LennardJones14Scale);

    if (CMMremoval_frequency > 0)
    {
        OpenMM::CMMotionRemover * cmmotionremover = new OpenMM::CMMotionRemover(CMMremoval_frequency);

        system_openmm->addForce(cmmotionremover);

        if (Debug)
            qDebug() << "\n\nWill remove Center of Mass motion every " << CMMremoval_frequency << " steps\n\n";
    }


    int num_exceptions = nonbond_openmm->getNumExceptions();

    if (Debug)
        qDebug() << "NUM EXCEPTIONS = " << num_exceptions << "\n";

    for (int i = 0; i < num_exceptions; i++)
    {

        int p1, p2;

        double charge_prod, sigma_avg, epsilon_avg;


        nonbond_openmm->getExceptionParameters(i, p1, p2, charge_prod, sigma_avg, epsilon_avg);

        if (Debug)
            qDebug() << "Exception = " << i << " p1 = " << p1 << " p2 = " << p2 << " charge prod = " << charge_prod << " sigma avg = " << sigma_avg << " epsilon_avg = " << epsilon_avg << "\n";

        if (!(charge_prod == 0 && sigma_avg == 1 && epsilon_avg == 0))
        {//1-4 interactions

            QVector<double> perturbed_14_tmp(13);

            std::vector<double> p1_params(10);
            std::vector<double> p2_params(10);

            custom_force_field->getParticleParameters(p1, p1_params);
            custom_force_field->getParticleParameters(p2, p2_params);

            double Qstart_p1 = p1_params[0];
            double Qend_p1 = p1_params[1];
            double Epstart_p1 = p1_params[2];
            double Epend_p1 = p1_params[3];
            double Sigstart_p1 = p1_params[4];
            double Sigend_p1 = p1_params[5];
            double isHard_p1 = p1_params[6];
            double isTodummy_p1 = p1_params[7];
            double isFromdummy_p1 = p1_params[8];

            double Qstart_p2 = p2_params[0];
            double Qend_p2 = p2_params[1];
            double Epstart_p2 = p2_params[2];
            double Epend_p2 = p2_params[3];
            double Sigstart_p2 = p2_params[4];
            double Sigend_p2 = p2_params[5];
            double isHard_p2 = p2_params[6];
            double isTodummy_p2 = p2_params[7];
            double isFromdummy_p2 = p2_params[8];

            double charge_prod_start, charge_prod_end, charge_prod_mix;
            double sigma_avg_start, sigma_avg_end;
            double epsilon_avg_start, epsilon_avg_end, epsilon_avg_mix;

            double Coulomb14Scale_tmp = Coulomb14Scale;
            double LennardJones14Scale_tmp = LennardJones14Scale;

            if (special_14)
            {

                QPair<double, double> sc_factors;

                QPair<int, int> indices_pair(p1, p2);
                QHash< QPair<int, int>, QPair<double, double> >::const_iterator i_pair = custom14pairs.find(indices_pair);


                if (i_pair != custom14pairs.end())
                {

                    sc_factors = i_pair.value();
                    Coulomb14Scale_tmp = sc_factors.first;
                    LennardJones14Scale_tmp = sc_factors.second;

                    if (Debug)
                        qDebug() << "The pair ( " << p1 << ", " << p2 << " ) is 14 special no swap pair";
                }
                else
                {

                    QPair<int, int> indices_swap_pair(p2, p1);
                    QHash< QPair<int, int>, QPair<double, double> >::const_iterator i_swap_pair = custom14pairs.find(indices_swap_pair);

                    if (i_swap_pair != custom14pairs.end())
                    {

                        sc_factors = i_swap_pair.value();
                        Coulomb14Scale_tmp = sc_factors.first;
                        LennardJones14Scale_tmp = sc_factors.second;

                        if (Debug)
                            qDebug() << "The pair ( " << p2 << ", " << p1 << " ) is 14 special swap pair";

                    }
                }
            }

            charge_prod_start = Qstart_p1 * Qstart_p2 * Coulomb14Scale_tmp;
            charge_prod_end = Qend_p1 * Qend_p2 * Coulomb14Scale_tmp;
            charge_prod_mix = (Qend_p1 * Qstart_p2 + Qstart_p1 * Qend_p2) * Coulomb14Scale_tmp;

            sigma_avg_start = (Sigstart_p1 + Sigstart_p2) / 2.0;
            sigma_avg_end = (Sigend_p1 + Sigend_p2) / 2.0;

            epsilon_avg_start = Epstart_p1 * Epstart_p2 * LennardJones14Scale_tmp * LennardJones14Scale_tmp;
            epsilon_avg_end = Epend_p1 * Epend_p2 * LennardJones14Scale_tmp * LennardJones14Scale_tmp;
            epsilon_avg_mix = (Epend_p1 * Epstart_p2 + Epstart_p1 * Epend_p2) * LennardJones14Scale_tmp * LennardJones14Scale_tmp;


            std::vector<double> params(8);

            params[0] = charge_prod_start;
            params[1] = charge_prod_end;
            params[2] = charge_prod_mix;
            params[3] = epsilon_avg_start;
            params[4] = epsilon_avg_end;
            params[5] = epsilon_avg_mix;
            params[6] = sigma_avg_start;
            params[7] = sigma_avg_end;


            if (Debug)
            {

                qDebug() << "Particle p1 = " << p1 << "\nQstart = " << Qstart_p1 << "\nQend = " << Qend_p1
                    << "\nEpstart = " << Epstart_p1 << "\nEpend = " << Epend_p1
                    << "\nSgstart = " << Sigstart_p1 << "\nSgend = " << Sigend_p1
                    << "\nisHard = " << isHard_p1 << "\nisTodummy = " << isTodummy_p1 << "\nisFromdummy = " << isFromdummy_p1 << "\n";
                qDebug() << "Particle p2 = " << p2 << "\nQstart = " << Qstart_p2 << "\nQend = " << Qend_p2
                    << "\nEpstart = " << Epstart_p2 << "\nEpend = " << Epend_p2
                    << "\nSgstart = " << Sigstart_p2 << "\nSgend = " << Sigend_p2
                    << "\nisHard = " << isHard_p2 << "\nisTodummy = " << isTodummy_p2 << "\nisFromdummy = " << isFromdummy_p2 << "\n";



                qDebug() << "Product Charge start = " << charge_prod_start << "\nProduct Charge end = " << charge_prod_end << "\nProduct Chrage mixed = " << charge_prod_mix
                    << "\nEpsilon average start = " << epsilon_avg_start << "\nEpsilon average end = " << epsilon_avg_end << "\nEpsilon average mixed = " << charge_prod_mix
                    << "\nSigma average start = " << sigma_avg_start << "\nSigma average end = " << sigma_avg_end;
                qDebug() << "Columbic Scale Factor = " << Coulomb14Scale_tmp << " Lennard-Jones Scale Factor = " << LennardJones14Scale_tmp << "\n";
            }

            if ((isHard_p1 == 1.0 && isHard_p2 == 1.0))
            {

                custom_intra_14_clj->addBond(p1, p2, params);


                if (Debug)
                    qDebug() << "Added clj Hard 1-4\n";
            }
            else if ((isTodummy_p1 == 1.0 && isTodummy_p2 == 1.0) || (isHard_p1 == 1.0 && isTodummy_p2 == 1.0) || (isHard_p2 == 1.0 && isTodummy_p1 == 1.0))
            {

                custom_intra_14_todummy->addBond(p1, p2, params);

                if (Debug)
                    qDebug() << "Added soft TO dummy 1-4\n";
            }

            else if ((isFromdummy_p1 == 1.0 && isFromdummy_p2 == 1.0) || (isHard_p1 == 1.0 && isFromdummy_p2 == 1.0) || (isHard_p2 == 1.0 && isFromdummy_p1 == 1.0))
            {

                custom_intra_14_fromdummy->addBond(p1, p2, params);

                if (Debug)
                    qDebug() << "Added soft FROM dummy 1-4\n";

            }

            else if ((isFromdummy_p1 == 1.0 && isTodummy_p2 == 1.0) || (isFromdummy_p2 == 1.0 && isTodummy_p1 == 1.0))
            {

                custom_intra_14_fromdummy_todummy->addBond(p1, p2, params);

                if (Debug)
                    qDebug() << "Added soft FROM dummy TO dummy 1-4\n";

            }


        }//end if 1-4 interactions

        custom_force_field->addExclusion(p1, p2);

    }

    /*************************************NON BONDED INTERACTIONS*********************************************/

    int npairs = (custom_force_field->getNumParticles() * (custom_force_field->getNumParticles() - 1)) / 2;

    if (Debug)
    {
        qDebug() << "Num pairs  = " << npairs;
        qDebug() << "Num bonds 1-4 Hard = " << custom_intra_14_clj->getNumBonds();
        qDebug() << "Num bonds 1-4 To Dummy = " << custom_intra_14_todummy->getNumBonds();
        qDebug() << "Num bonds 1-4 From Dummy = " << custom_intra_14_fromdummy->getNumBonds();
        qDebug() << "Num bonds 1-4 From Dummy To Dummy = " << custom_intra_14_fromdummy_todummy->getNumBonds();
    }

    if (npairs != num_exceptions)
    {
        custom_force_field->setForceGroup(0);
        system_openmm->addForce(custom_force_field);
        perturbed_energies_tmp[0] = true; //Custom non bonded 1-5 is added to the system
        if (Debug)
            qDebug() << "Added 1-5";
    }

    if (custom_intra_14_clj->getNumBonds() != 0)
    {
        custom_intra_14_clj->setForceGroup(0);
        system_openmm->addForce(custom_intra_14_clj);
        perturbed_energies_tmp[1] = true; //Custom non bonded 1-4 is added to the system
        if (Debug)
            qDebug() << "Added 1-4 CLJ";
    }

    if (custom_intra_14_todummy->getNumBonds() != 0)
    {
        custom_intra_14_todummy->setForceGroup(0);
        system_openmm->addForce(custom_intra_14_todummy);
        perturbed_energies_tmp[2] = true; //Custom non bonded 1-4 is added to the system
        if (Debug)
            qDebug() << "Added 1-4 To Dummy";
    }


    if (custom_intra_14_fromdummy->getNumBonds() != 0)
    {
        custom_intra_14_fromdummy->setForceGroup(0);
        system_openmm->addForce(custom_intra_14_fromdummy);
        perturbed_energies_tmp[3] = true; //Custom non bonded 1-4 is added to the system
        if (Debug)
            qDebug() << "Added 1-4 From Dummy";

    }
    if (custom_intra_14_fromdummy_todummy->getNumBonds() != 0)
    {
        custom_intra_14_fromdummy_todummy->setForceGroup(0);
        system_openmm->addForce(custom_intra_14_fromdummy_todummy);
        perturbed_energies_tmp[4] = true; //Custom non bonded 1-4 is added to the system
        if (Debug)
            qDebug() << "Added 1-4 From Dummy To Dummy";
    }

    /*****************************************BONDED INTERACTIONS***********************************************/

    if (bondStretch_openmm->getNumBonds() != 0)
    {
        bondStretch_openmm->setForceGroup(1);
        system_openmm->addForce(bondStretch_openmm);
        if (Debug)
            qDebug() << "Added Internal Bond energy term";
    }

    if (bondBend_openmm->getNumAngles() != 0)
    {
        bondBend_openmm->setForceGroup(1);
        system_openmm->addForce(bondBend_openmm);
        if (Debug)
            qDebug() << "Added Internal Angle energy term";
    }

    if (bondTorsion_openmm->getNumTorsions() != 0)
    {
        bondTorsion_openmm->setForceGroup(1);
        system_openmm->addForce(bondTorsion_openmm);
        if (Debug)
            qDebug() << "Added Internal Torsion energy term";
    }

    if (solute_bond_perturbation->getNumBonds() != 0)
    {
        solute_bond_perturbation->setForceGroup(0);
        system_openmm->addForce(solute_bond_perturbation);
        perturbed_energies_tmp[5] = true; //Custom bonded is added to the system
        if (Debug)
            qDebug() << "Added Perturbed Internal Bond energy term";
    }

    if (solute_angle_perturbation->getNumAngles() != 0)
    {
        solute_angle_perturbation->setForceGroup(0);
        system_openmm->addForce(solute_angle_perturbation);
        perturbed_energies_tmp[6] = true; //Custom bonded is added to the system
        if (Debug)
            qDebug() << "Added Perturbed Internal Angle energy term";
    }


    perturbed_energies = perturbed_energies_tmp;

    //IMPORTANT: PERTURBED ENERGY TORSIONS ARE ADDED ABOVE
    bool UseLink_flag = true;

    //Distance Restaint. All the information are stored in the first molecule only.

    if (UseLink_flag == true)
    {

        Molecule molecule = moleculegroup.moleculeAt(0).molecule();

        bool haslinkinfo = molecule.hasProperty("linkbonds");

        if (haslinkinfo)
        {

            std::vector<double> custom_bond_link_par(3);

            Properties linkprop = molecule.property("linkbonds").asA<Properties>();

            int nlinks = linkprop.property(QString("nbondlinks")).asA<VariantProperty>().toInt();

            if (Debug)
                qDebug() << "Number of constraint links = " << nlinks;

            for (int i = 0; i < nlinks; i++)
            {

                int atomnum0 = linkprop.property(QString("AtomNum0(%1)").arg(i)).asA<VariantProperty>().toInt();
                int atomnum1 = linkprop.property(QString("AtomNum1(%1)").arg(i)).asA<VariantProperty>().toInt();
                double reql = linkprop.property(QString("reql(%1)").arg(i)).asA<VariantProperty>().toDouble();
                double kl = linkprop.property(QString("kl(%1)").arg(i)).asA<VariantProperty>().toDouble();
                double dl = linkprop.property(QString("dl(%1)").arg(i)).asA<VariantProperty>().toDouble();

                int openmmindex0 = AtomNumToOpenMMIndex[atomnum0];
                int openmmindex1 = AtomNumToOpenMMIndex[atomnum1];

                custom_bond_link_par[0] = reql * OpenMM::NmPerAngstrom; //req
                custom_bond_link_par[1] = kl * (OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm); //k
                custom_bond_link_par[2] = dl * OpenMM::NmPerAngstrom; //dl

                if (Debug)
                {
                    qDebug() << "atomnum0 = " << atomnum0 << " openmmindex0 =" << openmmindex0;
                    qDebug() << "atomnum1 = " << atomnum1 << " openmmindex1 =" << openmmindex1;
                    qDebug() << "Req = " << reql << " kl = " << kl << " dl = " << dl;
                }

                custom_link_bond->addBond(openmmindex0, openmmindex1, custom_bond_link_par);

            }

            system_openmm->addForce(custom_link_bond);
        }

    }//end of bond link flag



    this->openmm_system = system_openmm;
    this->isSystemInitialised = true;


}

/**
 *
 * @param workspace
 * @param timestep
 */

void OpenMMFrEnergyST::createContext(IntegratorWorkspace &workspace, SireUnits::Dimension::Time timestep)
{
    bool Debug = false;

    if (Debug)
    {
        qDebug() << "In OpenMMFrEnergyST::createContext()\n\n";
        qDebug() << isContextInitialised;
        qDebug() << reinitialise_context;
    }

    // Check that the openmm system has been initialised
    // !! Should check that the workspace is compatible with molgroup
    if (not this->isSystemInitialised)
    {

        qDebug() << "Not initialised ! ";
        throw SireError::program_bug(QObject::tr(
                                                 "OpenMMFrEnergyST should have been initialised before calling integrate."), CODELOC);
    }

    OpenMM::System *system_openmm = openmm_system;

    int nats = system_openmm->getNumParticles();

    if (Debug)
        qDebug() << " openmm nats " << nats;


    // Integrator

    const double dt = convertTo(timestep.value(), picosecond);
    const double converted_Temperature = convertTo(Temperature.value(), kelvin);
    const double converted_friction = convertTo(friction.value(), picosecond);

    if (!isContextInitialised || (isContextInitialised && reinitialise_context))
    {
        OpenMM::Integrator * integrator_openmm = NULL;

        if (Integrator_type == "leapfrogverlet")
            integrator_openmm = new OpenMM::VerletIntegrator(dt); //dt in picosecond
        else if (Integrator_type == "variableleapfrogverlet")
            integrator_openmm = new OpenMM::VariableVerletIntegrator(integration_tol); //integration tolerance error unitless
        else if (Integrator_type == "langevin")
            integrator_openmm = new OpenMM::LangevinIntegrator(converted_Temperature, converted_friction, dt);
        else if (Integrator_type == "variablelangevin")
            integrator_openmm = new OpenMM::VariableLangevinIntegrator(converted_Temperature, converted_friction, integration_tol);
        else if (Integrator_type == "brownian")
            integrator_openmm = new OpenMM::BrownianIntegrator(converted_Temperature, converted_friction, dt);
        else
            throw SireError::program_bug(QObject::tr("The user defined Integrator type is not supported. Available types are leapfrogverlet, variableleapfrogverlet, langevin, variablelangevin, brownian"), CODELOC);

        if (Debug)
        {
            qDebug() << "Using Integrator: " << Integrator_type;

            qDebug() << "Integration step = " << dt << " ps";

            if (Integrator_type == "variablelangevin" || Integrator_type == "variableleapfrogverlet")
            {
                qDebug() << "Integration Tol = " << integration_tol;
            }
            if (Integrator_type == "langevin" || Integrator_type == "variablelangevin" || Integrator_type == "brownian")
            {
                qDebug() << "Converted Friction = " << converted_friction << "1/ps";
            }
        }

        OpenMM::Platform& platform_openmm = OpenMM::Platform::getPlatformByName(platform_type.toStdString());

        if (platform_type == "OpenCL")
        {

            const std::string prop = std::string("OpenCLDeviceIndex");
            const std::string prec = std::string("OpenCLPrecision");

            platform_openmm.setPropertyDefaultValue(prop, device_index.toStdString());
            platform_openmm.setPropertyDefaultValue(prec, precision.toStdString());

            if (Debug)
            {
                qDebug() << "Setting up OpenCL default Index to " << device_index;
                qDebug() << "Setting up OpenCL precision to" << precision;
            }
        }
        else if (platform_type == "CUDA")
        {

            const std::string prop = std::string("CudaDeviceIndex");
            const std::string prec = std::string("CudaPrecision");

            platform_openmm.setPropertyDefaultValue(prop, device_index.toStdString());
            platform_openmm.setPropertyDefaultValue(prec, precision.toStdString());

            if (Debug)
            {
                qDebug() << "Setting up CUDA default Index to " << device_index;
                qDebug() << "Setting up CUDA precision to" << precision;
            }

        }



        delete openmm_context;
        openmm_context = new OpenMM::Context(*system_openmm, *integrator_openmm, platform_openmm);
        this->isContextInitialised = true;

    }


    if (Debug)
        qDebug() << "\n Using OpenMM platform = " << openmm_context->getPlatform().getName().c_str() << "\n";

    // Now update coordinates / velocities / dimensions with sire data
    AtomicVelocityWorkspace &ws = workspace.asA<AtomicVelocityWorkspace>();

    if (CutoffType == "cutoffperiodic")
    {

        const System & ptr_sys = ws.system();
        const PropertyName &space_property = PropertyName("space");
        const PeriodicBox &space = ptr_sys.property(space_property).asA<PeriodicBox>();

        const double Box_x_Edge_Length = space.dimensions()[0] * OpenMM::NmPerAngstrom; //units in nm
        const double Box_y_Edge_Length = space.dimensions()[1] * OpenMM::NmPerAngstrom; //units in nm
        const double Box_z_Edge_Length = space.dimensions()[2] * OpenMM::NmPerAngstrom; //units in nm

        if (Debug)
            qDebug() << "\nBOX SIZE [A] = (" << space.dimensions()[0] << " , " << space.dimensions()[1] << " ,  " << space.dimensions()[2] << ")\n\n";

        //Set Periodic Box Condition

        system_openmm->setDefaultPeriodicBoxVectors(OpenMM::Vec3(Box_x_Edge_Length, 0, 0), OpenMM::Vec3(0, Box_y_Edge_Length, 0), OpenMM::Vec3(0, 0, Box_z_Edge_Length));
        openmm_context->setPeriodicBoxVectors(OpenMM::Vec3(Box_x_Edge_Length, 0, 0), OpenMM::Vec3(0, Box_y_Edge_Length, 0), OpenMM::Vec3(0, 0, Box_z_Edge_Length));
        openmm_context->reinitialize();
    }

    //OpenMM vector coordinate
    std::vector<OpenMM::Vec3> positions_openmm(nats);

    //OpenMM vector momenta
    std::vector<OpenMM::Vec3> velocities_openmm(nats);

    // Conversion factor because sire units of time are in AKMA, whereas OpenMM uses picoseconds

    double AKMAPerPs = 0.04888821;
    double PsPerAKMA = 1.0 / AKMAPerPs;

    const int nmols = ws.nMolecules();

    int system_index = 0;

    for (int i = 0; i < nmols; ++i)
    {

        const int nats_mol = ws.nAtoms(i);

        Vector *c = ws.coordsArray(i);
        Vector *p = ws.momentaArray(i);
        const double *m = ws.massArray(i);

        for (int j = 0; j < nats_mol; ++j)
        {

            positions_openmm[system_index] = OpenMM::Vec3(c[j].x() * (OpenMM::NmPerAngstrom), c[j].y() * (OpenMM::NmPerAngstrom), c[j].z() * (OpenMM::NmPerAngstrom));

            if (m[j] == 0.0)
                qDebug() << "\nWARNING - THE MASS OF PARTICLE " << system_index << " is ZERO\n";

            if (m[j] > SireMaths::small)
            {
                velocities_openmm[system_index] = OpenMM::Vec3(p[j].x() / m[j] * (OpenMM::NmPerAngstrom) * PsPerAKMA, p[j].y() / m[j] * (OpenMM::NmPerAngstrom) * PsPerAKMA, p[j].z() / m[j] * (OpenMM::NmPerAngstrom) * PsPerAKMA);
            }
            else
            {
                velocities_openmm[system_index] = OpenMM::Vec3(0.0, 0.0, 0.0);
            }

            if (Debug)
            {
                qDebug() << "Particle num = " << system_index;
                qDebug() << "Particle mass = " << m[j];
                qDebug() << "X = " << positions_openmm[system_index][0] * OpenMM::AngstromsPerNm << " A" <<
                    " Y = " << positions_openmm[system_index][1] * OpenMM::AngstromsPerNm << " A" <<
                    " Z = " << positions_openmm[system_index][2] * OpenMM::AngstromsPerNm << " A";
                qDebug() << "Vx = " << velocities_openmm[system_index][0] << " Vy = " << velocities_openmm[system_index][1] << " Vz = " << velocities_openmm[system_index][2] << "\n";
            }
            system_index++;
        }
    }

    if (system_index != nats)
    {
        if (Debug)
            qDebug() << " system_index " << system_index << " nats " << nats;
        throw SireError::program_bug(QObject::tr("The number of atoms in the openmm system does not match the number of atoms in the sire workspace"), CODELOC);
    }

    openmm_context->setPositions(positions_openmm);
    openmm_context->setVelocities(velocities_openmm);


}

void OpenMMFrEnergyST::destroyContext()
{
    if (this->isContextInitialised)
    {
        delete openmm_context;
        openmm_context = 0;
        this->isContextInitialised = false;
    }
}

MolarEnergy OpenMMFrEnergyST::getPotentialEnergy(const System &system)
{
    cout << "Energy function, does this ever get called?" << endl;
    IntegratorWorkspacePtr ws = this->createWorkspace(molgroup);
    ws.edit().setSystem(system);

    createContext(ws.edit(), 2 * femtosecond);

    int infoMask = 0;
    infoMask = infoMask + OpenMM::State::Energy;
    OpenMM::State state_openmm = openmm_context->getState(infoMask);

    MolarEnergy nrg = state_openmm.getPotentialEnergy() * kJ_per_mol;

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
System OpenMMFrEnergyST::minimiseEnergy(System &system, double tolerance = 1.0e-10, int max_iteration = 1)
{
    bool Debug = false;
    const MoleculeGroup moleculegroup = this->molgroup.read();
    IntegratorWorkspacePtr workspace = this->createWorkspace(moleculegroup);
    if (system.nMolecules() != moleculegroup.nMolecules())
    {
        std::cerr << "Number of molecules do not agree!";
        exit(1);
    }
    workspace.edit().setSystem(system);
    // Use helper function to create a Context
    SireUnits::Dimension::Time timestep = 0.0 * picosecond;
    createContext(workspace.edit(), timestep);
    // Step 2 minimise
    OpenMM::LocalEnergyMinimizer::minimize(*openmm_context, tolerance, max_iteration);
    // Step 3 update the positions in the system
    int infoMask = OpenMM::State::Positions;
    OpenMM::State state_openmm = openmm_context->getState(infoMask);
    std::vector<OpenMM::Vec3> positions_openmm = state_openmm.getPositions();
    // Recast to atomicvelocityworkspace because want to use commitCoordinates() method to update system
    AtomicVelocityWorkspace &ws = workspace.edit().asA<AtomicVelocityWorkspace>();
    const int nmols = ws.nMolecules();
    int k = 0;

    for (int i = 0; i < nmols; i++)
    {
        Vector *sire_coords = ws.coordsArray(i);
        for (int j = 0; j < ws.nAtoms(i); j++)
        {
            sire_coords[j] = Vector(positions_openmm[j + k][0] * (OpenMM::AngstromsPerNm),
                                    positions_openmm[j + k][1] * (OpenMM::AngstromsPerNm),
                                    positions_openmm[j + k][2] * (OpenMM::AngstromsPerNm));
            if (Debug)
            {
                std::cout << "X = " << positions_openmm[j + k][0] * OpenMM::AngstromsPerNm << " A" <<
                    " Y = " << positions_openmm[j + k][1] * OpenMM::AngstromsPerNm << " A" <<
                    " Z = " << positions_openmm[j + k][2] * OpenMM::AngstromsPerNm << " A";
            }
        }
        k = k + ws.nAtoms(i);
    }

    // This causes the workspace to update the system coordinates with the
    // contents of *sire_coords. Note that velocities aren't touched.
    ws.commitCoordinates();
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

System OpenMMFrEnergyST::annealSystemToLambda(System &system,
                                      SireUnits::Dimension::Time anneal_step_size,
                                      int annealing_steps)
{
    bool Debug = false;
    const double AKMAPerPs = 0.04888821;

    const MoleculeGroup moleculegroup = this->molgroup.read();
    IntegratorWorkspacePtr workspace = this->createWorkspace(moleculegroup);
    //TODO: Add some sanity checks here.
    if (system.nMolecules() != moleculegroup.nMolecules())
    {
        std::cerr << "Number of molecules in do not agree!";
        exit(1);
    }

    workspace.edit().setSystem(system);
    //SireUnits::Dimension::Time timestep = stepSize * picosecond;
    createContext(workspace.edit(), anneal_step_size);

    int max = ceil(Alchemical_value / 0.1);

    double lam = 0.0;

    for (int i = 0; i < max + 1; i++)
    {
        updateOpenMMContextLambda(lam);
        (openmm_context->getIntegrator()).step(annealing_steps);

        if (i == max - 1)
            lam = Alchemical_value;
        else
            lam = lam + 0.1;
    }
    int infoMask = OpenMM::State::Positions;
    infoMask = infoMask + OpenMM::State::Velocities;
    OpenMM::State state_openmm = openmm_context->getState(infoMask);
    std::vector<OpenMM::Vec3> positions_openmm = state_openmm.getPositions();
    std::vector<OpenMM::Vec3> velocities_openmm = state_openmm.getVelocities();

    // Recast to atomicvelocityworkspace because want to use commitCoordinates() method to update system
    AtomicVelocityWorkspace &ws = workspace.edit().asA<AtomicVelocityWorkspace>();

    const int nmols = ws.nMolecules();
    int k = 0;

    for (int i = 0; i < nmols; i++)
    {

        Vector *sire_coords = ws.coordsArray(i);
        Vector *sire_momenta = ws.momentaArray(i);
        const double *m = ws.massArray(i);

        for (int j = 0; j < ws.nAtoms(i); j++)
        {

            sire_coords[j] = Vector(positions_openmm[j + k][0] * (OpenMM::AngstromsPerNm),
                                    positions_openmm[j + k][1] * (OpenMM::AngstromsPerNm),
                                    positions_openmm[j + k][2] * (OpenMM::AngstromsPerNm));

            sire_momenta[j] = Vector(velocities_openmm[j + k][0] * m[j] * (OpenMM::AngstromsPerNm) * AKMAPerPs,
                                     velocities_openmm[j + k][1] * m[j] * (OpenMM::AngstromsPerNm) * AKMAPerPs,
                                     velocities_openmm[j + k][2] * m[j] * (OpenMM::AngstromsPerNm) * AKMAPerPs);

        }
        k = k + ws.nAtoms(i);
    }


    ws.commitCoordinatesAndVelocities();
    //Now we also want to update the systems box in case it is very different!
    if(MCBarostat_flag)
    {
        // dummy buffered dimensions vector, maybe there is better solution
        //to this than just passing an empty vector
        QVector< Vector> dimensions;
        updateBoxDimensions(state_openmm, dimensions, Debug, ws);
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
void OpenMMFrEnergyST::integrate(IntegratorWorkspace &workspace,
                                 const Symbol &nrg_component, SireUnits::Dimension::Time timestep,
                                 int nmoves, bool record_stats)
{


    createContext(workspace, timestep);
    bool Debug = false;
    const int nats = openmm_system->getNumParticles();

    AtomicVelocityWorkspace &ws = workspace.asA<AtomicVelocityWorkspace>();

    const int nmols = ws.nMolecules();

    const double AKMAPerPs = 0.04888821;

    const double dt = convertTo(timestep.value(), picosecond);

    if (Debug)
        qDebug() << " Doing " << nmoves << " steps of dynamics ";


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

    if (coord_freq > 0)
    {
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
    if (coord_freq > 0)
    {
        if (nframes > MAXFRAMES)
        {
            throw SireError::program_bug(QObject::tr("You are requesting to "
                                                     "buffer %1 frames, which is above the hardcoded limit "
                                                     "of %2.").arg(n_samples, MAXFRAMES), CODELOC);
        }
    }
    else
    {
        nframes = 0;
    }

    QVector< std::vector<OpenMM::Vec3> > buffered_positions;
    QVector< Vector > buffered_dimensions;

    OpenMM::Vec3 a;
    OpenMM::Vec3 b;
    OpenMM::Vec3 c;

    const double beta = 1.0 / (0.0083144621 * convertTo(Temperature.value(), kelvin)); //mol/kJ

    int infoMask = 0;

    infoMask = OpenMM::State::Positions;
    infoMask = infoMask + OpenMM::State::Velocities;
    infoMask = infoMask + OpenMM::State::Energy;
    //infoMask = infoMask + OpenMM::State::Parameters;

    OpenMM::State state_openmm; //OpenMM State


    int sample_count = 1;

    if (coord_freq > 0 && Debug)
        qDebug() << "Saving atom coordinates every " << coord_freq << "\n";


    if (Debug)
    {
        for (int i = 0; i < perturbed_energies.size(); i++)
            qDebug() << "Perturbed energy flag index" << i << " Value = " << perturbed_energies[i];
    }

    //OpenMM vector coordinate
    std::vector<OpenMM::Vec3> positions_openmm(nats);
    //OpenMM vector momenta
    std::vector<OpenMM::Vec3> velocities_openmm(nats);
    //Time skipping
    const double time_skip = convertTo(timeskip.value(), picosecond);

    if (time_skip != 0.0)
    {

        if (Debug)
            qDebug() << "Time to Skip = " << time_skip << "ps";

        int new_nmoves = time_skip / dt;

        if (new_nmoves >= nmoves)
        {
            throw SireError::program_bug(QObject::tr("Time to Skip is greater than the simulation time"), CODELOC);
            exit(-1);
        }

        (openmm_context->getIntegrator()).step(new_nmoves);

        n_samples = (nmoves - new_nmoves) / energy_frequency;

        if (coord_freq > 0)
            nframes = (nmoves - new_nmoves) / coord_freq;

    }


    bool IsFiniteNumber = true;


    double increment = delta_alchemical;
    double incr_plus = Alchemical_value + increment;
    double incr_minus = Alchemical_value - increment;


    double actual_gradient = 0.0;
    emptyContainers();
    while (sample_count <= n_samples)
    {
        //*********************MD STEPS****************************
        (openmm_context->getIntegrator()).step(energy_frequency);
        state_openmm = openmm_context->getState(infoMask, false, 0x01);
        double p_energy_lambda = state_openmm.getPotentialEnergy();
        if (Debug)
        {
            printf("Lambda = %f Potential energy = %.5f kcal/mol\n", Alchemical_value, p_energy_lambda * OpenMM::KcalPerKJ);
            //exit(-1);
        }
        IsFiniteNumber = (p_energy_lambda <= DBL_MAX && p_energy_lambda >= -DBL_MAX);

        if (!IsFiniteNumber)
        {
            qDebug() << "NaN or Inf has been generated along the simulation";
            exit(-1);
        }
        pot_energies.append(p_energy_lambda * OpenMM::KcalPerKJ);
        
        if (perturbed_energies[0])
        {
            openmm_context->setParameter("SPOnOff", 1.0); //Solvent-Solvent and Protein Protein Non Bonded OFF
        }
        state_openmm = openmm_context->getState(infoMask, false, 0x01);

        if (Debug)
            qDebug() << "Total Time = " << state_openmm.getTime() << " ps";

        // Because looping from 1 to n_samples
        int modulo = 1;

        if (coord_freq > 0)
            modulo = sample_count % ((coord_freq / energy_frequency));

        if (Debug)
            qDebug() << "modulo is " << modulo;

        if (coord_freq > 0 and modulo == 0)
        {
            if (Debug)
                qDebug() << "buffering coordinates and dimensions";

            positions_openmm = state_openmm.getPositions();
            buffered_positions.append(positions_openmm);

            if (MCBarostat_flag == true)
            {
                state_openmm.getPeriodicBoxVectors(a, b, c);
                Vector dims = Vector(a[0] * OpenMM::AngstromsPerNm, b[1] * OpenMM::AngstromsPerNm, c[2] * OpenMM::AngstromsPerNm);
                buffered_dimensions.append(dims);
            }
        }

        //Computing the potential energies and gradients
        p_energy_lambda = state_openmm.getPotentialEnergy();


        //Let's calculate the gradients
        double m_forward, m_backward;
        boost::tuples::tie(actual_gradient, m_forward, m_backward) = calculateGradient(incr_plus, 
                             incr_minus, p_energy_lambda, beta);

        if (alchemical_array.size()>1)
        {
            //Let's calculate the biased energies
            reduced_perturbed_energies.append(computeReducedPerturbedEnergies(beta));
        }

        //Now we append all the calculated information to the useful accumulation arrays
        finite_diff_gradients.append(actual_gradient * beta);
        forward_Metropolis.append(m_forward);
        backward_Metropolis.append(m_backward);


        //RESET coupling parameter to its original value
        if (perturbed_energies[0])
        {
            openmm_context->setParameter("SPOnOff", 0.0); //Solvent-Solvent and Protein Protein Non Bonded ON
        }
        updateOpenMMContextLambda(Alchemical_value);
        sample_count = sample_count + 1.0;

    }//end while
    if (time_skip != 0)
    {
        timeskip = SireUnits::Dimension::Time(0.0);
    }


    state_openmm = openmm_context->getState(infoMask);
    positions_openmm = state_openmm.getPositions();
    velocities_openmm = state_openmm.getVelocities();

    // Vector of Vector of molecules that are vector of atomic coordinates...
    QVector< QVector< QVector< Vector > > > buffered_workspace(nframes);
    for (int i = 0; i < buffered_workspace.size(); i++)
    {

        buffered_workspace[i].resize(nmols);

        for (int j = 0; j < nmols; j++)
        {
            int nats = ws.nAtoms(j);
            buffered_workspace[i][j].resize(nats);
        }
    }


    int k = 0;

    for (int i = 0; i < nmols; i++)
    {

        Vector *sire_coords = ws.coordsArray(i);
        Vector *sire_momenta = ws.momentaArray(i);
        const double *m = ws.massArray(i);

        for (int j = 0; j < ws.nAtoms(i); j++)
        {

            sire_coords[j] = Vector(positions_openmm[j + k][0] * (OpenMM::AngstromsPerNm),
                                    positions_openmm[j + k][1] * (OpenMM::AngstromsPerNm),
                                    positions_openmm[j + k][2] * (OpenMM::AngstromsPerNm));

            if (Debug)
                qDebug() << "X = " << positions_openmm[j + k][0] * OpenMM::AngstromsPerNm << " A" <<
                " Y = " << positions_openmm[j + k][1] * OpenMM::AngstromsPerNm << " A" <<
                " Z = " << positions_openmm[j + k][2] * OpenMM::AngstromsPerNm << " A";

            for (int l = 0; l < nframes; l++)
            {
                //qDebug() << " i " << i << " j " << j << " k " << k << " l " << l;
                Vector buffered_atcoord = Vector(buffered_positions[l][j + k][0] * (OpenMM::AngstromsPerNm),
                                                 buffered_positions[l][j + k][1] * (OpenMM::AngstromsPerNm),
                                                 buffered_positions[l][j + k][2] * (OpenMM::AngstromsPerNm));

                buffered_workspace[l][i][j] = buffered_atcoord;
            }

            sire_momenta[j] = Vector(velocities_openmm[j + k][0] * m[j] * (OpenMM::AngstromsPerNm) * AKMAPerPs,
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
    if (MCBarostat_flag == true)
    {

        updateBoxDimensions(state_openmm, buffered_dimensions, Debug, ws);
    }
    // Clear all buffers

    buffered_workspace.clear();
    buffered_dimensions.clear();
    System & ptr_sys = ws.nonConstsystem();
    ptr_sys.mustNowRecalculateFromScratch();

    return;
}

double OpenMMFrEnergyST::getPotentialEnergyAtLambda(double lambda)
{
    double curr_potential_energy = 0.0;
    int infoMask = 0;
    infoMask = infoMask + OpenMM::State::Energy;
    updateOpenMMContextLambda(lambda);
    OpenMM::State state_openmm = openmm_context->getState(infoMask);
    state_openmm = openmm_context->getState(infoMask, false, 0x01);
    curr_potential_energy = state_openmm.getPotentialEnergy();
    return curr_potential_energy;
}

void OpenMMFrEnergyST::updateOpenMMContextLambda(double lambda)
{
    //NON BONDED TERMS
    if (perturbed_energies[0])
        openmm_context->setParameter("lam", lambda); //1-5 HD
    //1-4 Interactions
    if (perturbed_energies[1])
        openmm_context->setParameter("lamhd", lambda); //1-4 HD
    if (perturbed_energies[2])
        openmm_context->setParameter("lamtd", 1.0 - lambda); //1-4 To Dummy
    if (perturbed_energies[3])
        openmm_context->setParameter("lamfd", lambda); //1-4 From Dummy
    if (perturbed_energies[4])
        openmm_context->setParameter("lamftd", lambda); //1-4 From Dummy to Dummy

    //BONDED PERTURBED TERMS
    if (perturbed_energies[5])
        openmm_context->setParameter("lambond", lambda); //Bonds
    if (perturbed_energies[6])
        openmm_context->setParameter("lamangle", lambda); //Angles
    if (perturbed_energies[7])
        openmm_context->setParameter("lamdih", lambda); //Torsions
}

boost::tuples::tuple<double, double, double> OpenMMFrEnergyST::calculateGradient(
    double incr_plus, double incr_minus, double p_energy_lambda, double beta)
{
    double double_increment = incr_plus - incr_minus;
    double gradient = 0;
    double potential_energy_lambda_plus_delta;
    double potential_energy_lambda_minus_delta;
    double forward_m;
    double backward_m;
    if (incr_plus < 1.0)
    {
        potential_energy_lambda_plus_delta = getPotentialEnergyAtLambda(incr_plus);
    }
    if (incr_minus > 0.0)
    {
        potential_energy_lambda_minus_delta = getPotentialEnergyAtLambda(incr_minus);
    }
    if (incr_minus < 0.0)
    {
        gradient = (potential_energy_lambda_plus_delta-p_energy_lambda)*2/double_increment;
        backward_m = exp(beta * (potential_energy_lambda_plus_delta - p_energy_lambda));
        forward_m = exp(-beta * (potential_energy_lambda_plus_delta - p_energy_lambda));
    }
    else if(incr_plus > 1.0)
    {
        gradient = -(potential_energy_lambda_minus_delta-p_energy_lambda)*2/double_increment;
        backward_m = exp(-beta * (potential_energy_lambda_minus_delta - p_energy_lambda));
        forward_m = exp(beta * (potential_energy_lambda_minus_delta - p_energy_lambda));
    }
    else
    {
        gradient = (potential_energy_lambda_plus_delta-potential_energy_lambda_minus_delta)/double_increment;

        backward_m = exp(-beta * (potential_energy_lambda_minus_delta - p_energy_lambda));
        forward_m = exp(-beta * (potential_energy_lambda_plus_delta - p_energy_lambda));
    }
    return boost::tuples::make_tuple(gradient, forward_m, backward_m);
}

QVector<double> OpenMMFrEnergyST::computeReducedPerturbedEnergies(double beta)
{
    bool Debug = false;
    QVector<double> perturbed;
    QVector<double>::iterator i;
    for (i=alchemical_array.begin(); i!=alchemical_array.end(); i++)
    {
        perturbed.append(getPotentialEnergyAtLambda(*i)*beta);
    }
    if (Debug)
    {
        for (i=perturbed.begin(); i!=perturbed.end(); i++)
        {
            qDebug() <<"bias is: "<<*i<<endl;
        }
    }
    return perturbed;
}

void OpenMMFrEnergyST::emptyContainers()
{
    finite_diff_gradients.clear();
    pot_energies.clear();
    forward_Metropolis.clear();
    backward_Metropolis.clear();
    reduced_perturbed_energies.clear();
}

void OpenMMFrEnergyST::updateBoxDimensions(OpenMM::State &state_openmm, 
                                             QVector< Vector> &buffered_dimensions, 
                                             bool Debug, AtomicVelocityWorkspace &ws)
{
    OpenMM::Vec3 a;
    OpenMM::Vec3 b;
    OpenMM::Vec3 c;

    state_openmm.getPeriodicBoxVectors(a, b, c);
    Vector new_dims = Vector(a[0] * OpenMM::AngstromsPerNm, b[1] * OpenMM::AngstromsPerNm, c[2] * OpenMM::AngstromsPerNm);
    if (Debug)
        qDebug() << " a " << a[0] << " b " << b[1] << " c " << c[2];

    System & ptr_sys = ws.nonConstsystem();
    PeriodicBox sp = ptr_sys.property("space").asA<PeriodicBox>();

    sp.setDimensions(new_dims);
    const QString string = "space";
    ptr_sys.setProperty(string, sp);

    /** Buffer dimensions if necessary */
    for (int k = 0; k < buffered_dimensions.size(); k++)
    {
        const QString buffered_space = "buffered_space_" + QString::number(k);
        PeriodicBox buff_space = PeriodicBox(buffered_dimensions[k]);
        ptr_sys.setProperty(buffered_space, buff_space);
    }
}

/** Get the cutoff type: nocutoff, cutoffnonperiodic, cutoffperiodic */
QString OpenMMFrEnergyST::getCutoffType(void)
{
    return CutoffType;
}

/** Set the cutoff type: nocutoff, cutoffnonperiodic, cutoffperiodic */
void OpenMMFrEnergyST::setCutoffType(QString cutoff_type)
{
    CutoffType = cutoff_type;
}

/** Get the cutoff distance in A */
SireUnits::Dimension::Length OpenMMFrEnergyST::getCutoffDistance(void)
{
    return cutoff_distance;
}

/** Set the cutoff distance in A */
void OpenMMFrEnergyST::setCutoffDistance(SireUnits::Dimension::Length distance)
{
    cutoff_distance = distance;
}

/** Get the dielectric constant */
double OpenMMFrEnergyST::getFieldDielectric(void)
{
    return field_dielectric;
}

/** Set the dielectric constant */
void OpenMMFrEnergyST::setFieldDielectric(double dielectric)
{
    field_dielectric = dielectric;
}

/** Set Andersen thermostat */

void OpenMMFrEnergyST::setAndersen(bool andersen)
{
    Andersen_flag = andersen;
}

/** Get Andersen thermostat status on/off */
bool OpenMMFrEnergyST::getAndersen(void)
{

    return Andersen_flag;

}

/** Get the Andersen Thermostat frequency collision */
double OpenMMFrEnergyST::getAndersenFrequency(void)
{
    return Andersen_frequency;
}

/** Set the Andersen Thermostat frequency collision */
void OpenMMFrEnergyST::setAndersenFrequency(double freq)
{
    Andersen_frequency = freq;
}

/** Get the Integrator random seed */
int OpenMMFrEnergyST::getRandomSeed(void)
{
    return random_seed;
}

/** Set the Integrator random seed */
void OpenMMFrEnergyST::setRandomSeed(int seed)
{
    random_seed = seed;
}

/** Get the bath Temperature */
SireUnits::Dimension::Temperature OpenMMFrEnergyST::getTemperature(void)
{
    return Temperature;
}

/** Set the Temperature */
void OpenMMFrEnergyST::setTemperature(SireUnits::Dimension::Temperature temperature)
{
    Temperature = temperature;
}

/** Set Monte Carlo Barostat on/off */

void OpenMMFrEnergyST::setMCBarostat(bool MCBarostat)
{
    MCBarostat_flag = MCBarostat;
}

/** Get Andersen thermostat status on/off */
bool OpenMMFrEnergyST::getMCBarostat(void)
{
    return MCBarostat_flag;
}

/** Get the Monte Carlo Barostat frequency in time speps */
int OpenMMFrEnergyST::getMCBarostatFrequency(void)
{
    return MCBarostat_frequency;
}

/** Set the Monte Carlo Barostat frequency in time speps */
void OpenMMFrEnergyST::setMCBarostatFrequency(int freq)
{
    MCBarostat_frequency = freq;

}

/** Get the Presure */
SireUnits::Dimension::Pressure OpenMMFrEnergyST::getPressure(void)
{
    return Pressure;
}

/** Set the Pressure */
void OpenMMFrEnergyST::setPressure(SireUnits::Dimension::Pressure pressure)
{
    Pressure = pressure;
}

/** Get the Constraint type: none, hbonds, allbonds, hangles */
QString OpenMMFrEnergyST::getConstraintType(void)
{
    return ConstraintType;
}

/** Set the Constraint type: none, hbonds, allbonds, hangles */
void OpenMMFrEnergyST::setConstraintType(QString constrain)
{
    ConstraintType = constrain;
}

/** Get the OpenMMMD Platform: CUDA, OpenCL, CPU */
QString OpenMMFrEnergyST::getPlatform(void)
{
    return platform_type;
}

/** Set the OpenMM Platform: CUDA, OpenCL, CPU */
void OpenMMFrEnergyST::setPlatform(QString platform)
{
    platform_type = platform;
}

/** Get the OpenMMMD Platform: CUDA, OpenCL, CPU */
QString OpenMMFrEnergyST::getDeviceIndex(void)
{
    return device_index;
}

/** Set the OpenMM Precision */
void OpenMMFrEnergyST::setPrecision(QString prec)
{
    precision = prec;
}

/** Get the OpenMMMD Precision */
QString OpenMMFrEnergyST::getPrecision(void)
{
    return precision;
}

/** Set the OpenMM Platform: CUDA, OpenCL, CPU */
void OpenMMFrEnergyST::setDeviceIndex(QString deviceidx)
{
    device_index = deviceidx;
}

/** Get the Restaint mode*/
bool OpenMMFrEnergyST::getRestraint(void)
{
    return Restraint_flag;
}

/** Set the Retraint mode */
void OpenMMFrEnergyST::setRestraint(bool Restraint)
{
    Restraint_flag = Restraint;
}

/** Get the Center of Mass motion removal frequency */
int OpenMMFrEnergyST::getCMMremovalFrequency(void)
{
    return CMMremoval_frequency;
}

/** Set the Center of Mass motion removal frequency */
void OpenMMFrEnergyST::setCMMremovalFrequency(int frequency)
{
    CMMremoval_frequency = frequency;
}

/** Get the frequency of buffering coordinates */
int OpenMMFrEnergyST::getBufferFrequency()
{
    return buffer_frequency;
}

/** Set the Center of Mass motion removal frequency */
void OpenMMFrEnergyST::setBufferFrequency(int frequency)
{
    buffer_frequency = frequency;
}

/** Get the frequency of buffering coordinates */
int OpenMMFrEnergyST::getEnergyFrequency()
{
    return energy_frequency;
}

/** Set the Center of Mass motion removal frequency */
void OpenMMFrEnergyST::setEnergyFrequency(int frequency)
{
    energy_frequency = frequency;
}

/** Get the alchemical value used to calculate the free energy change via TI method*/
double OpenMMFrEnergyST::getAlchemicalValue(void)
{
    return Alchemical_value;
}

/** Set the alchemical value used to calculate the free energy change via TI method*/
void OpenMMFrEnergyST::setAlchemicalValue(double lambda_value)
{
    Alchemical_value = max(0.0, min(1.0, lambda_value));
}

void OpenMMFrEnergyST::setAlchemicalArray(QVector<double> lambda_array)
{
    for (int i =0; i< lambda_array.size(); i++)
    {
        alchemical_array.append(max(0.0, min(1.0, lambda_array[i])));
    }
}

/** Get the coulomb power used in the soft core potential*/
//int OpenMMFrEnergyST::getCoulomb_power(void)

float OpenMMFrEnergyST::getCoulombPower(void)
{
    return coulomb_power;
}

/** Set the coulomb power used in the soft core potential*/
//void OpenMMFrEnergyST::setCoulomb_power(int coulomb)

void OpenMMFrEnergyST::setCoulombPower(float coulomb)
{
    coulomb_power = coulomb;
}

/** Get the shift used in the soft core potential*/
double OpenMMFrEnergyST::getShiftDelta(void)
{
    return shift_delta;
}

/**
 * <Set the shift used in the soft core potential>
 * @param shiftdelta
 */
void OpenMMFrEnergyST::setShiftDelta(double shiftdelta)
{

    shift_delta = shiftdelta;

}

/** Get the delta alchemical used in the FEP method*/
double OpenMMFrEnergyST::getDeltaAlchemical(void)
{
    return delta_alchemical;
}

/**
 * Set the delta alchemical used in the FEP method
 * @param deltaalchemical
 */
void OpenMMFrEnergyST::setDeltatAlchemical(double deltaalchemical)
{
    delta_alchemical = deltaalchemical;
}

/** Calculated Gradients*/
QVector<double> OpenMMFrEnergyST::getGradients(void)
{
    return finite_diff_gradients;
}

/** Average energies*/
QVector<double> OpenMMFrEnergyST::getEnergies(void)
{
    return pot_energies;
}
/** Average energies*/
QVector<double> OpenMMFrEnergyST::getForwardMetropolis(void)
{
    return forward_Metropolis;
}
/** Average energies*/
QVector<double> OpenMMFrEnergyST::getBackwardMetropolis(void)
{
    return backward_Metropolis;
}


QVector<QVector <double> > OpenMMFrEnergyST::getReducedPerturbedEnergies(void)
{
    return reduced_perturbed_energies;
}

/** Get the Integrator type*/
QString OpenMMFrEnergyST::getIntegrator(void)
{
    return Integrator_type;
}

/** Set the Integrator type*/
void OpenMMFrEnergyST::setIntegrator(QString intgrator)
{
    Integrator_type = intgrator;
}

/** Get the friction used in specific Integrator type*/
SireUnits::Dimension::Time OpenMMFrEnergyST::getFriction(void)
{
    return friction;
}

/** Set the friction used in specific Integrator type*/
void OpenMMFrEnergyST::setFriction(SireUnits::Dimension::Time thefriction)
{
    friction = thefriction;
}

/** Get the integration tolerance */
double OpenMMFrEnergyST::getIntegrationTolerance(void)
{
    return integration_tol;
}

/** Set the integration tolerance*/
void OpenMMFrEnergyST::setIntegrationTolerance(double tolerance)
{
    integration_tol = tolerance;
}

/** Get total time to skip*/
SireUnits::Dimension::Time OpenMMFrEnergyST::getTimetoSkip(void)
{
    return timeskip;
}

/** Get total time to skip*/
void OpenMMFrEnergyST::setTimetoSkip(SireUnits::Dimension::Time skip)
{
    timeskip = skip;
}

/** Set the flag to reinitialise the context*/
void OpenMMFrEnergyST::setReinitialiseContext(bool reinitialise)
{
    reinitialise_context = reinitialise;
}

/** Create an empty workspace */
IntegratorWorkspacePtr OpenMMFrEnergyST::createWorkspace(const PropertyMap &map) const
{
    return IntegratorWorkspacePtr(new AtomicVelocityWorkspace(map));
}

/** Return the ensemble of this integrator */
Ensemble OpenMMFrEnergyST::ensemble() const
{
    return Ensemble::NVE();
}

/** Return whether or not this integrator is time-reversible */
bool OpenMMFrEnergyST::isTimeReversible() const
{
    return true;
}

/** Create a workspace for this integrator for the molecule group 'molgroup' */
IntegratorWorkspacePtr OpenMMFrEnergyST::createWorkspace(const MoleculeGroup &molgroup, const PropertyMap &map) const
{
    return IntegratorWorkspacePtr(new AtomicVelocityWorkspace(molgroup, map));
}

const char* OpenMMFrEnergyST::typeName()
{
    return QMetaType::typeName(qMetaTypeId<OpenMMFrEnergyST>());
}
