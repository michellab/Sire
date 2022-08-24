/********************************************\
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

#include "openmmmdintegrator.h"
#include "ensemble.h"

#include "SireMol/moleculegroup.h"
#include "SireMol/partialmolecule.h"
#include "SireMol/molecule.h"
#include "SireMol/atommasses.h"
#include "SireMol/atomcoords.h"
#include "SireMol/moleditor.h"
#include "SireMol/core.h"

#include "SireBase/variantproperty.h"

#include "SireMol/amberparameters.h"

#include "SireSystem/system.h"

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
#include <iostream>
/* include <QElapsedTimer> */
#include <QDebug>
#include <QTime>

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

//ADDED BY JM
using namespace SireMM;
using namespace SireIO;
//using namespace SireMaths;

//ADDED BY GAC
using namespace std;

enum
{
    NOCUTOFF = 0,
    CUTOFFNONPERIODIC = 1,
    CUTOFFPERIODIC = 2,
    EWALD = 3,
    PME = 4
};

enum
{
    NONE = 0,
    HBONDS = 1,
    ALLBONDS = 2,
    HANGLES = 3

};

static const RegisterMetaType<OpenMMMDIntegrator> r_openmmint;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const OpenMMMDIntegrator &velver)
{
    writeHeader(ds, r_openmmint, 2);

    SharedDataStream sds(ds);

    sds << velver.frequent_save_velocities << velver.molgroup
        << velver.Integrator_type << velver.friction
        << velver.CutoffType << velver.cutoff_distance << velver.field_dielectric
        << velver.tolerance_ewald_pme
        << velver.Andersen_flag << velver.Andersen_frequency
        << velver.MCBarostat_flag << velver.MCBarostat_frequency << velver.ConstraintType
        << velver.Pressure << velver.Temperature
        << velver.platform_type << velver.Restraint_flag << velver.CMMremoval_frequency
        << velver.buffer_frequency
        << velver.device_index << velver.LJ_dispersion << velver.precision << velver.integration_tol
        << velver.timeskip
        << velver.reinetialise_context
        << velver.is_periodic
        << static_cast<const Integrator&> (velver);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, OpenMMMDIntegrator &velver)
{

    VersionID v = readHeader(ds, r_openmmint);

    if (v == 2)
    {
        SharedDataStream sds(ds);

        velver.is_periodic = false;

        sds >> velver.frequent_save_velocities >> velver.molgroup
            >> velver.Integrator_type >> velver.friction
            >> velver.CutoffType >> velver.cutoff_distance >> velver.field_dielectric
            >> velver.tolerance_ewald_pme
            >> velver.Andersen_flag >> velver.Andersen_frequency
            >> velver.MCBarostat_flag >> velver.MCBarostat_frequency >> velver.ConstraintType
            >> velver.Pressure >> velver.Temperature
            >> velver.platform_type >> velver.Restraint_flag >> velver.CMMremoval_frequency
            >> velver.buffer_frequency
            >> velver.device_index >> velver.LJ_dispersion >> velver.precision
            >> velver.integration_tol
            >> velver.timeskip
            >> velver.reinetialise_context
            >> velver.is_periodic
            >> static_cast<Integrator&> (velver);

        // Maybe....need to reinitialise from molgroup because openmm system was not serialised...
        velver.isSystemInitialised = false;
        velver.isContextInitialised = false;

        //qDebug() << " Re-initialisation of openmmmdintegrator from datastream";

        velver.initialise();
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);

        velver.is_periodic = false;

        sds >> velver.frequent_save_velocities >> velver.molgroup
            >> velver.Integrator_type >> velver.friction
            >> velver.CutoffType >> velver.cutoff_distance >> velver.field_dielectric
            >> velver.tolerance_ewald_pme
            >> velver.Andersen_flag >> velver.Andersen_frequency
            >> velver.MCBarostat_flag >> velver.MCBarostat_frequency >> velver.ConstraintType
            >> velver.Pressure >> velver.Temperature
            >> velver.platform_type >> velver.Restraint_flag >> velver.CMMremoval_frequency
            >> velver.buffer_frequency
            >> velver.device_index >> velver.LJ_dispersion >> velver.precision
            >> velver.integration_tol
            >> velver.timeskip
            >> velver.reinetialise_context
            >> static_cast<Integrator&> (velver);

        // Maybe....need to reinitialise from molgroup because openmm system was not serialised...
        velver.isSystemInitialised = false;
        velver.isContextInitialised = false;

        //qDebug() << " Re-initialisation of openmmmdintegrator from datastream";

        velver.initialise();
    }
    else
        throw version_error(v, "1,2", r_openmmint, CODELOC);

    return ds;
}

/** Constructor*/
OpenMMMDIntegrator::OpenMMMDIntegrator(bool frequent_save)
: ConcreteProperty<OpenMMMDIntegrator, Integrator>(),
frequent_save_velocities(frequent_save), molgroup(MoleculeGroup()),
openmm_system(0), openmm_context(0), isSystemInitialised(false),
isContextInitialised(false),
Integrator_type("leapfrogverlet"), friction(1.0 / picosecond),
CutoffType("nocutoff"), cutoff_distance(1.0 * nanometer), field_dielectric(78.3),
tolerance_ewald_pme(0.0001),
Andersen_flag(false), Andersen_frequency(90.0), MCBarostat_flag(false),
MCBarostat_frequency(25), ConstraintType("none"),
Pressure(1.0 * bar), Temperature(300.0 * kelvin), platform_type("Reference"),
Restraint_flag(false),
CMMremoval_frequency(0), buffer_frequency(0), device_index("0"),
LJ_dispersion(true), precision("single"),
reinetialise_context(false), integration_tol(0.001), timeskip(0.0 * picosecond),
//minimise(false), minimise_tol(1.0), minimise_iterations(0),
//equilib_iterations(5000), equilib_time_step(0.0005 * picosecond),
is_periodic(false)
{
}

/** Constructor using the passed molecule group */
OpenMMMDIntegrator::OpenMMMDIntegrator(const MoleculeGroup &molecule_group, bool frequent_save)
: ConcreteProperty<OpenMMMDIntegrator, Integrator>(),
frequent_save_velocities(frequent_save), molgroup(molecule_group),
openmm_system(0), openmm_context(0), isSystemInitialised(false),
isContextInitialised(false),
Integrator_type("leapfrogverlet"), friction(1.0 / picosecond),
CutoffType("nocutoff"), cutoff_distance(1.0 * nanometer), field_dielectric(78.3),
tolerance_ewald_pme(0.0001),
Andersen_flag(false), Andersen_frequency(90.0), MCBarostat_flag(false),
MCBarostat_frequency(25), ConstraintType("none"),
Pressure(1.0 * bar), Temperature(300.0 * kelvin), platform_type("Reference"),
Restraint_flag(false),
CMMremoval_frequency(0), buffer_frequency(0), device_index("0"),
LJ_dispersion(true), precision("single"),
reinetialise_context(false), integration_tol(0.001), timeskip(0.0 * picosecond),
//minimise(false), minimise_tol(1.0), minimise_iterations(0),
//equilib_iterations(5000), equilib_time_step(0.0005 * picosecond),
is_periodic(false)
{
}

/** Copy constructor */
OpenMMMDIntegrator::OpenMMMDIntegrator(const OpenMMMDIntegrator &other)
: ConcreteProperty<OpenMMMDIntegrator, Integrator>(other),
frequent_save_velocities(other.frequent_save_velocities),
molgroup(other.molgroup),
openmm_system(other.openmm_system), openmm_context(other.openmm_context),
isSystemInitialised(other.isSystemInitialised),
isContextInitialised(other.isContextInitialised),
Integrator_type(other.Integrator_type), friction(other.friction),
CutoffType(other.CutoffType), cutoff_distance(other.cutoff_distance),
field_dielectric(other.field_dielectric),
tolerance_ewald_pme(other.tolerance_ewald_pme),
Andersen_flag(other.Andersen_flag),
Andersen_frequency(other.Andersen_frequency),
MCBarostat_flag(other.MCBarostat_flag),
MCBarostat_frequency(other.MCBarostat_frequency),
ConstraintType(other.ConstraintType),
Pressure(other.Pressure), Temperature(other.Temperature),
platform_type(other.platform_type),
Restraint_flag(other.Restraint_flag),
CMMremoval_frequency(other.CMMremoval_frequency),
buffer_frequency(other.buffer_frequency), device_index(other.device_index),
LJ_dispersion(other.LJ_dispersion), precision(other.precision),
reinetialise_context(other.reinetialise_context),
integration_tol(other.integration_tol), timeskip(other.timeskip),
is_periodic(other.is_periodic)
{
}

/** Destructor */
OpenMMMDIntegrator::~OpenMMMDIntegrator()
{
    //delete openmm_system;
}

/** Copy assignment operator */
OpenMMMDIntegrator& OpenMMMDIntegrator::operator=(const OpenMMMDIntegrator &other)
{
    Integrator::operator=(other);
    frequent_save_velocities = other.frequent_save_velocities;
    molgroup = other.molgroup;
    openmm_system = other.openmm_system;
    openmm_context = other.openmm_context;
    isSystemInitialised = other.isSystemInitialised;
    isContextInitialised = other.isContextInitialised;
    Integrator_type = other.Integrator_type;
    friction = other.friction;
    CutoffType = other.CutoffType;
    cutoff_distance = other.cutoff_distance;
    field_dielectric = other.field_dielectric;
    tolerance_ewald_pme = other.tolerance_ewald_pme;
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
    device_index = other.device_index;
    LJ_dispersion = other.LJ_dispersion;
    precision = other.precision;
    reinetialise_context = other.reinetialise_context;
    integration_tol = other.integration_tol;
    timeskip = other.timeskip;
    is_periodic = other.is_periodic;

    return *this;
}

/** Comparison operator */
bool OpenMMMDIntegrator::operator==(const OpenMMMDIntegrator &other) const
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
        and device_index == other.device_index
        and precision == other.precision
        and Integrator_type == other.Integrator_type
        and friction == other.friction
        and integration_tol == other.integration_tol
        and timeskip == other.timeskip
        and reinetialise_context == other.reinetialise_context
        and is_periodic == other.is_periodic
        and Integrator::operator==(other);
}

/** Comparison operator */
bool OpenMMMDIntegrator::operator!=(const OpenMMMDIntegrator &other) const
{
    return not OpenMMMDIntegrator::operator==(other);
}

/** Return a string representation of this integrator */
QString OpenMMMDIntegrator::toString() const
{
    return QObject::tr("OpenMMMDIntegrator()");
}

/** Integrate the coordinates of the atoms in the molecules in 'molgroup'
    using the forces in 'forcetable', using the optionally supplied
    property map to find the necessary molecular properties

    \throw SireMol::missing_molecule
    \throw SireBase::missing_property
    \throw SireError:invalid_cast
    \throw SireError::incompatible_error
 */

void OpenMMMDIntegrator::initialise()
{
    bool Debug = false;

    if (Debug)
        qDebug() << " initialising OpenMMMDIntegrator";

    // Create a workspace using the stored molgroup


    const MoleculeGroup moleculegroup = this->molgroup.read();

    if (moleculegroup.isEmpty())
    {
        throw SireError::program_bug(QObject::tr(
                                                 "Cannot initialise OpenMMMDIntegrator because molgroup has not been defined"), CODELOC);
    }

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
    else if (CutoffType == "ewald")
        flag_cutoff = EWALD;
    else if (CutoffType == "PME")
        flag_cutoff = PME;
    else
        throw SireError::program_bug(QObject::tr(
                                                 "The CutOff method has not been specified. Possible choices: nocutoff, cutoffnonperiodic, cutoffperiodic,ewal,PME"), CODELOC);

    if (Debug)
        qDebug() << "\nCutoffType = " << CutoffType << "\n";

    if (ConstraintType == "none")
        flag_constraint = NONE;
    else if (ConstraintType == "hbonds")
        flag_constraint = HBONDS;
    else if (ConstraintType == "allbonds")
        flag_constraint = ALLBONDS;
    else if (ConstraintType == "hangles")
        flag_constraint = HANGLES;
    else
        throw SireError::program_bug(QObject::tr(
                                                 "The Constraints method has not been specified. Possible choices: none, hbonds, allbonds, hangles"), CODELOC);

    if (Debug)
        qDebug() << "\nConstraint Type = " << ConstraintType << "\n";

    //Load Plugins from the OpenMM standard Plugin Directory
    OpenMM::Platform::loadPluginsFromDirectory(OpenMM::Platform::getDefaultPluginsDirectory());

    OpenMM::System *system_openmm = new OpenMM::System();

    system_openmm->setDefaultPeriodicBoxVectors(OpenMM::Vec3(4, 0, 0),
                                                OpenMM::Vec3(0, 5, 0),
                                                OpenMM::Vec3(0, 0, 6));

    OpenMM::NonbondedForce * nonbond_openmm = new OpenMM::NonbondedForce();

    system_openmm->addForce(nonbond_openmm);

    if (flag_cutoff == NOCUTOFF)
    {
        nonbond_openmm->setNonbondedMethod(OpenMM::NonbondedForce::NoCutoff);

        if (Debug)
            qDebug() << "\nCut off type = " << CutoffType << "\n";
    }
    else
    {
        const double converted_cutoff_distance = convertTo(cutoff_distance.value(), nanometer);

        if (flag_cutoff == CUTOFFNONPERIODIC)
        {
            nonbond_openmm->setNonbondedMethod(OpenMM::NonbondedForce::CutoffNonPeriodic);
            //Set Dielectric constant media
            nonbond_openmm->setReactionFieldDielectric(field_dielectric);
        }
        else if (flag_cutoff == CUTOFFPERIODIC)
        {
            nonbond_openmm->setNonbondedMethod(OpenMM::NonbondedForce::CutoffPeriodic);
            //Set Dielectric constant media
            nonbond_openmm->setReactionFieldDielectric(field_dielectric);
        }
        else if (flag_cutoff == EWALD)
        {
            nonbond_openmm->setNonbondedMethod(OpenMM::NonbondedForce::Ewald);
            nonbond_openmm->setEwaldErrorTolerance(tolerance_ewald_pme);
        }
        else if (flag_cutoff == PME)
        {
            nonbond_openmm->setNonbondedMethod(OpenMM::NonbondedForce::PME);
            nonbond_openmm->setEwaldErrorTolerance(tolerance_ewald_pme);
        }

        nonbond_openmm->setCutoffDistance(converted_cutoff_distance);

        if (Debug)
        {
            qDebug() << "\nCut off type = " << CutoffType << "\n";
            qDebug() << "CutOff distance = " << converted_cutoff_distance << " Nm" << "\n";
            if (flag_cutoff == CUTOFFNONPERIODIC || flag_cutoff == CUTOFFPERIODIC)
                qDebug() << "Dielectric constant = " << field_dielectric << "\n\n";
            else if (flag_cutoff == EWALD || flag_cutoff == PME)
                qDebug() << "Tolerance EWALD/PME = " << tolerance_ewald_pme << "\n\n";
        }
    }


    // Andersen thermostat. Complain if NOT using Verlet
    if (Andersen_flag == true)
    {
        if (Debug)
            qDebug() << " Integrator_type " << Integrator_type;

        if (Integrator_type != "leapfrogverlet" and Integrator_type != "variableleapfrogverlet")
        {
            throw SireError::program_bug(QObject::tr(
                                                     "The Andersen thermostat can only be used with the leapfrogverlet or variableleapfrogverlet integrators"), CODELOC);
        }

        const double converted_Temperature = convertTo(Temperature.value(), kelvin);

        OpenMM::AndersenThermostat * thermostat = new OpenMM::AndersenThermostat(converted_Temperature, Andersen_frequency);

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
        system_openmm->addForce(barostat);

        if (Debug)
        {
            qDebug() << "\nMonte Carlo Barostat set\n";
            qDebug() << "Temperature = " << converted_Temperature << " K\n";
            qDebug() << "Pressure = " << converted_Pressure << " bar\n";
            qDebug() << "Frequency every " << MCBarostat_frequency << " steps\n";
            qDebug() << "Lennard Jones Dispersion term is set to " << LJ_dispersion << "\n";
        }
    }


    //Setting Lennard Jones dispersion globally, since its default is set to true!
    nonbond_openmm->setUseDispersionCorrection(LJ_dispersion);

    //OpenMM Bonded Forces

    OpenMM::HarmonicBondForce * bondStretch_openmm = new OpenMM::HarmonicBondForce();

    OpenMM::HarmonicAngleForce * bondBend_openmm = new OpenMM::HarmonicAngleForce();

    OpenMM::PeriodicTorsionForce * bondTorsion_openmm = new OpenMM::PeriodicTorsionForce();

    system_openmm->addForce(bondStretch_openmm);

    system_openmm->addForce(bondBend_openmm);

    system_openmm->addForce(bondTorsion_openmm);

    // Check whether positional restraints have been defined for a set of atoms in that molecule.
    // You can get the information out by getting the property and casting to VariantProperty
    //From VariantProperty you have the QVariant, so you can call .toDouble() and .toInt() there
    //so VariantProperty num = mol.property(QString("AtomNum(%1)").arg(i)).asA<VariantProperty>();
    //AtomNum atomnum( num.toInt() );
    //double x = mol.property(QString("x(%1)").arg(i)).asA<VariantProperty>().toDouble();
    //double y = ...; double z = ...;
    //QVector< QPair<AtomNum,Vector> > vals;
    //vals.append( QPair<AtomNum,Vector>(AtomNum(num.toInt()), Vector(x,y,z) ) );

    OpenMM::CustomExternalForce * positionalRestraints_openmm = NULL;

    if (Restraint_flag == true)
    {
        positionalRestraints_openmm = new OpenMM::CustomExternalForce("k*( (x-xref)^2 + (y-yref)^2 + (z-zref)^2 )");
        positionalRestraints_openmm->addPerParticleParameter("xref");
        positionalRestraints_openmm->addPerParticleParameter("yref");
        positionalRestraints_openmm->addPerParticleParameter("zref");
        positionalRestraints_openmm->addPerParticleParameter("k");

        system_openmm->addForce(positionalRestraints_openmm);

        if (Debug)
            qDebug() << "\nRestraint = ON\n\n";

    }

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

    // Conversion factor because sire units of time are in AKMA, whereas OpenMM uses picoseconds
    //double AKMAPerPs = 0.04888821;
    //double PsPerAKMA = 1 / AKMAPerPs;

    // The default 1,4 scaling factors
    const double Coulomb14Scale = 1.0 / 1.2;
    const double LennardJones14Scale = 1.0 / 2.0;

    // A list of 1,4 atom pairs with non default scale factors
    // for each entry, first pair has pair of indices, second has pair of scale factors
    QList< QPair< QPair<int, int>, QPair<double, double > > > custom14pairs;

    for (int i = 0; i < nmols; ++i)
    {
        const int nats_mol = ws.nAtoms(i);

        //Vector *c = ws.coordsArray(i);
        //Vector *p = ws.momentaArray(i);
        const double *m = ws.massArray(i);

        MolNum molnum = moleculegroup.molNumAt(i);
        const ViewsOfMol &molview = moleculegroup[molnum].data();
        const Molecule &mol = molview.molecule();
        Selector<Atom> molatoms = mol.atoms();

        for (int j = 0; j < nats_mol; ++j)
        {

            system_openmm->addParticle(m[j]);

            Atom at = molatoms(j);
            AtomNum atnum = at.number();

            //if (Debug)
                //qDebug() << " openMM_index " << system_index << " Sire Atom Number " << atnum.toString();

            AtomNumToOpenMMIndex[atnum.value()] = system_index;

            // JM Nov 12
            // The code below implements a ThreeParticleAverageSite for virtual sites for EPW atoms present in a WAT residue
            // This is very AMBER specific.

            AtomName atname = at.name();

            if (Debug)
                //qDebug() << " atname " << atname.value() << " mol " << i;

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

                        if ((at0name == AtomName("O") and at1name == AtomName("H1")) or
                            ( at0name == AtomName("H1") and at1name == AtomName("O")))
                        {
                            distoh = r0;
                        }
                        else if ((at0name == AtomName("H1") and at1name == AtomName("H2")) or
                                 ( at0name == AtomName("H2") and at1name == AtomName("H1")))
                        {
                            disthh = r0;
                        }
                        else if ((at0name == AtomName("EPW") and at1name == AtomName("O")) or
                                 ( at0name == AtomName("O") and at1name == AtomName("EPW")))
                        {
                            distoe = r0;
                        }
                    }

                    if (distoh < 0 or disthh < 0 or distoe < 0)
                    {
                        throw SireError::program_bug(QObject::tr("Could not find expected atoms in TIP4P water molecule."), CODELOC);
                    }

                    //qDebug() << " distoe " << distoe << " distoh " << distoh << " disthh " << disthh;

                    double weightH = distoe / sqrt((distoh * distoh) - (0.25 * disthh * disthh));

                    int o_index = AtomNumToOpenMMIndex[oatom.number().value()];
                    int h1_index = AtomNumToOpenMMIndex[h1atom.number().value()];
                    int h2_index = AtomNumToOpenMMIndex[h2atom.number().value()];

                    //if (Debug)
                        //qDebug() << "virtual site " << system_index <<
                        //" o " << o_index << " h1 " << h1_index <<
                        //" h2 " << h2_index << " 1 - weightH " << 1 - weightH <<
                        //" weightH/2 " << weightH / 2;

                    OpenMM::ThreeParticleAverageSite * vsite = new OpenMM::ThreeParticleAverageSite(o_index, h1_index, h2_index, 1 - weightH, weightH / 2, weightH / 2);

                    system_openmm->setVirtualSite(system_index, vsite);

                }
            }
            system_index = system_index + 1;
        }// end of loop on atoms in molecule
    }//end of loop on molecules in workspace

    int num_atoms_till_i = 0;

    for (int i = 0; i < nmols; i++)
    {
        const Vector *c = ws.coordsArray(i);

        Molecule molecule = moleculegroup.moleculeAt(i).molecule();
        int num_atoms_molecule = molecule.nAtoms();

        // The atomic parameters
        AtomLJs atomvdws = molecule.property("LJ").asA<AtomLJs>();
        AtomCharges atomcharges = molecule.property("charge").asA<AtomCharges>();
        QVector<SireMM::LJParameter> ljparameters = atomvdws.toVector();
        QVector<SireUnits::Dimension::Charge> charges = atomcharges.toVector();

        for (int j = 0; j < ljparameters.size(); j++)
        {
            double sigma = ljparameters[j].sigma();
            double epsilon = ljparameters[j].epsilon();
            double charge = charges[j].value();
            nonbond_openmm->addParticle(charge, sigma * OpenMM::NmPerAngstrom, epsilon * OpenMM::KJPerKcal);
        }


        if (Restraint_flag == true)
        {
            bool hasRestrainedAtoms = molecule.hasProperty("restrainedatoms");

            if (hasRestrainedAtoms)
            {
                Properties restrainedAtoms = molecule.property("restrainedatoms").asA<Properties>();

                int nrestrainedatoms = restrainedAtoms.property(QString("nrestrainedatoms")).asA<VariantProperty>().toInt();
                //if (Debug)
                    //qDebug() << " nrestrainedatoms " << nrestrainedatoms;

                for (int i = 0; i < nrestrainedatoms; i++)
                {
                    int atomnum = restrainedAtoms.property(QString("AtomNum(%1)").arg(i)).asA<VariantProperty>().toInt();
                    double xref = restrainedAtoms.property(QString("x(%1)").arg(i)).asA<VariantProperty>().toDouble();
                    double yref = restrainedAtoms.property(QString("y(%1)").arg(i)).asA<VariantProperty>().toDouble();
                    double zref = restrainedAtoms.property(QString("z(%1)").arg(i)).asA<VariantProperty>().toDouble();
                    double k = restrainedAtoms.property(QString("k(%1)").arg(i)).asA<VariantProperty>().toDouble();

                    int openmmindex = AtomNumToOpenMMIndex[atomnum];

                    //if (Debug)
                       // qDebug() << " atomnum " << atomnum << " openmmindex " << openmmindex << " x " << xref << " y " << yref << " z " << zref << " k " << k;

                    int posrestrdim = 4;
                    std::vector<double> params(posrestrdim);

                    params[0] = xref * OpenMM::NmPerAngstrom;
                    params[1] = yref * OpenMM::NmPerAngstrom;
                    params[2] = zref * OpenMM::NmPerAngstrom;
                    params[3] = k * (OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);

                    positionalRestraints_openmm->addParticle(openmmindex, params);
                }
            }
        }//end of restraint flag



        // The bonded parameters
        bool hasConnectivity = molecule.hasProperty("connectivity");

        if (!hasConnectivity)
        {
            num_atoms_till_i = num_atoms_till_i + num_atoms_molecule;
           // if (Debug)
            //{
               // qDebug() << "\nAtoms = " << num_atoms_molecule << " Num atoms till i =" << num_atoms_till_i << "\n";
               // qDebug() << "\n*********************MONOATOMIC MOLECULE DETECTED**************************\n";
            //}
            continue;
        }

        // The bonded parameters are stored in "amberparameters"
        AmberParameters amber_params = molecule.property("amberparameters").asA<AmberParameters>();

        //Bonds

        QList<BondID> bonds_ff = amber_params.getAllBonds();
        QVector<BondID> bonds = bonds_ff.toVector();

        for (int j = 0; j < bonds_ff.length(); j++)
        {
            BondID bond_ff = bonds_ff[j];
            QList<double> bond_params = amber_params.getParams(bond_ff);

            double k = bond_params[0];
            double r0 = bond_params[1];

            int idx0 = bonds[j].atom0().asA<AtomIdx>().value();
            int idx1 = bonds[j].atom1().asA<AtomIdx>().value();

            //Select the atom type
            QString atom0 = molecule.atom(AtomIdx(idx0)).toString();
            QString atom1 = molecule.atom(AtomIdx(idx1)).toString();

            // JM 04/14 Skip constraints involving atom named EPW. Necessary for EPW support in T4P in OpenMM >6.0
            AtomName atname0 = molecule.atom(bonds[j].atom0()).name();
            AtomName atname1 = molecule.atom(bonds[j].atom1()).name();
            // A more general fix may be to skip constraints involving an atom with mass 0 g.mol-1
            // But need to figure out how to get this info from Sire and should forbid changes to masses between
            // initialisation and integration...
            //const PropertyName &mass_property = PropertyName("mass");
            //const double mass0 = molecule.atom(AtomIdx(idx0)).property(mass_property).asA<AtomMasses>.value();
            //const double mass1 = molecule.atom(AtomIdx(idx1)).property(mass_property).asA<AtomMasses>.value();

            idx0 = idx0 + num_atoms_till_i;
            idx1 = idx1 + num_atoms_till_i;

            if (flag_constraint == NONE)
            {
                bondStretch_openmm->addBond(idx0, idx1, r0 * OpenMM::NmPerAngstrom,
                                            k * 2.0 * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
                //if (Debug)
                   // qDebug() << "\nBOND ADDED TO " << atom0 << " AND " << atom1 << "\n";
            }
            else if (flag_constraint == ALLBONDS || flag_constraint == HANGLES)
            {
                //if (Debug)
                   // qDebug() << " atname0 " << atname0.toString() << " atname1 " << atname1.toString();
                if (atname0 != AtomName("EPW") && atname1 != AtomName("EPW"))
                {
                    system_openmm->addConstraint(idx0, idx1, r0 * OpenMM::NmPerAngstrom);
                    //if (Debug)
                    //    qDebug() << "\nALLBONDS or HANGLES ADDED BOND CONSTRAINT TO " << atom0 << " AND " << atom1 << "\n";
                }
            }
            else if (flag_constraint == HBONDS)
            {
                if ((atom0[6] == 'H') || (atom1[6] == 'H'))
                {
                    system_openmm->addConstraint(idx0, idx1, r0 * OpenMM::NmPerAngstrom);
                    //if (Debug)
                     //   qDebug() << "\nHBONDS ADDED BOND CONSTRAINT TO " << atom0 << " AND " << atom1 << "\n";
                }
                else
                {
                    bondStretch_openmm->addBond(idx0, idx1, r0 * OpenMM::NmPerAngstrom,
                                                k * 2.0 * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
                    //if (Debug)
                    //    qDebug() << "\nHBONDS ADDED BOND TO " << atom0 << " AND " << atom1 << "\n";
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
                    //if (Debug)
                    //    qDebug() << "\nHANGLES ANGLE CONSTRAINT TO " << atom0 << " AND " << atom2 << "\n";
                }
                else if (((atom0[6] == 'H') && (atom1[6] == 'O')) || ((atom1[6] == 'O') && (atom2[6] == 'H')))
                {
                    system_openmm->addConstraint(idx0, idx2, diff.length() * OpenMM::NmPerAngstrom);
                    //if (Debug)
                    //    qDebug() << "\n ANGLE CONSTRAINT TO " << atom0 << " AND " << atom2 << "\n";
                }
                else
                {
                    bondBend_openmm->addAngle(idx0, idx1, idx2, theta0, k * 2.0 * OpenMM::KJPerKcal);
                    //if (Debug)
                    //    qDebug() << "\n ANGLE TO " << atom0 << " AND " << atom1 << " AND " << atom2 << "\n";
                }
            }
            else
            {
                bondBend_openmm->addAngle(idx0, idx1, idx2, theta0, k * 2.0 * OpenMM::KJPerKcal);
                //if (Debug)
                //    qDebug() << "\n ANGLE TO " << atom0 << " AND " << atom1 << " AND " << atom2 << "\n";
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

            // Variable number of parameters
            for (int k = 0; k < dihedral_params.length(); k = k + 3)
            {
                double v = dihedral_params[ k ];
                int periodicity = dihedral_params[ k + 1 ];
                double phase = dihedral_params[ k + 2 ];

                bondTorsion_openmm->addTorsion(idx0, idx1, idx2, idx3, periodicity, phase, v * OpenMM::KJPerKcal);

                /*cout << "Dihedral between atom global index " << idx0 << " and " << idx1 << " and " << idx2 << " and " << idx3<<"\n";
                cout << "Amplitude_dih = " << v << " periodicity " << periodicity << " phase " << phase<<"\n";*/
                //cout << "Dihedral local" << dihedral_ff.toString() << " v " << v << " periodicity " << periodicity << " phase " << phase;
                //cout << "\n";
            }
        }// end of dihedrals

        //Improper Dihedrals
        QList<ImproperID> impropers_ff = amber_params.getAllImpropers();
        QVector<ImproperID> impropers = impropers_ff.toVector();

        for (int j = 0; j < impropers_ff.length(); j++)
        {
            ImproperID improper_ff = impropers_ff[j];
            QList<double> improper_params = amber_params.getParams(improper_ff);
            // Variable number of parameters
            int idx0 = impropers[j].atom0().asA<AtomIdx>().value() + num_atoms_till_i;
            int idx1 = impropers[j].atom1().asA<AtomIdx>().value() + num_atoms_till_i;
            int idx2 = impropers[j].atom2().asA<AtomIdx>().value() + num_atoms_till_i;
            int idx3 = impropers[j].atom3().asA<AtomIdx>().value() + num_atoms_till_i;

            for (int k = 0; k < improper_params.length(); k = k + 3)
            {
                double v = improper_params[ k ];
                int periodicity = improper_params[ k + 1 ];
                double phase = improper_params[ k + 2 ];

                bondTorsion_openmm->addTorsion(idx0, idx1, idx2, idx3, periodicity, phase, v * OpenMM::KJPerKcal);
                /*cout << "Improper Dihedral between atom global index " << idx0 << " and " << idx1 << " and " << idx2 << " and " << idx3<<"\n";
                cout << "Amplitude_impr = " << v << " periodicity " << periodicity << " phase " << phase <<"\n";*/
                //cout << "\n";
            }
        }//end of impropers

        // Variable 1,4 scaling factors
        QList<BondID> pairs14_ff = amber_params.getAll14Pairs();
        QVector<BondID> pairs14 = pairs14_ff.toVector();
        //
        for (int j = 0; j < pairs14_ff.length(); j++)
        {
            BondID pair14_ff = pairs14_ff[j];
            QList<double> pair14_params = amber_params.get14PairParams(pair14_ff);
            double cscl = pair14_params[0];
            double ljscl = pair14_params[1];


            // Add to custom pairs if scale factor differs from default
            if (abs(cscl - Coulomb14Scale) > 0.0001 or abs(ljscl - LennardJones14Scale) > 0.0001)
            {
                int idx0 = pair14_ff.atom0().asA<AtomIdx>().value() + num_atoms_till_i;
                int idx1 = pair14_ff.atom1().asA<AtomIdx>().value() + num_atoms_till_i;
                QPair<int, int> indices_pair(idx0, idx1);
                QPair<double, double> scl_pair(cscl, ljscl);
                QPair< QPair<int, int>, QPair<double, double> > custom14pair(indices_pair, scl_pair);
                custom14pairs.append(custom14pair);
            }
        }// end of variable 1,4 scaling factors

        num_atoms_till_i = num_atoms_till_i + num_atoms_molecule;
    }// end of loop over molecules

    //Exclude the 1-2, 1-3 bonded atoms from nonbonded forces, and scale down 1-4 bonded atoms
    // using default scaling factors for the 1,4 interactions
    nonbond_openmm->createExceptionsFromBonds(bondPairs, Coulomb14Scale, LennardJones14Scale);

    // Now manually override exceptions for 1,4 pairs that do not have a standard 1,4 scale factor
    for (int i = 0; i < custom14pairs.length(); i++)
    {
        QPair< QPair<int, int>, QPair<double, double> > custom14pair = custom14pairs[ i ];
        int p1 = custom14pair.first.first;
        int p2 = custom14pair.first.second;
        double cscl = custom14pair.second.first;
        double ljscl = custom14pair.second.second;
        // get the particle parameters
        double q1;
        double sig1;
        double eps1;
        double q2;
        double sig2;
        double eps2;
        nonbond_openmm->getParticleParameters(p1, q1, sig1, eps1);
        nonbond_openmm->getParticleParameters(p2, q2, sig2, eps2);
        // now find chargeprod sigma epsilon
        double chargeprod = cscl * q1*q2;
        double sigma = 0.5 * (sig1 + sig2);
        double epsilon = ljscl * sqrt(eps1 * eps2);
        nonbond_openmm->addException(p1, p2, chargeprod, sigma, epsilon, true);
        //if (Debug)
        //    qDebug() << "manual exception " << p1 << " " << p2 << " cscl " << cscl << " ljscl " << ljscl << "\n";
    }

    if (CMMremoval_frequency > 0)
    {
        OpenMM::CMMotionRemover * cmmotionremover = new OpenMM::CMMotionRemover(CMMremoval_frequency);

        system_openmm->addForce(cmmotionremover);

        //if (Debug)
        //    qDebug() << "\n\nWill remove Center of Mass motion every " << CMMremoval_frequency << " steps\n\n";
    }

    this->openmm_system = system_openmm;

    this->isSystemInitialised = true;
    this->is_periodic = false;
}

void OpenMMMDIntegrator::createContext(IntegratorWorkspace &workspace,
                                       SireUnits::Dimension::Time timestep)
{
    bool Debug = false;

    if (Debug)
        qDebug() << "In OpenMMMDIntegrator::integrate()\n\n";

    // Check that the openmm system has been initialised
    // !! Should check that the workspace is compatible with molgroup
    if (not this->isSystemInitialised)
    {
        qDebug() << "Not initialised ! ";
        throw SireError::program_bug(QObject::tr(
                                                 "OpenMMMDintegrator should have been initialised before calling integrate."), CODELOC);
    }

    OpenMM::System *system_openmm = openmm_system;

    int nats = system_openmm->getNumParticles();

    if (Debug)
        qDebug() << " openmm nats " << nats;


    const double converted_Temperature = convertTo(Temperature.value(), kelvin);
    const double converted_friction = convertTo(friction.value(), picosecond);
    const double dt = convertTo(timestep.value(), picosecond);

    if (!isContextInitialised || (isContextInitialised && reinetialise_context))
    {
        OpenMM::Integrator * integrator_openmm = NULL;

        if (Integrator_type == "leapfrogverlet")
            integrator_openmm = new OpenMM::VerletIntegrator(dt); //dt in picosecond
        else if (Integrator_type == "variableleapfrogverlet")
            integrator_openmm = new OpenMM::VariableVerletIntegrator(integration_tol); //integration tolerance error unitless
        else if (Integrator_type == "langevin")
            integrator_openmm = new OpenMM::LangevinIntegrator(converted_Temperature, converted_friction, dt);
        else if (Integrator_type == "langevinmiddle")
            integrator_openmm = new OpenMM::LangevinMiddleIntegrator(converted_Temperature, converted_friction, dt);
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

    //OpenMM vector coordinate
    std::vector<OpenMM::Vec3> positions_openmm(nats);
    //OpenMM vector momenta
    std::vector<OpenMM::Vec3> velocities_openmm(nats);

    AtomicVelocityWorkspace &ws = workspace.asA<AtomicVelocityWorkspace>();

    is_periodic = false;

    if (CutoffType == "cutoffperiodic" || CutoffType == "ewald" || CutoffType == "pme")
    {
        const System & ptr_sys = ws.system();
        const PropertyName &space_property = PropertyName("space");

        // PeriodicBox.
        if (ptr_sys.property(space_property).isA<PeriodicBox>())
        {
            const PeriodicBox &space = ptr_sys.property(space_property).asA<PeriodicBox>();

            const double Box_x_Edge_Length = space.dimensions()[0] * OpenMM::NmPerAngstrom; //units in nm
            const double Box_y_Edge_Length = space.dimensions()[1] * OpenMM::NmPerAngstrom; //units in nm
            const double Box_z_Edge_Length = space.dimensions()[2] * OpenMM::NmPerAngstrom; //units in nm

            if (Debug)
                qDebug() << "\nBOX SIZE [A] = (" << space.dimensions()[0] << " , " << space.dimensions()[1] << " ,  " << space.dimensions()[2] << ")\n\n";

            //Set Periodic Box Condition

            openmm_context->setPeriodicBoxVectors(OpenMM::Vec3(Box_x_Edge_Length, 0, 0),
                                                  OpenMM::Vec3(0, Box_y_Edge_Length, 0),
                                                  OpenMM::Vec3(0, 0, Box_z_Edge_Length));
        }
        // TriclinicBox.
        else if (ptr_sys.property(space_property).isA<TriclinicBox>())
        {
            const TriclinicBox &space = ptr_sys.property(space_property).asA<TriclinicBox>();

            // Get the three triclinic box vectors.
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

            openmm_context->setPeriodicBoxVectors(OpenMM::Vec3(xx, xy, xz),
                                                  OpenMM::Vec3(yx, yy, yz),
                                                  OpenMM::Vec3(zx, zy, zz));
        }

        is_periodic = true;
    }


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

            positions_openmm[system_index] = OpenMM::Vec3(c[j].x() * (OpenMM::NmPerAngstrom),
                                                          c[j].y() * (OpenMM::NmPerAngstrom), c[j].z() * (OpenMM::NmPerAngstrom));

            if (m[j] == 0.0)
                qDebug() << "\nWARNING - THE MASS OF PARTICLE " << system_index << " is ZERO\n";

            if (m[j] > SireMaths::small)
            {
                velocities_openmm[system_index] = OpenMM::Vec3(p[j].x() / m[j] *
                                                               (OpenMM::NmPerAngstrom) * PsPerAKMA, p[j].y() / m[j] *
                                                               (OpenMM::NmPerAngstrom) * PsPerAKMA, p[j].z() / m[j] *
                                                               (OpenMM::NmPerAngstrom) * PsPerAKMA);
            }
            else
            {
                velocities_openmm[system_index] = OpenMM::Vec3(0.0, 0.0, 0.0);
            }

            if (Debug)
            {
                qDebug() << "Particle num = " << system_index;
                qDebug() << "Particle mass = " << m[j];
                qDebug() << "X = " << positions_openmm[system_index][0] * OpenMM::AngstromsPerNm
                    << " A" <<
                    " Y = " << positions_openmm[system_index][1] * OpenMM::AngstromsPerNm
                    << " A" <<
                    " Z = " << positions_openmm[system_index][2] * OpenMM::AngstromsPerNm
                    << " A";
                qDebug() << "Vx = " << velocities_openmm[system_index][0] << " Vy = "
                    << velocities_openmm[system_index][1] << " Vz = "
                    << velocities_openmm[system_index][2] << "\n";
            }
            system_index++;
        }
    }
    //openmmKineticEnergy = state_openmm.getKineticEnergy()* OpenMM::KcalPerKJ*kcal_per_mol;

    if (system_index != nats)
    {
        if (Debug)
            qDebug() << " system_index " << system_index << " nats " << nats;
        throw SireError::program_bug(QObject::tr("The number of atoms in the openmm system does not match the number of atoms in the sire workspace"), CODELOC);
    }

    openmm_context->setPositions(positions_openmm);
    openmm_context->setVelocities(velocities_openmm);


}

void OpenMMMDIntegrator::destroyContext()
{
    if (this->isContextInitialised)
    {
        delete openmm_context;
        openmm_context = 0;
        this->isContextInitialised = false;
    }
}

MolarEnergy OpenMMMDIntegrator::getPotentialEnergy(const System &system)
{
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
 * <Returns the kinetic energy of the OpenMM system>
 * minimizeEnergy will find the nearest local potential energy minimum,
 * given the current Sire::System. It calls the
 * LocalEnergyMinimizer :: minimize() function of OpenMM.
 * @param system                Sire System including molegroup, forcefield
 *                              positions etc
 * @return                      Kinetic energy computed with OpenMM.
 */
MolarEnergy OpenMMMDIntegrator::getKineticEnergy()
{
    //We need to compute the OpenMM kinetic energy because of the Verlet half
    //step algorithm. Sire kinetic energies will not be the same as the OpenMM
    //ones.
    if (!isSystemInitialised)
    {
        throw SireError::program_bug(QObject::tr(
                                                 "System was not initialised! Initialise first before requesting energies!"), CODELOC);
    }
    else
    {
        return openmmKineticEnergy;
    }
}


/**
 * <Runs an energy Minimisation on the current system.>
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
System OpenMMMDIntegrator::minimiseEnergy(System &system, double tolerance = 1, int max_iteration = 100)
{
    bool Debug = false;
    const MoleculeGroup moleculegroup = this->molgroup.read();
    IntegratorWorkspacePtr workspace = this->createWorkspace(moleculegroup);
    if (system.nMolecules() != moleculegroup.nMolecules())
    {
        std::cerr << "Number of molecules in do not agree!";
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
 * annealLambda will equilibrate the system to the current alchemical lambda
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

System OpenMMMDIntegrator::equilibrateSystem(System &system,
                                      SireUnits::Dimension::Time equib_time_step,
                                      int equib_steps)
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
    if (Debug)
    {
        PeriodicBox sp = system.property("space").asA<PeriodicBox>();
        cout << "Box dimensions are: "<< sp.dimensions()[0]<< " "<<
            sp.dimensions()[1]<<" " << sp.dimensions()[2]<<endl;
    }
    workspace.edit().setSystem(system);
    createContext(workspace.edit(), equib_time_step);
    (openmm_context->getIntegrator()).step(equib_steps);

    int infoMask = OpenMM::State::Positions;
    infoMask = infoMask + OpenMM::State::Velocities+OpenMM::State::Energy;
    OpenMM::State state_openmm = openmm_context->getState(infoMask);
    std::vector<OpenMM::Vec3> positions_openmm = state_openmm.getPositions();
    std::vector<OpenMM::Vec3> velocities_openmm = state_openmm.getVelocities();

    openmmKineticEnergy = state_openmm.getKineticEnergy()* OpenMM::KcalPerKJ*kcal_per_mol;

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
            if(Debug)
            {
                qDebug() <<"Sire momenta before update "<<sire_momenta[j].toString();
            }
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
    // update periodic box
    if(is_periodic)
    {
        // dummy buffered dimensions vector, maybe there is better solution
        //to this than just passing an empty vector
        QVector<QVector<Vector>> dimensions;
        updateBoxDimensions(state_openmm, dimensions, Debug, ws);
    }

    this->destroyContext();
    const System & ptr_sys = ws.system();
    return ptr_sys;
}


void OpenMMMDIntegrator::integrate(IntegratorWorkspace &workspace, const Symbol &nrg_component,
                                   SireUnits::Dimension::Time timestep,
                                   int nmoves, bool record_stats)
{
    bool Debug = false;

    createContext(workspace, timestep);

    const int nats = openmm_system->getNumParticles();

    AtomicVelocityWorkspace &ws = workspace.asA<AtomicVelocityWorkspace>();

    const double AKMAPerPs = 0.04888821;

    const int nmols = ws.nMolecules();

    QVector< std::vector<OpenMM::Vec3> > buffered_positions;
    QVector<QVector<Vector>> buffered_dimensions;

    OpenMM::State state_openmm;

    OpenMM::Vec3 a;
    OpenMM::Vec3 b;
    OpenMM::Vec3 c;

    int infoMask = 0;
    infoMask = OpenMM::State::Positions;
    infoMask = infoMask + OpenMM::State::Velocities;
    infoMask = infoMask + OpenMM::State::Energy;
    state_openmm = openmm_context->getState(infoMask);
    double EnergyInKcal = (state_openmm.getPotentialEnergy() + state_openmm.getKineticEnergy())* OpenMM::KcalPerKJ;

    const double dt = convertTo(timestep.value(), picosecond);

    if (Debug)
    {
        qDebug() << " Total energy from Openmm is "<<EnergyInKcal;
    }
      if (Debug)
    {
        qDebug() <<"\n";
        qDebug() <<"-------------info before ---------------";
        qDebug() << "Kinetic energy Sire is "<<ws.kineticEnergy();
        qDebug() << "Kinetic energy OpenMM is "<<state_openmm.getKineticEnergy()*OpenMM::KcalPerKJ;
        //qDebug() << "Total energy from Sire after step is "<<ws.energy()+ws.kineticEnergy();
        qDebug() <<"---------------------------------";
        qDebug() <<"\n";
    }

    // Coordinates are buffered every coord_freq steps
    int coord_freq = buffer_frequency;

    int nframes;
    int MAXFRAMES = 1000;


    if (Debug)
        qDebug() << " nmoves " << nmoves << " coord_freq " << coord_freq;

    // Limit excessive internal buffering
    if (coord_freq > 0)
    {
        nframes = (nmoves / coord_freq);

        if (nframes > MAXFRAMES)
        {
            throw SireError::program_bug(QObject::tr(
                                                     "You are requesting to buffer %1 frames, which is above the hardcoded limit of %2.").arg(nframes, MAXFRAMES), CODELOC);
        }
    }
    else
    {
        nframes = 0;
    }

    //OpenMM vector coordinate
    std::vector<OpenMM::Vec3> positions_openmm(nats);
    //OpenMM vector momenta
    std::vector<OpenMM::Vec3> velocities_openmm(nats);

    state_openmm = openmm_context->getState(infoMask);
    double mypotential_energy = state_openmm.getPotentialEnergy();

    if (Debug)
    {
        qDebug() << " pot nrg bef dyn " << mypotential_energy*0.239006;
        qDebug() << "pot nrg sire bef dyn "<<workspace.nonConstsystem().energy();
    }

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
        //OpenMM::integrator = (openmm_context->getIntegrator());

        (openmm_context->getIntegrator()).step(new_nmoves);


        nmoves = (nmoves - new_nmoves);

        if (coord_freq > 0)
            nframes = nmoves / coord_freq;

    }
    if (coord_freq > 0)
    {/** Break nmoves in several steps to buffer coordinates*/
        for (int i = 0; i < nmoves; i = i + coord_freq)
        {

            if (Debug)
                qDebug() << " about to step ";

            (openmm_context->getIntegrator()).step(coord_freq);

            if (Debug)
                qDebug() << " i now " << i;

            state_openmm = openmm_context->getState(infoMask);
            positions_openmm = state_openmm.getPositions();
            buffered_positions.append(positions_openmm);

            state_openmm.getPeriodicBoxVectors(a, b, c);

            Vector v0 = Vector(a[0] * OpenMM::AngstromsPerNm, a[1] * OpenMM::AngstromsPerNm, a[2] * OpenMM::AngstromsPerNm);
            Vector v1 = Vector(b[0] * OpenMM::AngstromsPerNm, b[1] * OpenMM::AngstromsPerNm, b[2] * OpenMM::AngstromsPerNm);
            Vector v2 = Vector(c[0] * OpenMM::AngstromsPerNm, c[1] * OpenMM::AngstromsPerNm, c[2] * OpenMM::AngstromsPerNm);

            QVector<Vector> lattice_vectors{v0, v1, v2};

            buffered_dimensions.append(lattice_vectors);
        }
    }
    else
    {/** No buffering*/
        (openmm_context->getIntegrator()).step(nmoves);
    }

    if (time_skip != 0)
    {
        timeskip = SireUnits::Dimension::Time(0.0);
    }
    bool IsFiniteNumber = true;

    state_openmm = openmm_context->getState(infoMask);
    positions_openmm = state_openmm.getPositions();
    velocities_openmm = state_openmm.getVelocities();
    double potential_energy = state_openmm.getPotentialEnergy();


    IsFiniteNumber = (potential_energy <= DBL_MAX && potential_energy >= -DBL_MAX);

    if (!IsFiniteNumber)
    {
        qDebug() << "NaN or Inf has been generated along the simulation";
        exit(-1);
    }


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
    double Ekin_openmm=0;

    for (int i = 0; i < nmols; i++)
    {
        Vector *sire_coords = ws.coordsArray(i);
        Vector *sire_momenta = ws.momentaArray(i);
        const double *m = ws.massArray(i);

        for (int j = 0; j < ws.nAtoms(i); j++)
        {
            if(Debug)
            {
                qDebug() <<"Sire momenta before update "<<sire_momenta[j].toString();
            }

            sire_coords[j] = Vector(positions_openmm[j + k][0] * (OpenMM::AngstromsPerNm),
                                    positions_openmm[j + k][1] * (OpenMM::AngstromsPerNm),
                                    positions_openmm[j + k][2] * (OpenMM::AngstromsPerNm));


            for (int l = 0; l < nframes; l++)
            {
                //qDebug() << " i " << i << " j " << j << " k " << k << " l " << l;

                Vector buffered_atcoord = Vector(buffered_positions[l][j + k][0] * (OpenMM::AngstromsPerNm),
                                                 buffered_positions[l][j + k][1] * (OpenMM::AngstromsPerNm),
                                                 buffered_positions[l][j + k][2] * (OpenMM::AngstromsPerNm));
                buffered_workspace[l][i][j] = buffered_atcoord;
            }

            sire_momenta[j] = Vector(
                                     velocities_openmm[j + k][0] * m[j] * (OpenMM::AngstromsPerNm) * AKMAPerPs,
                                     velocities_openmm[j + k][1] * m[j] * (OpenMM::AngstromsPerNm) * AKMAPerPs,
                                     velocities_openmm[j + k][2] * m[j] * (OpenMM::AngstromsPerNm) * AKMAPerPs);
            if (Debug)
            {
                double vsqrt = velocities_openmm[j + k][0]*velocities_openmm[j + k][0];
                vsqrt += velocities_openmm[j + k][1]*velocities_openmm[j + k][1];
                vsqrt += velocities_openmm[j + k][2]*velocities_openmm[j + k][2];
                Ekin_openmm+=0.5*vsqrt*m[j];
                OpenMM::System *system_openmm = openmm_system;
                double mass = system_openmm->getParticleMass(j);
                qDebug()<<"particle mass "<<mass;
                qDebug()<<"j is : "<<j<<" k is ";
;
            }
        }

        k = k + ws.nAtoms(i);
    }
    //if (Debug)
    //    qDebug()<< "Kinetic energy calculated by hand opennmm "<<Ekin_openmm*OpenMM::KcalPerKJ;
    openmmKineticEnergy = state_openmm.getKineticEnergy()*OpenMM::KcalPerKJ*kcal_per_mol;

    System & ptr_sys = ws.nonConstsystem();
    ptr_sys.mustNowRecalculateFromScratch();

    if (nframes <= 0)
        ws.commitCoordinatesAndVelocities();
    else
        ws.commitBufferedCoordinatesAndVelocities(buffered_workspace);

    /** Now the box dimensions (if the simulation used a periodic space) */
    if (is_periodic)
    {
        updateBoxDimensions(state_openmm, buffered_dimensions, true, ws);
    }

    if (ws.system().contains(molgroup.read().number()))
    {
        molgroup = ws.system()[molgroup.read().number()];
    }
    else
    {
        molgroup.edit().update(ws.system().molecules());
    }

    EnergyInKcal = (state_openmm.getPotentialEnergy() + state_openmm.getKineticEnergy())* OpenMM::KcalPerKJ;


    if (Debug)
    {
        qDebug() <<"-------------OpenMM infor ---------------";
        qDebug() << "Potential energy Openmm is "<<state_openmm.getPotentialEnergy()*OpenMM::KcalPerKJ;
        qDebug() << "Kinteic energy Openmm is "<<state_openmm.getKineticEnergy()*OpenMM::KcalPerKJ;
        qDebug() << "Total energy from Openmm  after step is "<<EnergyInKcal;
        qDebug() <<"-------------Sire info ---------------";
        qDebug() << "Potential energy Sire is"<<ptr_sys.energy();
        qDebug() << "Kinteic energy Sire is "<<ws.kineticEnergy();
        qDebug()<< "Kinetic energy calculated by hand opennmm "<<Ekin_openmm*OpenMM::KcalPerKJ;

        qDebug() << "Total energy from Sire after step is "<<ptr_sys.energy()+ws.kineticEnergy();
        qDebug() <<"---------------------------------";
    }

    /** Clear all buffers */
    buffered_workspace.clear();
    buffered_dimensions.clear();

    return;
}

void OpenMMMDIntegrator::updateBoxDimensions(OpenMM::State &state_openmm,
                                             QVector<QVector<Vector>> &buffered_dimensions,
                                             bool Debug, AtomicVelocityWorkspace &ws)
{
    Debug = false;
    OpenMM::Vec3 a;
    OpenMM::Vec3 b;
    OpenMM::Vec3 c;

    System & ptr_sys = ws.nonConstsystem();

    // TriclinicBox.
    if (ptr_sys.property("space").isA<TriclinicBox>())
    {
        state_openmm.getPeriodicBoxVectors(a, b, c);
        Vector v0 = Vector(a[0] * OpenMM::AngstromsPerNm, a[1] * OpenMM::AngstromsPerNm, a[2] * OpenMM::AngstromsPerNm);
        Vector v1 = Vector(b[0] * OpenMM::AngstromsPerNm, b[1] * OpenMM::AngstromsPerNm, b[2] * OpenMM::AngstromsPerNm);
        Vector v2 = Vector(c[0] * OpenMM::AngstromsPerNm, c[1] * OpenMM::AngstromsPerNm, c[2] * OpenMM::AngstromsPerNm);

        System & ptr_sys = ws.nonConstsystem();
        TriclinicBox sp(v0, v1, v2);

        const QString string = "space";
        ptr_sys.setProperty(string, sp);

        /** Buffer dimensions if necessary */
        for (int k = 0; k < buffered_dimensions.size(); k++)
        {
            const QString buffered_space = "buffered_space_" + QString::number(k);
            TriclinicBox buff_space = TriclinicBox(buffered_dimensions[k][0],
                                                   buffered_dimensions[k][1],
                                                   buffered_dimensions[k][2]);
            ptr_sys.setProperty(buffered_space, buff_space);
        }
    }
    // PeriodicBox.
    else
    {
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
            Vector dims(buffered_dimensions[k][0].x(),
                        buffered_dimensions[k][1].y(),
                        buffered_dimensions[k][2].z());
            PeriodicBox buff_space = PeriodicBox(dims);
            ptr_sys.setProperty(buffered_space, buff_space);
        }
    }
}

QString OpenMMMDIntegrator::getIntegrator(void)
{
    return Integrator_type;
}

void OpenMMMDIntegrator::setIntegrator(QString intgrator)
{
    Integrator_type = intgrator;
}

SireUnits::Dimension::Time OpenMMMDIntegrator::getFriction(void)
{
    return friction;
}

void OpenMMMDIntegrator::setFriction(SireUnits::Dimension::Time thefriction)
{
    friction = thefriction;
}

/** Get the cufott type: nocutoff, cutoffnonperiodic, cutoffperiodic */
QString OpenMMMDIntegrator::getCutoffType(void)
{
    return CutoffType;
}

/** Set the cufott type: nocutoff, cutoffnonperiodic, cutoffperiodic */
void OpenMMMDIntegrator::setCutoffType(QString cutoff_type)
{
    CutoffType = cutoff_type;
}

/** Get the cutoff distance in A */
SireUnits::Dimension::Length OpenMMMDIntegrator::getCutoffDistance(void)
{
    return cutoff_distance;
}

/** Set the cutoff distance in A */
void OpenMMMDIntegrator::setCutoffDistance(SireUnits::Dimension::Length distance)
{
    cutoff_distance = distance;
}

/** Get the dielectric constant */
double OpenMMMDIntegrator::getFieldDielectric(void)
{
    return field_dielectric;
}

/** Set the dielectric constant */
void OpenMMMDIntegrator::setFieldDielectric(double dielectric)
{
    field_dielectric = dielectric;
}

/** Set the Ewald or PME tolerance */
double OpenMMMDIntegrator::getToleranceEwaldPME(void)
{
    return tolerance_ewald_pme;
}

void OpenMMMDIntegrator::setToleranceEwaldPME(double toll)
{
    tolerance_ewald_pme = toll;
}

/** Set Andersen thermostat */

void OpenMMMDIntegrator::setAndersen(bool andersen)
{
    Andersen_flag = andersen;
}

/** Get Andersen thermostat status on/off */
bool OpenMMMDIntegrator::getAndersen(void)
{
    return Andersen_flag;
}

/** Get the Andersen Thermostat frequency collision */
double OpenMMMDIntegrator::getAndersenFrequency(void)
{
    return Andersen_frequency;
}

/** Set the Andersen Thermostat frequency collision */
void OpenMMMDIntegrator::setAndersenFrequency(double freq)
{
    Andersen_frequency = freq;
}

/** Get the bath Temperature */
SireUnits::Dimension::Temperature OpenMMMDIntegrator::getTemperature(void)
{
    return Temperature;
}

/** Set the Temperature */
void OpenMMMDIntegrator::setTemperature(SireUnits::Dimension::Temperature temperature)
{
    Temperature = temperature;
}

/** Set Monte Carlo Barostat on/off */
void OpenMMMDIntegrator::setMCBarostat(bool MCBarostat)
{
    MCBarostat_flag = MCBarostat;
}

/** Get Andersen thermostat status on/off */
bool OpenMMMDIntegrator::getMCBarostat(void)
{
    return MCBarostat_flag;
}

/** Get the Monte Carlo Barostat frequency in time speps */
int OpenMMMDIntegrator::getMCBarostatFrequency(void)
{
    return MCBarostat_frequency;
}

/** Set the Monte Carlo Barostat frequency in time speps */
void OpenMMMDIntegrator::setMCBarostatFrequency(int freq)
{
    MCBarostat_frequency = freq;
}

/** Get the Presure */
SireUnits::Dimension::Pressure OpenMMMDIntegrator::getPressure(void)
{
    return Pressure;
}

/** Set the Pressure */
void OpenMMMDIntegrator::setPressure(SireUnits::Dimension::Pressure pressure)
{
    Pressure = pressure;
}

/** Get the Constraint type: none, hbonds, allbonds, hangles */
QString OpenMMMDIntegrator::getConstraintType(void)
{
    return ConstraintType;
}

/** Set the Constraint type: none, hbonds, allbonds, hangles */
void OpenMMMDIntegrator::setConstraintType(QString constrain)
{
    ConstraintType = constrain;
}

/** Get the OpenMMMD Platform: CUDA, OpenCL, CPU */
QString OpenMMMDIntegrator::getPlatform(void)
{
    return platform_type;
}

/** Set the OpenMM Platform: CUDA, OpenCL, CPU */
void OpenMMMDIntegrator::setPlatform(QString platform)
{
    platform_type = platform;
}

/** Get the OpenMMMD Platform: CUDA, OpenCL, CPU */
QString OpenMMMDIntegrator::getDeviceIndex(void)
{
    return device_index;
}

/** Set the OpenMM Platform: CUDA, OpenCL, CPU */
void OpenMMMDIntegrator::setDeviceIndex(QString deviceidx)
{
    device_index = deviceidx;
}

/** Set the OpenMM Precision */
void OpenMMMDIntegrator::setPrecision(QString prec)
{

    precision = prec;

}

/** Get the OpenMMMD Precision */
QString OpenMMMDIntegrator::getPrecision(void)
{

    return precision;

}

/** Get the Restaint mode*/
bool OpenMMMDIntegrator::getLJDispersion(void)
{

    return LJ_dispersion;

}

/** Set the Retraint mode */
void OpenMMMDIntegrator::setLJDispersion(bool LJ_disp)
{

    LJ_dispersion = LJ_disp;
}

/** Get the Restaint mode*/
bool OpenMMMDIntegrator::getRestraint(void)
{

    return Restraint_flag;

}

/** Set the Retraint mode */
void OpenMMMDIntegrator::setRestraint(bool Restraint)
{

    Restraint_flag = Restraint;

}

/** Get the Center of Mass motion removal frequency */
int OpenMMMDIntegrator::getCMMremovalFrequency(void)
{

    return CMMremoval_frequency;

}

/** Set the Center of Mass motion removal frequency */
void OpenMMMDIntegrator::setCMMremovalFrequency(int frequency)
{

    CMMremoval_frequency = frequency;

}

/** Get the frequency of buffering coordinates */
int OpenMMMDIntegrator::getBufferFrequency()
{

    return buffer_frequency;

}

/** Set the Center of Mass motion removal frequency */
void OpenMMMDIntegrator::setBufferFrequency(int frequency)
{

    buffer_frequency = frequency;

}

/** Set the flag to reinitialise the context*/
void OpenMMMDIntegrator::setReinitialiseContext(bool reinitialise)
{

    reinetialise_context = reinitialise;
}

/** Get the integration tolerance */
double OpenMMMDIntegrator::getIntegrationTolerance(void)
{

    return integration_tol;

}

/** Set the integration tolerance*/
void OpenMMMDIntegrator::setIntegrationTolerance(double tolerance)
{
    integration_tol = tolerance;
}

/** Get total time to skip*/
SireUnits::Dimension::Time OpenMMMDIntegrator::getTimetoSkip(void)
{
    return timeskip;
}

/** Get total time to skip*/
void OpenMMMDIntegrator::setTimetoSkip(SireUnits::Dimension::Time skip)
{
    timeskip = skip;
}

/** Create an empty workspace */
IntegratorWorkspacePtr OpenMMMDIntegrator::createWorkspace(const PropertyMap &map) const
{

    return IntegratorWorkspacePtr(new AtomicVelocityWorkspace(map));

}

/** Return the ensemble of this integrator */
Ensemble OpenMMMDIntegrator::ensemble() const
{

    return Ensemble::NVE();

}

/** Return whether or not this integrator is time-reversible */
bool OpenMMMDIntegrator::isTimeReversible() const
{

    return true;

}

/** Create a workspace for this integrator for the molecule group 'molgroup' */
IntegratorWorkspacePtr OpenMMMDIntegrator::createWorkspace(const MoleculeGroup &molgroup, const PropertyMap &map) const
{

    return IntegratorWorkspacePtr(new AtomicVelocityWorkspace(molgroup, map));

}

const char* OpenMMMDIntegrator::typeName()
{

    return QMetaType::typeName(qMetaTypeId<OpenMMMDIntegrator>());

}
