/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2010  Christopher Woods
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

#include "hybridmc.h"

#include "SireMol/atomelements.h"
#include "SireMol/atommasses.h"
#include "SireMol/atomvelocities.h"

#include "SireUnits/units.h"
#include "SireUnits/temperature.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireMove;
using namespace SireMol;
using namespace SireBase;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

//////////
////////// Implementation of HMCVelGen
//////////

static const RegisterMetaType<HMCVelGen> r_hmcvelgen( MAGIC_ONLY,
                                                      HMCVelGen::typeName() );

QDataStream &operator<<(QDataStream &ds, const HMCVelGen &hmcvelgen)
{
    writeHeader(ds, r_hmcvelgen, 1);
    
    ds << hmcvelgen.rangen << hmcvelgen.temp 
       << static_cast<const VelocityGenerator&>(hmcvelgen);
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds, HMCVelGen &hmcvelgen)
{
    VersionID v = readHeader(ds, r_hmcvelgen);
    
    if (v == 1)
    {
        ds >> hmcvelgen.rangen >> hmcvelgen.temp 
           >> static_cast<VelocityGenerator&>(hmcvelgen);
    }
    else
        throw version_error(v, "1", r_hmcvelgen, CODELOC);

    return ds;
}

/** Constructor */
HMCVelGen::HMCVelGen() : VelocityGenerator(), temp(25*celsius)
{}

/** Copy constructor */
HMCVelGen::HMCVelGen(const HMCVelGen &other) 
          : VelocityGenerator(other), rangen(other.rangen), temp(other.temp)
{}

/** Destructor */
HMCVelGen::~HMCVelGen()
{}

/** Copy assignment operator */
HMCVelGen& HMCVelGen::operator=(const HMCVelGen &other)
{
    if (this != &other)
    {
        temp = other.temp;
        rangen = other.rangen;
        VelocityGenerator::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool HMCVelGen::operator==(const HMCVelGen &other) const
{
    return temp == other.temp and VelocityGenerator::operator==(other);
}

/** Comparison operator */
bool HMCVelGen::operator!=(const HMCVelGen &other) const
{
    return not HMCVelGen::operator==(other);
}

/** Set the random number generator used to generate the random 
    numbers needed by this generator */
void HMCVelGen::setGenerator(const RanGenerator &generator)
{
    rangen = generator;
}

/** Return the random number generator used by this generator */
const RanGenerator& HMCVelGen::rand() const
{
    return rangen;
}

/** Return the temperature of the velocities to be generated */
Temperature HMCVelGen::temperature() const
{
    return temp;
}

/** Set the temperature of the velocities to be generated */
void HMCVelGen::setTemperature(Temperature temperature)
{
    temp = temperature;
}

//////////
////////// Implementation of HMCGenerator
//////////

static const RegisterMetaType<HMCGenerator> r_hmcgen;

QDataStream &operator<<(QDataStream &ds, const HMCGenerator &hmcgen)
{
    writeHeader(ds, r_hmcgen, 1);
    
    ds << static_cast<const HMCVelGen&>(hmcgen);
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds, HMCGenerator &hmcgen)
{
    VersionID v = readHeader(ds, r_hmcgen);
    
    if (v == 1)
    {
        ds >> static_cast<HMCVelGen&>(hmcgen);
    }
    else
        throw version_error(v, "1", r_hmcgen, CODELOC);
        
    return ds;
}

/** Constructor */
HMCGenerator::HMCGenerator() : ConcreteProperty<HMCGenerator,HMCVelGen>()
{}

/** Copy constructor */
HMCGenerator::HMCGenerator(const HMCGenerator &other)
             : ConcreteProperty<HMCGenerator,HMCVelGen>(other)
{}

/** Destructor */
HMCGenerator::~HMCGenerator()
{}

const char* HMCGenerator::typeName()
{
    return QMetaType::typeName( qMetaTypeId<HMCGenerator>() );
}

/** Copy assignment operator */
HMCGenerator& HMCGenerator::operator=(const HMCGenerator &other)
{
    HMCVelGen::operator=(other);
    return *this;
}

/** Comparison operator */
bool HMCGenerator::operator==(const HMCGenerator &other) const
{
    return HMCVelGen::operator==(other);
}

/** Comparison operator */
bool HMCGenerator::operator!=(const HMCGenerator &other) const
{
    return HMCVelGen::operator!=(other);
}

/** Return the bias of the current state of the passed MD move */
double HMCGenerator::getBias(const MolecularDynamics &md)
{
    return std::exp( -md.kineticEnergy().value() / (k_boltz*temperature().to(kelvin)) );
}

/** Generate the velocities for the passed MD move, returning
    the biasing factor */
double HMCGenerator::generate(const System &system, MolecularDynamics &md)
{
    md.regenerateVelocities(system, *this);
    return HMCGenerator::getBias(md);
}

static QVector<MolarMass> getMassesFromElements(const MoleculeView &molview,
                                                const PropertyMap &map)
{
    const AtomElements &elements = molview.data().property(map["element"])
                                                 .asA<AtomElements>();
                                                 
    QVector<Element> elms;
    
    if (molview.selectedAll())
        elms = elements.toVector();
    else
        elms = elements.toVector(molview.selection());
        
    int sz = elms.count();
        
    QVector<MolarMass> masses(sz);
    
    const Element *elm_array = elms.constData();
    MolarMass *masses_array = masses.data();
    
    for (int i=0; i<sz; ++i)
    {
        masses_array[i] = elm_array[i].mass();
    }
    
    return masses;
}

static QVector<MolarMass> getMasses(const MoleculeView &molview, 
                                    const PropertyMap &map)
{
    PropertyName mass_property = map["mass"];
    
    if (not molview.data().hasProperty(mass_property))
        return ::getMassesFromElements(molview, map);
        
    const AtomMasses &masses = molview.data().property(mass_property)
                                             .asA<AtomMasses>();
                                             
    if (molview.selectedAll())
        return masses.toVector();
    else
        return masses.toVector(molview.selection());
}

/** Generate the velocities for the passed molecule view */
AtomVelocities HMCGenerator::generate(const MoleculeView &molview, 
                                      const PropertyMap &map) const
{
    const double kT = temperature().to(kelvin) * k_boltz;

    AtomVelocities molvels( molview.data().info() );
    
    QVector<Velocity3D> vels;
    
    if (molview.selectedAll())
        vels = molvels.toVector();
    else
        vels = molvels.toVector(molview.selection());
        
    QVector<MolarMass> masses = ::getMasses(molview, map);
    
    int sz = vels.count();
    BOOST_ASSERT( masses.count() == sz );
    
    Velocity3D *vels_array = vels.data();
    const MolarMass *masses_array = masses.constData();
    
    for (int i=0; i<sz; ++i)
    {
        //generate random momenta from the following gaussian
        //distribution (see documents/hybridmc.pdf)
        //
        //  prob(P) proportional exp{ -beta * P^2 / 2m }
        //
        // The normal distribution is proportional to;
        //
        //   exp{ -(x - mu)^2 / 2 sigma^2 }
        //
        // where mu is the mean, and sigma the variance
        //
        // To generate the desired momenta, we thus need
        // to generate normally distributed values, with mean 0
        // and variance sqrt( m / beta ) = sqrt( kT m ) 
        //
        // However, as we work with velocities (P = mV), we need
        //
        // prob(V) proportional exp{ -beta * m^2 V^2 / 2m }  = exp{ -beta * m V^2 / 2m }
        //
        // so need variance sqrt( kT / m )
        //
        // First generate a random vector from a normal distribution
        // with 0 mean and unit variance
        Vector norm_rand( rand().randNorm(0,1),
                          rand().randNorm(0,1),
                          rand().randNorm(0,1) );
        
        // the velocity is this, multiplied by sqrt( kT / m )
        if (masses_array[i].value() == 0)
            vels_array[i] = Velocity3D(0*miles_per_hour);
        else
            vels_array[i] = Velocity3D( std::sqrt(kT / masses_array[i].value()) 
                                                            * norm_rand );
    }

    if (molview.selectedAll())
        molvels.copyFrom(vels);
    else
        molvels.copyFrom(vels, molview.selection());
        
    return molvels;
}

//////////
////////// Implementation of HybridMC
//////////

static const RegisterMetaType<HybridMC> r_hmc;

QDataStream &operator<<(QDataStream &ds, const HybridMC &hmc)
{
    writeHeader(ds, r_hmc, 1);
    
    SharedDataStream sds(ds);
    
    sds << hmc.md << hmc.nsteps << static_cast<const MonteCarlo&>(hmc);
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds, HybridMC &hmc)
{
    VersionID v = readHeader(ds, r_hmc);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> hmc.md >> hmc.nsteps >> static_cast<MonteCarlo&>(hmc);
    }
    else
        throw version_error(v, "1", r_hmc, CODELOC);
    
    
    return ds;
}

/** Empty constructor */
HybridMC::HybridMC(const PropertyMap &map)
         : ConcreteProperty<HybridMC,MonteCarlo>(map),
           md(map), velgen(HMCGenerator()), nsteps(100)
{
    BOOST_ASSERT( md.ensemble() == Ensemble::NVE() );
    BOOST_ASSERT( md.integrator().isTimeReversible() );

    MonteCarlo::setEnsemble( Ensemble::NVT(25*celsius) );
}

/** Construct to perform hybrid MC moves on the passed molecule group */
HybridMC::HybridMC(const MoleculeGroup &molgroup, const PropertyMap &map)
         : ConcreteProperty<HybridMC,MonteCarlo>(map),
           md(molgroup, map), velgen(HMCGenerator()), nsteps(100)
{
    BOOST_ASSERT( md.ensemble() == Ensemble::NVE() );
    BOOST_ASSERT( md.integrator().isTimeReversible() );

    MonteCarlo::setEnsemble( Ensemble::NVT(25*celsius) );
}

void HybridMC::assertValidIntegrator(const Integrator &integrator) const
{
    if (integrator.ensemble() != Ensemble::NVE())
        throw SireError::incompatible_error( QObject::tr(
                "You can only use an integrator that samples from the "
                "microcanonical (NVE) ensemble with a Hybrid MC move. "
                "The integrator you've passed (%1) uses the ensemble %2. "
                "Try an integrator like VelocityVerlet().")
                    .arg(integrator.toString())
                    .arg(integrator.ensemble().toString()), CODELOC );
                    
    if (not integrator.isTimeReversible())
        throw SireError::incompatible_error( QObject::tr(
                "You can only use a time-reversible integrator with "
                "a Hybrid MC move. The integrator you've passed (%1) "
                "is not time-reversible. Type an integrator like "
                "VelocityVerlet().")
                    .arg(integrator.toString()), CODELOC );
}

/** Construct to perform hybrid MC moves on the passed molecule group
    using the passed integrator */
HybridMC::HybridMC(const MoleculeGroup &molgroup,
                   const Integrator &integrator,
                   const PropertyMap &map)
         : ConcreteProperty<HybridMC,MonteCarlo>(map),
           md(molgroup, integrator, map), velgen(HMCGenerator()), nsteps(100)
{
    assertValidIntegrator(integrator);
    MonteCarlo::setEnsemble( Ensemble::NVT(25*celsius) );
}

/** Construct to perform hybrid MC moves on the passed molecule group
    with a timestep of 'timestep' */
HybridMC::HybridMC(const MoleculeGroup &molgroup, Time timestep, const PropertyMap &map)
         : ConcreteProperty<HybridMC,MonteCarlo>(map),
           md(molgroup, timestep, map), velgen(HMCGenerator()), nsteps(100)
{
    BOOST_ASSERT( md.ensemble() == Ensemble::NVE() );
    BOOST_ASSERT( md.integrator().isTimeReversible() );

    MonteCarlo::setEnsemble( Ensemble::NVT(25*celsius) );
}

/** Construct to perform hybrid MC moves on the passed molecule group
    using the passed integrator and a timestep of 'timestep' */
HybridMC::HybridMC(const MoleculeGroup &molgroup, const Integrator &integrator,
                   Time timestep, const PropertyMap &map)
         : ConcreteProperty<HybridMC,MonteCarlo>(map),
           md(molgroup, integrator, timestep, map), velgen(HMCGenerator()), nsteps(100)
{
    assertValidIntegrator(integrator);
    MonteCarlo::setEnsemble( Ensemble::NVT(25*celsius) );
}

/** Construct to perform hybrid MC moves on the passed molecule group,
    running 'nmoves' MD moves per MC move */
HybridMC::HybridMC(const MoleculeGroup &molgroup, int n, const PropertyMap &map)
         : ConcreteProperty<HybridMC,MonteCarlo>(map),
           md(molgroup, map), velgen(HMCGenerator()), nsteps(0)
{
    if (n > 0)
        nsteps = n;

    BOOST_ASSERT( md.ensemble() == Ensemble::NVE() );
    BOOST_ASSERT( md.integrator().isTimeReversible() );

    MonteCarlo::setEnsemble( Ensemble::NVT(25*celsius) );
}

/** Construct to perform hybrid MC moves on the passed molecule group,
    using the passed integrator, running 'nmoves' MD moves per MC move */
HybridMC::HybridMC(const MoleculeGroup &molgroup, const Integrator &integrator, 
                   int n, const PropertyMap &map)
         : ConcreteProperty<HybridMC,MonteCarlo>(map),
           md(molgroup, integrator, map), velgen(HMCGenerator()), nsteps(0)
{
    if (n > 0)
        nsteps = n;

    assertValidIntegrator(integrator);
    MonteCarlo::setEnsemble( Ensemble::NVT(25*celsius) );
}

/** Construct to perform hybrid MC moves on the passed molecule group,
    with a timestep of 'timestep', running 'nmoves' MD moves per MC move */
HybridMC::HybridMC(const MoleculeGroup &molgroup, Time timestep, int n,
                   const PropertyMap &map)
         : ConcreteProperty<HybridMC,MonteCarlo>(map),
           md(molgroup, timestep, map), velgen(HMCGenerator()), nsteps(0)
{
    if (n > 0)
        nsteps = n;

    BOOST_ASSERT( md.ensemble() == Ensemble::NVE() );
    BOOST_ASSERT( md.integrator().isTimeReversible() );

    MonteCarlo::setEnsemble( Ensemble::NVT(25*celsius) );
}

/** Construct to perform hybrid MC moves on the passed molecule group,
    using the passed integrator, with a timestep of 'timestep', 
    running 'nmoves' MD moves per MC move */
HybridMC::HybridMC(const MoleculeGroup &molgroup, const Integrator &integrator,
                   Time timestep, int n, const PropertyMap &map)
         : ConcreteProperty<HybridMC,MonteCarlo>(map),
           md(molgroup, integrator, timestep, map), velgen(HMCGenerator()), nsteps(0)
{
    if (n > 0)
        nsteps = n;

    assertValidIntegrator(integrator);
    MonteCarlo::setEnsemble( Ensemble::NVT(25*celsius) );
}

/** Copy constructor */
HybridMC::HybridMC(const HybridMC &other)
         : ConcreteProperty<HybridMC,MonteCarlo>(other),
           md(other.md), velgen(other.velgen), nsteps(other.nsteps)
{}

/** Destructor */
HybridMC::~HybridMC()
{}

const char* HybridMC::typeName()
{
    return QMetaType::typeName( qMetaTypeId<HybridMC>() );
}

/** Copy assignment operator */
HybridMC& HybridMC::operator=(const HybridMC &other)
{
    if (this != &other)
    {
        md = other.md;
        nsteps = other.nsteps;
        velgen = other.velgen;
        MonteCarlo::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool HybridMC::operator==(const HybridMC &other) const
{
    return this == &other or
           (md == other.md and nsteps == other.nsteps and
            velgen == other.velgen and MonteCarlo::operator==(other));
}

/** Comparison operator */
bool HybridMC::operator!=(const HybridMC &other) const
{
    return not HybridMC::operator==(other);
}

QString HybridMC::toString() const
{
    return QObject::tr("HybridMC( MD = %1, nDynamicsSteps() = %2, "
                       "velocityGenerator() = %3 )")
                .arg(md.toString()).arg(nsteps)
                .arg(velgen.read().toString());
}

/** Return the molecule group that is moved by this move */
const MoleculeGroup& HybridMC::moleculeGroup() const
{
    return md.moleculeGroup();
}

/** Set the location of the coordinates that are affected by this move */
void HybridMC::setCoordinatesProperty(const PropertyName &coords_property)
{
    md.setCoordinatesProperty(coords_property);
    MonteCarlo::setCoordinatesProperty(coords_property);
}

/** Set the location of the space property used by this move */
void HybridMC::setSpaceProperty(const PropertyName &space_property)
{
    md.setSpaceProperty(space_property);
    MonteCarlo::setSpaceProperty(space_property);
}

/** Set the random number generator that is used by this move */
void HybridMC::setGenerator(const RanGenerator &rangenerator)
{
    md.setGenerator(rangenerator);
    MonteCarlo::setGenerator(rangenerator);
}

/** Set the velocity generator used to generate velocities that
    are compatible with the Hybrid MC move */
void HybridMC::setVelocityGenerator(const HMCVelGen &generator)
{
    velgen = generator;
    velgen.edit().asA<HMCVelGen>().setTemperature( this->temperature() );
}

/** Return the velocity generator used to generate the velocities
    that are used by the hybrid MC move */
const HMCVelGen& HybridMC::velocityGenerator() const
{
    return velgen.read().asA<HMCVelGen>();
}

/** Set the timestep of the dynamics part of the move */
void HybridMC::setTimeStep(Time timestep)
{
    md.setTimeStep(timestep);
}

/** Return the timestep of the dynamics part of the move */
Time HybridMC::timeStep() const
{
    return md.timeStep();
}

/** Return the number of MD steps to perform per hybrid MC move */
int HybridMC::nDynamicsSteps()
{
    return nsteps;
}

/** Set the number of MD steps to perform per hybrid MC move */
void HybridMC::setNDynamicsSteps(int n)
{
    if (n > 0)
        nsteps = n;
    else
        nsteps = 0;
}

/** Perform 'nmoves' hybrid MC moves on the passed system, recording
    statistics is 'record_stats' is true */
void HybridMC::move(System &system, int nmoves, bool record_stats)
{
    if (nmoves <= 0 or nsteps == 0)
        return;

    //save our, and the system's, current state
    HybridMC old_state(*this);

    System old_system_state(system);
    
    HMCVelGen &gen = velgen.edit().asA<HMCVelGen>();
    
    try
    {
        gen.setTemperature( this->temperature() );
    
        for (int i=0; i<nmoves; ++i)
        {
            //save the old dynamics move
            MolecularDynamics old_md = md;

            //get the old total energy of the system
            double old_nrg = system.energy( this->energyComponent() );

            //save the old system
            System old_system(system);
        
            //regenerate random velocities for this temperature
            double old_bias = gen.generate(system, md);
        
            qDebug() << old_nrg << md.kineticEnergy().value()
                     << (old_nrg + md.kineticEnergy().value());

            qDebug() << old_bias;

            //run the dynamics moves (we don't record statistics
            //during the MD moves)
            md.move(system, nsteps, false);
    
            //calculate the energy of the system
            double new_nrg = system.energy( this->energyComponent() );
            
            //calculate the new kinetic energy
            double new_bias = gen.getBias(md);

            qDebug() << new_nrg << md.kineticEnergy().value()
                     << (new_nrg + md.kineticEnergy().value());

            qDebug() << new_bias;

            //accept or reject the move
            if (not this->test(new_nrg, old_nrg, new_bias, old_bias))
            {
                //the move has been rejected - reset the state
                system = old_system;
                md = old_md;
            }

            if (record_stats)
            {
                system.collectStats();
            }
        }
    }
    catch(...)
    {
        system = old_system_state;
        this->operator=(old_state);

        throw;
    }
}
