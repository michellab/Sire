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

#include "velocitygenerator.h"

#include "SireMol/moleculeview.h"
#include "SireMol/moleculedata.h"
#include "SireMol/moleculeinfodata.h"
#include "SireMol/atomelements.h"
#include "SireMol/atommasses.h"
#include "SireMol/atomvelocities.h"

#include "SireCAS/symbol.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/units.h"
#include "SireUnits/temperature.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMove;
using namespace SireMol;
using namespace SireCAS;
using namespace SireMaths;
using namespace SireBase;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

/////////
///////// Implementation of VelocityGenerator
/////////

static const RegisterMetaType<VelocityGenerator> r_velgen( MAGIC_ONLY,
                                                    VelocityGenerator::typeName() );
                                                    
/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                        const VelocityGenerator &velgen)
{
    writeHeader(ds, r_velgen, 1);
    
    ds << static_cast<const Property&>(velgen);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, VelocityGenerator &velgen)
{
    VersionID v = readHeader(ds, r_velgen);
    
    if (v == 1)
    {
        ds >> static_cast<Property&>(velgen);
    }
    else
        throw version_error( v, "1", r_velgen, CODELOC );
        
    return ds;
}

/** Constructor */
VelocityGenerator::VelocityGenerator() : Property()
{}

/** Copy constructor */
VelocityGenerator::VelocityGenerator(const VelocityGenerator &other)
                  : Property(other)
{}

/** Destructor */
VelocityGenerator::~VelocityGenerator()
{}

/** Copy assignment operator */
VelocityGenerator& VelocityGenerator::operator=(const VelocityGenerator &other)
{
    Property::operator=(other);
    return *this;
}

/** Comparison operator */
bool VelocityGenerator::operator==(const VelocityGenerator&) const
{
    return true;
}

/** Comparison operator */
bool VelocityGenerator::operator!=(const VelocityGenerator&) const
{
    return false;
}

/** Set the random number generator that may be used to help
    generate the initial velocities */
void VelocityGenerator::setGenerator(const RanGenerator&)
{}

Q_GLOBAL_STATIC( NullVelocityGenerator, nullVelocityGenerator )

/** Return the global null generator */
const NullVelocityGenerator& VelocityGenerator::null()
{
    return *(nullVelocityGenerator());
}

/////////
///////// Implementation of NullVelocityGenerator
/////////

static const RegisterMetaType<NullVelocityGenerator> r_nullvelgen;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, 
                                        const NullVelocityGenerator &nullvelgen)
{
    writeHeader(ds, r_nullvelgen, 1);
    
    ds << static_cast<const VelocityGenerator&>(nullvelgen);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, 
                                        NullVelocityGenerator &nullvelgen)
{
    VersionID v = readHeader(ds, r_nullvelgen);
    
    if (v == 1)
    {
        ds >> static_cast<VelocityGenerator&>(nullvelgen);
    }
    else
        throw version_error(v, "1", r_nullvelgen, CODELOC);
        
    return ds;
}

/** Constructor */
NullVelocityGenerator::NullVelocityGenerator()
                      : ConcreteProperty<NullVelocityGenerator,VelocityGenerator>()
{}

/** Copy constructor */
NullVelocityGenerator::NullVelocityGenerator(const NullVelocityGenerator &other)
                      : ConcreteProperty<NullVelocityGenerator,VelocityGenerator>(other)
{}

/** Destructor */
NullVelocityGenerator::~NullVelocityGenerator()
{}

/** Copy assignment operator */
NullVelocityGenerator& 
NullVelocityGenerator::operator=(const NullVelocityGenerator &other)
{
    VelocityGenerator::operator=(other);
    return *this;
}

/** Comparison operator */
bool NullVelocityGenerator::operator==(const NullVelocityGenerator &other) const
{
    return VelocityGenerator::operator==(other);
}

/** Comparison operator */
bool NullVelocityGenerator::operator!=(const NullVelocityGenerator &other) const
{
    return VelocityGenerator::operator!=(other);
}

/** Zero velocities are generated */
AtomVelocities NullVelocityGenerator::generate(const MoleculeView &molview,
                                               const PropertyMap &map) const
{
    return AtomVelocities( molview.data().info() );
}

const char* NullVelocityGenerator::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullVelocityGenerator>() );
}

/////////
///////// Implementation of VelocitiesFromProperty
/////////

static const RegisterMetaType<VelocitiesFromProperty> r_velfromprop;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, 
                                        const VelocitiesFromProperty &velfromprop)
{
    writeHeader(ds, r_velfromprop, 1);
    
    SharedDataStream sds(ds);
    
    sds << velfromprop.vel_property
        << static_cast<const VelocityGenerator&>(velfromprop);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                        VelocitiesFromProperty &velfromprop)
{
    VersionID v = readHeader(ds, r_velfromprop);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> velfromprop.vel_property
            >> static_cast<VelocityGenerator&>(velfromprop);
    }
    else
        throw version_error(v, "1", r_velfromprop, CODELOC);
        
    return ds;
}

/** Constructor */
VelocitiesFromProperty::VelocitiesFromProperty()
                       : ConcreteProperty<VelocitiesFromProperty,VelocityGenerator>()
{}

/** Construct to get the velocities from the property 'property' */
VelocitiesFromProperty::VelocitiesFromProperty(const PropertyName &property)
                       : ConcreteProperty<VelocitiesFromProperty,VelocityGenerator>(),
                         vel_property(property)
{}

/** Copy constructor */
VelocitiesFromProperty::VelocitiesFromProperty(const VelocitiesFromProperty &other)
                : ConcreteProperty<VelocitiesFromProperty,VelocityGenerator>(other),
                  vel_property(other.vel_property)
{}

/** Destructor */
VelocitiesFromProperty::~VelocitiesFromProperty()
{}

/** Copy assignment operator */
VelocitiesFromProperty& 
VelocitiesFromProperty::operator=(const VelocitiesFromProperty &other)
{
    VelocityGenerator::operator=(other);
    vel_property = other.vel_property;
    
    return *this;
}

/** Comparison operator */
bool VelocitiesFromProperty::operator==(const VelocitiesFromProperty &other) const
{
    return vel_property == other.vel_property and
           VelocityGenerator::operator==(other);
}

/** Comparison operator */
bool VelocitiesFromProperty::operator!=(const VelocitiesFromProperty &other) const
{
    return vel_property != other.vel_property or
           VelocityGenerator::operator!=(other);
}

/** Return the velocities from the specified property */
AtomVelocities VelocitiesFromProperty::generate(const MoleculeView &molview,
                                                const PropertyMap&) const
{
    return molview.data().property(vel_property).asA<AtomVelocities>();
}

const char* VelocitiesFromProperty::typeName()
{
    return QMetaType::typeName( qMetaTypeId<VelocitiesFromProperty>() );
}

/////////
///////// Implementation of MaxwellBoltzmann
/////////

static const RegisterMetaType<MaxwellBoltzmann> r_maxboltz;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, 
                                        const MaxwellBoltzmann &maxboltz)
{
    writeHeader(ds, r_maxboltz, 1);
    
    SharedDataStream sds(ds);
    
    sds << maxboltz.ran_generator << maxboltz.temp.to(kelvin)
        << static_cast<const VelocityGenerator&>(maxboltz);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, MaxwellBoltzmann &maxboltz)
{
    VersionID v = readHeader(ds, r_maxboltz);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        double temp;
        
        sds >> maxboltz.ran_generator >> temp
            >> static_cast<VelocityGenerator&>(maxboltz);
            
        maxboltz.temp = temp*kelvin;
    }
    else
        throw version_error(v, "1", r_maxboltz);
        
    return ds;
}

/** Construct to generate velocities correct for a temperature of 25 C */
MaxwellBoltzmann::MaxwellBoltzmann()
                 : ConcreteProperty<MaxwellBoltzmann,VelocityGenerator>(),
                   temp(25*celsius)
{}

/** Construct to generate velocities correct for the passed temperature */
MaxwellBoltzmann::MaxwellBoltzmann(Temperature temperature)
                 : ConcreteProperty<MaxwellBoltzmann,VelocityGenerator>(),
                   temp(temperature)
{}

/** Copy constructor */
MaxwellBoltzmann::MaxwellBoltzmann(const MaxwellBoltzmann &other)
                 : ConcreteProperty<MaxwellBoltzmann,VelocityGenerator>(other),
                   ran_generator(other.ran_generator), temp(other.temp)
{}

/** Destructor */
MaxwellBoltzmann::~MaxwellBoltzmann()
{}

/** Copy assignment operator */
MaxwellBoltzmann& MaxwellBoltzmann::operator=(const MaxwellBoltzmann &other)
{
    if (this != &other)
    {
        VelocityGenerator::operator=(other);

        ran_generator = other.ran_generator;
        temp = other.temp;
    }
    
    return *this;
}

/** Comparison operator */
bool MaxwellBoltzmann::operator==(const MaxwellBoltzmann &other) const
{
    return temp == other.temp and
           VelocityGenerator::operator==(other);
}

/** Comparison operator */
bool MaxwellBoltzmann::operator!=(const MaxwellBoltzmann &other) const
{
    return not this->operator==(other);
}

/** Set the random number generator used to generate the random velocities */
void MaxwellBoltzmann::setGenerator(const RanGenerator &rangenerator)
{
    ran_generator = rangenerator;
}

/** Return the random number generator used to generate random velocities */
const RanGenerator& MaxwellBoltzmann::generator() const
{
    return ran_generator;
}

/** Return the temperature for which the velocities will be generated */
Temperature MaxwellBoltzmann::temperature() const
{
    return temp;
}

/** Set the temperature at which the velocities will be generated */
void MaxwellBoltzmann::setTemperature(Temperature temperature)
{
    temp = temperature;
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

/** Generate completely random velocities */
AtomVelocities MaxwellBoltzmann::generate(const MoleculeView &molview,
                                          const PropertyMap &map) const
{
    const double kT = temp.to(kelvin) * k_boltz;

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
        //generate random velocities from a Maxwell-Boltzmann distribution.
        //
        // First generate a random vector from a normal distribution
        // with 0 mean and unit variance
        Vector norm_rand( ran_generator.randNorm(0, 1),
                          ran_generator.randNorm(0, 1),
                          ran_generator.randNorm(0, 1) );
        
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

const char* MaxwellBoltzmann::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MaxwellBoltzmann>() );
}
