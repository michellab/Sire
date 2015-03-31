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

#include "molinserter.h"

#include "SireMol/molecule.h"
#include "SireMol/partialmolecule.h"

#include "SireVol/space.h"

#include "SireSystem/system.h"

#include "SireMaths/quaternion.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/units.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMove;
using namespace SireSystem;
using namespace SireMol;
using namespace SireVol;
using namespace SireMaths;
using namespace SireBase;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

//////////
////////// Implementation of MolInserter
//////////

static const RegisterMetaType<MolInserter> r_molinserter( MAGIC_ONLY,
                                                          MolInserter::typeName() );
                                                          
/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds,
                                        const MolInserter &molinserter)
{
    writeHeader(ds, r_molinserter, 1);
    
    SharedDataStream sds(ds);
    
    sds << molinserter.rangen << molinserter.mgids
        << static_cast<const Property&>(molinserter);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, MolInserter &molinserter)
{
    VersionID v = readHeader(ds, r_molinserter);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> molinserter.rangen >> molinserter.mgids
            >> static_cast<Property&>(molinserter);
            
        molinserter.rebuildCoordsProperties();
    }
    else
        throw version_error(v, "1", r_molinserter, CODELOC);
        
    return ds;
}

/** Constructor */
MolInserter::MolInserter() : Property()
{}

/** Copy constructor */
MolInserter::MolInserter(const MolInserter &other)
            : Property(other),
              rangen(other.rangen), mgids(other.mgids), 
              coords_properties(other.coords_properties)
{}

/** Destructor */
MolInserter::~MolInserter()
{}

/** Copy assignment operator */
MolInserter& MolInserter::operator=(const MolInserter &other)
{
    if (this != &other)
    {
        rangen = other.rangen;
        mgids = other.mgids;
        coords_properties = other.coords_properties;
        
        Property::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool MolInserter::operator==(const MolInserter &other) const
{
    return mgids == other.mgids and Property::operator==(other);
}

/** Comparison operator */
bool MolInserter::operator!=(const MolInserter &other) const
{
    return not MolInserter::operator==(other);
}

const char* MolInserter::typeName()
{
    return "SireMove::MolInserter";
}

/** Internal function used to get all of the coordinates properties
    that need to be used when inserting the molecule into the groups */
void MolInserter::rebuildCoordsProperties()
{
    int ngroups = mgids.count();
    const PropertyMap *maps_array = mgids.propertyMaps().constData();

    QSet<QString> coords_set;
    
    for (int i=0; i<ngroups; ++i)
    {
        const PropertyName &coords_property = maps_array[i]["coordinates"];
        
        if (coords_property.hasSource())
            coords_set.insert( coords_property.source() );
    }
    
    coords_properties = coords_set.toList();
}

/** Internal function used to return the different coordinates 
    properties needed for the insert */
const QStringList& MolInserter::coordsProperties() const
{
    return coords_properties;
}

/** Set the random number generator used to generate the random 
    numbers used to insert the molecule */
void MolInserter::setGenerator(const RanGenerator &generator)
{
    rangen = generator;
}

/** Set the groups (and associated properties) used when adding
    the molecule to the system */
void MolInserter::setGroups(const MGIDsAndMaps &group_ids)
{
    mgids = group_ids;
    this->rebuildCoordsProperties();
}

/** Return the random number generator used to generate the random  
    numbers used when inserting the molecule */
const RanGenerator& MolInserter::generator() const
{
    return rangen;
}

/** Return the group IDs (and associated) properties to which
    the molecule will be inserted */
const MGIDsAndMaps& MolInserter::groups() const
{
    return mgids;
}

static SharedPolyPointer<NullInserter> shared_null;

/** Return the global null MolInserter */
const NullInserter& MolInserter::null()
{
    if (shared_null.constData() == 0)
    {
        QMutexLocker lkr( SireBase::globalLock() );
        
        if (shared_null.constData() == 0)
            shared_null = new NullInserter();
    }
    
    return *(shared_null.constData());
}

//////////
////////// Implementation of NullInserter
//////////

static const RegisterMetaType<NullInserter> r_nullinserter;

/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds,
                                        const NullInserter &nullinserter)
{
    writeHeader(ds, r_nullinserter, 1);
    
    ds << static_cast<const MolInserter&>(nullinserter);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, NullInserter &nullinserter)
{
    VersionID v = readHeader(ds, r_nullinserter);
    
    if (v == 1)
    {
        ds >> static_cast<MolInserter&>(nullinserter);
    }
    else
        throw version_error(v, "1", r_nullinserter, CODELOC);
        
    return ds;
}

/** Constructor */
NullInserter::NullInserter() : ConcreteProperty<NullInserter,MolInserter>()
{}

/** Copy constructor */
NullInserter::NullInserter(const NullInserter &other)
             : ConcreteProperty<NullInserter,MolInserter>(other)
{}

/** Destructor */
NullInserter::~NullInserter()
{}

const char* NullInserter::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullInserter>() );
}

/** Copy assignment operator */
NullInserter& NullInserter::operator=(const NullInserter &other)
{
    MolInserter::operator=(other);
    return *this;
}

/** Comparison operator */
bool NullInserter::operator==(const NullInserter &other) const
{
    return MolInserter::operator==(other);
}

/** Comparison operator */
bool NullInserter::operator!=(const NullInserter &other) const
{
    return MolInserter::operator!=(other);
}

/** This does nothing! */
double NullInserter::insert(const Molecule&, System&, const Space&)
{
    return 0;
}
              
/** This does nothing! */
double NullInserter::insert(const PartialMolecule&, System&, const Space&)
{
    return 0;
}

//////////
////////// Implementation of UniformInserter
//////////

static const RegisterMetaType<UniformInserter> r_uniforminserter;

/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds,
                                        const UniformInserter &uniforminserter)
{
    writeHeader(ds, r_uniforminserter, 1);
    
    ds << static_cast<const MolInserter&>(uniforminserter);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds,
                                        UniformInserter &uniforminserter)
{
    VersionID v = readHeader(ds, r_uniforminserter);
    
    if (v == 1)
    {
        ds >> static_cast<MolInserter&>(uniforminserter);
    }
    else
        throw version_error(v, "1", r_uniforminserter, CODELOC);

    return ds;
}

/** Constructor */
UniformInserter::UniformInserter()
                : ConcreteProperty<UniformInserter,MolInserter>()
{}

/** Construct to insert molecules into the groups identified in 'mgids' 
    (using the associated property maps to find the properties needed
    for those insertions) */
UniformInserter::UniformInserter(const MGIDsAndMaps &mgids)
                : ConcreteProperty<UniformInserter,MolInserter>()
{
    UniformInserter::setGroups(mgids);
}

/** Copy constructor */
UniformInserter::UniformInserter(const UniformInserter &other)
                : ConcreteProperty<UniformInserter,MolInserter>(other)
{}

/** Destructor */
UniformInserter::~UniformInserter()
{}

const char* UniformInserter::typeName()
{
    return QMetaType::typeName( qMetaTypeId<UniformInserter>() );
}

/** Copy assignment operator */
UniformInserter& UniformInserter::operator=(const UniformInserter &other)
{
    MolInserter::operator=(other);
    return *this;
}

/** Comparison operator */
bool UniformInserter::operator==(const UniformInserter &other) const
{
    return MolInserter::operator==(other);
}

/** Comparison operator */
bool UniformInserter::operator!=(const UniformInserter &other) const
{
    return MolInserter::operator!=(other);
}

template<class T>
void UniformInserter::uniform_insert(const T &molecule, System &system, 
                                     const Space &space) const
{
    //pick a random point in the space
    Vector insertion_point = space.getRandomPoint( generator() );

    //now pick a random orientation - this is a random vector and 
    //random angle around which to rotate the molecule
    Vector orientation_vector = generator().vectorOnSphere();
    Angle orientation_angle = generator().rand(-two_pi, two_pi) * radians;

    //we now need to move the molecule to this point. This may
    //have to be done multiple times, as there may be multiple
    //coordinates properties (different coordinates properties
    //for different molecule groups)
    
    PropertyMap map;
    
    T moved_mol(molecule);
    
    foreach (const QString &coords_property, coordsProperties())
    {
        map.set("coordinates", coords_property);
        
        //we need to know the center of the molecule as we will 
        //rotate around that, and also as we need to translate
        //the molecule so that its center is at the insertion point
        Vector mol_center = moved_mol.evaluate().center(map);
        
        moved_mol = moved_mol.move()
                             .rotate( Quaternion(orientation_angle, orientation_vector),
                                      mol_center, map )
                             .translate( insertion_point - mol_center, map )
                             .commit();
    }

    //ok - the molecule now has all of the necessary coordinate properties
    // - lets add it to the required molecule groups
    int ngroups = groups().count();
    
    const MGIdentifier *mgids_array = groups().mgIDs().constData();
    const PropertyMap *maps_array = groups().propertyMaps().constData();
    
    for (int i=0; i<ngroups; ++i)
    {
        system.add( moved_mol, mgids_array[i], maps_array[i] );
    }
}

/** This funciton inserts the molecule 'molecule' into 'system' at 
    a random orientation and position within the space 'space' */
double UniformInserter::insert(const Molecule &molecule, System &system,
                               const Space &space)
{
    this->uniform_insert<Molecule>(molecule, system, space);

    return 1;
}

/** This funciton inserts the molecule 'molecule' into 'system' at 
    a random orientation and position within the space 'space' */
double UniformInserter::insert(const PartialMolecule &molecule, System &system,
                               const Space &space)
{
    this->uniform_insert<PartialMolecule>(molecule, system, space);

    return 1;
}
