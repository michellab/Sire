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

#include "volumechanger.h"

#include "SireMol/moleculegroup.h"
#include "SireMol/partialmolecule.h"

#include "SireSystem/system.h"

#include "SireVol/space.h"

#include "SireUnits/units.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "sire_config.h"

#ifdef HAVE_CUBEROOT
    #include "math.h"
#else
    #include <cmath>
#endif

using namespace SireMove;
using namespace SireSystem;
using namespace SireMol;
using namespace SireVol;
using namespace SireFF;
using namespace SireMaths;
using namespace SireBase;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

////////////
//////////// Implementation of VolumeChanger
////////////

static const RegisterMetaType<VolumeChanger> r_volchanger( MAGIC_ONLY,
                                                           VolumeChanger::typeName() );
                                                           
/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const VolumeChanger &volchanger)
{
    writeHeader(ds, r_volchanger, 1);
    
    SharedDataStream sds(ds);
    
    sds << volchanger.rangen << volchanger.mgid
        << static_cast<const Property&>(volchanger);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, VolumeChanger &volchanger)
{
    VersionID v = readHeader(ds, r_volchanger);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> volchanger.rangen >> volchanger.mgid
            >> static_cast<Property&>(volchanger);

    }
    else
        throw version_error( v, "1", r_volchanger, CODELOC );
        
    return ds;
}

/** Null constructor */
VolumeChanger::VolumeChanger() : Property()
{}

/** Construct to operate on the molecule group(s) that 
    match the ID 'mgid' */
VolumeChanger::VolumeChanger(const MGID &mg_id) : Property(), mgid(mg_id)
{}
  
/** Construct to operate on the passed molecule group
    (this matches the group based on its molecule group number) */            
VolumeChanger::VolumeChanger(const MoleculeGroup &molgroup)
              : Property(), mgid(molgroup.number())
{}

/** Copy constructor */
VolumeChanger::VolumeChanger(const VolumeChanger &other)
              : Property(other), mgid(other.mgid)
{}

/** Destructor */
VolumeChanger::~VolumeChanger()
{}

/** Copy assignment operator */
VolumeChanger& VolumeChanger::operator=(const VolumeChanger &other)
{
    if (this != &other)
    {
        rangen = other.rangen;
        mgid = other.mgid;
        Property::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool VolumeChanger::operator==(const VolumeChanger &other) const
{
    return mgid == other.mgid and Property::operator==(other);
}

/** Comparison operator */
bool VolumeChanger::operator!=(const VolumeChanger &other) const
{
    return mgid != other.mgid or Property::operator!=(other);
}

/** Set the random number generator that may be used to generate
    new volumes */
void VolumeChanger::setGenerator(const RanGenerator &generator)
{
    rangen = generator;
}

/** Return the random number generator that may be used to generate
    new volumes */
const RanGenerator& VolumeChanger::generator() const
{
    return rangen;
}

/** Return the ID of the molecule group(s) that will be affected
    by this volume changer */
const MGID& VolumeChanger::groupID() const
{
    return mgid.base();
}

/** Set the ID of the molecule group(s) that will be affected by
    this volume changer */
void VolumeChanger::setGroup(const MGID &new_mgid)
{
    mgid = new_mgid;
}

/** Set the molecule group that is affected by this volume changer
     - this will match the group based on its molecule group number */
void VolumeChanger::setGroup(const MoleculeGroup &molgroup)
{
    mgid = molgroup.number();
}

/** Change the volume of the passed system 'system' by 'delta', using
    the optionally supplied property map to find the names of the 
    necessary properties
    
    This returns the number of molecules which were involved in
    the volume change
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::unsupported
    \throw SireError::incompatible_error
*/
int VolumeChanger::changeVolume(System &system, const Volume &delta,
                                const PropertyMap &map) const
{
    Volume current_volume = system.property( map["space"] ).asA<Space>().volume();
    
    return this->setVolume(system, current_volume + delta, map);
}

/** Change the volume of the passed system 'system' by a random
    amount between -maxvolchange and maxvolchange, using
    the optionally supplied property map to find the names of the 
    necessary properties.

    If this move is biased, then this sets old_bias to the
    bias before the move, and new_bias to the bias afterwards.

    This returns the number of molecules which were involved in
    the volume change
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::unsupported
    \throw SireError::incompatible_error
*/
int VolumeChanger::randomChangeVolume(System &system, const Volume &maxvolchange,
                                      double &new_bias, double &old_bias,
                                      const PropertyMap &map) const
{
    new_bias = 1;
    old_bias = 1;

    return this->changeVolume(system, Volume(rangen.rand(-maxvolchange,maxvolchange)),
                              map);
}

Q_GLOBAL_STATIC( NullVolumeChanger, nullVolumeChanger );

/** Return the global null changer (which doesn't change anything!) */
const NullVolumeChanger& VolumeChanger::null()
{
    return *(nullVolumeChanger());
}

////////////
//////////// Implementation of NullVolumeChanger
////////////

static const RegisterMetaType<NullVolumeChanger> r_nullvolchanger;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, 
                                        const NullVolumeChanger &nullvolchanger)
{
    writeHeader(ds, r_nullvolchanger, 1);
    
    ds << static_cast<const VolumeChanger&>(nullvolchanger);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                        NullVolumeChanger &nullvolchanger)
{
    VersionID v = readHeader(ds, r_nullvolchanger);
    
    if (v == 1)
    {
        ds >> static_cast<VolumeChanger&>(nullvolchanger);
    }
    else 
        throw version_error( v, "1", r_nullvolchanger, CODELOC );
        
    return ds;
}

/** Constructor */
NullVolumeChanger::NullVolumeChanger() 
                  : ConcreteProperty<NullVolumeChanger,VolumeChanger>()
{}

/** Copy constructor */
NullVolumeChanger::NullVolumeChanger(const NullVolumeChanger &other)
                  : ConcreteProperty<NullVolumeChanger,VolumeChanger>(other)
{}

/** Destructor */
NullVolumeChanger::~NullVolumeChanger()
{}

/** Copy assignment operator */
NullVolumeChanger& NullVolumeChanger::operator=(const NullVolumeChanger &other)
{
    VolumeChanger::operator=(other);
    return *this;
}

/** Comparison operator */
bool NullVolumeChanger::operator==(const NullVolumeChanger &other) const
{
    return VolumeChanger::operator==(other);
}

/** Comparison operator */
bool NullVolumeChanger::operator!=(const NullVolumeChanger &other) const
{
    return VolumeChanger::operator!=(other);
}

const char* NullVolumeChanger::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullVolumeChanger>() );
}

/** The null volume changer doesn't change anything */
int NullVolumeChanger::setVolume(System &system, const Volume&, const PropertyMap&) const
{
    return system.molecules(this->groupID()).nMolecules();
}

////////////
//////////// Implementation of ScaleVolumeFromCenter
////////////

static const RegisterMetaType<ScaleVolumeFromCenter> r_scalevol;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                        const ScaleVolumeFromCenter &scalevol)
{
    writeHeader(ds, r_scalevol, 1);
    
    SharedDataStream sds(ds);
    
    sds << scalevol.scale_point
        << static_cast<const VolumeChanger&>(scalevol);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                        ScaleVolumeFromCenter &scalevol)
{
    VersionID v = readHeader(ds, r_scalevol);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> scalevol.scale_point 
            >> static_cast<VolumeChanger&>(scalevol);
    }
    else
        throw version_error( v, "1", r_scalevol, CODELOC );
    
    return ds;
}

/** Null constructor */
ScaleVolumeFromCenter::ScaleVolumeFromCenter()
                      : ConcreteProperty<ScaleVolumeFromCenter,VolumeChanger>()
{}

/** Construct to scale the molecules in the group(s) with ID 'mgid',
    scaling the molecules from the origin (0,0,0) */
ScaleVolumeFromCenter::ScaleVolumeFromCenter(const MGID &mgid)
                      : ConcreteProperty<ScaleVolumeFromCenter,VolumeChanger>(mgid),
                        scale_point( VectorPoint(Vector(0)) )
{}

/** Construct to scale the molecules in the group 'molgroup',
    scaling the molecules from the origin (0,0,0) */
ScaleVolumeFromCenter::ScaleVolumeFromCenter(const MoleculeGroup &molgroup)
                      : ConcreteProperty<ScaleVolumeFromCenter,VolumeChanger>(molgroup),
                        scale_point( VectorPoint(Vector(0)) )
{}

/** Construct to scale the molecules in the groups with ID 'mgid',
    scaling the molecules from the point 'point' */
ScaleVolumeFromCenter::ScaleVolumeFromCenter(const MGID &mgid, const PointRef &point)
                      : ConcreteProperty<ScaleVolumeFromCenter,VolumeChanger>(mgid),
                        scale_point(point)
{}

/** Construct to scale the molecules in the group 'molgroup',
    scaling the molecules from the point 'point' */
ScaleVolumeFromCenter::ScaleVolumeFromCenter(const MoleculeGroup &molgroup, 
                                             const PointRef &point)
                      : ConcreteProperty<ScaleVolumeFromCenter,VolumeChanger>(molgroup),
                        scale_point(point)
{}

/** Copy constructor */
ScaleVolumeFromCenter::ScaleVolumeFromCenter(const ScaleVolumeFromCenter &other)
                      : ConcreteProperty<ScaleVolumeFromCenter,VolumeChanger>(other),
                        scale_point(other.scale_point)
{}

/** Destructor */
ScaleVolumeFromCenter::~ScaleVolumeFromCenter()
{}

/** Copy assignment operator */
ScaleVolumeFromCenter& ScaleVolumeFromCenter::operator=(
                                                    const ScaleVolumeFromCenter &other)
{
    VolumeChanger::operator=(other);
    scale_point = other.scale_point;
    
    return *this;
}

/** Comparison operator */
bool ScaleVolumeFromCenter::operator==(const ScaleVolumeFromCenter &other) const
{
    return VolumeChanger::operator==(other) and 
           scale_point == other.scale_point;
}

/** Comparison operator */
bool ScaleVolumeFromCenter::operator!=(const ScaleVolumeFromCenter &other) const
{
    return VolumeChanger::operator!=(other) or 
           scale_point != other.scale_point;
}

const char* ScaleVolumeFromCenter::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ScaleVolumeFromCenter>() );
}

/** Return the center point from which the molecules are scaled */
const Point& ScaleVolumeFromCenter::center() const
{
    return scale_point.read();
}

/** Set the center point from which the molecules will be scaled */
void ScaleVolumeFromCenter::setCenter(const PointRef &center)
{
    scale_point = center;
}

/** Set the volume of the system 'system' to 'volume', using the 
    optionally supplied property map to find the names of the 
    properties needed to change the system volume.
    
    This returns the number of molecules involved in the volume change
    
    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireVol::incompatible_space
    \throw SireError::invalid_cast
    \throw SireError::invalid_state
    \throw SireError::unsupported
    \throw SireError::incompatible_error
*/
int ScaleVolumeFromCenter::setVolume(System &system, const Volume &volume,
                                     const PropertyMap &map) const
{
    const PropertyName &space_property = map["space"];
    
    //get the molecules to be scaled in the molecule groups
    Molecules old_molecules = system.molecules(this->groupID());

    //get the space
    const Space &old_space = system.property(space_property).asA<Space>();

    const Volume old_volume = old_space.volume();

    if (old_volume <= Volume(0))
        throw SireError::invalid_state( QObject::tr(
            "Cannot set the volume to %1 A^3 as the volume of the space "
            "%2 (%3 A^3) is less than or equal to zero!")
                .arg(volume.to(angstrom3))
                .arg(old_space.toString())
                .arg(old_volume.to(angstrom3)), CODELOC );
    
    else if (old_volume == volume)
        //there is nothing to change - but we still need to return the
        //number of molecules that are part of this volume
        return old_molecules.nMolecules();

    //get the center point from which we will be scaling
    Vector center;
    
    if (scale_point.read().usesMoleculesIn(system))
    {
        //need to do a update on a copy here as we are not allowed
        //to change this volume changer
        PointPtr new_point(scale_point);
        
        new_point.edit().update(system);
        new_point.edit().setSpace(old_space);    //(could this violate detailed balance,
                                                 // as the center could depend on the 
                                                 // current space, so would change,
                                                 // meaning that reverse move probability
                                                 // would be hard to compute...)
        
        center = new_point.read().point();
    }
    else
        center = scale_point.read().point();
	
    //set the new volume of this space
    SpacePtr space = old_space.setVolume( volume );

    //now work out by how much we have to scale the distance between
    //each molecule and the center point (we will scale x,y,z by an equal
    //amount, which is the cube-root of the ratio of the volume change)
    #ifdef HAVE_CUBEROOT
        const double scale_ratio = ::cbrt( volume / old_volume );
    #else
        const double scale_ratio = std::pow( volume / old_volume, 1.0/3.0 );
    #endif
    
    Molecules molecules = old_molecules;
    
    for (Molecules::const_iterator it = old_molecules.constBegin();
         it != old_molecules.constEnd();
         ++it)
    {
        PartialMolecule mol( *it );

		//all we need to do is scale the molecule from the center!
        Vector old_mol_center = mol.evaluate().center(map);

		Vector delta = (scale_ratio - 1) * (old_mol_center - center);
        
        if (not delta.isZero())
        {
            mol = mol.move().translate(delta, map).commit();
            molecules.update(mol);
        }
    }

	if (space_property.hasSource())
	    system.setProperty(space_property.source(), space);    

    system.update(molecules);
    
    return old_molecules.nMolecules();
}
