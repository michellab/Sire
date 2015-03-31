/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#include "angleid.h"

#include "atomcoords.h"
#include "moleculedata.h"
#include "moleculeinfodata.h"

#include "SireBase/property.h"

#include "SireMaths/vector.h"
#include "SireMaths/triangle.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireID;
using namespace SireMaths;
using namespace SireStream;

using SireUnits::Dimension::Angle;

static const RegisterMetaType<AngleID> r_angleid;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds,
                                       const AngleID &angleid)
{
    writeHeader(ds, r_angleid, 1);
    
    SharedDataStream sds(ds);
    
    sds << angleid.atm0 << angleid.atm1 << angleid.atm2;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds,
                                       AngleID &angleid)
{
    VersionID v = readHeader(ds, r_angleid);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> angleid.atm0 >> angleid.atm1 >> angleid.atm2;
    }
    else
        throw version_error(v, "1", r_angleid, CODELOC);
        
    return ds;
}

/** Null constructor */
AngleID::AngleID() : ID()
{}

/** Construct a angle between the two specified atoms. The order
    is important, as this angle may be between two different
    molecules */
AngleID::AngleID(const AtomID &atom0, const AtomID &atom1,
                 const AtomID &atom2)
       : ID(), atm0(atom0), atm1(atom1), atm2(atom2)
{}

/** Copy constructor */
AngleID::AngleID(const AngleID &other)
       : ID(other), atm0(other.atm0), atm1(other.atm1),
                    atm2(other.atm2)
{}

/** Destructor */
AngleID::~AngleID()
{}

/** Copy assignment operator */
AngleID& AngleID::operator=(const AngleID &other)
{
    atm0 = other.atm0;
    atm1 = other.atm1;
    atm2 = other.atm2;
    
    return *this;
}

/** Comparison operator - the order is important */
bool AngleID::operator==(const AngleID &other) const
{
    return atm0 == other.atm0 and atm1 == other.atm1 and
           atm2 == other.atm2;
}

/** Comparison operator - the order is important */
bool AngleID::operator!=(const AngleID &other) const
{
    return atm0 != other.atm0 or atm1 != other.atm1 or
           atm2 != other.atm2;
}

/** Return the mirror of this AngleID - i.e. if this is 
    AngleID(atom0, atom1, atom2), this returns 
    AngleID(atom2, atom1, atom0).
    
    This is useful if you know that AngleID(atom0,atom1,atom2) equals
    AngleID(atom2,atom1,atom0), e.g. you can now write;
    
    if (not (angles.contains(angle) or angles.contains(angle.mirror())) )
    {
        angles.insert(angle);
    }
    
    or
    
    if (angle == other_angle or angle.mirror() == other.angle())
    {
        //this is the same angle
    }
*/
AngleID AngleID::mirror() const
{
    return AngleID(atm2, atm1, atm0);
}

/** Return a hash for this ID */
uint AngleID::hash() const
{
    return (atm0.hash() << 16) | ( (atm1.hash()*atm2.hash()) & 0x0000FFFF);
}

/** Return a string representation of this ID */
QString AngleID::toString() const
{
    return QString("Angle( %1, %2, %3 )")
                .arg(atm0.toString(), atm1.toString(),
                     atm2.toString());
}

/** Return whether this is a null ID */
bool AngleID::isNull() const
{
    return atm0.isNull() and atm1.isNull() and
           atm2.isNull();
}

/** Comparison operator with another ID */
bool AngleID::operator==(const SireID::ID &other) const
{
    const AngleID *other_angle = dynamic_cast<const AngleID*>(&other);
    
    return other_angle and this->operator==(*other_angle);
}

/** Return the indicies of the three atoms in this angle - this returns 
    them in the order tuple(angle.atom0(),angle.atom1(),angle.atom2())
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
tuple<AtomIdx,AtomIdx,AtomIdx> 
AngleID::map(const MoleculeInfoData &molinfo) const
{
    return tuple<AtomIdx,AtomIdx,AtomIdx>( molinfo.atomIdx(atm0),
                                           molinfo.atomIdx(atm1),
                                           molinfo.atomIdx(atm2) );
}

/** Return the indicies of the three atoms of this angle, between the
    two molecules whose data is in 'mol0info' (containing angle.atom0()),
    'mol1info' (containing angle.atom1()) and 'mol2info' (containing
    angle.atom2())
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
tuple<AtomIdx,AtomIdx,AtomIdx>
AngleID::map(const MoleculeInfoData &mol0info,
             const MoleculeInfoData &mol1info,
             const MoleculeInfoData &mol2info) const
{
    return tuple<AtomIdx,AtomIdx,AtomIdx>( mol0info.atomIdx(atm0),
                                           mol1info.atomIdx(atm1),
                                           mol2info.atomIdx(atm2) );
}

/** Return the geometric triangle formed by the three atoms
    of this angle in the molecule whose data is in 'moldata',
    using 'map' to find the coordinates property.
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/                      
Triangle AngleID::triangle(const MoleculeData &moldata,
                           const PropertyMap &map) const
{
    const AtomCoords &coords = moldata.property(map["coordinates"])
                                    .asA<AtomCoords>();
                                    
    return Triangle( coords.at( moldata.info().cgAtomIdx(atm0) ),
                     coords.at( moldata.info().cgAtomIdx(atm1) ),
                     coords.at( moldata.info().cgAtomIdx(atm2) ) );
}

/** Return the geometric triangle formed by the three atoms,
    atom0() in the molecule whose data is in 'mol0data',
    atom1() from 'mol1data' and atom2() from 'mol2data',
    using 'map0' to find the coordinates property of mol0,
    'map1' to find the coordinates property of mol1 and
    'map2' to find the coordinates property of mol2.
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/                      
Triangle AngleID::triangle(const MoleculeData &mol0data,
                           const PropertyMap &map0,
                           const MoleculeData &mol1data,
                           const PropertyMap &map1,
                           const MoleculeData &mol2data,
                           const PropertyMap &map2) const
{
    const AtomCoords &coords0 = mol0data.property(map0["coordinates"])
                                    .asA<AtomCoords>();
    const AtomCoords &coords1 = mol1data.property(map1["coordinates"])
                                    .asA<AtomCoords>();
    const AtomCoords &coords2 = mol2data.property(map2["coordinates"])
                                    .asA<AtomCoords>();

    return Triangle( coords0.at( mol0data.info().cgAtomIdx(atm0) ),
                     coords1.at( mol1data.info().cgAtomIdx(atm1) ),
                     coords2.at( mol2data.info().cgAtomIdx(atm2) ) );
}
                  
/** Return the geometric triangle formed by the three atoms,
    atom0() in the molecule whose data is in 'mol0data',
    atom1() from 'mol1data' and atom2() from 'mol2data',
    using 'map' to find the coordinates property of the
    molecules
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/                      
Triangle AngleID::triangle(const MoleculeData &mol0data,
                           const MoleculeData &mol1data,
                           const MoleculeData &mol2data,
                           const PropertyMap &map) const
{
    return this->triangle(mol0data, map,
                          mol1data, map,
                          mol2data, map);
}
                           
/** Return the vector that is perpendicular to the plane
    formed by atoms atom0(), atom1() and atom2() in the 
    molecule whose data is in 'moldata', using the supplied
    property map to find the property that contains the 
    coordinates to be used
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Vector AngleID::vector(const MoleculeData &moldata,
                       const PropertyMap &map) const
{
    return this->triangle(moldata, map).vector();
}

/** Return the vector that is perpendicular to the plane
    formed by the atoms atom0() in 'mol0data', atom1() in
    'mol1data' and atom2() in 'mol2data', using map0 to find the 
    coordinates property of 'mol0', map1 to find the 
    coordinates property of 'mol1' and map2 to find the 
    coordinates property of 'mol2'
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Vector AngleID::vector(const MoleculeData &mol0data, 
                       const PropertyMap &map0,
                       const MoleculeData &mol1data,
                       const PropertyMap &map1,
                       const MoleculeData &mol2data,
                       const PropertyMap &map2) const
{
    return this->triangle(mol0data, map0, mol1data, map1,
                          mol2data, map2).vector();
}

/** Return the vector that is perpendicular to the plane
    formed by the atoms atom0() in 'mol0data', atom1() in
    'mol1data' and atom2() in 'mol2data', using 'map' to find the 
    coordinates property of the molecules
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Vector AngleID::vector(const MoleculeData &mol0data, 
                       const MoleculeData &mol1data,
                       const MoleculeData &mol2data,
                       const PropertyMap &map) const
{
    return this->triangle(mol0data, mol1data, mol2data, map).vector();
}
     
/** Return the size of this angle in the molecule whose data
    is in 'moldata', using 'map' to find the coordinates property
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Angle AngleID::size(const MoleculeData &moldata,
                    const PropertyMap &map) const
{
    return this->triangle(moldata,map).angle();
}

/** Return the size of the angle between atom0() in the 
    molecule whose data is in 'mol0data', atom1() in the 
    molecule whose data is in 'mol1data' and atom2() in 
    the molecule whose data is in 'mol2data', using 'map0'
    to the find the coordinates property of 'mol0',
    'map1' to find the coordinates property of 'mol1'
    and 'map2' to find the coordinates property of 'mol2'
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Angle AngleID::size(const MoleculeData &mol0data,
                    const PropertyMap &map0,
                    const MoleculeData &mol1data,
                    const PropertyMap &map1,
                    const MoleculeData &mol2data,
                    const PropertyMap &map2) const
{
    return this->triangle(mol0data, map0, 
                          mol1data, map1,
                          mol2data, map2).angle();
}

/** Return the size of the angle between atom0() in the 
    molecule whose data is in 'mol0data', atom1() in the 
    molecule whose data is in 'mol1data' and atom2() in 
    the molecule whose data is in 'mol2data', using 'map'
    to the find the coordinates property the molecules
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Angle AngleID::size(const MoleculeData &mol0data,
                    const MoleculeData &mol1data,
                    const MoleculeData &mol2data,
                    const PropertyMap &map) const
{
    return this->triangle(mol0data, mol1data, mol2data, map).angle();
}

/** Return the ID of the first atom of the angle */
const AtomID& AngleID::atom0() const
{
    return atm0.base();
}

/** Return the ID of the second atom of the angle */
const AtomID& AngleID::atom1() const
{
    return atm1.base();
}

/** Return the ID of the third atom of the angle */
const AtomID& AngleID::atom2() const
{
    return atm2.base();
}

const char* AngleID::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AngleID>() );
}

AngleID* AngleID::clone() const
{
    return new AngleID(*this);
}

