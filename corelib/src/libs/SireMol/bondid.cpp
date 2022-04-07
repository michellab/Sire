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

#include "bondid.h"

#include "atomcoords.h"
#include "moleculedata.h"
#include "moleculeinfodata.h"

#include "SireBase/property.h"

#include "SireMaths/vector.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireID;
using namespace SireMaths;
using namespace SireStream;

static const RegisterMetaType<BondID> r_bondid;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                       const BondID &bondid)
{
    writeHeader(ds, r_bondid, 1);

    SharedDataStream sds(ds);

    sds << bondid.atm0 << bondid.atm1;

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                       BondID &bondid)
{
    VersionID v = readHeader(ds, r_bondid);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> bondid.atm0 >> bondid.atm1;
    }
    else
        throw version_error(v, "1", r_bondid, CODELOC);

    return ds;
}

/** Null constructor */
BondID::BondID() : ID()
{}

/** Construct a bond between the two specified atoms. The order
    is important, as this bond may be between two different
    molecules */
BondID::BondID(const AtomID &atom0, const AtomID &atom1)
       : ID(), atm0(atom0), atm1(atom1)
{}

/** Copy constructor */
BondID::BondID(const BondID &other)
       : ID(other), atm0(other.atm0), atm1(other.atm1)
{}

/** Destructor */
BondID::~BondID()
{}

/** Copy assignment operator */
BondID& BondID::operator=(const BondID &other)
{
    atm0 = other.atm0;
    atm1 = other.atm1;

    return *this;
}

/** Comparison operator - the order is important */
bool BondID::operator==(const BondID &other) const
{
    return atm0 == other.atm0 and atm1 == other.atm1;
}

/** Comparison operator - the order is important */
bool BondID::operator!=(const BondID &other) const
{
    return atm0 != other.atm0 or atm1 != other.atm1;
}

/** Return the mirror of this BondID - i.e. if this is
    Bond(atom0, atom1), this returns Bond(atom1, atom0).

    This is useful if you know that Bond(atom0,atom1) equals
    Bond(atom1,atom0), e.g. you can now write;

    if (not (bonds.contains(bond) or bonds.contains(bond.mirror())) )
    {
        bonds.insert(bond);
    }

    or

    if (bond == other_bond or bond.mirror() == other.bond())
    {
        //this is the same bond
    }
*/
BondID BondID::mirror() const
{
    return BondID(atm1, atm0);
}

/** Return a hash for this ID */
uint BondID::hash() const
{
    return (atm0.hash() << 16) | (atm1.hash() & 0x0000FFFF);
}

/** Return a string representation of this ID */
QString BondID::toString() const
{
    return QString("Bond( %1, %2 )").arg(atm0.toString(), atm1.toString());
}

/** Return whether this is a null ID */
bool BondID::isNull() const
{
    return atm0.isNull() and atm1.isNull();
}

/** Comparison operator with another ID */
bool BondID::operator==(const SireID::ID &other) const
{
    const BondID *other_bond = dynamic_cast<const BondID*>(&other);

    return other_bond and this->operator==(*other_bond);
}

/** Return the indicies of the two atoms in this bond - this returns
    them in the order tuple(bond.atom0(),bond.atom1())

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
tuple<AtomIdx,AtomIdx> BondID::map(const MoleculeInfoData &molinfo) const
{
    return tuple<AtomIdx,AtomIdx>( molinfo.atomIdx(atm0),
                                   molinfo.atomIdx(atm1) );
}

/** Return the indicies of the two atoms of this bond, between the
    two molecules whose data is in 'mol0info' (containing bond.atom0())
    and 'mol1info' (containing bond.atom1())

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
tuple<AtomIdx,AtomIdx> BondID::map(const MoleculeInfoData &mol0info,
                                   const MoleculeInfoData &mol1info) const
{
    return tuple<AtomIdx,AtomIdx>( mol0info.atomIdx(atm0),
                                   mol1info.atomIdx(atm1) );
}

/** Return a BondID that comprises two AtomIdx IDs, in AtomIdx order
    (the lowest index atom is atom0)
*/
BondID BondID::mapToOrderedBondIdx(const MoleculeInfoData &molinfo) const
{
    auto atom0 = molinfo.atomIdx(atm0);
    auto atom1 = molinfo.atomIdx(atm1);

    if (atom0 <= atom1)
    {
        return BondID(atom0, atom1);
    }
    else
    {
        return BondID(atom1, atom0);
    }
}

/** Return the vector that goes from atom0() to atom1() in the
    molecule whose data is in 'moldata', using the supplied
    property map to find the property that contains the
    coordinates to be used

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Vector BondID::vector(const MoleculeData &moldata,
                      const PropertyMap &map) const
{
    const AtomCoords &coords = moldata.property(map["coordinates"])
                                             .asA<AtomCoords>();

    return coords.at( moldata.info().cgAtomIdx(atm1) ) -
           coords.at( moldata.info().cgAtomIdx(atm0) );
}

/** Return the vector that goes from atom0() in the molecule
    whose data is in 'mol0data' to atom1() in the molecule
    whose data is in 'mol1data', using map0 to find the
    coordinates property of 'mol0' and map1 to find the
    coordinates property of 'mol1'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Vector BondID::vector(const MoleculeData &mol0data,
                      const PropertyMap &map0,
                      const MoleculeData &mol1data,
                      const PropertyMap &map1) const
{
    const AtomCoords &coords0 = mol0data.property(map0["coordinates"])
                                      .asA<AtomCoords>();

    const AtomCoords &coords1 = mol1data.property(map1["coordinates"])
                                      .asA<AtomCoords>();

    return coords1.at( mol1data.info().cgAtomIdx(atm1) ) -
           coords0.at( mol0data.info().cgAtomIdx(atm0) );
}

/** Return the vector that goes from atom0() in the molecule
    whose data is in 'mol0data' to atom1() in the molecule
    whose data is in 'mol1data', using the supplied map
    to find the coordinates property in both molecules.

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Vector BondID::vector(const MoleculeData &mol0data,
                      const MoleculeData &mol1data,
                      const PropertyMap &map) const
{
    return this->vector(mol0data, map, mol1data, map);
}

/** Return the length of this bond in the molecule whose data
    is in 'moldata', using 'map' to find the coordinates property

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
double BondID::length(const MoleculeData &moldata,
                      const PropertyMap &map) const
{
    return this->vector(moldata,map).length();
}

/** Return the length of the bond from atom0() in the
    molecule whose data is in 'mol0data' to atom1() in the
    molecule whose data is in 'mol1data', using 'map0'
    to the find the coordinates property of 'mol0' and
    'map1' to find the coordinates property of 'mol1'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
double BondID::length(const MoleculeData &mol0data,
                      const PropertyMap &map0,
                      const MoleculeData &mol1data,
                      const PropertyMap &map1) const
{
    return this->vector(mol0data, map0, mol1data, map1).length();
}

/** Return the length of the bond from atom0() in the
    molecule whose data is in 'mol0data' to atom1() in the
    molecule whose data is in 'mol1data', using 'map'
    to the find the coordinates properties both molecules

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
double BondID::length(const MoleculeData &mol0data,
                      const MoleculeData &mol1data,
                      const PropertyMap &map) const
{
    return this->vector(mol0data, mol1data, map).length();
}

/** Synonym for BondID::length(const MoleculeData&, const PropertyMap&)

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
double BondID::size(const MoleculeData &moldata,
                    const PropertyMap &map) const
{
    return this->length(moldata, map);
}

/** Synonym for BondID::length(const MoleculeData&, const MoleculeData&,
                               const PropertyMap&)

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
double BondID::size(const MoleculeData &mol0data,
                    const MoleculeData &mol1data,
                    const PropertyMap &map) const
{
    return this->length(mol0data, mol1data, map);
}

/** Synonym for BondID::length(const MoleculeData&, const PropertyMap&,
                               const MoleculeData&, const PropertyMap&)

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
double BondID::size(const MoleculeData &mol0data,
                    const PropertyMap &map0,
                    const MoleculeData &mol1data,
                    const PropertyMap &map1) const
{
    return this->length(mol0data, map0, mol1data, map1);
}

/** Return the ID of the first atom of the bond */
const AtomID& BondID::atom0() const
{
    return atm0.base();
}

/** Return the ID of the second atom of the bond */
const AtomID& BondID::atom1() const
{
    return atm1.base();
}

const char* BondID::typeName()
{
    return QMetaType::typeName( qMetaTypeId<BondID>() );
}

BondID* BondID::clone() const
{
    return new BondID(*this);
}

