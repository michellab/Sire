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

#include "improperid.h"

#include "atomcoords.h"
#include "moleculedata.h"
#include "moleculeinfodata.h"

#include "SireBase/property.h"

#include "SireMaths/vector.h"
#include "SireMaths/torsion.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireID;
using namespace SireMaths;
using namespace SireStream;

using SireUnits::Dimension::Angle;

static const RegisterMetaType<ImproperID> r_improperid;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                       const ImproperID &improperid)
{
    writeHeader(ds, r_improperid, 1);
    
    SharedDataStream sds(ds);
    
    sds << improperid.atm0 << improperid.atm1 
        << improperid.atm2 << improperid.atm3;
         
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                       ImproperID &improperid)
{
    VersionID v = readHeader(ds, r_improperid);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> improperid.atm0 >> improperid.atm1 
            >> improperid.atm2 >> improperid.atm3;
    }
    else
        throw version_error(v, "1", r_improperid, CODELOC);
        
    return ds;
}

/** Null constructor */
ImproperID::ImproperID() : ID()
{}

/** Construct a improper between the two specified atoms. The order
    is important, as this improper may be between two different
    molecules */
ImproperID::ImproperID(const AtomID &atom0, const AtomID &atom1,
                       const AtomID &atom2, const AtomID &atom3)
       : ID(), atm0(atom0), atm1(atom1), atm2(atom2), atm3(atom3)
{}

/** Copy constructor */
ImproperID::ImproperID(const ImproperID &other)
       : ID(other), atm0(other.atm0), atm1(other.atm1),
                    atm2(other.atm2), atm3(other.atm3)
{}

/** Destructor */
ImproperID::~ImproperID()
{}

/** Copy assignment operator */
ImproperID& ImproperID::operator=(const ImproperID &other)
{
    atm0 = other.atm0;
    atm1 = other.atm1;
    atm2 = other.atm2;
    atm3 = other.atm3;
    
    return *this;
}

/** Comparison operator - the order is important */
bool ImproperID::operator==(const ImproperID &other) const
{
    return atm0 == other.atm0 and atm1 == other.atm1 and
           atm2 == other.atm2 and atm3 == other.atm3;
}

/** Comparison operator - the order is important */
bool ImproperID::operator!=(const ImproperID &other) const
{
    return atm0 != other.atm0 or atm1 != other.atm1 or
           atm2 != other.atm2 or atm3 != other.atm3;
}

/** Are these impropers generally equivalent, i.e. do they contain the same
    atom indices. This is useful since the ordering of improper atoms is
    inconsistent between different molecular topology formats.
*/
bool ImproperID::equivalent(const ImproperID &other) const
{
    return (atm0 == other.atm0 or atm0 == other.atm1 or atm0 == other.atm2 or atm0 == other.atm3) and
           (atm1 == other.atm0 or atm1 == other.atm1 or atm1 == other.atm2 or atm1 == other.atm3) and
           (atm2 == other.atm0 or atm2 == other.atm1 or atm2 == other.atm2 or atm2 == other.atm3) and
           (atm3 == other.atm0 or atm3 == other.atm1 or atm3 == other.atm2 or atm3 == other.atm3);
}

/** Return a hash for this ID */
uint ImproperID::hash() const
{
    return ( (atm0.hash()*atm1.hash()) << 16) | 
           ( (atm2.hash()*atm3.hash()) & 0x0000FFFF);
}

/** Return a string representation of this ID */
QString ImproperID::toString() const
{
    return QString("Improper( %1, %2, %3, %4 )")
                .arg(atm0.toString(), atm1.toString(),
                     atm2.toString(), atm3.toString());
}

/** Return whether this is a null ID */
bool ImproperID::isNull() const
{
    return atm0.isNull() and atm1.isNull() and
           atm2.isNull() and atm3.isNull();
}

/** Comparison operator with another ID */
bool ImproperID::operator==(const SireID::ID &other) const
{
    const ImproperID *other_improper = dynamic_cast<const ImproperID*>(&other);
    
    return other_improper and this->operator==(*other_improper);
}

/** Return the indicies of the four atoms in this improper - this returns 
    them in the order 
    tuple(improper.atom0(),improper.atom1(),improper.atom2(),improper.atom3())
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
tuple<AtomIdx,AtomIdx,AtomIdx,AtomIdx> 
ImproperID::map(const MoleculeInfoData &molinfo) const
{
    return tuple<AtomIdx,AtomIdx,AtomIdx,AtomIdx>( 
                        molinfo.atomIdx(atm0), molinfo.atomIdx(atm1),
                        molinfo.atomIdx(atm2), molinfo.atomIdx(atm3) );
}

/** Return the indicies of the four atoms of this improper, between the
    two molecules whose data is in 'mol0info' (containing improper.atom0()),
    'mol1info' (containing improper.atom1()), 'mol2info' (containing
    improper.atom2()) and 'mol3info' (containing improper.atom3())
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
tuple<AtomIdx,AtomIdx,AtomIdx,AtomIdx>
ImproperID::map(const MoleculeInfoData &mol0info, 
                const MoleculeInfoData &mol1info,
                const MoleculeInfoData &mol2info, 
                const MoleculeInfoData &mol3info) const
{
    return tuple<AtomIdx,AtomIdx,AtomIdx,AtomIdx>( 
                        mol0info.atomIdx(atm0), mol1info.atomIdx(atm1),
                        mol2info.atomIdx(atm2), mol3info.atomIdx(atm3) );
}

/** Return the geometric torsion formed by the four atoms
    of this improper in the molecule whose data is in 'moldata',
    using 'map' to find the coordinates property.
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/                      
Torsion ImproperID::torsion(const MoleculeData &moldata,
                            const PropertyMap &map) const
{
    const AtomCoords &coords = moldata.property(map["coordinates"])
                                      .asA<AtomCoords>();
                                    
    return Torsion( coords.at( moldata.info().cgAtomIdx(atm0) ),
                    coords.at( moldata.info().cgAtomIdx(atm1) ),
                    coords.at( moldata.info().cgAtomIdx(atm2) ),
                    coords.at( moldata.info().cgAtomIdx(atm3) ) );
}

/** Return the geometric torsion formed by the four atoms,
    atom0() in the molecule whose data is in 'mol0data',
    atom1() from 'mol1data', atom2() from 'mol2data', and
    atom3() from 'mol3data',
    using 'map0' to find the coordinates property of mol0,
    'map1' to find the coordinates property of mol1,
    'map2' to find the coordinates property of mol2 and
    'map3' to find the coordinates property of mol3.
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/                      
Torsion ImproperID::torsion(const MoleculeData &mol0data,
                            const PropertyMap &map0,
                            const MoleculeData &mol1data,
                            const PropertyMap &map1,
                            const MoleculeData &mol2data,
                            const PropertyMap &map2,
                            const MoleculeData &mol3data,
                            const PropertyMap &map3) const
{
    const AtomCoords &coords0 = mol0data.property(map0["coordinates"])
                                        .asA<AtomCoords>();
    const AtomCoords &coords1 = mol1data.property(map1["coordinates"])
                                        .asA<AtomCoords>();
    const AtomCoords &coords2 = mol2data.property(map2["coordinates"])
                                        .asA<AtomCoords>();
    const AtomCoords &coords3 = mol3data.property(map3["coordinates"])
                                        .asA<AtomCoords>();

    return Torsion( coords0.at( mol0data.info().cgAtomIdx(atm0) ),
                    coords1.at( mol1data.info().cgAtomIdx(atm1) ),
                    coords2.at( mol2data.info().cgAtomIdx(atm2) ),
                    coords3.at( mol3data.info().cgAtomIdx(atm3) ) );
}
                  
/** Return the geometric torsion formed by the four atoms,
    atom0() in the molecule whose data is in 'mol0data',
    atom1() from 'mol1data', atom2() from 'mol2data', and
    atom3() from 'mol3data',
    using 'map' to find the coordinates property of 
    the molecules
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/                      
Torsion ImproperID::torsion(const MoleculeData &mol0data,
                            const MoleculeData &mol1data,
                            const MoleculeData &mol2data,
                            const MoleculeData &mol3data,
                            const PropertyMap &map) const
{
    return this->torsion(mol0data, map,
                         mol1data, map,
                         mol2data, map,
                         mol3data, map);
}
     
/** Return the size of this improper in the molecule whose data
    is in 'moldata', using 'map' to find the coordinates property
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Angle ImproperID::size(const MoleculeData &moldata,
                       const PropertyMap &map) const
{
    return this->torsion(moldata,map).angle();
}

/** Return the size of the improper between atom0() in the 
    molecule whose data is in 'mol0data', atom1() in the 
    molecule whose data is in 'mol1data', atom2() in 
    the molecule whose data is in 'mol2data', and
    atom3() in the molecule whose data is in 'mol3data', using 'map0'
    to the find the coordinates property of 'mol0',
    'map1' to find the coordinates property of 'mol1',
    'map2' to find the coordinates property of 'mol2' and
    'map3' to find the coordinates property of 'mol3'
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Angle ImproperID::size(const MoleculeData &mol0data,
                       const PropertyMap &map0,
                       const MoleculeData &mol1data,
                       const PropertyMap &map1,
                       const MoleculeData &mol2data,
                       const PropertyMap &map2,
                       const MoleculeData &mol3data,
                       const PropertyMap &map3) const
{
    return this->torsion(mol0data, map0, 
                         mol1data, map1,
                         mol2data, map2,
                         mol3data, map3).angle();
}

/** Return the size of the improper between atom0() in the 
    molecule whose data is in 'mol0data', atom1() in the 
    molecule whose data is in 'mol1data', atom2() in 
    the molecule whose data is in 'mol2data', and
    atom3() in the molecule whose data is in 'mol3data', 
    using 'map' to find the coordinates property of the 
    molecules
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Angle ImproperID::size(const MoleculeData &mol0data,
                       const MoleculeData &mol1data,
                       const MoleculeData &mol2data,
                       const MoleculeData &mol3data,
                       const PropertyMap &map) const
{
    return this->torsion(mol0data, mol1data, 
                         mol2data, mol3data, map).angle();
}

/** Return the ID of the first atom of the improper */
const AtomID& ImproperID::atom0() const
{
    return atm0.base();
}

/** Return the ID of the second atom of the improper */
const AtomID& ImproperID::atom1() const
{
    return atm1.base();
}

/** Return the ID of the third atom of the improper */
const AtomID& ImproperID::atom2() const
{
    return atm2.base();
}

/** Return the ID of the fourth atom of the improper */
const AtomID& ImproperID::atom3() const
{
    return atm3.base();
}

const char* ImproperID::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ImproperID>() );
}

ImproperID* ImproperID::clone() const
{
    return new ImproperID(*this);
}

