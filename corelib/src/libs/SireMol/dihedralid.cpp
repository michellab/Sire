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

#include "dihedralid.h"

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

static const RegisterMetaType<DihedralID> r_dihedralid;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                       const DihedralID &dihedralid)
{
    writeHeader(ds, r_dihedralid, 1);
    
    SharedDataStream sds(ds);
    
    sds << dihedralid.atm0 << dihedralid.atm1 
        << dihedralid.atm2 << dihedralid.atm3;
         
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                       DihedralID &dihedralid)
{
    VersionID v = readHeader(ds, r_dihedralid);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> dihedralid.atm0 >> dihedralid.atm1 
            >> dihedralid.atm2 >> dihedralid.atm3;
    }
    else
        throw version_error(v, "1", r_dihedralid, CODELOC);
        
    return ds;
}

/** Null constructor */
DihedralID::DihedralID() : ID()
{}

/** Construct a dihedral between the two specified atoms. The order
    is important, as this dihedral may be between two different
    molecules */
DihedralID::DihedralID(const AtomID &atom0, const AtomID &atom1,
                       const AtomID &atom2, const AtomID &atom3)
       : ID(), atm0(atom0), atm1(atom1), atm2(atom2), atm3(atom3)
{}

/** Copy constructor */
DihedralID::DihedralID(const DihedralID &other)
       : ID(other), atm0(other.atm0), atm1(other.atm1),
                    atm2(other.atm2), atm3(other.atm3)
{}

/** Destructor */
DihedralID::~DihedralID()
{}

/** Copy assignment operator */
DihedralID& DihedralID::operator=(const DihedralID &other)
{
    atm0 = other.atm0;
    atm1 = other.atm1;
    atm2 = other.atm2;
    atm3 = other.atm3;
    
    return *this;
}

/** Comparison operator - the order is important */
bool DihedralID::operator==(const DihedralID &other) const
{
    return atm0 == other.atm0 and atm1 == other.atm1 and
           atm2 == other.atm2 and atm3 == other.atm3;
}

/** Comparison operator - the order is important */
bool DihedralID::operator!=(const DihedralID &other) const
{
    return atm0 != other.atm0 or atm1 != other.atm1 or
           atm2 != other.atm2 or atm3 != other.atm3;
}

/** Return the mirror of this DihedralID - i.e. if this is 
    DihedralID(atom0, atom1, atom2, atom3), this returns 
    DihedralID(atom3, atom2, atom1, atom0).
    
    This is useful if you know that DihedralID(atom0,atom1,atom2,atom3) equals
    DihedralID(atom3,atom2,atom1,atom0), e.g. you can now write;
    
    if (not (dihedrals.contains(dihedral) or dihedrals.contains(dihedral.mirror())) )
    {
        dihedrals.insert(dihedral);
    }
    
    or
    
    if (dihedral == other_dihedral or dihedral.mirror() == other.dihedral())
    {
        //this is the same dihedral
    }
*/
DihedralID DihedralID::mirror() const
{
    return DihedralID(atm3, atm2, atm1, atm0);
}

/** Return a hash for this ID */
uint DihedralID::hash() const
{
    return ( (atm0.hash()*atm1.hash()) << 16) | 
           ( (atm2.hash()*atm3.hash()) & 0x0000FFFF);
}

/** Return a string representation of this ID */
QString DihedralID::toString() const
{
    return QString("Dihedral( %1, %2, %3, %4 )")
                .arg(atm0.toString(), atm1.toString(),
                     atm2.toString(), atm3.toString());
}

/** Return whether this is a null ID */
bool DihedralID::isNull() const
{
    return atm0.isNull() and atm1.isNull() and
           atm2.isNull() and atm3.isNull();
}

/** Comparison operator with another ID */
bool DihedralID::operator==(const SireID::ID &other) const
{
    const DihedralID *other_dihedral = dynamic_cast<const DihedralID*>(&other);
    
    return other_dihedral and this->operator==(*other_dihedral);
}

/** Return the indicies of the four atoms in this dihedral - this returns 
    them in the order 
    tuple(dihedral.atom0(),dihedral.atom1(),dihedral.atom2(),dihedral.atom3())
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
tuple<AtomIdx,AtomIdx,AtomIdx,AtomIdx> 
DihedralID::map(const MoleculeInfoData &molinfo) const
{
    return tuple<AtomIdx,AtomIdx,AtomIdx,AtomIdx>( 
                        molinfo.atomIdx(atm0), molinfo.atomIdx(atm1),
                        molinfo.atomIdx(atm2), molinfo.atomIdx(atm3) );
}

/** Return the indicies of the four atoms of this dihedral, between the
    two molecules whose data is in 'mol0info' (containing dihedral.atom0()),
    'mol1info' (containing dihedral.atom1()), 'mol2info' (containing
    dihedral.atom2()) and 'mol3info' (containing dihedral.atom3())
    
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
tuple<AtomIdx,AtomIdx,AtomIdx,AtomIdx>
DihedralID::map(const MoleculeInfoData &mol0info, 
                const MoleculeInfoData &mol1info,
                const MoleculeInfoData &mol2info, 
                const MoleculeInfoData &mol3info) const
{
    return tuple<AtomIdx,AtomIdx,AtomIdx,AtomIdx>( 
                        mol0info.atomIdx(atm0), mol1info.atomIdx(atm1),
                        mol2info.atomIdx(atm2), mol3info.atomIdx(atm3) );
}

/** Return the geometric torsion formed by the four atoms
    of this dihedral in the molecule whose data is in 'moldata',
    using 'map' to find the coordinates property.
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/                      
Torsion DihedralID::torsion(const MoleculeData &moldata,
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
Torsion DihedralID::torsion(const MoleculeData &mol0data,
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
Torsion DihedralID::torsion(const MoleculeData &mol0data,
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
     
/** Return the size of this dihedral in the molecule whose data
    is in 'moldata', using 'map' to find the coordinates property
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Angle DihedralID::size(const MoleculeData &moldata,
                       const PropertyMap &map) const
{
    return this->torsion(moldata,map).angle();
}

/** Return the size of the dihedral between atom0() in the 
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
Angle DihedralID::size(const MoleculeData &mol0data,
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

/** Return the size of the dihedral between atom0() in the 
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
Angle DihedralID::size(const MoleculeData &mol0data,
                       const MoleculeData &mol1data,
                       const MoleculeData &mol2data,
                       const MoleculeData &mol3data,
                       const PropertyMap &map) const
{
    return this->torsion(mol0data, mol1data, 
                         mol2data, mol3data, map).angle();
}

/** Return the ID of the first atom of the dihedral */
const AtomID& DihedralID::atom0() const
{
    return atm0.base();
}

/** Return the ID of the second atom of the dihedral */
const AtomID& DihedralID::atom1() const
{
    return atm1.base();
}

/** Return the ID of the third atom of the dihedral */
const AtomID& DihedralID::atom2() const
{
    return atm2.base();
}

/** Return the ID of the fourth atom of the dihedral */
const AtomID& DihedralID::atom3() const
{
    return atm3.base();
}

const char* DihedralID::typeName()
{
    return QMetaType::typeName( qMetaTypeId<DihedralID>() );
}

DihedralID* DihedralID::clone() const
{
    return new DihedralID(*this);
}

