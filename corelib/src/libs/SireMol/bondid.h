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

#ifndef SIREMOL_BONDID_H
#define SIREMOL_BONDID_H

#include "atomidentifier.h"

#include "SireBase/propertymap.h"

#include <boost/tuple/tuple.hpp>

SIRE_BEGIN_HEADER

namespace SireMol
{
class BondID;
}

QDataStream& operator<<(QDataStream&, const SireMol::BondID&);
QDataStream& operator>>(QDataStream&, SireMol::BondID&);

namespace SireMaths
{
class Vector;
}

namespace SireMol
{

class MoleculeData;
class MoleculeInfoData;
class AtomIdx;

using SireMaths::Vector;

using SireBase::PropertyMap;

using boost::tuple;

/** This class provides a generic ID for a bond between
    two atoms 
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT BondID : public SireID::ID
{

friend QDataStream& ::operator<<(QDataStream&, const BondID&);
friend QDataStream& ::operator>>(QDataStream&, BondID&);

public:
    BondID();
    BondID(const AtomID &atom0, const AtomID &atom1);

    BondID(const BondID &other);
    
    ~BondID();
    
    static const char* typeName();
    
    const char* what() const
    {
        return BondID::typeName();
    }
    
    BondID* clone() const;
    
    uint hash() const;

    QString toString() const;
    
    bool isNull() const;
    
    BondID& operator=(const BondID &other);
    
    bool operator==(const SireID::ID &other) const;
    
    bool operator==(const BondID &other) const;
    bool operator!=(const BondID &other) const;
    
    BondID mirror() const;
    
    tuple<AtomIdx,AtomIdx> map(const MoleculeInfoData &molinfo) const;
    tuple<AtomIdx,AtomIdx> map(const MoleculeInfoData &mol0info,
                               const MoleculeInfoData &mol1info) const;
                               
    Vector vector(const MoleculeData &moldata,
                  const PropertyMap &map = PropertyMap()) const;
    
    Vector vector(const MoleculeData &mol0data, 
                  const MoleculeData &mol1data,
                  const PropertyMap &map = PropertyMap()) const;

    Vector vector(const MoleculeData &mol0data,
                  const PropertyMap &map0,
                  const MoleculeData &mol1data,
                  const PropertyMap &map1) const;

    double size(const MoleculeData &moldata,
                const PropertyMap &map = PropertyMap()) const;

    double size(const MoleculeData &mol0data, 
                const MoleculeData &mol1data,
                const PropertyMap &map = PropertyMap()) const;
                
    double size(const MoleculeData &mol0data,
                const PropertyMap &map0,
                const MoleculeData &mol1data,
                const PropertyMap &map1) const;
                
    double length(const MoleculeData &moldata,
                  const PropertyMap &map = PropertyMap()) const;

    double length(const MoleculeData &mol0data,
                  const MoleculeData &mol1data,
                  const PropertyMap &map = PropertyMap()) const;

    double length(const MoleculeData &mol0data,
                  const PropertyMap &map0,
                  const MoleculeData &mol1data,
                  const PropertyMap &map1) const;

    const AtomID& atom0() const;
    const AtomID& atom1() const;

private:
    /** The identifiers of the two atoms */
    AtomIdentifier atm0,atm1;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

inline uint qHash(const BondID &bondid)
{
    return bondid.hash();
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE(SireMol::BondID);

SIRE_EXPOSE_CLASS( SireMol::BondID )

SIRE_END_HEADER

#endif
