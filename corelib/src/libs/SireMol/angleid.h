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

#ifndef SIREMOL_ANGLEID_H
#define SIREMOL_ANGLEID_H

#include "atomidentifier.h"

#include "SireBase/propertymap.h"
#include "SireUnits/dimensions.h"

#include <boost/tuple/tuple.hpp>

SIRE_BEGIN_HEADER

namespace SireMol
{
class AngleID;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::AngleID&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::AngleID&);

namespace SireMaths
{
class Vector;
class Triangle;
}

namespace SireMol
{

class MoleculeData;
class MoleculeInfoData;
class AtomIdx;

using SireMaths::Vector;
using SireMaths::Triangle;

using SireBase::PropertyMap;

using boost::tuple;

/** This class provides a generic ID for an angle between
    three atoms
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT AngleID : public SireID::ID
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const AngleID&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, AngleID&);

public:
    AngleID();
    AngleID(const AtomID &atom0, const AtomID &atom1,
            const AtomID &atom2);

    AngleID(const AngleID &other);
    
    ~AngleID();
    
    static const char* typeName();
    
    const char* what() const
    {
        return AngleID::typeName();
    }
    
    AngleID* clone() const;
        
    uint hash() const;

    QString toString() const;
    
    bool isNull() const;
    
    AngleID& operator=(const AngleID &other);
    
    bool operator==(const SireID::ID &other) const;
    
    bool operator==(const AngleID &other) const;
    bool operator!=(const AngleID &other) const;
    
    AngleID mirror() const;
    
    tuple<AtomIdx,AtomIdx,AtomIdx> map(const MoleculeInfoData &molinfo) const;
    tuple<AtomIdx,AtomIdx,AtomIdx> map(const MoleculeInfoData &mol0info,
                                       const MoleculeInfoData &mol1info,
                                       const MoleculeInfoData &mol2info) const;
                               
    Vector vector(const MoleculeData &moldata,
                  const PropertyMap &map = PropertyMap()) const;
    
    Vector vector(const MoleculeData &mol0data, 
                  const MoleculeData &mol1data,
                  const MoleculeData &mol2data,
                  const PropertyMap &map = PropertyMap()) const;

    Vector vector(const MoleculeData &mol0data,
                  const PropertyMap &map0,
                  const MoleculeData &mol1data,
                  const PropertyMap &map1,
                  const MoleculeData &mol2data,
                  const PropertyMap &map2) const;

    Triangle triangle(const MoleculeData &moldata,
                      const PropertyMap &map = PropertyMap()) const;
                      
    Triangle triangle(const MoleculeData &mol0data,
                      const MoleculeData &mol1data,
                      const MoleculeData &mol2data,
                      const PropertyMap &map = PropertyMap()) const;
                      
    Triangle triangle(const MoleculeData &mol0data,
                      const PropertyMap &map0,
                      const MoleculeData &mol1data,
                      const PropertyMap &map1,
                      const MoleculeData &mol2data,
                      const PropertyMap &map2) const;

    SireUnits::Dimension::Angle size(const MoleculeData &moldata,
                                     const PropertyMap &map = PropertyMap()) const;

    SireUnits::Dimension::Angle size(const MoleculeData &mol0data, 
                                     const MoleculeData &mol1data,
                                     const MoleculeData &mol2data,
                                     const PropertyMap &map = PropertyMap()) const;
                
    SireUnits::Dimension::Angle size(const MoleculeData &mol0data,
                                     const PropertyMap &map0,
                                     const MoleculeData &mol1data,
                                     const PropertyMap &map1,
                                     const MoleculeData &mol2data,
                                     const PropertyMap &map2) const;

    const AtomID& atom0() const;
    const AtomID& atom1() const;
    const AtomID& atom2() const;

private:
    /** The identifiers of the three atoms */
    AtomIdentifier atm0,atm1,atm2;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

inline uint qHash(const AngleID &angleid)
{
    return angleid.hash();
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE(SireMol::AngleID);

SIRE_EXPOSE_CLASS( SireMol::AngleID )

SIRE_END_HEADER

#endif
