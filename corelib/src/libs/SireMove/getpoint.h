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

#ifndef SIREMOVE_GETPOINT_H
#define SIREMOVE_GETPOINT_H

#include "SireBase/property.h"
#include "SireBase/propertymap.h"

#include "SireID/idorset.hpp"

#include "SireMol/atomidx.h"
#include "SireMol/atomidentifier.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class GetPoint;

class NullGetPoint;
class GetCOMPoint;
class GetCOGPoint;
class GetAtomPoint;
class GetIntersectionPoint;
}

QDataStream& operator<<(QDataStream&, const SireMove::GetPoint&);
QDataStream& operator>>(QDataStream&, SireMove::GetPoint&);

QDataStream& operator<<(QDataStream&, const SireMove::NullGetPoint&);
QDataStream& operator>>(QDataStream&, SireMove::NullGetPoint&);

QDataStream& operator<<(QDataStream&, const SireMove::GetCOMPoint&);
QDataStream& operator>>(QDataStream&, SireMove::GetCOMPoint&);

QDataStream& operator<<(QDataStream&, const SireMove::GetCOGPoint&);
QDataStream& operator>>(QDataStream&, SireMove::GetCOGPoint&);

QDataStream& operator<<(QDataStream&, const SireMove::GetAtomPoint&);
QDataStream& operator>>(QDataStream&, SireMove::GetAtomPoint&);

QDataStream& operator<<(QDataStream&, const SireMove::GetIntersectionPoint&);
QDataStream& operator>>(QDataStream&, SireMove::GetIntersectionPoint&);

namespace SireMaths
{
class Vector;
}

namespace SireMol
{
class MoleculeView;
}

namespace SireMove
{

using SireBase::PropertyMap;

using SireMaths::Vector;

using SireMol::AtomID;
using SireMol::MoleculeView;

/** This is the base class of the function objects that are 
    used to return the coordinates of a point in space based
    on the passed PartialMolecule. This is used, for example,
    to find the center of rotation of a set of atoms
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT GetPoint : public SireBase::Property
{

friend QDataStream& ::operator<<(QDataStream&, const GetPoint&);
friend QDataStream& ::operator>>(QDataStream&, GetPoint&);

public:
    GetPoint();
    GetPoint(const GetPoint &other);
    
    virtual ~GetPoint();
    
    static const char* typeName();
    
    virtual GetPoint* clone() const=0;
    
    virtual Vector getPoint(const MoleculeView &molecule,
                            const PropertyMap &map = PropertyMap()) const=0;

    Vector operator()(const MoleculeView &molecule) const;
    Vector operator()(const MoleculeView &molecule, const PropertyMap &map) const;
    
    static NullGetPoint null();
    
protected:
    GetPoint& operator=(const GetPoint &other);
    
    bool operator==(const GetPoint &other) const;
    bool operator!=(const GetPoint &other) const;
};

/** This is the null GetPoint function - this just returns
    the point (0,0,0) for all passed molecule views
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT NullGetPoint 
            : public SireBase::ConcreteProperty<NullGetPoint,GetPoint>
{

friend QDataStream& ::operator<<(QDataStream&, const NullGetPoint&);
friend QDataStream& ::operator>>(QDataStream&, NullGetPoint&);

public:
    NullGetPoint();
    NullGetPoint(const NullGetPoint &other);
    
    ~NullGetPoint();
    
    static const char* typeName();
    
    NullGetPoint& operator=(const NullGetPoint &other);
    
    bool operator==(const NullGetPoint &other) const;
    bool operator!=(const NullGetPoint &other) const;
    
    Vector getPoint(const MoleculeView &molecule,
                    const PropertyMap &map = PropertyMap()) const;
};

/** This function returns the center of geometry (COG) of the 
    atoms in the passed view of the molecule
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT GetCOGPoint 
            : public SireBase::ConcreteProperty<GetCOGPoint,GetPoint>
{

friend QDataStream& ::operator<<(QDataStream&, const GetCOGPoint&);
friend QDataStream& ::operator>>(QDataStream&, GetCOGPoint&);

public:
    GetCOGPoint();

    GetCOGPoint(const AtomID &atomid);
    GetCOGPoint(const AtomID &atomid0, const AtomID &atomid1);
    GetCOGPoint(const QList<SireMol::AtomIdentifier> &atomids);
    
    GetCOGPoint(const GetCOGPoint &other);
    
    ~GetCOGPoint();
    
    static const char* typeName();
    
    GetCOGPoint& operator=(const GetCOGPoint &other);
    
    bool operator==(const GetCOGPoint &other) const;
    bool operator!=(const GetCOGPoint &other) const;

    const AtomID& atomID() const;
    
    Vector getPoint(const MoleculeView &molecule,
                    const PropertyMap &map = PropertyMap()) const;

private:
    /** The list of AtomIDs to use to limit the atoms over which the 
        COG is calculated */
    SireID::IDOrSet<AtomID> atomids;
};

/** This function returns the center of mass (COG) of the 
    atoms in the passed view of the molecule
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT GetCOMPoint 
            : public SireBase::ConcreteProperty<GetCOMPoint,GetPoint>
{

friend QDataStream& ::operator<<(QDataStream&, const GetCOMPoint&);
friend QDataStream& ::operator>>(QDataStream&, GetCOMPoint&);

public:
    GetCOMPoint();

    GetCOMPoint(const AtomID &atomid);
    GetCOMPoint(const AtomID &atomid0, const AtomID &atomid1);
    GetCOMPoint(const QList<SireMol::AtomIdentifier> &atomids);

    GetCOMPoint(const GetCOMPoint &other);
    
    ~GetCOMPoint();
    
    static const char* typeName();
    
    GetCOMPoint& operator=(const GetCOMPoint &other);
    
    bool operator==(const GetCOMPoint &other) const;
    bool operator!=(const GetCOMPoint &other) const;

    const AtomID& atomID() const;
   
    Vector getPoint(const MoleculeView &molecule,
                    const PropertyMap &map = PropertyMap()) const;

private:
    /** The list of AtomIDs to use to limit the atoms over which the 
        COG is calculated */
    SireID::IDOrSet<AtomID> atomids;
};

typedef SireBase::PropPtr<GetPoint> GetPointPtr;

}

Q_DECLARE_METATYPE( SireMove::NullGetPoint )
Q_DECLARE_METATYPE( SireMove::GetCOMPoint )
Q_DECLARE_METATYPE( SireMove::GetCOGPoint )

SIRE_EXPOSE_CLASS( SireMove::GetPoint )
SIRE_EXPOSE_CLASS( SireMove::NullGetPoint )
SIRE_EXPOSE_CLASS( SireMove::GetCOMPoint )
SIRE_EXPOSE_CLASS( SireMove::GetCOGPoint )

SIRE_EXPOSE_PROPERTY( SireMove::GetPointPtr, SireMove::GetPoint )

SIRE_END_HEADER

#endif
