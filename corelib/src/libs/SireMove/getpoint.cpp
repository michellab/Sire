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

#include "getpoint.h"

#include "SireMol/moleculeview.h"
#include "SireMol/evaluator.h"
#include "SireMol/atom.h"
#include "SireMol/atomcoords.h"
#include "SireMol/mover.hpp"
#include "SireMol/selector.hpp"

#include "SireMaths/vector.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireMove;
using namespace SireMol;
using namespace SireMaths;
using namespace SireBase;
using namespace SireStream;

///////////
/////////// Implementation of GetPoint
///////////

static const RegisterMetaType<GetPoint> r_getpoint( MAGIC_ONLY,
                                                    GetPoint::typeName() );

QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds, const GetPoint &getpoint)
{
    writeHeader(ds, r_getpoint, 1);
    
    ds << static_cast<const Property&>(getpoint);
    
    return ds;
}

QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, GetPoint &getpoint)
{
    VersionID v = readHeader(ds, r_getpoint);

    if (v == 1)
    {
        ds >> static_cast<Property&>(getpoint);
    }
    else
        throw version_error(v, "1", r_getpoint, CODELOC);
        
    return ds;
}

/** Constructor */
GetPoint::GetPoint() : Property()
{}

/** Copy constructor */
GetPoint::GetPoint(const GetPoint &other) : Property(other)
{}

/** Destructor */
GetPoint::~GetPoint()
{}

/** Copy assignment operator */
GetPoint& GetPoint::operator=(const GetPoint &other)
{
    Property::operator=(other);
    return *this;
}

/** Comparison operator */
bool GetPoint::operator==(const GetPoint &other) const
{
    return Property::operator==(other);
}

/** Comparison operator */
bool GetPoint::operator!=(const GetPoint &other) const
{
    return Property::operator!=(other);
}

const char* GetPoint::typeName()
{
    return "SireMove::GetPoint";
}

/** Return the point associated with the passed molecule view */
Vector GetPoint::operator()(const MoleculeView &molecule) const
{
    return this->getPoint(molecule, PropertyMap());
}

/** Return the point associated with the passed molecule view,
    using the passed property map to find any necessary properties */
Vector GetPoint::operator()(const MoleculeView &molecule, const PropertyMap &map) const
{
    return this->getPoint(molecule, map);
}

NullGetPoint GetPoint::null()
{
    return NullGetPoint();
}

///////////
/////////// Implementation of NullGetPoint
///////////

static const RegisterMetaType<NullGetPoint> r_nullgetpoint;

QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds,
                                        const NullGetPoint &getpoint)
{
    writeHeader(ds, r_nullgetpoint, 1);
    
    ds << static_cast<const GetPoint&>(getpoint);
    
    return ds;
}

QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, NullGetPoint &getpoint)
{
    VersionID v = readHeader(ds, r_nullgetpoint);
    
    if (v == 1)
    {
        ds >> static_cast<GetPoint&>(getpoint);
    }
    else
        throw version_error(v, "1", r_nullgetpoint, CODELOC);
        
    return ds;
}

/** Constructor */
NullGetPoint::NullGetPoint() : ConcreteProperty<NullGetPoint,GetPoint>()
{}

/** Copy constructor */
NullGetPoint::NullGetPoint(const NullGetPoint &other)
             : ConcreteProperty<NullGetPoint,GetPoint>(other)
{}

/** Destructor */
NullGetPoint::~NullGetPoint()
{}

const char* NullGetPoint::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullGetPoint>() );
}

/** Copy assignment operator */
NullGetPoint& NullGetPoint::operator=(const NullGetPoint &other)
{
    GetPoint::operator=(other);
    return *this;
}

/** Comparison operator */
bool NullGetPoint::operator==(const NullGetPoint &other) const
{
    return GetPoint::operator==(other);
}

/** Comparison operator */
bool NullGetPoint::operator!=(const NullGetPoint &other) const
{
    return GetPoint::operator!=(other);
}

Vector NullGetPoint::getPoint(const MoleculeView &molecule,
                              const PropertyMap &map) const
{
    return Vector(0);
}

///////////
/////////// Implementation of GetCOGPoint
///////////

static const RegisterMetaType<GetCOGPoint> r_getcogpoint;

QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds,
                                        const GetCOGPoint &getpoint)
{
    writeHeader(ds, r_getcogpoint, 1);
    
    SharedDataStream sds(ds);
    
    sds << getpoint.atomids << static_cast<const GetPoint&>(getpoint);
    
    return ds;
}

QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, GetCOGPoint &getpoint)
{
    VersionID v = readHeader(ds, r_getcogpoint);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> getpoint.atomids >> static_cast<GetPoint&>(getpoint);
    }
    else
        throw version_error(v, "1", r_getcogpoint, CODELOC);
        
    return ds;
}

/** Constructor */
GetCOGPoint::GetCOGPoint() : ConcreteProperty<GetCOGPoint,GetPoint>()
{}

/** Construct to get the COG of the atoms in the molecule that match
    the passed AtomID */
GetCOGPoint::GetCOGPoint(const AtomID &atomid)
            : ConcreteProperty<GetCOGPoint,GetPoint>(),
              atomids(atomid)
{}

/** Construct to get the COG of the atoms in the molecule that
    match either of the two passed AtomIDs */
GetCOGPoint::GetCOGPoint(const AtomID &atomid0, const AtomID &atomid1)
            : ConcreteProperty<GetCOGPoint,GetPoint>(),
              atomids(atomid0,atomid1)
{}

/** Construct to get the COG of the atoms in the molecule that
    match any of the passed AtomIDs */
GetCOGPoint::GetCOGPoint(const QList<AtomIdentifier> &ids)
            : ConcreteProperty<GetCOGPoint,GetPoint>(),
              atomids(ids)
{}

/** Copy constructor */
GetCOGPoint::GetCOGPoint(const GetCOGPoint &other)
            : ConcreteProperty<GetCOGPoint,GetPoint>(other),
              atomids(other.atomids)
{}

/** Destructor */
GetCOGPoint::~GetCOGPoint()
{}

const char* GetCOGPoint::typeName()
{
    return QMetaType::typeName( qMetaTypeId<GetCOGPoint>() );
}

/** Copy assignment operator */
GetCOGPoint& GetCOGPoint::operator=(const GetCOGPoint &other)
{
    GetPoint::operator=(other);
    atomids = other.atomids;
    
    return *this;
}

/** Comparison operator */
bool GetCOGPoint::operator==(const GetCOGPoint &other) const
{
    return atomids == other.atomids and GetPoint::operator==(other);
}

/** Comparison operator */
bool GetCOGPoint::operator!=(const GetCOGPoint &other) const
{
    return not GetCOGPoint::operator==(other);
}

/** Return the AtomID(s) used to limit the search for the point */
const AtomID& GetCOGPoint::atomID() const
{
    return atomids;
}

Vector GetCOGPoint::getPoint(const MoleculeView &molecule,
                             const PropertyMap &map) const
{
    if (atomids.IDs().isEmpty())
        return Evaluator(molecule).centerOfGeometry(map);

    else
        return molecule.atoms(atomids,map).evaluate().centerOfGeometry(map);
}

///////////
/////////// Implementation of GetCOMPoint
///////////

static const RegisterMetaType<GetCOMPoint> r_getcompoint;

QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds,
                                        const GetCOMPoint &getpoint)
{
    writeHeader(ds, r_getcompoint, 1);
    
    SharedDataStream sds(ds);
    
    sds << getpoint.atomids << static_cast<const GetPoint&>(getpoint);
    
    return ds;
}

QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, GetCOMPoint &getpoint)
{
    VersionID v = readHeader(ds, r_getcompoint);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> getpoint.atomids >> static_cast<GetPoint&>(getpoint);
    }
    else
        throw version_error(v, "1", r_getcompoint, CODELOC);
        
    return ds;
}

/** Constructor */
GetCOMPoint::GetCOMPoint() : ConcreteProperty<GetCOMPoint,GetPoint>()
{}

/** Construct to get the COM of the atoms in the molecule that match
    the passed AtomID */
GetCOMPoint::GetCOMPoint(const AtomID &atomid)
            : ConcreteProperty<GetCOMPoint,GetPoint>(),
              atomids(atomid)
{}

/** Construct to get the COM of the atoms in the molecule that
    match either of the two passed AtomIDs */
GetCOMPoint::GetCOMPoint(const AtomID &atomid0, const AtomID &atomid1)
            : ConcreteProperty<GetCOMPoint,GetPoint>(),
              atomids(atomid0,atomid1)
{}

/** Construct to get the COM of the atoms in the molecule that
    match any of the passed AtomIDs */
GetCOMPoint::GetCOMPoint(const QList<AtomIdentifier> &ids)
            : ConcreteProperty<GetCOMPoint,GetPoint>(),
              atomids(ids)
{}

/** Copy constructor */
GetCOMPoint::GetCOMPoint(const GetCOMPoint &other)
            : ConcreteProperty<GetCOMPoint,GetPoint>(other),
              atomids(other.atomids)
{}

/** Destructor */
GetCOMPoint::~GetCOMPoint()
{}

const char* GetCOMPoint::typeName()
{
    return QMetaType::typeName( qMetaTypeId<GetCOMPoint>() );
}

/** Copy assignment operator */
GetCOMPoint& GetCOMPoint::operator=(const GetCOMPoint &other)
{
    GetPoint::operator=(other);
    atomids = other.atomids;
    return *this;
}

/** Comparison operator */
bool GetCOMPoint::operator==(const GetCOMPoint &other) const
{
    return atomids == other.atomids and GetPoint::operator==(other);
}

/** Comparison operator */
bool GetCOMPoint::operator!=(const GetCOMPoint &other) const
{
    return not GetCOMPoint::operator==(other);
}

/** Return the AtomID(s) used to limit the search for the point */
const AtomID& GetCOMPoint::atomID() const
{
    return atomids;
}

Vector GetCOMPoint::getPoint(const MoleculeView &molecule,
                             const PropertyMap &map) const
{
    if (atomids.IDs().isEmpty())
        return Evaluator(molecule).centerOfMass(map);
    else
        return molecule.atoms(atomids,map).evaluate().centerOfMass(map);
}
