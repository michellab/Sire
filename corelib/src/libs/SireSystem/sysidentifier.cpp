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

#include "sysidentifier.h"
#include "sysidx.h"
#include "sysname.h"

#include "SireStream/datastream.h"
#include "SireStream/streampolypointer.hpp"

using namespace SireSystem;
using namespace SireID;
using namespace SireStream;

static const RegisterMetaType<SysIdentifier> r_sysid;

/** Serialise to a binary datastream */
QDataStream SIRESYSTEM_EXPORT &operator<<(QDataStream &ds, 
                                          const SysIdentifier &sysid)
{
    writeHeader(ds, r_sysid, 1);
    
    SireStream::savePolyPointer(ds, sysid.d);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIRESYSTEM_EXPORT &operator>>(QDataStream &ds, SysIdentifier &sysid)
{
    VersionID v = readHeader(ds, r_sysid);
    
    if (v == 1)
    {
        SireStream::loadPolyPointer(ds, sysid.d);
    }
    else
        throw version_error( v, "1", r_sysid, CODELOC );
        
    return ds;
}

/** Null constructor */
SysIdentifier::SysIdentifier() : SysID()
{}

/** Construct from the passed SysID */
SysIdentifier::SysIdentifier(const SysID &sysid)
              : SysID()
{
    if (sysid.isA<SysIdentifier>())
        d = sysid.asA<SysIdentifier>().d;
    else if (not sysid.isNull())
        d.reset( sysid.clone() );
}

/** Copy constructor */
SysIdentifier::SysIdentifier(const SysIdentifier &other)
              : SysID(other), d(other.d)
{}

/** Destructor */
SysIdentifier::~SysIdentifier()
{}

/** Is this selection null? */
bool SysIdentifier::isNull() const
{
    return d.get() == 0;
}

/** Return a hash of this identifier */
uint SysIdentifier::hash() const
{
    if (d.get() == 0)
        return 0;
    else
        return d->hash();
}
            
/** Return a string representatio of this ID */
QString SysIdentifier::toString() const
{
    if (d.get() == 0)
        return "null";
    else
        return d->toString();
}

/** Return the base type of this ID */
const SysID& SysIdentifier::base() const
{
    if (d.get() == 0)
        return *this;
    else
        return *d;
}

/** Copy assignment operator */
SysIdentifier& SysIdentifier::operator=(const SysIdentifier &other)
{
    d = other.d;
    return *this;
}

/** Copy assignment operator */
SysIdentifier& SysIdentifier::operator=(const SysID &other)
{
    if (other.isA<SysIdentifier>())
        d = other.asA<SysIdentifier>().d;
    else if (other.isNull())
        d.reset();
    else
        d.reset(other.clone());
    
    return *this;
}

/** Comparison operator */
bool SysIdentifier::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<SysIdentifier>(*this, other);
}

/** Comparison operator */
bool SysIdentifier::operator==(const SysIdentifier &other) const
{
    if (d.get() == 0 or other.d.get() == 0)
        return d.get() == other.d.get();
    else
        return d == other.d or *d == *(other.d);
}

/** Comparison operator */
bool SysIdentifier::operator!=(const SysIdentifier &other) const
{
    if (d.get() == 0 or other.d.get() == 0)
        return d.get() != other.d.get();
    else
        return d != other.d and *d != *(other.d);
}

/** Comparison operator */
bool SysIdentifier::operator==(const SysID &other) const
{
    if (d.get() == 0)
        return other.isNull();
    else if (other.isA<SysIdentifier>())
        return this->operator==(other.asA<SysIdentifier>());
    else
        return d->operator==(other);
}

/** Comparison operator */
bool SysIdentifier::operator!=(const SysID &other) const
{
    if (d.get() == 0)
        return not other.isNull();
    else if (other.isA<SysIdentifier>())
        return this->operator!=(other.asA<SysIdentifier>());
    else
        return d->operator!=(other);
}

/** Map this ID */
QList<SysIdx> SysIdentifier::map(const Systems &systems) const
{
    if (d.get() == 0)
        return systems.getSystems();
    else
        return d->map(systems);
}

const char* SysIdentifier::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SysIdentifier>() );
}

SysIdentifier* SysIdentifier::clone() const
{
    return new SysIdentifier(*this);
}
