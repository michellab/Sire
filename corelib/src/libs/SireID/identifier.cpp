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

#include "identifier.h"

#include "SireStream/datastream.h"
#include "SireStream/streampolypointer.hpp"

using namespace SireID;
using namespace SireStream;

static const RegisterMetaType<Identifier> r_id;

/** Serialise to a binary datastream */
QDataStream SIREID_EXPORT &operator<<(QDataStream &ds, const Identifier &id)
{
    writeHeader(ds, r_id, 1);
    
    SireStream::savePolyPointer(ds, id.d);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREID_EXPORT &operator>>(QDataStream &ds, Identifier &id)
{
    VersionID v = readHeader(ds, r_id);
    
    if (v == 1)
    {
        SireStream::loadPolyPointer(ds, id.d);
    }
    else
        throw version_error( v, "1", r_id, CODELOC );
        
    return ds;
}

/** Null constructor */
Identifier::Identifier() : ID()
{}

/** Construct from the passed ID */
Identifier::Identifier(const ID &id) : ID()
{
    if (id.isA<Identifier>())
        d = id.asA<Identifier>().d;
    else if (not id.isNull())
        d.reset( id.clone() );
}

/** Copy constructor */
Identifier::Identifier(const Identifier &other)
           : ID(other), d(other.d)
{}

/** Destructor */
Identifier::~Identifier()
{}

/** Is this ID null? */
bool Identifier::isNull() const
{
    return d.get() == 0;
}

/** Return a hash of this identifier */
uint Identifier::hash() const
{
    if (d.get() == 0)
        return 0;
    else
        return d->hash();
}
            
/** Return a string representatio of this ID */
QString Identifier::toString() const
{
    if (d.get() == 0)
        return "null";
    else
        return d->toString();
}

/** Return the base type of this ID */
const ID& Identifier::base() const
{
    if (d.get() == 0)
        return *this;
    else
        return *d;
}

/** Copy assignment operator */
Identifier& Identifier::operator=(const Identifier &other)
{
    d = other.d;
    return *this;
}

/** Copy assignment operator */
Identifier& Identifier::operator=(const ID &other)
{
    if (other.isA<Identifier>())
        d = other.asA<Identifier>().d;
    else if (other.isNull())
        d.reset();
    else
        d.reset(other.clone());
    
    return *this;
}

/** Comparison operator */
bool Identifier::operator==(const ID &other) const
{
    return ID::compare<Identifier>(*this, other);
}

/** Comparison operator */
bool Identifier::operator!=(const ID &other) const
{
    return ID::operator!=(other);
}

/** Comparison operator */
bool Identifier::operator==(const Identifier &other) const
{
    if (d.get() == 0 or other.d.get() == 0)
        return d.get() == other.d.get();
    else
        return d == other.d or *d == *(other.d);
}

/** Comparison operator */
bool Identifier::operator!=(const Identifier &other) const
{
    if (d.get() == 0 or other.d.get() == 0)
        return d.get() != other.d.get();
    else
        return d != other.d and *d != *(other.d);
}

const char* Identifier::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Identifier>() );
}

Identifier* Identifier::clone() const
{
    return new Identifier(*this);
}

