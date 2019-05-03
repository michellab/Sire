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

#include "ffidentifier.h"
#include "forcefields.h"
#include "ffidx.h"
#include "ffname.h"

#include "SireStream/datastream.h"
#include "SireStream/streampolypointer.hpp"

using namespace SireFF;
using namespace SireID;
using namespace SireStream;

static const RegisterMetaType<FFIdentifier> r_ffid;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const FFIdentifier &ffid)
{
    writeHeader(ds, r_ffid, 1);
    
    SireStream::savePolyPointer(ds, ffid.d);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, FFIdentifier &ffid)
{
    VersionID v = readHeader(ds, r_ffid);
    
    if (v == 1)
    {
        SireStream::loadPolyPointer(ds, ffid.d);
    }
    else
        throw version_error( v, "1", r_ffid, CODELOC );
        
    return ds;
}

/** Null constructor */
FFIdentifier::FFIdentifier() : FFID()
{}

/** Construct from the passed FFID */
FFIdentifier::FFIdentifier(const FFID &ffid)
             : FFID()
{
    if (ffid.isA<FFIdentifier>())
        d = ffid.asA<FFIdentifier>().d;
    else if (not ffid.isNull())
        d.reset( ffid.clone() );
}

/** Copy constructor */
FFIdentifier::FFIdentifier(const FFIdentifier &other)
             : FFID(other), d(other.d)
{}

/** Destructor */
FFIdentifier::~FFIdentifier()
{}

/** Is this selection null? */
bool FFIdentifier::isNull() const
{
    return d.get() == 0;
}

/** Return a hash of this identifier */
uint FFIdentifier::hash() const
{
    if (d.get() == 0)
        return 0;
    else
        return d->hash();
}
            
/** Return a string representatio of this ID */
QString FFIdentifier::toString() const
{
    if (d.get() == 0)
        return "null";
    else
        return d->toString();
}

/** Return the base type of this ID */
const FFID& FFIdentifier::base() const
{
    if (d.get() == 0)
        return *this;
    else
        return *d;
}

/** Copy assignment operator */
FFIdentifier& FFIdentifier::operator=(const FFIdentifier &other)
{
    d = other.d;
    return *this;
}

/** Copy assignment operator */
FFIdentifier& FFIdentifier::operator=(const FFID &other)
{
    if (other.isA<FFIdentifier>())
        d = other.asA<FFIdentifier>().d;
    else if (other.isNull())
        d.reset();
    else
        d.reset(other.clone());
    
    return *this;
}

/** Comparison operator */
bool FFIdentifier::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<FFIdentifier>(*this, other);
}

/** Comparison operator */
bool FFIdentifier::operator==(const FFIdentifier &other) const
{
    if (d.get() == 0 or other.d.get() == 0)
        return d.get() == other.d.get();
    else
        return d == other.d or *d == *(other.d);
}

/** Comparison operator */
bool FFIdentifier::operator!=(const FFIdentifier &other) const
{
    if (d.get() == 0 or other.d.get() == 0)
        return d.get() != other.d.get();
    else
        return d != other.d and *d != *(other.d);
}

/** Comparison operator */
bool FFIdentifier::operator==(const FFID &other) const
{
    if (d.get() == 0)
        return other.isNull();
    else if (other.isA<FFIdentifier>())
        return this->operator==(other.asA<FFIdentifier>());
    else
        return d->operator==(other);
}

/** Comparison operator */
bool FFIdentifier::operator!=(const FFID &other) const
{
    if (d.get() == 0)
        return not other.isNull();
    else if (other.isA<FFIdentifier>())
        return this->operator!=(other.asA<FFIdentifier>());
    else
        return d->operator!=(other);
}

static QList<FFIdx> allFFIdxs(int n)
{
    QList<FFIdx> ffidxs;
    
    for (int i=0; i<n; ++i)
    {
        ffidxs.append( FFIdx(i) );
    }
    
    return ffidxs;
}

/** Map this name to the index of the matching forcefield in the 
    passed ForceFields object 
    
    \throw SireFF::missing_forcefield
*/
QList<FFIdx> FFName::map(const ForceFields &ffields) const
{
    if (this->isNull())
        return ::allFFIdxs(ffields.nForceFields());

    return ffields.map(*this);
}

/** Short cut function to map this index to the index of the
    matching forcefield in the passed ForceFields object
    
    \throw SireError::invalid_index
*/
QList<FFIdx> FFIdx::map(const ForceFields &ffields) const
{
    if (this->isNull())
        return ::allFFIdxs(ffields.nForceFields());

    return ffields.map(*this);
}

/** Map this identifier to the indicies of the matching
    forcefields in the passed ForceFields object 
    
    \throw SireFF::missing_forcefield
    \throw SireError::invalid_index
*/
QList<FFIdx> FFIdentifier::map(const ForceFields &ffields) const
{
    if (d.get() == 0)
    {
        return ::allFFIdxs(ffields.nForceFields());
    }
    else
        return d->map(ffields);
}

const char* FFIdentifier::typeName()
{
    return QMetaType::typeName( qMetaTypeId<FFIdentifier>() );
}

FFIdentifier* FFIdentifier::clone() const
{
    return new FFIdentifier(*this);
}

