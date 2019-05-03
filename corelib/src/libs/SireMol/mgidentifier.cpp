/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#include "mgidentifier.h"
#include "mgnum.h"
#include "mgidx.h"
#include "mgname.h"
#include "moleculegroups.h"

#include "SireStream/datastream.h"
#include "SireStream/streampolypointer.hpp"

using namespace SireMol;
using namespace SireID;
using namespace SireStream;

//////////
////////// Implementation of MGIdentifier
//////////

static const RegisterMetaType<MGIdentifier> r_mgid;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const MGIdentifier &mgid)
{
    writeHeader(ds, r_mgid, 1);
    
    SireStream::savePolyPointer(ds, mgid.d);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, MGIdentifier &mgid)
{
    VersionID v = readHeader(ds, r_mgid);
    
    if (v == 1)
    {
        SireStream::loadPolyPointer(ds, mgid.d);
    }
    else
        throw version_error( v, "1", r_mgid, CODELOC );
        
    return ds;
}

/** Null constructor */
MGIdentifier::MGIdentifier() : MGID()
{}

/** Construct from the passed MGID */
MGIdentifier::MGIdentifier(const MGID &molid)
             : MGID()
{
    if (molid.isA<MGIdentifier>())
        d = molid.asA<MGIdentifier>().d;
    else if (not molid.isNull())
        d.reset( molid.clone() );
}

/** Copy constructor */
MGIdentifier::MGIdentifier(const MGIdentifier &other)
             : MGID(other), d(other.d)
{}

/** Destructor */
MGIdentifier::~MGIdentifier()
{}

/** Is this selection null? */
bool MGIdentifier::isNull() const
{
    return d.get() == 0;
}

/** Return a hash of this identifier */
uint MGIdentifier::hash() const
{
    if (d.get() == 0)
        return 0;
    else
        return d->hash();
}
            
/** Return a string representatio of this ID */
QString MGIdentifier::toString() const
{
    if (d.get() == 0)
        return "null";
    else
        return d->toString();
}

/** Return the base type of this ID */
const MGID& MGIdentifier::base() const
{
    if (d.get() == 0)
        return *this;
    else
        return *d;
}

QList<MGNum> MGIdentifier::map(const MolGroupsBase &molgroups) const
{
    if (d.get() == 0)
        return molgroups.groupNumbers();
    else
        return molgroups.map(*d);
}

/** Copy assignment operator */
MGIdentifier& MGIdentifier::operator=(const MGIdentifier &other)
{
    d = other.d;
    return *this;
}

/** Copy assignment operator */
MGIdentifier& MGIdentifier::operator=(const MGID &other)
{
    if (other.isA<MGIdentifier>())
        d = other.asA<MGIdentifier>().d;
    else if (other.isNull())
        d.reset();
    else
        d.reset(other.clone());
    
    return *this;
}

/** Comparison operator */
bool MGIdentifier::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<MGIdentifier>(*this, other);
}

/** Comparison operator */
bool MGIdentifier::operator==(const MGIdentifier &other) const
{
    if (d.get() == 0 or other.d.get() == 0)
        return d.get() == other.d.get();
    else
        return d == other.d or *d == *(other.d);
}

/** Comparison operator */
bool MGIdentifier::operator!=(const MGIdentifier &other) const
{
    if (d.get() == 0 or other.d.get() == 0)
        return d.get() != other.d.get();
    else
        return d != other.d and *d != *(other.d);
}

/** Comparison operator */
bool MGIdentifier::operator==(const MGID &other) const
{
    if (d.get() == 0)
        return other.isNull();
    else if (other.isA<MGIdentifier>())
        return this->operator==(other.asA<MGIdentifier>());
    else
        return d->operator==(other);
}

/** Comparison operator */
bool MGIdentifier::operator!=(const MGID &other) const
{
    if (d.get() == 0)
        return not other.isNull();
    else if (other.isA<MGIdentifier>())
        return this->operator!=(other.asA<MGIdentifier>());
    else
        return d->operator!=(other);
}

//////////
////////// Implementation of MGIdx
//////////

static const RegisterMetaType<MGIdx> r_mgidx;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const MGIdx &mgidx)
{
    writeHeader(ds, r_mgidx, 1);
    
    ds << static_cast<const SireID::Index_T_<MGIdx>&>(mgidx);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, MGIdx &mgidx)
{
    VersionID v = readHeader(ds, r_mgidx);
    
    if (v == 1)
    {
        ds >> static_cast<SireID::Index_T_<MGIdx>&>(mgidx);
    }
    else
        throw version_error( v, "1", r_mgidx, CODELOC );
        
    return ds;
}

//////////
////////// Implementation of MGNum
//////////

static const RegisterMetaType<MGNum> r_mgnum;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const MGNum &mgnum)
{
    writeHeader(ds, r_mgnum, 1);
    
    ds << static_cast<const SireID::Number&>(mgnum);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, MGNum &mgnum)
{
    VersionID v = readHeader(ds, r_mgnum);
    
    if (v == 1)
    {
        ds >> static_cast<SireID::Number&>(mgnum);
    }
    else
        throw version_error( v, "1", r_mgnum, CODELOC );
        
    return ds;
}

//////////
////////// Implementation of MGName
//////////

static const RegisterMetaType<MGName> r_mgname;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const MGName &mgname)
{
    writeHeader(ds, r_mgname, 1);
    
    ds << static_cast<const SireID::Name&>(mgname);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, MGName &mgname)
{
    VersionID v = readHeader(ds, r_mgname);
    
    if (v == 1)
    {
        ds >> static_cast<SireID::Name&>(mgname);
    }
    else
        throw version_error( v, "1", r_mgname, CODELOC );
        
    return ds;
}

const char* MGIdentifier::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MGIdentifier>() );
}

MGIdentifier* MGIdentifier::clone() const
{
    return new MGIdentifier(*this);
}
