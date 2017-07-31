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

#include "segidentifier.h"
#include "segidx.h"
#include "segname.h"
#include "molinfo.h"

#include "SireStream/datastream.h"
#include "SireStream/streampolypointer.hpp"

using namespace SireMol;
using namespace SireID;
using namespace SireStream;

////////
//////// Implementation of SegIdentifier
////////

static const RegisterMetaType<SegIdentifier> r_segid;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const SegIdentifier &segid)
{
    writeHeader(ds, r_segid, 1);
    
    SireStream::savePolyPointer(ds, segid.d);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, SegIdentifier &segid)
{
    VersionID v = readHeader(ds, r_segid);
    
    if (v == 1)
    {
        SireStream::loadPolyPointer(ds, segid.d);
    }
    else
        throw version_error( v, "1", r_segid, CODELOC );
        
    return ds;
}

/** Null constructor */
SegIdentifier::SegIdentifier() : SegID()
{}

/** Construct from the passed SegID */
SegIdentifier::SegIdentifier(const SegID &segid)
              : SegID()
{
    if (segid.isA<SegIdentifier>())
        d = segid.asA<SegIdentifier>().d;
    else if (not segid.isNull())
        d.reset( segid.clone() );
}

/** Copy constructor */
SegIdentifier::SegIdentifier(const SegIdentifier &other)
              : SegID(other), d(other.d)
{}

/** Destructor */
SegIdentifier::~SegIdentifier()
{}

/** Is this selection null? */
bool SegIdentifier::isNull() const
{
    return d.get() == 0;
}

/** Return a hash of this identifier */
uint SegIdentifier::hash() const
{
    if (d.get() == 0)
        return 0;
    else
        return d->hash();
}
            
/** Return a string representatio of this ID */
QString SegIdentifier::toString() const
{
    if (d.get() == 0)
        return "null";
    else
        return d->toString();
}

/** Return the base type of this ID */
const SegID& SegIdentifier::base() const
{
    if (d.get() == 0)
        return *this;
    else
        return *d;
}

/** Copy assignment operator */
SegIdentifier& SegIdentifier::operator=(const SegIdentifier &other)
{
    d = other.d;
    return *this;
}

/** Copy assignment operator */
SegIdentifier& SegIdentifier::operator=(const SegID &other)
{
    if (other.isA<SegIdentifier>())
        d = other.asA<SegIdentifier>().d;
    else if (other.isNull())
        d.reset();
    else
        d.reset(other.clone());
    
    return *this;
}

/** Comparison operator */
bool SegIdentifier::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<SegIdentifier>(*this, other);
}

/** Comparison operator */
bool SegIdentifier::operator==(const SegIdentifier &other) const
{
    if (d.get() == 0 or other.d.get() == 0)
        return d.get() == other.d.get();
    else
        return d == other.d or *d == *(other.d);
}

/** Comparison operator */
bool SegIdentifier::operator!=(const SegIdentifier &other) const
{
    if (d.get() == 0 or other.d.get() == 0)
        return d.get() != other.d.get();
    else
        return d != other.d and *d != *(other.d);
}

/** Comparison operator */
bool SegIdentifier::operator==(const SegID &other) const
{
    if (d.get() == 0)
        return other.isNull();
    else if (other.isA<SegIdentifier>())
        return this->operator==(other.asA<SegIdentifier>());
    else
        return d->operator==(other);
}

/** Comparison operator */
bool SegIdentifier::operator!=(const SegID &other) const
{
    if (d.get() == 0)
        return not other.isNull();
    else if (other.isA<SegIdentifier>())
        return this->operator!=(other.asA<SegIdentifier>());
    else
        return d->operator!=(other);
}

/** Map this ID to the list of indicies of segments that match this ID

    \throw SireMol::missing_segment
    \throw SireError::invalid_index
*/
QList<SegIdx> SegIdentifier::map(const MolInfo &molinfo) const
{
    if (d.get() == 0)
        return molinfo.getSegments();
    else
        return d->map(molinfo);
}

const char* SegIdentifier::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SegIdentifier>() );
}

///////
/////// Implementation of SegIdx
///////

static const RegisterMetaType<SegIdx> r_segidx;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const SegIdx &segidx)
{
    writeHeader(ds, r_segidx, 1);
    
    ds << static_cast<const SireID::Index_T_<SegIdx>&>(segidx);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, SegIdx &segidx)
{
    VersionID v = readHeader(ds, r_segidx);
    
    if (v == 1)
    {
        ds >> static_cast<SireID::Index_T_<SegIdx>&>(segidx);
    }
    else
        throw version_error( v, "1", r_segidx, CODELOC );
        
    return ds;
}

SegIdx::SegIdx() : SireID::Index_T_<SegIdx>(), SegID()
{}

SegIdx::SegIdx(quint32 idx) 
          : SireID::Index_T_<SegIdx>(idx), SegID()
{}

SegIdx::SegIdx(const SegIdx &other) 
          : SireID::Index_T_<SegIdx>(other), SegID(other)
{}

SegIdx::~SegIdx()
{}

SegIdx SegIdx::null()
{
    return SegIdx();
}

bool SegIdx::isNull() const
{
    return SireID::Index_T_<SegIdx>::isNull();
}

uint SegIdx::hash() const
{
    return SireID::Index_T_<SegIdx>::hash();
}

QString SegIdx::toString() const
{
    return QString("SegIdx(%1)").arg(_idx);
}

SegIdx& SegIdx::operator=(const SegIdx &other)
{
    SireID::IndexBase::operator=(other);
    SegID::operator=(other);
    return *this;
}

bool SegIdx::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<SegIdx>(*this, other);
}

QList<SegIdx> SegIdx::map(const MolInfo &molinfo) const
{
    return molinfo.map(*this);
}

const char* SegIdx::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SegIdx>() );
}

///////
/////// Implementation of SegName
///////

static const RegisterMetaType<SegName> r_segname;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const SegName &segname)
{
    writeHeader(ds, r_segname, 1);
    
    ds << static_cast<const SireID::Name&>(segname);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, SegName &segname)
{
    VersionID v = readHeader(ds, r_segname);
    
    if (v == 1)
    {
        ds >> static_cast<SireID::Name&>(segname);
    }
    else
        throw version_error( v, "1", r_segname, CODELOC );
        
    return ds;
}

SegName::SegName() : SireID::Name(), SegID()
{}

SegName::SegName(const QString &name) : SireID::Name(name), SegID()
{}

SegName::SegName(const QString &name, SireID::CaseSensitivity case_sensitivity)
        : SireID::Name(name, case_sensitivity), SegID()
{}

SegName::SegName(const SegName &other) : SireID::Name(other), SegID(other)
{}

SegName::~SegName()
{}

bool SegName::isNull() const
{
    return SireID::Name::isNull();
}

uint SegName::hash() const
{
    return qHash(_name);
}

QString SegName::toString() const
{
    if (case_sensitive)
        return QString("SegName('%1')").arg(_name);
    else
        return QString("SegName('%1', isCaseSensitive=False)").arg(_name);
}

SegName& SegName::operator=(const SegName &other)
{
    SireID::Name::operator=(other);
    SegID::operator=(other);
    return *this;
}

bool SegName::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<SegName>(*this, other);
}

bool SegName::operator==(const SegName &other) const
{
    return SireID::Name::operator==(other);
}

bool SegName::operator!=(const SegName &other) const
{
    return SireID::Name::operator!=(other);
}

QList<SegIdx> SegName::map(const MolInfo &molinfo) const
{
    return molinfo.map(*this);
}

const char* SegName::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SegName>() );
}

SegIdentifier* SegIdentifier::clone() const
{
    return new SegIdentifier(*this);
}

SegIdx* SegIdx::clone() const
{
    return new SegIdx(*this);
}

