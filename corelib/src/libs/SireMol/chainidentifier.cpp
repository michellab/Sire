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

#include "chainidentifier.h"
#include "chainidx.h"
#include "chainname.h"
#include "molinfo.h"

#include "SireStream/datastream.h"
#include "SireStream/streampolypointer.hpp"

using namespace SireMol;
using namespace SireID;
using namespace SireStream;

////////
//////// Implementation of ChainIdentifier
////////

static const RegisterMetaType<ChainIdentifier> r_chainid;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const ChainIdentifier &chainid)
{
    writeHeader(ds, r_chainid, 1);
    
    SireStream::savePolyPointer(ds, chainid.d);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, ChainIdentifier &chainid)
{
    VersionID v = readHeader(ds, r_chainid);
    
    if (v == 1)
    {
        SireStream::loadPolyPointer(ds, chainid.d);
    }
    else
        throw version_error( v, "1", r_chainid, CODELOC );
        
    return ds;
}

/** Null constructor */
ChainIdentifier::ChainIdentifier() : ChainID()
{}

/** Construct from the passed ChainID */
ChainIdentifier::ChainIdentifier(const ChainID &chainid)
               : ChainID()
{
    if (chainid.isA<ChainIdentifier>())
        d = chainid.asA<ChainIdentifier>().d;
    else if (not chainid.isNull())
        d.reset( chainid.clone() );
}

/** Copy constructor */
ChainIdentifier::ChainIdentifier(const ChainIdentifier &other)
               : ChainID(other), d(other.d)
{}

/** Destructor */
ChainIdentifier::~ChainIdentifier()
{}

/** Is this selection null? */
bool ChainIdentifier::isNull() const
{
    return d.get() == 0;
}

/** Return a hash of this identifier */
uint ChainIdentifier::hash() const
{
    if (d.get() == 0)
        return 0;
    else
        return d->hash();
}
            
/** Return a string representatio of this ID */
QString ChainIdentifier::toString() const
{
    if (d.get() == 0)
        return "null";
    else
        return d->toString();
}

/** Return the base type of this ID */
const ChainID& ChainIdentifier::base() const
{
    if (d.get() == 0)
        return *this;
    else
        return *d;
}

/** Copy assignment operator */
ChainIdentifier& ChainIdentifier::operator=(const ChainIdentifier &other)
{
    d = other.d;
    return *this;
}

/** Copy assignment operator */
ChainIdentifier& ChainIdentifier::operator=(const ChainID &other)
{
    if (other.isA<ChainIdentifier>())
        d = other.asA<ChainIdentifier>().d;
    else if (other.isNull())
        d.reset();
    else
        d.reset(other.clone());
    
    return *this;
}

/** Comparison operator */
bool ChainIdentifier::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<ChainIdentifier>(*this, other);
}

/** Comparison operator */
bool ChainIdentifier::operator==(const ChainIdentifier &other) const
{
    if (d.get() == 0 or other.d.get() == 0)
        return d.get() == other.d.get();
    else
        return d == other.d or *d == *(other.d);
}

/** Comparison operator */
bool ChainIdentifier::operator!=(const ChainIdentifier &other) const
{
    if (d.get() == 0 or other.d.get() == 0)
        return d.get() != other.d.get();
    else
        return d != other.d and *d != *(other.d);
}

/** Comparison operator */
bool ChainIdentifier::operator==(const ChainID &other) const
{
    if (d.get() == 0)
        return other.isNull();
    else if (other.isA<ChainIdentifier>())
        return this->operator==(other.asA<ChainIdentifier>());
    else
        return d->operator==(other);
}

/** Comparison operator */
bool ChainIdentifier::operator!=(const ChainID &other) const
{
    if (d.get() == 0)
        return not other.isNull();
    else if (other.isA<ChainIdentifier>())
        return this->operator!=(other.asA<ChainIdentifier>());
    else
        return d->operator!=(other);
}

/** Map this ID to the list of indicies of chains that match this ID

    \throw SireMol::missing_chain
    \throw SireError::invalid_index
*/
QList<ChainIdx> ChainIdentifier::map(const MolInfo &molinfo) const
{
    if (d.get() == 0)
        return molinfo.getChains();
    else
        return d->map(molinfo);
}

const char* ChainIdentifier::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ChainIdentifier>() );
}

////////
//////// Implementation of ChainIdx
////////

static const RegisterMetaType<ChainIdx> r_chainidx;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const ChainIdx &chainidx)
{
    writeHeader(ds, r_chainidx, 1);
    
    ds << static_cast<const SireID::Index_T_<ChainIdx>&>(chainidx);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, ChainIdx &chainidx)
{
    VersionID v = readHeader(ds, r_chainidx);
    
    if (v == 1)
    {
        ds >> static_cast<SireID::Index_T_<ChainIdx>&>(chainidx);
    }
    else
        throw version_error( v, "1", r_chainidx, CODELOC );
        
    return ds;
}

ChainIdx::ChainIdx() : SireID::Index_T_<ChainIdx>(), ChainID()
{}

ChainIdx::ChainIdx(qint32 idx) : SireID::Index_T_<ChainIdx>(idx), ChainID()
{}

ChainIdx::ChainIdx(const ChainIdx &other) : SireID::Index_T_<ChainIdx>(other), ChainID(other)
{}

ChainIdx::~ChainIdx()
{}

ChainIdx ChainIdx::null()
{
    return ChainIdx();
}

bool ChainIdx::isNull() const
{
    return SireID::Index_T_<ChainIdx>::isNull();
}

uint ChainIdx::hash() const
{
    return SireID::Index_T_<ChainIdx>::hash();
}

QString ChainIdx::toString() const
{
    return QString("ChainIdx(%1)").arg(_idx);
}

ChainIdx& ChainIdx::operator=(const ChainIdx &other)
{
    SireID::IndexBase::operator=(other);
    ChainID::operator=(other);
    return *this;
}

bool ChainIdx::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<ChainIdx>(*this, other);
}

QList<ChainIdx> ChainIdx::map(const MolInfo &molinfo) const
{
    return molinfo.map(*this);
}

const char* ChainIdx::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ChainIdx>() );
}

////////
//////// Implementation of ChainName
////////

static const RegisterMetaType<ChainName> r_chainname;

/** Serialise to a binary datastream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const ChainName &chainname)
{
    writeHeader(ds, r_chainname, 1);
    
    ds << static_cast<const SireID::Name&>(chainname);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, ChainName &chainname)
{
    VersionID v = readHeader(ds, r_chainname);
    
    if (v == 1)
    {
        ds >> static_cast<SireID::Name&>(chainname);
    }
    else
        throw version_error( v, "1", r_chainname, CODELOC );
        
    return ds;
}

ChainName::ChainName() : SireID::Name(), ChainID()
{}

ChainName::ChainName(const QString &name) : SireID::Name(name), ChainID()
{}

ChainName::ChainName(const QString &name, SireID::CaseSensitivity case_sensitivity)
          : SireID::Name(name, case_sensitivity), ChainID()
{}

ChainName::ChainName(const ChainName &other) : SireID::Name(other), ChainID(other)
{}

ChainName::~ChainName()
{}

bool ChainName::isNull() const
{
    return SireID::Name::isNull();
}

uint ChainName::hash() const
{
    return qHash(_name);
}

QString ChainName::toString() const
{
    if (case_sensitive)
        return QString("ChainName('%1')").arg(_name);
    else
        return QString("ChainName('%1', isCaseSensitive=False)").arg(_name);
}

ChainName& ChainName::operator=(const ChainName &other)
{
    SireID::Name::operator=(other);
    ChainID::operator=(other);
    return *this;
}

bool ChainName::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<ChainName>(*this, other);
}

bool ChainName::operator==(const ChainName &other) const
{
    return SireID::Name::operator==(other);
}

bool ChainName::operator!=(const ChainName &other) const
{
    return SireID::Name::operator!=(other);
}

QList<ChainIdx> ChainName::map(const MolInfo &molinfo) const
{
    return molinfo.map(*this);
}

const char* ChainName::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ChainName>() );
}

ChainIdentifier* ChainIdentifier::clone() const
{
    return new ChainIdentifier(*this);
}

ChainName* ChainName::clone() const
{
    return new ChainName(*this);
}


ChainIdx* ChainIdx::clone() const
{
    return new ChainIdx(*this);
}

