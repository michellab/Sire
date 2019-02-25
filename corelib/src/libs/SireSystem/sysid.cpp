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

#include "sysid.h"
#include "sysidx.h"
#include "sysname.h"

#include "SireSystem/errors.h"

#include "SireStream/datastream.h"

using namespace SireSystem;
using namespace SireStream;

///////
/////// Implementation of SysID
///////

SysID::SysID() : SireID::ID()
{}

SysID::SysID(const SysID &other) : SireID::ID(other)
{}

SysID::~SysID()
{}

Specify<SysID> SysID::operator[](int i) const
{
    return Specify<SysID>(*this, i);
}

Specify<SysID> SysID::operator()(int i) const
{
    return this->operator[](i);
}

Specify<SysID> SysID::operator()(int i, int j) const
{
    return Specify<SysID>(*this, i, j);
}

IDAndSet<SysID> SysID::operator+(const SysID &other) const
{
    return IDAndSet<SysID>(*this, other);
}

IDAndSet<SysID> SysID::operator&&(const SysID &other) const
{
    return this->operator+(other);
}

IDAndSet<SysID> SysID::operator&(const SysID &other) const
{
    return this->operator+(other);
}

IDOrSet<SysID> SysID::operator*(const SysID &other) const
{
    return IDOrSet<SysID>(*this, other);
}

IDOrSet<SysID> SysID::operator||(const SysID &other) const
{
    return this->operator*(other);
}

IDOrSet<SysID> SysID::operator|(const SysID &other) const
{
    return this->operator*(other);
}

QList<SysIdx> SysID::processMatches(QList<SysIdx> &matches,
                                    const Systems &systems) const
{
    if (matches.isEmpty())
        throw SireSystem::missing_system( QObject::tr(
            "There are no systems that match the ID \"%1\".")
                .arg(this->toString()), CODELOC );
                
    qSort(matches);
    
    return matches;
}

///////
/////// Implementation of the Systems stub class
///////

QList<SysIdx> Systems::getSystems() const
{
    return QList<SysIdx>();
}

QList<SysIdx> Systems::map(const SysID&) const
{
    return QList<SysIdx>();
}

///////
/////// Implementation of SysIdx
///////

static const RegisterMetaType<SysIdx> r_sysidx;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const SysIdx &sysidx)
{
    writeHeader(ds, r_sysidx, 1);
    
    ds << static_cast<const SireID::Index_T_<SysIdx>&>(sysidx);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, SysIdx &sysidx)
{
    VersionID v = readHeader(ds, r_sysidx);
    
    if (v == 1)
    {
        ds >> static_cast<SireID::Index_T_<SysIdx>&>(sysidx);
    }
    else
        throw version_error( v, "1", r_sysidx, CODELOC );
        
    return ds;
}

SysIdx::SysIdx() : SireID::Index_T_<SysIdx>(), SysID()
{}

SysIdx::SysIdx(qint32 idx) : SireID::Index_T_<SysIdx>(idx), SysID()
{}

SysIdx::SysIdx(const SysIdx &other) : SireID::Index_T_<SysIdx>(other), SysID(other)
{}

SysIdx::~SysIdx()
{}
  
SysIdx SysIdx::null()
{
    return SysIdx();
}

bool SysIdx::isNull() const
{
    return SireID::Index_T_<SysIdx>::isNull();
}

uint SysIdx::hash() const
{
    return SireID::Index_T_<SysIdx>::hash();
}

QString SysIdx::toString() const
{
    return QString("SysIdx(%1)").arg(_idx);
}

SysIdx& SysIdx::operator=(const SysIdx &other)
{
    SireID::IndexBase::operator=(other);
    SysID::operator=(other);
    return *this;
}

bool SysIdx::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<SysIdx>(*this, other);
}

QList<SysIdx> SysIdx::map(const Systems &systems) const
{
    return systems.map(*this);
}

const char* SysIdx::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SysIdx>() );
}

///////
/////// Implementation of SysName
///////

static const RegisterMetaType<SysName> r_sysname;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const SysName &sysname)
{
    writeHeader(ds, r_sysname, 1);
    
    ds << static_cast<const SireID::Name&>(sysname);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, SysName &sysname)
{
    VersionID v = readHeader(ds, r_sysname);
    
    if (v == 1)
    {
        ds >> static_cast<SireID::Name&>(sysname);
    }
    else
        throw version_error( v, "1", r_sysname, CODELOC );
        
    return ds;
}

SysName::SysName() : SireID::Name(), SysID()
{}

SysName::SysName(const QString &name) : SireID::Name(name), SysID()
{}

SysName::SysName(const SysName &other) : SireID::Name(other), SysID(other)
{}

SysName::~SysName()
{}
    
bool SysName::isNull() const
{
    return SireID::Name::isNull();
}

uint SysName::hash() const
{
    return qHash(_name);
}

QString SysName::toString() const
{
    return QString("SysName('%1')").arg(_name);
}

SysName& SysName::operator=(const SysName &other)
{
    SireID::Name::operator=(other);
    SysID::operator=(other);
    return *this;
}

bool SysName::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<SysName>(*this, other);
}

bool SysName::operator==(const SysName &other) const
{
    return _name == other._name;
}

bool SysName::operator!=(const SysName &other) const
{
    return _name != other._name;
}

QList<SysIdx> SysName::map(const Systems &systems) const
{
    return systems.map(*this);
}

const char* SysName::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SysName>() );
}

///////
///////

namespace SireID
{
    template class Specify<SysID>;
    template class IDAndSet<SysID>;
    template class IDOrSet<SysID>;
}

static const RegisterMetaType< Specify<SysID> > r_specify_sysid;
static const RegisterMetaType< IDAndSet<SysID> > r_idandset_sysid;
static const RegisterMetaType< IDOrSet<SysID> > r_idorset_sysid;

SysIdx* SysIdx::clone() const
{
    return new SysIdx(*this);
}


SysName* SysName::clone() const
{
    return new SysName(*this);
}

