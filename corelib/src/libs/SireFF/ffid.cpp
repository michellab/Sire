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

#include "ffid.h"
#include "ffidx.h"
#include "ffname.h"
#include "forcefields.h"

#include "SireFF/errors.h"

#include "SireStream/datastream.h"

using namespace SireFF;
using namespace SireID;
using namespace SireStream;

///////
/////// Implementation of FFID
///////

FFID::FFID() : SireID::ID()
{}

FFID::FFID(const FFID &other) : SireID::ID(other)
{}

FFID::~FFID()
{}

Specify<FFID> FFID::operator[](int i) const
{
    return Specify<FFID>(*this, i);
}

Specify<FFID> FFID::operator()(int i) const
{
    return Specify<FFID>(*this, i);
}

Specify<FFID> FFID::operator()(int i, int j) const
{
    return Specify<FFID>(*this, i, j);
}

IDAndSet<FFID> FFID::operator+(const FFID &other) const
{
    return IDAndSet<FFID>(*this, other);
}

IDAndSet<FFID> FFID::operator&&(const FFID &other) const
{
    return this->operator+(other);
}

IDAndSet<FFID> FFID::operator&(const FFID &other) const
{
    return this->operator+(other);
}

IDOrSet<FFID> FFID::operator*(const FFID &other) const
{
    return IDOrSet<FFID>(*this, other);
}

IDOrSet<FFID> FFID::operator||(const FFID &other) const
{
    return this->operator*(other);
}

IDOrSet<FFID> FFID::operator|(const FFID &other) const
{
    return this->operator*(other);
}

QList<FFIdx> FFID::processMatches(QList<FFIdx> &ffidxs, 
                                  const ForceFields &ffields) const
{
    if (ffidxs.isEmpty())
        throw SireFF::missing_forcefield( QObject::tr(
            "No forcefield in the passed forcefields object matches the ID "
            "\"%1\".")  
                .arg(this->toString()), CODELOC );
                
    qSort(ffidxs);

    return ffidxs;
}

///////
/////// Implementation of FFIdx
///////

static const RegisterMetaType<FFIdx> r_ffidx;

/** Serialise to a binary datastream */
QDataStream SIREFF_EXPORT &operator<<(QDataStream &ds, const FFIdx &ffidx)
{
    writeHeader(ds, r_ffidx, 1);
    
    ds << static_cast<const SireID::Index_T_<FFIdx>&>(ffidx);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREFF_EXPORT &operator>>(QDataStream &ds, FFIdx &ffidx)
{
    VersionID v = readHeader(ds, r_ffidx);
    
    if (v == 1)
    {
        ds >> static_cast<SireID::Index_T_<FFIdx>&>(ffidx);
    }
    else
        throw version_error( v, "1", r_ffidx, CODELOC );
        
    return ds;
}

FFIdx::FFIdx() : SireID::Index_T_<FFIdx>(), FFID()
{}

FFIdx::FFIdx(qint32 idx) : SireID::Index_T_<FFIdx>(idx), FFID()
{}

FFIdx::FFIdx(const FFIdx &other) : SireID::Index_T_<FFIdx>(other), FFID(other)
{}

FFIdx::~FFIdx()
{}

FFIdx FFIdx::null()
{
    return FFIdx();
}

bool FFIdx::isNull() const
{
    return SireID::Index_T_<FFIdx>::isNull();
}

uint FFIdx::hash() const
{
    return SireID::Index_T_<FFIdx>::hash();
}

QString FFIdx::toString() const
{
    return QString("FFIdx(%1)").arg(_idx);
}

FFIdx& FFIdx::operator=(const FFIdx &other)
{
    SireID::IndexBase::operator=(other);
    FFID::operator=(other);
    return *this;
}

bool FFIdx::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<FFIdx>(*this, other);
}

const char* FFIdx::typeName()
{
    return QMetaType::typeName( qMetaTypeId<FFIdx>() );
}

///////
/////// Implementation of FFName
///////

static const RegisterMetaType<FFName> r_ffname;

/** Serialise to a binary datastream */
QDataStream SIREFF_EXPORT &operator<<(QDataStream &ds, const FFName &ffname)
{
    writeHeader(ds, r_ffname, 1);
    
    ds << static_cast<const SireID::Name&>(ffname);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREFF_EXPORT &operator>>(QDataStream &ds, FFName &ffname)
{
    VersionID v = readHeader(ds, r_ffname);
    
    if (v == 1)
    {
        ds >> static_cast<SireID::Name&>(ffname);
    }
    else
        throw version_error( v, "1", r_ffname, CODELOC );
        
    return ds;
}

FFName::FFName() : SireID::Name(), FFID()
{}

FFName::FFName(const QString &name) : SireID::Name(name), FFID()
{}

FFName::FFName(const QString &name, SireID::CaseSensitivity case_sensitivity)
       : SireID::Name(name, case_sensitivity), FFID()
{}

FFName::FFName(const FFName &other) : SireID::Name(other), FFID(other)
{}

FFName::~FFName()
{}

bool FFName::isNull() const
{
    return SireID::Name::isNull();
}

uint FFName::hash() const
{
    return qHash(_name);
}

QString FFName::toString() const
{
    return QString("FFName('%1')").arg(_name);
}

FFName& FFName::operator=(const FFName &other)
{
    SireID::Name::operator=(other);
    FFID::operator=(other);
    return *this;
}

bool FFName::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<FFName>(*this, other);
}

bool FFName::operator==(const FFName &other) const
{
    return SireID::Name::operator==(other);
}

bool FFName::operator!=(const FFName &other) const
{
    return SireID::Name::operator!=(other);
}

const char* FFName::typeName()
{
    return QMetaType::typeName( qMetaTypeId<FFName>() );
}

///////
///////

namespace SireID
{
    template class Specify<FFID>;
    template class IDAndSet<FFID>;
    template class IDOrSet<FFID>;
}

static const RegisterMetaType< Specify<FFID> > r_specify_ffid;
static const RegisterMetaType< IDAndSet<FFID> > r_idandset_ffid;
static const RegisterMetaType< IDOrSet<FFID> > r_idorset_ffid;

FFIdx* FFIdx::clone() const
{
    return new FFIdx(*this);
}


FFName* FFName::clone() const
{
    return new FFName(*this);
}

