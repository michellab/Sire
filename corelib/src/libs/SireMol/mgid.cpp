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

#include "mgid.h"
#include "mgidx.h"
#include "mgname.h"
#include "mgnum.h"
#include "mgidentifier.h"
#include "moleculegroups.h"

#include "SireBase/incremint.h"

#include "SireMol/errors.h"
#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "tostring.h"

using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

////////
//////// Implementation of MGID
////////

MGID::MGID() : SireID::ID()
{}

MGID::MGID(const MGID &other) : SireID::ID(other)
{}

MGID::~MGID()
{}

Specify<MGID> MGID::operator[](int i) const
{
    return Specify<MGID>(*this, i);
}

Specify<MGID> MGID::operator()(int i) const
{
    return this->operator[](i);
}

Specify<MGID> MGID::operator()(int i, int j) const
{
    return Specify<MGID>(*this, i, j);
}

IDAndSet<MGID> MGID::operator+(const MGID &other) const
{
    return IDAndSet<MGID>(*this, other);
}

IDAndSet<MGID> MGID::operator&&(const MGID &other) const
{
    return this->operator+(other);
}

IDAndSet<MGID> MGID::operator&(const MGID &other) const
{
    return this->operator+(other);
}

IDOrSet<MGID> MGID::operator*(const MGID &other) const
{
    return IDOrSet<MGID>(*this, other);
}

IDOrSet<MGID> MGID::operator||(const MGID &other) const
{
    return this->operator*(other);
}

IDOrSet<MGID> MGID::operator|(const MGID &other) const
{
    return this->operator*(other);
}

void MGID::processMatches(QList<MGNum> &matches, const MolGroupsBase&) const
{
    if (matches.isEmpty())
        throw SireMol::missing_group( QObject::tr(
                "There is no group in the passed groups that matches "
                "the ID \"%1\".")
                    .arg(this->toString()), CODELOC );
}

////////
//////// Implementation of MGIdx
////////

MGIdx::MGIdx() : SireID::Index_T_<MGIdx>(), MGID()
{}

MGIdx::MGIdx(qint32 idx) : SireID::Index_T_<MGIdx>(idx), MGID()
{}

MGIdx::MGIdx(const MGIdx &other) : SireID::Index_T_<MGIdx>(other), MGID(other)
{}

MGIdx::~MGIdx()
{}

MGIdx MGIdx::null()
{
    return MGIdx();
}

bool MGIdx::isNull() const
{
    return SireID::Index_T_<MGIdx>::isNull();
}

uint MGIdx::hash() const
{
    return SireID::Index_T_<MGIdx>::hash();
}

QString MGIdx::toString() const
{
    return QString("MGIdx(%1)").arg(_idx);
}

MGIdx& MGIdx::operator=(const MGIdx &other)
{
    SireID::IndexBase::operator=(other);
    MGID::operator=(other);
    return *this;
}

bool MGIdx::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<MGIdx>(*this, other);
}

QList<MGNum> MGIdx::map(const MolGroupsBase &molgroups) const
{
    return molgroups.map(*this);
}

const char* MGIdx::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MGIdx>() );
}

////////
//////// Implementation of MGName
////////

MGName::MGName() : SireID::Name(), MGID()
{}

MGName::MGName(const QString &name) : SireID::Name(name), MGID()
{}

MGName::MGName(const QString &name, SireID::CaseSensitivity case_sensitivity)
       : SireID::Name(name, case_sensitivity), MGID()
{}

MGName::MGName(const MGName &other) : SireID::Name(other), MGID(other)
{}

MGName::~MGName()
{}

bool MGName::isNull() const
{
    return SireID::Name::isNull();
}

uint MGName::hash() const
{
    return qHash(_name);
}

QString MGName::toString() const
{
    if (case_sensitive)
        return QString("MGName('%1')").arg(_name);
    else
        return QString("MGName('%1', isCaseSensitive=False)").arg(_name);
}

MGName& MGName::operator=(const MGName &other)
{
    SireID::Name::operator=(other);
    MGID::operator=(other);
    return *this;
}

bool MGName::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<MGName>(*this, other);
}

bool MGName::operator==(const MGName &other) const
{
    return SireID::Name::operator==(other);
}

bool MGName::operator!=(const MGName &other) const
{
    return SireID::Name::operator!=(other);
}

QList<MGNum> MGName::map(const MolGroupsBase &molgroups) const
{
    return molgroups.map(*this);
}

const char* MGName::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MGName>() );
}

////////
//////// Implementation of MGNum
////////

MGNum::MGNum() : SireID::Number(), MGID()
{}

MGNum::MGNum(quint32 num) : SireID::Number(num), MGID()
{}

MGNum::MGNum(const MGNum &other) : SireID::Number(other), MGID(other)
{}

MGNum::~MGNum()
{}
    
bool MGNum::isNull() const
{
    return SireID::Number::isNull();
}

uint MGNum::hash() const
{
    return ::qHash( static_cast<const SireID::Number&>(*this) );
}

QString MGNum::toString() const
{
    return QString("MGNum(%1)").arg(_num);
}

MGNum& MGNum::operator=(const MGNum &other)
{
    SireID::Number::operator=(other);
    MGID::operator=(other);
    return *this;
}

bool MGNum::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<MGNum>(*this, other);
}

bool MGNum::operator==(const MGNum &other) const
{
    return _num == other._num;
}

bool MGNum::operator!=(const MGNum &other) const
{
    return _num != other._num;
}

bool MGNum::operator<(const MGNum &other) const
{
    return _num < other._num;
}

bool MGNum::operator<=(const MGNum &other) const
{
    return _num <= other._num;
}

bool MGNum::operator>(const MGNum &other) const
{
    return _num > other._num;
}

bool MGNum::operator>=(const MGNum &other) const
{
    return _num >= other._num;
}

QList<MGNum> MGNum::map(const MolGroupsBase &molgroups) const
{
    return molgroups.map(*this);
}

const char* MGNum::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MGNum>() );
}

///////

namespace SireID
{
    template class IDAndSet<MGID>;
    template class IDOrSet<MGID>;
    template class Specify<MGID>;
}

static const RegisterMetaType< IDAndSet<MGID> > r_idandset_mgid;
static const RegisterMetaType< IDOrSet<MGID> > r_idorset_mgid;
static const RegisterMetaType< Specify<MGID> > r_specify_mgid;


MGNum* MGNum::clone() const
{
    return new MGNum(*this);
}


MGIdx* MGIdx::clone() const
{
    return new MGIdx(*this);
}


MGName* MGName::clone() const
{
    return new MGName(*this);
}

