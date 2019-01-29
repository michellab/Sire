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

#include "id.h"
#include "name.h"
#include "number.h"

#include <QDebug>

using namespace SireID;

///////////
/////////// Implementation of ID
///////////

ID::ID()
{}

ID::ID(const ID&)
{}

ID::~ID()
{}

///////////
/////////// Implementation of Name
///////////

/** Serialise a Name class */
QDataStream &operator<<(QDataStream &ds, const SireID::Name &name)
{
    SireStream::SharedDataStream sds(ds);
    sds << name._name << name.case_sensitive;

    return ds;
}

/** Deserialise a Name class */
QDataStream &operator>>(QDataStream &ds, SireID::Name &name)
{
    SireStream::SharedDataStream sds(ds);
    sds >> name._name >> name.case_sensitive;
    
    return ds;
}

Name::Name(const QString &name, CaseSensitivity case_sensitivity) 
     : _name(name)
{
    switch (case_sensitivity)
    {
        case CaseSensitive:
            case_sensitive = true;
            break;
        case CaseInsensitive:
            case_sensitive = false;
            break;
    }
}

Name::Name(const Name &other) 
     : _name(other._name), case_sensitive(other.case_sensitive)
{}

Name::~Name()
{}

Name& Name::operator=(const Name &other)
{
    _name = other._name;
    case_sensitive = other.case_sensitive;
    
    return *this;
}

bool Name::operator==(const Name &other) const
{
    return _name == other._name and case_sensitive == other.case_sensitive;
}

bool Name::operator!=(const Name &other) const
{
    return _name != other._name or case_sensitive != other.case_sensitive;
}

Name::operator QString() const
{
    return _name;
}

bool Name::isNull() const
{
    return _name.isNull();
}

bool Name::isCaseSensitive() const
{
    return case_sensitive;
}

uint Name::hash() const
{
    return ::qHash(_name);
}

bool Name::isEmpty() const
{
    return _name.isEmpty();
}

const QString& Name::value() const
{
    return _name;
}

///////////
/////////// Implementation of Number
///////////


/** Serialise a Number class */
QDataStream &operator<<(QDataStream &ds, const SireID::Number &number)
{
    ds << number._num;
    return ds;
}

/** Deserialise a Number class */
QDataStream &operator>>(QDataStream &ds, SireID::Number &number)
{
    ds >> number._num;
    return ds;
}

Number::Number(qint32 num) : _num(num)
{}

Number::Number(const Number &other) : _num(other._num)
{}

Number::~Number()
{}

Number::operator qint32() const
{
    return _num;
}

qint32 Number::null()
{
    return std::numeric_limits<qint32>::min();
}

uint Number::hash() const
{
    return quint32(_num);
}

bool Number::isNull() const
{
    return _num == null();
}

qint32 Number::value() const
{
    return _num;
}
