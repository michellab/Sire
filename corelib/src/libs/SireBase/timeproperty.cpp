/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2017  Christopher Woods
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

#include "SireBase/timeproperty.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<TimeProperty> r_prop;

QDataStream SIREBASE_EXPORT &operator<<(QDataStream &ds, const TimeProperty &prop)
{
    writeHeader(ds, r_prop, 1);
    ds << prop.val.value();
    return ds;
}

QDataStream SIREBASE_EXPORT &operator>>(QDataStream &ds, TimeProperty &prop)
{
    VersionID v = readHeader(ds, r_prop);
    
    if (v == 1)
    {
        double val;
        ds >> val;
        prop.val = Time(val);
    }
    else
        throw version_error(v, "1", r_prop, CODELOC);
    
    return ds;
}

/** Constructor - this constructs the integer "0" */
TimeProperty::TimeProperty() : ConcreteProperty<TimeProperty,Property>(), val(0)
{}

/** Construct from the passed length */
TimeProperty::TimeProperty(Time value)
               : ConcreteProperty<TimeProperty,Property>(), val(value)
{}

/** Construct from a VariantProperty */
TimeProperty::TimeProperty(const VariantProperty &other)
               : ConcreteProperty<TimeProperty,Property>(other),
                 val(other.convertTo<Time>())
{}

/** Copy constructor */
TimeProperty::TimeProperty(const TimeProperty &other)
               : ConcreteProperty<TimeProperty,Property>(other), val(other.val)
{}

/** Destructor */
TimeProperty::~TimeProperty()
{}

const char* TimeProperty::typeName()
{
    return QMetaType::typeName( qMetaTypeId<TimeProperty>() );
}

/** Copy assignment operator */
TimeProperty& TimeProperty::operator=(const TimeProperty &other)
{
    if (this != &other)
    {
        val = other.val;
    }
    
    return *this;
}

/** Comparison operator */
bool TimeProperty::operator==(const TimeProperty &other) const
{
    return val == other.val;
}

/** Comparison operator */
bool TimeProperty::operator!=(const TimeProperty &other) const
{
    return not this->operator==(other);
}

/** Return this number cast as a Time */
Time TimeProperty::value() const
{
    return val;
}

QString TimeProperty::toString() const
{
    return val.toString();
}

/** Return this number cast as a double */
TimeProperty::operator Time() const
{
    return val;
}
