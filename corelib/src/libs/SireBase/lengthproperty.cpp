/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2014  Christopher Woods
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

#include "SireBase/lengthproperty.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<LengthProperty> r_prop;

QDataStream SIREBASE_EXPORT &operator<<(QDataStream &ds, const LengthProperty &prop)
{
    writeHeader(ds, r_prop, 1);
    ds << prop.val.value();
    return ds;
}

QDataStream SIREBASE_EXPORT &operator>>(QDataStream &ds, LengthProperty &prop)
{
    VersionID v = readHeader(ds, r_prop);
    
    if (v == 1)
    {
        double val;
        ds >> val;
        prop.val = Length(val);
    }
    else
        throw version_error(v, "1", r_prop, CODELOC);
    
    return ds;
}

/** Constructor - this constructs the integer "0" */
LengthProperty::LengthProperty() : ConcreteProperty<LengthProperty,Property>(), val(0)
{}

/** Construct from the passed length */
LengthProperty::LengthProperty(Length value)
               : ConcreteProperty<LengthProperty,Property>(), val(value)
{}

/** Construct from a VariantProperty */
LengthProperty::LengthProperty(const VariantProperty &other)
               : ConcreteProperty<LengthProperty,Property>(other),
                 val(other.convertTo<Length>())
{}

/** Copy constructor */
LengthProperty::LengthProperty(const LengthProperty &other)
               : ConcreteProperty<LengthProperty,Property>(other), val(other.val)
{}

/** Destructor */
LengthProperty::~LengthProperty()
{}

const char* LengthProperty::typeName()
{
    return QMetaType::typeName( qMetaTypeId<LengthProperty>() );
}

/** Copy assignment operator */
LengthProperty& LengthProperty::operator=(const LengthProperty &other)
{
    if (this != &other)
    {
        val = other.val;
    }
    
    return *this;
}

/** Comparison operator */
bool LengthProperty::operator==(const LengthProperty &other) const
{
    return val == other.val;
}

/** Comparison operator */
bool LengthProperty::operator!=(const LengthProperty &other) const
{
    return not this->operator==(other);
}

/** Return this number cast as a double */
Length LengthProperty::value() const
{
    return val;
}

QString LengthProperty::toString() const
{
    return val.toString();
}

/** Return this number cast as a double */
LengthProperty::operator Length() const
{
    return val;
}
