/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2013  Christopher Woods
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

#include "SireBase/numberproperty.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<NumberProperty> r_numprop;

QDataStream SIREBASE_EXPORT &operator<<(QDataStream &ds, const NumberProperty &numprop)
{
    writeHeader(ds, r_numprop, 1);
    
    ds << numprop.is_int;
    
    if (numprop.is_int)
        ds << numprop.val.ival;
    else
        ds << numprop.val.dval;
    
    return ds;
}

QDataStream SIREBASE_EXPORT &operator>>(QDataStream &ds, NumberProperty &numprop)
{
    VersionID v = readHeader(ds, r_numprop);
    
    if (v == 1)
    {
        ds >> numprop.is_int;
        
        if (numprop.is_int)
            ds >> numprop.val.ival;
        else
            ds >> numprop.val.dval;
    }
    else
        throw version_error(v, "1", r_numprop, CODELOC);
    
    return ds;
}

/** Constructor - this constructs the integer "0" */
NumberProperty::NumberProperty() : ConcreteProperty<NumberProperty,Property>(),
                                   is_int(true)
{
    val.ival = 0;
}

/** Construct from the passed double */
NumberProperty::NumberProperty(double value)
               : ConcreteProperty<NumberProperty,Property>(),
                 is_int(false)
{
    val.dval = value;
}

/** Construct from the passed integer */
NumberProperty::NumberProperty(int value)
               : ConcreteProperty<NumberProperty,Property>(),
                 is_int(true)
{
    val.ival = value;
}

/** Construct from the passed integer */
NumberProperty::NumberProperty(qint64 value)
               : ConcreteProperty<NumberProperty,Property>(),
                 is_int(true)
{
    val.ival = value;
}

/** Copy constructor */
NumberProperty::NumberProperty(const NumberProperty &other)
               : ConcreteProperty<NumberProperty,Property>(other),
                 is_int(other.is_int)
{
    if (is_int)
        val.ival = other.val.ival;
    else
        val.dval = other.val.dval;
}

/** Destructor */
NumberProperty::~NumberProperty()
{}

const char* NumberProperty::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NumberProperty>() );
}

/** Copy assignment operator */
NumberProperty& NumberProperty::operator=(const NumberProperty &other)
{
    if (this != &other)
    {
        is_int = other.is_int;
        
        if (is_int)
            val.ival = other.val.ival;
        else
            val.dval = other.val.dval;
    }
    
    return *this;
}

/** Comparison operator */
bool NumberProperty::operator==(const NumberProperty &other) const
{
    return this->value() == other.value();
}

/** Comparison operator */
bool NumberProperty::operator!=(const NumberProperty &other) const
{
    return not this->operator==(other);
}

/** Return this number cast as an integer */
qint64 NumberProperty::toInt() const
{
    if (is_int)
        return val.ival;
    else
        return qint64( val.dval );
}

/** Return this number cast as a double */
double NumberProperty::toDouble() const
{
    if (is_int)
        return double(val.ival);
    else
        return val.dval;
}

/** Return this number cast as a double */
double NumberProperty::value() const
{
    return this->toDouble();
}

QString NumberProperty::toString() const
{
    if (is_int)
        return QString::number(val.ival);
    else
        return QString::number(val.dval);
}

/** Return this number cast as a double */
NumberProperty::operator double() const
{
    return this->toDouble();
}

/** Return this number cast as an int */
NumberProperty::operator int() const
{
    return this->toInt();
}

/** Return this number cast as an int */
NumberProperty::operator qint64() const
{
    return this->toInt();
}
