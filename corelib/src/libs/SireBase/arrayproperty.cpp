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

#include "SireBase/propertylist.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireBase;
using namespace SireStream;

//////////
////////// Implementation of DoubleArrayProperty
//////////

static const RegisterMetaType<DoubleArrayProperty> r_doublearray;

QDataStream SIREBASE_EXPORT &operator<<(QDataStream &ds, const DoubleArrayProperty &array)
{
    writeHeader(ds, r_doublearray, 1);
    
    SharedDataStream sds(ds);
    sds << array.array();
    
    return ds;
}

QDataStream SIREBASE_EXPORT &operator>>(QDataStream &ds, DoubleArrayProperty &array)
{
    VersionID v = readHeader(ds, r_doublearray);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        QVector<double> a;
        sds >> a;
        array = DoubleArrayProperty(a);
    }
    else
        throw version_error(v, "1", r_doublearray, CODELOC);
    
    return ds;
}

DoubleArrayProperty::DoubleArrayProperty()
                    : ConcreteProperty<DoubleArrayProperty, ArrayProperty<double> >()
{}

DoubleArrayProperty::DoubleArrayProperty(const QList<double> &array)
                    : ConcreteProperty<DoubleArrayProperty, ArrayProperty<double> >(array)
{}

DoubleArrayProperty::DoubleArrayProperty(const QVector<double> &array)
                    : ConcreteProperty<DoubleArrayProperty, ArrayProperty<double> >(array)
{}

DoubleArrayProperty::DoubleArrayProperty(const DoubleArrayProperty &other)
                    : ConcreteProperty<DoubleArrayProperty, ArrayProperty<double> >(other)
{}

DoubleArrayProperty::~DoubleArrayProperty()
{}

const char* DoubleArrayProperty::typeName()
{
    return QMetaType::typeName(qMetaTypeId<DoubleArrayProperty>());
}

DoubleArrayProperty& DoubleArrayProperty::operator=(const DoubleArrayProperty &other)
{
    ArrayProperty<double>::operator=(other);
    return *this;
}

bool DoubleArrayProperty::operator==(const DoubleArrayProperty &other) const
{
    return ArrayProperty<double>::operator==(other);
}

bool DoubleArrayProperty::operator!=(const DoubleArrayProperty &other) const
{
    return ArrayProperty<double>::operator!=(other);
}

DoubleArrayProperty DoubleArrayProperty::operator+(const DoubleArrayProperty &other) const
{
    return DoubleArrayProperty( ArrayProperty<double>::operator+(other) );
}

DoubleArrayProperty& DoubleArrayProperty::operator+=(const DoubleArrayProperty &other)
{
    ArrayProperty<double>::operator+=(other);
    return *this;
}

//////////
////////// Implementation of IntegerArrayProperty
//////////

static const RegisterMetaType<IntegerArrayProperty> r_intarray;

QDataStream SIREBASE_EXPORT &operator<<(QDataStream &ds, const IntegerArrayProperty &array)
{
    writeHeader(ds, r_intarray, 1);
    
    SharedDataStream sds(ds);
    sds << array.array();
    
    return ds;
}

QDataStream SIREBASE_EXPORT &operator>>(QDataStream &ds, IntegerArrayProperty &array)
{
    VersionID v = readHeader(ds, r_intarray);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        QVector<qint64> a;
        sds >> a;
        array = IntegerArrayProperty(a);
    }
    else
        throw version_error(v, "1", r_intarray, CODELOC);
    
    return ds;
}

IntegerArrayProperty::IntegerArrayProperty()
                     : ConcreteProperty<IntegerArrayProperty, ArrayProperty<qint64> >()
{}

IntegerArrayProperty::IntegerArrayProperty(const QList<qint64> &array)
                    : ConcreteProperty<IntegerArrayProperty, ArrayProperty<qint64> >(array)
{}

IntegerArrayProperty::IntegerArrayProperty(const QVector<qint64> &array)
                    : ConcreteProperty<IntegerArrayProperty, ArrayProperty<qint64> >(array)
{}

IntegerArrayProperty::IntegerArrayProperty(const IntegerArrayProperty &other)
                    : ConcreteProperty<IntegerArrayProperty, ArrayProperty<qint64> >(other)
{}

IntegerArrayProperty::~IntegerArrayProperty()
{}

const char* IntegerArrayProperty::typeName()
{
    return QMetaType::typeName(qMetaTypeId<IntegerArrayProperty>());
}

IntegerArrayProperty& IntegerArrayProperty::operator=(const IntegerArrayProperty &other)
{
    ArrayProperty<qint64>::operator=(other);
    return *this;
}

bool IntegerArrayProperty::operator==(const IntegerArrayProperty &other) const
{
    return ArrayProperty<qint64>::operator==(other);
}

bool IntegerArrayProperty::operator!=(const IntegerArrayProperty &other) const
{
    return ArrayProperty<qint64>::operator!=(other);
}

IntegerArrayProperty IntegerArrayProperty::operator+(const IntegerArrayProperty &other) const
{
    return IntegerArrayProperty( ArrayProperty<qint64>::operator+(other) );
}

IntegerArrayProperty& IntegerArrayProperty::operator+=(const IntegerArrayProperty &other)
{
    ArrayProperty<qint64>::operator+=(other);
    return *this;
}

//////////
////////// Implementation of StringArrayProperty
//////////

static const RegisterMetaType<StringArrayProperty> r_stringarray;

QDataStream SIREBASE_EXPORT &operator<<(QDataStream &ds, const StringArrayProperty &array)
{
    writeHeader(ds, r_stringarray, 1);
    
    SharedDataStream sds(ds);
    sds << array.array();
    
    return ds;
}

QDataStream SIREBASE_EXPORT &operator>>(QDataStream &ds, StringArrayProperty &array)
{
    VersionID v = readHeader(ds, r_stringarray);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        QVector<QString> a;
        sds >> a;
        array = StringArrayProperty(a);
    }
    else
        throw version_error(v, "1", r_stringarray, CODELOC);
    
    return ds;
}

StringArrayProperty::StringArrayProperty()
                    : ConcreteProperty<StringArrayProperty, ArrayProperty<QString> >()
{}

StringArrayProperty::StringArrayProperty(const QList<QString> &array)
                    : ConcreteProperty<StringArrayProperty, ArrayProperty<QString> >(array)
{}

StringArrayProperty::StringArrayProperty(const QVector<QString> &array)
                    : ConcreteProperty<StringArrayProperty, ArrayProperty<QString> >(array)
{}

StringArrayProperty::StringArrayProperty(const StringArrayProperty &other)
                    : ConcreteProperty<StringArrayProperty, ArrayProperty<QString> >(other)
{}

StringArrayProperty::~StringArrayProperty()
{}

const char* StringArrayProperty::typeName()
{
    return QMetaType::typeName(qMetaTypeId<StringArrayProperty>());
}

StringArrayProperty& StringArrayProperty::operator=(const StringArrayProperty &other)
{
    ArrayProperty<QString>::operator=(other);
    return *this;
}

bool StringArrayProperty::operator==(const StringArrayProperty &other) const
{
    return ArrayProperty<QString>::operator==(other);
}

bool StringArrayProperty::operator!=(const StringArrayProperty &other) const
{
    return ArrayProperty<QString>::operator!=(other);
}

StringArrayProperty StringArrayProperty::operator+(const StringArrayProperty &other) const
{
    return StringArrayProperty( ArrayProperty<QString>::operator+(other) );
}

StringArrayProperty& StringArrayProperty::operator+=(const StringArrayProperty &other)
{
    ArrayProperty<QString>::operator+=(other);
    return *this;
}

