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

#include "SireMaths/vectorproperty.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireBase;
using namespace SireMaths;
using namespace SireStream;

namespace SireMaths
{
    PropertyPtr SIREMATHS_EXPORT wrap(const Vector &vector)
    {
        return PropertyPtr( VectorProperty(vector) );
    }
    
    PropertyPtr SIREMATHS_EXPORT wrap(const QVector<Vector> &vector)
    {
        return PropertyPtr( VectorArrayProperty(vector) );
    }
    
    PropertyPtr SIREMATHS_EXPORT wrap(const QList<Vector> &vector)
    {
        return PropertyPtr( VectorArrayProperty(vector) );
    }
}

//////////
////////// Implementation of VectorArrayProperty
//////////

static const RegisterMetaType<VectorArrayProperty> r_vectorarray;

QDataStream SIREMATHS_EXPORT &operator<<(QDataStream &ds, const VectorArrayProperty &vecarray)
{
    writeHeader(ds, r_vectorarray, 1);
    
    SharedDataStream sds(ds);
    sds << vecarray.array();
    
    return ds;
}

QDataStream SIREMATHS_EXPORT &operator>>(QDataStream &ds, VectorArrayProperty &vecarray)
{
    VersionID v = readHeader(ds, r_vectorarray);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        QVector<Vector> a;
        sds >> a;
        vecarray = VectorArrayProperty(a);
    }
    else
        throw version_error(v, "1", r_vectorarray, CODELOC);
    
    return ds;
}

VectorArrayProperty::VectorArrayProperty()
                    : ConcreteProperty<VectorArrayProperty, ArrayProperty<Vector> >()
{}

VectorArrayProperty::VectorArrayProperty(const QList<Vector> &array)
                    : ConcreteProperty<VectorArrayProperty, ArrayProperty<Vector> >(array)
{}

VectorArrayProperty::VectorArrayProperty(const QVector<Vector> &array)
                    : ConcreteProperty<VectorArrayProperty, ArrayProperty<Vector> >(array)
{}

VectorArrayProperty::VectorArrayProperty(const VectorArrayProperty &other)
                    : ConcreteProperty<VectorArrayProperty, ArrayProperty<Vector> >(other)
{}

VectorArrayProperty::~VectorArrayProperty()
{}

const char* VectorArrayProperty::typeName()
{
    return QMetaType::typeName(qMetaTypeId<VectorArrayProperty>());
}

VectorArrayProperty& VectorArrayProperty::operator=(const VectorArrayProperty &other)
{
    ArrayProperty<Vector>::operator=(other);
    return *this;
}

bool VectorArrayProperty::operator==(const VectorArrayProperty &other) const
{
    return ArrayProperty<Vector>::operator==(other);
}

bool VectorArrayProperty::operator!=(const VectorArrayProperty &other) const
{
    return ArrayProperty<Vector>::operator!=(other);
}

VectorArrayProperty VectorArrayProperty::operator+(const VectorArrayProperty &other) const
{
    return VectorArrayProperty( ArrayProperty<Vector>::operator+(other) );
}

VectorArrayProperty& VectorArrayProperty::operator+=(const VectorArrayProperty &other)
{
    ArrayProperty<Vector>::operator+=(other);
    return *this;
}

//////////
////////// Implementation of VectorProperty
//////////

static const RegisterMetaType<VectorProperty> r_vecprop;

QDataStream SIREMATHS_EXPORT &operator<<(QDataStream &ds, const VectorProperty &vecprop)
{
    writeHeader(ds, r_vecprop, 1);
    
    ds << static_cast<const Vector&>(vecprop);
    
    return ds;
}

QDataStream SIREMATHS_EXPORT &operator>>(QDataStream &ds, VectorProperty &vecprop)
{
    VersionID v = readHeader(ds, r_vecprop);
    
    if (v == 1)
    {
        ds >> static_cast<Vector&>(vecprop);
    }
    else
        throw version_error(v, "1", r_vecprop, CODELOC);
    
    return ds;
}

/** Constructor */
VectorProperty::VectorProperty() : ConcreteProperty<VectorProperty,Property>(), Vector()
{}

/** Construct a copy of the passed vector */
VectorProperty::VectorProperty(const Vector &value)
               : ConcreteProperty<VectorProperty,Property>(), Vector(value)
{}

/** Copy constructor */
VectorProperty::VectorProperty(const VectorProperty &other)
               : ConcreteProperty<VectorProperty,Property>(other), Vector(other)
{}

/** Destructor */
VectorProperty::~VectorProperty()
{}

const char* VectorProperty::typeName()
{
    return QMetaType::typeName( qMetaTypeId<VectorProperty>() );
}

const char* VectorProperty::what() const
{
    return VectorProperty::typeName();
}

VectorProperty* VectorProperty::clone() const
{
    return new VectorProperty(*this);
}

/** Copy assignment operator */
VectorProperty& VectorProperty::operator=(const VectorProperty &other)
{
    Vector::operator=(other);
    return *this;
}

/** Copy assignment operator */
VectorProperty& VectorProperty::operator=(const Vector &other)
{
    Vector::operator=(other);
    return *this;
}

/** Comparison operator */
bool VectorProperty::operator==(const VectorProperty &other) const
{
    return Vector::operator==(other);
}

/** Comparison operator */
bool VectorProperty::operator!=(const VectorProperty &other) const
{
    return Vector::operator!=(other);
}

/** Return the actual value of the vector */
Vector VectorProperty::value() const
{
    return Vector(*this);
}

QString VectorProperty::toString() const
{
    return Vector::toString();
}

