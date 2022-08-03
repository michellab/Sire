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

#include "SireMaths/sphereproperty.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireBase;
using namespace SireMaths;
using namespace SireStream;

namespace SireMaths
{
    PropertyPtr wrap(const Sphere &sphere)
    {
        return PropertyPtr( SphereProperty(sphere) );
    }

    PropertyPtr wrap(const QVector<Sphere> &sphere)
    {
        return PropertyPtr( SphereArrayProperty(sphere) );
    }

    PropertyPtr wrap(const QList<Sphere> &sphere)
    {
        return PropertyPtr( SphereArrayProperty(sphere) );
    }
}

//////////
////////// Implementation of SphereArrayProperty
//////////

static const RegisterMetaType<SphereArrayProperty> r_spherearray;

QDataStream &operator<<(QDataStream &ds, const SphereArrayProperty &array)
{
    writeHeader(ds, r_spherearray, 1);

    SharedDataStream sds(ds);
    sds << array.array();

    return ds;
}

QDataStream &operator>>(QDataStream &ds, SphereArrayProperty &array)
{
    VersionID v = readHeader(ds, r_spherearray);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        QVector<Sphere> a;
        sds >> a;
        array = SphereArrayProperty(a);
    }
    else
        throw version_error(v, "1", r_spherearray, CODELOC);

    return ds;
}

SphereArrayProperty::SphereArrayProperty()
                    : ConcreteProperty<SphereArrayProperty, ArrayProperty<Sphere> >()
{}

SphereArrayProperty::SphereArrayProperty(const QList<Sphere> &array)
                    : ConcreteProperty<SphereArrayProperty, ArrayProperty<Sphere> >(array)
{}

SphereArrayProperty::SphereArrayProperty(const QVector<Sphere> &array)
                    : ConcreteProperty<SphereArrayProperty, ArrayProperty<Sphere> >(array)
{}

SphereArrayProperty::SphereArrayProperty(const SphereArrayProperty &other)
                    : ConcreteProperty<SphereArrayProperty, ArrayProperty<Sphere> >(other)
{}

SphereArrayProperty::~SphereArrayProperty()
{}

const char* SphereArrayProperty::typeName()
{
    return QMetaType::typeName(qMetaTypeId<SphereArrayProperty>());
}

SphereArrayProperty& SphereArrayProperty::operator=(const SphereArrayProperty &other)
{
    ArrayProperty<Sphere>::operator=(other);
    return *this;
}

bool SphereArrayProperty::operator==(const SphereArrayProperty &other) const
{
    return ArrayProperty<Sphere>::operator==(other);
}

bool SphereArrayProperty::operator!=(const SphereArrayProperty &other) const
{
    return ArrayProperty<Sphere>::operator!=(other);
}

SphereArrayProperty SphereArrayProperty::operator+(const SphereArrayProperty &other) const
{
    return SphereArrayProperty( ArrayProperty<Sphere>::operator+(other) );
}

SphereArrayProperty& SphereArrayProperty::operator+=(const SphereArrayProperty &other)
{
    ArrayProperty<Sphere>::operator+=(other);
    return *this;
}

//////////
////////// Implementation of SphereProperty
//////////

static const RegisterMetaType<SphereProperty> r_sphereprop;

QDataStream &operator<<(QDataStream &ds, const SphereProperty &vecprop)
{
    writeHeader(ds, r_sphereprop, 1);

    ds << static_cast<const Sphere&>(vecprop);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, SphereProperty &vecprop)
{
    VersionID v = readHeader(ds, r_sphereprop);

    if (v == 1)
    {
        ds >> static_cast<Sphere&>(vecprop);
    }
    else
        throw version_error(v, "1", r_sphereprop, CODELOC);

    return ds;
}

/** Constructor */
SphereProperty::SphereProperty() : ConcreteProperty<SphereProperty,Property>(), Sphere()
{}

/** Construct a copy of the passed Sphere */
SphereProperty::SphereProperty(const Sphere &value)
               : ConcreteProperty<SphereProperty,Property>(), Sphere(value)
{}

/** Copy constructor */
SphereProperty::SphereProperty(const SphereProperty &other)
               : ConcreteProperty<SphereProperty,Property>(other), Sphere(other)
{}

/** Destructor */
SphereProperty::~SphereProperty()
{}

const char* SphereProperty::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SphereProperty>() );
}

const char* SphereProperty::what() const
{
    return SphereProperty::typeName();
}

SphereProperty* SphereProperty::clone() const
{
    return new SphereProperty(*this);
}

/** Copy assignment operator */
SphereProperty& SphereProperty::operator=(const SphereProperty &other)
{
    Sphere::operator=(other);
    return *this;
}

/** Copy assignment operator */
SphereProperty& SphereProperty::operator=(const Sphere &other)
{
    Sphere::operator=(other);
    return *this;
}

/** Comparison operator */
bool SphereProperty::operator==(const SphereProperty &other) const
{
    return Sphere::operator==(other);
}

/** Comparison operator */
bool SphereProperty::operator!=(const SphereProperty &other) const
{
    return Sphere::operator!=(other);
}

/** Return the actual value of the Sphere */
Sphere SphereProperty::value() const
{
    return Sphere(*this);
}

QString SphereProperty::toString() const
{
    return Sphere::toString();
}

