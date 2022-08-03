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

#ifndef SIREMATHS_SPHEREPROPERTY_H
#define SIREMATHS_SPHEREPROPERTY_H

#include "SireBase/property.h"
#include "SireBase/arrayproperty.hpp"

#include "SireMaths/sphere.h"

SIRE_BEGIN_HEADER

namespace SireMaths
{
class SphereProperty;
class SphereArrayProperty;
}

SIREMATHS_EXPORT QDataStream& operator<<(QDataStream&, const SireMaths::SphereProperty&);
SIREMATHS_EXPORT QDataStream& operator>>(QDataStream&, SireMaths::SphereProperty&);

SIREMATHS_EXPORT QDataStream& operator<<(QDataStream&, const SireMaths::SphereArrayProperty&);
SIREMATHS_EXPORT QDataStream& operator>>(QDataStream&, SireMaths::SphereArrayProperty&);

namespace SireMaths
{

/** This class provides a simple Property wrapper around a Vector, thereby
    allowing the vector to be stored as a Property, e.g. for the center
    of a molecule

    @author Christopher Woods
*/
class SIREMATHS_EXPORT SphereProperty
            : public SireBase::ConcreteProperty<SphereProperty,SireBase::Property>,
              public Sphere
{

friend SIREMATHS_EXPORT QDataStream& ::operator<<(QDataStream&, const SphereProperty&);
friend SIREMATHS_EXPORT QDataStream& ::operator>>(QDataStream&, SphereProperty&);

public:
    SphereProperty();
    SphereProperty(const Sphere &value);
    SphereProperty(const SphereProperty &value);

    ~SphereProperty();

    static const char* typeName();
    const char* what() const;

    SphereProperty& operator=(const SphereProperty &other);
    SphereProperty& operator=(const Sphere &other);

    bool operator==(const SphereProperty &other) const;
    bool operator!=(const SphereProperty &other) const;

    QString toString() const;

    Sphere value() const;

    SphereProperty* clone() const;
};

class SIREMATHS_EXPORT SphereArrayProperty
        : public SireBase::ConcreteProperty<SphereArrayProperty,SireBase::ArrayProperty<Sphere> >
{

friend SIREMATHS_EXPORT QDataStream& ::operator<<(QDataStream&, const SphereArrayProperty&);
friend SIREMATHS_EXPORT QDataStream& ::operator>>(QDataStream&, SphereArrayProperty&);

public:
    SphereArrayProperty();
    SphereArrayProperty(const QList<Sphere> &array);
    SphereArrayProperty(const QVector<Sphere> &array);
    SphereArrayProperty(const SphereArrayProperty &other);

    ~SphereArrayProperty();

    static const char* typeName();

    SphereArrayProperty& operator=(const SphereArrayProperty &other);

    bool operator==(const SphereArrayProperty &other) const;
    bool operator!=(const SphereArrayProperty &other) const;

    SphereArrayProperty operator+(const SphereArrayProperty &other) const;
    SphereArrayProperty& operator+=(const SphereArrayProperty &other);
};

SIREMATHS_EXPORT SireBase::PropertyPtr wrap(const Sphere &sphere);
SIREMATHS_EXPORT SireBase::PropertyPtr wrap(const QVector<Sphere> &sphere);
SIREMATHS_EXPORT SireBase::PropertyPtr wrap(const QList<Sphere> &sphere);

}

Q_DECLARE_METATYPE( SireMaths::SphereProperty )
Q_DECLARE_METATYPE( SireMaths::SphereArrayProperty )

SIRE_EXPOSE_FUNCTION( SireMaths::wrap )

SIRE_EXPOSE_CLASS( SireMaths::SphereProperty )
SIRE_EXPOSE_CLASS( SireMaths::SphereArrayProperty )

SIRE_EXPOSE_ALIAS( (SireBase::ArrayProperty<SireMaths::Sphere>),
                    SireBase::ArrayProperty_Sphere_ )

#ifdef SIRE_INSTANTIATE_TEMPLATES
template class SireBase::ArrayProperty<SireMaths::Sphere>;
#endif

SIRE_END_HEADER

#endif
