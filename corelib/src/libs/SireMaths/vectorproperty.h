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

#ifndef SIREMATHS_VECTORPROPERTY_H
#define SIREMATHS_VECTORPROPERTY_H

#include "SireBase/property.h"
#include "SireBase/arrayproperty.hpp"

#include "SireMaths/vector.h"

SIRE_BEGIN_HEADER

namespace SireMaths
{
class VectorProperty;
class VectorArrayProperty;
}

QDataStream& operator<<(QDataStream&, const SireMaths::VectorProperty&);
QDataStream& operator>>(QDataStream&, SireMaths::VectorProperty&);

QDataStream& operator<<(QDataStream&, const SireMaths::VectorArrayProperty&);
QDataStream& operator>>(QDataStream&, SireMaths::VectorArrayProperty&);

namespace SireMaths
{

/** This class provides a simple Property wrapper around a Vector, thereby
    allowing the vector to be stored as a Property, e.g. for the center
    of a molecule
    
    @author Christopher Woods
*/
class SIREMATHS_EXPORT VectorProperty
            : public SireBase::ConcreteProperty<VectorProperty,SireBase::Property>,
              public Vector
{

friend QDataStream& ::operator<<(QDataStream&, const VectorProperty&);
friend QDataStream& ::operator>>(QDataStream&, VectorProperty&);

public:
    VectorProperty();
    VectorProperty(const Vector &value);
    VectorProperty(const VectorProperty &value);
    
    ~VectorProperty();
    
    static const char* typeName();
    const char* what() const;
    
    VectorProperty& operator=(const VectorProperty &other);
    VectorProperty& operator=(const Vector &other);
    
    bool operator==(const VectorProperty &other) const;
    bool operator!=(const VectorProperty &other) const;
    
    QString toString() const;
    
    VectorProperty* clone() const;
};

class SIREMATHS_EXPORT VectorArrayProperty
        : public SireBase::ConcreteProperty<VectorArrayProperty,SireBase::ArrayProperty<Vector> >
{

friend QDataStream& ::operator<<(QDataStream&, const VectorArrayProperty&);
friend QDataStream& ::operator>>(QDataStream&, VectorArrayProperty&);

public:
    VectorArrayProperty();
    VectorArrayProperty(const QList<Vector> &array);
    VectorArrayProperty(const QVector<Vector> &array);
    VectorArrayProperty(const VectorArrayProperty &other);
    
    ~VectorArrayProperty();
    
    static const char* typeName();
    
    VectorArrayProperty& operator=(const VectorArrayProperty &other);
    
    bool operator==(const VectorArrayProperty &other) const;
    bool operator!=(const VectorArrayProperty &other) const;

    VectorArrayProperty operator+(const VectorArrayProperty &other) const;
    VectorArrayProperty& operator+=(const VectorArrayProperty &other);
};

SireBase::PropertyPtr wrap(const Vector &vector);
SireBase::PropertyPtr wrap(const QVector<Vector> &vector);
SireBase::PropertyPtr wrap(const QList<Vector> &vector);

}

Q_DECLARE_METATYPE( SireMaths::VectorProperty )
Q_DECLARE_METATYPE( SireMaths::VectorArrayProperty )

SIRE_EXPOSE_FUNCTION( SireMaths::wrap )

SIRE_EXPOSE_CLASS( SireMaths::VectorProperty )
SIRE_EXPOSE_CLASS( SireMaths::VectorArrayProperty )

SIRE_END_HEADER

#endif
