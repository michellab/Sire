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

#ifndef SIREBASE_BOOLEANPROPERTY_H
#define SIREBASE_BOOLEANPROPERTY_H

#include "SireBase/property.h"

SIRE_BEGIN_HEADER

namespace SireBase
{
class BooleanProperty;
}

QDataStream& operator<<(QDataStream&, const SireBase::BooleanProperty&);
QDataStream& operator>>(QDataStream&, SireBase::BooleanProperty&);

namespace SireBase
{

/** This class provides a thin Property wrapper around bools

    @author Christopher Woods
*/
class SIREBASE_EXPORT BooleanProperty : public ConcreteProperty<BooleanProperty,Property>
{

friend QDataStream& ::operator<<(QDataStream&, const BooleanProperty&);
friend QDataStream& ::operator>>(QDataStream&, BooleanProperty&);

public:
    BooleanProperty();
    BooleanProperty(const QString &value);
    BooleanProperty(bool value);

    BooleanProperty(const Property &other);
    BooleanProperty(const BooleanProperty &other);
    
    ~BooleanProperty();
    
    static const char* typeName();
    
    BooleanProperty& operator=(const BooleanProperty &other);
    
    bool operator==(const BooleanProperty &other) const;
    bool operator!=(const BooleanProperty &other) const;
    
    bool value() const;
    
    QString toString() const;
    
    bool isAString() const;
    bool isADouble() const;
    bool isAnInteger() const;
    bool isABoolean() const;
    
    QString asAString() const;
    double asADouble() const;
    int asAnInteger() const;
    bool asABoolean() const;
    
private:
    bool val;
};

}

Q_DECLARE_METATYPE( SireBase::BooleanProperty )

SIRE_EXPOSE_CLASS( SireBase::BooleanProperty )

SIRE_END_HEADER

#endif
