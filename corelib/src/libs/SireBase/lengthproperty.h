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

#ifndef SIREBASE_LENGTHPROPERTY_H
#define SIREBASE_LENGTHPROPERTY_H

#include "SireBase/property.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireBase
{
class LengthProperty;
}

QDataStream& operator<<(QDataStream&, const SireBase::LengthProperty&);
QDataStream& operator>>(QDataStream&, SireBase::LengthProperty&);

namespace SireBase
{

using SireUnits::Dimension::Length;

/** This class provides a thin Property wrapper around lengths

    @author Christopher Woods
*/
class SIREBASE_EXPORT LengthProperty : public ConcreteProperty<LengthProperty,Property>
{

friend QDataStream& ::operator<<(QDataStream&, const LengthProperty&);
friend QDataStream& ::operator>>(QDataStream&, LengthProperty&);

public:
    LengthProperty();
    LengthProperty(Length value);

    LengthProperty(const LengthProperty &other);
    LengthProperty(const Property &other);
    
    ~LengthProperty();
    
    static const char* typeName();
    
    LengthProperty& operator=(const LengthProperty &other);
    
    bool operator==(const LengthProperty &other) const;
    bool operator!=(const LengthProperty &other) const;
    
    Length value() const;
    
    QString toString() const;
    
private:
    Length val;
};

}

Q_DECLARE_METATYPE( SireBase::LengthProperty )

SIRE_EXPOSE_CLASS( SireBase::LengthProperty )

SIRE_END_HEADER

#endif
