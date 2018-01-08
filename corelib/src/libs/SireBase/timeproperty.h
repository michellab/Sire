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

#ifndef SIREBASE_TIMEPROPERTY_H
#define SIREBASE_TIMEPROPERTY_H

#include "SireBase/property.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireBase
{
class TimeProperty;
}

QDataStream& operator<<(QDataStream&, const SireBase::TimeProperty&);
QDataStream& operator>>(QDataStream&, SireBase::TimeProperty&);

namespace SireBase
{

using SireUnits::Dimension::Time;

/** This class provides a thin Property wrapper around times

    @author Christopher Woods
*/
class SIREBASE_EXPORT TimeProperty : public ConcreteProperty<TimeProperty,Property>
{

friend QDataStream& ::operator<<(QDataStream&, const TimeProperty&);
friend QDataStream& ::operator>>(QDataStream&, TimeProperty&);

public:
    TimeProperty();
    TimeProperty(Time value);

    TimeProperty(const TimeProperty &other);
    TimeProperty(const Property &other);    
    
    ~TimeProperty();
    
    static const char* typeName();
    
    TimeProperty& operator=(const TimeProperty &other);
    
    bool operator==(const TimeProperty &other) const;
    bool operator!=(const TimeProperty &other) const;
    
    Time value() const;
    
    QString toString() const;
    
    operator Time() const;
    
private:
    Time val;
};

}

Q_DECLARE_METATYPE( SireBase::TimeProperty )

SIRE_EXPOSE_CLASS( SireBase::TimeProperty )

SIRE_END_HEADER

#endif
