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

#ifndef SIREBASE_NUMBERPROPERTY_H
#define SIREBASE_NUMBERPROPERTY_H

#include "SireBase/property.h"

SIRE_BEGIN_HEADER

namespace SireBase
{
class NumberProperty;
class VariantProperty;
}

QDataStream& operator<<(QDataStream&, const SireBase::NumberProperty&);
QDataStream& operator>>(QDataStream&, SireBase::NumberProperty&);

namespace SireBase
{

/** This class provides a thin Property wrapper around numbers (doubles and ints)

    @author Christopher Woods
*/
class SIREBASE_EXPORT NumberProperty : public ConcreteProperty<NumberProperty,Property>
{

friend QDataStream& ::operator<<(QDataStream&, const NumberProperty&);
friend QDataStream& ::operator>>(QDataStream&, NumberProperty&);

public:
    NumberProperty();
    NumberProperty(double value);
    NumberProperty(int value);
    NumberProperty(qint64 value);

    NumberProperty(const VariantProperty &other);    
    NumberProperty(const NumberProperty &other);
    
    ~NumberProperty();
    
    static const char* typeName();
    
    NumberProperty& operator=(const NumberProperty &other);
    
    bool operator==(const NumberProperty &other) const;
    bool operator!=(const NumberProperty &other) const;
    
    qint64 toInt() const;
    double toDouble() const;
    
    double value() const;
    
    QString toString() const;
    
private:
    union
    {
        double dval;
        qint64 ival;
    } val;
    
    bool is_int;
};

}

Q_DECLARE_METATYPE( SireBase::NumberProperty )

SIRE_EXPOSE_CLASS( SireBase::NumberProperty )

SIRE_END_HEADER

#endif
