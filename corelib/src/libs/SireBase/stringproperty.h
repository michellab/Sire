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

#ifndef SIREBASE_STRINGPROPERTY_H
#define SIREBASE_STRINGPROPERTY_H

#include "property.h"

SIRE_BEGIN_HEADER

namespace SireBase
{
class StringProperty;
class VariantProperty;
}

QDataStream& operator<<(QDataStream&, const SireBase::StringProperty&);
QDataStream& operator>>(QDataStream&, SireBase::StringProperty&);

namespace SireBase
{

/** This class provides a thin Property wrapper around a QString 

    @author Christopher Woods
*/
class SIREBASE_EXPORT StringProperty : public ConcreteProperty<StringProperty,Property>
{

friend QDataStream& ::operator<<(QDataStream&, const StringProperty&);
friend QDataStream& ::operator>>(QDataStream&, StringProperty&);

public:
    StringProperty();
    StringProperty(const QString &s);
    StringProperty(const VariantProperty &other);
    StringProperty(const StringProperty &other);
    
    ~StringProperty();
    
    StringProperty& operator=(const StringProperty &other);
    
    bool operator==(const StringProperty &other) const;
    bool operator!=(const StringProperty &other) const;
    
    static const char* typeName();
    
    QString toString() const;
    
    operator QString() const;

private:
    /** The actual string */
    QString s;
};

}

Q_DECLARE_METATYPE( SireBase::StringProperty )

SIRE_EXPOSE_CLASS( SireBase::StringProperty )

SIRE_END_HEADER

#endif


