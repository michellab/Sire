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

#ifndef SIRECAS_EXPRESSIONPROPERTY_H
#define SIRECAS_EXPRESSIONPROPERTY_H

#include "SireBase/property.h"
#include "SireCAS/expression.h"

SIRE_BEGIN_HEADER

namespace SireCAS
{
class ExpressionProperty;
}

SIRECAS_EXPORT QDataStream& operator<<(QDataStream&, const SireCAS::ExpressionProperty&);
SIRECAS_EXPORT QDataStream& operator>>(QDataStream&, SireCAS::ExpressionProperty&);

namespace SireCAS
{

/** This class provides a thin Property wrapper around SireCAS objects

    @author Christopher Woods
*/
class SIRECAS_EXPORT ExpressionProperty
    : public SireBase::ConcreteProperty<ExpressionProperty,SireBase::Property>
{

friend SIRECAS_EXPORT QDataStream& ::operator<<(QDataStream&, const ExpressionProperty&);
friend SIRECAS_EXPORT QDataStream& ::operator>>(QDataStream&, ExpressionProperty&);

public:
    ExpressionProperty();
    ExpressionProperty(const ExBase &exbase);
    ExpressionProperty(const Expression &expression);
    ExpressionProperty(const ExpressionProperty &other);
    
    ~ExpressionProperty();
    
    static const char* typeName();
    
    ExpressionProperty& operator=(const ExpressionProperty &other);
    
    bool operator==(const ExpressionProperty &other) const;
    bool operator!=(const ExpressionProperty &other) const;
    
    QString toString() const;
    
    Expression value() const;
    
    bool isADouble() const;
    bool isAnInteger() const;
    bool isABoolean() const;
    
    double asADouble() const;
    int asAnInteger() const;
    bool asABoolean() const;
    
private:
    Expression _val;
};

SIRECAS_EXPORT SireBase::PropertyPtr wrap(const ExBase &expression);
SIRECAS_EXPORT SireBase::PropertyPtr wrap(const Expression &expression);

}

Q_DECLARE_METATYPE( SireCAS::ExpressionProperty )

SIRE_EXPOSE_CLASS( SireCAS::ExpressionProperty )

SIRE_EXPOSE_FUNCTION( SireCAS::wrap )

SIRE_END_HEADER

#endif
