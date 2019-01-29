/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIRECAS_INTEGRATIONCONSTANT_H
#define SIRECAS_INTEGRATIONCONSTANT_H

#include "symbol.h"

SIRE_BEGIN_HEADER

namespace SireCAS
{
class IntegrationConstant;
}

SIRECAS_EXPORT QDataStream& operator<<(QDataStream&, const SireCAS::IntegrationConstant&);
SIRECAS_EXPORT QDataStream& operator>>(QDataStream&, SireCAS::IntegrationConstant&);

namespace SireCAS
{

/**
This class represents a constant of integration. This is not a number or 
function in the normal sense, but rather a placeholder that is created 
during indefinite integration.

@author Christopher Woods
*/
class SIRECAS_EXPORT IntegrationConstant : public Symbol
{

friend QDataStream& ::operator<<(QDataStream&, const IntegrationConstant&);
friend QDataStream& ::operator>>(QDataStream&, IntegrationConstant&);

public:
    IntegrationConstant();

    IntegrationConstant(const IntegrationConstant &other);

    ~IntegrationConstant();

    bool operator==(const ExBase &other) const;

    uint hash() const;

    static const char* typeName();

    const char* what() const
    {
        return IntegrationConstant::typeName();
    }

    IntegrationConstant* clone() const;

    Expression integrate(const Symbol &symbol) const;
};

}

Q_DECLARE_METATYPE(SireCAS::IntegrationConstant)

SIRE_EXPOSE_CLASS( SireCAS::IntegrationConstant )

SIRE_END_HEADER

#endif
