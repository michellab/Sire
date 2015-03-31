/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#ifndef SIRECAS_ABS_H
#define SIRECAS_ABS_H

#include "singlefunc.h"

SIRE_BEGIN_HEADER

namespace SireCAS
{
class Abs;
}

QDataStream& operator<<(QDataStream&, const SireCAS::Abs&);
QDataStream& operator>>(QDataStream&, SireCAS::Abs&);

namespace SireCAS
{

/**
This is the absolute value, abs.

abs(x) returns x if x >= 0, else -x

For complex values, this returns abs(x) + abs(y) i

@author Christopher Woods
*/
class SIRECAS_EXPORT Abs : public SingleFunc
{

friend QDataStream& ::operator<<(QDataStream&, const Abs&);
friend QDataStream& ::operator>>(QDataStream&, Abs&);

public:
    Abs();
    Abs(const Expression &power);

    Abs(const Abs &other);

    ~Abs();

    bool operator==(const ExBase &other) const;

    static const char* typeName();

    const char* what() const
    {
        return Abs::typeName();
    }

    Abs* clone() const;

    double evaluate(const Values &values) const;
    Complex evaluate(const ComplexValues &values) const;

protected:

    Expression functionOf(const Expression &arg) const
    {
        if (arg == argument())
            return Expression(*this);
        else
            return Expression(Abs(arg));
    }

    QString stringRep() const
    {
        return "abs";
    }

    uint magic() const;

    Expression diff() const;
    Expression integ() const;
};

}

Q_DECLARE_METATYPE(SireCAS::Abs)

SIRE_EXPOSE_CLASS( SireCAS::Abs )

SIRE_END_HEADER

#endif
