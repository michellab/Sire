/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#ifndef SIRECAS_MINMAX_H
#define SIRECAS_MINMAX_H

#include "doublefunc.h"

SIRE_BEGIN_HEADER

namespace SireCAS
{
class Min;
class Max;
}

SIRECAS_EXPORT QDataStream& operator<<(QDataStream&, const SireCAS::Min&);
SIRECAS_EXPORT QDataStream& operator>>(QDataStream&, SireCAS::Min&);

SIRECAS_EXPORT QDataStream& operator<<(QDataStream&, const SireCAS::Max&);
SIRECAS_EXPORT QDataStream& operator>>(QDataStream&, SireCAS::Max&);

namespace SireCAS
{

/** Minimum value of two expressions */
class SIRECAS_EXPORT Min : public DoubleFunc
{

friend SIRECAS_EXPORT QDataStream& ::operator<<(QDataStream&, const Min&);
friend SIRECAS_EXPORT QDataStream& ::operator>>(QDataStream&, Min&);

public:
    Min();
    Min(const Expression &x, const Expression &y);

    Min(const Min &other);

    ~Min();

    bool operator==(const ExBase &other) const;

    static const char* typeName();

    const char* what() const
    {
        return Min::typeName();
    }

    Min* clone() const;

    double evaluate(const Values &values) const;
    Complex evaluate(const ComplexValues &values) const;

protected:
    Expression functionOf(const Expression &x, const Expression &y) const;

    QString stringRep() const
    {
        return "min";
    }

    uint magic() const;
};

/** Maximum value of two expressions */
class SIRECAS_EXPORT Max : public DoubleFunc
{

friend SIRECAS_EXPORT QDataStream& ::operator<<(QDataStream&, const Max&);
friend SIRECAS_EXPORT QDataStream& ::operator>>(QDataStream&, Max&);

public:
    Max();
    Max(const Expression &x, const Expression &y);

    Max(const Max &other);

    ~Max();

    bool operator==(const ExBase &other) const;

    static const char* typeName();

    const char* what() const
    {
        return Max::typeName();
    }

    Max* clone() const;

    double evaluate(const Values &values) const;
    Complex evaluate(const ComplexValues &values) const;

protected:
    Expression functionOf(const Expression &x, const Expression &y) const;

    QString stringRep() const
    {
        return "max";
    }

    uint magic() const;
};

}

Q_DECLARE_METATYPE( SireCAS::Min )
Q_DECLARE_METATYPE( SireCAS::Max )

SIRE_EXPOSE_CLASS( SireCAS::Min )
SIRE_EXPOSE_CLASS( SireCAS::Max )

SIRE_END_HEADER

#endif
