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

#ifndef SIRECAS_EXP_H
#define SIRECAS_EXP_H

#include "power.h"
#include "singlefunc.h"

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireCAS
{
class Exp;
class Ln;
}

SIRECAS_EXPORT QDataStream& operator<<(QDataStream&, const SireCAS::Exp&);
SIRECAS_EXPORT QDataStream& operator>>(QDataStream&, SireCAS::Exp&);

SIRECAS_EXPORT QDataStream& operator<<(QDataStream&, const SireCAS::Ln&);
SIRECAS_EXPORT QDataStream& operator>>(QDataStream&, SireCAS::Ln&);

namespace SireCAS
{

/**
This is the exponential function, e^x

@author Christopher Woods
*/
class SIRECAS_EXPORT Exp : public PowerFunction
{

friend SIRECAS_EXPORT QDataStream& ::operator<<(QDataStream&, const Exp&);
friend SIRECAS_EXPORT QDataStream& ::operator>>(QDataStream&, Exp&);

public:
    Exp();
    Exp(const Expression &power);

    Exp(const Exp &other);

    ~Exp();

    bool operator==(const ExBase &other) const;

    uint hash() const;

    static const char* typeName();

    const char* what() const
    {
        return Exp::typeName();
    }

    Exp* clone() const;

    QString toString() const;

    double evaluate(const Values &values) const;
    Complex evaluate(const ComplexValues &values) const;

    Expression differentiate(const Symbol &symbol) const;
    Expression integrate(const Symbol &symbol) const;

    Expression core() const;
    Expression power() const;

private:

    /** The expression to which 'e' is raised to */
    Expression pwr;
};

/** This is the natural logarithm (ln) function

@author Christopher Woods
*/
class SIRECAS_EXPORT Ln : public SingleFunc
{

friend SIRECAS_EXPORT QDataStream& ::operator<<(QDataStream&, const Ln&);
friend SIRECAS_EXPORT QDataStream& ::operator>>(QDataStream&, Ln&);

public:
    Ln();
    Ln(const Expression &expression);

    Ln(const Ln &other);
    ~Ln();

    bool operator==(const ExBase &other) const;

    static const char* typeName();

    const char* what() const
    {
        return Ln::typeName();
    }

    Ln* clone() const;

    double evaluate(const Values &values) const;
    Complex evaluate(const ComplexValues &values) const;

protected:

    Expression functionOf(const Expression &arg) const
    {
        if (arg == argument())
            return Expression(*this);
        else
            return Expression(Ln(arg));
    }

    QString stringRep() const
    {
        return "ln";
    }

    uint magic() const;

    Expression diff() const;
    Expression integ() const;
};

}

Q_DECLARE_METATYPE(SireCAS::Exp)
Q_DECLARE_METATYPE(SireCAS::Ln)

SIRE_EXPOSE_CLASS( SireCAS::Exp )
SIRE_EXPOSE_CLASS( SireCAS::Ln )

SIRE_END_HEADER

#endif
