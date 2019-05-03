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

#ifndef SIRECAS_POWER_H
#define SIRECAS_POWER_H

#include "exbase.h"
#include "expression.h"
#include "symbols.h"
#include "functions.h"

SIRE_BEGIN_HEADER

namespace SireCAS
{
class PowerFunction;
class Power;
}

SIRECAS_EXPORT QDataStream& operator<<(QDataStream&, const SireCAS::PowerFunction&);
SIRECAS_EXPORT QDataStream& operator>>(QDataStream&, SireCAS::PowerFunction&);

SIRECAS_EXPORT QDataStream& operator<<(QDataStream&, const SireCAS::Power&);
SIRECAS_EXPORT QDataStream& operator>>(QDataStream&, SireCAS::Power&);

namespace SireCAS
{

/**
This is the base class of all power expressions, e.g. x^y (all of the form core^power). There are several sub-classes that depend on exactly what is being raised to which power, e.g. Exp is e^y, Power is x^y, PowerConstant is c^y and ConstantPower is x^c (with ConstantPower further derived into RationalPower and RealPower based on whether the constant is rational). All of these can be constructed transparently by creating a Power and then calling 'reduce' on the resulting object.

@author Christopher Woods
*/
class SIRECAS_EXPORT PowerFunction : public ExBase
{

friend SIRECAS_EXPORT QDataStream& ::operator<<(QDataStream&, const PowerFunction&);
friend SIRECAS_EXPORT QDataStream& ::operator>>(QDataStream&, PowerFunction&);

public:
    PowerFunction() : ExBase()
    {}

    ~PowerFunction()
    {}

    virtual Expression core() const=0;
    virtual Expression power() const=0;

    static const char* typeName()
    {
        return "SireCAS::PowerFunction";
    }

    Symbols symbols() const
    {
        Symbols s = core().symbols();
        s.insert(power().symbols());
        return s;
    }

    Functions functions() const
    {
        Functions f = core().functions();
        f.insert(power().functions());
        return f;
    }

    bool isCompound() const
    {
        return false;
    }

    QString toString() const;
    QString toOpenMMString() const;

    Expression substitute(const Identities &identities) const;

    Expression differentiate(const Symbol &symbol) const;
    Expression integrate(const Symbol &symbol) const;

    Expressions children() const;

    bool isFunction(const Symbol &symbol) const;
    bool isConstant() const;

    Expression reduce() const;
    
    QList<Factor> expand(const Symbol &symbol) const;
};

/**
This class represents an expression raised to a generic power (e.g. x^y). This is also the route to raising expressions to real-number powers, and the base of the implementation of the exp() and invlog_10() functions.

@author Christopher Woods
*/
class SIRECAS_EXPORT Power : public PowerFunction
{

friend SIRECAS_EXPORT QDataStream& ::operator<<(QDataStream&, const Power&);
friend SIRECAS_EXPORT QDataStream& ::operator>>(QDataStream&, Power&);

public:
    Power();
    Power(const Expression &base, const Expression &power);

    Power(const Power &other);

    ~Power();

    bool operator==(const ExBase &other) const;

    uint hash() const;

    static const char* typeName();

    const char* what() const
    {
        return Power::typeName();
    }

    double evaluate(const Values &values) const;
    Complex evaluate(const ComplexValues &values) const;

    Expression core() const;
    Expression power() const;

    bool isCompound() const
    {
        return pwr.isCompound();
    }

    Power* clone() const;

private:

    /** The core expression */
    Expression ex;

    /** The power */
    Expression pwr;

};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the core of the power (e.g. return x for x^y) */
inline Expression Power::core() const
{
    return ex;
}

/** Return the power of the power (e.g. return y for x^y) */
inline Expression Power::power() const
{
    return pwr;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE(SireCAS::Power)

SIRE_EXPOSE_CLASS( SireCAS::PowerFunction )
SIRE_EXPOSE_CLASS( SireCAS::Power )

SIRE_END_HEADER

#endif
