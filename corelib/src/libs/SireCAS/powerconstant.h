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

#ifndef SIRECAS_POWERCONSTANT_H
#define SIRECAS_POWERCONSTANT_H

#include "power.h"

#include "SireMaths/complex.h"

SIRE_BEGIN_HEADER

namespace SireCAS
{
class PowerConstant;
class ConstantPower;
class IntegerPower;
class RationalPower;
class RealPower;
class ComplexPower;
}

QDataStream& operator<<(QDataStream&, const SireCAS::PowerConstant&);
QDataStream& operator>>(QDataStream&, SireCAS::PowerConstant&);

QDataStream& operator<<(QDataStream&, const SireCAS::ConstantPower&);
QDataStream& operator>>(QDataStream&, SireCAS::ConstantPower&);

QDataStream& operator<<(QDataStream&, const SireCAS::IntegerPower&);
QDataStream& operator>>(QDataStream&, SireCAS::IntegerPower&);

QDataStream& operator<<(QDataStream&, const SireCAS::RationalPower&);
QDataStream& operator>>(QDataStream&, SireCAS::RationalPower&);

QDataStream& operator<<(QDataStream&, const SireCAS::RealPower&);
QDataStream& operator>>(QDataStream&, SireCAS::RealPower&);

QDataStream& operator<<(QDataStream&, const SireCAS::ComplexPower&);
QDataStream& operator>>(QDataStream&, SireCAS::ComplexPower&);

namespace SireCAS
{

using SireMaths::Complex;

/**
This class represents a constant raised to a generic power, e.g. 10^x

@author Christopher Woods
*/
class SIRECAS_EXPORT PowerConstant : public PowerFunction
{

friend QDataStream& ::operator<<(QDataStream&, const PowerConstant&);
friend QDataStream& ::operator>>(QDataStream&, PowerConstant&);

public:
    PowerConstant();
    PowerConstant(double val, const Expression &power);

    PowerConstant(const PowerConstant &other);

    ~PowerConstant();

    bool operator==(const ExBase &other) const;

    uint hash() const;

    static const char* typeName();

    const char* what() const
    {
        return PowerConstant::typeName();
    }

    PowerConstant* clone() const;

    double evaluate(const Values &values) const;
    Complex evaluate(const ComplexValues &values) const;

    Expression core() const
    {
        return Expression(cre);
    }

    Expression power() const
    {
        return pwr;
    }

private:

    /** The constant value */
    double cre;

    /** The expression by which the constant is raised to the power */
    Expression pwr;

};

/**
This class represents an expression raised to a constant power. This is the
base class of RationalPower (expression raised to a rational power) and
RealPower (expression raised to a non-rational power).

@author Christopher Woods
*/
class SIRECAS_EXPORT ConstantPower : public PowerFunction
{

friend QDataStream& ::operator<<(QDataStream&, const ConstantPower&);
friend QDataStream& ::operator>>(QDataStream&, ConstantPower&);

public:
    ConstantPower() : PowerFunction()
    {}

    ConstantPower(const Expression &expression) : PowerFunction(), ex(expression)
    {}

    ConstantPower(const ConstantPower &other) : PowerFunction(), ex(other.ex)
    {}

    ~ConstantPower()
    {}

    Expression core() const
    {
        return ex;
    }

    uint hash() const;

protected:

    /** The expression that is raised to a power */
    Expression ex;
};

/** This class represents an expression raised to a constant integer power */
class SIRECAS_EXPORT IntegerPower : public ConstantPower
{

friend QDataStream& ::operator<<(QDataStream&, const IntegerPower&);
friend QDataStream& ::operator>>(QDataStream&, IntegerPower&);

public:
    IntegerPower();
    IntegerPower(const Expression &expression, int power);

    IntegerPower(const IntegerPower &other);

    ~IntegerPower();

    bool operator==(const ExBase &other) const;

    uint hash() const;

    static const char* typeName();

    const char* what() const
    {
        return IntegerPower::typeName();
    }

    double evaluate(const Values &values) const;
    Complex evaluate(const ComplexValues &values) const;

    Expression power() const
    {
        return Expression(pwr);
    }

    IntegerPower* clone() const;

private:

    /** The integer power */
    int pwr;
};


/** This class represents an expression raised to a rational power */
class SIRECAS_EXPORT RationalPower : public ConstantPower
{

friend QDataStream& ::operator<<(QDataStream&, const RationalPower&);
friend QDataStream& ::operator>>(QDataStream&, RationalPower&);

public:
    RationalPower();
    RationalPower(const Expression &expression, const Rational &power);

    RationalPower(const RationalPower &other);

    ~RationalPower();

    bool operator==(const ExBase &other) const;

    uint hash() const;

    static const char* typeName();

    const char* what() const
    {
        return RationalPower::typeName();
    }

    double evaluate(const Values &values) const;
    Complex evaluate(const ComplexValues &values) const;

    Expression power() const
    {
        return Expression(pwr);
    }

    RationalPower* clone() const;

private:

    /** The rational power */
    Rational pwr;
};

/** This class represents an expression raised to a real power */
class SIRECAS_EXPORT RealPower : public ConstantPower
{

friend QDataStream& ::operator<<(QDataStream&, const RealPower&);
friend QDataStream& ::operator>>(QDataStream&, RealPower&);

public:
    RealPower();
    RealPower(const Expression &expression, double power);

    RealPower(const RealPower &other);

    ~RealPower();

    bool operator==(const ExBase &other) const;

    uint hash() const;

    static const char* typeName();

    const char* what() const
    {
        return RealPower::typeName();
    }

    double evaluate(const Values &values) const;
    Complex evaluate(const ComplexValues &values) const;

    Expression power() const
    {
        return Expression(pwr);
    }

    RealPower* clone() const;

private:

    /** The real power */
    double pwr;
};

/** This class represents an expression raised to a complex power */
class SIRECAS_EXPORT ComplexPower : public ConstantPower
{

friend QDataStream& ::operator<<(QDataStream&, const ComplexPower&);
friend QDataStream& ::operator>>(QDataStream&, ComplexPower&);

public:
    ComplexPower();
    ComplexPower(const Expression &expression, const Complex &power);

    ComplexPower(const ComplexPower &other);

    ~ComplexPower();

    bool operator==(const ExBase &other) const;

    uint hash() const;

    static const char* typeName();

    const char* what() const
    {
        return ComplexPower::typeName();
    }

    double evaluate(const Values &values) const;
    Complex evaluate(const ComplexValues &values) const;

    Expression power() const
    {
        return Expression(pwr);
    }

    bool isComplex() const
    {
        return true;
    }

    ComplexPower* clone() const;

private:

    /** The complex power */
    Complex pwr;
};

}

Q_DECLARE_METATYPE(SireCAS::PowerConstant)
Q_DECLARE_METATYPE(SireCAS::IntegerPower)
Q_DECLARE_METATYPE(SireCAS::RationalPower)
Q_DECLARE_METATYPE(SireCAS::RealPower)
Q_DECLARE_METATYPE(SireCAS::ComplexPower)

SIRE_EXPOSE_CLASS( SireCAS::PowerConstant )
SIRE_EXPOSE_CLASS( SireCAS::IntegerPower )
SIRE_EXPOSE_CLASS( SireCAS::RationalPower )
SIRE_EXPOSE_CLASS( SireCAS::RealPower )
SIRE_EXPOSE_CLASS( SireCAS::ComplexPower )

SIRE_END_HEADER

#endif
