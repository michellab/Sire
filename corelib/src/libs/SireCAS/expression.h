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

#ifndef SIRECAS_EXPRESSION_H
#define SIRECAS_EXPRESSION_H

#include "SireMaths/rational.h"

#include "expressionbase.h"
#include "expressions.h"

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireCAS
{
class Expression;
}

QDataStream& operator<<(QDataStream&, const SireCAS::Expression&);
QDataStream& operator>>(QDataStream&, SireCAS::Expression&);

namespace SireCAS
{

using SireMaths::Rational;

class Identities;
class Symbols;
class Functions;
class Factor;

/**
An Expression is the base class of all algebraic entities.

@author Christopher Woods
*/
class SIRECAS_EXPORT Expression
{

friend QDataStream& ::operator<<(QDataStream&, const Expression&);
friend QDataStream& ::operator>>(QDataStream&, Expression&);

public:
    Expression();

    Expression(int constant);
    Expression(const Rational &constant);
    Expression(double constant);
    Expression(const Complex &constant);

    Expression(const ExpressionBase &base);
    Expression(const ExBase &base);

    Expression(const Expression &other);

    ~Expression();

    static const char* typeName();
    
    const char* what() const
    {
        return Expression::typeName();
    }

    bool operator==(const Expression &other) const;
    bool operator!=(const Expression &other) const;

    Expression& operator=(const Expression &other);

    Expression& operator+=(const Expression &other);
    Expression& operator-=(const Expression &other);
    Expression& operator*=(const Expression &other);
    Expression& operator/=(const Expression &other);

    double operator()(const Values &values) const;
    Complex operator()(const ComplexValues &values) const;

    Expression operator-() const;

    Expression add(const Expression &ex) const;
    Expression add(double val) const;
    Expression add(const Complex &val) const;

    Expression subtract(const Expression &ex) const;
    Expression subtract(double val) const;
    Expression subtract(const Complex &val) const;

    Expression multiply(const Expression &ex) const;
    Expression multiply(double val) const;
    Expression multiply(const Complex &val) const;

    Expression divide(const Expression &ex) const;
    Expression divide(double val) const;
    Expression divide(const Complex &val) const;

    Expression negate() const;
    Expression invert() const;
    Expression conjugate() const;

    Expression pow(int n) const;
    Expression squared() const;
    Expression cubed() const;

    Expression pow(double n) const;
    Expression pow(const Complex &n) const;
    Expression pow(const Rational &n) const;
    Expression pow(const Expression &n) const;

    Expression root(int n) const;

    Expression substitute(const Identities &identities) const;
    Expression simplify(int options=0) const;

    QList<Factor> expand(const Symbol &symbol) const;

    double factor() const;

    double evaluate(const Values &values) const;
    Complex evaluate(const ComplexValues &values) const;

    Expression differentiate(const Symbol &symbol, int level=1) const;
    Expression integrate(const Symbol &symbol) const;

    Expression diff(const Symbol &symbol, int level=1) const;
    Expression integ(const Symbol &symbol) const;

    Expression series(const Symbol &symbol, int order) const;

    bool isZero() const;
    bool isConstant() const;
    bool isFunction(const Symbol &symbol) const;
    bool isCompound() const;
    bool isComplex() const;

    QString toString() const;
    QString toOpenMMString() const;

    uint hash() const
    {
        return exbase.hash();
    }

    const ExpressionBase& base() const;

    Symbols symbols() const;
    Functions functions() const;
    Expressions children() const;

    template<class T>
    QList<T> children() const;

private:

    /** The base of this expression */
    ExpressionBase exbase;

    /** The factor of the expression */
    double fac;
};

Expression operator+(const Expression &ex0, const Expression &ex1);
Expression operator+(const Expression &ex, double val);
Expression operator+(double val, const Expression &ex);
Expression operator+(const Expression &ex, const Complex &val);
Expression operator+(const Complex &val, const Expression &ex);
Expression operator-(const Expression &ex0, const Expression &ex1);
Expression operator-(const Expression &ex, double val);
Expression operator-(double val, const Expression &ex);
Expression operator*(const Expression &ex0, const Expression &ex1);
Expression operator*(double val, const Expression &ex);
Expression operator*(const Expression &ex, double val);
Expression operator*(const Complex &val, const Expression &ex);
Expression operator*(const Expression &ex, const Complex &val);
Expression operator/(const Expression &ex0, const Expression &ex1);
Expression operator/(const Expression &ex, double val);
Expression operator/(double val, const Expression &ex);
Expression operator/(const Expression &ex, const Complex &val);
Expression operator/(const Complex &val, const Expression &ex);
Expression pow(const Expression &ex0, int n);
Expression pow(const Expression &ex0, double n);
Expression pow(const Expression &ex0, const Expression &n);
Expression pow(const Expression &ex0, const Complex &n);
Expression pow(const Expression &ex0, const Rational &n);
Expression root(const Expression &ex0, int n);
Expression sqrt(const Expression &ex0);
Expression cbrt(const Expression &ex0);

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return a hash for an expression */
inline uint qHash(const Expression &ex)
{
    return ex.hash();
}

/** Return a list of all children of type 'T' in this expression */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QList<T> Expression::children() const
{
    Expressions exs = this->children();

    QList<T> children_t;

    for (Expressions::const_iterator it = exs.constBegin();
         it != exs.constEnd();
         ++it)
    {
        const ExpressionBase &base = it->base();

        //gccxml doesn't like this section, so remove it
        //when we are generating the python wrappers
        #ifndef SKIP_BROKEN_GCCXML_PARTS
        if (base.isA<T>())
            children_t.append( base.asA<T>() );
        #endif
    }

    return children_t;
}

/** Comparison operator */
inline bool Expression::operator==(const Expression &other) const
{
    return fac == other.fac and exbase == other.exbase;
}

/** Comparison operator */
inline bool Expression::operator!=(const Expression &other) const
{
    return fac != other.fac or exbase != other.exbase;
}

/** Addition operator */
inline Expression operator+(const Expression &ex0,
                            const Expression &ex1)
{
    return ex0.add(ex1);
}

/** Addition operator */
inline Expression operator+(const Expression &ex,
                            double val)
{
    return ex.add(val);
}

/** Addition operator */
inline Expression operator+(double val,
                            const Expression &ex)
{
    return ex.add(val);
}

/** Addition operator */
inline Expression operator+(const Expression &ex,
                            const Complex &val)
{
    return ex.add(val);
}

/** Addition operator */
inline Expression operator+(const Complex &val,
                            const Expression &ex)
{
    return ex.add(val);
}

/** Subtraction operator */
inline Expression operator-(const Expression &ex0,
                            const Expression &ex1)
{
    return ex0.subtract(ex1);
}

/** Subtraction operator */
inline Expression operator-(const Expression &ex,
                            double val)
{
    return ex.subtract(val);
}

/** Subtraction operator */
inline Expression operator-(double val,
                            const Expression &ex)
{
    return ex.negate().add(val);
}


/** Multiplication operator */
inline Expression operator*(const Expression &ex0,
                            const Expression &ex1)
{
    return ex0.multiply(ex1);
}

/** Multiplication operator */
inline Expression operator*(double val, const Expression &ex)
{
    return ex.multiply(val);
}

/** Multiplication operator */
inline Expression operator*(const Expression &ex, double val)
{
    return ex.multiply(val);
}

/** Multiplication operator */
inline Expression operator*(const Complex &val, const Expression &ex)
{
    return ex.multiply(val);
}

/** Multiplication operator */
inline Expression operator*(const Expression &ex, const Complex &val)
{
    return ex.multiply(val);
}

/** Division operator */
inline Expression operator/(const Expression &ex0,
                            const Expression &ex1)
{
    return ex0.divide(ex1);
}

/** Division operator */
inline Expression operator/(const Expression &ex,
                            double val)
{
    return ex.divide(val);
}

/** Division operator */
inline Expression operator/(double val,
                            const Expression &ex)
{
    return ex.invert().multiply(val);
}

/** Division operator */
inline Expression operator/(const Expression &ex,
                            const Complex &val)
{
    return ex.divide(val);
}

/** Division operator */
inline Expression operator/(const Complex &val,
                            const Expression &ex)
{
    return ex.invert().multiply(val);
}

/** Raise an expression to the nth power */
inline Expression pow(const Expression &ex0, int n)
{
    return ex0.pow(n);
}

/** Raise an expression to a real power */
inline Expression pow(const Expression &ex0, double n)
{
    return ex0.pow(n);
}

/** Raise an expression to a functional power */
inline Expression pow(const Expression &ex0,
                      const Expression &n)
{
    return ex0.pow(n);
}

/** Raise an expression to a complex power */
inline Expression pow(const Expression &ex0, const Complex &n)
{
    return ex0.pow(n);
}

/** Raise an expression to a rational power */
inline Expression pow(const Expression &ex0, const Rational &n)
{
    return ex0.pow(n);
}

/** Take the nth root of an expression */
inline Expression root(const Expression &ex0, int n)
{
    return ex0.root(n);
}

/** Take the square root of an expression */
inline Expression sqrt(const Expression &ex0)
{
    return ex0.root(2);
}

/** Take the cube root of an expression */
inline Expression cbrt(const Expression &ex0)
{
    return ex0.root(3);
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE(SireCAS::Expression)

SIRE_EXPOSE_CLASS( SireCAS::Expression )
SIRE_EXPOSE_FUNCTION( SireCAS::pow )
SIRE_EXPOSE_FUNCTION( SireCAS::sqrt )
SIRE_EXPOSE_FUNCTION( SireCAS::cbrt )

SIRE_END_HEADER

#endif
