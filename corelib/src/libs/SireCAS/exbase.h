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

#ifndef SIRECAS_EXBASE_H
#define SIRECAS_EXBASE_H

#include <QString>

#include "SireMaths/rational.h"
#include "SireMaths/complex.h"

#include "SireBase/refcountdata.h"

SIRE_BEGIN_HEADER

namespace SireCAS
{
class ExBase;
}

SIRECAS_EXPORT QDataStream& operator<<(QDataStream&, const SireCAS::ExBase&);
SIRECAS_EXPORT QDataStream& operator>>(QDataStream&, SireCAS::ExBase&);

namespace SireCAS
{

class Values;
class ComplexValues;
class Symbol;
class Symbols;
class Functions;
class Expression;
class Expressions;
class Identities;
class Factor;

class ExpressionBase;
class Expression;

using SireMaths::Complex;
using SireMaths::Rational;

/**
Pure-virtual base class of all of the parts of mathematical expressions.

This class provides the 'atom' of SireCAS. ExBase objects are combined together
to form complete expressions. All constants, functions and symbols are derived
from this object.

This class is an example of an implicitly shared, self-managed object, and
is designed so that it can be held by SharedPolyPointer (indeed,
ExpressionBase is just a proxy for SharedPolyPointer<ExBase>).

@author Christopher Woods
*/
class SIRECAS_EXPORT ExBase : public SireBase::RefCountData
{

friend QDataStream& ::operator<<(QDataStream&, const ExBase&);
friend QDataStream& ::operator>>(QDataStream&, ExBase&);

public:
    typedef ExBase ROOT;

    ExBase();

    ExBase(const ExBase &other);

    virtual ~ExBase();

    ///////
    /////// Non-virtual functions
    ///////

    /** Assignment operator */
    ExBase& operator=(const ExBase&)
    {
        return *this;
    }

    bool operator!=(const ExBase &other) const;
    Expression operator-() const;

    /** Return whether or not this is an ExBase of type 'T' */
    template<class T>
    bool isA() const
    {
        return dynamic_cast<const T*>(this) != 0;
    }

    /** Return this ExBase cast as a const reference to type 'T'.
        Note that this is only safe if 'isA<T>()' returns true. */
    template<class T>
    const T& asA() const
    {
        return dynamic_cast<const T&>(*this);
    }

    ///////
    /////// Virtual functions - you may wish to override these
    /////// in your derived class
    ///////

    virtual Expression differentiate(const Symbol &symbol) const;
    virtual Expression integrate(const Symbol &symbol) const;

    virtual Expression series(const Symbol &symbol, int n) const;

    virtual Expression simplify(int options=0) const;

    virtual Expression conjugate() const;

    virtual bool isFunction(const Symbol&) const;
    virtual bool isConstant() const;
    virtual bool isComplex() const;
    virtual bool isCompound() const;

    /** Return a string representation of this object in the OpenMM syntax*/
    virtual QString toOpenMMString() const;

    ///////
    /////// Pure-virtual functions - these must be overridden
    /////// in your derived class
    ///////

    /** Comparison operator - only return true if these are the same
        class and contain the same data. */
    virtual bool operator==(const ExBase &other) const=0;

    /** Return a hash of this object - return a combination of the
        identifying magic for the class and a hash for its contents. */
    virtual uint hash() const=0;

    /** Return the name of the type of this ExBase object */
    virtual const char* what() const=0;

    /** Return the name of this class type */
    static const char* typeName()
    {
        return "SireCAS::ExBase";
    }

    /** Return a string representation of this object */
    virtual QString toString() const=0;


    /** Return a clone of this object */
    virtual ExBase* clone() const=0;

    /** Evaluate this ExBase using values 'values'. Any
        missing symbols are assumed to equal zero.

        Note that an exception will be thrown if the result of the
        evaluation of this, or one of its children, is complex.

        \throw SireMaths::domain_error
    */
    virtual double evaluate(const Values &values) const=0;

    /** Evaluate this ExBase using the complex values 'values'.
        Any missing symbols are assumed to equal zero. */
    virtual Complex evaluate(const ComplexValues &values) const=0;

    /** Return an expression that has the identities in 'identities'
        substituted into this expression */
    virtual Expression substitute(const Identities &identities) const=0;

    /** Return the set of Symbols that appear in this ExBase */
    virtual Symbols symbols() const=0;

    /** Return the set of Functions that appear in this ExBase */
    virtual Functions functions() const=0;

    /** Return the child expressions of this Expression */
    virtual Expressions children() const=0;
    
    /** Rearrange this expression into the form
        m x^i + n x^j + ... + constant
        and return the factors and powers of x
        
        \throw SireCAS::rearrangement_error
    */
    virtual QList<Factor> expand(const Symbol &symbol) const=0;

};

/** Return a hash of an ExBase object */
inline uint qHash(const ExBase &ex)
{
    return ex.hash();
}

SIRECAS_EXPORT Expression operator+(const ExBase &base0, const ExBase &base1);
SIRECAS_EXPORT Expression operator+(const ExBase &base, const Expression &ex);
SIRECAS_EXPORT Expression operator+(const Expression &ex, const ExBase &base);
SIRECAS_EXPORT Expression operator+(const ExBase &base, double val);
SIRECAS_EXPORT Expression operator+(double val, const ExBase &base);
SIRECAS_EXPORT Expression operator+(const ExBase &base, const Complex &val);
SIRECAS_EXPORT Expression operator+(const Complex &val, const ExBase &base);

SIRECAS_EXPORT Expression operator-(const ExBase &base0, const ExBase &base1);
SIRECAS_EXPORT Expression operator-(const ExBase &base, const Expression &ex);
SIRECAS_EXPORT Expression operator-(const Expression &ex, const ExBase &base);
SIRECAS_EXPORT Expression operator-(const ExBase &base, double val);
SIRECAS_EXPORT Expression operator-(double val, const ExBase &base);
SIRECAS_EXPORT Expression operator-(const ExBase &base, const Complex &val);
SIRECAS_EXPORT Expression operator-(const Complex &val, const ExBase &base);

SIRECAS_EXPORT Expression operator*(const ExBase &base0, const ExBase &base1);
SIRECAS_EXPORT Expression operator*(const ExBase &base, const Expression &ex);
SIRECAS_EXPORT Expression operator*(const Expression &ex, const ExBase &base);
SIRECAS_EXPORT Expression operator*(const ExBase &base, double val);
SIRECAS_EXPORT Expression operator*(double val, const ExBase &base);
SIRECAS_EXPORT Expression operator*(const ExBase &base, const Complex &val);
SIRECAS_EXPORT Expression operator*(const Complex &val, const ExBase &base);

SIRECAS_EXPORT Expression operator/(const ExBase &base0, const ExBase &base1);
SIRECAS_EXPORT Expression operator/(const ExBase &base, const Expression &ex);
SIRECAS_EXPORT Expression operator/(const Expression &ex, const ExBase &base);
SIRECAS_EXPORT Expression operator/(const ExBase &base, double val);
SIRECAS_EXPORT Expression operator/(double val, const ExBase &base);
SIRECAS_EXPORT Expression operator/(const ExBase &base, const Complex &val);
SIRECAS_EXPORT Expression operator/(const Complex &val, const ExBase &base);

SIRECAS_EXPORT Expression pow(const ExBase &base, int n);
SIRECAS_EXPORT Expression pow(const ExBase &base, const Rational &n);
SIRECAS_EXPORT Expression pow(const ExBase &base, double n);
SIRECAS_EXPORT Expression pow(const ExBase &base, const Complex &n);
SIRECAS_EXPORT Expression pow(const ExBase &base, const Expression &n);
SIRECAS_EXPORT Expression pow(const ExBase &base, const ExBase &n);

}

SIRE_EXPOSE_CLASS( SireCAS::ExBase )
SIRE_EXPOSE_FUNCTION( SireCAS::pow )

SIRE_END_HEADER

#endif
