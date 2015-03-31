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

#ifndef SIRECAS_FUNCTION_H
#define SIRECAS_FUNCTION_H

#include <QHash>
#include <QSet>

#include "symbol.h"
#include "functionsignature.h"

SIRE_BEGIN_HEADER

namespace SireCAS
{
class FunctionPvt;
class Function;
}

QDataStream& operator<<(QDataStream&, const SireCAS::FunctionPvt&);
QDataStream& operator>>(QDataStream&, SireCAS::FunctionPvt&);

QDataStream& operator<<(QDataStream&, const SireCAS::Function&);
QDataStream& operator>>(QDataStream&, SireCAS::Function&);

namespace SireCAS
{

class Symbols;

/** This is a private class that is used to hold the data of
    Function objects.

    \author Christopher Woods
*/
class FunctionPvt
{

friend QDataStream& ::operator<<(QDataStream&, const FunctionPvt&);
friend QDataStream& ::operator>>(QDataStream&, FunctionPvt&);

public:
    FunctionPvt();

    FunctionPvt(const QString &name);

    FunctionPvt(const FunctionPvt &other);

    ~FunctionPvt();

    bool operator==(const FunctionPvt &other) const;
    bool operator!=(const FunctionPvt &other) const;

    void add(const Symbol &symbol);

    QString toString() const;

    QString name() const
    {
        return sig.name();
    }

    const FunctionSignature& signature() const
    {
        return sig;
    }

    const QHash<Symbol,int>& symbols() const
    {
        return syms;
    }

    const QSet<Function>& functions() const
    {
        return funcs;
    }

    void differentiate(const Symbol &symbol);
    void integrate(const Symbol &symbol);

private:

    /** The signature of the function */
    FunctionSignature sig;

    /** All of the symbols in the function */
    QHash<Symbol, int> syms;

    /** All of the functions of the function */
    QSet<Function> funcs;
};

/**
This class represents a placeholder for a generic function (e.g. f(x)).
This is meant to be used when the exact form of the function is unknown, and
the values of the function are supplied during normal evaluation, e.g.

\code
Symbol x("x");

Function f("f");

Expression ex = 3*x + f(x);

double result = ex.evaluate( Values( x==0.4, f(x)==0.2 ) );
\endcode

In this example, the expression equals 3x + f(x), and the result
for x==0.4 and f(x)==0.2 is 1.4.

In many ways, a function is just like a symbol. However, there is one
very important difference. The differential of a symbol is a constant,
while the differential of a function is another function, e.g.;

\code
Symbol x("x");
Function f("f");

Expression ex = 3*x + f(x);

Expression dx = ex.differentiate(x);

//dx == 3 + f(x')

//evaluate the differential
double val = dx.evaluate( Values( x==0.4, f(x).diff(x) == 0.2 ) );

Expression intx = ex.integrate(x);

//dx == (3/2)x^2 + f(x`) + C

//evaluate the integral for C == 0.4
double val = intx.evaluate( Values( x==0.4, f(x).integ(x)==0.2, C==0.4 ) );

\endcode

The differential of f(x) has returned f(x') (the prime indicates that the
function has been differentiated with respect to x). Similarly, the integration
of f(x) has resulted in f(x`) (the mark indicates that the function has been
integrated with respect to x). Now the differentials of integrals of this function
can be substituted directly when evaluating the expression.

The use of primes and marks on the symbols is a bit weird, but it in response
to the need to be able to represent complex patterns of differentiation and
integration, e.g. f(x,y') implies f(x,y) has been differentiated with respect
to y, f(x',y') imples differentiation with respect to x and y, f(x'',y) implies
double differentiation with respect to x etc. Similarly, f(x`,y) implies
integration with respect to x, f(x`,y`) implies integration with respect to
x and y, and f(x``,y) implies double integration with respect to x.

Just as symbols can be substituted in an expression, so to can functions. Functions
however provide an additional guard that the substituting expression is truly
an expression of the symbols of the function, and also provides automatic
substitution of the differentials and integrals, e.g.

\code
Symbol x("x");
Function f("f");

Expression ex = 3*x + f(x);

Expression ex2 = ex.substitute( f(x) = Sin(x) );

//ex2 now equals 3x + sin(x)

Expression dx = ex.differentiate(x);

//dx now equals 3 + f(x')

Expression dx2 = dx.substitute( f(x) = Sin(x) );

// dx2 now equals 3 + cos(x)  - note how f(x) has been automatically
// differentiated upon substitution.
\endcode

You can also try to make a function of a function. The code will, however,
expand this out to a new function, e.g.

\code
Symbol x("x");
Function f("f");
Function g("g");

Expression ex = 3*x + f( g(x) );

// ex now equals 3x + f_g(x)   - note how f( g(x) ) has been expanded
//                               into a single combined function.

Expression dx = ex.differentiate(x);

// dx now equals 3 + f_g(x')

double val = ex.evalute( Values( x==0.3, f(g(x))==5.4 ) );

// val == 6.3

Expression dx2 = dx.substitute( f(g(x)) == Sin(x) );

// dx2 now equals 3 + cos(x)

\endcode

Note that it makes no sense to substitute for only g(x) in the above example,
as it cannot be known what the functional form of f(g(x)) will be. Also note that
the code allow you to create repeated functions (e.g. f( g( f(x) ) ), == f_g_f(x) ).

Functions can also have their dependent symbols extended, e.g.

\code
Symbol x("x"), y("y");

Function f("f");

Function fx = f(x);

//fx == f(x)

Function fxy = fx(y);

//fxy == f(x,y)
\endcode

@author Christopher Woods
*/
class SIRECAS_EXPORT Function : public Symbol
{

friend class FunctionPvt;
friend QDataStream& ::operator<<(QDataStream&, const Function&);
friend QDataStream& ::operator>>(QDataStream&, Function&);

public:
    Function();
    Function(const QString &name);

    Function(const QString &name, const Symbols &symbols);

    Function(const QString &name,
             const Symbol &sym0);
    Function(const QString &name,
             const Symbol &sym0, const Symbol &sym1);
    Function(const QString &name,
             const Symbol &sym0, const Symbol &sym1, const Symbol &sym2);
    Function(const QString &name,
             const Symbol &sym0, const Symbol &sym1, const Symbol &sym2,
             const Symbol &sym3);
    Function(const QString &name,
             const Symbol &sym0, const Symbol &sym1, const Symbol &sym2,
             const Symbol &sym3, const Symbol &sym4);
    Function(const QString &name,
             const Symbol &sym0, const Symbol &sym1, const Symbol &sym2,
             const Symbol &sym3, const Symbol &sym4, const Symbol &sym5);
    Function(const QString &name,
             const Symbol &sym0, const Symbol &sym1, const Symbol &sym2,
             const Symbol &sym3, const Symbol &sym4, const Symbol &sym5,
             const Symbol &sym6);
    Function(const QString &name,
             const Symbol &sym0, const Symbol &sym1, const Symbol &sym2,
             const Symbol &sym3, const Symbol &sym4, const Symbol &sym5,
             const Symbol &sym6, const Symbol &sym7);
    Function(const QString &name,
             const Symbol &sym0, const Symbol &sym1, const Symbol &sym2,
             const Symbol &sym3, const Symbol &sym4, const Symbol &sym5,
             const Symbol &sym6, const Symbol &sym7, const Symbol &sym8);
    Function(const QString &name,
             const Symbol &sym0, const Symbol &sym1, const Symbol &sym2,
             const Symbol &sym3, const Symbol &sym4, const Symbol &sym5,
             const Symbol &sym6, const Symbol &sym7, const Symbol &sym8,
             const Symbol &sym9);

    Function(const Function &other);

    ~Function();

    bool operator==(const ExBase &other) const;

    Function& operator=(const Function &other);

    /** Convienient operator used to combine a function with a value */
    SymbolValue operator==(double val) const
    {
        return SymbolValue(ID(), val);
    }

    /** Convienient operator used to combine a function with a value */
    SymbolValue operator==(int val) const
    {
        return SymbolValue(ID(), val);
    }

    /** Convienient operator used to combine a function with a Complex */
    SymbolComplex operator==(const Complex &val) const
    {
        return SymbolComplex(ID(), val);
    }

    /** Convieient operator used to combine a function with an equivalent
        expression */
    SymbolExpression operator==(const Expression &ex) const
    {
        return SymbolExpression(*this, ex);
    }

    Expression differentiate(const Symbol &symbol) const;
    Expression integrate(const Symbol &symbol) const;

    bool isFunction(const Symbol&) const;

    /** Return the name of this function */
    QString name() const
    {
        return d.name();
    }

    /** Return the signature of this function */
    FunctionSignature signature() const
    {
        return d.signature();
    }

    uint hash() const;

    static const char* typeName();

    const char* what() const
    {
        return Function::typeName();
    }

    Function* clone() const;

    Expression substitute(const Identities &identities) const;

    Symbols symbols() const;
    Functions functions() const;

    Function operator()(const Symbols &symbols) const;

    Function operator()(const Symbol &sym0) const;
    Function operator()(const Symbol &sym0, const Symbol &sym1) const;
    Function operator()(const Symbol &sym0, const Symbol &sym1, const Symbol &sym2) const;
    Function operator()(const Symbol &sym0, const Symbol &sym1, const Symbol &sym2,
                        const Symbol &sym3) const;
    Function operator()(const Symbol &sym0, const Symbol &sym1, const Symbol &sym2,
                        const Symbol &sym3, const Symbol &sym4) const;
    Function operator()(const Symbol &sym0, const Symbol &sym1, const Symbol &sym2,
                        const Symbol &sym3, const Symbol &sym4, const Symbol &sym5) const;
    Function operator()(const Symbol &sym0, const Symbol &sym1, const Symbol &sym2,
                        const Symbol &sym3, const Symbol &sym4, const Symbol &sym5,
                        const Symbol &sym6) const;
    Function operator()(const Symbol &sym0, const Symbol &sym1, const Symbol &sym2,
                        const Symbol &sym3, const Symbol &sym4, const Symbol &sym5,
                        const Symbol &sym6, const Symbol &sym7) const;
    Function operator()(const Symbol &sym0, const Symbol &sym1, const Symbol &sym2,
                        const Symbol &sym3, const Symbol &sym4, const Symbol &sym5,
                        const Symbol &sym6, const Symbol &sym7, const Symbol &sym8) const;
    Function operator()(const Symbol &sym0, const Symbol &sym1, const Symbol &sym2,
                        const Symbol &sym3, const Symbol &sym4, const Symbol &sym5,
                        const Symbol &sym6, const Symbol &sym7, const Symbol &sym8,
                        const Symbol &sym9) const;

    QList<Factor> factorise(const Symbol &symbol) const;

private:

    Function(const FunctionPvt &other);

    Expression match(const Function &func, const Expression &ex) const;

    const FunctionPvt& data() const
    {
        return d;
    }

    /** Class holding the metadata for this function */
    FunctionPvt d;
};

inline uint qHash(const Function &func)
{
    return func.hash();
}

}

Q_DECLARE_METATYPE(SireCAS::Function)

SIRE_END_HEADER

#endif
