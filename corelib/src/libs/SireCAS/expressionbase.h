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

#ifndef SIRECAS_EXPRESSIONBASE_H
#define SIRECAS_EXPRESSIONBASE_H

#include <QString>

#include "exbase.h"

#include "SireBase/sharedpolypointer.hpp"

SIRE_BEGIN_HEADER

namespace SireMaths
{
class Complex;
}

namespace SireCAS
{
class ExpressionBase;
}

class QDataStream;
SIRECAS_EXPORT QDataStream& operator<<(QDataStream&, const SireCAS::ExpressionBase&);
SIRECAS_EXPORT QDataStream& operator>>(QDataStream&, SireCAS::ExpressionBase&);

namespace SireCAS
{

class ExBase;
class Expression;
class Expressions;
class Identities;
class Symbol;
class Symbols;
class Functions;
class Values;
class ComplexValues;

using SireMaths::Complex;
using SireBase::SharedPolyPointer;

/** This class provides implicitly shared access to ExBase objects, which
    are a virtual class hierarchy that provide all of the functionality of
    the expressions.

    @author Christopher Woods
*/
class SIRECAS_EXPORT ExpressionBase
{

friend SIRECAS_EXPORT QDataStream& ::operator<<(QDataStream&, const ExpressionBase&);
friend SIRECAS_EXPORT QDataStream& ::operator>>(QDataStream&, ExpressionBase&);

public:
    ExpressionBase();

    ExpressionBase(const ExBase &ex);

    ExpressionBase(const ExpressionBase &other);

    ~ExpressionBase();

    bool operator==(const ExpressionBase &other) const;
    bool operator!=(const ExpressionBase &other) const;

    ExpressionBase& operator=(const ExpressionBase &other);
    ExpressionBase& operator=(const ExBase &other);

    Expression differentiate(const Symbol &symbol) const;
    Expression integrate(const Symbol &symbol) const;

    Expression series(const Symbol &symbol, int n) const;

    Expression simplify(int options=0) const;

    QList<Factor> expand(const Symbol &symbol) const;

    Expression conjugate() const;

    bool isFunction(const Symbol&) const;
    bool isConstant() const;
    bool isComplex() const;
    bool isCompound() const;

    uint hash() const;

    const char* what() const;
    QString toString() const;
    QString toOpenMMString() const;

    double evaluate(const Values &values) const;
    Complex evaluate(const ComplexValues &values) const;

    Expression substitute(const Identities &identities) const;

    Symbols symbols() const;
    Functions functions() const;
    Expressions children() const;

    template<class T>
    bool isA() const
    {
        return d->isA<T>();
    }

    template<class T>
    const T& asA() const
    {
        return d->asA<T>();
    }

private:
    /** Shared pointer to the expression */
    SharedPolyPointer<ExBase> d;
};

/** Return a hash for an ExpressionBase */
inline uint qHash(const ExpressionBase &ex)
{
    return ex.hash();
}

}

Q_DECLARE_METATYPE(SireCAS::ExpressionBase)

SIRE_EXPOSE_CLASS( SireCAS::ExpressionBase )

SIRE_END_HEADER

#endif
