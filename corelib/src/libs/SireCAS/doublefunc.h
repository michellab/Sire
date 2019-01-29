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

#ifndef SIRECAS_DOUBLEFUNC_H
#define SIRECAS_DOUBLEFUNC_H

#include "exbase.h"
#include "symbol.h"
#include "symbols.h"
#include "functions.h"
#include "expression.h"
#include "expressions.h"

SIRE_BEGIN_HEADER

namespace SireCAS
{
class DoubleFunc;
}

SIRECAS_EXPORT QDataStream& operator<<(QDataStream&, const SireCAS::DoubleFunc&);
SIRECAS_EXPORT QDataStream& operator>>(QDataStream&, SireCAS::DoubleFunc&);

namespace SireCAS
{

/** Base class of all double-expression functions ( e.g. f(x(), y()) )

    @author Christopher Woods
*/
class SIRECAS_EXPORT DoubleFunc : public ExBase
{

friend QDataStream& ::operator<<(QDataStream&, const DoubleFunc&);
friend QDataStream& ::operator>>(QDataStream&, DoubleFunc&);

public:
    DoubleFunc();
    DoubleFunc(const Expression &x, const Expression &y);

    DoubleFunc(const DoubleFunc &other);

    ~DoubleFunc();

    DoubleFunc& operator=(const DoubleFunc &other);

    uint hash() const;

    const Expression& x() const;
    const Expression& y() const;

    Expression conjugate() const;

    bool isFunction(const Symbol &symbol) const;
    bool isConstant() const;
    bool isComplex() const;
    bool isCompound() const;

    QString toString() const;

    Expression substitute(const Identities &identities) const;
    Symbols symbols() const;
    Functions functions() const;
    Expressions children() const;

    Expression differentiate(const Symbol &symbol) const;
    Expression integrate(const Symbol &symbol) const;

    QList<Factor> expand(const Symbol &symbol) const;

protected:

    virtual Expression functionOf(const Expression &x, 
                                  const Expression &y) const=0;

    virtual QString stringRep() const=0;

    virtual uint magic() const=0;

    /** The two expressions that this function operates on */
    Expression ex0, ex1;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the first argument - viewed as this is f( x(), y() ) */
inline const Expression& DoubleFunc::x() const
{
    return ex0;
}

/** Return the second argument - viewed as this is f( x(), y() ) */
inline const Expression& DoubleFunc::y() const
{
    return ex1;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

/** To declare a new function copy Min or Max (in minmax.h) */

}

SIRE_EXPOSE_CLASS( SireCAS::DoubleFunc )

SIRE_END_HEADER

#endif
