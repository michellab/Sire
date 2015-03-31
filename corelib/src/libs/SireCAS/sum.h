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

#ifndef SIRECAS_SUM_H
#define SIRECAS_SUM_H

#include <QHash>

#include "exbase.h"
#include "expression.h"

SIRE_BEGIN_HEADER

namespace SireCAS
{
class Sum;
}

QDataStream& operator<<(QDataStream&, const SireCAS::Sum&);
QDataStream& operator>>(QDataStream&, SireCAS::Sum&);

namespace SireCAS
{

class Product;

/**
This class holds a collection of expressions that are to be added (or subtracted) from one another

@author Christopher Woods
*/
class SIRECAS_EXPORT Sum : public ExBase
{

friend QDataStream& ::operator<<(QDataStream&, const Sum&);
friend QDataStream& ::operator>>(QDataStream&, Sum&);

public:
    Sum();
    Sum(const Expression &ex0, const Expression &ex1);
    Sum(const Expressions &expressions);

    Sum(const Sum &other);

    ~Sum();

    Expression differentiate(const Symbol &symbol) const;
    Expression integrate(const Symbol &symbol) const;

    Expression series(const Symbol &symbol, int n) const;

    Expression simplify(int options=0) const;

    Expression conjugate() const;

    bool isFunction(const Symbol&) const;
    bool isConstant() const;
    bool isComplex() const;
    bool isCompound() const;

    bool operator==(const ExBase &other) const;
    uint hash() const;

    static const char* typeName();

    const char* what() const
    {
        return Sum::typeName();
    }

    QString toString() const;
    QString toOpenMMString() const;

    double evaluate(const Values &values) const;
    Complex evaluate(const ComplexValues &values) const;

    Expression substitute(const Identities &identities) const;

    Symbols symbols() const;
    Functions functions() const;
    Expressions children() const;

    Expression reduce() const;

    Sum* clone() const;

    QList<Factor> expand(const Symbol &symbol) const;

private:
    friend class Product;

    void add(const Expression &ex);

    Expression take(const ExpressionBase &ex);
    void add(double fac, const ExpressionBase &ex);

    /** The positive parts of the sum */
    QHash<ExpressionBase, Expression> posparts;

    /** The negative parts of the sum */
    QHash<ExpressionBase, Expression> negparts;

    /** The start value of the sum */
    double strtval;
};

}

Q_DECLARE_METATYPE(SireCAS::Sum)

SIRE_EXPOSE_CLASS( SireCAS::Sum )

SIRE_END_HEADER

#endif
