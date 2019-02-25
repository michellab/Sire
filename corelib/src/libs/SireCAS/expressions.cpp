/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#include "expression.h"
#include "expressions.h"

using namespace SireCAS;

Expressions::Expressions() : QList<Expression>()
{}

Expressions::Expressions(const Expression &expression) : QList<Expression>()
{
    this->append(expression);
}

Expressions::Expressions(const QList<Expression> &expressions)
            : QList<Expression>(expressions)
{}

Expressions::~Expressions()
{}

Expressions Expressions::differentiate(const Symbol &symbol) const
{
    Expressions diffs;

    int sz = count();

    for (int i=0; i<sz; ++i)
    {
        Expression diff = at(i).differentiate(symbol);

        if (not diff.isZero())
            diffs.append(diff);
    }

    return diffs;
}

Expressions Expressions::integrate(const Symbol &symbol) const
{
    Expressions ints;

    int sz = count();

    for (int i=0; i<sz; ++i)
    {
        Expression integ = at(i).integrate(symbol);

        if (not integ.isZero())
            ints.append(integ);
    }

    return ints;
}
