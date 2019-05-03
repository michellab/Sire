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

#include "constant.h"
#include "expression.h"
#include "expressions.h"
#include "symbol.h"
#include "symbols.h"
#include "values.h"
#include "complexvalues.h"
#include "identities.h"
#include "functions.h"
#include "integrationconstant.h"

#include "SireStream/datastream.h"

using namespace SireStream;
using namespace SireCAS;

static const RegisterMetaType<Constant> r_constant;

/** Serialise a constant to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const Constant &constant)
{
    writeHeader(ds, r_constant, 1)
            << static_cast<const ExBase&>(constant);

    return ds;
}

/** Deserialise a constant from a binary datastream */
QDataStream &operator>>(QDataStream &ds, Constant &constant)
{
    VersionID v = readHeader(ds, r_constant);

    if (v == 1)
    {
        ds >> static_cast<ExBase&>(constant);
    }
    else
        throw version_error(v, "1", r_constant, CODELOC);

    return ds;
}

/** Construct a constant */
Constant::Constant() : ExBase()
{}

/** Copy constructor */
Constant::Constant(const Constant&) : ExBase()
{}

/** Destructor */
Constant::~Constant()
{}

/** Differential of a constant is zero */
Expression Constant::differentiate(const Symbol&) const
{
    return 0;
}

/** Integral of a constant is = constant*symbol + C */
Expression Constant::integrate(const Symbol &symbol) const
{
    return *this * symbol;
}

/** Comparison operator */
bool Constant::operator==(const ExBase &other) const
{
    return typeid(other).name() == typeid(*this).name();
}

/** Hash a constant */
uint Constant::hash() const
{
    return (r_constant.magicID() << 16) | (r_constant.magicID() & 0x0000FFFF);
}

/** Return a string representation of this constant (actually an empty string!) */
QString Constant::toString() const
{
    return QString("1.0");
}

/** Evaluation of a constant is 1 */
double Constant::evaluate(const Values&) const
{
    return 1;
}

/** Evaluation of a constant is 1 */
Complex Constant::evaluate(const ComplexValues&) const
{
    return Complex(1,0);
}

/** Can't substitute into a constant */
Expression Constant::substitute(const Identities&) const
{
    return Expression(*this);
}

/** No symbols in a constant */
Symbols Constant::symbols() const
{
    return Symbols();
}

/** No functions in a constant */
Functions Constant::functions() const
{
    return Functions();
}

/** No children in a constant */
Expressions Constant::children() const
{
    return Expressions();
}

QList<Factor> Constant::expand(const Symbol &symbol) const
{
    QList<Factor> ret;
    ret.append( Factor( symbol, 1, 0 ) );
    return ret;
}

const char* Constant::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Constant>() );
}

Constant* Constant::clone() const
{
    return new Constant(*this);
}

