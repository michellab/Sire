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

#include "exp.h"
#include "values.h"
#include "complexvalues.h"
#include "integrationconstant.h"

#include "SireStream/datastream.h"

using namespace SireStream;
using namespace SireCAS;

/////////////
///////////// Implementation of Exp
/////////////

static const RegisterMetaType<Exp> r_exp;

/** Serialise an Exp to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const Exp &exp)
{
    writeHeader(ds, r_exp, 1)
          << exp.pwr << static_cast<const PowerFunction&>(exp);

    return ds;
}

/** Deserialise an Exp from a binary datastream */
QDataStream &operator>>(QDataStream &ds, Exp &exp)
{
    VersionID v = readHeader(ds, r_exp);

    if (v == 1)
    {
        ds >> exp.pwr >> static_cast<PowerFunction&>(exp);
    }
    else
        throw version_error(v, "1", r_exp, CODELOC);

    return ds;
}

/** Construct an empty Exp (e^0) */
Exp::Exp() : PowerFunction()
{}

/** Construct e^power */
Exp::Exp(const Expression &power) : PowerFunction(), pwr(power)
{}

/** Copy constructor */
Exp::Exp(const Exp &other) : PowerFunction(), pwr(other.pwr)
{}

/** Destructor */
Exp::~Exp()
{}

/** Comparison operator */
bool Exp::operator==(const ExBase &other) const
{
    const Exp *other_exp = dynamic_cast<const Exp*>(&other);

    return other_exp != 0 and typeid(other).name() == typeid(*this).name()
             and pwr == other_exp->pwr;
}

/** Return a hash of this Exp function */
uint Exp::hash() const
{
    return ( r_exp.magicID() << 16) | (pwr.hash() & 0x0000FFFF);
}

/** Return a string representation of this function */
QString Exp::toString() const
{
    if (pwr.isZero())
        return "1";
    else
        return QString("exp(%1)").arg(pwr.toString());
}

/** Evaluate this expression at 'values' */
double Exp::evaluate(const Values &values) const
{
    return std::exp( pwr.evaluate(values) );
}

/** Complex evaluate */
Complex Exp::evaluate(const ComplexValues &values) const
{
    return SireMaths::exp( pwr.evaluate(values) );
}

/** Differentiatial of e^x is e^x */
Expression Exp::differentiate(const Symbol &symbol) const
{
    // d e^f(x) / dx = f'(x) * e^f(x)
    return pwr.differentiate(symbol) * Expression(*this);
}

/** Integral of e^x is e^x + c */
Expression Exp::integrate(const Symbol &symbol) const
{
    return Exp(pwr) / pwr.diff(symbol);
}

/** Return the core of this power */
Expression Exp::core() const
{
    return Expression( SireMaths::e );
}

/** Return the power of this power */
Expression Exp::power() const
{
    return pwr;
}

const char* Exp::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Exp>() );
}

/////////////
///////////// Implementation of Ln
/////////////

static const RegisterMetaType<Ln> r_ln;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const Ln &ln)
{
    writeHeader(ds, r_ln, 1) << static_cast<const SingleFunc&>(ln);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds, Ln &ln)
{
    VersionID v = readHeader(ds, r_ln);

    if (v == 1)
    {
        ds >> static_cast<SingleFunc&>(ln);
    }
    else
        throw version_error(v, "1", r_ln, CODELOC);

    return ds;
}

/** Construct an empty ln (ln(0)) */
Ln::Ln() : SingleFunc()
{}

/** Construct ln(expression) */
Ln::Ln(const Expression &expression) : SingleFunc(expression)
{}

/** Copy constructor */
Ln::Ln(const Ln &other) : SingleFunc(other)
{}

/** Destructor */
Ln::~Ln()
{}

/** Comparison operator */
bool Ln::operator==(const ExBase &other) const
{
    const Ln *other_ln = dynamic_cast<const Ln*>(&other);
    return other_ln != 0 and typeid(other).name() == typeid(*this).name()
             and this->argument() == other_ln->argument();
}

/** Return the magic for this function */
uint Ln::magic() const
{
    return r_ln.magicID();
}

/** Evaluate this function */
double Ln::evaluate(const Values &values) const
{
    return std::log( x().evaluate(values) );
}

/** Complex evaluation */
Complex Ln::evaluate(const ComplexValues &values) const
{
    return SireMaths::log( x().evaluate(values) );
}

/** differential of ln(x) = 1/x */
Expression Ln::diff() const
{
    // d ln(f(x)) / dx = f'(x) / f(x)
    return 1 / x();
}

/** Integral of ln(x) = x ln(x) - x + C */
Expression Ln::integ() const
{
    return x()*Ln(x()) - x();
}

const char* Ln::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Ln>() );
}

Ln* Ln::clone() const
{
    return new Ln(*this);
}


Exp* Exp::clone() const
{
    return new Exp(*this);
}

