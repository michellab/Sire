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

#include "exbase.h"
#include "values.h"
#include "identities.h"
#include "complexvalues.h"
#include "functions.h"
#include "expressionbase.h"
#include "expression.h"

#include "SireCAS/errors.h"

#include "SireStream/datastream.h"

using namespace SireCAS;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<ExBase> r_exbase(MAGIC_ONLY, "SireCAS::ExBase");

/** Serialise to a binary datastream */
QDataStream SIRECAS_EXPORT &operator<<(QDataStream &ds, const ExBase&)
{
    writeHeader(ds, r_exbase, 0);
    return ds;
}

/** Deserialise from a binary datastream */
QDataStream SIRECAS_EXPORT &operator>>(QDataStream &ds, ExBase&)
{
    VersionID v = readHeader(ds, r_exbase);

    if (v != 0)
        throw version_error(v, "0", r_exbase, CODELOC);

    return ds;
}

/** Constructor */
ExBase::ExBase() : RefCountData()
{}

/** Copy constructor */
ExBase::ExBase(const ExBase &) : RefCountData()
{}

/** Destructor */
ExBase::~ExBase()
{}

/** Comparison operator */
bool ExBase::operator!=(const ExBase &other) const
{
    return not operator==(other);
}

/** Return the negative of this ExBase */
Expression ExBase::operator-() const
{
    return -(Expression(*this));
}

/** Return an expression that the differential of this ExBase
    with respect to 'symbol'. Note an exception may
    be thrown if this ExBase cannot be differentiated.

    \throw SireCAS::unavailable_differential
*/
Expression ExBase::differentiate(const Symbol &symbol) const
{
    throw SireCAS::unavailable_differential(QObject::tr(
        "The differential of \"%1\" with respect to \"%2\" is not available.")
            .arg(toString(), symbol.toString()), CODELOC);
}

/** Return the indefinite integral of this 'ExBase' with respect to
    symbol. This is not guaranteed to work(!) and will return an
    expression of the form Sum( integral(exbase) + integral_constant ).
    If it doesn't work then an exception will be throw.

    \throw SireCAS::unavailable_integral
*/
Expression ExBase::integrate(const Symbol &symbol) const
{
    throw SireCAS::unavailable_integral(QObject::tr(
        "The integral of \"%1\" with respect to \"%2\" is not available.")
            .arg(toString(), symbol.toString()), CODELOC);
}

/** Return a series expansion of this expression with respect to
    'symbol', if possible, to order
    'n'. This is not guaranteed to work, and will return this expression
    unchanged if it doesn't work. If it is expanded, then a series
    will be returned, together with an estimate of the error (e.g. O(x^2)) */
Expression ExBase::series(const Symbol&, int) const
{
    return Expression(*this);
}

/** Try to simplify this expression. This will try to use known mathematical
    identities to convert complex expressions down to more simple ones.
    If SireCAS::UNSAFE_COMPLEX_SIMPLIFICATIONS is true, then identities
    that are not safe for complex math are used, e.g. z = sin(arcsin(z)). */
Expression ExBase::simplify(int) const
{
    return Expression(*this);
}

/** Return the complex conjugate of this expression */
Expression ExBase::conjugate() const
{
    return Expression(*this);
}

/** Return whether or not this is a function of the passed Symbol */
bool ExBase::isFunction(const Symbol&) const
{
    return false;
}

/** Return whether or not this is a constant expression (does not
    depend on any symbols) */
bool ExBase::isConstant() const
{
    return true;
}

/** Return whether or not this expression contains any complex (imaginary)
    parts */
bool ExBase::isComplex() const
{
    return false;
}

/** Return whether or not this is a compound expression, and thus as such
    requires brackets placed around it when it is printed. Examples include
    Sum, Product and Power. For most other functions it is safe to leave
    this as false. */
bool ExBase::isCompound() const
{
    return false;
}


QString ExBase::toOpenMMString() const{

    return this->toString();
}




namespace SireCAS
{

Expression SIRECAS_EXPORT operator+(const ExBase &base0, const ExBase &base1)
{
    return Expression(base0) + Expression(base1);
}

Expression SIRECAS_EXPORT operator+(const Expression &ex, const ExBase &base)
{
    return ex.add(base);
}

Expression SIRECAS_EXPORT operator+(const ExBase &base, const Expression &ex)
{
    return ex.add(base);
}

Expression SIRECAS_EXPORT operator+(double val, const ExBase &base)
{
    return Expression(base).add(val);
}

Expression SIRECAS_EXPORT operator+(const ExBase &base, double val)
{
    return Expression(base).add(val);
}

Expression SIRECAS_EXPORT operator+(const Complex &val, const ExBase &base)
{
    return Expression(base).add(val);
}

Expression SIRECAS_EXPORT operator+(const ExBase &base, const Complex &val)
{
    return Expression(base).add(val);
}

Expression SIRECAS_EXPORT operator-(const ExBase &base0, const ExBase &base1)
{
    return Expression(base0) - Expression(base1);
}

Expression SIRECAS_EXPORT operator-(const Expression &ex, const ExBase &base)
{
    return ex.subtract(base);
}

Expression SIRECAS_EXPORT operator-(const ExBase &base, const Expression &ex)
{
    return Expression(base).subtract(ex);
}

Expression SIRECAS_EXPORT operator-(double val, const ExBase &base)
{
    return Expression(val).subtract(base);
}

Expression SIRECAS_EXPORT operator-(const ExBase &base, double val)
{
    return Expression(base).subtract(val);
}

Expression SIRECAS_EXPORT operator-(const Complex &val, const ExBase &base)
{
    return Expression(val).subtract(base);
}

Expression SIRECAS_EXPORT operator-(const ExBase &base, const Complex &val)
{
    return Expression(base).subtract(val);
}

Expression SIRECAS_EXPORT operator*(const ExBase &base0, const ExBase &base1)
{
    return Expression(base0) * Expression(base1);
}

Expression SIRECAS_EXPORT operator*(const Expression &ex, const ExBase &base)
{
    return ex.multiply(base);
}

Expression SIRECAS_EXPORT operator*(const ExBase &base, const Expression &ex)
{
    return ex.multiply(base);
}

Expression SIRECAS_EXPORT operator*(double val, const ExBase &base)
{
    return Expression(base).multiply(val);
}

Expression SIRECAS_EXPORT operator*(const ExBase &base, double val)
{
    return Expression(base).multiply(val);
}

Expression SIRECAS_EXPORT operator*(const Complex &val, const ExBase &base)
{
    return Expression(base).multiply(val);
}

Expression SIRECAS_EXPORT operator*(const ExBase &base, const Complex &val)
{
    return Expression(base).multiply(val);
}

Expression SIRECAS_EXPORT operator/(const ExBase &base0, const ExBase &base1)
{
    return Expression(base0) / Expression(base1);
}

Expression SIRECAS_EXPORT operator/(const Expression &ex, const ExBase &base)
{
    return ex.divide(base);
}

Expression SIRECAS_EXPORT operator/(const ExBase &base, const Expression &ex)
{
    return Expression(base).divide(ex);
}

Expression SIRECAS_EXPORT operator/(double val, const ExBase &base)
{
    return Expression(val).divide(base);
}

Expression SIRECAS_EXPORT operator/(const ExBase &base, double val)
{
    return Expression(base).divide(val);
}

Expression SIRECAS_EXPORT operator/(const Complex &val, const ExBase &base)
{
    return Expression(val).divide(base);
}

Expression SIRECAS_EXPORT operator/(const ExBase &base, const Complex &val)
{
    return Expression(base).divide(val);
}

Expression SIRECAS_EXPORT pow(const ExBase &base, int n)
{
    return Expression(base).pow(n);
}

Expression SIRECAS_EXPORT pow(const ExBase &base, const Rational &n)
{
    return Expression(base).pow(n);
}

Expression SIRECAS_EXPORT pow(const ExBase &base, double n)
{
    return Expression(base).pow(n);
}

Expression SIRECAS_EXPORT pow(const ExBase &base, const Complex &n)
{
    return Expression(base).pow(n);
}

Expression SIRECAS_EXPORT pow(const ExBase &base, const Expression &n)
{
    return Expression(base).pow(n);
}

Expression SIRECAS_EXPORT pow(const ExBase &base, const ExBase &n)
{
    return Expression(base).pow(Expression(n));
}

}
