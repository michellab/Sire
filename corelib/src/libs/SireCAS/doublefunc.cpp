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

#include "doublefunc.h"
#include "values.h"
#include "complexvalues.h"
#include "identities.h"

#include "SireCAS/errors.h"

#include "SireStream/datastream.h"

#include <QDebug>

using namespace SireStream;
using namespace SireCAS;

//register the pure virtual base class
static const RegisterMetaType<DoubleFunc> r_doublefunc(MAGIC_ONLY, "SireCAS::DoubleFunc");

/** Serialise to a binary datastream */
QDataStream SIRECAS_EXPORT &operator<<(QDataStream &ds, const DoubleFunc &func)
{
    writeHeader(ds, r_doublefunc, 1) 
            << func.ex0 << func.ex1 << static_cast<const ExBase&>(func);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream SIRECAS_EXPORT &operator>>(QDataStream &ds, DoubleFunc &func)
{
    VersionID v = readHeader(ds, r_doublefunc);

    if (v == 1)
    {
        ds >> func.ex0 >> func.ex1 >> static_cast<ExBase&>(func);
    }
    else
        throw version_error(v, "1", r_doublefunc, CODELOC);

    return ds;
}

/** Null constructor */
DoubleFunc::DoubleFunc() : ExBase()
{}

/** Construct a function that operates on the expressions 'x' and 'y' */
DoubleFunc::DoubleFunc(const Expression &x, const Expression &y) 
           : ExBase(), ex0(x), ex1(y)
{}

/** Copy constructor */
DoubleFunc::DoubleFunc(const DoubleFunc &other)
           : ExBase(), ex0(other.ex0), ex1(other.ex1)
{}

/** Destructor */
DoubleFunc::~DoubleFunc()
{}

/** Copy assignment */
DoubleFunc& DoubleFunc::operator=(const DoubleFunc &other)
{
    ExBase::operator=(other);
    ex0 = other.ex0;
    ex1 = other.ex1;

    return *this;
}

/** Return the conjugate of this function */
Expression DoubleFunc::conjugate() const
{
    return functionOf( ex0.conjugate(), ex1.conjugate() );
}

/** Return if this is a function of 'symbol' */
bool DoubleFunc::isFunction(const Symbol &symbol) const
{
    return ex0.isFunction(symbol) or ex1.isFunction(symbol);
}

/** Return whether or not this is constant */
bool DoubleFunc::isConstant() const
{
    return ex0.isConstant() and ex1.isConstant();
}

/** Return whether this is complex */
bool DoubleFunc::isComplex() const
{
    return ex0.isComplex() or ex1.isComplex();
}

/** A function is not compound */
bool DoubleFunc::isCompound() const
{
    return false;
}

/** Return a string representation of this function */
QString DoubleFunc::toString() const
{
    return QString("%1(%2,%3)").arg(stringRep(), ex0.toString(), ex1.toString());
}

/** Substitute into this expression */
Expression DoubleFunc::substitute(const Identities &identities) const
{
    return functionOf( ex0.substitute(identities),
                       ex1.substitute(identities) );
}

/** Return the symbols used in this function */
Symbols DoubleFunc::symbols() const
{
    return ex0.symbols() + ex1.symbols();
}

/** Return the functions used in this function */
Functions DoubleFunc::functions() const
{
    return ex0.functions() + ex1.functions();
}

/** Return the child expression of this function */
Expressions DoubleFunc::children() const
{
    return Expressions(ex0) + Expressions(ex1);
}

/** Return the differential of this function with respect to 'symbol' */
Expression DoubleFunc::differentiate(const Symbol &symbol) const
{
    if (ex0.isFunction(symbol) or ex1.isFunction(symbol))
        throw SireCAS::unavailable_differential( QObject::tr(
            "SireCAS does not know how to differentiate functions derived "
            "from \"DoubleFunc\" - so it can't differentiate %1 with respect "
            "to %2.")
                .arg(this->toString(), symbol.toString()), CODELOC );
                
    return Expression(0);
}

/** Return the integral of this function with respect to 'symbol' */
Expression DoubleFunc::integrate(const Symbol &symbol) const
{
    if (ex0.isFunction(symbol) or ex1.isFunction(symbol))
        throw SireCAS::unavailable_integral( QObject::tr(
            "SireCAS does not know how to integrate functions derived "
            "from \"DoubleFunc\" - so it can't integrate %1 with respect "
            "to %2.")
                .arg(this->toString(), symbol.toString()), CODELOC );

    return *this * symbol;
}

QList<Factor> DoubleFunc::expand(const Symbol &symbol) const
{
    if (this->isFunction(symbol))
    {
        //we cannot expand a function of this symbol...
        throw SireCAS::rearrangement_error( QObject::tr(
            "You cannot expand the function %1 in terms of the symbol %2.")
                .arg(this->toString(), symbol.toString()), CODELOC );
    }

    QList<Factor> ret;
    ret.append( Factor(symbol, *this, 0) );
    
    return ret;
}

/** Return a has for the function */
uint DoubleFunc::hash() const
{
    qDebug() << CODELOC;
    return (magic() << 16) | (ex0.hash() & 0x0000FF00) | (ex1.hash() & 0x000000FF);
}
