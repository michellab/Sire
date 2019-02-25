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

#include "singlefunc.h"
#include "values.h"
#include "complexvalues.h"
#include "identities.h"

#include "SireCAS/errors.h"

#include "SireStream/datastream.h"

using namespace SireStream;
using namespace SireCAS;

//register the pure virtual base class
static const RegisterMetaType<SingleFunc> r_singlefunc(MAGIC_ONLY, "SireCAS::SingleFunc");

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const SingleFunc &func)
{
    writeHeader(ds, r_singlefunc, 1) << func.ex << static_cast<const ExBase&>(func);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds, SingleFunc &func)
{
    VersionID v = readHeader(ds, r_singlefunc);

    if (v == 1)
    {
        ds >> func.ex >> static_cast<ExBase&>(func);
    }
    else
        throw version_error(v, "1", r_singlefunc, CODELOC);

    return ds;
}

/** Null constructor */
SingleFunc::SingleFunc() : ExBase()
{}

/** Construct a function that operates on the expression 'expression' */
SingleFunc::SingleFunc(const Expression &expression) : ExBase(), ex(expression)
{}

/** Copy constructor */
SingleFunc::SingleFunc(const SingleFunc &other)
           : ExBase(), ex(other.ex)
{}

/** Destructor */
SingleFunc::~SingleFunc()
{}

/** Copy assignment */
SingleFunc& SingleFunc::operator=(const SingleFunc &other)
{
    ExBase::operator=(other);
    ex = other.ex;

    return *this;
}

/** Return the conjugate of this function */
Expression SingleFunc::conjugate() const
{
    return functionOf( ex.conjugate() );
}

/** Return if this is a function of 'symbol' */
bool SingleFunc::isFunction(const Symbol &symbol) const
{
    return ex.isFunction(symbol);
}

/** Return whether or not this is constant */
bool SingleFunc::isConstant() const
{
    return ex.isConstant();
}

/** Return whether this is complex */
bool SingleFunc::isComplex() const
{
    return ex.isComplex();
}

/** A function is not compound */
bool SingleFunc::isCompound() const
{
    return false;
}

/** Return a string representation of this function */
QString SingleFunc::toString() const
{
    return QString("%1(%2)").arg(stringRep(), ex.toString());
}


/** Return a string representation of this function in the OpenMM syntax*/
QString SingleFunc::toOpenMMString() const
{
    return QString("%1(%2)").arg(stringRep(), ex.toOpenMMString());
}


/** Substitute into this expression */
Expression SingleFunc::substitute(const Identities &identities) const
{
    return functionOf( ex.substitute(identities) );
}

/** Return the symbols used in this function */
Symbols SingleFunc::symbols() const
{
    return ex.symbols();
}

/** Return the functions used in this function */
Functions SingleFunc::functions() const
{
    return ex.functions();
}

/** Return the child expression of this function */
Expressions SingleFunc::children() const
{
    return Expressions(ex);
}

/** Return the differential of this function with respect to 'symbol' */
Expression SingleFunc::differentiate(const Symbol &symbol) const
{
    if (ex.isFunction(symbol))
        return ex.differentiate(symbol) * this->diff();
    else
        return 0;
}

/** Return the integral of this function with respect to 'symbol' */
Expression SingleFunc::integrate(const Symbol &symbol) const
{
    if (ex.isFunction(symbol))
        return this->integ() / ex.differentiate(symbol);
    else
        return *this * symbol;
}

/** Default 'differentiate' function */
Expression SingleFunc::diff() const
{
    throw SireCAS::unavailable_differential(QObject::tr(
        "The function \"%1\" does not have an available differential. "
        "(full function = \"%2\")").arg(stringRep()).arg(toString()), CODELOC);

    return 0;
}

/** Default 'integrate' function */
Expression SingleFunc::integ() const
{
    throw SireCAS::unavailable_integral(QObject::tr(
        "The function \"%1\" does not have an available integral. "
        "(full function = \"%2\")").arg(stringRep()).arg(toString()), CODELOC);

    return 0;
}

QList<Factor> SingleFunc::expand(const Symbol &symbol) const
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
