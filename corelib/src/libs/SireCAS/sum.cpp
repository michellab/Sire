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

#include <QMap>

#include "sum.h"
#include "expression.h"
#include "symbols.h"
#include "functions.h"
#include "complexvalues.h"
#include "values.h"
#include "identities.h"
#include "i.h"
#include "integrationconstant.h"

#include "SireStream/datastream.h"

#include <QDebug>

using namespace SireStream;
using namespace SireCAS;

static const RegisterMetaType<Sum> r_sum;

/** Return a hash for the Sum */
uint Sum::hash() const
{
    return ( r_sum.magicID() << 16 ) | ( (posparts.count() << 8) & 0x0000FF00 )
                                     | ( negparts.count() & 0x000000FF );
}

/** Serialise a Sum to a binary datastream */
QDataStream SIRECAS_EXPORT &operator<<(QDataStream &ds, const Sum &sum)
{
    writeHeader(ds, r_sum, 1) << sum.posparts.values()
                              << sum.negparts.values() << sum.strtval
                              << static_cast<const ExBase&>(sum);
    return ds;
}

/** Deserialise a Sum from a binary datastream */
QDataStream SIRECAS_EXPORT &operator>>(QDataStream &ds, Sum &sum)
{
    VersionID v = readHeader(ds, r_sum);

    if (v == 1)
    {
        QList<Expression> posparts, negparts;
        
        sum.posparts.clear();
        sum.negparts.clear();
        
        ds >> posparts >> negparts >> sum.strtval
           >> static_cast<ExBase&>(sum);

        foreach (Expression ex, posparts)
        {
            sum.posparts.insert( ex.base(), ex );
        }

        foreach (Expression ex, negparts)
        {
            sum.negparts.insert( ex.base(), ex );
        }
    }
    else
        throw version_error(v, "1", r_sum, CODELOC);

    return ds;
}

/** Construct an empty (zero) sum */
Sum::Sum() : ExBase(), strtval(0)
{}

/** Construct the sum of two expressions */
Sum::Sum(const Expression &ex0, const Expression &ex1) : ExBase(), strtval(0)
{
    this->add(ex0);
    this->add(ex1);
}

/** Construct the sum of the expressions in 'expressions' */
Sum::Sum(const Expressions &expressions) : ExBase(), strtval(0)
{
    int n = expressions.count();
    for (int i=0; i<n; ++i)
        this->add(expressions.at(i));
}

/** Copy constructor */
Sum::Sum(const Sum &other)
    : ExBase(), posparts(other.posparts), negparts(other.negparts),
      strtval(other.strtval)
{}

/** Destructor */
Sum::~Sum()
{}

/** Comparison operator */
bool Sum::operator==(const ExBase &other) const
{
    const Sum *othersum = dynamic_cast<const Sum*>(&other);

    return othersum != 0 and typeid(*this).name() == typeid(other).name() and
               othersum->posparts == posparts and othersum->negparts == negparts and
                   strtval == othersum->strtval;
}

/** Return a string representation of the sum */
QString Sum::toString() const
{
    if (posparts.count() == 0 and negparts.count() == 0)
        return QString::number(strtval);

    QString ret;

    int i = 0;

    for (QHash<ExpressionBase,Expression>::const_iterator it = posparts.begin();
         it != posparts.end();
         ++it)
    {
        if (i == 0)
            ret = it->toString();
        else
            ret = QString("%1 + %2").arg(ret,it->toString());

        ++i;
    }

    for (QHash<ExpressionBase,Expression>::const_iterator it = negparts.begin();
         it != negparts.end();
         ++it)
    {
        if (i == 0)
            ret = QString("-%1").arg(it->toString());
        else
            ret = QString("%1 - %2").arg(ret,it->toString());

        ++i;
    }

    if (not SireMaths::isZero(strtval))
    {
        if (i == 0)
            ret = QString::number(strtval);
        else if (strtval < 0)
            ret = QString("%1 - %2").arg(ret).arg(-strtval);
        else
            ret = QString("%1 + %2").arg(ret).arg(strtval);
    }

    return ret;
}

/** Return a string representation of the sum */
QString Sum::toOpenMMString() const
{
    if (posparts.count() == 0 and negparts.count() == 0)
        return QString::number(strtval);

    QString ret;

    int i = 0;

    for (QHash<ExpressionBase,Expression>::const_iterator it = posparts.begin();
         it != posparts.end();
         ++it)
    {
        if (i == 0)
            ret = it->toOpenMMString();
        else
            ret = QString("%1 + %2").arg(ret,it->toOpenMMString());

        ++i;
    }

    for (QHash<ExpressionBase,Expression>::const_iterator it = negparts.begin();
         it != negparts.end();
         ++it)
    {
        if (i == 0)
            ret = QString("-%1").arg(it->toOpenMMString());
        else
            ret = QString("%1 - %2").arg(ret,it->toOpenMMString());

        ++i;
    }

    if (not SireMaths::isZero(strtval))
    {
        if (i == 0)
            ret = QString::number(strtval);
        else if (strtval < 0)
            ret = QString("%1 - %2").arg(ret).arg(-strtval);
        else
            ret = QString("%1 + %2").arg(ret).arg(strtval);
    }

    return ret;
}


/** Remove the current version of 'ex', which will contain its current factor.
    Expression(0) will be returned if this expression is not in this Sum. */
Expression Sum::take(const ExpressionBase &ex)
{
    if (posparts.contains(ex))
        return posparts.take(ex);
    else if (negparts.contains(ex))
    {
        Expression negex = negparts.take(ex);
        return negex.negate();
    }
    else
        return Expression(0);
}

/** Add fac * ex */
void Sum::add(double fac, const ExpressionBase &ex)
{
    if (ex.isA<IntegrationConstant>())
    {
        //do not add these together - remove the current one...
        this->take(ex);

        //add 1 IntegrationConstant to the sum
        posparts.insert( ex, ex );
    }
    else
    {
        //get any expression that is currently in this Sum with this base
        Expression current_ex = this->take(ex);

        //get the current factor of this expression
        double current_factor = current_ex.factor();

        //calculate the new factor
        double new_factor = current_factor + fac;

        if (new_factor > 0)
            posparts.insert( ex, new_factor * ex );
        else if (new_factor < 0)
            negparts.insert( ex, -new_factor * ex );
    }
}

/** Add the expression 'ex' to this sum */
void Sum::add(const Expression &ex)
{
    //no need to add zero onto a sum!
    if (ex.isZero())
        return;

    else if (ex.base().isA<I>())
    {
        //this is a multiple of the complex number 'i'
        add(ex.factor(), ex.base());
    }
    else if (ex.isConstant())
    {
        Complex exval = ex.evaluate(ComplexValues());

        strtval += exval.real();

        if (not exval.isReal())
            this->add( exval.imag() * I() );
    }
    else if (ex.base().isA<Sum>())
    {
        //add the elements of the sum individually
        const Sum &sum = ex.base().asA<Sum>();

        if (posparts.count() == 0 and negparts.count() == 0)
        {
            if (ex.factor() == 1)
            {
                posparts = sum.posparts;
                negparts = sum.negparts;
                strtval += sum.strtval;
                
                return;
            }
            else if (ex.factor() == -1)
            {
                posparts = sum.negparts;
                negparts = sum.posparts;
                strtval -= sum.strtval;
                
                return;
            }
        }
            
        for (QHash<ExpressionBase,Expression>::const_iterator it = sum.posparts.begin();
             it != sum.posparts.end();
             ++it)
        {
           this->add( ex.factor()*it->factor(), it->base() );
        }

        for (QHash<ExpressionBase,Expression>::const_iterator it = sum.negparts.begin();
             it != sum.negparts.end();
             ++it)
        {
            this->add( -(ex.factor()*it->factor()), it->base() );
        }

        //add the start value to the sum
        strtval += ex.factor() * sum.strtval;
    }
    else
    {
        this->add(ex.factor(), ex.base());
    }
}

/** Reduce a Sum down to a simple form. This replaces the Sum with a single expression
    or a constant if this is no longer a Sum. It does not collapse together common
    factors - use 'collapse()' if you want to do this */
Expression Sum::reduce() const
{
    if (negparts.count() == 0 and posparts.count() == 0)
        return Expression(strtval);
    else if (SireMaths::isZero(strtval))
    {
        if (negparts.count() == 0 and posparts.count() == 1)
        {
            //sum contains only the positive part
            return posparts.values()[0];
        }
        else if (posparts.count() == 0 and negparts.count() == 1)
        {
            //sum contains only the negative part
            return negparts.values()[0].negate();
        }
    }

    return *this;
}

/** Simplify this sum */
Expression Sum::simplify(int options) const
{
    Sum ret;
    ret.strtval = strtval;

    //simplify the positive parts...
    for (QHash<ExpressionBase,Expression>::const_iterator it = posparts.begin();
         it != posparts.end();
         ++it)
    {
        ret.add( it->simplify(options) );
    }

    //now simplify the negative parts...
    for (QHash<ExpressionBase,Expression>::const_iterator it = negparts.begin();
         it != negparts.end();
         ++it)
    {
        ret.add( -(it->simplify(options)) );
    }

    return ret;
}

/** Return the conjugate of this sum */
Expression Sum::conjugate() const
{
    Sum ret;
    ret.strtval = strtval;

    //simplify the positive parts...
    for (QHash<ExpressionBase,Expression>::const_iterator it = posparts.begin();
         it != posparts.end();
         ++it)
    {
        ret.add( it->conjugate() );
    }

    //now simplify the negative parts...
    for (QHash<ExpressionBase,Expression>::const_iterator it = negparts.begin();
         it != negparts.end();
         ++it)
    {
        ret.add( -(it->conjugate()) );
    }

    return ret;
}

/** Evaluate the Sum for the values 'values'. Any missing values are assumed to be
    equal to zero. */
double Sum::evaluate(const Values &values) const
{
    double result = strtval;

    for ( QHash<ExpressionBase,Expression>::const_iterator it = posparts.begin();
          it != posparts.end();
          ++it )
    {
        result += it->evaluate(values);
    }

    for ( QHash<ExpressionBase,Expression>::const_iterator it = negparts.begin();
          it != negparts.end();
          ++it )
    {
        result -= it->evaluate(values);
    }

    return result;
}

/** Evaluate the Sum for the values 'values'. Any missing values are assumed to be
    equal to zero. */
Complex Sum::evaluate(const ComplexValues &values) const
{
    Complex result(strtval);

    for ( QHash<ExpressionBase,Expression>::const_iterator it = posparts.begin();
          it != posparts.end();
          ++it )
    {
        result += it->evaluate(values);
    }

    for ( QHash<ExpressionBase,Expression>::const_iterator it = negparts.begin();
          it != negparts.end();
          ++it )
    {
        result -= it->evaluate(values);
    }

    return result;
}

/** Return the differential of this Sum with respect to 'symbol'. */
Expression Sum::differentiate(const Symbol &symbol) const
{
    Sum diff;

    //add the differentials of all of the positive parts...
    for (QHash<ExpressionBase,Expression>::const_iterator it = posparts.begin();
         it != posparts.end();
         ++it)
    {
        diff.add( it->differentiate(symbol) );
    }

    //now add the differentials of all of the negative parts...
    for (QHash<ExpressionBase,Expression>::const_iterator it = negparts.begin();
         it != negparts.end();
         ++it)
    {
        diff.add( -(it->differentiate(symbol)) );
    }

    return diff.reduce();
}

/** Return the integral of this Sum with respect to 'symbol'. The integral of
    a sum is the sum of the integrals */
Expression Sum::integrate(const Symbol &symbol) const
{
    Sum integ;

    //add the integrals of all of the positive parts...
    for (QHash<ExpressionBase,Expression>::const_iterator it = posparts.begin();
         it != posparts.end();
         ++it)
    {
        integ.add( it->integrate(symbol) );
    }

    //now add the integrals of all of the negative parts...
    for (QHash<ExpressionBase,Expression>::const_iterator it = negparts.begin();
         it != negparts.end();
         ++it)
    {
        integ.add( -(it->integrate(symbol)) );
    }

    //add strtval * symbol
    if (not SireMaths::isZero(strtval))
        integ.add( strtval * symbol );

    return integ.reduce();
}

/** Return a series expansion of this sum about 'symbol' to order 'n' */
Expression Sum::series(const Symbol &symbol, int n) const
{
    Sum s;

    //add the expansions of all of the positive parts...
    for (QHash<ExpressionBase,Expression>::const_iterator it = posparts.begin();
         it != posparts.end();
         ++it)
    {
        s.add( it->series(symbol,n) );
    }

    //now add the expansions of all of the negative parts...
    for (QHash<ExpressionBase,Expression>::const_iterator it = negparts.begin();
         it != negparts.end();
         ++it)
    {
        s.add( -(it->series(symbol,n)) );
    }

    s.strtval = strtval;

    return s.reduce();
}

/** Return an expression that is this expression with 'identities' substituted in */
Expression Sum::substitute(const Identities &identities) const
{
    if (isConstant())
        return Expression(*this);
    else
    {
        Sum subsum;
        subsum.strtval = strtval;

        for (QHash<ExpressionBase,Expression>::const_iterator it = posparts.begin();
             it != posparts.end();
             ++it)
        {
            subsum.add( it->substitute(identities) );
        }

        for (QHash<ExpressionBase,Expression>::const_iterator it = negparts.begin();
             it != negparts.end();
             ++it)
        {
            subsum.add( it->negate().substitute(identities) );
        }

        return subsum.reduce();
    }
}

/** Return whether or not this is constant */
bool Sum::isConstant() const
{
    for (QHash<ExpressionBase,Expression>::const_iterator it = posparts.begin();
         it != posparts.end();
         ++it)
    {
        if (not it->isConstant())
            return false;
    }

    for (QHash<ExpressionBase,Expression>::const_iterator it = negparts.begin();
         it != negparts.end();
         ++it)
    {
        if (not it->isConstant())
            return false;
    }

    return true;
}

/** Return whether or not this is a function of 'symbol' */
bool Sum::isFunction(const Symbol &symbol) const
{
    for (QHash<ExpressionBase,Expression>::const_iterator it = posparts.begin();
         it != posparts.end();
         ++it)
    {
        if (it->isFunction(symbol))
            return true;
    }

    for (QHash<ExpressionBase,Expression>::const_iterator it = negparts.begin();
         it != negparts.end();
         ++it)
    {
        if (it->isFunction(symbol))
            return true;
    }

    return false;
}

/** Return whether or not this function contains any complex parts */
bool Sum::isComplex() const
{
    for (QHash<ExpressionBase,Expression>::const_iterator it = posparts.begin();
         it != posparts.end();
         ++it)
    {
        if (it->isComplex())
            return true;
    }

    for (QHash<ExpressionBase,Expression>::const_iterator it = negparts.begin();
         it != negparts.end();
         ++it)
    {
        if (it->isComplex())
            return true;
    }

    return false;
}

/** Return all of the symbols involved in this sum (and all expressions in this sum) */
Symbols Sum::symbols() const
{
    Symbols syms;

    for (QHash<ExpressionBase,Expression>::const_iterator it = posparts.begin();
         it != posparts.end();
         ++it)
    {
        syms.insert( it->symbols() );
    }

    for (QHash<ExpressionBase,Expression>::const_iterator it = negparts.begin();
         it != negparts.end();
         ++it)
    {
        syms.insert( it->symbols() );
    }

    return syms;
}

/** Return all of the functions involved in this sum (and all expressions in this sum) */
Functions Sum::functions() const
{
    Functions funcs;

    for (QHash<ExpressionBase,Expression>::const_iterator it = posparts.begin();
         it != posparts.end();
         ++it)
    {
        funcs.insert( it->functions() );
    }

    for (QHash<ExpressionBase,Expression>::const_iterator it = negparts.begin();
         it != negparts.end();
         ++it)
    {
        funcs.insert( it->functions() );
    }

    return funcs;
}

/** Return all of the child expressions in this Sum */
Expressions Sum::children() const
{
    Expressions exps;

    for (QHash<ExpressionBase,Expression>::const_iterator it = posparts.begin();
         it != posparts.end();
         ++it)
    {
        exps.append( *it );
    }

    for (QHash<ExpressionBase,Expression>::const_iterator it = negparts.begin();
         it != negparts.end();
         ++it)
    {
        exps.append( it->negate() );
    }

    return exps;
}

/** Return whether or not this is compound (needs brakets when printed) */
bool Sum::isCompound() const
{
    return posparts.count() + negparts.count() >= ( 2 - (not SireMaths::isZero(strtval)) );
}

QList<Factor> Sum::expand(const Symbol &symbol) const
{
    QHash<Expression, Expression> factors;
    
    for (QHash<ExpressionBase, Expression>::const_iterator it = posparts.constBegin();
         it != posparts.constEnd();
         ++it)
    {
        QList<Factor> facs = it->expand(symbol);

        foreach (const Factor &fac, facs)
        {
            factors[fac.power()] += fac.factor();
        }
    }

    for (QHash<ExpressionBase, Expression>::const_iterator it = negparts.constBegin();
         it != negparts.constEnd();
         ++it)
    {
        QList<Factor> facs = it->expand(symbol);

        foreach (const Factor &fac, facs)
        {
            factors[fac.power()] -= fac.factor();
        }
    }
    
    QList<Factor> ret;
    
    for (QHash<Expression,Expression>::const_iterator it = factors.constBegin();
         it != factors.constEnd();
         ++it)
    {
        ret.append( Factor( symbol, it.value(), it.key() ) );
    }
    
    if (strtval != 0)
    {
        ret.append( Factor(symbol, strtval, 0) );
    }
    
    return ret;
}

const char* Sum::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Sum>() );
}

Sum* Sum::clone() const
{
    return new Sum(*this);
}

