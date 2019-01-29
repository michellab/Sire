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

#include "product.h"
#include "sum.h"
#include "expression.h"
#include "expressions.h"
#include "symbol.h"
#include "symbols.h"
#include "functions.h"
#include "power.h"
#include "integrationconstant.h"
#include "complexvalues.h"
#include "identities.h"
#include "i.h"

#include "SireCAS/errors.h"

#include "SireStream/datastream.h"

#include <boost/assert.hpp>

#include <QDebug>

using namespace SireStream;
using namespace SireCAS;

static const RegisterMetaType<Product> r_product;

/** Return a hash for this product */
uint Product::hash() const
{
    return ( r_product.magicID() << 16) | ( (numparts.count()<<8) & 0x0000FF00 )
                                        | ( denomparts.count() & 0x000000FF );
}

/** Serialise a Product to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const Product &product)
{
    writeHeader(ds, r_product, 1)
          << product.powers << product.strtval << static_cast<const ExBase&>(product);

    return ds;
}

/** Deserialise a Product from a binary datastream */
QDataStream &operator>>(QDataStream &ds, Product &product)
{
    VersionID v = readHeader(ds, r_product);

    if (v == 1)
    {
        ds >> product.powers >> product.strtval
           >> static_cast<ExBase&>(product);

        //rebuild numparts and denomparts from powers
        product.rebuild();
    }
    else
        throw version_error(v, "1", r_product, CODELOC);

    return ds;
}

/** Construct an empty, zero Product */
Product::Product() : ExBase(), strtval(0)
{}

/** Construct the product of two expressions */
Product::Product(const Expression &ex0, const Expression &ex1)
        : ExBase(), strtval(1)
{
    multiply(ex0);
    multiply(ex1);
}

/** Construct the product of all of the expressions in 'expressions' */
Product::Product(const Expressions &expressions)
        : ExBase(), strtval(1)
{
    int n = expressions.count();

    for (int i=0; i<n; ++i)
        multiply(expressions.at(i));
}

/** Copy constructor */
Product::Product(const Product &other)
        : ExBase(), powers(other.powers), numparts(other.numparts),
          denomparts(other.denomparts), strtval(other.strtval)
{}

/** Destructor */
Product::~Product()
{}

/** Comparison operator */
bool Product::operator==(const ExBase &other) const
{
    const Product *otherprod = dynamic_cast<const Product*>(&other);

    return otherprod != 0 and typeid(*this).name() == typeid(other).name() and
               strtval == otherprod->strtval and
                  powers == otherprod->powers;

}

/** Evaluate this product */
double Product::evaluate(const Values &values) const
{
    if (SireMaths::isZero(strtval))
        return 0;
    else
    {
        //evaluate the numerator
        double numerator = strtval;

        for ( QHash<Expression,Expression>::const_iterator it = numparts.begin();
              it != numparts.end();
              ++it)
        {
            numerator *= it->evaluate(values);
        }

        if (denomparts.count() == 0 or SireMaths::isZero(numerator))
            return numerator;

        //evaluate the denominator
        double denominator = 1.0;
        for ( QHash<Expression,Expression>::const_iterator it = denomparts.begin();
              it != denomparts.end();
              ++it)
        {
            denominator *= it->evaluate(values);
        }

        //return the ratio (eventually we could try and be more clever about this...)
        return numerator / denominator;
    }
}

/** Evaluate this product */
Complex Product::evaluate(const ComplexValues &values) const
{
    if (SireMaths::isZero(strtval))
        return Complex(0);
    else if (numparts.count() == 0 and denomparts.count() == 0)
        return Complex(strtval);
    else
    {
        //evaluate the numerator
        Complex numerator(strtval);

        for ( QHash<Expression,Expression>::const_iterator it = numparts.begin();
              it != numparts.end();
              ++it)
        {
            numerator *= it->evaluate(values);
        }

        if (numerator.isZero() or denomparts.count() == 0)
            return numerator;

        //evaluate the denominator
        Complex denominator(1);
        for ( QHash<Expression,Expression>::const_iterator it = denomparts.begin();
              it != denomparts.end();
              ++it)
        {
            denominator *= it->evaluate(values);
        }

        //return the ratio (eventually we could try and be more clever about this...)
        return numerator / denominator;
    }
}

/** Need to write integral of a product... */
Expression Product::integrate(const Symbol &symbol) const
{
    if (this->isFunction(symbol))
        throw SireCAS::unavailable_integral(QObject::tr(
            "Integral of a product has yet to be coded. Cannot integrate \"%1\" with "
            "respect to %2").arg(toString(),symbol.toString()), CODELOC);

    return *this * symbol;
}

/** Return the product with the identities in 'identities' substituted in */
Expression Product::substitute(const Identities &identities) const
{
    if (this->isConstant())
        return Expression(*this);
    else
    {
        Product subproduct;
        subproduct.strtval = strtval;

        for (QHash<Expression,Expression>::const_iterator it = numparts.begin();
             it != numparts.end();
             ++it)
        {
            subproduct.multiply( it->substitute(identities) );
        }

        for (QHash<Expression,Expression>::const_iterator it = denomparts.begin();
             it != denomparts.end();
             ++it)
        {
            subproduct.multiply( it->invert().substitute(identities) );
        }

        return subproduct.reduce();
    }
}

/** Reduce a Product down to a simple form. This will not attempt to collapse
    common factors - if you want to do this then call the 'collapse' function. */
Expression Product::reduce() const
{
    //return nothing if this is zero
    if (isConstant())
        return Expression(evaluate(ComplexValues()));
    else if (numparts.count() == 1 and denomparts.count() == 0)
    {
        //product only contains a single numerator
        return numparts.values()[0].multiply(strtval);
    }
    else if (denomparts.count() == 1 and numparts.count() == 0)
    {
        //product only contains the denominator
        return denomparts.values()[0].invert().multiply(strtval);
    }

    //move the 'strtval' into the expression
    Product p(*this);
    p.strtval = 1;

    return p * strtval;
}

/** Internal function: Multiply this Product by a constant */
void Product::multiply(double val)
{
    strtval *= val;

    if ( SireMaths::isZero(strtval) )
    {
        numparts.clear();
        denomparts.clear();
        powers.clear();
    }
}

/** Remove the current value of core^power, and return 'power', or return
    '0' if there is no expression of the form core^power in this Product */
Expression Product::take(const Expression &core)
{
    //remove this core from numparts and denomparts (if it exists)
    numparts.remove(core);
    denomparts.remove(core);

    //remove and return the current expression from powers - this returns
    //Expression(0) if an expression of the form 'core^power' is
    //not in this Product
    return powers.take(core);
}

/** Internal function: Multiply by ex^power. No expansion or inspection of 'ex'
    is performed by this function */
void Product::multiplyPvt(const Expression &ex, const Expression &power)
{
    //get the old power
    Expression oldpower = this->take(ex);

    Expression newpower;

    if (oldpower.isZero())
        newpower = power;
    else
        //calculate the new power
        newpower = oldpower + power;

    if (newpower.isConstant())
    {
        if (ex.base().isA<IntegrationConstant>())
        {
            //we only want C^1
            powers.insert( ex.base(), Expression(1) );
            numparts.insert( ex.base(), ex.base() );
        }
        else
        {
            //calculate the power
            Complex powerval = newpower.evaluate(ComplexValues());

            //is this zero?
            if (powerval.isZero())
            {
                //x^0 == 1, so leaves the product
                return;
            }
            //is this i?
            else if (ex.base().isA<I>())
            {
                Complex z = SireMaths::pow( Complex(0,1), powerval );

                if (z.isReal())
                    multiply(z.real());
                else if (SireMaths::isZero(z.real()))
                {
                    multiply(z.imag());
                    powers.insert( ex.base(), Expression(1) );
                    numparts.insert( ex.base(), ex.base() );
                }
                else
                {
                    Expression cex = z.real() + z.imag() * I();
                    powers.insert( cex.base(), Expression(1) );
                    numparts.insert( cex.base(), Expression(1) );
                }
            }
            else if (powerval.isReal())
            {
                double realpower = powerval.real();

                powers.insert( ex, Expression(realpower) );

                if (realpower > 0)
                    numparts.insert( ex, ex.pow(realpower) );
                else
                    denomparts.insert( ex, ex.pow(-realpower) );
            }
            else
            {
                powers.insert( ex, Expression(powerval) );
                numparts.insert( ex, ex.pow(powerval) );
            }
        }
    }
    else
    {
        //just put the power on the top
        numparts.insert( ex, ex.pow(newpower) );

        //save the expression and its power
        powers.insert(ex, newpower);
    }
}

/** Internal function: Multiply this Product by a Complex */
void Product::multiply(const Complex &complex)
{
    if (complex.isReal())
        multiply(complex.real());
    else if ( SireMaths::isZero(complex.real()) )
        multiplyPvt( complex.imag() * I(), Expression(1) );
    else
        multiplyPvt( complex.real() + complex.imag()*I(), Expression(1) );
}

/** Internal function: Multiply by a number raised to a power */
void Product::multiply(double val, const Expression &power)
{
    if (SireMaths::areEqual(val,1.0))
        //nothing needs doing as 1^n == 1, even if n is complex
        return;
    else if (power.isConstant())
    {
        Complex powerval = power.evaluate(ComplexValues());
        multiply( SireMaths::pow(val, powerval) );
    }
    else
        multiplyPvt(Expression(val), power);
}

/** Internal function: Multiply by a complex raised to a power */
void Product::multiply(const Complex &complex, const Expression &power)
{
    if (complex.isReal())
        multiply(complex.real(), power);
    else if (power.isConstant())
    {
        multiply( SireMaths::pow(complex, power.evaluate(ComplexValues())) );
    }
    else
        multiplyPvt( complex.real() + complex.imag()*I(), power );
}

/** Internal function: Multiply this Product by a Product */
void Product::multiply(const Product &product, const Expression &power)
{
    //first, the strtval of the product
    multiply( product.strtval, power );

    //now the product's numerator...
    for (QHash<Expression,Expression>::const_iterator it = product.numparts.begin();
         it != product.numparts.end();
         ++it)
    {
        multiply( it->factor(), power );
        multiply( it->base(), power );
    }

    //finally the denominator
    for (QHash<Expression,Expression>::const_iterator it = product.denomparts.begin();
         it != product.denomparts.end();
         ++it)
    {
        multiply( it->factor(), -power );
        multiply( it->base(), -power );
    }
}

/** Internal function: Multiply by base^power */
void Product::multiply(const ExpressionBase &base, const Expression &power)
{
    if (base.isConstant())
    {
        Complex baseval = base.evaluate(ComplexValues());
        this->multiply(baseval, power);
    }
    else if (base.isA<PowerFunction>())
    {
        const PowerFunction &powerfunc = base.asA<PowerFunction>();

        Expression core = powerfunc.core();
        Expression combined_power = power * powerfunc.power();

        multiply( core.factor(), combined_power );
        multiply( core.base(), combined_power );
    }
    else if (base.isA<Product>())
    {
        multiply( base.asA<Product>(), power );
    }
    else
    {
        multiplyPvt( Expression(base), power );
    }
}

/** Multiply this product by the expression 'ex' */
void Product::multiply(const Expression &ex)
{
    if (SireMaths::isZero(strtval))
        //no change if we are already equal to zero!
        return;
    else if (ex.base().isA<Product>())
    {
        //multiply by the factor on the product
        multiply(ex.factor());

        const Product &product = ex.base().asA<Product>();

        if (numparts.count() == 0 and denomparts.count() == 0)
        {
            numparts = product.numparts;
            denomparts = product.denomparts;
            powers = product.powers;
            multiply(product.strtval);
        }
        else
            multiply(product, Expression(1));
    }
    else if (ex.isConstant())
    {
        Complex constval = ex.evaluate(ComplexValues());
        this->multiply( constval, Expression(1) );
    }
    else
    {
        this->multiply( ex.factor() );
        this->multiply( ex.base(), Expression(1) );
    }
}

/** Rebuild numparts and denomparts from the values in powers */
void Product::rebuild()
{
    if (SireMaths::isZero(strtval))
    {
        powers.clear();
        numparts.clear();
        denomparts.clear();
        return;
    }

    //clear numparts and denomparts...
    numparts.clear();
    denomparts.clear();

    //rebuild them from 'powers'
    for (QHash<Expression,Expression>::const_iterator it = powers.begin();
         it != powers.end();
         ++it)
    {
        Expression core = it.key();
        Expression power = it.value();

        if (power.isZero())
            continue;
        else if (power.isConstant())
        {
            Complex powerval = power.evaluate(ComplexValues());

            if (powerval.isReal())
            {
                if (powerval.real() > 0)
                    numparts.insert( core, core.pow(powerval.real()) );
                else
                    denomparts.insert( core, core.pow(-powerval.real()) );
            }
            else
                numparts.insert( core, core.pow(powerval) );
        }
        else
            numparts.insert( core, core.pow(power) );
    }
}

/** Return a string representation of this Product */
QString Product::toString() const
{
    if (isConstant())
        return evaluate(ComplexValues()).toString();

    QString top("");

    bool strtval_is_one = SireMaths::areEqual(strtval,1.0);
    bool strtval_is_minus_one = SireMaths::areEqual(strtval,-1.0);

    if (numparts.count() == 0)
    {
        top = QString::number(strtval);
    }
    else if (numparts.count() == 1)
    {
        Expression toppart = numparts.values()[0];

        if ( not toppart.isCompound() or
             (denomparts.count() == 0 and (strtval_is_one or strtval_is_minus_one)) )
            top = toppart.toString();
        else
            top = QString("[%1]").arg(toppart.toString());
    }
    else
    {
        for (QHash<Expression,Expression>::const_iterator it = numparts.begin();
            it != numparts.end();
            ++it)
        {
            if (it->base().isCompound())
                top = QString("%1 [%2]").arg(top,it->toString());
            else
                top = QString("%1 %2").arg(top,it->toString());
        }
    }

    if (strtval_is_minus_one)
        top = QString("-%1").arg(top);
    else if ( not (strtval_is_one or numparts.count() == 0) )
        top = QString("%1 %2").arg(strtval).arg(top);

    if (denomparts.count() == 0)
        return top.simplified();
    else if (denomparts.count() == 1)
    {
        Expression bottom = denomparts.values()[0];
        if (bottom.base().isCompound())
            return QString("%1 / [%2]").arg(top.simplified(),bottom.toString());
        else
            return QString("%1 / %2").arg(top.simplified(),bottom.toString());
    }
    else
    {
        QString bottom("");

        for (QHash<Expression,Expression>::const_iterator it = denomparts.begin();
             it != denomparts.end();
             ++it)
        {
            if (it->base().isCompound())
                bottom = QString("%1 [%2]").arg(bottom,it->toString());
            else
                bottom = QString("%1 %2").arg(bottom,it->toString());
        }

        return QString("%1 / [%2]").arg(top.simplified(),bottom.simplified());
    }
}


/** Return a string representation of this Product in the OpenMM syntax*/
QString Product::toOpenMMString() const
{
    if (isConstant())
        return evaluate(ComplexValues()).toString();

    QString top("");

    bool strtval_is_one = SireMaths::areEqual(strtval,1.0);
    bool strtval_is_minus_one = SireMaths::areEqual(strtval,-1.0);

    if (numparts.count() == 0)
    {
        top = QString::number(strtval);
    }
    else if (numparts.count() == 1)
    {
        Expression toppart = numparts.values()[0];

        if ( not toppart.isCompound() or
             (denomparts.count() == 0 and (strtval_is_one or strtval_is_minus_one)) )
            top = toppart.toOpenMMString();
        else
            top = QString("(%1)").arg(toppart.toOpenMMString());
    }
    else
    {
        for (QHash<Expression,Expression>::const_iterator it = numparts.begin();
            it != numparts.end();
            ++it)
        {
            if(top.isEmpty()){
                 if (it->base().isCompound())
                    top = QString("(%1)").arg(it->toOpenMMString());
                else
                    top = it->toOpenMMString();
            }
            else{
                if (it->base().isCompound())
                    top = QString("%1 * (%2)").arg(top,it->toOpenMMString());
                else
                    top = QString("%1 * %2").arg(top,it->toOpenMMString());
            }
        }
    }

    if (strtval_is_minus_one)
        top = QString("-%1").arg(top);
    else if ( not (strtval_is_one or numparts.count() == 0) )
        top = QString("%1 * %2").arg(strtval).arg(top);

    if (denomparts.count() == 0)
        return top.simplified();
    else if (denomparts.count() == 1)
    {
        Expression bottom = denomparts.values()[0];
        if (bottom.base().isCompound())
            return QString("%1 / (%2)").arg(top.simplified(),bottom.toOpenMMString());
        else
            return QString("%1 / %2").arg(top.simplified(),bottom.toOpenMMString());
    }
    else
    {
        QString bottom("");

        for (QHash<Expression,Expression>::const_iterator it = denomparts.begin();
             it != denomparts.end();
             ++it)
        {
            if (it->base().isCompound())
                bottom = QString("%1 * (%2)").arg(bottom,it->toOpenMMString());
            else
                bottom = QString("%1 * %2").arg(bottom,it->toOpenMMString());
        }

        return QString("%1 / (%2)").arg(top.simplified(),bottom.simplified());
    }
}


/** Return whether this is a constant */
bool Product::isConstant() const
{
    for (QHash<Expression,Expression>::const_iterator it = numparts.begin();
         it != numparts.end();
         ++it)
    {
        if (not it->isConstant())
            return false;
    }

    for (QHash<Expression,Expression>::const_iterator it = denomparts.begin();
         it != denomparts.end();
         ++it)
    {
        if (not it->isConstant())
            return false;
    }

    return true;
}

/** Return whether or not this Product contains any complex terms */
bool Product::isComplex() const
{
    if (isConstant())
    {
        Complex cval = evaluate(ComplexValues());
        return not cval.isReal();
    }
    else
    {
        for (QHash<Expression,Expression>::const_iterator it = numparts.begin();
             it != numparts.end();
             ++it)
        {
            if (it->isComplex())
                return true;
        }

        for (QHash<Expression,Expression>::const_iterator it = denomparts.begin();
             it != denomparts.end();
             ++it)
        {
            if (it->isComplex())
                return true;
        }

        return false;
    }
}

/** Return whether or not this is a function of 'symbol' */
bool Product::isFunction(const Symbol &symbol) const
{
    if (isConstant())
        return false;
    else
    {
        for (QHash<Expression,Expression>::const_iterator it = numparts.begin();
             it != numparts.end();
             ++it)
        {
            if (it->isFunction(symbol))
                return true;
        }

        for (QHash<Expression,Expression>::const_iterator it = denomparts.begin();
             it != denomparts.end();
             ++it)
        {
            if (it->isFunction(symbol))
                return true;
        }

        return false;
    }
}

/** Return whether or not this is compound (requires brackets when printing) */
bool Product::isCompound() const
{
    return (not isConstant())
        and (denomparts.count() > 0 or numparts.count() > 1
                or not SireMaths::areEqual(strtval,1.0));
}

/** Return all of the symbols used in this product */
Symbols Product::symbols() const
{
    if (isConstant())
        return Symbols();
    else
    {
        Symbols syms;

        for (QHash<Expression,Expression>::const_iterator it = numparts.begin();
            it != numparts.end();
            ++it)
        {
            syms.insert( it->symbols() );
        }

        for (QHash<Expression,Expression>::const_iterator it = denomparts.begin();
            it != denomparts.end();
            ++it)
        {
            syms.insert( it->symbols() );
        }

        return syms;
    }
}

/** Return all of the functions used in this product */
Functions Product::functions() const
{
    if (isConstant())
        return Functions();
    else
    {
        Functions funcs;

        for (QHash<Expression,Expression>::const_iterator it = numparts.begin();
            it != numparts.end();
            ++it)
        {
            funcs.insert( it->functions() );
        }

        for (QHash<Expression,Expression>::const_iterator it = denomparts.begin();
            it != denomparts.end();
            ++it)
        {
            funcs.insert( it->functions() );
        }

        return funcs;
    }
}

/** Return the child expressions of this product */
Expressions Product::children() const
{
    Expressions exps;

    for (QHash<Expression,Expression>::const_iterator it = numparts.begin();
         it != numparts.end();
         ++it)
    {
        exps.append( *it );
    }

    for (QHash<Expression,Expression>::const_iterator it = denomparts.begin();
         it != denomparts.end();
         ++it)
    {
        exps.append( it->invert() );
    }

    return exps;
}

/** Return whether or not this is a pure product (has no denominator) */
bool Product::isPureProduct() const
{
    return denomparts.count() == 0;
}

/** Return the Product of expressions on the numerator of this Product */
Product Product::numerator() const
{
    if (isPureProduct())
        return *this;
    else
        return Product( Expressions(numparts.values()) );
}

/** Return the Product of expressions on the denominator of this Product */
Product Product::denominator() const
{
    return Product( Expressions(denomparts.values()) );
}

/** Remove the first expression in the product that depends on 'symbol' and
    return that expression. If there are no expressions with this symbol,
    then Expression(1) is returned. Note that this scans the expressions
    in the numerator before it moves to the expressions in the denominator */
Expression Product::takeFirst(const Symbol &symbol)
{
    QHash<Expression,Expression>::iterator it = numparts.begin();

    while( it != numparts.end() )
    {
        if (it.value().isFunction(symbol))
        {
            Expression f_symbol = it.value();

            //remove this expression from 'powers'
            powers.remove(it.key());

            //remove this expression from numparts
            it = numparts.erase(it);

            //return this expression
            return f_symbol;
        }
        else
            ++it;
    }

    it = denomparts.begin();
    while( it != denomparts.end() )
    {
        if (it.value().isFunction(symbol))
        {
            Expression f_symbol = it.value();

            //remove this expression from 'powers'
            powers.remove(it.key());

            //remove this expression from denomparts
            it = denomparts.erase(it);

            //return this expression (inverted, as it is on the denominator)
            return f_symbol.invert();
        }
        else
            ++it;
    }

    //no expression has been found - return Expression(1)
    return Expression(1);
}

/** Internal function: Return the differential of 'product' using the product rule. Note that
    'product' must be a pure product */
Expression Product::productRule(Product product, const Symbol &symbol)
{
    BOOST_ASSERT(product.isPureProduct());

    //is this product a function of 'symbol'? If not then return 0
    if (not product.isFunction(symbol))
        return Expression(0);

    //ok - this could be a product of many individual functions...
    // break this down a function at a time...

    //remove the first function of symbol - f(x)
    Expression f = product.takeFirst(symbol);

    BOOST_ASSERT(f.isFunction(symbol));  //has to be true, as product.isFunction(symbol)

    if (product.isFunction(symbol))
    {
        //product = g(x)
        Expression g = product.reduce();

        //product rule is d f(x)g(x) / dx  = f'(x)g(x) + f(x)g'(x)
        Expression df = f.differentiate(symbol);
        Expression dg = g.differentiate(symbol);

        return df*g + f*dg;
    }
    else
    {
        //this is differential of constant * f(x) = constant * f'(x)
        return product * f.differentiate(symbol);
    }
}

/** Internal function: Return the differential of 'f(x)/g(x)'. Note
    that both of these products must be pure products */
Expression Product::quotientRule(const Product &f, const Product &g,
                                 const Symbol &symbol)
{
    BOOST_ASSERT(f.isPureProduct());
    BOOST_ASSERT(f.isPureProduct());

    //quotient rule for f(x)/g(x),  calculate d/dx
    //   = [g(x)f'(x) - f(x)g'(x)] / g(x)^2

    Expression g2 = pow(g, 2);
    BOOST_ASSERT( not g2.isZero() );

    Expression dg = productRule(g, symbol);
    Expression df = productRule(f, symbol);

    return ( (g*df) - (f*dg) ) / g2;
}

/** Return the differential of this product... */
Expression Product::differentiate(const Symbol &symbol) const
{
    if (isConstant() or not isFunction(symbol))
        return Expression(0);
    else if (isPureProduct())
    {
        //pure product, so can use product rule
        return productRule(*this, symbol);
    }
    else
    {
        //this is really a quotient. My need to use the quotient rule.
        //split into numerator and denominator
        Product num = this->numerator();
        Product denom = this->denominator();

        bool numerator_is_variable = num.isFunction(symbol);
        bool denominator_is_variable = denom.isFunction(symbol);

        if (numerator_is_variable and denominator_is_variable)
        {
            //use the quotient rule to get the differentials of these two parts
            return quotientRule(num, denom, symbol);
        }
        else if (numerator_is_variable)
        {
            //the result is numerator' / denominator (since denominator is a constant)
            return num.differentiate(symbol) / denom;
        }
        else
        {
            //the result is numerator * [ (denominator)^1 ]'
            //  = -numerator * denominator' * (denominator)^-2
            return -( num * denom.differentiate(symbol) * pow(denom,-2) );
        }
    }
}

/** Return the series expansion of this product with respect to 'symbol', to order 'n'*/
Expression Product::series(const Symbol &symbol, int n) const
{
    if (not isFunction(symbol))
        return *this;
    else
    {
        Product ret;
        ret.strtval = strtval;

        for (QHash<Expression,Expression>::const_iterator it = numparts.begin();
            it != numparts.end();
            ++it)
        {
            ret.multiply( it->series(symbol,n) );
        }

        for (QHash<Expression,Expression>::const_iterator it = denomparts.begin();
            it != denomparts.end();
            ++it)
        {
            ret.multiply( pow(it->series(symbol,n),-1) );
        }

        return ret.reduce();
    }
}

/** Try to simplify this product */
Expression Product::simplify(int options) const
{
    Expression ret(strtval);

    for (QHash<Expression,Expression>::const_iterator it = numparts.begin();
        it != numparts.end();
        ++it)
    {
        if (it->base().isA<Sum>())
        {
            const Sum &sum = it->base().asA<Sum>();
            
            Expression part = ret;
            ret = ((it->factor() * sum.strtval) 
                        * part.simplify(options))
                                .simplify(options);
            
            for (QHash<ExpressionBase, Expression>::const_iterator
                                                        it2 = sum.posparts.constBegin();
                 it2 != sum.posparts.constEnd();
                 ++it2)
            {
                ret += (part.simplify(options) * 
                            ( it->factor() * it2.value().simplify(options) ))
                                    .simplify(options);
            }

            for (QHash<ExpressionBase, Expression>::const_iterator
                                                        it2 = sum.negparts.constBegin();
                 it2 != sum.negparts.constEnd();
                 ++it2)
            {
                ret -= (part.simplify(options) * 
                            ( it->factor() * it2.value().simplify(options) ))
                                    .simplify(options);
            }
        }
        else if (ret.base().isA<Sum>())
        {
            Expression old_ret = ret;
            const Sum &sum = old_ret.base().asA<Sum>();
            
            Expression part = it.value();
            ret = ((old_ret.factor() * sum.strtval) 
                        * part.simplify(options))
                                .simplify(options);
            
            for (QHash<ExpressionBase, Expression>::const_iterator
                                                        it2 = sum.posparts.constBegin();
                 it2 != sum.posparts.constEnd();
                 ++it2)
            {
                ret += (part.simplify(options) * 
                            ( it->factor() * it2.value().simplify(options) ))
                                    .simplify(options);
            }

            for (QHash<ExpressionBase, Expression>::const_iterator
                                                        it2 = sum.negparts.constBegin();
                 it2 != sum.negparts.constEnd();
                 ++it2)
            {
                ret -= (part.simplify(options) * 
                            ( it->factor() * it2.value().simplify(options) ))
                                    .simplify(options);
            }
        }
        else
        {
            ret *= it->simplify(options);
        }
    }

    for (QHash<Expression,Expression>::const_iterator it = denomparts.begin();
        it != denomparts.end();
        ++it)
    {
        ret *= pow(it->simplify(options),-1);
    }

    if (ret.base().isA<Product>())
        return ret.factor() * ret.base().asA<Product>().reduce();
    else
        return ret;
}

/** Return the complex conjugate of this product */
Expression Product::conjugate() const
{
    Product ret;
    ret.strtval = strtval;

    for (QHash<Expression,Expression>::const_iterator it = numparts.begin();
        it != numparts.end();
        ++it)
    {
        ret.multiply( it->conjugate() );
    }

    for (QHash<Expression,Expression>::const_iterator it = denomparts.begin();
        it != denomparts.end();
        ++it)
    {
        ret.multiply( pow(it->conjugate(), -1) );
    }

    return ret.reduce();
}

static QList<Factor> multiply(const QList<Factor> &f0s, const QList<Factor> &f1s)
{
    if (f0s.isEmpty())
        return f1s;
    else if (f1s.isEmpty())
        return f0s;

    QHash<Expression,Expression> factors;
    
    foreach (const Factor &f0, f0s)
    {
        foreach (const Factor &f1, f1s)
        {
            factors[ f0.power() + f1.power() ] += (f0.factor() * f1.factor());
        }
    }
    
    QList<Factor> ret;
    
    for (QHash<Expression,Expression>::const_iterator it = factors.constBegin();
         it != factors.constEnd();
         ++it)
    {
        ret.append( Factor( f0s.at(0).symbol(), it.value(), it.key() ) );
    }
    
    return ret;
}

static QList<Factor> expand(const QHash<Expression,Expression> &product,
                            const Symbol &symbol)
{
    QList<Factor> ret;

    for (QHash<Expression,Expression>::const_iterator it = product.constBegin();
         it != product.constEnd();
         ++it)
    {
        QList<Factor> factors = it->expand(symbol);
        
        ret = multiply(ret, factors);
    }
    
    return ret;
}

QList<Factor> Product::expand(const Symbol &symbol) const
{
    if (strtval == 0)
        return QList<Factor>();

    QList<Factor> denom_factors = ::expand(denomparts, symbol);

    if (denom_factors.count() > 1)
    {
        //we cannot expand functions like ( mx^i + nx^j )^-1
        throw SireCAS::rearrangement_error( QObject::tr(
            "The product %1 cannot be expanded in terms of the symbol %2 "
            "as the denominator contains more than one power of %2.")
                .arg(this->toString(), symbol.toString()), CODELOC );
    }
    
    QList<Factor> num_factors = ::expand(numparts, symbol);
    
    //multiply the numerator by the constant...
    if (strtval != 1)
    {
        for (QList<Factor>::iterator it = num_factors.begin();
             it != num_factors.end();
             ++it)
        {
            *it = Factor( symbol, strtval * it->factor(), it->power() );
        }
    }
    
    //now invert the power of the denominator...
    for (QList<Factor>::iterator it = denom_factors.begin();
         it != denom_factors.end();
         ++it)
    {
        *it = Factor( it->symbol(), it->factor(), -(it->power()) );
    }
    
    return ::multiply(num_factors, denom_factors);
}

const char* Product::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Product>() );
}

Product* Product::clone() const
{
    return new Product(*this);
}

