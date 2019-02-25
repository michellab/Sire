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

#include "power.h"
#include "exp.h"
#include "powerconstant.h"
#include "values.h"
#include "complexvalues.h"
#include "identities.h"
#include "integrationconstant.h"

#include "SireMaths/complex.h"

#include "SireError/errors.h"
#include "SireCAS/errors.h"

#include "SireStream/datastream.h"

#include <QDebug>

using namespace SireStream;
using namespace SireCAS;

////////////
//////////// Implementation of PowerFunction
////////////

static const RegisterMetaType<PowerFunction> r_powerfunc(MAGIC_ONLY,
                                                         "SireCAS::PowerFunction");

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const PowerFunction &powerfunc)
{
    writeHeader(ds, r_powerfunc, 1) << static_cast<const ExBase&>(powerfunc);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds, PowerFunction &powerfunc)
{
    VersionID v = readHeader(ds, r_powerfunc);

    if (v == 1)
    {
        ds >> static_cast<ExBase&>(powerfunc);
    }
    else
        throw version_error(v, "1", r_powerfunc, CODELOC);

    return ds;
}

/** Return a string representation of this power */
QString PowerFunction::toString() const
{
    Expression pwr = power();

    if (pwr.isZero())
        return "1";
    else
    {
        Expression ex = core();

        if (pwr.isConstant() and ex.isConstant())
            return this->evaluate(ComplexValues()).toString();
        else if (pwr.isConstant())
        {
            Complex pwrval = pwr.evaluate(ComplexValues());

            if (pwrval.isReal())
            {
                double realpwr = pwrval.real();

                if ( SireMaths::areEqual(realpwr,1.0) )
                    return ex.toString();
                else if (ex.isCompound())
                    return QString("[%1]^%2").arg(ex.toString()).arg(realpwr);
                else
                    return QString("%1^%2").arg(ex.toString()).arg(realpwr);
            }
            else if ( SireMaths::isZero(pwrval.real()) )
            {
                if ( ex.isCompound() )
                    return QString("[%1]^%2").arg(ex.toString(),pwrval.toString());
                else
                    return QString("%1^%2").arg(ex.toString(),pwrval.toString());
            }
            else
            {
                if ( ex.isCompound() )
                    return QString("[%1]^[%2]").arg(ex.toString(),pwrval.toString());
                else
                    return QString("%1^[%2]").arg(ex.toString(),pwrval.toString());
            }
        }
        else
        {
            QString ret;
            if (ex.isConstant())
                ret = ex.evaluate(ComplexValues()).toString();
            else if (ex.isCompound())
                ret = QString("[%1]").arg(ex.toString());
            else
                ret = ex.toString();

            if (pwr.isCompound())
                return QString("%1^{%2}").arg(ret, pwr.toString());
            else
                return QString("%1^%2").arg(ret, pwr.toString());
        }
    }
}



/** Return a string representation of this power in the OpenMM syntax */
QString PowerFunction::toOpenMMString() const
{
    Expression pwr = power();

    if (pwr.isZero())
        return "1";
    else
    {
        Expression ex = core();

        if (pwr.isConstant() and ex.isConstant())
            return this->evaluate(ComplexValues()).toString();
        else if (pwr.isConstant())
        {
            Complex pwrval = pwr.evaluate(ComplexValues());

            if (pwrval.isReal())
            {
                double realpwr = pwrval.real();

                if ( SireMaths::areEqual(realpwr,1.0) )
                    return ex.toOpenMMString();
                else if (ex.isCompound())
                    return QString("(%1)^%2").arg(ex.toOpenMMString()).arg(realpwr);
                else
                    return QString("%1^%2").arg(ex.toOpenMMString()).arg(realpwr);
            }
            else if ( SireMaths::isZero(pwrval.real()) )
            {
                if ( ex.isCompound() )
                    return QString("(%1)^%2").arg(ex.toOpenMMString(),pwrval.toString());
                else
                    return QString("%1^%2").arg(ex.toOpenMMString(),pwrval.toString());
            }
            else
            {
                if ( ex.isCompound() )
                    return QString("(%1)^(%2)").arg(ex.toOpenMMString(),pwrval.toString());
                else
                    return QString("%1^(%2)").arg(ex.toOpenMMString(),pwrval.toString());
            }
        }
        else
        {
            QString ret;
            if (ex.isConstant())
                ret = ex.evaluate(ComplexValues()).toString();
            else if (ex.isCompound())
                ret = QString("(%1)").arg(ex.toOpenMMString());
            else
                ret = ex.toString();

            if (pwr.isCompound())
                return QString("%1^(%2)").arg(ret, pwr.toOpenMMString());
            else
                return QString("%1^%2").arg(ret, pwr.toOpenMMString());
        }
    }
}



/** Return this expression with the supplied substitutions */
Expression PowerFunction::substitute(const Identities &identities) const
{
    return Power( core().substitute(identities), power().substitute(identities) ).reduce();
}

/** Return the differential of this expression with respect to 'symbol' */
Expression PowerFunction::differentiate(const Symbol &symbol) const
{
    Expression ex = core();
    Expression pwr = power();

    //f(x)^0 == 1, do d f(x)^0 / dx = 0
    if (pwr.isZero())
        return Expression(0);
    //is the power a function of this symbol?
    else if (pwr.isFunction(symbol))
    {
        // d (a^x) / dx = a^x ln(a)

        if (ex.isFunction(symbol))
        {
            // use identity;
            // f(x)^g(x) = exp( g(x) * ln(f(x)) )

            return Exp( pwr * Ln(ex) ).differentiate(symbol);
        }
        else
            // d a^(f(x)) / dx = a^(f(x)) ln(a) * d f(x) / dx
            return Expression(*this) * Ln(ex) * pwr.differentiate(symbol);
    }
    else
    {
        //differentiate the base
        Expression diff = ex.differentiate(symbol);

        if (diff.isZero())
            return Expression(0);

        //constant power, so d (x^n) / dx = n x^(n-1)
        return diff * pwr * ex.pow( pwr - 1 );
    }
}

/** Return the integral of this power with respect to 'symbol' */
Expression PowerFunction::integrate(const Symbol &symbol) const
{
    Expression pwr = power();
    Expression cre = core();

    bool cre_is_func = cre.isFunction(symbol);
    bool pwr_is_func = pwr.isFunction(symbol);

    if (cre_is_func and not pwr_is_func)
    {
        //power of form f(x)^n,  integral is 1/(f'(x)*(n+1)) * f(x)^n+1
        Expression newpower = pwr + 1;
        return Power(cre,newpower).reduce() / (newpower * cre.diff(symbol));
    }
    else if (pwr_is_func and not cre_is_func)
    {
        //what is integral of n^f(x)?
        return ExBase::integrate(symbol);
    }
    else if (cre_is_func and pwr_is_func)
    {
        //what is the integral of f(x)^g(x)?
        return ExBase::integrate(symbol);
    }
    else
        //this expression is constant with respect to symbol
        return *this * symbol;
}

/** Return the child expressions of this Power - this contains the core() and the power() */
Expressions PowerFunction::children() const
{
    Expressions expressions(core());
    expressions.append(power());
    return expressions;
}

/** Return whether this is a function of 'symbol' */
bool PowerFunction::isFunction(const Symbol &symbol) const
{
    return core().isFunction(symbol) or power().isFunction(symbol);
}

/** Return whether or not this is a constant */
bool PowerFunction::isConstant() const
{
    Expression pwr = power();

    return pwr.isConstant() and (pwr.isZero() or core().isConstant());
}

/** Reduce this Power to a simplified expression (if possible) */
Expression PowerFunction::reduce() const
{
    Expression ex = core();
    Expression pwr = power();

    if (ex.base().isA<PowerFunction>())
    {
        const PowerFunction &powerfunc = ex.base().asA<PowerFunction>();

        if (ex.factor() == 1)
        {
            //this is a power of a power. We can combine the powers together
            return Power( powerfunc.core(), pwr * powerfunc.power() ).reduce();
        }
        else
        {
            //this is ex.factor()^pwr * (power of a power)
            return PowerConstant(ex.factor(),pwr).reduce() *
                   Power(powerfunc.core(), pwr * powerfunc.power() ).reduce();
        }
    }
    else if (ex.base().isA<IntegrationConstant>())
    {
        if (pwr.isConstant())
            //IntegrationConstant ^ const_power  =  IntegrationConstant
            return ex.base();
        else
            return *this;
    }
    else
    {
        bool core_is_constant = ex.isConstant();
        bool power_is_constant = pwr.isConstant();

        if (core_is_constant and power_is_constant)
        {
            //this is a constant
            return Expression( this->evaluate(ComplexValues()) );
        }
        else if (power_is_constant)
        {
            //get the power as an irrational number
            Complex pwrval = pwr.evaluate(ComplexValues());

            if (pwrval.isZero())
            {
                return Expression(1);
            }
            else if (not pwrval.isReal())
            {
                return ComplexPower(ex, pwrval);
            }
            else
            {
                double realpower = pwrval.real();
                if (SireMaths::areEqual(realpower,1.0))
                {
                    return ex;
                }
                else if (SireMaths::isInteger(realpower))
                {
                    return IntegerPower(ex, int(realpower));
                }
                else if (SireMaths::isRational(realpower))
                {
                    return RationalPower(ex, SireMaths::toRational(realpower));
                }
                else
                {
                    return RealPower(ex, realpower);
                }
            }
        }
        else if (core_is_constant)
        {
            Complex val = ex.evaluate(ComplexValues());

            if (val.isReal())
            {
                double realval = val.real();

                if (SireMaths::areEqual(realval, SireMaths::e))
                    return Expression( Exp(pwr) );
                else
                    return Expression( PowerConstant(realval, pwr) );
            }
            else
                return *this;
        }
        else
        {
            return *this;
        }
    }
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

QList<Factor> PowerFunction::expand(const Symbol &symbol) const
{
    QList<Factor> ret;

    if (not this->isFunction(symbol))
    {
        ret.append( Factor(symbol, *this, 0) );
        return ret;
    }

    //get the core
    Expression ex = core();
    
    //expand the core in powers of this symbol
    QList<Factor> core_factors = ex.expand(symbol);

    //get the power...
    Expression pwr = power();

    //we cannot expand if the power is a function of the symbol
    if (pwr.isFunction(symbol))
        throw SireCAS::rearrangement_error( QObject::tr(
            "The expression %1 cannot be rearranged in terms of powers of %2 "
            "as this power is itself a function of %2.")
                .arg(this->toString(), symbol.toString()), CODELOC );

    //if there is only a single power of the symbol, then it doesn't
    //matter what the power is
    if (core_factors.count() == 1)
    {
        core_factors[0] = Factor( core_factors[0].symbol(),
                                  core_factors[0].factor(),
                                  pwr * core_factors[0].power() );
                                  
        return core_factors;
    }

    //we can only expand if the power is a constant integer...
    Complex pwrval = pwr.evaluate(ComplexValues());

    if (pwrval.isZero())
    {
        //this expression is just equal to 1
        ret.append( Factor(symbol,1,0) );
        return ret;
    }
    else if (pwrval.isReal() and SireMaths::isInteger(pwrval.real()))
    {
        int int_power = int(pwrval.real());

        if (int_power == 1)
            return core_factors;
        else if (int_power == 0)
        {
            QList<Factor> ret;
            ret.append( Factor(symbol,1,0) );
            return ret;
        }
        else if (int_power < 0)
        {
            if (core_factors.count() > 1)
                throw SireCAS::rearrangement_error( QObject::tr(
                    "The expression %1 cannot be rearranged in terms of powers of %2 "
                    "because multiple powers of %2 are raised to a negative power.")
                        .arg(this->toString(), symbol.toString()), CODELOC );
                        
            core_factors[0] = Factor(  core_factors[0].symbol(),
                                       core_factors[0].factor(), 
                                       int_power * core_factors[0].power() );
                                      
            return core_factors;
        }
        else
        {
            QList<Factor> ret = core_factors;
            
            for (int i=1; i<int_power; ++i)
            {
                ret = ::multiply(ret, core_factors);
            }
            
            return ret;
        }
    }
    else
        throw SireCAS::rearrangement_error( QObject::tr(
            "The expression %1 cannot be rearranged in terms of powers of %2 "
            "because it involves a non-integer power of a compound expression "
            "involving %2.")
                .arg(this->toString(), symbol.toString()), CODELOC );
    
    return ret;
}

////////////
//////////// Implementation of Power
////////////

static const RegisterMetaType<Power> r_power;

/** Serialise a Power to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const Power &power)
{
    writeHeader(ds, r_power, 1)
          << power.ex << power.pwr << static_cast<const PowerFunction&>(power);

    return ds;
}

/** Deserialise a Power from a binary datastream */
QDataStream &operator>>(QDataStream &ds, Power &power)
{
    VersionID v = readHeader(ds, r_power);

    if (v == 1)
    {
        ds >> power.ex >> power.pwr >> static_cast<PowerFunction&>(power);
    }
    else
        throw version_error(v, "1", r_power, CODELOC);

    return ds;
}

/** Construct a null power */
Power::Power() : PowerFunction()
{}

/** Construct a power that represents core^power */
Power::Power(const Expression &core, const Expression &power)
      : PowerFunction(), ex(core), pwr(power)
{}

/** Copy constructor */
Power::Power(const Power &other)
      : PowerFunction(), ex(other.ex), pwr(other.pwr)
{}

/** Destructor */
Power::~Power()
{}

/** Comparison operator */
bool Power::operator==(const ExBase &other) const
{
    const Power *other_power = dynamic_cast<const Power*>(&other);

    return other_power != 0 and typeid(other).name() == typeid(*this).name()
               and ex == other_power->ex and pwr == other_power->pwr;
}

/** Return a hash for this power */
uint Power::hash() const
{
    return ( r_power.magicID() << 16 ) | ( ex.hash() & 0x0000FF00 )
                | ( pwr.hash() & 0x000000FF );
}

/** Evaluate this power - this could be dodgy for negative bases with
    non-integer powers */
double Power::evaluate(const Values &values) const
{
    double pwrval = pwr.evaluate(values);
    if (SireMaths::isZero(pwrval))
        return 1.0;
    else
    {
        double val = ex.evaluate(values);
        if (SireMaths::isZero(val))
            return val;
        else
            return SireMaths::pow(val, pwrval);
    }
}

/** Evaluate this power - this could be dodgy for negative bases with
    non-integer powers */
Complex Power::evaluate(const ComplexValues &values) const
{
    Complex pwrval = pwr.evaluate(values);

    if (pwrval.isZero())
        return Complex(1);
    else
    {
        Complex val = ex.evaluate(values);
        if (val.isZero())
            return val;
        else
            return SireMaths::pow(val, pwrval);
    }
}

const char* Power::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Power>() );
}

Power* Power::clone() const
{
    return new Power(*this);
}

