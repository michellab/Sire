/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#include "conditional.h"
#include "functions.h"
#include "symbols.h"
#include "identities.h"
#include "expressions.h"
#include "values.h"
#include "complexvalues.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireError/errors.h"
#include "SireCAS/errors.h"
#include "SireMaths/errors.h"

using namespace SireCAS;
using namespace SireStream;

///////////
/////////// Implementation of Condition
///////////

static const RegisterMetaType<Condition> r_condition(MAGIC_ONLY,Condition::typeName());

/** Serialise to a binary datastream */
QDataStream SIRECAS_EXPORT &operator<<(QDataStream &ds, const Condition &condition)
{
    writeHeader(ds, r_condition, 1);
    
    SharedDataStream sds(ds);
    
    sds << condition.lhs << condition.rhs 
        << static_cast<const ExBase&>(condition);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIRECAS_EXPORT &operator>>(QDataStream &ds, Condition &condition)
{
    VersionID v = readHeader(ds, r_condition);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> condition.lhs >> condition.rhs
            >> static_cast<ExBase&>(condition);
    }
    else
        throw version_error( v, "1", r_condition, CODELOC );
        
    return ds;
}

/** Null constructor */
Condition::Condition() : ExBase()
{}

/** Construct with the passed left and right hand side of the condition */
Condition::Condition(const Expression &lhs_, const Expression &rhs_)
          : ExBase(), lhs(lhs_), rhs(rhs_)
{}

/** Copy constructor */
Condition::Condition(const Condition &other)
          : ExBase(other), lhs(other.lhs), rhs(other.rhs)
{}

/** Destructor */
Condition::~Condition()
{}

/** Copy assignment operator */
Condition& Condition::operator=(const Condition &other)
{
    lhs = other.lhs;
    rhs = other.rhs;
    ExBase::operator=(other);
    
    return *this;
}

/** Comparison operator */
bool Condition::operator==(const Condition &other) const
{
    return this == &other or
           (lhs == other.lhs and rhs == other.rhs);
}

/** Differentiate this expression */
Expression Condition::differentiate(const Symbol &symbol) const
{
    throw SireCAS::unavailable_differential( QObject::tr(
              "No differential is available for the condition \"%1\"")
                   .arg(this->toString()), CODELOC );

    return Expression();
}

/** Integrate this expression */
Expression Condition::integrate(const Symbol &symbol) const
{
    throw SireCAS::unavailable_integral( QObject::tr(
            "No integral is available for the condition \"%1\"")
                .arg(this->toString()), CODELOC );
                    
    return Expression();
}

/** Simplify this condition */
Expression Condition::simplify(int options) const
{
    SharedPolyPointer<Condition> ret = *this;
    
    ret->lhs = lhs.simplify(options);
    ret->rhs = rhs.simplify(options);
    
    if (ret->alwaysTrue())
        return AlwaysTrue();
    else if (ret->alwaysFalse())
        return AlwaysFalse();
    else
        return *ret;
}

/** The complex conjugate of a condition is the condition */
Expression Condition::conjugate() const
{
    return *this;
}

/** Return whether or not this function is a function of 'symbol' */
bool Condition::isFunction(const Symbol &symbol) const
{
    return lhs.isFunction(symbol) or rhs.isFunction(symbol);
}

/** Return whether or not this is a constant */
bool Condition::isConstant() const
{
    return (lhs.isConstant() and rhs.isConstant()) or (lhs == rhs);
}

/** Return whether or not this is a complex condition */
bool Condition::isComplex() const
{
    return lhs.isComplex() or rhs.isComplex();
}

/** Return whether or not this is a compound expression */
bool Condition::isCompound() const
{
    return not this->isConstant();
}

/** Return whether or not this is null */
bool Condition::isNull() const
{
    return false;
}

/** Hash this condition */
uint Condition::hash() const
{
    return lhs.hash() + rhs.hash();
}

/** Substitute in the passed identities */
Expression Condition::substitute(const Identities &identities) const
{
    SharedPolyPointer<Condition> ret( *this );
    
    ret->lhs = lhs.substitute(identities);
    ret->rhs = rhs.substitute(identities);
    
    return *ret;
}

/** Return the symbols used in this expression */
Symbols Condition::symbols() const
{
    if (this->alwaysTrue())
        return lhs.symbols();
    else if (this->alwaysFalse())
        return rhs.symbols();
    else
        return lhs.symbols() + rhs.symbols();
}

/** Return the functions in this expression */
Functions Condition::functions() const
{
    return lhs.functions() + rhs.functions();
}

/** Return all of the child expressions in this condition */
Expressions Condition::children() const
{
    Expressions exps;
    exps.append(lhs);
    exps.append(rhs);
    
    return exps;
}

/** Expand this condition into factors of the passed symbol */
QList<Factor> Condition::expand(const Symbol &symbol) const
{
    throw SireCAS::rearrangement_error( QObject::tr(
        "The conditional expression \"%1\" cannot be expanded in terms "
        "of %2.")
            .arg(this->toString(), symbol.toString()), CODELOC );
                
    return QList<Factor>();
}

/** Return the left hand side of this condition */
const Expression& Condition::leftHandSide() const
{
    return lhs;
}

/** Return the right hand side of this condition */
const Expression& Condition::rightHandSide() const
{
    return rhs;
}

/** Return a string representation of this expression */
QString Condition::toString() const
{
    return QString("%1 %2 %3")
                .arg(lhs.toString(), this->operatorString(), rhs.toString());
}

/** Evalute this condition, returning whether or not it is true or false */
bool Condition::evaluateCondition(const Values &values) const
{
    double lhs_value = lhs.evaluate(values);
    double rhs_value = rhs.evaluate(values);
    
    return this->compareValues(lhs_value, rhs_value);
}

/** Evalute this condition, returning whether or not it is true or false */
bool Condition::evaluateCondition(const ComplexValues &values) const
{
    Complex lhs_value = lhs.evaluate(values);
    Complex rhs_value = rhs.evaluate(values);
    
    return this->compareValues(lhs_value, rhs_value);
}

/** Evaluate this expression - this returns '1' if it is true, else it
    returns '0' if it is false */
double Condition::evaluate(const Values &values) const
{
    if (this->evaluateCondition(values))
        return 1;
    else
        return 0;
}

/** Evaluate this expression - this returns '1' if it is true, else
    it returns '0' if it is false */
Complex Condition::evaluate(const ComplexValues &values) const
{
    if (this->evaluateCondition(values))
        return Complex(1);
    else
        return Complex(0);
}

///////////
/////////// Implementation of GreaterThan
///////////

static const RegisterMetaType<GreaterThan> r_greaterthan;

/** Serialise to a binary datastream */
QDataStream SIRECAS_EXPORT &operator<<(QDataStream &ds, const GreaterThan &greaterthan)
{
    writeHeader(ds, r_greaterthan, 1);
    
    ds << static_cast<const Condition&>(greaterthan);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIRECAS_EXPORT &operator>>(QDataStream &ds, GreaterThan &greaterthan)
{
    VersionID v = readHeader(ds, r_greaterthan);
    
    if (v == 1)
    {
        ds >> static_cast<Condition&>(greaterthan);
    }
    else
        throw version_error( v, "1", r_greaterthan, CODELOC );
        
    return ds;
}

/** Constructor */
GreaterThan::GreaterThan() : Condition()
{}

/** Construct to compare 'left_hand_side' with 'right_hand_side' */
GreaterThan::GreaterThan(const Expression &left_hand_side,
                         const Expression &right_hand_side)
            : Condition(left_hand_side, right_hand_side)
{}

/** Copy constructor */
GreaterThan::GreaterThan(const GreaterThan &other)
            : Condition(other)
{}

/** Destructor */
GreaterThan::~GreaterThan()
{}

/** Copy assignment operator */
GreaterThan& GreaterThan::operator=(const GreaterThan &other)
{
    Condition::operator=(other);
    return *this;
}

/** Comparison operator */
bool GreaterThan::operator==(const GreaterThan &other) const
{
    return Condition::operator==(other);
}

/** Comparison operator */
bool GreaterThan::operator==(const ExBase &other) const
{
    const GreaterThan *other_t  = dynamic_cast<const GreaterThan*>(&other);
    
    if (other_t)
        return GreaterThan::operator==(*other_t);
    else
        return false;
}

const char* GreaterThan::typeName()
{
    return QMetaType::typeName( qMetaTypeId<GreaterThan>() );
}

/** Compare the values 'val0' and 'val1' */
bool GreaterThan::compareValues(double val0, double val1) const
{
    return val0 > val1;
}

/** Compare the values 'val0' and 'val1' 

    \throw SireMaths::domain_error
*/
bool GreaterThan::compareValues(const Complex &val0, const Complex &val1) const
{
    if (val0 == val1)
        return false;
    
    if (val0.isReal() and val1.isReal())
        return compareValues( val0.real(), val1.real() );
    
    else if (val0.isPurelyComplex() and val1.isPurelyComplex())
        return compareValues( val0.imag(), val1.imag() );
        
    else
        throw SireMaths::domain_error( QObject::tr(
                "Cannot compare the values of complex numbers %1 and %2.")
                    .arg(val0.toString(), val1.toString()), CODELOC );
    
    return false;
}

/** Return the string representation of this operator */
QString GreaterThan::operatorString() const
{
    return QString(">");
}

/** Return whether or not we can be absolutely sure that this
    condition will always be true. Note that this doesn't try
    too hard, so some things that are always true may not
    be reported here as being always true, e.g. x + 1 > x */
bool GreaterThan::alwaysTrue() const
{
    if (leftHandSide() == rightHandSide())
        return false;
    else if (this->isConstant())
        return this->evaluateCondition( ComplexValues() );
    else
        return false;
}

/** Return whether or not we can be absolutely sure that this
    condition will always be false. Note that this doesn't try
    too hard, so some things that are always false may not
    be reported here as being always false, e.g. x > x + 1 */
bool GreaterThan::alwaysFalse() const
{
    if (leftHandSide() == rightHandSide())
        return true;
    else if (this->isConstant())
        return not (this->evaluateCondition( ComplexValues() ));
    else
        return false;
}

///////////
/////////// Implementation of LessThan
///////////

static const RegisterMetaType<LessThan> r_lessthan;

/** Serialise to a binary datastream */
QDataStream SIRECAS_EXPORT &operator<<(QDataStream &ds, const LessThan &lessthan)
{
    writeHeader(ds, r_lessthan, 1);
    
    ds << static_cast<const Condition&>(lessthan);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIRECAS_EXPORT &operator>>(QDataStream &ds, LessThan &lessthan)
{
    VersionID v = readHeader(ds, r_lessthan);
    
    if (v == 1)
    {
        ds >> static_cast<Condition&>(lessthan);
    }
    else
        throw version_error( v, "1", r_lessthan, CODELOC );
        
    return ds;
}

/** Constructor */
LessThan::LessThan() : Condition()
{}

/** Construct to compare 'left_hand_side' with 'right_hand_side' */
LessThan::LessThan(const Expression &left_hand_side,
                   const Expression &right_hand_side)
         : Condition(left_hand_side, right_hand_side)
{}

/** Copy constructor */
LessThan::LessThan(const LessThan &other)
         : Condition(other)
{}

/** Destructor */
LessThan::~LessThan()
{}

/** Copy assignment operator */
LessThan& LessThan::operator=(const LessThan &other)
{
    Condition::operator=(other);
    return *this;
}

/** Comparison operator */
bool LessThan::operator==(const LessThan &other) const
{
    return Condition::operator==(other);
}

/** Comparison operator */
bool LessThan::operator==(const ExBase &other) const
{
    const LessThan *other_t  = dynamic_cast<const LessThan*>(&other);
    
    if (other_t)
        return LessThan::operator==(*other_t);
    else
        return false;
}

const char* LessThan::typeName()
{
    return QMetaType::typeName( qMetaTypeId<LessThan>() );
}

/** Compare the values 'val0' and 'val1' */
bool LessThan::compareValues(double val0, double val1) const
{
    return val0 < val1;
}

/** Compare the values 'val0' and 'val1' 

    \throw SireMaths::domain_error
*/
bool LessThan::compareValues(const Complex &val0, const Complex &val1) const
{
    if (val0 == val1)
        return false;
    
    if (val0.isReal() and val1.isReal())
        return compareValues( val0.real(), val1.real() );
    
    else if (val0.isPurelyComplex() and val1.isPurelyComplex())
        return compareValues( val0.imag(), val1.imag() );
        
    else
        throw SireMaths::domain_error( QObject::tr(
                "Cannot compare the values of complex numbers %1 and %2.")
                    .arg(val0.toString(), val1.toString()), CODELOC );
    
    return false;
}

/** Return the string representation of this operator */
QString LessThan::operatorString() const
{
    return QString("<");
}

/** Return whether or not we can be absolutely sure that this
    condition will always be true. Note that this doesn't try
    too hard, so some things that are always true may not
    be reported here as being always true, e.g. x < x + 1 */
bool LessThan::alwaysTrue() const
{
    if (leftHandSide() == rightHandSide())
        return false;
    else if (this->isConstant())
        return this->evaluateCondition( ComplexValues() );
    else
        return false;
}

/** Return whether or not we can be absolutely sure that this
    condition will always be false. Note that this doesn't try
    too hard, so some things that are always false may not
    be reported here as being always false, e.g. x + 1 < x */
bool LessThan::alwaysFalse() const
{
    if (leftHandSide() == rightHandSide())
        return true;
    else if (this->isConstant())
        return not (this->evaluateCondition( ComplexValues() ));
    else
        return false;
}

///////////
/////////// Implementation of GreaterOrEqualThan
///////////

static const RegisterMetaType<GreaterOrEqualThan> r_gethan;

/** Serialise to a binary datastream */
QDataStream SIRECAS_EXPORT &operator<<(QDataStream &ds, const GreaterOrEqualThan &gethan)
{
    writeHeader(ds, r_gethan, 1);
    
    ds << static_cast<const Condition&>(gethan);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIRECAS_EXPORT &operator>>(QDataStream &ds, GreaterOrEqualThan &gethan)
{
    VersionID v = readHeader(ds, r_gethan);
    
    if (v == 1)
    {
        ds >> static_cast<Condition&>(gethan);
    }
    else
        throw version_error( v, "1", r_gethan, CODELOC );
        
    return ds;
}

/** Constructor */
GreaterOrEqualThan::GreaterOrEqualThan() : Condition()
{}

/** Construct to compare 'left_hand_side' with 'right_hand_side' */
GreaterOrEqualThan::GreaterOrEqualThan(const Expression &left_hand_side,
                                       const Expression &right_hand_side)
                   : Condition(left_hand_side, right_hand_side)
{}

/** Copy constructor */
GreaterOrEqualThan::GreaterOrEqualThan(const GreaterOrEqualThan &other)
                   : Condition(other)
{}

/** Destructor */
GreaterOrEqualThan::~GreaterOrEqualThan()
{}

/** Copy assignment operator */
GreaterOrEqualThan& GreaterOrEqualThan::operator=(const GreaterOrEqualThan &other)
{
    Condition::operator=(other);
    return *this;
}

/** Comparison operator */
bool GreaterOrEqualThan::operator==(const GreaterOrEqualThan &other) const
{
    return Condition::operator==(other);
}

/** Comparison operator */
bool GreaterOrEqualThan::operator==(const ExBase &other) const
{
    const GreaterOrEqualThan *other_t  = dynamic_cast<const GreaterOrEqualThan*>(&other);
    
    if (other_t)
        return GreaterOrEqualThan::operator==(*other_t);
    else
        return false;
}

const char* GreaterOrEqualThan::typeName()
{
    return QMetaType::typeName( qMetaTypeId<GreaterOrEqualThan>() );
}

/** Compare the values 'val0' and 'val1' */
bool GreaterOrEqualThan::compareValues(double val0, double val1) const
{
    return val0 >= val1;
}

/** Compare the values 'val0' and 'val1' 

    \throw SireMaths::domain_error
*/
bool GreaterOrEqualThan::compareValues(const Complex &val0, const Complex &val1) const
{
    if (val0 == val1)
        return true;
    
    if (val0.isReal() and val1.isReal())
        return compareValues( val0.real(), val1.real() );
    
    else if (val0.isPurelyComplex() and val1.isPurelyComplex())
        return compareValues( val0.imag(), val1.imag() );
        
    else
        throw SireMaths::domain_error( QObject::tr(
                "Cannot compare the values of complex numbers %1 and %2.")
                    .arg(val0.toString(), val1.toString()), CODELOC );
    
    return false;
}

/** Return the string representation of this operator */
QString GreaterOrEqualThan::operatorString() const
{
    return QString(">=");
}

/** Return whether or not we can be absolutely sure that this
    condition will always be true. Note that this doesn't try
    too hard, so some things that are always true may not
    be reported here as being always true, e.g. x + 1 > x */
bool GreaterOrEqualThan::alwaysTrue() const
{
    if (leftHandSide() == rightHandSide())
        return true;
    else if (this->isConstant())
        return this->evaluateCondition( ComplexValues() );
    else
        return false;
}

/** Return whether or not we can be absolutely sure that this
    condition will always be false. Note that this doesn't try
    too hard, so some things that are always false may not
    be reported here as being always false, e.g. x > x + 1 */
bool GreaterOrEqualThan::alwaysFalse() const
{
    if (leftHandSide() == rightHandSide())
        return false;
    else if (this->isConstant())
        return not (this->evaluateCondition( ComplexValues() ));
    else
        return false;
}

///////////
/////////// Implementation of LessOrEqualThan
///////////

static const RegisterMetaType<LessOrEqualThan> r_lethan;

/** Serialise to a binary datastream */
QDataStream SIRECAS_EXPORT &operator<<(QDataStream &ds, const LessOrEqualThan &lethan)
{
    writeHeader(ds, r_lethan, 1);
    
    ds << static_cast<const Condition&>(lethan);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIRECAS_EXPORT &operator>>(QDataStream &ds, LessOrEqualThan &lethan)
{
    VersionID v = readHeader(ds, r_lethan);
    
    if (v == 1)
    {
        ds >> static_cast<Condition&>(lethan);
    }
    else
        throw version_error( v, "1", r_lethan, CODELOC );
        
    return ds;
}

/** Constructor */
LessOrEqualThan::LessOrEqualThan() : Condition()
{}

/** Construct to compare 'left_hand_side' with 'right_hand_side' */
LessOrEqualThan::LessOrEqualThan(const Expression &left_hand_side,
                                 const Expression &right_hand_side)
                : Condition(left_hand_side, right_hand_side)
{}

/** Copy constructor */
LessOrEqualThan::LessOrEqualThan(const LessOrEqualThan &other)
                : Condition(other)
{}

/** Destructor */
LessOrEqualThan::~LessOrEqualThan()
{}

/** Copy assignment operator */
LessOrEqualThan& LessOrEqualThan::operator=(const LessOrEqualThan &other)
{
    Condition::operator=(other);
    return *this;
}

/** Comparison operator */
bool LessOrEqualThan::operator==(const LessOrEqualThan &other) const
{
    return Condition::operator==(other);
}

/** Comparison operator */
bool LessOrEqualThan::operator==(const ExBase &other) const
{
    const LessOrEqualThan *other_t  = dynamic_cast<const LessOrEqualThan*>(&other);
    
    if (other_t)
        return LessOrEqualThan::operator==(*other_t);
    else
        return false;
}

const char* LessOrEqualThan::typeName()
{
    return QMetaType::typeName( qMetaTypeId<LessOrEqualThan>() );
}

/** Compare the values 'val0' and 'val1' */
bool LessOrEqualThan::compareValues(double val0, double val1) const
{
    return val0 <= val1;
}

/** Compare the values 'val0' and 'val1' 

    \throw SireMaths::domain_error
*/
bool LessOrEqualThan::compareValues(const Complex &val0, const Complex &val1) const
{
    if (val0 == val1)
        return true;
    
    if (val0.isReal() and val1.isReal())
        return compareValues( val0.real(), val1.real() );
    
    else if (val0.isPurelyComplex() and val1.isPurelyComplex())
        return compareValues( val0.imag(), val1.imag() );
        
    else
        throw SireMaths::domain_error( QObject::tr(
                "Cannot compare the values of complex numbers %1 and %2.")
                    .arg(val0.toString(), val1.toString()), CODELOC );
    
    return false;
}

/** Return the string representation of this operator */
QString LessOrEqualThan::operatorString() const
{
    return QString("<=");
}

/** Return whether or not we can be absolutely sure that this
    condition will always be true. Note that this doesn't try
    too hard, so some things that are always true may not
    be reported here as being always true, e.g. x + 1 > x */
bool LessOrEqualThan::alwaysTrue() const
{
    if (leftHandSide() == rightHandSide())
        return true;
    else if (this->isConstant())
        return this->evaluateCondition( ComplexValues() );
    else
        return false;
}

/** Return whether or not we can be absolutely sure that this
    condition will always be false. Note that this doesn't try
    too hard, so some things that are always false may not
    be reported here as being always false, e.g. x > x + 1 */
bool LessOrEqualThan::alwaysFalse() const
{
    if (leftHandSide() == rightHandSide())
        return false;
    else if (this->isConstant())
        return not (this->evaluateCondition( ComplexValues() ));
    else
        return false;
}

///////////
/////////// Implementation of EqualTo
///////////

static const RegisterMetaType<EqualTo> r_equalto;

/** Serialise to a binary datastream */
QDataStream SIRECAS_EXPORT &operator<<(QDataStream &ds, const EqualTo &equalto)
{
    writeHeader(ds, r_equalto, 1);
    
    ds << static_cast<const Condition&>(equalto);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIRECAS_EXPORT &operator>>(QDataStream &ds, EqualTo &equalto)
{
    VersionID v = readHeader(ds, r_equalto);
    
    if (v == 1)
    {
        ds >> static_cast<Condition&>(equalto);
    }
    else
        throw version_error( v, "1", r_equalto, CODELOC );
        
    return ds;
}

/** Constructor */
EqualTo::EqualTo() : Condition()
{}

/** Construct to compare 'left_hand_side' with 'right_hand_side' */
EqualTo::EqualTo(const Expression &left_hand_side,
                 const Expression &right_hand_side)
        : Condition(left_hand_side, right_hand_side)
{}

/** Copy constructor */
EqualTo::EqualTo(const EqualTo &other)
        : Condition(other)
{}

/** Destructor */
EqualTo::~EqualTo()
{}

/** Copy assignment operator */
EqualTo& EqualTo::operator=(const EqualTo &other)
{
    Condition::operator=(other);
    return *this;
}

/** Comparison operator */
bool EqualTo::operator==(const EqualTo &other) const
{
    return Condition::operator==(other);
}

/** Comparison operator */
bool EqualTo::operator==(const ExBase &other) const
{
    const EqualTo *other_t  = dynamic_cast<const EqualTo*>(&other);
    
    if (other_t)
        return EqualTo::operator==(*other_t);
    else
        return false;
}

const char* EqualTo::typeName()
{
    return QMetaType::typeName( qMetaTypeId<EqualTo>() );
}

/** Compare the values 'val0' and 'val1' */
bool EqualTo::compareValues(double val0, double val1) const
{
    return val0 == val1;
}

/** Compare the values 'val0' and 'val1' 

    \throw SireMaths::domain_error
*/
bool EqualTo::compareValues(const Complex &val0, const Complex &val1) const
{
    if (val0 == val1)
        return true;
    
    if (val0.isReal() and val1.isReal())
        return compareValues( val0.real(), val1.real() );
    
    else if (val0.isPurelyComplex() and val1.isPurelyComplex())
        return compareValues( val0.imag(), val1.imag() );
        
    else
        throw SireMaths::domain_error( QObject::tr(
                "Cannot compare the values of complex numbers %1 and %2.")
                    .arg(val0.toString(), val1.toString()), CODELOC );
    
    return false;
}

/** Return the string representation of this operator */
QString EqualTo::operatorString() const
{
    return QString("==");
}

/** Return whether or not we can be absolutely sure that this
    condition will always be true. Note that this doesn't try
    too hard, so some things that are always true may not
    be reported here as being always true, e.g. x + 1 > x */
bool EqualTo::alwaysTrue() const
{
    if (leftHandSide() == rightHandSide())
        return true;
    else if (this->isConstant())
        return this->evaluateCondition( ComplexValues() );
    else
        return false;
}

/** Return whether or not we can be absolutely sure that this
    condition will always be false. Note that this doesn't try
    too hard, so some things that are always false may not
    be reported here as being always false, e.g. x > x + 1 */
bool EqualTo::alwaysFalse() const
{
    if (leftHandSide() == rightHandSide())
        return false;
    else if (this->isConstant())
        return not (this->evaluateCondition( ComplexValues() ));
    else
        return false;
}

///////////
/////////// Implementation of NotEqualTo
///////////

static const RegisterMetaType<NotEqualTo> r_notequalto;

/** Serialise to a binary datastream */
QDataStream SIRECAS_EXPORT &operator<<(QDataStream &ds, const NotEqualTo &notequalto)
{
    writeHeader(ds, r_notequalto, 1);
    
    ds << static_cast<const Condition&>(notequalto);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIRECAS_EXPORT &operator>>(QDataStream &ds, NotEqualTo &notequalto)
{
    VersionID v = readHeader(ds, r_notequalto);
    
    if (v == 1)
    {
        ds >> static_cast<Condition&>(notequalto);
    }
    else
        throw version_error( v, "1", r_notequalto, CODELOC );
        
    return ds;
}

/** Constructor */
NotEqualTo::NotEqualTo() : Condition()
{}

/** Construct to compare 'left_hand_side' with 'right_hand_side' */
NotEqualTo::NotEqualTo(const Expression &left_hand_side,
                       const Expression &right_hand_side)
           : Condition(left_hand_side, right_hand_side)
{}

/** Copy constructor */
NotEqualTo::NotEqualTo(const NotEqualTo &other)
           : Condition(other)
{}

/** Destructor */
NotEqualTo::~NotEqualTo()
{}

/** Copy assignment operator */
NotEqualTo& NotEqualTo::operator=(const NotEqualTo &other)
{
    Condition::operator=(other);
    return *this;
}

/** Comparison operator */
bool NotEqualTo::operator==(const NotEqualTo &other) const
{
    return Condition::operator==(other);
}

/** Comparison operator */
bool NotEqualTo::operator==(const ExBase &other) const
{
    const NotEqualTo *other_t  = dynamic_cast<const NotEqualTo*>(&other);
    
    if (other_t)
        return NotEqualTo::operator==(*other_t);
    else
        return false;
}

const char* NotEqualTo::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NotEqualTo>() );
}

/** Compare the values 'val0' and 'val1' */
bool NotEqualTo::compareValues(double val0, double val1) const
{
    return val0 != val1;
}

/** Compare the values 'val0' and 'val1' 

    \throw SireMaths::domain_error
*/
bool NotEqualTo::compareValues(const Complex &val0, const Complex &val1) const
{
    if (val0 == val1)
        return false;
    
    if (val0.isReal() and val1.isReal())
        return compareValues( val0.real(), val1.real() );
    
    else if (val0.isPurelyComplex() and val1.isPurelyComplex())
        return compareValues( val0.imag(), val1.imag() );
        
    else
        throw SireMaths::domain_error( QObject::tr(
                "Cannot compare the values of complex numbers %1 and %2.")
                    .arg(val0.toString(), val1.toString()), CODELOC );
    
    return false;
}

/** Return the string representation of this operator */
QString NotEqualTo::operatorString() const
{
    return QString("!=");
}

/** Return whether or not we can be absolutely sure that this
    condition will always be true. Note that this doesn't try
    too hard, so some things that are always true may not
    be reported here as being always true, e.g. x + 1 > x */
bool NotEqualTo::alwaysTrue() const
{
    if (leftHandSide() == rightHandSide())
        return false;
    else if (this->isConstant())
        return this->evaluateCondition( ComplexValues() );
    else
        return false;
}

/** Return whether or not we can be absolutely sure that this
    condition will always be false. Note that this doesn't try
    too hard, so some things that are always false may not
    be reported here as being always false, e.g. x > x + 1 */
bool NotEqualTo::alwaysFalse() const
{
    if (leftHandSide() == rightHandSide())
        return true;
    else if (this->isConstant())
        return not (this->evaluateCondition( ComplexValues() ));
    else
        return false;
}

///////////
/////////// Implementation of AlwaysTrue
///////////

static const RegisterMetaType<AlwaysTrue> r_alwaystrue;

/** Serialise to a binary datastream */
QDataStream SIRECAS_EXPORT &operator<<(QDataStream &ds, const AlwaysTrue &alwaystrue)
{
    writeHeader(ds, r_alwaystrue, 1);
    
    ds << static_cast<const Condition&>(alwaystrue);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIRECAS_EXPORT &operator>>(QDataStream &ds, AlwaysTrue &alwaystrue)
{
    VersionID v = readHeader(ds, r_alwaystrue);
    
    if (v == 1)
    {
        ds >> static_cast<Condition&>(alwaystrue);
    }
    else
        throw version_error( v, "1", r_alwaystrue, CODELOC );
        
    return ds;
}

/** Constructor */
AlwaysTrue::AlwaysTrue() : Condition()
{}

/** Copy constructor */
AlwaysTrue::AlwaysTrue(const AlwaysTrue &other) : Condition(other)
{}

/** Destructor */
AlwaysTrue::~AlwaysTrue()
{}

/** Copy assignment operator */
AlwaysTrue& AlwaysTrue::operator=(const AlwaysTrue &other)
{
    Condition::operator=(other);
    return *this;
}

/** Comparison operator */
bool AlwaysTrue::operator==(const AlwaysTrue&) const
{
    return true;
}

/** Comparison operator */
bool AlwaysTrue::operator==(const ExBase &other) const
{
    return dynamic_cast<const AlwaysTrue*>(&other) != 0;
}

const char* AlwaysTrue::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AlwaysTrue>() );
}

/** Return a string representation of truth */
QString AlwaysTrue::toString() const
{
    return QObject::tr( "True" );
}
    
/** This cannot be further simplified */
Expression AlwaysTrue::simplify(int) const
{
    return *this;
}

/** AlwaysTrue is always true */
bool AlwaysTrue::alwaysTrue() const
{
    return true;
}

/** AlwaysTrue is never false */
bool AlwaysTrue::alwaysFalse() const
{
    return false;
}

/** This is not a function of anything */
bool AlwaysTrue::isFunction(const Symbol&) const
{
    return false;
}

/** Truth is always constant */
bool AlwaysTrue::isConstant() const
{
    return true;
}

/** Truth is never complex */
bool AlwaysTrue::isComplex() const
{
    return false;
}

/** Truth is always simple */
bool AlwaysTrue::isCompound() const
{
    return false;
}

/** Truth is never empty */
bool AlwaysTrue::isNull() const
{
    return false;
}

/** Hash truth */
uint AlwaysTrue::hash() const
{
    return 1;
}

/** There is no substituting the truth */
Expression AlwaysTrue::substitute(const Identities&) const
{
    return *this;
}

/** There are no symbols in truth */
Symbols AlwaysTrue::symbols() const
{
    return Symbols();
}

/** There are no functions in truth */
Functions AlwaysTrue::functions() const
{
    return Functions();
}

/** The truth has no children */
Expressions AlwaysTrue::children() const
{
    return Expressions();
}

/** The truth cannot be expanded */
QList<Factor> AlwaysTrue::expand(const Symbol &symbol) const
{
    QList<Factor> ret;
    ret.append( Factor(symbol, *this, 0) );
    
    return ret;
}

/** Truth is always true */
bool AlwaysTrue::evaluateCondition(const Values&) const
{
    return true;
}

/** Truth is always true */
bool AlwaysTrue::evaluateCondition(const ComplexValues&) const
{
    return true;
}

/** Truth is always true */
double AlwaysTrue::evaluate(const Values&) const
{
    return 1;
}

/** Truth is always true */
Complex AlwaysTrue::evaluate(const ComplexValues&) const
{
    return Complex(1);
}

QString AlwaysTrue::operatorString() const
{
    return QString();
}

bool AlwaysTrue::compareValues(double, double) const
{
    return true;
}

bool AlwaysTrue::compareValues(const Complex&, const Complex&) const
{
    return true;
}

///////////
/////////// Implementation of AlwaysFalse
///////////

static const RegisterMetaType<AlwaysFalse> r_alwaysfalse;

/** Serialise to a binary datastream */
QDataStream SIRECAS_EXPORT &operator<<(QDataStream &ds, const AlwaysFalse &alwaysfalse)
{
    writeHeader(ds, r_alwaysfalse, 1);
    
    ds << static_cast<const Condition&>(alwaysfalse);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIRECAS_EXPORT &operator>>(QDataStream &ds, AlwaysFalse &alwaysfalse)
{
    VersionID v = readHeader(ds, r_alwaysfalse);
    
    if (v == 1)
    {
        ds >> static_cast<Condition&>(alwaysfalse);
    }
    else
        throw version_error( v, "1", r_alwaysfalse, CODELOC );
        
    return ds;
}

/** Constructor */
AlwaysFalse::AlwaysFalse() : Condition()
{}

/** Copy constructor */
AlwaysFalse::AlwaysFalse(const AlwaysFalse &other) : Condition(other)
{}

/** Destructor */
AlwaysFalse::~AlwaysFalse()
{}

/** Copy assignment operator */
AlwaysFalse& AlwaysFalse::operator=(const AlwaysFalse &other)
{
    Condition::operator=(other);
    return *this;
}

/** Comparison operator */
bool AlwaysFalse::operator==(const AlwaysFalse&) const
{
    return true;
}

/** Comparison operator */
bool AlwaysFalse::operator==(const ExBase &other) const
{
    return dynamic_cast<const AlwaysFalse*>(&other) != 0;
}

const char* AlwaysFalse::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AlwaysFalse>() );
}

/** Return a string representation of truth */
QString AlwaysFalse::toString() const
{
    return QObject::tr( "False" );
}
    
/** This cannot be further simplified */
Expression AlwaysFalse::simplify(int) const
{
    return *this;
}

/** AlwaysFalse is never true */
bool AlwaysFalse::alwaysTrue() const
{
    return false;
}

/** AlwaysFalse is always false */
bool AlwaysFalse::alwaysFalse() const
{
    return true;
}

/** This is not a function of anything */
bool AlwaysFalse::isFunction(const Symbol&) const
{
    return false;
}

/** Truth is always constant */
bool AlwaysFalse::isConstant() const
{
    return true;
}

/** False is never complex */
bool AlwaysFalse::isComplex() const
{
    return false;
}

/** False is always simple */
bool AlwaysFalse::isCompound() const
{
    return false;
}

/** False is never empty */
bool AlwaysFalse::isNull() const
{
    return false;
}

/** Hash false */
uint AlwaysFalse::hash() const
{
    return 0;
}

/** There is no substituting false */
Expression AlwaysFalse::substitute(const Identities&) const
{
    return *this;
}

/** There are no symbols in false */
Symbols AlwaysFalse::symbols() const
{
    return Symbols();
}

/** There are no functions in false */
Functions AlwaysFalse::functions() const
{
    return Functions();
}

/** False has no children */
Expressions AlwaysFalse::children() const
{
    return Expressions();
}

/** False cannot be expanded */
QList<Factor> AlwaysFalse::expand(const Symbol &symbol) const
{
    QList<Factor> ret;
    ret.append( Factor(symbol, *this, 0) );
    
    return ret;
}

/** False is never true */
bool AlwaysFalse::evaluateCondition(const Values&) const
{
    return false;
}

/** False is never true */
bool AlwaysFalse::evaluateCondition(const ComplexValues&) const
{
    return false;
}

/** False is never true */
double AlwaysFalse::evaluate(const Values&) const
{
    return 0;
}

/** False is never true */
Complex AlwaysFalse::evaluate(const ComplexValues&) const
{
    return Complex(0);
}

QString AlwaysFalse::operatorString() const
{
    return QString();
}

bool AlwaysFalse::compareValues(double, double) const
{
    return false;
}

bool AlwaysFalse::compareValues(const Complex&, const Complex&) const
{
    return false;
}

///////////
/////////// Implementation of Conditional
///////////

static const RegisterMetaType<Conditional> r_conditional;

/** Serialise to a binary datastream */
QDataStream SIRECAS_EXPORT &operator<<(QDataStream &ds, const Conditional &conditional)
{
    writeHeader(ds, r_conditional, 1);
    
    SharedDataStream sds(ds);
    
    sds << conditional.cond
        << conditional.true_expression << conditional.false_expression
        << static_cast<const ExBase&>(conditional);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIRECAS_EXPORT &operator>>(QDataStream &ds, Conditional &conditional)
{
    VersionID v = readHeader(ds, r_conditional);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> conditional.cond
            >> conditional.true_expression >> conditional.false_expression
            >> static_cast<ExBase&>(conditional);
    }
    else
        throw version_error( v, "1", r_conditional, CODELOC );
        
    return ds;
}

/** Null constructor */
Conditional::Conditional() : ExBase(), cond(AlwaysTrue())
{}

/** Construct a conditional where if 'condition' is true, then 
    'true_expression' is evaluated, while if 'condition' is false,
    then 'false_expression' is evaluated */
Conditional::Conditional(const Condition &condition, 
                         const Expression &t_expression,
                         const Expression &f_expression)
            : ExBase(),
              cond(condition),
              true_expression(t_expression),
              false_expression(f_expression)
{
    if (condition.isConstant())
    {
        cond = AlwaysTrue();
    
        if (condition.alwaysTrue())
        {
            false_expression = Expression();
        }
        else
        {
            true_expression = false_expression;
            false_expression = Expression();
        }
    }
}
            
/** Copy constructor */
Conditional::Conditional(const Conditional &other)
            : ExBase(other), cond(other.cond),
              true_expression(other.true_expression),
              false_expression(other.false_expression)
{}

/** Destructor */
Conditional::~Conditional()
{}

/** Copy assignment operator */
Conditional& Conditional::operator=(const Conditional &other)
{
    if (this != &other)
    {
        cond = other.cond;
        true_expression = other.true_expression;
        false_expression = other.false_expression;
    }

    return *this;
}

/** Comparison operator */
bool Conditional::operator==(const Conditional &other) const
{
    return this == &other or
           ( (cond == other.cond or 
                *(static_cast<const ExBase*>(cond.constData())) == 
                *(static_cast<const ExBase*>(other.cond.constData()))) and 
            true_expression == other.true_expression and
            false_expression == other.false_expression );
}

/** Comparison operator */
bool Conditional::operator==(const ExBase &other) const
{
    const Conditional *othercond = dynamic_cast<const Conditional*>(&other);
    
    if (othercond)
        return Conditional::operator==(*othercond);
    else
        return false;
}

const char* Conditional::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Conditional>() );
}

/** Return the condition */
const Condition& Conditional::condition() const
{
    BOOST_ASSERT( cond.constData() != 0 );
    return *cond;
}

/** Return the expression to be evaluated if the condition is true */
const Expression& Conditional::trueExpression() const
{
    return true_expression;
}

/** Return the expression to be evaluated if the condition is false */
const Expression& Conditional::falseExpression() const
{
    return false_expression;
}

/** Return the series expansion of this product with respect to 'symbol', to order 'n'*/
Expression Conditional::series(const Symbol &symbol, int n) const
{
    if (condition().alwaysTrue())
        return true_expression.series(symbol, n);
    else
        return *this;
}

/** Try to simplify this condition */
Expression Conditional::simplify(int options) const
{
    if (condition().isConstant())
    {
        if (condition().alwaysTrue())
            return true_expression.simplify(options);
        else
            return false_expression.simplify(options);
    }
    else
    {
        Conditional ret;
        
        Expression simple_cond = cond->simplify(options);
        
        if (simple_cond.isConstant())
        {
            if (simple_cond.isZero())
                return false_expression.simplify(options);
            else
                return true_expression.simplify(options);
        }
        
        if (not simple_cond.base().isA<Condition>())
            throw SireError::program_bug( QObject::tr(
                    "How have we lost the condition %1 and got %2 instead???")
                        .arg(cond->toString(), simple_cond.toString()), CODELOC );
        
        ret.cond = simple_cond.base().asA<Condition>();
        ret.true_expression = true_expression.simplify(options);
        ret.false_expression = false_expression.simplify(options);
        
        return ret;
    }
}

/** Return the complex conjugate of this expression */
Expression Conditional::conjugate() const
{
    if (condition().alwaysTrue())
        return true_expression.conjugate();
    else
    {
        Conditional ret;
        
        ret.cond = cond;
        ret.true_expression = true_expression.conjugate();
        ret.false_expression = false_expression.conjugate();
        
        return ret;
    }
}

/** Return the differential of this expression */
Expression Conditional::differentiate(const Symbol &symbol) const
{
    if (condition().alwaysTrue())
        return true_expression.differentiate(symbol);
    else
    {
        if (true_expression.isFunction(symbol) or 
            false_expression.isFunction(symbol))
        {
            throw SireCAS::unavailable_differential( QObject::tr(
                    "Cannot differentiate the condition %1 with respect to %2.")
                        .arg(this->toString(), symbol.toString()), CODELOC );
        }
        
        return 0;
    }
}

/** Return the integral of this expression */
Expression Conditional::integrate(const Symbol &symbol) const
{
    if (condition().alwaysTrue())
        return true_expression.integrate(symbol);
    else
    {
        if (true_expression.isFunction(symbol) or 
            false_expression.isFunction(symbol))
        {
            throw SireCAS::unavailable_integral( QObject::tr(
                    "Cannot integrate the condition %1 with respect to %2.")
                        .arg(this->toString(), symbol.toString()), CODELOC );
        }
        
        return *this * symbol;
    }
}

/** Return whether or not this is a function of 'symbol' */
bool Conditional::isFunction(const Symbol &symbol) const
{
    return true_expression.isFunction(symbol) or 
           false_expression.isFunction(symbol);
}

/** Return whether or not this is constant */
bool Conditional::isConstant() const
{
    return cond->isConstant() and 
           true_expression.isConstant() and false_expression.isConstant();
}

/** Is this a complex expression? */
bool Conditional::isComplex() const
{
    return condition().isComplex() or 
           true_expression.isComplex() or false_expression.isComplex();
}

/** Is this a compound expression? */
bool Conditional::isCompound() const
{
    if (condition().alwaysTrue())
        return true_expression.isCompound();
    else
        return true;
}

/** Hash this conditional */
uint Conditional::hash() const
{
    return cond->hash() + true_expression.hash() + false_expression.hash();
}

/** Return a string representation of this conditional */
QString Conditional::toString() const
{
    if (condition().alwaysTrue())
    {
        return true_expression.toString();
    }
    else
    {
        return QObject::tr( "if ( %1 ){ %2 ? %3 }" )
                    .arg(cond->toString(),
                         true_expression.toString(),
                         false_expression.toString());
    }
}

/** Return whether or not this is null */
bool Conditional::isNull() const
{
    return condition().alwaysTrue() and true_expression.isZero();
}

/** Evaluate this expression for the passed values */
double Conditional::evaluate(const Values &values) const
{
    if (cond->evaluateCondition(values))
    {
        return true_expression.evaluate(values);
    }
    else
    {
        return false_expression.evaluate(values);
    }
}

/** Evaluate this expresion for the passed values */
Complex Conditional::evaluate(const ComplexValues &values) const
{
    if (cond->evaluateCondition(values))
    {
        return true_expression.evaluate(values);
    }
    else
    {
        return false_expression.evaluate(values);
    }
}

/** Substitute 'identities' into this expression */
Expression Conditional::substitute(const Identities &identities) const
{
    Conditional ret;

    Expression new_cond = cond->substitute(identities);
    
    if (new_cond.base().isA<Condition>())
        ret.cond = new_cond.base().asA<Condition>();
    else
    {
        if (new_cond.isConstant())
        {
            if (new_cond.isZero())
                ret.cond = AlwaysFalse();
            else
                ret.cond = AlwaysTrue();
        }
        else
            throw SireError::program_bug( QObject::tr(
                    "How did the condition %1 turn into %2???")
                        .arg(cond->toString(), new_cond.toString()), CODELOC );
    }

    if (ret.cond->isConstant())
    {
        if (ret.cond->alwaysTrue())
        {
            ret.true_expression = true_expression.substitute(identities);
            ret.false_expression = Expression();
        }
        else
        {
            ret.true_expression = false_expression.substitute(identities);
            ret.false_expression = Expression();
        }
        
        ret.cond = AlwaysTrue();
    }
    else
    {
        ret.true_expression = true_expression.substitute(identities);
        ret.false_expression = false_expression.substitute(identities);
    }
    
    return ret;
}

/** Return the symbols used in this expression */
Symbols Conditional::symbols() const
{
    if (condition().alwaysTrue())
        return true_expression.symbols();
    else
    {
        return condition().symbols() + 
               true_expression.symbols() + false_expression.symbols();
    }
}

/** Return the functions used in this expression */
Functions Conditional::functions() const
{
    if (condition().alwaysTrue())
        return true_expression.functions();
    else
    {
        return condition().functions() + 
               true_expression.functions() + false_expression.functions();
    }
}

/** Return the children of this expression */
Expressions Conditional::children() const
{
    if (condition().alwaysTrue())
        return true_expression.children();
    else
    {
        Expressions ret = condition().children();
        ret.append(true_expression);
        ret.append(false_expression);
        
        return ret;
    }
}

/** Expand this expression in terms of 'symbol' */
QList<Factor> Conditional::expand(const Symbol &symbol) const
{
    if (condition().alwaysTrue())
        return true_expression.expand(symbol);
        
    else if (this->isFunction(symbol))
        throw SireCAS::rearrangement_error( QObject::tr(
                 "Cannot expand the conditional %1 in terms of %2.")
                    .arg(this->toString(), symbol.toString()), CODELOC );
                    
    QList<Factor> ret;
    
    ret.append( Factor(symbol, *this, 0) );
    
    return ret;
}

Conditional* Conditional::clone() const
{
    return new Conditional(*this);
}


LessOrEqualThan* LessOrEqualThan::clone() const
{
    return new LessOrEqualThan(*this);
}


GreaterOrEqualThan* GreaterOrEqualThan::clone() const
{
    return new GreaterOrEqualThan(*this);
}


GreaterThan* GreaterThan::clone() const
{
    return new GreaterThan(*this);
}


EqualTo* EqualTo::clone() const
{
    return new EqualTo(*this);
}


AlwaysTrue* AlwaysTrue::clone() const
{
    return new AlwaysTrue(*this);
}


LessThan* LessThan::clone() const
{
    return new LessThan(*this);
}


AlwaysFalse* AlwaysFalse::clone() const
{
    return new AlwaysFalse(*this);
}


NotEqualTo* NotEqualTo::clone() const
{
    return new NotEqualTo(*this);
}

