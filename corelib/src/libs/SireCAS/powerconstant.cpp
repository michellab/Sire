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

#include "powerconstant.h"
#include "expression.h"
#include "expressions.h"
#include "symbol.h"
#include "symbols.h"
#include "values.h"
#include "complexvalues.h"
#include "identities.h"
#include "integrationconstant.h"

#include "SireStream/datastream.h"

#include <QDebug>

using namespace SireStream;
using namespace SireCAS;

//////////
////////// Implementation of PowerConstant
//////////

static const RegisterMetaType<PowerConstant> r_powerconstant;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const PowerConstant &power)
{
    writeHeader(ds, r_powerconstant, 1)
          << power.cre << power.pwr << static_cast<const PowerFunction&>(power);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds, PowerConstant &power)
{
    VersionID v = readHeader(ds, r_powerconstant);

    if (v == 1)
    {
        ds >> power.cre >> power.pwr >> static_cast<PowerFunction&>(power);
    }
    else
        throw version_error(v, "1", r_powerconstant, CODELOC);

    return ds;
}

/** Create a null PowerConstant (0^0) */
PowerConstant::PowerConstant() : PowerFunction(), cre(0), pwr(1)
{}

/** Construct the PowerConstant val^power */
PowerConstant::PowerConstant(double val, const Expression &power)
              : PowerFunction(), cre(val), pwr(power)
{}

/** Copy constructor */
PowerConstant::PowerConstant(const PowerConstant &other)
              : PowerFunction(), cre(other.cre), pwr(other.pwr)
{}

/** Destructor */
PowerConstant::~PowerConstant()
{}

/** Comparison operator */
bool PowerConstant::operator==(const ExBase &other) const
{
    const PowerConstant *other_power = dynamic_cast<const PowerConstant*>(&other);

    return other_power != 0 and typeid(other).name() == typeid(*this).name()
                 and cre == other_power->cre and pwr == other_power->pwr;
}

/** Return a hash of this power */
uint PowerConstant::hash() const
{
    return ( r_powerconstant.magicID() <<16 ) | ( pwr.hash() & 0x0000FFFF );
}

/** Evaluate this function */
double PowerConstant::evaluate(const Values &values) const
{
    return SireMaths::pow( cre, pwr.evaluate(values) );
}

/** Evaluate this function */
Complex PowerConstant::evaluate(const ComplexValues &values) const
{
    return SireMaths::pow( cre, pwr.evaluate(values) );
}

//////////
////////// Implementation of ConstantPower
//////////

//register a pure virtual class
static const RegisterMetaType<ConstantPower> r_constantpower(MAGIC_ONLY,
                                                             "SireCAS::ConstantPower");

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const ConstantPower &power)
{
    writeHeader(ds, r_constantpower, 1)
          << power.ex << static_cast<const PowerFunction&>(power);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds, ConstantPower &power)
{
    VersionID v = readHeader(ds, r_constantpower);

    if (v == 1)
    {
        ds >> power.ex >> static_cast<PowerFunction&>(power);
    }
    else
        throw version_error(v, "1", r_constantpower, CODELOC);

    return ds;
}

/** Return a hash of this power */
uint ConstantPower::hash() const
{
    return ( r_constantpower.magicID() <<16 ) | ( ex.hash() & 0x0000FFFF );
}

//////////
////////// Implementation of IntegerPower
//////////

static const RegisterMetaType<IntegerPower> r_integerpower;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const IntegerPower &power)
{
    writeHeader(ds, r_integerpower, 1) << power.pwr
                                       << static_cast<const ConstantPower&>(power);
    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds, IntegerPower &power)
{
    VersionID v = readHeader(ds, r_integerpower);

    if (v == 1)
    {
        ds >> power.pwr >> static_cast<ConstantPower&>(power);
    }
    else
        throw version_error(v, "1", r_integerpower, CODELOC);

    return ds;
}

/** Null constructor */
IntegerPower::IntegerPower() : ConstantPower(), pwr(0)
{}

/** Construct expression^power */
IntegerPower::IntegerPower(const Expression &expression, int power)
             : ConstantPower(expression), pwr(power)
{}

/** Copy constructor */
IntegerPower::IntegerPower(const IntegerPower &other)
             : ConstantPower(other), pwr(other.pwr)
{}

/** Destructor */
IntegerPower::~IntegerPower()
{}

/** Comparison operator */
bool IntegerPower::operator==(const ExBase &other) const
{
    const IntegerPower *other_power = dynamic_cast<const IntegerPower*>(&other);

    return other_power != 0 and typeid(other).name() == typeid(*this).name()
                 and pwr == other_power->pwr and ex == other_power->ex;
}

/** Return a hash of this power */
uint IntegerPower::hash() const
{
    return ( r_integerpower.magicID() <<16 ) | ( ex.hash() & 0x0000FFFF );
}

/** Evaluate this power */
double IntegerPower::evaluate(const Values &values) const
{
    return SireMaths::pow( ex.evaluate(values), pwr );
}

/** Evaluate this power */
Complex IntegerPower::evaluate(const ComplexValues &values) const
{
    return SireMaths::pow( ex.evaluate(values), pwr );
}

//////////
////////// Implementation of RationalPower
//////////

static const RegisterMetaType<RationalPower> r_rationalpower;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const RationalPower &power)
{
    writeHeader(ds, r_rationalpower, 1) << power.pwr
                                        << static_cast<const ConstantPower&>(power);
    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds, RationalPower &power)
{
    VersionID v = readHeader(ds, r_rationalpower);

    if (v == 1)
    {
        ds >> power.pwr >> static_cast<ConstantPower&>(power);
    }
    else
        throw version_error(v, "1", r_rationalpower, CODELOC);

    return ds;
}

/** Null constructor */
RationalPower::RationalPower() : ConstantPower(), pwr(0)
{}

/** Construct expression^power */
RationalPower::RationalPower(const Expression &expression, const Rational &power)
              : ConstantPower(expression), pwr(power)
{}

/** Copy constructor */
RationalPower::RationalPower(const RationalPower &other)
              : ConstantPower(other), pwr(other.pwr)
{}

/** Destructor */
RationalPower::~RationalPower()
{}

/** Comparison operator */
bool RationalPower::operator==(const ExBase &other) const
{
    const RationalPower *other_power = dynamic_cast<const RationalPower*>(&other);

    return other_power != 0 and typeid(other).name() == typeid(*this).name()
                 and pwr == other_power->pwr and ex == other_power->ex;
}

/** Return a hash of this power */
uint RationalPower::hash() const
{
    return ( r_rationalpower.magicID() <<16 ) | ( ex.hash() & 0x0000FFFF );
}

/** Evaluate this power */
double RationalPower::evaluate(const Values &values) const
{
    return SireMaths::pow( ex.evaluate(values), pwr );
}

/** Evaluate this power */
Complex RationalPower::evaluate(const ComplexValues &values) const
{
    return SireMaths::pow( ex.evaluate(values), pwr );
}

//////////
////////// Implementation of RealPower
//////////

static const RegisterMetaType<RealPower> r_realpower;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const RealPower &power)
{
    writeHeader(ds, r_realpower, 1) << power.pwr
                                    << static_cast<const ConstantPower&>(power);
    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds, RealPower &power)
{
    VersionID v = readHeader(ds, r_realpower);

    if (v == 1)
    {
        ds >> power.pwr >> static_cast<ConstantPower&>(power);
    }
    else
        throw version_error(v, "1", r_realpower, CODELOC);

    return ds;
}

/** Null constructor */
RealPower::RealPower() : ConstantPower(), pwr(0)
{}

/** Construct expression^power */
RealPower::RealPower(const Expression &expression, double power)
          : ConstantPower(expression), pwr(power)
{}

/** Copy constructor */
RealPower::RealPower(const RealPower &other)
          : ConstantPower(other), pwr(other.pwr)
{}

/** Destructor */
RealPower::~RealPower()
{}

/** Comparison operator */
bool RealPower::operator==(const ExBase &other) const
{
    const RealPower *other_power = dynamic_cast<const RealPower*>(&other);

    return other_power != 0 and typeid(other).name() == typeid(*this).name()
                 and pwr == other_power->pwr and ex == other_power->ex;
}

/** Return a hash of this power */
uint RealPower::hash() const
{
    return ( r_realpower.magicID() <<16 ) | ( ex.hash() & 0x0000FFFF );
}

/** Evaluate this power */
double RealPower::evaluate(const Values &values) const
{
    return SireMaths::pow( ex.evaluate(values), pwr );
}

/** Evaluate this power */
Complex RealPower::evaluate(const ComplexValues &values) const
{
    return SireMaths::pow( ex.evaluate(values), pwr );
}

//////////
////////// Implementation of ComplexPower
//////////

static const RegisterMetaType<ComplexPower> r_complexpower;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const ComplexPower &power)
{
    writeHeader(ds, r_complexpower, 1) << power.pwr
                                       << static_cast<const ConstantPower&>(power);
    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds, ComplexPower &power)
{
    VersionID v = readHeader(ds, r_complexpower);

    if (v == 1)
    {
        ds >> power.pwr >> static_cast<ConstantPower&>(power);
    }
    else
        throw version_error(v, "1", r_complexpower, CODELOC);

    return ds;
}

/** Null constructor */
ComplexPower::ComplexPower() : ConstantPower(), pwr(0)
{}

/** Construct expression^power */
ComplexPower::ComplexPower(const Expression &expression, const Complex &power)
             : ConstantPower(expression), pwr(power)
{}

/** Copy constructor */
ComplexPower::ComplexPower(const ComplexPower &other)
             : ConstantPower(other), pwr(other.pwr)
{}

/** Destructor */
ComplexPower::~ComplexPower()
{}

/** Comparison operator */
bool ComplexPower::operator==(const ExBase &other) const
{
    const ComplexPower *other_power = dynamic_cast<const ComplexPower*>(&other);

    return other_power != 0 and typeid(other).name() == typeid(*this).name()
                 and pwr == other_power->pwr and ex == other_power->ex;
}

/** Return a hash of this power */
uint ComplexPower::hash() const
{
    return ( r_complexpower.magicID() <<16 ) | ( ex.hash() & 0x0000FFFF );
}

/** Evaluate this power */
double ComplexPower::evaluate(const Values &values) const
{
    //calculate the result...
    Complex val = SireMaths::pow( ex.evaluate(values), pwr );

    if (not val.isReal())
        throw SireMaths::domain_error(QObject::tr(
            "Raising the expression \"%1\" to the complex power \"%2\" has "
            "resulted in a complex value, \"%3\"")
                .arg(ex.toString(), pwr.toString(), val.toString()), CODELOC);

    return val.real();
}

/** Evaluate this power */
Complex ComplexPower::evaluate(const ComplexValues &values) const
{
    return SireMaths::pow( ex.evaluate(values), pwr );
}

const char* PowerConstant::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PowerConstant>() );
}

const char* IntegerPower::typeName()
{
    return QMetaType::typeName( qMetaTypeId<IntegerPower>() );
}

const char* RationalPower::typeName()
{
    return QMetaType::typeName( qMetaTypeId<RationalPower>() );
}

const char* RealPower::typeName()
{
    return QMetaType::typeName( qMetaTypeId<RealPower>() );
}

const char* ComplexPower::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ComplexPower>() );
}

RealPower* RealPower::clone() const
{
    return new RealPower(*this);
}


ComplexPower* ComplexPower::clone() const
{
    return new ComplexPower(*this);
}


IntegerPower* IntegerPower::clone() const
{
    return new IntegerPower(*this);
}


RationalPower* RationalPower::clone() const
{
    return new RationalPower(*this);
}


PowerConstant* PowerConstant::clone() const
{
    return new PowerConstant(*this);
}

