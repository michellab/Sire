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

#include "hyperbolicfuncs.h"
#include "trigfuncs.h"
#include "exp.h"
#include "identities.h"
#include "expression.h"
#include "complexvalues.h"

#include "SireStream/datastream.h"

using namespace SireStream;
using namespace SireCAS;

////////////
//////////// Register the trig functions
////////////

static const RegisterMetaType<Cosh> r_cosh;
static const RegisterMetaType<Sinh> r_sinh;
static const RegisterMetaType<Tanh> r_tanh;
static const RegisterMetaType<Csch> r_csch;
static const RegisterMetaType<Sech> r_sech;
static const RegisterMetaType<Coth> r_coth;

////////////
//////////// Stream the trig functions
////////////

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const Cosh &cosh)
{
    writeHeader(ds, r_cosh, 1) << static_cast<const SingleFunc&>(cosh);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds, Cosh &cosh)
{
    VersionID v = readHeader(ds, r_cosh);

    if (v == 1)
    {
        ds >> static_cast<SingleFunc&>(cosh);
    }
    else
        throw version_error(v, "1", r_cosh, CODELOC);

    return ds;
}

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const Sinh &sinh)
{
    writeHeader(ds, r_sinh, 1) << static_cast<const SingleFunc&>(sinh);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds, Sinh &sinh)
{
    VersionID v = readHeader(ds, r_sinh);

    if (v == 1)
    {
        ds >> static_cast<SingleFunc&>(sinh);
    }
    else
        throw version_error(v, "1", r_sinh, CODELOC);

    return ds;
}

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const Tanh &tanh)
{
    writeHeader(ds, r_tanh, 1) << static_cast<const SingleFunc&>(tanh);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds, Tanh &tanh)
{
    VersionID v = readHeader(ds, r_tanh);

    if (v == 1)
    {
        ds >> static_cast<SingleFunc&>(tanh);
    }
    else
        throw version_error(v, "1", r_tanh, CODELOC);

    return ds;
}

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const Csch &csch)
{
    writeHeader(ds, r_csch, 1) << static_cast<const SingleFunc&>(csch);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds, Csch &csch)
{
    VersionID v = readHeader(ds, r_csch);

    if (v == 1)
    {
        ds >> static_cast<SingleFunc&>(csch);
    }
    else
        throw version_error(v, "1", r_csch, CODELOC);

    return ds;
}

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const Sech &sech)
{
    writeHeader(ds, r_sech, 1) << static_cast<const SingleFunc&>(sech);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds, Sech &sech)
{
    VersionID v = readHeader(ds, r_sech);

    if (v == 1)
    {
        ds >> static_cast<SingleFunc&>(sech);
    }
    else
        throw version_error(v, "1", r_sech, CODELOC);

    return ds;
}

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const Coth &coth)
{
    writeHeader(ds, r_coth, 1) << static_cast<const SingleFunc&>(coth);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds, Coth &coth)
{
    VersionID v = readHeader(ds, r_coth);

    if (v == 1)
    {
        ds >> static_cast<SingleFunc&>(coth);
    }
    else
        throw version_error(v, "1", r_coth, CODELOC);

    return ds;
}

////////////
//////////// Implementation of Hyperbolic Cosine
////////////

/** Null constructor */
Cosh::Cosh() : SingleFunc()
{}

/** Construct cos(expression) */
Cosh::Cosh(const Expression &expression) : SingleFunc(expression)
{}

/** Create cos(cos(expression)) */
Cosh::Cosh(const Cosh &other) : SingleFunc(other)
{}

/** Destructor */
Cosh::~Cosh()
{}

/** Return the magic */
uint Cosh::magic() const
{
    return r_cosh.magicID();
}

/** Comparison operator */
bool Cosh::operator==(const ExBase &other) const
{
    const Cosh *other_cos = dynamic_cast<const Cosh*>(&other);

    return other_cos != 0 and typeid(other).name() == typeid(*this).name()
                 and this->argument() == other_cos->argument();
}

/** Evaluate this function */
double Cosh::evaluate(const Values &values) const
{
    return std::cosh( x().evaluate(values) );
}

/** Complex evaluation */
Complex Cosh::evaluate(const ComplexValues &values) const
{
    return SireMaths::cosh( x().evaluate(values) );
}

/** The differential of cosh(x) = sinh(x) */
Expression Cosh::diff() const
{
    return Sinh(x());
}

/** Integral of cosh(x) = sinh(x) */
Expression Cosh::integ() const
{
    return Sinh(x());
}

////////////
//////////// Implementation of hyperbolic sine
////////////

/** Null constructor */
Sinh::Sinh() : SingleFunc()
{}

/** Construct cos(expression) */
Sinh::Sinh(const Expression &expression) : SingleFunc(expression)
{}

/** Copy constructor */
Sinh::Sinh(const Sinh &other) : SingleFunc(other)
{}

/** Destructor */
Sinh::~Sinh()
{}

/** Return the magic */
uint Sinh::magic() const
{
    return r_sinh.magicID();
}

/** Comparison operator */
bool Sinh::operator==(const ExBase &other) const
{
    const Sinh *other_cos = dynamic_cast<const Sinh*>(&other);

    return other_cos != 0 and typeid(other).name() == typeid(*this).name()
                 and this->argument() == other_cos->argument();
}

/** Evaluate this function */
double Sinh::evaluate(const Values &values) const
{
    return std::sinh( x().evaluate(values) );
}

/** Complex evaluation */
Complex Sinh::evaluate(const ComplexValues &values) const
{
    return SireMaths::sinh( x().evaluate(values) );
}

/** The differential of sinh(x) = cosh(x) */
Expression Sinh::diff() const
{
    return Cosh(x());
}

/** Integral of sinh(x) = cosh(x) */
Expression Sinh::integ() const
{
    return  -Cosh(x());
}

////////////
//////////// Implementation of hyperbolic tangent
////////////

/** Null constructor */
Tanh::Tanh() : SingleFunc()
{}

/** Construct cos(expression) */
Tanh::Tanh(const Expression &expression) : SingleFunc(expression)
{}

/** Copy constructor */
Tanh::Tanh(const Tanh &other) : SingleFunc(other)
{}

/** Destructor */
Tanh::~Tanh()
{}

/** Return the magic */
uint Tanh::magic() const
{
    return r_tanh.magicID();
}

/** Comparison operator */
bool Tanh::operator==(const ExBase &other) const
{
    const Tanh *other_cos = dynamic_cast<const Tanh*>(&other);

    return other_cos != 0 and typeid(other).name() == typeid(*this).name()
                 and this->argument() == other_cos->argument();
}

/** Evaluate this function */
double Tanh::evaluate(const Values &values) const
{
    return std::tanh( x().evaluate(values) );
}

/** Complex evaluation */
Complex Tanh::evaluate(const ComplexValues &values) const
{
    return SireMaths::tanh( x().evaluate(values) );
}

/** The differential of tanh(x) = sech^2(x) */
Expression Tanh::diff() const
{
    return pow( Sech(x()), 2 );
}

/** Integral of tan(x) = ln [ cosh(x) ] + C */
Expression Tanh::integ() const
{
    return Ln( Cosh(x()) );
}

////////////
//////////// Implementation of hyperbolic cosecant
////////////

/** Null constructor */
Csch::Csch() : SingleFunc()
{}

/** Construct cos(expression) */
Csch::Csch(const Expression &expression) : SingleFunc(expression)
{}

/** Copy constructor */
Csch::Csch(const Csch &other) : SingleFunc(other)
{}

/** Destructor */
Csch::~Csch()
{}

/** Return the magic */
uint Csch::magic() const
{
    return r_csch.magicID();
}

/** Comparison operator */
bool Csch::operator==(const ExBase &other) const
{
    const Csch *other_cos = dynamic_cast<const Csch*>(&other);

    return other_cos != 0 and typeid(other).name() == typeid(*this).name()
                 and this->argument() == other_cos->argument();
}

/** Evaluate this function */
double Csch::evaluate(const Values &values) const
{
    //csch = 1 / sinh
    return double(1.0) / std::sinh( x().evaluate(values) );
}

/** Complex evaluation */
Complex Csch::evaluate(const ComplexValues &values) const
{
    return SireMaths::csch( x().evaluate(values) );
}

/** The differential of csc(x) = -coth(x) csch(x) */
Expression Csch::diff() const
{
    return -Csch(x()) * Coth(x());
}

/** Integral of csc(x) = ln( sinh(x/2) ) - ln( cosh(x/2) )  */
Expression Csch::integ() const
{
    return Ln( Sinh(x()/2) ) - Ln( Cosh(x()/2) );
}

////////////
//////////// Implementation of hyperbolic secant
////////////

/** Null constructor */
Sech::Sech() : SingleFunc()
{}

/** Construct cos(expression) */
Sech::Sech(const Expression &expression) : SingleFunc(expression)
{}

/** Copy constructor */
Sech::Sech(const Sech &other) : SingleFunc(other)
{}

/** Destructor */
Sech::~Sech()
{}

/** Return the magic */
uint Sech::magic() const
{
    return r_sech.magicID();
}

/** Comparison operator */
bool Sech::operator==(const ExBase &other) const
{
    const Sech *other_cos = dynamic_cast<const Sech*>(&other);

    return other_cos != 0 and typeid(other).name() == typeid(*this).name()
                 and this->argument() == other_cos->argument();
}

/** Evaluate this function */
double Sech::evaluate(const Values &values) const
{
    //sech = 1 / cosh
    return double(1.0) / std::cosh( x().evaluate(values) );
}

/** Complex evaluation */
Complex Sech::evaluate(const ComplexValues &values) const
{
    return SireMaths::sech( x().evaluate(values) );
}

/** The differential of sec(x) = -sech(x) tanh(x) */
Expression Sech::diff() const
{
    return -Sech(x()) * Tanh(x());
}

/** Integral of sech(x) = 2 cot( tanh(x/2) ) */
Expression Sech::integ() const
{
    return 2 * Cot( Tanh(x()/2) );
}

////////////
//////////// Implementation of hyperbolic cotangent
////////////

/** Null constructor */
Coth::Coth() : SingleFunc()
{}

/** Construct cos(expression) */
Coth::Coth(const Expression &expression) : SingleFunc(expression)
{}

/** Copy constructor */
Coth::Coth(const Coth &other) : SingleFunc(other)
{}

/** Destructor */
Coth::~Coth()
{}

/** Return the magic */
uint Coth::magic() const
{
    return r_coth.magicID();
}

/** Comparison operator */
bool Coth::operator==(const ExBase &other) const
{
    const Coth *other_cos = dynamic_cast<const Coth*>(&other);

    return other_cos != 0 and typeid(other).name() == typeid(*this).name()
                 and this->argument() == other_cos->argument();
}

/** Evaluate this function */
double Coth::evaluate(const Values &values) const
{
    //coth = 1 / tanh
    return double(1.0) / std::tanh( x().evaluate(values) );
}

/** Complex evaluation */
Complex Coth::evaluate(const ComplexValues &values) const
{
    return SireMaths::coth( x().evaluate(values) );
}

/** The differential of coth(x) = -csch^2(x) */
Expression Coth::diff() const
{
    return -(pow( Csch(x()), 2 ));
}

/** Integral of coth(x) = ln( sinh(x) ) */
Expression Coth::integ() const
{
    return Ln( Sinh(x()) );
}


const char* Cosh::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Cosh>() );
}

const char* Sinh::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Sinh>() );
}

const char* Tanh::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Tanh>() );
}

const char* Csch::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Csch>() );
}

const char* Sech::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Sech>() );
}

const char* Coth::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Coth>() );
}

Cosh* Cosh::clone() const
{
    return new Cosh(*this);
}


Tanh* Tanh::clone() const
{
    return new Tanh(*this);
}


Sinh* Sinh::clone() const
{
    return new Sinh(*this);
}


Coth* Coth::clone() const
{
    return new Coth(*this);
}


Csch* Csch::clone() const
{
    return new Csch(*this);
}


Sech* Sech::clone() const
{
    return new Sech(*this);
}

