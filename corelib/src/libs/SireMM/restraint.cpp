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

#include "restraint.h"

#include "SireMol/moleculedata.h"
#include "SireMol/molecules.h"
#include "SireMol/molid.h"
#include "SireMol/molnum.h"

#include "SireCAS/symbols.h"
#include "SireCAS/values.h"
#include "SireCAS/expression.h"

#include "SireFF/forcetable.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireCAS/errors.h"
#include "SireError/errors.h"

using namespace SireMM;
using namespace SireFF;
using namespace SireMol;
using namespace SireCAS;
using namespace SireBase;
using namespace SireUnits::Dimension;
using namespace SireStream;

////////////
//////////// Implementation of Restraint
////////////

static const RegisterMetaType<Restraint> r_restraint( MAGIC_ONLY, 
                                                      Restraint::typeName() );

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const Restraint &restraint)
{
    writeHeader(ds, r_restraint, 1);
    
    ds << static_cast<const Property&>(restraint);
       
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, Restraint &restraint)
{
    VersionID v = readHeader(ds, r_restraint);
    
    if (v == 1)
    {
        ds >> static_cast<Property&>(restraint);
    }
    else
        throw version_error( v, "1", r_restraint, CODELOC );
        
    return ds;
}

/** Constructor */
Restraint::Restraint() : Property()
{}

/** Copy constructor */
Restraint::Restraint(const Restraint &other) : Property(other)
{}

/** Destructor */
Restraint::~Restraint()
{}

/** Copy assignment operator */
Restraint& Restraint::operator=(const Restraint &other)
{
    if (this != &other)
    {
        Property::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool Restraint::operator==(const Restraint &other) const
{
    return Property::operator==(other);
}

/** Comparison operator */
bool Restraint::operator!=(const Restraint &other) const
{
    return Property::operator!=(other);
}

Q_GLOBAL_STATIC( NullRestraint, nullRestraint )

/** Return the global null restraint */
const NullRestraint& Restraint::null()
{
    return *(nullRestraint());
}

////////////
//////////// Implementation of Restraint3D
////////////

static const RegisterMetaType<Restraint3D> r_restraint3d( MAGIC_ONLY,
                                                          Restraint3D::typeName() );

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds,
                                      const Restraint3D &restraint3d)
{
    writeHeader(ds, r_restraint3d, 1);

    SharedDataStream sds(ds);

    sds << restraint3d.spce
        << static_cast<const Restraint&>(restraint3d);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, Restraint3D &restraint3d)
{
    VersionID v = readHeader(ds, r_restraint3d);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> restraint3d.spce
            >> static_cast<Restraint&>(restraint3d);
    }
    else
        throw version_error( v, "1", r_restraint3d, CODELOC );
        
    return ds;
}

/** Constructor */
Restraint3D::Restraint3D() : Restraint()
{}

/** Copy constructor */
Restraint3D::Restraint3D(const Restraint3D &other)
            : Restraint(), spce(other.spce)
{}

/** Destructor */
Restraint3D::~Restraint3D()
{}

/** Copy assignment operator */
Restraint3D& Restraint3D::operator=(const Restraint3D &other)
{
    Restraint::operator=(other);
    spce = other.spce;
    
    return *this;
}

/** Comparison operator */
bool Restraint3D::operator==(const Restraint3D &other) const
{
    return spce == other.spce and Restraint::operator==(other);
}

/** Comparison operator */
bool Restraint3D::operator!=(const Restraint3D &other) const
{
    return spce != other.spce or Restraint::operator!=(other);
}

/** Return the 3D space in which this restraint operates */
const Space& Restraint3D::space() const
{
    return spce.read();
}

/** Set the 3D space in which this restraint operates */
void Restraint3D::setSpace(const Space &space)
{
    spce = space;
}

///////////
//////////// Implementation of ExpressionRestraint3D
////////////

static const RegisterMetaType<ExpressionRestraint3D> r_exprestraint3d( MAGIC_ONLY,
                                                   ExpressionRestraint3D::typeName() );

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, 
                                      const ExpressionRestraint3D &exprestraint3d)
{
    writeHeader(ds, r_exprestraint3d, 1);
    
    SharedDataStream sds(ds);
    
    sds << exprestraint3d.nrg_expression << exprestraint3d.vals
        << static_cast<const Restraint3D&>(exprestraint3d);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, 
                                      ExpressionRestraint3D &exprestraint3d)
{
    VersionID v = readHeader(ds, r_exprestraint3d);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> exprestraint3d.nrg_expression >> exprestraint3d.vals
            >> static_cast<Restraint3D&>(exprestraint3d);
    }
    else
        throw version_error( v, "1", r_exprestraint3d, CODELOC );
        
    return ds;
}

/** Constructor */
ExpressionRestraint3D::ExpressionRestraint3D() : Restraint3D()
{}

/** Construct to use the passed energy expression, with the supplied  
    user values */
ExpressionRestraint3D::ExpressionRestraint3D(const Expression &expression,
                                             const Values &values)
                      : Restraint3D(), nrg_expression(expression)
{
    if (expression.isConstant())
    {
        nrg_expression = expression.evaluate(Values());
    }
    else
    {
        if (not values.isEmpty())
        {
            //put in any missing values into 'user_vals'
            foreach (Symbol symbol, nrg_expression.symbols())
            {
                if (values.contains(symbol))
                    vals.set( symbol, values[symbol] );
            }
        }
    }
}

/** Copy constructor */
ExpressionRestraint3D::ExpressionRestraint3D(const ExpressionRestraint3D &other)
                      : Restraint3D(other), nrg_expression(other.nrg_expression),
                        vals(other.vals)
{}

/** Destructor */
ExpressionRestraint3D::~ExpressionRestraint3D()
{}

/** Copy assignment operator */
ExpressionRestraint3D& ExpressionRestraint3D::operator=(
                                                const ExpressionRestraint3D &other)
{
    if (this != &other)
    {
        Restraint3D::operator=(other);
        nrg_expression = other.nrg_expression;
        vals = other.vals;
    }
    
    return *this;
}

/** Comparison operator */
bool ExpressionRestraint3D::operator==(const ExpressionRestraint3D &other) const
{
    return this == &other or
           (Restraint3D::operator==(other) and
            nrg_expression == other.nrg_expression and
            vals == other.vals);
}

/** Comparison operator */
bool ExpressionRestraint3D::operator!=(const ExpressionRestraint3D &other) const
{
    return not ExpressionRestraint3D::operator==(other);
}

/** Return a string representation of this restraint */
QString ExpressionRestraint3D::toString() const
{
    if (vals.isEmpty())
        return QString("%1( %2 )").arg(this->what()).arg(nrg_expression.toString());
    else
        return QString("%1( %2 ; %3 )").arg( this->what() )
                                       .arg( nrg_expression.toString() )
                                       .arg( vals.toString() );
}

/** Return the function used to evaluate the restraint */
const Expression& ExpressionRestraint3D::restraintFunction() const
{
    return nrg_expression;
}

/** Internal function used to set the value of the passed symbol to 'value' */
void ExpressionRestraint3D::_pvt_setValue(const Symbol &symbol, double value)
{
    vals.set(symbol, value);
}

/** Set the value of 'symbol' to 'value'. Nothing is done if this
    symbol is not in the passed expression. 
    
    \throw SireError::invalid_arg
*/
void ExpressionRestraint3D::setValue(const Symbol &symbol, double value)
{
    if (nrg_expression.isFunction(symbol))
    {
        if (this->builtinSymbols().contains(symbol))
            throw SireError::invalid_arg( QObject::tr(
                "You cannot set the value of the symbol %1 to %2 in "
                "the restraint %3 as this symbol is one of the built-in "
                "unchangable symbols of this restraint (%4).")
                    .arg(symbol.toString()).arg(value)
                    .arg(this->toString())
                    .arg( Sire::toString(this->builtinSymbols()) ), CODELOC );

        vals.set( symbol, value );
    }
}

/** Return the value of the symbol 'symbol'

    \throw SireCAS::missing_symbol
*/
double ExpressionRestraint3D::getValue(const Symbol &symbol) const
{
    if (not nrg_expression.isFunction(symbol))
    {
        throw SireCAS::missing_symbol( QObject::tr(
                "There is no symbol %1 in the restraint %2. Available symbols "
                "are %3.")
                    .arg(symbol.toString(), this->toString())
                    .arg(Sire::toString(this->symbols())), CODELOC );
    }
    
    return vals[symbol];
}

/** Return whether or not this restraint has a value for the symbol 'symbol' */
bool ExpressionRestraint3D::hasValue(const Symbol &symbol) const
{
    return nrg_expression.isFunction(symbol);
}

/** Return all of the symbols available in this restraint */
Symbols ExpressionRestraint3D::symbols() const
{
    return nrg_expression.symbols();
}

/** Return all of the symbols in this expression that can be set by the user */
Symbols ExpressionRestraint3D::userSymbols() const
{
    return this->symbols() - this->builtinSymbols();
}

/** Return all of the values set in the restraint expression */
Values ExpressionRestraint3D::values() const
{
    return vals;
}

/** Return all of the values set by the user in this restraint expression */
Values ExpressionRestraint3D::userValues() const
{
    Values ret;
    
    foreach (Symbol symbol, this->userSymbols())
    {
        ret.set( symbol, vals[symbol] );
    }
    
    return ret;
}

/** Return the current energy of this restraint */
MolarEnergy ExpressionRestraint3D::energy() const
{
    return MolarEnergy( nrg_expression.evaluate(vals) );
}

////////////
//////////// Implementation of NullRestraint
////////////

static const RegisterMetaType<NullRestraint> r_nullrestraint;

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, 
                                      const NullRestraint &nullrestraint)
{
    writeHeader(ds, r_nullrestraint, 1);
    
    ds << static_cast<const Restraint3D&>(nullrestraint);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, NullRestraint &nullrestraint)
{
    VersionID v = readHeader(ds, r_nullrestraint);
    
    if (v == 1)
    {
        ds >> static_cast<Restraint3D&>(nullrestraint);
    }
    else
        throw version_error( v, "1", r_nullrestraint, CODELOC );
        
    return ds;
}

/** Constructor */
NullRestraint::NullRestraint() : ConcreteProperty<NullRestraint,Restraint3D>()
{}

/** Copy constructor */
NullRestraint::NullRestraint(const NullRestraint &other)
              : ConcreteProperty<NullRestraint,Restraint3D>(other)
{}

/** Destructor */
NullRestraint::~NullRestraint()
{}

/** Copy assignment operator */
NullRestraint& NullRestraint::operator=(const NullRestraint &other)
{
    Restraint3D::operator=(other);
    return *this;
}

/** Comparison operator */
bool NullRestraint::operator==(const NullRestraint &other) const
{
    return Restraint3D::operator==(other);
}

/** Comparison operator */
bool NullRestraint::operator!=(const NullRestraint &other) const
{
    return Restraint3D::operator!=(other);
}

const char* NullRestraint::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullRestraint>() );
}

/** Return a string representation of this restraint */
QString NullRestraint::toString() const
{
    return QObject::tr("NullRestraint()");
}

/** The null restraint has no energy */
MolarEnergy NullRestraint::energy() const
{
    return MolarEnergy(0);
}

/** The null restraint will not change the force */
void NullRestraint::force(MolForceTable&, double) const
{}

/** The null restraint will not change the force */
void NullRestraint::force(ForceTable&, double) const
{}

/** The null restraint cannot be updated */
void NullRestraint::update(const MoleculeData&)
{}

/** The null restraint cannot be updated */
void NullRestraint::update(const Molecules&)
{}

/** There are no molecules in the NullRestraint */
Molecules NullRestraint::molecules() const
{
    return Molecules();
}

/** There are no molecules in the NullRestraint */
bool NullRestraint::contains(MolNum) const
{
    return false;
}

/** There are no molecules in the NullRestraint */
bool NullRestraint::contains(const MolID&) const
{
    return false;
}
    
/** There are no molecules in the NullRestraint */
bool NullRestraint::usesMoleculesIn(const ForceTable &forcetable) const
{
    return false;
}

/** There are no molecules in the NullRestraint */
bool NullRestraint::usesMoleculesIn(const Molecules &molecules) const
{
    return false;
}

/** Set the value of the symbol 'symbol' in this restraint to 'value'.
    This does nothing if this symbol is not used in this restraint */
void NullRestraint::setValue(const Symbol&, double)
{}

/** Return the value of the symbol 'symbol' in this restraint. This
    raises an exception if this symbol is not used
    
    \throw SireCAS::missing_symbol
*/
double NullRestraint::getValue(const Symbol &symbol) const
{
    throw SireCAS::missing_symbol( QObject::tr(
            "The NullRestraint class does not use the symbol %1.")
                .arg(symbol.toString()), CODELOC );
                
    return 0;
}

/** Return whether or not this restraint has a value for the symbol 'symbol' */
bool NullRestraint::hasValue(const Symbol&) const
{
    return false;
}

/** Return all of the symbols uses by this restraint */
Symbols NullRestraint::symbols() const
{
    return Symbols();
}

/** Return the symbols that can be set by the user */
Symbols NullRestraint::userSymbols() const
{
    return Symbols();
}

/** Return the symbols that are built into this restraint */
Symbols NullRestraint::builtinSymbols() const
{
    return Symbols();
}

/** Return all of the values of all of the symbols used in this restraint */
Values NullRestraint::values() const
{
    return Values();
}

/** Return the values that have been supplied by the user */
Values NullRestraint::userValues() const
{
    return Values();
}

/** Return the values that are built into this restraint */
Values NullRestraint::builtinValues() const
{
    return Values();
}

/** Return the differential of this restraint with respect to the 
    symbol 'symbol'
    
    \throw SireCAS::unavailable_differential
*/
RestraintPtr NullRestraint::differentiate(const Symbol&) const
{
    return *this;
}
