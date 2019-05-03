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

#include <QMutex>
#include <QHash>

#include "symbol.h"
#include "symbols.h"
#include "functions.h"
#include "expressions.h"
#include "identities.h"
#include "complexvalues.h"
#include "values.h"

#include "SireCAS/errors.h"

#include "SireStream/datastream.h"

#include <QDebug>

using namespace SireStream;
using namespace SireCAS;

////////
//////// Implementation of Factor
////////

Factor::Factor() : f(1), p(1)
{}

Factor::Factor(const Symbol &symbol,
               const Expression &factor, const Expression &power)
       : s(symbol), f(factor), p(power)
{}

Factor::Factor(const Symbol &symbol, 
               double factor, double power) 
       : s(symbol), f(factor), p(power)
{}

Factor::Factor(const Factor &other) : s(other.s), f(other.f), p(other.p)
{}

Factor::~Factor()
{}

Factor& Factor::operator=(const Factor &other)
{
    f = other.f;
    p = other.p;
    s = other.s;
    
    return *this;
}

bool Factor::operator==(const Factor &other) const
{
    return f == other.f and p == other.p and s == other.s;
}

bool Factor::operator!=(const Factor &other) const
{
    return f != other.f or p != other.p or s != other.s;
}

static QString get_string(const Expression &expression)
{
    if (expression.isCompound())
    {
        return QString("[%1]").arg(expression.toString());
    }
    else
        return expression.toString();
}

QString Factor::toString() const
{
    return QString("%1 times { %2 }^%3")
                    .arg(::get_string(f), s.toString(), ::get_string(p));
}

////////
//////// Implementation of Symbol
////////

typedef struct
{
QHash<QString,SymbolID> name2id;
SymbolID lastid;
} SymbolReg;

static SymbolReg *registry = 0;

static QMutex global_reg_mutex;

static SymbolReg& getRegistry()
{
    if (registry == 0)
    {
        registry = new SymbolReg();
        registry->lastid = 0;
    }

    return *registry;
}

/** Return an ID for the symbol with representation 'rep'. This
    creates a new ID if there is no symbol currently registered with
    this ID. */
SymbolID Symbol::getNewID(const QString &rep)
{
    if (rep.isNull() or rep.isEmpty())
        return 0;

    QMutexLocker lkr(&global_reg_mutex);

    SymbolReg &registry = getRegistry();

    if (registry.name2id.contains(rep))
        return registry.name2id.value(rep);
    else
    {
        registry.lastid++;
        registry.name2id.insert( rep, registry.lastid );
        return registry.lastid;
    }
}

/** Return the name of the symbol with ID == symid

    \throw SireCAS::invalid_symbol
*/
QString Symbol::getName(SymbolID symid)
{
    QMutexLocker lkr(&global_reg_mutex);

    SymbolReg &registry = getRegistry();

    for (QHash<QString,SymbolID>::const_iterator it = registry.name2id.constBegin();
         it != registry.name2id.constEnd();
         ++it)
    {
        if (*it == symid)
            return it.key();
    }

    throw SireCAS::invalid_symbol( QObject::tr(
              "There is no symbol with ID == %1.")
                  .arg(symid), CODELOC );
}

static const RegisterMetaType<Symbol> r_symbol;

/** Hash a symbol */
uint Symbol::hash() const
{
    //return (r_symbol.magicID() << 16) | (id & 0x0000FFFF);
    return id;
}

/** Serialise a Symbol to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const Symbol &sym)
{
    writeHeader(ds, r_symbol, 1) << sym.stringrep
                                 << static_cast<const ExBase&>(sym);

    return ds;
}

/** Deserialise a Symbol from a binary datastream */
QDataStream &operator>>(QDataStream &ds, Symbol &sym)
{
    VersionID v = readHeader(ds, r_symbol);

    if (v == 1)
    {
        ds >> sym.stringrep >> static_cast<ExBase&>(sym);

        //get the ID number for this symbol
        sym.id = Symbol::getNewID(sym.stringrep);
    }
    else
        throw version_error(v, "1", r_symbol, CODELOC);

    return ds;
}

/** Null constructor */
Symbol::Symbol() : ExBase(), id(0), stringrep(QString::null)
{}

/** Construct a symbol from the passed ID number */
Symbol::Symbol(SymbolID symid)
       : ExBase(), id(symid), stringrep( Symbol::getName(symid) )
{}

/** Construct a new symbol, with string representation 'rep' */
Symbol::Symbol(const QString &rep)
       : ExBase(), id( Symbol::getNewID(rep) ), stringrep(rep)
{}

/** Copy constructor */
Symbol::Symbol(const Symbol &other)
       : ExBase(), id(other.id), stringrep(other.stringrep)
{}

/** Destructor */
Symbol::~Symbol()
{}

/** Return whether or not the symbol is null */
bool Symbol::isNull() const
{
    return id == 0;
}

/** Assignment operator */
Symbol& Symbol::operator=(const Symbol &other)
{
    id = other.id;
    stringrep = other.stringrep;
    return *this;
}

/** Assignment operator */
Symbol& Symbol::operator=(SymbolID symid)
{
    return this->operator=( Symbol(symid) );
}

/** Comparison operator */
bool Symbol::operator==(const ExBase &other) const
{
    const Symbol *sym = dynamic_cast<const Symbol*>(&other);

    return sym != 0 and sym->ID() == this->ID();
}

/** Comparison operator - a Symbol is greater than another
    symbol if it's string representation is greater - this
    allows lists of symbols to be sorted alphabetically */
bool Symbol::operator<(const Symbol &other) const
{
    return stringrep < other.stringrep;
}

/** Comparison operator - a Symbol is greater than another
    symbol if it's string representation is greater - this
    allows lists of symbols to be sorted alphabetically */
bool Symbol::operator>(const Symbol &other) const
{
    return stringrep > other.stringrep;
}

/** Comparison operator - a Symbol is greater than another
    symbol if it's string representation is greater - this
    allows lists of symbols to be sorted alphabetically */
bool Symbol::operator<=(const Symbol &other) const
{
    return stringrep <= other.stringrep;
}

/** Comparison operator - a Symbol is greater than another
    symbol if it's string representation is greater - this
    allows lists of symbols to be sorted alphabetically */
bool Symbol::operator>=(const Symbol &other) const
{
    return stringrep >= other.stringrep;
}

/** Return a string representation of this symbol */
QString Symbol::toString() const
{
    return stringrep;
}

/** There are no child expressions in a symbol */
Expressions Symbol::children() const
{
    return Expressions();
}

/** Evaluate this symbol - returns the value of the symbol in 'values' if
    it is present, else it returns 0.0 */
double Symbol::evaluate(const Values &values) const
{
    return values.value(*this);
}

/** Evaluate this symbol - returns the value of the symbol in 'values' if
    it is present, else it returns 0 */
Complex Symbol::evaluate(const ComplexValues &values) const
{
    return values.value(*this);
}

/** Differentiate this symbol with respect to 'sym'. This returns 1.0 if this
    is 'sym', else it returns 0.0 */
Expression Symbol::differentiate(const Symbol &sym) const
{
    //assume int(true) == 1 and int(false) == 0
    return Expression( int( sym.ID() == ID() ) );
}

/** Integrate this symbol with respect to 'sym'. If 'sym' == this, then
    return 0.5 sym^2, else return *this * sym */
Expression Symbol::integrate(const Symbol &sym) const
{
    if (sym.ID() == ID())
        return 0.5 * pow(sym,2);
    else
        return *this * sym;
}

/** Return the expression that matches this symbol in 'identities' - or return
    an expression holding only this symbol if it does no exist in 'identities' */
Expression Symbol::substitute(const Identities &identities) const
{
    return identities.expression(*this);
}

/** Return this symbol */
Symbols Symbol::symbols() const
{
    return Symbols(*this);
}

/** This is not a function */
Functions Symbol::functions() const
{
    return Functions();
}

/** Is this a function of 'symbol' */
bool Symbol::isFunction(const Symbol &sym) const
{
    return (id != 0) and (sym.ID() == ID());
}

/** A symbol is by definition not constant */
bool Symbol::isConstant() const
{
    return false;
}

QList<Factor> Symbol::expand(const Symbol &symbol) const
{
    QList<Factor> factors;

    if ( *this == symbol )
    {
        factors.append( Factor(symbol, 1,1) );
    }
    else
    {
        factors.append( Factor(symbol, *this,0) );
    }
    
    return factors;
}

const char* Symbol::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Symbol>() );
}

Symbol* Symbol::clone() const
{
    return new Symbol(*this);
}

