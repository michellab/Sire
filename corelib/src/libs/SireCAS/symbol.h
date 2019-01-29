/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIRECAS_SYMBOL_H
#define SIRECAS_SYMBOL_H

#include <QString>

#include "exbase.h"
#include "symbolvalue.h"
#include "symbolcomplex.h"
#include "symbolexpression.h"

SIRE_BEGIN_HEADER

namespace SireCAS
{
class Symbol;
}

SIRECAS_EXPORT QDataStream& operator<<(QDataStream&, const SireCAS::Symbol&);
SIRECAS_EXPORT QDataStream& operator>>(QDataStream&, SireCAS::Symbol&);

namespace SireCAS
{

class Symbols;
class Values;

/** This class represents an algebraic symbol in the equation (e.g. 'x' or 'y')

    @author Christopher Woods
*/
class SIRECAS_EXPORT Symbol : public ExBase
{

friend QDataStream& ::operator<<(QDataStream&, const Symbol&);
friend QDataStream& ::operator>>(QDataStream&, Symbol&);

public:
    Symbol();
    Symbol(SymbolID symid);
    Symbol(const QString &rep);

    Symbol(const Symbol &other);

    ~Symbol();

    Symbol& operator=(const Symbol &other);
    Symbol& operator=(SymbolID symid);

    /** Return the unique ID number of the symbol */
    SymbolID ID() const
    {
        return id;
    }

    /** Convienient operator used to combine a symbol with a value */
    SymbolValue operator==(double val) const
    {
        return SymbolValue(ID(), val);
    }

    /** Convienient operator used to combine a symbol with a value */
    SymbolValue operator==(int val) const
    {
        return SymbolValue(ID(), val);
    }

    /** Convienient operator used to combine a symbol with a Complex */
    SymbolComplex operator==(const Complex &val) const
    {
        return SymbolComplex(ID(), val);
    }

    /** Convieient operator used to combine a symbol with an equivalent
        expression */
    SymbolExpression operator==(const Expression &ex) const
    {
        return SymbolExpression(*this, ex);
    }

    Expression differentiate(const Symbol &symbol) const;
    Expression integrate(const Symbol &symbol) const;

    bool isFunction(const Symbol&) const;
    bool isConstant() const;

    bool operator==(const ExBase &other) const;

    bool operator<(const Symbol &other) const;
    bool operator>(const Symbol &other) const;
    bool operator<=(const Symbol &other) const;
    bool operator>=(const Symbol &other) const;

    uint hash() const;

    static const char* typeName();

    const char* what() const
    {
        return Symbol::typeName();
    }

    Symbol* clone() const;

    QString toString() const;

    bool isNull() const;

    double evaluate(const Values &values) const;
    Complex evaluate(const ComplexValues &values) const;

    Expression substitute(const Identities &identities) const;

    Symbols symbols() const;
    Functions functions() const;
    Expressions children() const;

    QList<Factor> expand(const Symbol &symbol) const;

protected:

    static SymbolID getNewID(const QString &symbol);
    static QString getName(SymbolID symid);

    /** Unique ID number that is given to every symbol */
    SymbolID id;

    /** String representation of this symbol */
    QString stringrep;
};

/** This is a small class that can hold the factor and power of a symbol

    @author Christopher Woods
*/
class SIRECAS_EXPORT Factor
{
public:
    Factor();
    
    Factor(const Symbol &symbol, double factor, double power);
    Factor(const Symbol &symbol,
           const Expression &factor, const Expression &power);
    
    Factor(const Factor &other);
    
    ~Factor();
    
    Factor& operator=(const Factor &other);
    
    bool operator==(const Factor &other) const;
    bool operator!=(const Factor &other) const;
    
    QString toString() const;
    
    const Symbol& symbol() const
    {
        return s;
    }
    
    const Expression& factor() const
    {
        return f;
    }
    
    const Expression& power() const
    {
        return p;
    }

private:
    /** The symbol for the factor */
    Symbol s;

    /** The factor and power */
    Expression f, p;
};

inline uint qHash(const Symbol &symbol)
{
    return symbol.hash();
}

}

Q_DECLARE_METATYPE(SireCAS::Symbol)

SIRE_EXPOSE_CLASS( SireCAS::Symbol )
SIRE_EXPOSE_CLASS( SireCAS::Factor )

SIRE_END_HEADER

#endif
