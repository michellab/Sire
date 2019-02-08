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

#ifndef SIRECAS_IDENTITIES_H
#define SIRECAS_IDENTITIES_H

#include <QHash>
#include <QSet>
#include <QList>

#include "symbolexpression.h"
#include "functionsignature.h"

SIRE_BEGIN_HEADER

namespace SireCAS
{
class Identities;
}

class QDataStream;
SIRECAS_EXPORT QDataStream& operator<<(QDataStream&, const SireCAS::Identities&);
SIRECAS_EXPORT QDataStream& operator>>(QDataStream&, SireCAS::Identities&);

namespace SireCAS
{

/** This class holds a collection of identities. 
    These are symbols, and the expressions that they are equal to, 
    e.g. x == y.squared() + sin(z). These identities can then be 
    substituted into expressions, with all instances of the variable 
    being replaced by the equivalent expression.

    @author Christopher Woods
*/
class SIRECAS_EXPORT Identities
{

friend SIRECAS_EXPORT QDataStream& ::operator<<(QDataStream&, const Identities&);
friend SIRECAS_EXPORT QDataStream& ::operator>>(QDataStream&, Identities&);

public:
    Identities();
    Identities(const QList<SymbolExpression> &expressions);
    Identities(const QHash<Symbol,Expression> &expressions);
    Identities(const SymbolExpression &symex0);
    Identities(const SymbolExpression &symex0, const SymbolExpression &symex1);
    Identities(const SymbolExpression &symex0, const SymbolExpression &symex1,
               const SymbolExpression &symex2);
    Identities(const SymbolExpression &symex0, const SymbolExpression &symex1,
               const SymbolExpression &symex2, const SymbolExpression &symex3);
    Identities(const SymbolExpression &symex0, const SymbolExpression &symex1,
               const SymbolExpression &symex2, const SymbolExpression &symex3,
               const SymbolExpression &symex4);
    Identities(const SymbolExpression &symex0, const SymbolExpression &symex1,
               const SymbolExpression &symex2, const SymbolExpression &symex3,
               const SymbolExpression &symex4, const SymbolExpression &symex5);
    Identities(const SymbolExpression &symex0, const SymbolExpression &symex1,
               const SymbolExpression &symex2, const SymbolExpression &symex3,
               const SymbolExpression &symex4, const SymbolExpression &symex5,
               const SymbolExpression &symex6);
    Identities(const SymbolExpression &symex0, const SymbolExpression &symex1,
               const SymbolExpression &symex2, const SymbolExpression &symex3,
               const SymbolExpression &symex4, const SymbolExpression &symex5,
               const SymbolExpression &symex6, const SymbolExpression &symex7);
    Identities(const SymbolExpression &symex0, const SymbolExpression &symex1,
               const SymbolExpression &symex2, const SymbolExpression &symex3,
               const SymbolExpression &symex4, const SymbolExpression &symex5,
               const SymbolExpression &symex6, const SymbolExpression &symex7,
               const SymbolExpression &symex8);
    Identities(const SymbolExpression &symex0, const SymbolExpression &symex1,
               const SymbolExpression &symex2, const SymbolExpression &symex3,
               const SymbolExpression &symex4, const SymbolExpression &symex5,
               const SymbolExpression &symex6, const SymbolExpression &symex7,
               const SymbolExpression &symex8, const SymbolExpression &symex9);

    Identities(const Identities &other);

    ~Identities();

    static const char* typeName();
    
    const char* what() const
    {
        return Identities::typeName();
    }

    bool operator==(const Identities &other) const;
    bool operator!=(const Identities &other) const;
    
    void add(const SymbolExpression &symex0);
    void add(const SymbolExpression &symex0, const SymbolExpression &symex1);
    void add(const SymbolExpression &symex0, const SymbolExpression &symex1,
             const SymbolExpression &symex2);
    void add(const SymbolExpression &symex0, const SymbolExpression &symex1,
             const SymbolExpression &symex2, const SymbolExpression &symex3);
    void add(const SymbolExpression &symex0, const SymbolExpression &symex1,
             const SymbolExpression &symex2, const SymbolExpression &symex3,
             const SymbolExpression &symex4);
    void add(const SymbolExpression &symex0, const SymbolExpression &symex1,
             const SymbolExpression &symex2, const SymbolExpression &symex3,
             const SymbolExpression &symex4, const SymbolExpression &symex5);
    void add(const SymbolExpression &symex0, const SymbolExpression &symex1,
             const SymbolExpression &symex2, const SymbolExpression &symex3,
             const SymbolExpression &symex4, const SymbolExpression &symex5,
             const SymbolExpression &symex6);
    void add(const SymbolExpression &symex0, const SymbolExpression &symex1,
             const SymbolExpression &symex2, const SymbolExpression &symex3,
             const SymbolExpression &symex4, const SymbolExpression &symex5,
             const SymbolExpression &symex6, const SymbolExpression &symex7);
    void add(const SymbolExpression &symex0, const SymbolExpression &symex1,
             const SymbolExpression &symex2, const SymbolExpression &symex3,
             const SymbolExpression &symex4, const SymbolExpression &symex5,
             const SymbolExpression &symex6, const SymbolExpression &symex7,
             const SymbolExpression &symex8);
    void add(const SymbolExpression &symex0, const SymbolExpression &symex1,
             const SymbolExpression &symex2, const SymbolExpression &symex3,
             const SymbolExpression &symex4, const SymbolExpression &symex5,
             const SymbolExpression &symex6, const SymbolExpression &symex7,
             const SymbolExpression &symex8, const SymbolExpression &symex9);

    bool contains(const Symbol &sym) const;
    Expression expression(const Symbol &sym) const;

    Expression operator[](const Symbol &sym) const;

    bool contains(const Function &func) const;
    Function function(const Function &func) const;

    QString toString() const;

    QList<Symbol> symbols() const;

    void set( const Symbol &symbol, const Expression &expression );

private:

    /** Hash mapping a symbol to an expression */
    QHash<SymbolID, Expression> idhash;

    /** Hash mapping the signatures of functions to the actual function
        stored in this collection of identities */
    QHash<FunctionSignature, Expression> funchash;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Add the SymbolExpression to this set of identities */
inline void Identities::add(const SymbolExpression &symex0)
{
    set( symex0.symbol(), symex0.expression() );
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE(SireCAS::Identities)

SIRE_EXPOSE_CLASS( SireCAS::Identities )

SIRE_END_HEADER

#endif
