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

#ifndef SIRECAS_SYMBOLEXPRESSION_H
#define SIRECAS_SYMBOLEXPRESSION_H

#include <qglobal.h>

#include "symbolvalue.h"
#include "symbolcomplex.h"
#include "expression.h"

SIRE_BEGIN_HEADER

namespace SireCAS
{

class Symbol;
class Function;

/** Small class that holds a SymbolID number and an associated expression */
class SIRECAS_EXPORT SymbolExpression
{
public:
    SymbolExpression(const ExpressionBase &symbol, const Expression &expression) 
                : _ex(expression), _sym(symbol)
    {}
    
    ~SymbolExpression()
    {}
    
    const Symbol& symbol() const;
    
    bool isFunction() const;
    const Function& function() const;
    
    const Expression& expression() const
    {
        return _ex;
    }
private:

    Expression _ex;
    ExpressionBase _sym;
};

}

SIRE_EXPOSE_CLASS( SireCAS::SymbolExpression )

SIRE_END_HEADER

#endif
