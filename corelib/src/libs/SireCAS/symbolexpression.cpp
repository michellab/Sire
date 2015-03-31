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

#include "symbolexpression.h"
#include "symbol.h"
#include "function.h"

#include <boost/assert.hpp>

using namespace SireCAS;

/** Return the symbol */
const Symbol& SymbolExpression::symbol() const
{
    BOOST_ASSERT( _sym.isA<Symbol>() );
    return _sym.asA<Symbol>();
}

/** Return whether or not this is a function */
bool SymbolExpression::isFunction() const
{
    return _sym.isA<Function>();
}

/** Return the symbol as a function - this is only safe
    if 'isFunction()' returns true */
const Function& SymbolExpression::function() const
{
    BOOST_ASSERT( _sym.isA<Function>() );
    return _sym.asA<Function>();
}
