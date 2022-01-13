/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007   Christopher Woods
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

#include <Python.h>
#include <boost/python.hpp>

#include <QVector>
#include <QSet>

#include <boost/tuple/tuple.hpp>

#include "Helpers/convertlist.hpp"
#include "Helpers/convertdict.hpp"
#include "Helpers/convertset.hpp"
#include "Helpers/tuples.hpp"

#include "SireCAS/symbol.h"
#include "SireCAS/expression.h"
#include "SireCAS/symbolvalue.h"
#include "SireCAS/complexvalues.h"

using namespace SireCAS;

using boost::python::register_tuple;

void register_SireCAS_containers()
{
    register_list< QList<Expression> >();
    register_list< QVector<Expression> >();

    register_list< QList<Factor> >();

    register_list< QList<Symbol> >();
    register_list< QVector<Symbol> >();

    register_list< QList<SymbolValue> >();
    register_list< QList<SymbolComplex> >();
    register_list< QList<SymbolExpression> >();

    #if QT_VERSION >= QT_VERSION_CHECK(4, 2, 0)
    register_set< QSet<SymbolID> >();
    register_set< QSet<Symbol> >();

    register_set< Symbols >();

    register_dict< QHash<SymbolID,Complex> >();
    register_dict< QHash<SymbolID,Expression> >();
    register_dict< QHash<SymbolID,double> >();

    register_dict< QHash<Symbol,double> >();
    register_dict< QHash<Symbol,Complex> >();
    register_dict< QHash<Symbol,Expression> >();

    #else
    register_set< QSet<SymbolID>, SymbolID >();
    register_set< QSet<Symbol>, Symbol >();

    register_set< Symbols, Symbol >();

    register_dict< QHash<SymbolID,Complex>, SymbolID, Complex >();
    register_dict< QHash<SymbolID,Expression>, SymbolID, Expression >();
    register_dict< QHash<SymbolID,double>, SymbolID, double >();

    register_dict< QHash<Symbol,double>, Symbol, double >();
    register_dict< QHash<Symbol,Complex>, Symbol, Complex >();
    register_dict< QHash<Symbol,Expression>, Symbol, Expression >();

    #endif
}
