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

#ifndef SIRECAS_SYMBOLCOMPLEX_H
#define SIRECAS_SYMBOLCOMPLEX_H

#include "symbolvalue.h"

#include "SireMaths/complex.h"

SIRE_BEGIN_HEADER

namespace SireCAS
{

using SireMaths::Complex;

/** Small class that holds a SymbolID number and an associated complex value */
class SIRECAS_EXPORT SymbolComplex
{
public:
    SymbolComplex(SymbolID id, const Complex &val) : _val(val), _id(id)
    {}
    
    SymbolComplex(SymbolID id, double val) : _val(val), _id(id)
    {}
    
    SymbolComplex(const SymbolValue &val) : _val(val.value()), _id(val.ID())
    {}
    
    ~SymbolComplex()
    {}
    
    SymbolID ID() const
    {
        return _id;
    }
    
    const Complex& value() const
    {
        return _val;
    }
private:

    Complex _val;
    SymbolID _id;
};

}

SIRE_EXPOSE_CLASS( SireCAS::SymbolComplex )

SIRE_END_HEADER

#endif
