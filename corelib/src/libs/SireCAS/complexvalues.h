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

#ifndef SIRECAS_COMPLEXVALUES_H
#define SIRECAS_COMPLEXVALUES_H

#include <QHash>

#include "symbolcomplex.h"
#include "values.h"
#include "symbols.h"

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireCAS
{
class ComplexValues;
}

class QDataStream;
SIRECAS_EXPORT QDataStream& operator<<(QDataStream&, const SireCAS::ComplexValues&);
SIRECAS_EXPORT QDataStream& operator>>(QDataStream&, SireCAS::ComplexValues&);

namespace SireCAS
{

class Symbol;

/**
This class holds a set of Symbols and their associated complex values. 
This is used when numerically evaluating an equation using complex maths.

@author Christopher Woods
*/
class SIRECAS_EXPORT ComplexValues
{

friend SIRECAS_EXPORT QDataStream& ::operator<<(QDataStream&, const ComplexValues&);
friend SIRECAS_EXPORT QDataStream& ::operator>>(QDataStream&, ComplexValues&);

public:
    ComplexValues();

    ComplexValues(const QList<SymbolComplex> &values);
    ComplexValues(const QHash<Symbol,Complex> &values);

    ComplexValues(const SymbolComplex &symval0);
    ComplexValues(const SymbolComplex &symval0, const SymbolComplex &symval1);
    ComplexValues(const SymbolComplex &symval0, const SymbolComplex &symval1,
                  const SymbolComplex &symval2);
    ComplexValues(const SymbolComplex &symval0, const SymbolComplex &symval1,
                  const SymbolComplex &symval2, const SymbolComplex &symval3);
    ComplexValues(const SymbolComplex &symval0, const SymbolComplex &symval1,
                  const SymbolComplex &symval2, const SymbolComplex &symval3,
                  const SymbolComplex &symval4);
    ComplexValues(const SymbolComplex &symval0, const SymbolComplex &symval1,
                  const SymbolComplex &symval2, const SymbolComplex &symval3,
                  const SymbolComplex &symval4, const SymbolComplex &symval5);
    ComplexValues(const SymbolComplex &symval0, const SymbolComplex &symval1,
                  const SymbolComplex &symval2, const SymbolComplex &symval3,
                  const SymbolComplex &symval4, const SymbolComplex &symval5,
                  const SymbolComplex &symval6);
    ComplexValues(const SymbolComplex &symval0, const SymbolComplex &symval1,
                  const SymbolComplex &symval2, const SymbolComplex &symval3,
                  const SymbolComplex &symval4, const SymbolComplex &symval5,
                  const SymbolComplex &symval6, const SymbolComplex &symval7);
    ComplexValues(const SymbolComplex &symval0, const SymbolComplex &symval1,
                  const SymbolComplex &symval2, const SymbolComplex &symval3,
                  const SymbolComplex &symval4, const SymbolComplex &symval5,
                  const SymbolComplex &symval6, const SymbolComplex &symval7,
                  const SymbolComplex &symval8);
    ComplexValues(const SymbolComplex &symval0, const SymbolComplex &symval1,
                  const SymbolComplex &symval2, const SymbolComplex &symval3,
                  const SymbolComplex &symval4, const SymbolComplex &symval5,
                  const SymbolComplex &symval6, const SymbolComplex &symval7,
                  const SymbolComplex &symval8, const SymbolComplex &symval9);

    ComplexValues(const Values &other);

    ComplexValues(const ComplexValues &other);

    ~ComplexValues();

    static const char* typeName();

    const char* what() const
    {
        return ComplexValues::typeName();
    }

    void set(const Symbol &symbol, const Complex &value);

    void add(const SymbolComplex &symval0);
    void add(const SymbolComplex &symval0, const SymbolComplex &symval1);
    void add(const SymbolComplex &symval0, const SymbolComplex &symval1,
             const SymbolComplex &symval2);
    void add(const SymbolComplex &symval0, const SymbolComplex &symval1,
             const SymbolComplex &symval2, const SymbolComplex &symval3);
    void add(const SymbolComplex &symval0, const SymbolComplex &symval1,
             const SymbolComplex &symval2, const SymbolComplex &symval3,
             const SymbolComplex &symval4);
    void add(const SymbolComplex &symval0, const SymbolComplex &symval1,
             const SymbolComplex &symval2, const SymbolComplex &symval3,
             const SymbolComplex &symval4, const SymbolComplex &symval5);
    void add(const SymbolComplex &symval0, const SymbolComplex &symval1,
             const SymbolComplex &symval2, const SymbolComplex &symval3,
             const SymbolComplex &symval4, const SymbolComplex &symval5,
             const SymbolComplex &symval6);
    void add(const SymbolComplex &symval0, const SymbolComplex &symval1,
             const SymbolComplex &symval2, const SymbolComplex &symval3,
             const SymbolComplex &symval4, const SymbolComplex &symval5,
             const SymbolComplex &symval6, const SymbolComplex &symval7);
    void add(const SymbolComplex &symval0, const SymbolComplex &symval1,
             const SymbolComplex &symval2, const SymbolComplex &symval3,
             const SymbolComplex &symval4, const SymbolComplex &symval5,
             const SymbolComplex &symval6, const SymbolComplex &symval7,
             const SymbolComplex &symval8);
    void add(const SymbolComplex &symval0, const SymbolComplex &symval1,
             const SymbolComplex &symval2, const SymbolComplex &symval3,
             const SymbolComplex &symval4, const SymbolComplex &symval5,
             const SymbolComplex &symval6, const SymbolComplex &symval7,
             const SymbolComplex &symval8, const SymbolComplex &symval9);

    Complex value(const Symbol &sym) const;

    const QHash<SymbolID,Complex>& values() const;

private:

    /** Hash mapping Symbol IDs to actual numerical values */
    QHash<SymbolID, Complex> vals;

};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Add a SymbolComplex to the set of values */
SIRE_ALWAYS_INLINE void ComplexValues::add(const SymbolComplex &val0)
{
    vals.insert(val0.ID(), val0.value());
}

/** Set the Symbol 'symbol' equal to 'value' */
SIRE_ALWAYS_INLINE void ComplexValues::set(const Symbol &symbol, const Complex &value)
{
    vals.insert(symbol.ID(), value);
}

/** Return the hash mapping Symbol IDs to complex values */
SIRE_ALWAYS_INLINE const QHash<SymbolID,Complex>& ComplexValues::values() const
{
    return vals;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE(SireCAS::ComplexValues)

SIRE_EXPOSE_CLASS( SireCAS::ComplexValues )

SIRE_END_HEADER

#endif
