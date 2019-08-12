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

#ifndef SIRECAS_VALUES_H
#define SIRECAS_VALUES_H

#include <QHash>
#include <QList>

#include "symbolvalue.h"
#include "symbols.h"

SIRE_BEGIN_HEADER

namespace SireCAS
{
class Values;
}

class QDataStream;
SIRECAS_EXPORT QDataStream& operator<<(QDataStream&, const SireCAS::Values&);
SIRECAS_EXPORT QDataStream& operator>>(QDataStream&, SireCAS::Values&);

namespace SireCAS
{

class Symbol;

/**
    This class holds a set of Symbols and their associated values. This is used
    when numerically evaluating an equation.

    @author Christopher Woods
*/
class SIRECAS_EXPORT Values
{

friend SIRECAS_EXPORT QDataStream& ::operator<<(QDataStream&, const Values&);
friend SIRECAS_EXPORT QDataStream& ::operator>>(QDataStream&, Values&);

public:
    typedef QHash<SymbolID,double>::const_iterator const_iterator;

    Values();
    Values(const QList<SymbolValue> &values);
    Values(const QHash<Symbol,double> &values);
    Values(const SymbolValue &symval0);
    Values(const SymbolValue &symval0, const SymbolValue &symval1);
    Values(const SymbolValue &symval0, const SymbolValue &symval1, const SymbolValue &symval2);
    Values(const SymbolValue &symval0, const SymbolValue &symval1, const SymbolValue &symval2,
           const SymbolValue &symval3);
    Values(const SymbolValue &symval0, const SymbolValue &symval1, const SymbolValue &symval2,
           const SymbolValue &symval3, const SymbolValue &symval4);
    Values(const SymbolValue &symval0, const SymbolValue &symval1, const SymbolValue &symval2,
           const SymbolValue &symval3, const SymbolValue &symval4, const SymbolValue &symval5);
    Values(const SymbolValue &symval0, const SymbolValue &symval1, const SymbolValue &symval2,
           const SymbolValue &symval3, const SymbolValue &symval4, const SymbolValue &symval5,
           const SymbolValue &symval6);
    Values(const SymbolValue &symval0, const SymbolValue &symval1, const SymbolValue &symval2,
           const SymbolValue &symval3, const SymbolValue &symval4, const SymbolValue &symval5,
           const SymbolValue &symval6, const SymbolValue &symval7);
    Values(const SymbolValue &symval0, const SymbolValue &symval1, const SymbolValue &symval2,
           const SymbolValue &symval3, const SymbolValue &symval4, const SymbolValue &symval5,
           const SymbolValue &symval6, const SymbolValue &symval7, const SymbolValue &symval8);
    Values(const SymbolValue &symval0, const SymbolValue &symval1, const SymbolValue &symval2,
           const SymbolValue &symval3, const SymbolValue &symval4, const SymbolValue &symval5,
           const SymbolValue &symval6, const SymbolValue &symval7, const SymbolValue &symval8,
           const SymbolValue &symval9);

    Values(const Values &other);

    ~Values();

    static const char* typeName();
    
    const char* what() const
    {
        return Values::typeName();
    }

    bool operator==(const Values &other) const;
    bool operator!=(const Values &other) const;

    const_iterator begin() const;
    const_iterator end() const;

    const_iterator constBegin() const;
    const_iterator constEnd() const;
    
    const_iterator constFind(SymbolID symbolid) const;
    const_iterator constFind(const Symbol &symbol) const;

    void set(const Symbol &symbol, double value);

    void set(const const_iterator &it);

    void remove(const Symbol &symbol);
    void remove(const SymbolID &symbolid);

    void add(const SymbolValue &symval0);
    void add(const SymbolValue &symval0, const SymbolValue &symval1);
    void add(const SymbolValue &symval0, const SymbolValue &symval1, const SymbolValue &symval2);
    void add(const SymbolValue &symval0, const SymbolValue &symval1, const SymbolValue &symval2,
             const SymbolValue &symval3);
    void add(const SymbolValue &symval0, const SymbolValue &symval1, const SymbolValue &symval2,
             const SymbolValue &symval3, const SymbolValue &symval4);
    void add(const SymbolValue &symval0, const SymbolValue &symval1, const SymbolValue &symval2,
             const SymbolValue &symval3, const SymbolValue &symval4, const SymbolValue &symval5);
    void add(const SymbolValue &symval0, const SymbolValue &symval1, const SymbolValue &symval2,
             const SymbolValue &symval3, const SymbolValue &symval4, const SymbolValue &symval5,
             const SymbolValue &symval6);
    void add(const SymbolValue &symval0, const SymbolValue &symval1, const SymbolValue &symval2,
             const SymbolValue &symval3, const SymbolValue &symval4, const SymbolValue &symval5,
             const SymbolValue &symval6, const SymbolValue &symval7);
    void add(const SymbolValue &symval0, const SymbolValue &symval1, const SymbolValue &symval2,
             const SymbolValue &symval3, const SymbolValue &symval4, const SymbolValue &symval5,
             const SymbolValue &symval6, const SymbolValue &symval7, const SymbolValue &symval8);
    void add(const SymbolValue &symval0, const SymbolValue &symval1, const SymbolValue &symval2,
             const SymbolValue &symval3, const SymbolValue &symval4, const SymbolValue &symval5,
             const SymbolValue &symval6, const SymbolValue &symval7, const SymbolValue &symval8,
             const SymbolValue &symval9);

    double value(const Symbol &sym) const;

    double operator[](const Symbol &sym) const;
    double operator()(const Symbol &sym) const;

    const QHash<SymbolID, double>& values() const;

    QList<Symbol> keys() const;
    QList<Symbol> symbols() const;

    bool isEmpty() const;

    int count() const;

    bool contains(const Symbol &symbol) const;

    void reserve(int n);

    Values& operator+=(const SymbolValue &val);
    Values& operator+=(const Values &other);

    QString toString() const;

private:

    /** Hash mapping Symbol IDs to actual numerical values */
    QHash<SymbolID, double> vals;

};

SIRECAS_EXPORT Values operator+(const SymbolValue &val0, const SymbolValue &val1);

SIRECAS_EXPORT Values operator+(const Values &vals, const SymbolValue &val);
SIRECAS_EXPORT Values operator+(const SymbolValue &val, const Values &vals);

SIRECAS_EXPORT Values operator+(const Values &vals0, const Values &vals1);

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Reserve space for at least 'n' items */
SIRE_ALWAYS_INLINE void Values::reserve(int n)
{
    vals.reserve(n);
}

/** Return whether or not this set of values is empty */
SIRE_ALWAYS_INLINE bool Values::isEmpty() const
{
    return vals.isEmpty();
}

/** Return the number of specified values in this set */
SIRE_ALWAYS_INLINE int Values::count() const
{
    return vals.count();
}

/** Add a SymbolValue to the set of values */
SIRE_ALWAYS_INLINE void Values::add(const SymbolValue &val0)
{
    vals.insert(val0.ID(), val0.value());
}

/** Set the Symbol 'symbol' equal to 'value' */
SIRE_ALWAYS_INLINE void Values::set(const Symbol &symbol, double value)
{
    vals.insert(symbol.ID(), value);
}

/** Return the hash mapping the symbol ID to a value */
SIRE_ALWAYS_INLINE const QHash<SymbolID,double>& Values::values() const
{
    return vals;
}

/** Return whether or not a value for the symbol 'symbol' has been set */
SIRE_ALWAYS_INLINE bool Values::contains(const Symbol &symbol) const
{
    return vals.contains(symbol.ID());
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE(SireCAS::Values)

SIRE_EXPOSE_CLASS( SireCAS::Values )

SIRE_END_HEADER

#endif
