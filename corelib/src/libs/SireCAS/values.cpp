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

#include "values.h"
#include "symbol.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireStream;
using namespace SireCAS;

static const RegisterMetaType<Values> r_values(NO_ROOT);

/** Serialise to a binary data stream */
QDataStream &operator<<(QDataStream &ds, const Values &values)
{
    writeHeader(ds, r_values, 1);

    SharedDataStream sds(ds);

    sds << quint32( values.vals.count() );

    for (QHash<SymbolID,double>::const_iterator it = values.vals.constBegin();
         it != values.vals.constEnd();
         ++it)
    {
        sds << Symbol(it.key()) << it.value();
    }

    return ds;
}

/** Deserialise from a binary data stream */
QDataStream &operator>>(QDataStream &ds, Values &values)
{
    VersionID v = readHeader(ds, r_values);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        quint32 nvals;
        sds >> nvals;

        QHash<SymbolID,double> vals;
        vals.reserve(nvals);

        for (quint32 i=0; i<nvals; ++i)
        {
            Symbol symbol;
            double value;

            sds >> symbol >> value;

            vals.insert( symbol.ID(), value );
        }

        values.vals = vals;
    }
    else
        throw version_error(v, "1", r_values, CODELOC);

    return ds;
}

/** Construct an empty set of values */
Values::Values()
{}

/** Construct from a list of values */
Values::Values(const QList<SymbolValue> &values)
{
    for (QList<SymbolValue>::const_iterator it = values.begin();
         it != values.end();
         ++it)
    {
        add(*it);
    }
}

/** Construct from a hash of values indexed by symbols */
Values::Values(const QHash<Symbol,double> &values)
{
    for (QHash<Symbol,double>::const_iterator it = values.begin();
         it != values.end();
         ++it)
    {
        vals.insert( it.key().ID(), it.value() );
    }
}

/** Copy constructor */
Values::Values(const Values &other) : vals(other.vals)
{}

/** Comparison operator */
bool Values::operator==(const Values &other) const
{
    return vals == other.vals;
}

/** Comparison operator */
bool Values::operator!=(const Values &other) const
{
    return vals != other.vals;
}

/** Return a string representation of these values */
QString Values::toString() const
{
    QStringList words;
    QStringList lines;

    QList<Symbol> syms = this->symbols();

    std::sort(syms.begin(), syms.end());

    foreach (const Symbol &sym, syms)
    {
        words.append( QString("%1 == %2").arg(sym.toString())
                                         .arg(this->value(sym)) );

        if (words.count() == 4)
        {
            lines.append( words.join(", ") );
            words.clear();
        }
    }

    if (not words.isEmpty())
    {
        lines.append( words.join(", ") );
    }

    return QString("{ %1 }").arg( lines.join("\n  ") );
}

/** Return a list of the symbols that are present in this set */
QList<Symbol> Values::symbols() const
{
    QList<Symbol> s;

    for (QHash<SymbolID,double>::const_iterator it = vals.constBegin();
         it != vals.constEnd();
         ++it)
    {
        s.append( Symbol(it.key()) );
    }

    return s;
}

/** Return a list of the symbols that are present in this set */
QList<Symbol> Values::keys() const
{
    return this->symbols();
}

/** Construct from the passed values */
Values::Values(const SymbolValue &val0)
{
    add(val0);
}

/** Add the passed values */
void Values::add(const SymbolValue &val0, const SymbolValue &val1)
{
    add(val0);
    add(val1);
}

/** Add the passed values */
void Values::add(const SymbolValue &val0, const SymbolValue &val1, const SymbolValue &val2)
{
    add(val0);
    add(val1);
    add(val2);
}

/** Add the passed values */
void Values::add(const SymbolValue &val0, const SymbolValue &val1, const SymbolValue &val2,
                 const SymbolValue &val3)
{
    add(val0);
    add(val1);
    add(val2);
    add(val3);
}

/** Add the passed values */
void Values::add(const SymbolValue &val0, const SymbolValue &val1, const SymbolValue &val2,
                 const SymbolValue &val3, const SymbolValue &val4)
{
    add(val0);
    add(val1);
    add(val2);
    add(val3);
    add(val4);
}

/** Add the passed values */
void Values::add(const SymbolValue &val0, const SymbolValue &val1, const SymbolValue &val2,
                 const SymbolValue &val3, const SymbolValue &val4, const SymbolValue &val5)
{
    add(val0);
    add(val1);
    add(val2);
    add(val3);
    add(val4);
    add(val5);
}

/** Add the passed values */
void Values::add(const SymbolValue &val0, const SymbolValue &val1, const SymbolValue &val2,
                 const SymbolValue &val3, const SymbolValue &val4, const SymbolValue &val5,
                 const SymbolValue &val6)
{
    add(val0);
    add(val1);
    add(val2);
    add(val3);
    add(val4);
    add(val5);
    add(val6);
}

/** Add the passed values */
void Values::add(const SymbolValue &val0, const SymbolValue &val1, const SymbolValue &val2,
                 const SymbolValue &val3, const SymbolValue &val4, const SymbolValue &val5,
                 const SymbolValue &val6, const SymbolValue &val7)
{
    add(val0);
    add(val1);
    add(val2);
    add(val3);
    add(val4);
    add(val5);
    add(val6);
    add(val7);
}

/** Add the passed values */
void Values::add(const SymbolValue &val0, const SymbolValue &val1, const SymbolValue &val2,
                 const SymbolValue &val3, const SymbolValue &val4, const SymbolValue &val5,
                 const SymbolValue &val6, const SymbolValue &val7, const SymbolValue &val8)
{
    add(val0);
    add(val1);
    add(val2);
    add(val3);
    add(val4);
    add(val5);
    add(val6);
    add(val7);
    add(val8);
}

/** Add the passed values */
void Values::add(const SymbolValue &val0, const SymbolValue &val1, const SymbolValue &val2,
                 const SymbolValue &val3, const SymbolValue &val4, const SymbolValue &val5,
                 const SymbolValue &val6, const SymbolValue &val7, const SymbolValue &val8,
                 const SymbolValue &val9)
{
    add(val0);
    add(val1);
    add(val2);
    add(val3);
    add(val4);
    add(val5);
    add(val6);
    add(val7);
    add(val8);
    add(val9);
}

/** Construct from the passed values */
Values::Values(const SymbolValue &val0, const SymbolValue &val1)
{
    add(val0,val1);
}

/** Construct from the passed values */
Values::Values(const SymbolValue &val0, const SymbolValue &val1, const SymbolValue &val2)
{
    add(val0,val1,val2);
}

/** Construct from the passed values */
Values::Values(const SymbolValue &val0, const SymbolValue &val1, const SymbolValue &val2,
               const SymbolValue &val3)
{
    add(val0,val1,val2,val3);
}

/** Construct from the passed values */
Values::Values(const SymbolValue &val0, const SymbolValue &val1, const SymbolValue &val2,
               const SymbolValue &val3, const SymbolValue &val4)
{
    add(val0,val1,val2,val3,val4);
}

/** Construct from the passed values */
Values::Values(const SymbolValue &val0, const SymbolValue &val1, const SymbolValue &val2,
               const SymbolValue &val3, const SymbolValue &val4, const SymbolValue &val5)
{
    add(val0,val1,val2,val3,val4,val5);
}

/** Construct from the passed values */
Values::Values(const SymbolValue &val0, const SymbolValue &val1, const SymbolValue &val2,
               const SymbolValue &val3, const SymbolValue &val4, const SymbolValue &val5,
               const SymbolValue &val6)
{
    add(val0,val1,val2,val3,val4,val5,val6);
}

/** Construct from the passed values */
Values::Values(const SymbolValue &val0, const SymbolValue &val1, const SymbolValue &val2,
               const SymbolValue &val3, const SymbolValue &val4, const SymbolValue &val5,
               const SymbolValue &val6, const SymbolValue &val7)
{
    add(val0,val1,val2,val3,val4,val5,val6,val7);
}

/** Construct from the passed values */
Values::Values(const SymbolValue &val0, const SymbolValue &val1, const SymbolValue &val2,
               const SymbolValue &val3, const SymbolValue &val4, const SymbolValue &val5,
               const SymbolValue &val6, const SymbolValue &val7, const SymbolValue &val8)
{
    add(val0,val1,val2,val3,val4,val5,val6,val7,val8);
}

/** Construct from the passed values */
Values::Values(const SymbolValue &val0, const SymbolValue &val1, const SymbolValue &val2,
               const SymbolValue &val3, const SymbolValue &val4, const SymbolValue &val5,
               const SymbolValue &val6, const SymbolValue &val7, const SymbolValue &val8,
               const SymbolValue &val9)
{
    add(val0,val1,val2,val3,val4,val5,val6,val7,val8,val9);
}

/** Destructor */
Values::~Values()
{}

/** Return the value of the Symbol with ID 'id', or 0.0 if there is no such symbol */
double Values::value(const Symbol &sym) const
{
    return vals.value(sym.ID(),0.0);
}

/** Return the value of the Symbol with ID 'id', or 0.0 if there is no such symbol */
double Values::operator[](const Symbol &sym) const
{
    return this->value(sym);
}

/** Return the value of the Symbol with ID 'id', or 0.0 if there is no such symbol */
double Values::operator()(const Symbol &sym) const
{
    return this->value(sym);
}

/** Add the value 'val' to this set */
Values& Values::operator+=(const SymbolValue &val)
{
    this->add(val);
    return *this;
}

/** Add the contents of 'other' to this set - this overwrites any
    existing values that are also in 'other' */
Values& Values::operator+=(const Values &other)
{
    if (other.vals.isEmpty())
        return *this;
    else if (vals.isEmpty())
    {
        vals = other.vals;
        return *this;
    }
    else
    {
        vals.reserve( vals.count() + other.vals.count() );

        for (QHash<SymbolID,double>::const_iterator it = other.vals.begin();
             it != other.vals.end();
             ++it)
        {
            vals.insert( it.key(), it.value() );
        }

        return *this;
    }
}

/** Return an iterator to the first symbol/value pair */
Values::const_iterator Values::begin() const
{
    return vals.constBegin();
}

/** Return an iterator pointing to one past the last
    symbol/value pair */
Values::const_iterator Values::end() const
{
    return vals.constEnd();
}

/** Return an iterator to the first symbol/value pair */
Values::const_iterator Values::constBegin() const
{
    return vals.constBegin();
}

/** Return an iterator pointing to one past the last
    symbol/value pair */
Values::const_iterator Values::constEnd() const
{
    return vals.constEnd();
}

/** Return an iterator pointing to the symbol/value pair
    with symbol with ID 'symbolid' - or Values::end() if
    there is no matching symbol */
Values::const_iterator Values::constFind(SymbolID symbolid) const
{
    return vals.constFind(symbolid);
}

/** Return an iterator pointing to the symbol/value pair
    with symbol 'symbol' - or Values::end() if
    there is no matching symbol */
Values::const_iterator Values::constFind(const Symbol &symbol) const
{
    return vals.constFind(symbol.ID());
}

/** Set the value of the symbol/value pair pointed to by the
    iterator 'it' in this set */
void Values::set(const Values::const_iterator &it)
{
    if (it != vals.constEnd())
        vals.insert( it.key(), it.value() );
}

/** Remove the value for the symbol 'symbol' */
void Values::remove(const Symbol &symbol)
{
    vals.remove(symbol.ID());
}

/** Remove the value for the symbol with ID 'symbolid' */
void Values::remove(const SymbolID &symbolid)
{
    vals.remove(symbolid);
}

const char* Values::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Values>() );
}

namespace SireCAS
{
    Values operator+(const SymbolValue &val0, const SymbolValue &val1)
    {
        Values vals(val0);
        vals += val1;

        return vals;
    }

    Values operator+(const Values &vals, const SymbolValue &val)
    {
        Values new_vals(vals);
        new_vals += val;
        return new_vals;
    }

    Values operator+(const SymbolValue &val, const Values &vals)
    {
        Values new_vals(vals);
        new_vals += val;
        return new_vals;
    }

    Values operator+(const Values &vals0, const Values &vals1)
    {
        Values new_vals(vals0);
        new_vals += vals1;
        return new_vals;
    }
}
