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

#include "complexvalues.h"
#include "symbol.h"

#include "SireStream/datastream.h"

using namespace SireStream;
using namespace SireCAS;

static const RegisterMetaType<ComplexValues> r_complexvals(NO_ROOT);

/** Serialise to a binary data stream */
QDataStream SIRECAS_EXPORT &operator<<(QDataStream &ds, const ComplexValues &vals)
{
    writeHeader(ds, r_complexvals, 1) << vals.vals;

    return ds;
}

/** Deserialise from a binary data stream */
QDataStream SIRECAS_EXPORT &operator>>(QDataStream &ds, ComplexValues &vals)
{
    VersionID v = readHeader(ds, r_complexvals);

    if (v == 1)
    {
        ds >> vals.vals;
    }
    else
        throw version_error(v, "1", r_complexvals, CODELOC);

    return ds;
}

/** Construct an empty set of values */
ComplexValues::ComplexValues()
{}

/** Copy constructor */
ComplexValues::ComplexValues(const ComplexValues &other) : vals(other.vals)
{}

/** Construct from a list of values */
ComplexValues::ComplexValues(const QList<SymbolComplex> &values)
{
    for (QList<SymbolComplex>::const_iterator it = values.begin();
         it != values.end();
         ++it)
    {
        add(*it);
    }
}

/** Construct from a hash of values indexed by symbol */
ComplexValues::ComplexValues(const QHash<Symbol,Complex> &values)
{
    for (QHash<Symbol,Complex>::const_iterator it = values.begin();
         it != values.end();
         ++it)
    {
        vals.insert(it.key().ID(), it.value());
    }
}

/** Construct from Values */
ComplexValues::ComplexValues(const Values &other)
{
    for (QHash<SymbolID,double>::const_iterator it = other.values().begin();
         it != other.values().end();
         ++it)
    {
        add( SymbolComplex(it.key(), it.value()) );
    }
}

/** Construct from the passed values */
ComplexValues::ComplexValues(const SymbolComplex &val0)
{
    add(val0);
}

/** Add the passed values */
void ComplexValues::add(const SymbolComplex &val0, const SymbolComplex &val1)
{
    add(val0);
    add(val1);
}

/** Add the passed values */
void ComplexValues::add(const SymbolComplex &val0, const SymbolComplex &val1,
                        const SymbolComplex &val2)
{
    add(val0);
    add(val1);
    add(val2);
}

/** Add the passed values */
void ComplexValues::add(const SymbolComplex &val0, const SymbolComplex &val1,
                        const SymbolComplex &val2, const SymbolComplex &val3)
{
    add(val0);
    add(val1);
    add(val2);
    add(val3);
}

/** Add the passed values */
void ComplexValues::add(const SymbolComplex &val0, const SymbolComplex &val1,
                        const SymbolComplex &val2, const SymbolComplex &val3,
                        const SymbolComplex &val4)
{
    add(val0);
    add(val1);
    add(val2);
    add(val3);
    add(val4);
}

/** Add the passed values */
void ComplexValues::add(const SymbolComplex &val0, const SymbolComplex &val1,
                        const SymbolComplex &val2, const SymbolComplex &val3,
                        const SymbolComplex &val4, const SymbolComplex &val5)
{
    add(val0);
    add(val1);
    add(val2);
    add(val3);
    add(val4);
    add(val5);
}

/** Add the passed values */
void ComplexValues::add(const SymbolComplex &val0, const SymbolComplex &val1,
                        const SymbolComplex &val2, const SymbolComplex &val3,
                        const SymbolComplex &val4, const SymbolComplex &val5,
                        const SymbolComplex &val6)
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
void ComplexValues::add(const SymbolComplex &val0, const SymbolComplex &val1,
                        const SymbolComplex &val2, const SymbolComplex &val3,
                        const SymbolComplex &val4, const SymbolComplex &val5,
                        const SymbolComplex &val6, const SymbolComplex &val7)
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
void ComplexValues::add(const SymbolComplex &val0, const SymbolComplex &val1,
                        const SymbolComplex &val2, const SymbolComplex &val3,
                        const SymbolComplex &val4, const SymbolComplex &val5,
                        const SymbolComplex &val6, const SymbolComplex &val7,
                        const SymbolComplex &val8)
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
void ComplexValues::add(const SymbolComplex &val0, const SymbolComplex &val1,
                        const SymbolComplex &val2, const SymbolComplex &val3,
                        const SymbolComplex &val4, const SymbolComplex &val5,
                        const SymbolComplex &val6, const SymbolComplex &val7,
                        const SymbolComplex &val8, const SymbolComplex &val9)
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
ComplexValues::ComplexValues(const SymbolComplex &val0, const SymbolComplex &val1)
{
    add(val0,val1);
}

/** Construct from the passed values */
ComplexValues::ComplexValues(const SymbolComplex &val0, const SymbolComplex &val1,
                             const SymbolComplex &val2)
{
    add(val0,val1,val2);
}

/** Construct from the passed values */
ComplexValues::ComplexValues(const SymbolComplex &val0, const SymbolComplex &val1,
                             const SymbolComplex &val2, const SymbolComplex &val3)
{
    add(val0,val1,val2,val3);
}

/** Construct from the passed values */
ComplexValues::ComplexValues(const SymbolComplex &val0, const SymbolComplex &val1,
                             const SymbolComplex &val2, const SymbolComplex &val3,
                             const SymbolComplex &val4)
{
    add(val0,val1,val2,val3,val4);
}

/** Construct from the passed values */
ComplexValues::ComplexValues(const SymbolComplex &val0, const SymbolComplex &val1,
                             const SymbolComplex &val2, const SymbolComplex &val3,
                             const SymbolComplex &val4, const SymbolComplex &val5)
{
    add(val0,val1,val2,val3,val4,val5);
}

/** Construct from the passed values */
ComplexValues::ComplexValues(const SymbolComplex &val0, const SymbolComplex &val1,
                             const SymbolComplex &val2, const SymbolComplex &val3,
                             const SymbolComplex &val4, const SymbolComplex &val5,
                             const SymbolComplex &val6)
{
    add(val0,val1,val2,val3,val4,val5,val6);
}

/** Construct from the passed values */
ComplexValues::ComplexValues(const SymbolComplex &val0, const SymbolComplex &val1,
                             const SymbolComplex &val2, const SymbolComplex &val3,
                             const SymbolComplex &val4, const SymbolComplex &val5,
                             const SymbolComplex &val6, const SymbolComplex &val7)
{
    add(val0,val1,val2,val3,val4,val5,val6,val7);
}

/** Construct from the passed values */
ComplexValues::ComplexValues(const SymbolComplex &val0, const SymbolComplex &val1,
                             const SymbolComplex &val2, const SymbolComplex &val3,
                             const SymbolComplex &val4, const SymbolComplex &val5,
                             const SymbolComplex &val6, const SymbolComplex &val7,
                             const SymbolComplex &val8)
{
    add(val0,val1,val2,val3,val4,val5,val6,val7,val8);
}

/** Construct from the passed values */
ComplexValues::ComplexValues(const SymbolComplex &val0, const SymbolComplex &val1,
                             const SymbolComplex &val2, const SymbolComplex &val3,
                             const SymbolComplex &val4, const SymbolComplex &val5,
                             const SymbolComplex &val6, const SymbolComplex &val7,
                             const SymbolComplex &val8, const SymbolComplex &val9)
{
    add(val0,val1,val2,val3,val4,val5,val6,val7,val8,val9);
}

/** Destructor */
ComplexValues::~ComplexValues()
{}

/** Return the value of the Symbol with ID 'id', or 0.0 if there is no such symbol */
Complex ComplexValues::value(const Symbol &sym) const
{
    return vals.value(sym.ID(),Complex(0));
}

const char* ComplexValues::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ComplexValues>() );
}
