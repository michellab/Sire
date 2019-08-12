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

#ifndef SIRECAS_SYMBOLS_H
#define SIRECAS_SYMBOLS_H

#include <QSet>

#include "symbol.h"

#include "tostring.h"

SIRE_BEGIN_HEADER

namespace SireCAS
{

/** Trival derivation of QSet<Symbol> that adds a constructor that
    automatically adds the passed Symbol 
    
    @author Christopher Woods
*/
class SIRECAS_EXPORT Symbols : public QSet<Symbol>
{
public:
    Symbols();
    
    Symbols(const Symbol &sym);
    
    Symbols(const QSet<Symbol> &other);
    
    Symbols(const QList<Symbol> &other);
    
    ~Symbols();
    
    Symbols& operator+=(const Symbols &other);
    Symbols& operator-=(const Symbols &other);
    
    Symbols operator+(const Symbols &other) const;
    Symbols operator-(const Symbols &other) const;
    
    Symbols& add(const Symbols &other);
    Symbols& subtract(const Symbols &other);
        
    QString toString() const;
    
    void insert(const Symbol &symbol);
    
    void insert(const Symbols &symbols);
    
    bool intersects(const QSet<Symbol> &other) const;
    
    bool contains(const Symbol &symbol) const;
    bool contains(const QSet<Symbol> &other) const;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

SIRE_ALWAYS_INLINE Symbols::Symbols() : QSet<Symbol>()
{}

SIRE_ALWAYS_INLINE Symbols::Symbols(const Symbol &sym) : QSet<Symbol>()
{
    this->insert(sym);
}

SIRE_ALWAYS_INLINE Symbols::Symbols(const QSet<Symbol> &other) : QSet<Symbol>(other)
{}

SIRE_ALWAYS_INLINE Symbols::Symbols(const QList<Symbol> &other) : QSet<Symbol>()
{
    int n = other.count();
    for (int i=0; i<n; ++i)
        insert( other.at(i) );
}

SIRE_ALWAYS_INLINE Symbols::~Symbols()
{}

SIRE_ALWAYS_INLINE Symbols& Symbols::add(const Symbols &other)
{
    QSet<Symbol>::unite(other);
    return *this;
}

SIRE_ALWAYS_INLINE Symbols& Symbols::subtract(const Symbols &other)
{
    QSet<Symbol>::subtract(other);
    return *this;
}

SIRE_ALWAYS_INLINE Symbols& Symbols::operator+=(const Symbols &other)
{
    return Symbols::add(other);
}

SIRE_ALWAYS_INLINE Symbols& Symbols::operator-=(const Symbols &other)
{
    return Symbols::subtract(other);
}

SIRE_ALWAYS_INLINE Symbols Symbols::operator+(const Symbols &other) const
{
    Symbols ret(*this);
    ret += other;
    
    return ret;
}

SIRE_ALWAYS_INLINE Symbols Symbols::operator-(const Symbols &other) const
{
    Symbols ret(*this);
    ret -= other;
    
    return ret;
}

SIRE_ALWAYS_INLINE QString Symbols::toString() const
{
    return Sire::toString( QSet<Symbol>(*this) );
}

SIRE_ALWAYS_INLINE void Symbols::insert(const Symbol &symbol)
{
    QSet<Symbol>::insert(symbol);
}

SIRE_ALWAYS_INLINE void Symbols::insert(const Symbols &symbols)
{
    for (Symbols::const_iterator it = symbols.begin();
         it != symbols.end();
         ++it)
    {
        this->insert(*it);
    }
}

SIRE_ALWAYS_INLINE bool Symbols::intersects(const QSet<Symbol> &other) const
{
    if (other.count() < this->count())
    {
        for (QSet<Symbol>::const_iterator it = other.constBegin();
             it != other.constEnd();
             ++it)
        {
            if (this->contains(*it))
                return true;
        }
    }
    else
    {
        for (QSet<Symbol>::const_iterator it = this->constBegin();
             it != this->constEnd();
             ++it)
        {
            if (other.contains(*it))
                return true;
        }
    }
    
    return false;
}

SIRE_ALWAYS_INLINE bool Symbols::contains(const Symbol &symbol) const
{
    return QSet<Symbol>::contains(symbol);
}

SIRE_ALWAYS_INLINE bool Symbols::contains(const QSet<Symbol> &other) const
{
    if (other.count() > this->count())
        return false;

    for (QSet<Symbol>::const_iterator it = other.constBegin();
         it != other.constEnd();
         ++it)
    {
        if (not this->contains(*it))
            return false;
    }
    
    return true;
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

SIRE_END_HEADER

#endif
