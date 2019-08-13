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

#ifndef SIRECAS_FUNCTIONSIGNATURE_H
#define SIRECAS_FUNCTIONSIGNATURE_H

#include <QString>
#include <QSet>
#include <QDataStream>

#include "symbolvalue.h"

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireCAS
{
class FunctionSignature;
}

QDataStream& operator<<(QDataStream&, const SireCAS::FunctionSignature&);
QDataStream& operator>>(QDataStream&, SireCAS::FunctionSignature&);

namespace SireCAS
{

class Function;

/** This small class holds a signature for a function. This is ID that will uniquely
    ID a type of function, but not its current differentiation level.
    
    \author Christopher Woods
*/
class SIRECAS_EXPORT FunctionSignature
{
public:
    FunctionSignature()
    {}
    
    FunctionSignature(const QString &name)
             : _name(name)
    {}
    
    FunctionSignature(const QString &name, const QSet<SymbolID> &args)
             : _name(name), _args(args)
    {}
    
    FunctionSignature(const FunctionSignature &other)
             : _name(other._name), _args(other._args)
    {}
    
    ~FunctionSignature()
    {}
    
    static const char* typeName();
    
    const char* what() const
    {
        return FunctionSignature::typeName();
    }
    
    void add(SymbolID id)
    {
        _args.insert(id);
    }
    
    void setName(const QString &name)
    {
        _name = name;
    }
    
    bool operator==(const FunctionSignature &other) const
    {
        return _name == other._name and _args == other._args;
    }
    
    bool operator!=(const FunctionSignature &other) const
    {
        return _name != other._name or _args != other._args;
    }

    bool contains(SymbolID id) const
    {
        return _args.contains(id);
    }

    QString name() const
    {
        return _name;
    }

    const QSet<SymbolID>& args() const
    {
        return _args;
    }

private:
    
    QString _name;
    QSet<SymbolID> _args;
};

/** Return a hash for this signature */
SIRE_ALWAYS_INLINE uint qHash(const FunctionSignature &sig)
{
    return qHash(sig.name());
}

}

Q_DECLARE_METATYPE(SireCAS::FunctionSignature)

SIRE_END_HEADER

#endif
