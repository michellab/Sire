/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#ifndef SIREMM_ATOMFUNCTIONS_H
#define SIREMM_ATOMFUNCTIONS_H

#include <QSet>

#include "SireBase/shareddatapointer.hpp"

#include "SireMol/cgatomidx.h"
#include "SireMol/molviewproperty.h"

#include "SireCAS/expression.h"
#include "SireCAS/symbols.h"
#include "SireCAS/identities.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class AtomFunction;
class AtomFunctions;
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::AtomFunction&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::AtomFunction&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::AtomFunctions&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::AtomFunctions&);

namespace SireMol
{
class MoleculeData;
}

namespace SireMM
{

using SireCAS::Symbol;
using SireCAS::Identities;
using SireCAS::Expression;

using SireMol::MoleculeData;
using SireMol::MoleculeInfoData;

/** This is the base class of all objects that hold the raw data
    for an AtomFunction (a function that acts between
    a set number of atoms)
    
    @author Christopher Woods
*/
class SIREMM_EXPORT AtomFunction
{

friend SIREMM_EXPORT QDataStream& ::operator<<(QDataStream&, const AtomFunction&);
friend SIREMM_EXPORT QDataStream& ::operator>>(QDataStream&, AtomFunction&);

public:
    ~AtomFunction();
    
    const SireCAS::Expression& function() const;
    
protected:
    AtomFunction();
    AtomFunction(const SireCAS::Expression &function);
    
    AtomFunction(const AtomFunction &other);
    
    AtomFunction& operator=(const AtomFunction &other);

    bool operator==(const AtomFunction &other) const;
    bool operator!=(const AtomFunction &other) const;
    
private:
    /** The function acting between the atoms */
    SireCAS::Expression func;
};

/** This is the base class of all of the AtomFunctions molecular
    properties (these are properties that contain lots of AtomFunction
    functions for the atoms in a molecule)
    
    @author Christopher Woods
*/
class SIREMM_EXPORT AtomFunctions : public SireMol::MoleculeProperty
{

friend SIREMM_EXPORT QDataStream& ::operator<<(QDataStream&, const AtomFunctions&);
friend SIREMM_EXPORT QDataStream& ::operator>>(QDataStream&, AtomFunctions&);

public:
    ~AtomFunctions();

    const QSet<Symbol>& symbols() const;

    bool isCompatibleWith(const MoleculeInfoData &molinfo) const;

protected:
    AtomFunctions();
    AtomFunctions(const MoleculeData &moldata);
    AtomFunctions(const MoleculeInfoData &molinfo);
    
    AtomFunctions(const AtomFunctions &other);
    
    AtomFunctions& operator=(const AtomFunctions &other);
    
    bool operator==(const AtomFunctions &other) const;
    bool operator!=(const AtomFunctions &other) const;

    const MoleculeInfoData& info() const
    {
        return *molinfo;
    }

    void addSymbols(const QSet<Symbol> &symbols);
    
    void removeSymbols(const QSet<Symbol> &symbols);
    void removeSymbols();

private:
    /** The metadata for the molecule that contains these functions */
    SireBase::SharedDataPointer<MoleculeInfoData> molinfo;

    /** The set of symbols used by the atom functions */
    QSet<Symbol> symbls;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

////////
//////// Inline functions of AtomFunction
////////

/** Return the function that acts between the atoms */
SIRE_ALWAYS_INLINE const Expression& AtomFunction::function() const
{
    return func;
}

////////
//////// Inline functions of AtomFunctions
////////

/** Return all of the symbols used by all of the AtomFunctions in this
    set */
SIRE_ALWAYS_INLINE const QSet<Symbol>& AtomFunctions::symbols() const
{
    return symbls;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

SIRE_EXPOSE_CLASS( SireMM::AtomFunction )
SIRE_EXPOSE_CLASS( SireMM::AtomFunctions )

SIRE_END_HEADER

#endif
