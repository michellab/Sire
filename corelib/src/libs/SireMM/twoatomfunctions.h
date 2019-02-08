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

#ifndef SIREMM_TWOATOMFUNCTIONS_H
#define SIREMM_TWOATOMFUNCTIONS_H

#include <QHash>

#include "atomfunctions.h"

#include "SireMol/cgatomidx.h"
#include "SireMol/atomidx.h"
#include "SireMol/atomid.h"

#include "SireMol/bondid.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class TwoAtomFunction;
class TwoAtomFunctions;
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::TwoAtomFunction&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::TwoAtomFunction&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::TwoAtomFunctions&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::TwoAtomFunctions&);

namespace SireMol
{
class AtomSelection;
}

namespace SireMM
{

using SireMol::CGAtomIdx;
using SireMol::AtomIdx;
using SireMol::AtomID;
using SireMol::BondID;
using SireMol::AtomSelection;
using SireMol::AtomMatcher;

/** This class holds a function that acts using the 
    coordinate information of just two atoms */
class SIREMM_EXPORT TwoAtomFunction : public AtomFunction
{

friend SIREMM_EXPORT QDataStream& ::operator<<(QDataStream&, const TwoAtomFunction&);
friend SIREMM_EXPORT QDataStream& ::operator>>(QDataStream&, TwoAtomFunction&);

public:
    TwoAtomFunction();
    TwoAtomFunction(const CGAtomIdx &atom0, const CGAtomIdx &atom1,
                    const SireCAS::Expression &function);
                  
    TwoAtomFunction(const TwoAtomFunction &other);
    
    ~TwoAtomFunction();

    TwoAtomFunction& operator=(const TwoAtomFunction &other);
    
    bool operator==(const TwoAtomFunction &other) const;
    bool operator!=(const TwoAtomFunction &other) const;

    QString toString() const;
    
    const CGAtomIdx& atom0() const;
    const CGAtomIdx& atom1() const;

private:
    /** The indicies of the two atoms */
    CGAtomIdx atm0, atm1;
};

namespace detail
{

class IDPair
{
public:
    IDPair(quint32 atom0=0, quint32 atom1=0);
    IDPair(const IDPair &other);
    
    ~IDPair();
    
    IDPair& operator=(const IDPair &other);
    
    bool operator==(const IDPair &other) const;
    bool operator!=(const IDPair &other) const;
    
    quint32 atom0;
    quint32 atom1;
};

inline uint qHash(const IDPair &idpair)
{
    return (idpair.atom0 << 16) | (idpair.atom1 & 0x0000FFFF);
}

}

/** This class holds the set of TwoAtomFunction potentials that
    act between the atoms in a molecule
    
    @author Christopher Woods
*/
class SIREMM_EXPORT TwoAtomFunctions
        : public SireBase::ConcreteProperty<TwoAtomFunctions,AtomFunctions>
{

friend SIREMM_EXPORT QDataStream& ::operator<<(QDataStream&, const TwoAtomFunctions&);
friend SIREMM_EXPORT QDataStream& ::operator>>(QDataStream&, TwoAtomFunctions&);

public:
    TwoAtomFunctions();
    
    TwoAtomFunctions(const MoleculeData &moldata);
    TwoAtomFunctions(const MoleculeInfoData &molinfo);
    
    TwoAtomFunctions(const TwoAtomFunctions &other);
    
    ~TwoAtomFunctions();
    
    static const char* typeName();
    
    TwoAtomFunctions& operator=(const TwoAtomFunctions &other);
    
    bool operator==(const TwoAtomFunctions &other) const;
    bool operator!=(const TwoAtomFunctions &other) const;

    QString toString() const;

    void set(AtomIdx atom0, AtomIdx atom1,
             const Expression &expression);
             
    void set(const AtomID &atom0, const AtomID &atom1,
             const Expression &expression);
             
    void set(const BondID &bondid, const Expression &expression);

    void clear(AtomIdx atom);
    void clear(const AtomID &atom);

    void clear(AtomIdx atom0, AtomIdx atom1);
    void clear(const AtomID &atom0, const AtomID &atom1);
    void clear(const BondID &bondid);
    
    void clear();

    void substitute(const Identities &identities);

    bool isEmpty() const;
    
    int nFunctions() const;

    Expression potential(AtomIdx atom0, AtomIdx atom1) const;
    Expression potential(const AtomID &atom0, const AtomID &atom1) const;
    Expression potential(const BondID &bondid) const;
    
    Expression force(AtomIdx atom0, AtomIdx atom1, const Symbol &symbol) const;
    Expression force(const AtomID &atom0, const AtomID &atom1,
                     const Symbol &symbol) const;
    Expression force(const BondID &bondid, const Symbol &symbol) const;

    QVector<TwoAtomFunction> potentials() const;
    QVector<TwoAtomFunction> forces(const Symbol &symbol) const;
    
    TwoAtomFunctions includeOnly(const AtomSelection &selection,
                                 bool isstrict=true) const;

protected:
    SireBase::PropertyPtr _pvt_makeCompatibleWith(const MoleculeInfoData &molinfo,
                                                  const AtomMatcher &atommatcher) const;
    SireBase::PropertyPtr _pvt_makeCompatibleWith(const MoleculeInfoData &molinfo,
                                                  const QHash<AtomIdx,AtomIdx> &map) const;
    
private:
    void removeSymbols(QSet<Symbol> symbols);

    /** All of the potential functions, identified by the atom pair
        that contains that function */
    QHash<detail::IDPair,Expression> potentials_by_atoms;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

//////
////// Inline functions of TwoAtomFunction
//////

/** Return the first atom of the pair */
inline const CGAtomIdx& TwoAtomFunction::atom0() const
{
    return atm0;
}

/** Return the second atom of the pair */
inline const CGAtomIdx& TwoAtomFunction::atom1() const
{
    return atm1;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE( SireMM::TwoAtomFunctions );

SIRE_EXPOSE_CLASS( SireMM::TwoAtomFunction )
SIRE_EXPOSE_CLASS( SireMM::TwoAtomFunctions )

SIRE_END_HEADER

#endif
