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

#ifndef SIREMM_FOURATOMFUNCTIONS_H
#define SIREMM_FOURATOMFUNCTIONS_H

#include <QHash>

#include "atomfunctions.h"

#include "SireMol/cgatomidx.h"
#include "SireMol/atomidx.h"
#include "SireMol/atomid.h"

#include "SireMol/dihedralid.h"
#include "SireMol/improperid.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class FourAtomFunction;
class FourAtomFunctions;
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::FourAtomFunction&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::FourAtomFunction&);

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::FourAtomFunctions&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::FourAtomFunctions&);

namespace SireMol
{
class AtomSelection;
}

namespace SireMM
{

using SireMol::CGAtomIdx;
using SireMol::AtomIdx;
using SireMol::AtomID;
using SireMol::DihedralID;
using SireMol::ImproperID;
using SireMol::AtomSelection;
using SireMol::AtomMatcher;

/** This class holds a function that acts using the
    coordinate information of just four atoms */
class SIREMM_EXPORT FourAtomFunction : public AtomFunction
{

friend SIREMM_EXPORT QDataStream& ::operator<<(QDataStream&, const FourAtomFunction&);
friend SIREMM_EXPORT QDataStream& ::operator>>(QDataStream&, FourAtomFunction&);

public:
    FourAtomFunction();
    FourAtomFunction(const CGAtomIdx &atom0, const CGAtomIdx &atom1,
                     const CGAtomIdx &atom2, const CGAtomIdx &atom3,
                     const SireCAS::Expression &function);

    FourAtomFunction(const FourAtomFunction &other);

    ~FourAtomFunction();

    FourAtomFunction& operator=(const FourAtomFunction &other);

    bool operator==(const FourAtomFunction &other) const;
    bool operator!=(const FourAtomFunction &other) const;

    QString toString() const;

    const CGAtomIdx& atom0() const;
    const CGAtomIdx& atom1() const;
    const CGAtomIdx& atom2() const;
    const CGAtomIdx& atom3() const;

private:
    /** The indicies of the four atoms */
    CGAtomIdx atm0, atm1, atm2, atm3;
};

namespace detail
{

class IDQuad
{
public:
    IDQuad(quint32 atom0=0, quint32 atom1=0,
           quint32 atom2=0, quint32 atom3=0);

    IDQuad(const IDQuad &other);

    ~IDQuad();

    IDQuad& operator=(const IDQuad &other);

    bool operator==(const IDQuad &other) const;
    bool operator!=(const IDQuad &other) const;

    bool operator<(const IDQuad &other) const;
    bool operator<=(const IDQuad &other) const;
    bool operator>(const IDQuad &other) const;
    bool operator>=(const IDQuad &other) const;

    quint32 operator[](int i) const
    {
        switch(i)
        {
        case 0:
            return atom0;
        case 1:
            return atom1;
        case 2:
            return atom2;
        default:
            return atom3;
        }
    }

    quint32 atom0;
    quint32 atom1;
    quint32 atom2;
    quint32 atom3;
};

SIRE_ALWAYS_INLINE uint qHash(const IDQuad &idquad)
{
    return (idquad.atom0 << 24) |
           ( (idquad.atom1 << 16) & 0x00FF0000) |
           ( (idquad.atom2 << 8)  & 0x0000FF00) |
           (idquad.atom3 & 0x000000FF);
}

}

/** This class holds the set of FourAtomFunction potentials that
    act between the atoms in a molecule

    @author Christopher Woods
*/
class SIREMM_EXPORT FourAtomFunctions
        : public SireBase::ConcreteProperty<FourAtomFunctions,AtomFunctions>
{

friend SIREMM_EXPORT QDataStream& ::operator<<(QDataStream&, const FourAtomFunctions&);
friend SIREMM_EXPORT QDataStream& ::operator>>(QDataStream&, FourAtomFunctions&);

public:
    FourAtomFunctions();

    FourAtomFunctions(const MoleculeData &moldata);
    FourAtomFunctions(const MoleculeInfoData &molinfo);

    FourAtomFunctions(const FourAtomFunctions &other);

    ~FourAtomFunctions();

    static const char* typeName();

    FourAtomFunctions& operator=(const FourAtomFunctions &other);

    bool operator==(const FourAtomFunctions &other) const;
    bool operator!=(const FourAtomFunctions &other) const;

    QString toString() const;

    int nFunctions() const;

    void set(AtomIdx atom0, AtomIdx atom1,
             AtomIdx atom2, AtomIdx atom3,
             const Expression &expression);

    void set(const AtomID &atom0, const AtomID &atom1,
             const AtomID &atom2, const AtomID &atom3,
             const Expression &expression);

    void set(const DihedralID &dihedralid, const Expression &expression);
    void set(const ImproperID &improperid, const Expression &expression);

    void clear(AtomIdx atom);
    void clear(const AtomID &atom);

    void clear(AtomIdx atom0, AtomIdx atom1,
               AtomIdx atom2, AtomIdx atom3);
    void clear(const AtomID &atom0, const AtomID &atom1,
               const AtomID &atom2, const AtomID &atom3);

    void clear(const DihedralID &dihedralid);
    void clear(const ImproperID &improperid);

    void clear();

    void substitute(const Identities &identities);

    bool isEmpty() const;

    Expression potential(AtomIdx atom0, AtomIdx atom1,
                         AtomIdx atom2, AtomIdx atom3) const;
    Expression potential(const AtomID &atom0, const AtomID &atom1,
                         const AtomID &atom2, const AtomID &atom3) const;

    Expression potential(const DihedralID &dihedralid) const;
    Expression potential(const ImproperID &improperid) const;

    Expression force(AtomIdx atom0, AtomIdx atom1,
                     AtomIdx atom2, AtomIdx atom3,
                     const Symbol &symbol) const;
    Expression force(const AtomID &atom0, const AtomID &atom1,
                     const AtomID &atom2, const AtomID &atom3,
                     const Symbol &symbol) const;

    Expression force(const DihedralID &dihedralid, const Symbol &symbol) const;
    Expression force(const ImproperID &improperid, const Symbol &symbol) const;

    QVector<FourAtomFunction> potentials() const;
    QVector<FourAtomFunction> forces(const Symbol &symbol) const;

    FourAtomFunctions includeOnly(const AtomSelection &selected_atoms,
                                  bool isstrict=true) const;

protected:
    SireBase::PropertyPtr _pvt_makeCompatibleWith(const MoleculeInfoData &molinfo,
                                                  const AtomMatcher &atommatcher) const;
    SireBase::PropertyPtr _pvt_makeCompatibleWith(const MoleculeInfoData &molinfo,
                                                  const QHash<AtomIdx,AtomIdx> &map) const;

private:
    void removeSymbols(QSet<Symbol> symbols);

    /** All of the potential functions, identified by the atom quad
        that contains that function */
    QHash<detail::IDQuad,Expression> potentials_by_atoms;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

//////
////// Inline functions of FourAtomFunction
//////

/** Return the first atom of the quad */
SIRE_ALWAYS_INLINE const CGAtomIdx& FourAtomFunction::atom0() const
{
    return atm0;
}

/** Return the second atom of the quad */
SIRE_ALWAYS_INLINE const CGAtomIdx& FourAtomFunction::atom1() const
{
    return atm1;
}

/** Return the third atom of the quad */
SIRE_ALWAYS_INLINE const CGAtomIdx& FourAtomFunction::atom2() const
{
    return atm2;
}

/** Return the fourth atom of the quad */
SIRE_ALWAYS_INLINE const CGAtomIdx& FourAtomFunction::atom3() const
{
    return atm3;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE( SireMM::FourAtomFunctions );

SIRE_EXPOSE_CLASS( SireMM::FourAtomFunction )
SIRE_EXPOSE_CLASS( SireMM::FourAtomFunctions )

SIRE_END_HEADER

#endif
