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

#ifndef SIREMOL_ATOMID_H
#define SIREMOL_ATOMID_H

#include "SireBase/propertymap.h"

#include "SireID/id.h"

#include "SireID/specify.hpp"
#include "SireID/idorset.hpp"
#include "SireID/idandset.hpp"
#include "SireID/invertmatch.hpp"
#include "SireID/matchall.hpp"

SIRE_BEGIN_HEADER

namespace SireMol
{

using SireBase::PropertyMap;

using SireID::Specify;
using SireID::IDAndSet;
using SireID::IDOrSet;

template<class T>
class Selector;

template<class GROUPID, class ATOMID>
class GroupAtomID;

class Atom;

class MolInfo;

class AtomIdx;
class AtomIdentifier;

class Molecules;
class MoleculeGroup;
class MolGroupsBase;
class MolNum;

class MolID;
class MolAtomID;

class MoleculeView;

class CGID;
class ResID;
class ChainID;
class SegID;

class ResWithAtoms;
class CGsWithAtoms;
class ChainsWithAtoms;
class SegsWithAtoms;

/** This is the base class of all identifiers that are used
    to identify an atom

    @author Christopher Woods
*/
class SIREMOL_EXPORT AtomID : public SireID::ID
{
public:
    typedef AtomIdx Index;
    typedef AtomIdentifier Identifier;
    typedef MolInfo SearchObject;

    AtomID();
    AtomID(const AtomID &other);

    virtual ~AtomID();

    Specify<AtomID> operator[](qint64 i) const;
    Specify<AtomID> operator[](const SireBase::Range &range) const;
    Specify<AtomID> operator()(const SireBase::Range &range) const;
    Specify<AtomID> operator()(qint64 i) const;
    Specify<AtomID> operator()(qint64 start, qint64 end) const;
    Specify<AtomID> operator()(qint64 start, qint64 end, qint64 increment) const;

    IDAndSet<AtomID> operator+(const AtomID &other) const;

    GroupAtomID<CGID,AtomID> operator+(const CGID &other) const;
    GroupAtomID<ResID,AtomID> operator+(const ResID &other) const;
    GroupAtomID<ChainID,AtomID> operator+(const ChainID &other) const;
    GroupAtomID<SegID,AtomID> operator+(const SegID &other) const;
    MolAtomID operator+(const MolID &other) const;

    SireID::InvertMatch<AtomID> operator-() const;

    IDAndSet<AtomID> operator-(const AtomID &other) const;

    GroupAtomID<CGID,AtomID> operator-(const CGID &other) const;
    GroupAtomID<ResID,AtomID> operator-(const ResID &other) const;
    GroupAtomID<ChainID,AtomID> operator-(const ChainID &other) const;
    GroupAtomID<SegID,AtomID> operator-(const SegID &other) const;

    IDOrSet<AtomID> operator*(const AtomID &other) const;

    IDOrSet<AtomID> operator*(const CGID &other) const;
    IDOrSet<AtomID> operator*(const ResID &other) const;
    IDOrSet<AtomID> operator*(const ChainID &other) const;
    IDOrSet<AtomID> operator*(const SegID &other) const;
    IDOrSet<AtomID> operator*(const MolID &other) const;

    IDAndSet<AtomID> operator&&(const AtomID &other) const;
    IDAndSet<AtomID> operator&(const AtomID &other) const;

    GroupAtomID<CGID,AtomID> operator&&(const CGID &other) const;
    GroupAtomID<ResID,AtomID> operator&&(const ResID &other) const;
    GroupAtomID<ChainID,AtomID> operator&&(const ChainID &other) const;
    GroupAtomID<SegID,AtomID> operator&&(const SegID &other) const;
    MolAtomID operator&&(const MolID &other) const;

    GroupAtomID<CGID,AtomID> operator&(const CGID &other) const;
    GroupAtomID<ResID,AtomID> operator&(const ResID &other) const;
    GroupAtomID<ChainID,AtomID> operator&(const ChainID &other) const;
    GroupAtomID<SegID,AtomID> operator&(const SegID &other) const;
    MolAtomID operator&(const MolID &other) const;

    IDOrSet<AtomID> operator||(const AtomID &other) const;
    IDOrSet<AtomID> operator|(const AtomID &other) const;

    IDOrSet<AtomID> operator|(const CGID &other) const;
    IDOrSet<AtomID> operator|(const ResID &other) const;
    IDOrSet<AtomID> operator|(const ChainID &other) const;
    IDOrSet<AtomID> operator|(const SegID &other) const;
    IDOrSet<AtomID> operator|(const MolID &other) const;

    IDOrSet<AtomID> operator||(const CGID &other) const;
    IDOrSet<AtomID> operator||(const ResID &other) const;
    IDOrSet<AtomID> operator||(const ChainID &other) const;
    IDOrSet<AtomID> operator||(const SegID &other) const;
    IDOrSet<AtomID> operator||(const MolID &other) const;

    SireID::InvertMatch<AtomID> operator!() const;

    ResWithAtoms residues() const;
    CGsWithAtoms cutGroups() const;
    ChainsWithAtoms chains() const;
    SegsWithAtoms segments() const;

    static const char* typeName()
    {
        return "SireMol::AtomID";
    }

    virtual AtomID* clone() const=0;

    /** Map this ID back to the indicies of the matching atoms in the molecule,
        using the passed MoleculeInfo to do the mapping */
    virtual QList<AtomIdx> map(const MolInfo &molinfo) const=0;

    virtual QList<AtomIdx> map(const MoleculeView &molview,
                               const PropertyMap &map = PropertyMap()) const;

    virtual Atom selectFrom(const MoleculeView &molview,
                            const PropertyMap &map = PropertyMap()) const;

    virtual Selector<Atom> selectAllFrom(const MoleculeView &molview,
                                         const PropertyMap &map = PropertyMap()) const;

    virtual Atom selectFrom(const Molecules &molecules,
                            const PropertyMap &map = PropertyMap()) const;

    virtual QHash< MolNum,Selector<Atom> >
                selectAllFrom(const Molecules &molecules,
                              const PropertyMap &map = PropertyMap()) const;

    virtual Atom selectFrom(const MoleculeGroup &molgroup,
                            const PropertyMap &map = PropertyMap()) const;

    virtual QHash< MolNum,Selector<Atom> >
                selectAllFrom(const MoleculeGroup &molgroup,
                              const PropertyMap &map = PropertyMap()) const;

    virtual Atom selectFrom(const MolGroupsBase &molgroups,
                            const PropertyMap &map = PropertyMap()) const;
    virtual QHash< MolNum,Selector<Atom> >
                selectAllFrom(const MolGroupsBase &molgroups,
                              const PropertyMap &map = PropertyMap()) const;

    SireID::InvertMatch<AtomID> invert() const;
    SireID::InvertMatch<AtomID> inverse() const;

    static SireID::MatchAll<AtomID> any();

    static AtomIdentifier fromString(const QString &id);

protected:
    static QList<AtomIdx> matchAll(const MolInfo &molinfo);

    void processMatches(QList<AtomIdx> &matches, const MolInfo &molinfo) const;
};

}

#include "atomidentifier.h"

SIRE_EXPOSE_CLASS( SireMol::AtomID )
SIRE_EXPOSE_ALIAS( SireID::Specify<SireMol::AtomID>, SireMol::Specify_AtomID_ )
SIRE_EXPOSE_ALIAS( SireID::IDAndSet<SireMol::AtomID>, SireMol::IDAndSet_AtomID_ )
SIRE_EXPOSE_ALIAS( SireID::IDOrSet<SireMol::AtomID>, SireMol::IDOrSet_AtomID_ )
SIRE_EXPOSE_ALIAS( SireID::MatchAll<SireMol::AtomID>, SireMol::MatchAll_AtomID_ )
SIRE_EXPOSE_ALIAS( SireID::InvertMatch<SireMol::AtomID>, SireMol::InvertMatch_AtomID_ )

#ifdef SIRE_INSTANTIATE_TEMPLATES
template class SireID::Specify<SireMol::AtomID>;
template class SireID::IDAndSet<SireMol::AtomID>;
template class SireID::IDOrSet<SireMol::AtomID>;
template class SireID::MatchAll<SireMol::AtomID>;
template class SireID::InvertMatch<SireMol::AtomID>;
#endif

SIRE_END_HEADER

#endif
