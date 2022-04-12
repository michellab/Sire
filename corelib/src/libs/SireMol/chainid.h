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

#ifndef SIREMOL_CHAINID_H
#define SIREMOL_CHAINID_H

#include "SireID/id.h"

#include "atomsin.hpp"
#include "resin.hpp"

#include "SireID/idandset.hpp"
#include "SireID/idorset.hpp"
#include "SireID/specify.hpp"
#include "SireID/matchall.hpp"
#include "SireID/invertmatch.hpp"

SIRE_BEGIN_HEADER

namespace SireMol
{

using SireID::IDAndSet;
using SireID::IDOrSet;
using SireID::Specify;

template<class T>
class Selector;

template<class GROUPID, class ATOMID>
class GroupAtomID;

template<class G0, class G1>
class GroupGroupID;

class ChainResID;

class MolInfo;

class ChainIdx;
class ChainIdentifier;

class Chain;

class Molecules;
class MoleculeGroup;
class MolGroupsBase;
class MolNum;

/** This is the base class of all identifiers that are used
    to identify a chain within a molecule

    @author Christopher Woods
*/
class SIREMOL_EXPORT ChainID : public SireID::ID
{

public:
    typedef ChainIdx Index;
    typedef ChainIdentifier Identifier;
    typedef MolInfo SearchObject;

    ChainID();
    ChainID(const ChainID &other);

    virtual ~ChainID();

    Specify<ChainID> operator[](qint64 i) const;
    Specify<ChainID> operator[](const SireBase::Range &range) const;
    Specify<ChainID> operator()(const SireBase::Range &range) const;
    Specify<ChainID> operator()(qint64 i) const;
    Specify<ChainID> operator()(qint64 start, qint64 end) const;
    Specify<ChainID> operator()(qint64 start, qint64 end, qint64 increment) const;

    IDAndSet<ChainID> operator+(const ChainID &other) const;
    ChainResID operator+(const ResID &other) const;
    GroupAtomID<ChainID,AtomID> operator+(const AtomID &other) const;
    GroupGroupID<SegID,ChainID> operator+(const SegID &other) const;
    GroupGroupID<CGID,ChainID> operator+(const CGID &other) const;

    IDAndSet<ChainID> operator-(const ChainID &other) const;
    ChainResID operator-(const ResID &other) const;
    GroupAtomID<ChainID,AtomID> operator-(const AtomID &other) const;
    GroupGroupID<SegID,ChainID> operator-(const SegID &other) const;
    GroupGroupID<CGID,ChainID> operator-(const CGID &other) const;

    SireID::InvertMatch<ChainID> operator-() const;

    IDAndSet<ChainID> operator&&(const ChainID &other) const;
    ChainResID operator&&(const ResID &other) const;
    GroupAtomID<ChainID,AtomID> operator&&(const AtomID &other) const;
    GroupGroupID<SegID,ChainID> operator&&(const SegID &other) const;
    GroupGroupID<CGID,ChainID> operator&&(const CGID &other) const;

    IDAndSet<ChainID> operator&(const ChainID &other) const;
    ChainResID operator&(const ResID &other) const;
    GroupAtomID<ChainID,AtomID> operator&(const AtomID &other) const;
    GroupGroupID<SegID,ChainID> operator&(const SegID &other) const;
    GroupGroupID<CGID,ChainID> operator&(const CGID &other) const;

    IDOrSet<ChainID> operator*(const ChainID &other) const;
    IDOrSet<ChainID> operator||(const ChainID &other) const;
    IDOrSet<ChainID> operator|(const ChainID &other) const;

    IDOrSet<ResID> operator*(const ResID &other) const;
    IDOrSet<ResID> operator||(const ResID &other) const;
    IDOrSet<ResID> operator|(const ResID &other) const;

    IDOrSet<AtomID> operator*(const AtomID &other) const;
    IDOrSet<AtomID> operator||(const AtomID &other) const;
    IDOrSet<AtomID> operator|(const AtomID &other) const;

    SireID::InvertMatch<ChainID> operator!() const;

    AtomsIn<ChainID> atoms() const;
    AtomsIn<ChainID> atom(int i) const;
    AtomsIn<ChainID> atoms(int i, int j) const;

    ResIn<ChainID> residues() const;
    ResIn<ChainID> residue(int i) const;
    ResIn<ChainID> residues(int i, int j) const;

    static ChainIdentifier fromString(const QString &id);

    static SireID::MatchAll<ChainID> any();

    SireID::InvertMatch<ChainID> invert() const;
    SireID::InvertMatch<ChainID> inverse() const;

    static const char* typeName()
    {
        return "SireMol::ChainID";
    }

    virtual ChainID* clone() const=0;

    /** Map this ID back to the indicies of the chains in the molecule,
        using the passed MoleculeInfo to do the mapping */
    virtual QList<ChainIdx> map(const MolInfo &molinfo) const=0;


    virtual QList<ChainIdx> map(const MoleculeView &molview,
                                const PropertyMap &map = PropertyMap()) const;

    virtual Chain selectFrom(const MoleculeView &molview,
                             const PropertyMap &map = PropertyMap()) const;

    virtual Selector<Chain> selectAllFrom(const MoleculeView &molview,
                                          const PropertyMap &map = PropertyMap()) const;

    virtual Chain selectFrom(const Molecules &molecules,
                             const PropertyMap &map = PropertyMap()) const;

    virtual QHash< MolNum,Selector<Chain> >
                selectAllFrom(const Molecules &molecules,
                              const PropertyMap &map = PropertyMap()) const;

    virtual Chain selectFrom(const MoleculeGroup &molgroup,
                             const PropertyMap &map = PropertyMap()) const;

    virtual QHash< MolNum,Selector<Chain> >
                selectAllFrom(const MoleculeGroup &molgroup,
                              const PropertyMap &map = PropertyMap()) const;

    virtual Chain selectFrom(const MolGroupsBase &molgroups,
                             const PropertyMap &map = PropertyMap()) const;
    virtual QHash< MolNum,Selector<Chain> >
                selectAllFrom(const MolGroupsBase &molgroups,
                              const PropertyMap &map = PropertyMap()) const;

protected:
    static QList<ChainIdx> matchAll(const MolInfo &molinfo);

    void processMatches(QList<ChainIdx> &matches, const MolInfo &molinfo) const;

};

}

#include "chainidentifier.h"

SIRE_EXPOSE_CLASS( SireMol::ChainID )
SIRE_EXPOSE_ALIAS( SireID::Specify<SireMol::ChainID>, SireMol::Specify_ChainID_ )
SIRE_EXPOSE_ALIAS( SireMol::AtomsIn<SireMol::ChainID>, SireMol::AtomsIn_ChainID_ )
SIRE_EXPOSE_ALIAS( SireMol::ResIn<SireMol::ChainID>, SireMol::ResIn_ChainID_ )
SIRE_EXPOSE_ALIAS( SireID::IDAndSet<SireMol::ChainID>, SireMol::IDAndSet_ChainID_ )
SIRE_EXPOSE_ALIAS( SireID::IDOrSet<SireMol::ChainID>, SireMol::IDOrSet_ChainID_ )
SIRE_EXPOSE_ALIAS( SireID::MatchAll<SireMol::ChainID>, SireMol::MatchAll_ChainID_ )
SIRE_EXPOSE_ALIAS( SireID::InvertMatch<SireMol::ChainID>, SireMol::InvertMatch_ChainID_ )

#ifdef SIRE_INSTANTIATE_TEMPLATES
template class SireID::Specify<SireMol::ChainID>;
template class SireMol::AtomsIn<SireMol::ChainID>;
template class SireMol::ResIn<SireMol::ChainID>;
template class SireID::IDAndSet<SireMol::ChainID>;
template class SireID::IDOrSet<SireMol::ChainID>;
template class SireID::MatchAll<SireMol::ChainID>;
template class SireID::InvertMatch<SireMol::ChainID>;
#endif

SIRE_END_HEADER

#endif
