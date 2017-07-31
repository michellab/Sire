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

#ifndef SIREMOL_CGID_H
#define SIREMOL_CGID_H

#include "SireID/id.h"

#include "atomsin.hpp"

#include "SireID/idandset.hpp"
#include "SireID/idorset.hpp"
#include "SireID/specify.hpp"
#include "SireID/invertmatch.hpp"
#include "SireID/matchall.hpp"

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

template<class GROUP0, class GROUP1>
class GroupGroupID;

class MolInfo;

class CGIdx;
class CGIdentifier;

class CutGroup;

class Molecules;
class MoleculeGroup;
class MolGroupsBase;
class MolNum;

/** This is the base class of all identifiers that are used 
    to identify a CutGroup

    @author Christopher Woods
*/
class SIREMOL_EXPORT CGID : public SireID::ID
{

public:
    typedef CGIdx Index;
    typedef CGIdentifier Identifier;
    typedef MolInfo SearchObject;

    CGID();
    CGID(const CGID &other);

    virtual ~CGID();

    Specify<CGID> operator[](int i) const;
    Specify<CGID> operator()(int i) const;
    Specify<CGID> operator()(int i, int j) const;
    
    IDAndSet<CGID> operator+(const CGID &other) const;
    GroupAtomID<CGID,AtomID> operator+(const AtomID &other) const;
    GroupGroupID<SegID,CGID> operator+(const SegID &other) const;
    GroupGroupID<CGID,ChainID> operator+(const ChainID &other) const;
    GroupGroupID<CGID,ResID> operator+(const ResID &other) const;

    IDAndSet<CGID> operator&&(const CGID &other) const;
    GroupAtomID<CGID,AtomID> operator&&(const AtomID &other) const;
    GroupGroupID<SegID,CGID> operator&&(const SegID &other) const;
    GroupGroupID<CGID,ChainID> operator&&(const ChainID &other) const;
    GroupGroupID<CGID,ResID> operator&&(const ResID &other) const;

    IDAndSet<CGID> operator&(const CGID &other) const;
    GroupAtomID<CGID,AtomID> operator&(const AtomID &other) const;
    GroupGroupID<SegID,CGID> operator&(const SegID &other) const;
    GroupGroupID<CGID,ChainID> operator&(const ChainID &other) const;
    GroupGroupID<CGID,ResID> operator&(const ResID &other) const;
    
    IDOrSet<CGID> operator*(const CGID &other) const;
    IDOrSet<AtomID> operator*(const AtomID &other) const;

    IDOrSet<CGID> operator||(const CGID &other) const;
    IDOrSet<AtomID> operator||(const AtomID &other) const;

    IDOrSet<CGID> operator|(const CGID &other) const;
    IDOrSet<AtomID> operator|(const AtomID &other) const;
    
    SireID::InvertMatch<CGID> operator!() const;
    
    AtomsIn<CGID> atoms() const;
    AtomsIn<CGID> atom(int i) const;
    AtomsIn<CGID> atoms(int i, int j) const;
    
    static SireID::MatchAll<CGID> any();
    SireID::InvertMatch<CGID> invert() const;
    SireID::InvertMatch<CGID> inverse() const;
    
    static const char* typeName()
    {
        return "SireMol::CGID";
    }
    
    virtual CGID* clone() const=0;

    /** Map this ID back to the indicies of the CutGroups
        within the molecule described by the info in 'molinfo' */
    virtual QList<CGIdx> map(const MolInfo &molinfo) const=0;


    virtual QList<CGIdx> map(const MoleculeView &molview,
                             const PropertyMap &map = PropertyMap()) const;
    
    virtual CutGroup selectFrom(const MoleculeView &molview,
                                const PropertyMap &map = PropertyMap()) const;
    
    virtual Selector<CutGroup> selectAllFrom(const MoleculeView &molview,
                                         const PropertyMap &map = PropertyMap()) const;
    
    virtual CutGroup selectFrom(const Molecules &molecules,
                                const PropertyMap &map = PropertyMap()) const;
                            
    virtual QHash< MolNum,Selector<CutGroup> >
                selectAllFrom(const Molecules &molecules,
                              const PropertyMap &map = PropertyMap()) const;

    virtual CutGroup selectFrom(const MoleculeGroup &molgroup,
                                const PropertyMap &map = PropertyMap()) const;
                            
    virtual QHash< MolNum,Selector<CutGroup> >
                selectAllFrom(const MoleculeGroup &molgroup,
                              const PropertyMap &map = PropertyMap()) const;
    
    virtual CutGroup selectFrom(const MolGroupsBase &molgroups,
                               const PropertyMap &map = PropertyMap()) const;
    virtual QHash< MolNum,Selector<CutGroup> > 
                selectAllFrom(const MolGroupsBase &molgroups,
                              const PropertyMap &map = PropertyMap()) const;

protected:
    static QList<CGIdx> matchAll(const MolInfo &molinfo);

    void processMatches(QList<CGIdx> &matches, const MolInfo &molinfo) const;

};

}

#include "cgidentifier.h"

SIRE_EXPOSE_CLASS( SireMol::CGID )
SIRE_EXPOSE_ALIAS( (SireID::Specify<SireMol::CGID>), SireMol::Specify_CGID_ )
SIRE_EXPOSE_ALIAS( (SireMol::AtomsIn<SireMol::CGID>), SireMol::AtomsIn_CGID_ )
SIRE_EXPOSE_ALIAS( SireID::IDAndSet<SireMol::CGID>, SireMol::IDAndSet_CGID_ )
SIRE_EXPOSE_ALIAS( SireID::IDOrSet<SireMol::CGID>, SireMol::IDOrSet_CGID_ )
SIRE_EXPOSE_ALIAS( SireID::MatchAll<SireMol::CGID>, SireMol::MatchAll_CGID_ )
SIRE_EXPOSE_ALIAS( SireID::InvertMatch<SireMol::CGID>, SireMol::InvertMatch_CGID_ )

#ifdef SIRE_INSTANTIATE_TEMPLATES
template class SireID::Specify<SireMol::CGID>;
template class SireMol::AtomsIn<SireMol::CGID>;
template class SireID::IDAndSet<SireMol::CGID>;
template class SireID::IDOrSet<SireMol::CGID>;
template class SireID::InvertMatch<SireMol::CGID>;
template class SireID::MatchAll<SireMol::CGID>;
#endif

SIRE_END_HEADER

#endif
