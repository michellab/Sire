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

#ifndef SIREMOL_RESID_H
#define SIREMOL_RESID_H

#include "SireID/id.h"

#include "atomsin.hpp"

#include "SireID/idandset.hpp"
#include "SireID/idorset.hpp"
#include "SireID/specify.hpp"

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

class ResIdx;
class ResIdentifier;

class Residue;

class Molecules;
class MoleculeGroup;
class MolGroupsBase;
class MolNum;

class ChainsWithRes;

/** This is the base class of all identifiers that are used 
    to identify a residue within a molecule

    @author Christopher Woods
*/
class SIREMOL_EXPORT ResID : public SireID::ID
{

public:
    typedef ResIdx Index;
    typedef ResIdentifier Identifier;
    typedef MolInfo SearchObject;
    
    ResID();

    ResID(const ResID &other);

    virtual ~ResID();
    
    Specify<ResID> operator[](int i) const;
    Specify<ResID> operator()(int i) const;
    Specify<ResID> operator()(int i, int j) const;

    IDAndSet<ResID> operator+(const ResID &other) const;
    ChainResID operator+(const ChainID &other) const;
    GroupAtomID<ResID,AtomID> operator+(const AtomID &other) const;
    GroupGroupID<SegID,ResID> operator+(const SegID &other) const;
    GroupGroupID<CGID,ResID> operator+(const CGID &other) const;

    IDAndSet<ResID> operator&&(const ResID &other) const;
    ChainResID operator&&(const ChainID &other) const;
    GroupAtomID<ResID,AtomID> operator&&(const AtomID &other) const;
    GroupGroupID<SegID,ResID> operator&&(const SegID &other) const;
    GroupGroupID<CGID,ResID> operator&&(const CGID &other) const;

    IDAndSet<ResID> operator&(const ResID &other) const;
    ChainResID operator&(const ChainID &other) const;
    GroupAtomID<ResID,AtomID> operator&(const AtomID &other) const;
    GroupGroupID<SegID,ResID> operator&(const SegID &other) const;
    GroupGroupID<CGID,ResID> operator&(const CGID &other) const;

    IDOrSet<ResID> operator*(const ResID &other) const;
    IDOrSet<ResID> operator||(const ResID &other) const;
    IDOrSet<ResID> operator|(const ResID &other) const;

    AtomsIn<ResID> atoms() const;
    AtomsIn<ResID> atom(int i) const;
    AtomsIn<ResID> atoms(int i, int j) const;
    
    ChainsWithRes chains() const;
    
    static const char* typeName()
    {
        return "SireMol::ResID";
    }
    
    virtual ResID* clone() const=0;

    /** Map this ID back to the indicies of the residues in the molecule, 
        using the passed MoleculeInfo to do the mapping */
    virtual QList<ResIdx> map(const MolInfo &molinfo) const=0;

    virtual QList<ResIdx> map(const MoleculeView &molview,
                              const PropertyMap &map = PropertyMap()) const;
    
    virtual Residue selectFrom(const MoleculeView &molview,
                               const PropertyMap &map = PropertyMap()) const;
    
    virtual Selector<Residue> selectAllFrom(const MoleculeView &molview,
                                         const PropertyMap &map = PropertyMap()) const;
    
    virtual Residue selectFrom(const Molecules &molecules,
                               const PropertyMap &map = PropertyMap()) const;
                            
    virtual QHash< MolNum,Selector<Residue> >
                selectAllFrom(const Molecules &molecules,
                              const PropertyMap &map = PropertyMap()) const;

    virtual Residue selectFrom(const MoleculeGroup &molgroup,
                               const PropertyMap &map = PropertyMap()) const;
                            
    virtual QHash< MolNum,Selector<Residue> >
                selectAllFrom(const MoleculeGroup &molgroup,
                              const PropertyMap &map = PropertyMap()) const;
    
    virtual Residue selectFrom(const MolGroupsBase &molgroups,
                               const PropertyMap &map = PropertyMap()) const;
    virtual QHash< MolNum,Selector<Residue> > 
                selectAllFrom(const MolGroupsBase &molgroups,
                              const PropertyMap &map = PropertyMap()) const;

protected:
    void processMatches(QList<ResIdx> &matches, const MolInfo &molinfo) const;
};

}

#include "residentifier.h"

SIRE_EXPOSE_CLASS( SireMol::ResID )
SIRE_EXPOSE_ALIAS( SireID::Specify<SireMol::ResID>, SireMol::Specify_ResID_ )
SIRE_EXPOSE_ALIAS( SireMol::AtomsIn<SireMol::ResID>, SireMol::AtomsIn_ResID_ )
SIRE_EXPOSE_ALIAS( SireID::IDAndSet<SireMol::ResID>, SireMol::IDAndSet_ResID_ )
SIRE_EXPOSE_ALIAS( SireID::IDOrSet<SireMol::ResID>, SireMol::IDOrSet_ResID_ )

#ifdef SIRE_INSTANTIATE_TEMPLATES
template class SireID::Specify<SireMol::ResID>;
template class SireMol::AtomsIn<SireMol::ResID>;
template class SireID::IDAndSet<SireMol::ResID>;
template class SireID::IDOrSet<SireMol::ResID>;
#endif

SIRE_END_HEADER

#endif
