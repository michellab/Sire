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

#include "cgid.h"
#include "cgidentifier.h"

#include "groupatomids.h"
#include "groupgroupids.h"

#include "atom.h"
#include "selector.hpp"

#include "mover.hpp"
#include "editor.hpp"

#include "molinfo.h"

#include "partialmolecule.h"
#include "segment.h"
#include "chain.h"
#include "residue.h"
#include "cutgroup.h"
#include "atom.h"

#include "molecules.h"
#include "moleculegroup.h"
#include "moleculegroups.h"

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

#include "tostring.h"

using namespace SireMol;
using namespace SireID;
using namespace SireStream;

/** Constructor */
CGID::CGID() : ID()
{}

/** Copy constructor */
CGID::CGID(const CGID &other) : ID(other)
{}

/** Destructor */
CGID::~CGID()
{}
  
/** Return a specific atom that matches this ID */
Specify<CGID> CGID::operator[](int i) const
{
    return Specify<CGID>(*this, i);
}

/** Return a specific atom that matches this ID */
Specify<CGID> CGID::operator()(int i) const
{
    return this->operator[](i);
}

/** Return a range of atoms that match this ID */
Specify<CGID> CGID::operator()(int i, int j) const
{
    return Specify<CGID>(*this, i, j);
}

/** Combine with another ID */
IDAndSet<CGID> CGID::operator+(const CGID &other) const
{
    return IDAndSet<CGID>(*this, other);
}

/** Syntactic sugar for operator+ */
IDAndSet<CGID> CGID::operator&&(const CGID &other) const
{
    return this->operator+(other);
}

/** Syntactic sugar for operator+ */
IDAndSet<CGID> CGID::operator&(const CGID &other) const
{
    return this->operator+(other);
}

/** Combine with another ID */
GroupAtomID<CGID,AtomID> CGID::operator+(const AtomID &other) const
{
    return GroupAtomID<CGID,AtomID>(*this, other);
}

/** Syntactic sugar for operator+ */
GroupAtomID<CGID,AtomID> CGID::operator&&(const AtomID &other) const
{
    return this->operator+(other);
}

/** Syntactic sugar for operator+ */
GroupAtomID<CGID,AtomID> CGID::operator&(const AtomID &other) const
{
    return this->operator+(other);
}

/** Combine with another ID */
GroupGroupID<SegID,CGID> CGID::operator+(const SegID &other) const
{
    return GroupGroupID<SegID,CGID>(other, *this);
}

/** Syntactic sugar for operator+ */
GroupGroupID<SegID,CGID> CGID::operator&&(const SegID &other) const
{
    return this->operator+(other);
}

/** Syntactic sugar for operator+ */
GroupGroupID<SegID,CGID> CGID::operator&(const SegID &other) const
{
    return this->operator+(other);
}

/** Combine with another ID */
GroupGroupID<CGID,ChainID> CGID::operator+(const ChainID &other) const
{
    return GroupGroupID<CGID,ChainID>(*this, other);
}

/** Syntactic sugar for operator+ */
GroupGroupID<CGID,ChainID> CGID::operator&&(const ChainID &other) const
{
    return this->operator+(other);
}

/** Syntactic sugar for operator+ */
GroupGroupID<CGID,ChainID> CGID::operator&(const ChainID &other) const
{
    return this->operator+(other);
}

/** Combine with another ID */
GroupGroupID<CGID,ResID> CGID::operator+(const ResID &other) const
{
    return GroupGroupID<CGID,ResID>(*this, other);
}

/** Syntactic sugar for operator+ */
GroupGroupID<CGID,ResID> CGID::operator&&(const ResID &other) const
{
    return this->operator+(other);
}

/** Syntactic sugar for operator+ */
GroupGroupID<CGID,ResID> CGID::operator&(const ResID &other) const
{
    return this->operator+(other);
}

/** Return the combination of this ID or 'other' */
IDOrSet<CGID> CGID::operator*(const CGID &other) const
{
    return IDOrSet<CGID>(*this, other);
}

/** Syntactic sugar for operator* */
IDOrSet<CGID> CGID::operator||(const CGID &other) const
{
    return this->operator*(other);
}

/** Syntactic sugar for operator* */
IDOrSet<CGID> CGID::operator|(const CGID &other) const
{
    return this->operator*(other);
}

/** Return the combination of this ID or 'other' */
IDOrSet<AtomID> CGID::operator*(const AtomID &other) const
{
    return other * *this;
}

/** Syntactic sugar for operator* */
IDOrSet<AtomID> CGID::operator||(const AtomID &other) const
{
    return this->operator*(other);
}

/** Syntactic sugar for operator* */
IDOrSet<AtomID> CGID::operator|(const AtomID &other) const
{
    return this->operator*(other);
}

/** Return a match for any cutgroup */
MatchAll<CGID> CGID::any()
{
    return MatchAll<CGID>();
}

/** Return the inverse of this match */
InvertMatch<CGID> CGID::invert() const
{
    return InvertMatch<CGID>(*this);
}

/** Return the inverse of this match */
InvertMatch<CGID> CGID::inverse() const
{
    return this->invert();
}

/** Return the inverse of this match */
InvertMatch<CGID> CGID::operator!() const
{
    return this->invert();
}

/** Return all of the matching objects */
QList<CGIdx> CGID::matchAll(const MolInfo &molinfo)
{
    return molinfo.getCutGroups();
}

/** Return the atoms in the matching residues */
AtomsIn<CGID> CGID::atoms() const
{
    return AtomsIn<CGID>(*this);
}

/** Return a specific atom in the matching residues */
AtomsIn<CGID> CGID::atom(int i) const
{
    return AtomsIn<CGID>(*this, i);
}

/** Return a range of atoms in the matching residues */
AtomsIn<CGID> CGID::atoms(int i, int j) const
{
    return AtomsIn<CGID>(*this, i, j);
}

void CGID::processMatches(QList<CGIdx> &matches, const MolInfo&) const
{
    if (matches.isEmpty())
        throw SireMol::missing_cutgroup( QObject::tr(
            "There are no CutGroups that match the ID \"%1\" in the passed molecule.")
                .arg(this->toString()), CODELOC );

    qSort(matches);
}

/** Map this CGICutGroupD to the CutGroups in the passed molecule view

    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
*/
QList<CGIdx> CGID::map(const MoleculeView &molview, const PropertyMap&) const
{
    QList<CGIdx> cgidxs = this->map( molview.data().info() );
    
    if (molview.selectedAll())
        return cgidxs;
    else 
    {
        QMutableListIterator<CGIdx> it(cgidxs);
        
        const AtomSelection selected_cgs = molview.selection();
        
        while (it.hasNext())
        {
            it.next();
            
            if (not selected_cgs.selected(it.value()))
                it.remove();
        }
        
        if (cgidxs.isEmpty())
            throw SireMol::missing_cutgroup( QObject::tr(
                    "No CutGroups matching %1 can be found in the passed molecule.")
                        .arg(this->toString()), CODELOC );
                        
        return cgidxs;
    }
}

/** Select the CutGroup from the passed view that matches this ID

    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
    \throw SireMol::duplicate_cutgroup
*/
CutGroup CGID::selectFrom(const MoleculeView &molview, const PropertyMap &map) const
{
    QList<CGIdx> cgidxs = this->map(molview, map);
    
    if (cgidxs.count() > 1)
        throw SireMol::duplicate_cutgroup( QObject::tr(
                "More than one CutGroup matches the ID %1 (CutGroups %2).")
                    .arg(this->toString()).arg(Sire::toString(cgidxs)),
                        CODELOC );
                        
    return CutGroup(molview.data(), cgidxs.at(0)); 
}

/** Select all the CutGroups from the passed view that match this ID

    \throw SireMol::missing_cutgroup
    \throw SireError::invalid_index
    \throw SireMol::duplicate_cutgroup
*/
Selector<CutGroup> CGID::selectAllFrom(const MoleculeView &molview, 
                                       const PropertyMap &map) const
{
    QList<CGIdx> cgidxs = this->map(molview, map);
                        
    return Selector<CutGroup>(molview.data(), cgidxs); 
}

/** Return all of the CutGroups from the 'molecules' that match
    this ID
    
    \throw SireMol::missing_cutgroup
*/
QHash< MolNum,Selector<CutGroup> >
CGID::selectAllFrom(const Molecules &molecules, const PropertyMap &map) const
{
    QHash< MolNum,Selector<CutGroup> > selected_cgs;
    
    //loop over all molecules...
    for (Molecules::const_iterator it = molecules.constBegin();
         it != molecules.constEnd();
         ++it)
    {
        try
        {
            //try to find this CutGroup in this molecule
            selected_cgs.insert( it.key(), this->selectAllFrom(*it,map) );
        }
        catch(...)
        {}
    }
    
    if (selected_cgs.isEmpty())
        throw SireMol::missing_cutgroup( QObject::tr(
            "There was no atom matching the ID \"%1\" in "
            "the set of molecules.")
                .arg(this->toString()), CODELOC );

    return selected_cgs;
}

/** Return the atom from the molecules 'molecules' that matches
    this ID
    
    \throw SireMol::missing_cutgroup
    \throw SireMol::duplicate_cutgroup
*/
CutGroup CGID::selectFrom(const Molecules &molecules, const PropertyMap &map) const
{
    QHash< MolNum,Selector<CutGroup> > mols = this->selectAllFrom(molecules, map);
    
    if (mols.count() > 1)
        throw SireMol::duplicate_cutgroup( QObject::tr(
            "More than one molecule contains an atom that "
            "matches this ID (%1). These molecules have numbers %2.")
                .arg(this->toString()).arg(Sire::toString(mols.keys())),
                    CODELOC );
                    
    const Selector<CutGroup> &atoms = *(mols.constBegin());
    
    if (atoms.count() > 1)
        throw SireMol::duplicate_cutgroup( QObject::tr(
            "While only one molecule (MolNum == %1) "
            "contains an atom that matches this ID (%2), it contains "
            "more than one atom that matches.")
                .arg(atoms.data().number()).arg(this->toString()),
                    CODELOC );
                    
    return atoms(0);
}

/** Return the atom from the molecule group 'molgroup' that matches
    this ID
    
    \throw SireMol::missing_cutgroup
    \throw SireMol::duplicate_cutgroup
*/
CutGroup CGID::selectFrom(const MoleculeGroup &molgroup, 
                        const PropertyMap &map) const
{
    return this->selectFrom(molgroup.molecules(), map);
}

/** Return the atoms from the molecule group 'molgroup' that match
    this ID
    
    \throw SireMol::missing_cutgroup
*/
QHash< MolNum,Selector<CutGroup> >
CGID::selectAllFrom(const MoleculeGroup &molgroup,
                      const PropertyMap &map) const
{
    return this->selectAllFrom(molgroup.molecules(), map);
}

/** Return the atom from the molecule groups 'molgroups' that matches 
    this ID
    
    \throw SireMol::missing_cutgroup
    \throw SireMol::duplicate_cutgroup
*/
CutGroup CGID::selectFrom(const MolGroupsBase &molgroups,
                        const PropertyMap &map) const
{
    return this->selectFrom(molgroups.molecules(), map);
}

/** Return the set of atoms that match this ID in the molecule groups
    set 'molgroups' 
    
    \throw SireMol::missing_cutgroup
*/
QHash< MolNum,Selector<CutGroup> >
CGID::selectAllFrom(const MolGroupsBase &molgroups,
                      const PropertyMap &map) const
{
    return this->selectAllFrom(molgroups.molecules(), map);
}

//fully instantiate Specify<CGID> and AtomsIn<CGID>
namespace SireID
{
    template class Specify<CGID>;
    template class IDAndSet<CGID>;
    template class IDOrSet<CGID>;
}

namespace SireMol
{
    template class AtomsIn<CGID>;
}

static const RegisterMetaType< Specify<CGID> > r_specify_cgid;
static const RegisterMetaType< AtomsIn<CGID> > r_atomsin_cgid;
static const RegisterMetaType< IDAndSet<CGID> > r_idandset_cgid;
static const RegisterMetaType< IDOrSet<CGID> > r_idorset_cgid;
